#!/usr/bin/env snakemake
# Cleaned, executable version: removed duplicate/ambiguous rules (finalize & second all),
# keep a single final product results/done produced by raw_stats_full.

import os
import sys
import concurrent.futures
import pandas as pd
import seaborn as sns
import matplotlib
matplotlib.use("agg")
import matplotlib.pyplot as plt
from snakemake.utils import min_version
min_version("7.0")

########################################
# Final target
########################################
rule all:
    input:
        "results/done"

########################################
# Helper functions
########################################
def parse(stats_file):
    if os.path.exists(stats_file):
        try:
            df = pd.read_csv(stats_file, sep="\t")
        except pd.errors.EmptyDataError:
            print(f"{stats_file} is empty, please check")
            return None
        if not df.empty:
            return df
        return None
    else:
        print(f"{stats_file} is not exists")
        return None


def merge(input_list, func, workers, **kwargs):
    df_list = []
    with concurrent.futures.ProcessPoolExecutor(max_workers=workers) as executor:
        for df in executor.map(func, input_list):
            if df is not None:
                df_list.append(df)
    df_ = pd.concat(df_list)
    if "output" in kwargs:
        df_.to_csv(kwargs["output"], sep="\t", index=False)
    return df_


def parse_samples(samples_tsv):
    samples_df = pd.read_csv(
        samples_tsv, sep="\t",
        dtype={"sample_id": str, "fq1": str, "fq2": str}
    ).set_index(["sample_id"])

    cancel = False
    if "fq1" in samples_df.columns:
        for sample_id in samples_df.index.unique():
            sample_id = str(sample_id)
            if "." in sample_id:
                print(f"{sample_id} contain '.', please remove '.', now quitting :)")
                cancel = True
            fq1_list = samples_df.loc[[sample_id]]["fq1"].dropna().tolist()
            fq2_list = samples_df.loc[[sample_id]]["fq2"].dropna().tolist()
            for fq_file in fq1_list:
                if not fq_file.endswith(".gz"):
                    print(f"{fq_file} need gzip format")
                    cancel = True
                # (Optional existence checks can be re-enabled)
    else:
        print(f"wrong header: {samples_df.columns}")
        cancel = True

    if cancel:
        sys.exit(-1)
    return samples_df


def update_qcstats_row(row):
    rows = row["file"].split(".")
    sample_id = rows[0]
    step = rows[1]
    fq_type = "pe"
    if (rows[2] != "1") and (rows[2] != "2"):
        fq_type = "se"
    reads = "fq1"
    if ".2.fq.gz" in row["file"]:
        reads = "fq2"
    return pd.Series([sample_id, reads, step, fq_type])


def update_qcstats_df(df):
    stats_df = df.copy()
    stats_df[["sample_id", "reads", "step", "fq_type"]] = stats_df.apply(
        lambda row: update_qcstats_row(row), axis=1, result_type="expand")
    return stats_df


def qc_summary_merge(index_df, summary_df, md5=False, host_rate=True, **kwargs):
    if host_rate:
        df_host = summary_df.loc[:, ["sample_id", "host_rate"]].drop_duplicates()

    df_len = summary_df.loc[:, ["sample_id", "reads", "sum_len", "step"]] \
        .pivot(index=["sample_id", "step"], columns="reads", values="sum_len") \
        .reset_index() \
        .rename(columns={"fq1": "sumlen_fq1", "fq2": "sumlen_fq2"})

    if md5:
        df_md5 = summary_df.loc[:, ["sample_id", "reads", "md5", "step"]] \
            .pivot(index=["sample_id", "step"], columns="reads", values="md5") \
            .reset_index() \
            .rename(columns={"fq1": "md5_fq1", "fq2": "md5_fq2"})

    df_gc = summary_df.loc[:, ["sample_id", "reads", "GC(%)", "step"]] \
        .pivot(index=["sample_id", "step"], columns="reads", values="GC(%)") \
        .reset_index() \
        .rename(columns={"fq1": "GC_fq1(%)", "fq2": "GC_fq2(%)"})

    df_stats = summary_df.loc[:, [
        "sequencing_batch", "sample_id", "format", "type", "step", "fq_type",
        "num_seqs", "sum_len", "min_len", "avg_len", "max_len",
        "Q1", "Q2", "Q3", "sum_gap", "Q20(%)", "Q30(%)"
    ]]

    df_stats_ = df_stats.groupby(
        ["sequencing_batch", "sample_id", "format", "type", "step", "fq_type"]
    ).agg(
        num_seqs=("num_seqs", "sum"),
        sum_len=("sum_len", "sum"),
        min_len=("min_len", "min"),
        avg_len=("avg_len", "mean"),
        max_len=("max_len", "max"),
        Q1=("Q1", "mean"),
        Q2=("Q2", "mean"),
        Q3=("Q3", "mean"),
        Q20_per=("Q20(%)", "mean"),
        Q30_per=("Q30(%)", "mean")
    ).reset_index().rename(columns={"Q20_per": "Q20(%)", "Q30_per": "Q30(%)"})

    for c in ["avg_len", "Q1", "Q2", "Q3", "Q20(%)", "Q30(%)"]:
        df_stats_[c] = df_stats_[c].round(2)

    df_summary = pd.merge(df_stats_, df_len, how="inner", on=["sample_id", "step"])
    if md5:
        df_summary = pd.merge(df_summary, df_md5, how="inner", on=["sample_id", "step"])
    df_summary = pd.merge(df_summary, df_gc, how="inner", on=["sample_id", "step"])
    if host_rate:
        df_summary = pd.merge(df_summary, df_host, how="inner", on=["sample_id"])
    df_summary = pd.merge(df_summary, index_df, how="outer", on=["sample_id"])

    if "Lane" in list(df_summary.columns):
        df_summary = df_summary.sort_values(["Lane", "Sample_Project", "sum_len", "sample_id"])
    else:
        df_summary = df_summary.sort_values(["Sample_Project", "sum_len", "sample_id"])

    df_summary_undetermined = df_summary.query('sample_id=="Undetermined"')
    unclassified_base_number = 0
    if not df_summary_undetermined.empty:
        unclassified_base_number = float(df_summary_undetermined["sum_len"])

    if "output_tsv" in kwargs:
        df_summary.to_csv(kwargs["output_tsv"], sep="\t", index=False)

    if "output_xlsx" in kwargs:
        writer = pd.ExcelWriter(kwargs["output_xlsx"], engine="xlsxwriter")
        workbook = writer.book
        df_summary.to_excel(writer, sheet_name="Raw_Stats", index=False)
        red_format = workbook.add_format(
            {"bg_color": "#FFC7CE", "font_color": "#9C0006", "bold": True}
        )
        red_format.set_font_size(18)
        red_format.set_underline()
        ws = writer.sheets["Raw_Stats"]
        ws.set_tab_color("#fb8072")
        ws.set_zoom(120)
        rows = len(df_summary) + 1
        ws.conditional_format(
            f"B2:B{rows}",
            {"type": "cell", "criteria": "<=", "value": 5000000000, "format": red_format}
        )
        ws.conditional_format(
            f"H2:H{rows}",
            {"type": "cell", "criteria": "<=", "value": 5000000000, "format": red_format}
        )
        ws.set_column("A:A", 30)
        ws.set_column("G:G", 20)
        ws.write(f'G{rows+2}', 'Unclassified_base_ratio')
        ws.write_formula(f'H{rows+2}', f'={unclassified_base_number}/SUM(H2:H{rows})')
        writer.close()


def qc_bar_plot(df, engine, stacked=False, **kwargs):
    if engine == "seaborn":
        f, ax = plt.subplots(figsize=(10, 7))
        df_ = df.query('reads=="fq1"')
        sns.barplot(x="sample_id", y="num_seqs", hue="step", data=df_)
        ax.set_xticklabels(ax.get_xticklabels(), rotation=-90)
    elif engine == "pandas":
        if not stacked:
            df_ = (
                df.query('reads=="fq1"')
                .pivot(index="sample_id", columns="step", values="num_seqs")
                .loc[:, ["raw", "trimming", "rmhost"]]
            )
            df_.plot(kind="bar", figsize=(10, 7))
        else:
            dict_ = {"sample_id": [], "clean": [], "rmhost": [], "trim": []}
            df = df.set_index("id")
            for i in df.index.unique():
                reads_total = reads_trimmed = reads_host = reads_clean = 0
                if not df.loc[[i],].query('reads=="fq1" and step=="raw"').empty:
                    reads_total = df.loc[[i],].query('reads=="fq1" and step=="raw"')["num_seqs"][0]
                if not df.loc[[i],].query('reads=="fq1" and step=="trimming"').empty:
                    reads_trim = df.loc[[i],].query('reads=="fq1" and step=="trimming"')["num_seqs"][0]
                if not df.loc[[i],].query('reads=="fq1" and step=="rmhost"').empty:
                    reads_clean = df.loc[[i],].query('reads=="fq1" and step=="rmhost"')["num_seqs"][0]
                reads_trimmed = reads_total - reads_trim
                reads_host = reads_trim - reads_clean
                dict_["sample_id"].append(i)
                dict_["trim"].append(reads_trimmed)
                dict_["rmhost"].append(reads_host)
                dict_["clean"].append(reads_clean)
            df_ = pd.DataFrame(dict_).sort_values("sample_id").set_index("sample_id")
            df_.plot(kind="bar", stacked=True, color=["#2ca02c", "#ff7f0e", "#1f77b4"], figsize=(10, 7))
    plt.xlabel("Sample ID")
    plt.ylabel("The number of reads(-pair)")
    plt.title("Fastq quality control barplot", fontsize=11)
    if "output" in kwargs:
        plt.savefig(kwargs["output"])


########################################
# Load samples
########################################
SAMPLES = parse_samples(config["samples"])

########################################
# Rules
########################################
rule prepare_reads:
    input:
        r1 = lambda wc: SAMPLES.loc[wc.sample, "fq1"],
        r2 = lambda wc: SAMPLES.loc[wc.sample, "fq2"]
    output:
        r1 = "results/00.raw/reads/{sample}/{sample}.raw.1.fq.gz",
        r2 = "results/00.raw/reads/{sample}/{sample}.raw.2.fq.gz"
    shell:
        r'''
        R1=$(realpath {input.r1})
        R2=$(realpath {input.r2})
        mkdir -p $(dirname {output.r1})
        ln -s $R1 {output.r1}
        ln -s $R2 {output.r2}
        '''

rule raw_fastqc:
    input:
        r1 = "results/00.raw/reads/{sample}/{sample}.raw.1.fq.gz",
        r2 = "results/00.raw/reads/{sample}/{sample}.raw.2.fq.gz"
    output:
        directory("results/00.raw/fastqc/{sample}")
    log:
        "logs/00.raw_fastqc/{sample}.raw_fastqc.log"
    benchmark:
        "benchmark/00.raw_fastqc/{sample}.raw_fastqc.benchmark.txt"
    threads: 8
    shell:
        r'''
        mkdir -p {output}
        fastqc \
          --outdir {output} \
          --threads {threads} \
          --format fastq \
          {input.r1} {input.r2} \
          >{log} 2>&1
        '''

rule raw_fastqc_multiqc:
    input:
        expand("results/00.raw/fastqc/{sample}", sample=SAMPLES.index.unique())
    output:
        expand("results/02.multiqc/multiqc_fastqc/{batch}-fastqc_multiqc_report.html",
               batch=[config["batch"]])
    log:
        "logs/02.raw_fastqc_multiqc/raw_fastqc_multiqc.log"
    benchmark:
        "benchmark/02.raw_fastqc_multiqc/raw_fastqc_multiqc.benchmark.txt"
    params:
        batch = config["batch"]
    threads: 4
    shell:
        r'''
        OUTDIR=$(dirname {output})
        rm -rf $OUTDIR
        multiqc \
          --outdir $OUTDIR \
          --title {params.batch}-fastqc \
          --module fastqc \
          {input} \
          >{log} 2>&1
        '''

rule raw_fastp:
    input:
        r1 = "results/00.raw/reads/{sample}/{sample}.raw.1.fq.gz",
        r2 = "results/00.raw/reads/{sample}/{sample}.raw.2.fq.gz"
    output:
        json = "results/00.raw/fastp/{sample}.fastp.json",
        html = "results/00.raw/fastp/{sample}.fastp.html"
    log:
        "logs/00.raw_fastp/{sample}.raw_fastp.log"
    benchmark:
        "benchmark/00.raw_fastp/{sample}.raw_fastp.benchmark.txt"
    threads: 8
    shell:
        r'''
        fastp \
          --in1 {input.r1} \
          --in2 {input.r2} \
          --thread {threads} \
          --json {output.json} \
          --html {output.html} \
          >{log} 2>&1
        '''

rule raw_fastp_multiqc:
    input:
        expand("results/00.raw/fastp/{sample}.fastp.json", sample=SAMPLES.index.unique())
    output:
        expand("results/02.multiqc/multiqc_fastp/{batch}-fastp_multiqc_report.html",
               batch=[config["batch"]])
    log:
        "logs/02.raw_fastp_multiqc/raw_fastp_multiqc.log"
    benchmark:
        "benchmark/02.raw_fastp_multiqc/raw_fastp_multiqc.benchmark.txt"
    params:
        batch = config["batch"]
    threads: 4
    shell:
        r'''
        OUTDIR=$(dirname {output})
        rm -rf $OUTDIR
        multiqc \
          --outdir $OUTDIR \
          --title {params.batch}-fastp \
          --module fastp \
          {input} \
          >{log} 2>&1
        '''

rule raw_report_stats:
    input:
        r1 = "results/00.raw/reads/{sample}/{sample}.raw.1.fq.gz",
        r2 = "results/00.raw/reads/{sample}/{sample}.raw.2.fq.gz"
    output:
        "results/00.raw/stats/{sample}_raw_stats.tsv"
    log:
        "logs/00.raw_report_stats/{sample}.raw_report_stats.log"
    benchmark:
        "benchmark/00.raw_report_stats/{sample}.raw_report_stats.benchmark.txt"
    params:
        fq_encoding = "sanger"
    threads: 4
    shell:
        r'''
        seqkit stats \
          --all \
          --basename \
          --tabular \
          --fq-encoding {params.fq_encoding} \
          --out-file {output} \
          --threads {threads} \
          {input.r1} {input.r2} \
          >{log} 2>&1
        '''

rule raw_report_stats_merge:
    input:
        expand("results/00.raw/stats/{sample}_raw_stats.tsv", sample=SAMPLES.index.unique())
    output:
        "results/03.qcreport/raw_stats.tsv"
    log:
        "logs/03.raw_report_stats_merge/raw_report_stats_merge.log"
    benchmark:
        "benchmark/03.raw_report_stats_merge/raw_report_stats_merge.txt"
    threads: 4
    shell:
        r'''
        head -1 {input[0]} > {output}
        for report in {input}; do
          tail -q -n +2 "$report" >> {output} 2>>{log}
        done
        '''

rule raw_report_md5:
    input:
        r1 = "results/00.raw/reads/{sample}/{sample}.raw.1.fq.gz",
        r2 = "results/00.raw/reads/{sample}/{sample}.raw.2.fq.gz"
    output:
        "results/00.raw/md5/{sample}.md5"
    log:
        "logs/00.raw_report_md5/{sample}.log"
    benchmark:
        "benchmark/00.raw_report_md5/{sample}.txt"
    shell:
        r'''
        echo -e "md5\tfile" > {output}
        md5sum {input.r1} | awk -F'[ /]' '{{print $1 "\t" $NF}}' >> {output}
        md5sum {input.r2} | awk -F'[ /]' '{{print $1 "\t" $NF}}' >> {output}
        '''

rule raw_report_md5_merge:
    input:
        expand("results/00.raw/md5/{sample}.md5", sample=SAMPLES.index.unique())
    output:
        "results/03.qcreport/raw_md5.tsv"
    log:
        "logs/03.raw_report_md5_merge/raw_report_md5_merge.log"
    benchmark:
        "benchmark/03.raw_report_md5_merge/raw_report_md5_merge.txt"
    threads: 4
    shell:
        r'''
        head -1 {input[0]} > {output}
        for report in {input}; do
          tail -q -n +2 "$report" >> {output}
        done
        '''

rule raw_report_merge:
    input:
        stats = "results/03.qcreport/raw_stats.tsv",
        md5 = "results/03.qcreport/raw_md5.tsv"
    output:
        summary = "results/03.qcreport/raw_seqstats_summary_l.tsv"
    params:
        sequencing_batch = config["batch"]
    run:
        stats_df = pd.read_csv(input.stats, sep="\t")
        md5_df = pd.read_csv(input.md5, sep="\t")
        df = pd.merge(stats_df, md5_df, on=["file"])
        df = update_qcstats_df(df)
        df["sequencing_batch"] = params.sequencing_batch
        df.to_csv(output.summary, sep="\t", index=False)

rule raw_report_summary:
    input:
        index_file = config["index"],
        summary = "results/03.qcreport/raw_seqstats_summary_l.tsv"
    output:
        summary_tsv = "results/03.qcreport/raw_seqstats_summary_w.tsv",
        summary_xlsx = "results/03.qcreport/raw_seqstats_summary_w.xlsx"
    run:
        index_df = pd.read_csv(input.index_file, skiprows=18).rename(columns={"Sample_ID": "sample_id"})
        summary_df = pd.read_csv(input.summary, sep="\t")
        qc_summary_merge(
            index_df,
            summary_df,
            md5=True,
            host_rate=False,
            output_tsv=output.summary_tsv,
            output_xlsx=output.summary_xlsx
        )

rule raw_report_summary_plot:
    input:
        summary = "results/03.qcreport/raw_seqstats_summary_l.tsv"
    output:
        plot = "results/03.qcreport/raw_seqstats_summary.png"
    run:
        df = pd.read_csv(input.summary, sep="\t")
        qc_bar_plot(df, "seaborn", stacked=False, output=output.plot)

rule raw_stats_basic:
    input:
        summary = "results/03.qcreport/raw_seqstats_summary_l.tsv",
        summary_tsv = "results/03.qcreport/raw_seqstats_summary_w.tsv",
        summary_xlsx = "results/03.qcreport/raw_seqstats_summary_w.xlsx",
        plot = "results/03.qcreport/raw_seqstats_summary.png"
    output:
        summary = expand("results/03.qcreport/{batch}-raw_seqstats_summary.xlsx", batch=[config["batch"]]),
        done = "results/raw_stats_basic_done"
    shell:
        r'''
        cp {input.summary_xlsx} {output.summary}
        touch {output.done}
        '''

rule raw_stats_full:
    input:
        raw_stats_basic_done = "results/raw_stats_basic_done",
        multiqc_fastqc = expand("results/02.multiqc/multiqc_fastqc/{batch}-fastqc_multiqc_report.html",
                                batch=[config["batch"]]),
        multiqc_fastp = expand("results/02.multiqc/multiqc_fastp/{batch}-fastp_multiqc_report.html",
                               batch=[config["batch"]])
    output:
        "results/done"
    shell:
        r'''
        touch {output}
        '''