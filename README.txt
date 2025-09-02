# bcl2fastq Demultiplex + FASTQ 质控工作流

本仓库提供一个基于 **Snakemake** 的简单工作流：  
1. 依赖外部完成 `bcl2fastq` 拆分 (demultiplex) 得到原始 FASTQ。  
2. 标准化软链接整理输入文件。  
3. 对原始测序数据进行基础质控 (FastQC / fastp 统计 / seqkit 汇总 / MD5 校验 / MultiQC 汇总)。  
4. 生成多种 TSV / Excel / 图表汇总报表。  

> 当前版本只对“原始(raw)”数据做统计展示，未实际执行剪切(trimming)或去宿主过滤(host removal)。代码中为后续扩展（如 `trimming`, `rmhost`）预留了列结构。

---

## 目录
- [运行环境与依赖](#运行环境与依赖)
- [输入准备流程总览](#输入准备流程总览)
- [样本表 (`samples.tsv`)](#样本表-samplestsv)
- [配置文件 (`workflow_config.yaml`)](#配置文件-workflow_configyaml)
- [Illumina Index / SampleSheet 说明 (`index` 文件)](#illumina-index--samplesheet-说明-index-文件)
- [主要规则与数据流](#主要规则与数据流)
- [运行方式](#运行方式)
- [输出文件结构](#输出文件结构)
- [结果文件说明](#结果文件说明)
- [可视化与报表](#可视化与报表)
- [FAQ / 常见问题](#faq--常见问题)
- [后续可扩展方向](#后续可扩展方向)
- [License](#license)

---

## 运行环境与依赖

推荐使用 Conda / Mamba 管理环境：

```bash
mamba create -n demux_qc -c conda-forge -c bioconda \
  snakemake fastqc fastp seqkit multiqc \
  pandas seaborn matplotlib xlsxwriter
conda activate demux_qc
```

必要软件：
- Snakemake ≥ 7.0
- fastqc
- fastp（当前仅获取 JSON / HTML 指标，不产出过滤后 FASTQ）
- seqkit
- multiqc
- Python 库：pandas, seaborn, matplotlib, xlsxwriter
- 其他：GNU coreutils, awk, md5sum

---

## 输入准备流程总览

1. 使用 Illumina Experiment Manager 生成 `SampleSheet.csv`。  
2. 运行 `bcl2fastq`（或 `bcl-convert` 也可自行适配），例如：
   ```bash
   nohup bcl2fastq \
     -i /PATH/Run/Data/Intensities/BaseCalls \
     -R /PATH/Run \
     -o /PATH/fastq_output \
     --sample-sheet /PATH/SampleSheet.csv \
     --tiles s_2 > demultiplex.log 2>&1 &
   ```
3. 使用脚本 `scripts_make_samples_tsv.sh`（自行编辑其中 `FASTQ_DIR`）生成 `samples.tsv`：
   ```bash
   bash scripts_make_samples_tsv.sh
   ```
4. 准备 `workflow_config.yaml`（指定 `samples`, `batch`, `index` 等）。  
5. 运行 Snakemake 工作流。  

---

## 样本表 (`samples.tsv`)

制表符分隔，最少包含列：

| sample_id | fq1 | fq2 |
|-----------|-----|-----|

示例：

```text
sample_id	fq1	fq2
S01	/abs/path/FASTQ/S01_R1.fastq.gz	/abs/path/FASTQ/S01_R2.fastq.gz
S02	/abs/path/FASTQ/S02_R1.fastq.gz	/abs/path/FASTQ/S02_R2.fastq.gz
```

约束与校验：
- `sample_id` 不允许包含 `.`（脚本中会检测并退出）。
- `fq1/fq2` 必须以 `.gz` 结尾（需为 gzip 压缩）。
- 若为单端数据，可在脚本中扩展；当前规则假设双端（PE）。

---

## 配置文件 (`workflow_config.yaml`)

最小示例：

```yaml
samples: "config/samples.tsv"          # 指向上述样本表
batch: "Run2025Q3_42"                  # 一批次/项目标识；用于命名报告
index: "config/SampleSheet.csv"        # Illumina SampleSheet 或包含 Sample 信息的表
```

可放在仓库 `config/` 目录下。运行时通过 `--configfile` 指定。

---

## Illumina Index / SampleSheet 说明 (`index` 文件)

- 在规则 `raw_report_summary` 中读取：`pd.read_csv(index_file, skiprows=18)`  
  - 说明：默认 Illumina 标准 `SampleSheet.csv` 前 18 行为 `[Header]`、`[Reads]`、`[Settings]` 等区块，真正的样本表格从第 19 行开始。
- 需要至少存在列 `Sample_ID`（会被重命名为 `sample_id`）。
- 若包含 `Lane`, `Sample_Project` 等列，将在最终汇总表中用于排序和分组。
- `Undetermined` 条目（如果存在）将用于计算未分类碱基比例（写入 Excel 末尾公式）。

若你的 `SampleSheet` 结构不同，需手动调整：
- `skiprows=18` 参数；
- 列名映射。

---

## 主要规则与数据流

按执行顺序（部分合并/汇总类规则可能并行）：

1. `prepare_reads`  
   - 软链接原始 FASTQ 到标准命名：`{sample}.raw.1.fq.gz / .raw.2.fq.gz`
2. `raw_fastqc` → `raw_fastqc_multiqc`  
   - 每个样本 FastQC；MultiQC 汇总。
3. `raw_fastp` → `raw_fastp_multiqc`  
   - fastp 仅生成 JSON/HTML 质控报告（未输出修剪 FASTQ）。
4. `raw_report_stats` → `raw_report_stats_merge`  
   - 用 `seqkit stats` 统计每个样本双端文件 → 合并。
5. `raw_report_md5` → `raw_report_md5_merge`  
   - 计算每个样本 R1/R2 MD5 → 合并。
6. `raw_report_merge`  
   - 合并 stats + md5，生成长格式 `raw_seqstats_summary_l.tsv`，添加推断列：
     - `sample_id`, `reads (fq1/fq2)`, `step (raw)`, `fq_type (pe)`
7. `raw_report_summary`  
   - 读取 `index` 文件并宽表透视：生成
     - `raw_seqstats_summary_w.tsv`
     - `raw_seqstats_summary_w.xlsx`
8. `raw_report_summary_plot`  
   - 基于长表画柱状图 `raw_seqstats_summary.png`
9. `raw_stats_basic`  
   - 复制生成 `{batch}-raw_seqstats_summary.xlsx`（批次命名）。
10. `raw_stats_full`  
    - 聚合最终 MultiQC 报告 + 汇总 Excel，生成标志文件 `results/done`。
11. `all`  
    - Snakemake 终极目标：`results/done`

---

## 运行方式

基础运行：

```bash
snakemake -j 12 \
  --configfile workflow_config.yaml \
  --printshellcmds \
  --rerun-incomplete \
  --latency-wait 60
```

后台（Linux）：

```bash
nohup snakemake -j 12 \
  --configfile workflow_config.yaml \
  --printshellcmds \
  --rerun-incomplete \
  --latency-wait 60 > result.log 2>&1 &
```

启用 Conda（若在 Snakefile 中未来为各规则添加 `conda:` 环境）：

```bash
snakemake -j 12 --use-conda --conda-frontend mamba ...
```

Dry-run 查看计划：

```bash
snakemake -n
```

查看 DAG：

```bash
snakemake --dag | dot -Tpng > dag.png
```

---

## 输出文件结构

示例（两例样本：S01, S02；`batch = Run2025Q3_42`）：

```
results/
├── 00.raw
│   ├── reads
│   │   ├── S01
│   │   │   ├── S01.raw.1.fq.gz -> /abs/FASTQ/S01_R1.fastq.gz
│   │   │   └── S01.raw.2.fq.gz -> /abs/FASTQ/S01_R2.fastq.gz
│   │   └── S02
│   │       ├── S02.raw.1.fq.gz
│   │       └── S02.raw.2.fq.gz
│   ├── fastqc
│   │   ├── S01
│   │   │   ├── S01.raw.1_fastqc.html
│   │   │   ├── S01.raw.1_fastqc.zip
│   │   │   ├── S01.raw.2_fastqc.html
│   │   │   └── S01.raw.2_fastqc.zip
│   │   └── S02/...
│   ├── fastp
│   │   ├── S01.fastp.html
│   │   ├── S01.fastp.json
│   │   ├── S02.fastp.html
│   │   └── S02.fastp.json
│   ├── stats
│   │   ├── S01_raw_stats.tsv
│   │   └── S02_raw_stats.tsv
│   └── md5
│       ├── S01.md5
│       └── S02.md5
├── 02.multiqc
│   ├── multiqc_fastqc
│   │   ├── Run2025Q3_42-fastqc_multiqc_report.html
│   │   └── multiqc_data/...
│   └── multiqc_fastp
│       ├── Run2025Q3_42-fastp_multiqc_report.html
│       └── multiqc_data/...
├── 03.qcreport
│   ├── raw_stats.tsv
│   ├── raw_md5.tsv
│   ├── raw_seqstats_summary_l.tsv
│   ├── raw_seqstats_summary_w.tsv
│   ├── raw_seqstats_summary_w.xlsx
│   ├── raw_seqstats_summary.png
│   └── Run2025Q3_42-raw_seqstats_summary.xlsx
├── raw_stats_basic_done
└── done
logs/
│   ├── 00.raw_fastqc/*.log
│   ├── 00.raw_fastp/*.log
│   ├── 00.raw_report_stats/*.log
│   ├── 00.raw_report_md5/*.log
│   ├── 02.raw_fastqc_multiqc/*.log
│   ├── 02.raw_fastp_multiqc/*.log
│   └── 03.raw_report_stats_merge/*.log
benchmark/
    ├── 00.raw_fastqc/*.benchmark.txt
    ├── 00.raw_fastp/*.benchmark.txt
    ├── 00.raw_report_stats/*.benchmark.txt
    ├── 00.raw_report_md5/*.txt
    ├── 02.raw_fastqc_multiqc/*.benchmark.txt
    ├── 02.raw_fastp_multiqc/*.benchmark.txt
    └── 03.raw_report_stats_merge/*.benchmark.txt
```

> 注：`fastqc` 下每个样本目录内容随 FastQC 版本可能有细微差异。`multiqc_data/` 包含解析后的 JSON / TSV。

---

## 结果文件说明

| 文件 / 目录 | 说明 |
|-------------|------|
| `results/00.raw/reads/*` | 标准化命名的双端 FASTQ 软链接。 |
| `results/00.raw/fastqc/{sample}` | FastQC 原始质量报告（HTML + ZIP）。 |
| `results/00.raw/fastp/*.json/html` | fastp 指标报告（未输出修剪 FASTQ）。 |
| `results/00.raw/stats/{sample}_raw_stats.tsv` | `seqkit stats` 统计。 |
| `results/00.raw/md5/{sample}.md5` | 每个样本双端文件 MD5 值。 |
| `results/03.qcreport/raw_stats.tsv` | 所有样本 `seqkit` 合并表。 |
| `results/03.qcreport/raw_md5.tsv` | 所有样本 MD5 合并表。 |
| `results/03.qcreport/raw_seqstats_summary_l.tsv` | 长格式（含 reads=fq1/fq2、step、fq_type 等）。 |
| `results/03.qcreport/raw_seqstats_summary_w.tsv` | 宽格式（fq1/fq2 各列透视）。 |
| `results/03.qcreport/raw_seqstats_summary_w.xlsx` | 同上 Excel，含条件格式 & 未分类比例公式。 |
| `results/03.qcreport/raw_seqstats_summary.png` | 柱状图（当前仅 raw 层）。 |
| `results/03.qcreport/{batch}-raw_seqstats_summary.xlsx` | 以 batch 命名的最终交付汇总。 |
| `results/02.multiqc/*_multiqc_report.html` | FastQC / fastp MultiQC 汇总页面。 |
| `results/done` | 工作流最终完成标志。 |

Excel 中主要字段：
- 基础统计：`num_seqs`, `sum_len`, `min_len`, `avg_len`, `max_len`
- 质量分位：`Q1`, `Q2`, `Q3`
- 质量比率：`Q20(%)`, `Q30(%)`
- GC 比例：`GC_fq1(%)`, `GC_fq2(%)`
- 长度汇总：`sumlen_fq1`, `sumlen_fq2`
- MD5：`md5_fq1`, `md5_fq2`
- 可能来自 SampleSheet：`Lane`, `Sample_Project` 等

---

## 可视化与报表

当前仅提供：
- `raw_seqstats_summary.png`：每个样本 raw 层 reads (fq1) 的 `num_seqs` 分布（Seaborn 柱状图）。
- Excel 条件格式：对低于设定阈值（示例中 5,000,000,000）标红，可按需调整 Snakefile 中对应 `conditional_format` 的条件值。

如需：
- 堆叠图展示修剪 / 去宿主损失情况 → 需先实现对应 `trimming` / `rmhost` 规则并在文件命名中包含步骤关键字。
- 箱线图 / GC 分布图 → 可在 `raw_report_summary_plot` 中扩展。

---

## FAQ / 常见问题

1. Q: 为什么只有 `raw` 步骤，图里没有 trimming / rmhost？
   A: 当前工作流尚未添加实际剪切或去宿主规则，命名与函数仅为将来扩展预留。

2. Q: 运行时报 `sample_id contain '.'`？
   A: 请移除样本名中的点号（`.`），否则下游解析 `sample.step.reads` 命名结构会冲突。

3. Q: Excel 最后两行的未分类比例是什么？
   A: 若 `sample_id == "Undetermined"` 存在，则计算 `Undetermined sum_len / 所有样本 sum_len`。如果不需要，可在 `qc_summary_merge` 中删除相关代码。

4. Q: 如何添加修剪（如 fastp 输出过滤后 FASTQ）？
   A: 增加新规则输出 `{sample}.trimming.1.fq.gz` / `.2.fq.gz`，并更新统计/MD5 规则输入，再在合并逻辑中识别 `file` 名的第二段为 `trimming`。

5. Q: 多线程设置？
   A: 各规则 `threads:` 指定，实际并发由 `snakemake -j` 控制。

---

## 后续可扩展方向

- 添加 fastp 真实剪切输出并替换后续统计输入。
- 加入宿主过滤（如 bowtie2 / kraken2）并生成 `rmhost` 层数据。
- 增加 FastQC/fastp 指标二次解析进总表（如 per base quality）。
- 支持单端数据（调整 `parse_samples` + 命名约定）。
- 通过 `--cluster` / `--profile` 集成 HPC / SLURM / LSF。
- 提供 `conda:` / `container:` 字段保证可重复环境。
- 将统计绘图升级为交互式（Plotly / R markdown）。

---

## License

根据项目实际需要选择（如 MIT / Apache-2.0）。若未指定，建议添加 `LICENSE` 文件。

---

## 致谢

- Illumina 工具链（bcl2fastq）
- 开源社区：Snakemake / fastqc / fastp / seqkit / MultiQC / pandas / seaborn

如有问题或改进建议，欢迎提 Issue 或 PR（若公开仓库）。

---

（完）