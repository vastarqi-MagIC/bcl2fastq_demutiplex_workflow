# bcl2fastq_demutiplex_workflow

This is the demultiplex workflow to generate fastq file from BCL file and perform QC for the fastq

First use Illumina Experiment Manager to set up the samplesheet.

Use Illumina bcl2fastq to demultiplex

 nohup bcl2fastq -i /sequence_file_path/Data/Intensities/BaseCalls/ -R sequence_file_path -o /output_path/result --sample-sheet /samplesheet_path/samples_sheet.csv --tiles s_2 > ./demultiplex.log 2>&1 &

Use make scripts_make_samples_tsv.sh to generate the samples.tsv
bash scripts_make_samples_tsv.sh
note that change the FASTQ_DIR in the scripts_make_samples_tsv.sh

and change the parameters in the workflow_config.yaml

then

nohp snakemake --cores 12 --configfile workflow_config.yaml --printshellcmds --rerun-incomplete --latency-wait 60 > ./result.log 2>&1 &

the output should look like:

