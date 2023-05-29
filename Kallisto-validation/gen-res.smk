#!/bin/python3
# Author: Haoliang Xue

# nohup command:
#       nohup snakemake --snakefile gen-res.smk --cluster "qsub -q common -l nodes=node13:ppn=3" --jobs 5 -p --latency-wait 60 --rerun-incomplete >> log-gen-res.txt &

# Involved programs
SRATOOLS_BIN = "/home/haoliang.xue/.conda/envs/sra-tools/bin/"
CUTADAPT = "/home/haoliang.xue/.conda/envs/cutadapt/bin/cutadapt"
FASTQC = "/home/haoliang.xue/.conda/envs/multiqc/bin/fastqc"
MULTIQC = "/home/haoliang.xue/.conda/envs/multiqc/bin/multiqc"
RSCRIPT = "/home/haoliang.xue/.conda/envs/dekupl-hlx/bin/Rscript"
REINDEER_IMG = "/store/EQUIPES/SSFA/haoliang-shared/Reindeer.sif"
COUNTTAGS = "/home/haoliang.xue/tools/countTags/bin/countTags"

# Dataset
GTF_PATH = "/store/EQUIPES/SSFA/MEMBERS/haoliang.xue/reindeer_appli/data/Homo_sapiens.GRCh38.99.gtf"
DATA_DIR = "/store/EQUIPES/SSFA/Data/CCLE/"
KALLISTO_RES = DATA_DIR + "new-rmAdapt/kallisto-quant/" # Provided by Chloé & Benoit
KMERATOR_DIR = DATA_DIR + "kmerator-spec-ctgs/" # Provided by Chloé & Benoit

# Outputs
RES_DIR = "/store/EQUIPES/SSFA/MEMBERS/haoliang.xue/reindeer-appli/CCLE12-gene/"
RAW_FQ_DIR = RES_DIR + "raw-fq/"
TRIM_FQ_DIR = RES_DIR + "trimmed-fq/"
TXIMPORT_RES = RES_DIR + "gene-counts/"
REINDEER_QRY = RES_DIR + "reindeer-qres/"
COUNTTAGS_RES = RES_DIR + "countTags-res/"

# Wildcards
SMP_LIST = ["SRR8615893", "SRR8615897", "SRR8615898", "SRR8615899", "SRR8615900", "SRR8615901", 
            "SRR8615904", "SRR8615905", "SRR8615944", "SRR8616205", "SRR8616206", "SRR8616217"]
K_LIST = [31]
P_LIST = [0, 40, 70, 100]

# ===== Workflow ===== #
rule all:
    input: 
        TRIM_DIR + "multiqc_report.html",
        TXIMPORT_RES + "gene-abundance-tximport.tsv",
        TXIMPORT_RES + "gene-counts-tximport.tsv",
        expand(REINDEER_QRY + "k{k}/query_results/out_query_Reindeer_P{p}_krator-1000-random-genes_contigs-{k}_0.out", k = K_LIST, p = P_LIST)

rule fetch_sra:
    output:
        RAW_FQ_DIR + "{smp}_1.fastq.gz",
        RAW_FQ_DIR + "{smp}_2.fastq.gz"
    params:
        smp = "{smp}"
    log:
        RAW_FQ_DIR + "log-fetch_sra.{smp}.txt"
    threads: 1
    shell:
        """
        {SRATOOLS_BIN}prefetch --output-directory {RAW_FQ_DIR} {params.smp} &> {log}
        {SRATOOLS_BIN}fastq-dump --outdir {RAW_FQ_DIR} --gzip --split-3 --defline-qual '+' {FQ_DIR}{params.smp}.sra &>> {log}
        """

rule cutadapt:
    input:
        raw1 = RAW_FQ_DIR + "{smp}_1.fastq.gz",
        raw2 = RAW_FQ_DIR + "{smp}_2.fastq.gz",
    output:
        trim1 = TRIM_FQ_DIR + "{smp}_1.trim.fastq.gz",
        trim2 = TRIM_FQ_DIR + "{smp}_2.trim.fastq.gz"
    log:
        TRIM_FQ_DIR + "log-cutadapt.{smp}.txt"
    threads: 3
    shell:
        """
        {CUTADAPT} -j 3 -q 10,10 -m 31 -o {output.trim1} -p {output.trim2} \
                   -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
                   -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
                   {input.raw1} {input.raw2} &> {log}
        """

rule fastqc:
    input:
        trim1 = TRIM_FQ_DIR + "{smp}_1.trim.fastq.gz",
        trim2 = TRIM_FQ_DIR + "{smp}_2.trim.fastq.gz"
    output:
        qc1 = TRIM_FQ_DIR + "{smp}_1.trim_fastqc.html",
        qc2 = TRIM_FQ_DIR + "{smp}_2.trim_fastqc.html"
    log:
        TRIM_FQ_DIR + "log-fastqc.{smp}.txt"
    threads: 3
    shell:
        "{FASTQC} --noextract --outdir {TRIM_FQ_DIR} --threads 3 {input.trim1} {input.trim2} &> {log}"

rule multiqc:
    input:
        expand(TRIM_FQ_DIR + "{smp}_{num}.trim_fastqc.html", smp = SMP_LIST, num = [1, 2])
    output:
        TRIM_FQ_DIR + "multiqc_report.html"
    log:
        TRIM_FQ_DIR + "log-multiqc.txt"
    threads: 1
    shell:
        "{MULTIQC} --outdir {TRIM_FQ_DIR} {TRIM_FQ_DIR} &> {log}"

rule tximport:
    input:
        gtf = GTF_PATH,
        kallisto_dir = directory(KALLISTO_RES)
    output:
        TXIMPORT_RES + "gene-abundance-tximport.tsv",
        TXIMPORT_RES + "gene-counts-tximport.tsv"
    log:
        TXIMPORT_RES + "log-tximport.txt"
    shell:
        """
        {RSCRIPT} a_tximport.R {input.gtf} {input.kallisto_dir} {TXIMPORT_RES} {SMP_LIST}
        """

rule reindeer_query:
    input:
        reindeer_idx = DATA_DIR + "new-rmAdapt/reindeer-index-k{k}",
        kmerator_ctgs = KMERATOR_DIR + "krator-1000-random-genes_contigs-{k}.fa"
    output:
        reindeer_qry = REINDEER_QRY + "k{k}/query_results/out_query_Reindeer_P{p}_krator-1000-random-genes_contigs-{k}_0.out"
    log:
        REINDEER_QRY + "log-reindeer-qry-k{k}-P{p}.txt"
    params:
        p = "{p}",
        qry_dir = REINDEER_QRY + "k{k}"
    shell:
        """
        singularity exec --bind /store:/store {REINDEER_IMG} Reindeer --query -P {params} -q {input.kmerator_ctgs} -l {input.reindeer_idx} -o {params.qry_dir}
        """

rule countTags:
    input:
        kmerator_ctgs = KMERATOR_DIR + "krator-1000-random-genes_contigs-{k}.fa"
    output:
        COUNTTAGS_RES + "countTags-k{k}.trim.tsv"
    log:
        COUNTTAGS_RES + "log-countTags-k{k}.txt"
    params:
        "{k}"
    shell:
        """
        {COUNTTAGS} -i {input.kmerator_ctgs} -k {params} --alltags --nostranded {TRIM_FQ_DIR}*.trim.fastq.gz > {output} 2> {log}
        """