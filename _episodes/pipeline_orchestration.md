---
layout: page
title: Pipeline Orchestration — Slurm dependencies vs Snakemake vs Nextflow
published: true
---

> ## 🗺️ Where this fits
> The HPC tutorials in this course chain analysis steps using `sbatch --dependency=afterok:<JID>` — a "background" pipeline that runs entirely on Slurm. That's one of three common ways to orchestrate a bioinformatics workflow on a cluster. This page lays the three side-by-side using the same problem we already solved in [HPC RNA-Seq](../HPC_RNA_SEQ/index.html), so you can pick the right tool when you write your own pipeline.
{: .callout}

> ## Learning Objectives
> By the end of this page you will be able to:
> - Explain what a *pipeline orchestrator* does and why you need one
> - Read a Slurm dependency chain, a Snakefile, and a Nextflow `.nf` script for the same task
> - Compare the three approaches on real criteria: error recovery, partial reruns, portability, parallelism, learning curve
> - Decide which orchestrator fits a given project
{: .objectives}

---

## What "orchestration" means

You have N samples. For each sample you run **download → trim → align → count**. Each step:

1. Depends on the previous step finishing
2. Should run in parallel across samples
3. Should not re-run if its inputs haven't changed
4. Should fail loudly when a sample is broken, not silently produce an empty BAM

An *orchestrator* tracks dependencies, parallelism, and re-execution rules so you don't have to. The three popular ones in genomics are:

| Approach | What it is | Where it runs |
|----------|-----------|---------------|
| **Slurm `sbatch --dependency`** | Native Slurm primitive — submit a job with a "wait until job X finishes" flag | Cluster only |
| **Snakemake** | Python DSL inspired by GNU `make`. Rules say "I produce file Y from file X by running command Z". | Any machine; submits to Slurm/SGE/PBS/cloud as a backend |
| **Nextflow** | Groovy/Java DSL built on dataflow channels. Processes consume channels and emit channels. | Any machine; submits to Slurm/SGE/AWS Batch/Kubernetes/etc. as a backend |

The same N-sample pipeline below, three ways:

---

## Approach 1: Slurm `--dependency` chain (what this course teaches)

Each step is a separate Slurm script. You submit them in order with `--dependency=afterok:<previous JID>`. Slurm holds each later job in `PD (Dependency)` state until its parents complete with exit code 0.

```bash
#!/bin/bash
# run_all.sh — orchestration via Slurm only
set -euo pipefail
cd ~/scratch/rnaseq/ATH
micromamba activate RNASEQ_bch709

DUMP_JID=$(sbatch --parsable fastq-dump.sh)
IDX_JID=$(cd reference && sbatch --parsable index.sh)
TRIM_JID=$(sbatch --parsable --dependency=afterok:${DUMP_JID} trim.sh)
ALIGN_JID=$(sbatch --parsable --dependency=afterok:${TRIM_JID}:${IDX_JID} align.sh)
MQC_JID=$(sbatch --parsable --dependency=afterany:${ALIGN_JID} multiqc.sh)

echo "Submitted: dump=$DUMP_JID idx=$IDX_JID trim=$TRIM_JID align=$ALIGN_JID multiqc=$MQC_JID"
```

Each `*.sh` is a normal Slurm batch script with `#SBATCH` directives and the actual tool command (fastp, STAR, featureCounts, ...).

**Strengths**
- **Zero new tools** — if you can `sbatch`, you're done. Works on every cluster.
- **Transparent** — `squeue` / `sacct` show exactly what's happening; failures debug like any Slurm job.
- **No DSL learning curve** — it's just bash + Slurm flags students already know.

**Weaknesses**
- **No partial reruns** — if step 3 fails, you re-submit step 3 *manually*; the orchestrator doesn't track which sample succeeded.
- **No file-level dependency tracking** — Slurm only knows job IDs, not "this BAM is up to date relative to its FASTQ input." Re-running means re-doing every sample even if some finished.
- **Cluster-only** — the same script does *not* work on AWS, Google Cloud, your laptop, or a different scheduler. Moving institutions = rewriting the orchestration.
- **Manual fan-out** — N samples = N rows in `samples.tsv` + array jobs. The fan-out logic lives in the script, not the orchestrator.

---

## Approach 2: Snakemake

Snakemake is a Python-based workflow language. You write *rules* — "I produce target file Y from input file X by running this shell command." Snakemake walks the rule graph backwards from your final target and figures out everything that needs to happen.

### Install Snakemake on Pronghorn

```bash
# Create a dedicated env for the orchestrator (separate from your bioinformatics env)
micromamba create -n smk -c conda-forge -c bioconda \
    'snakemake>=8' snakemake-executor-plugin-slurm -y
micromamba activate smk
snakemake --version    # → 8.x
```

The `snakemake-executor-plugin-slurm` plugin is what lets `snakemake --executor slurm` translate `resources:` blocks into `sbatch` flags.

```python
# Snakefile — orchestration via Snakemake
SAMPLES = ["SRR1761506", "SRR1761507", "SRR1761508",
           "SRR1761509", "SRR1761510", "SRR1761511"]

rule all:
    input:
        "qc/ATH_report.html",
        expand("counts/{sample}.tsv", sample=SAMPLES)

rule fastq_dump:
    output:
        r1 = "raw_data/{sample}_1.fastq.gz",
        r2 = "raw_data/{sample}_2.fastq.gz"
    shell:
        "wget -q -O {output.r1} https://ftp.sra.ebi.ac.uk/vol1/fastq/{wildcards.sample[:6]}/{wildcards.sample}/{wildcards.sample}_1.fastq.gz && "
        "wget -q -O {output.r2} https://ftp.sra.ebi.ac.uk/vol1/fastq/{wildcards.sample[:6]}/{wildcards.sample}/{wildcards.sample}_2.fastq.gz"

rule trim:
    input:
        r1 = "raw_data/{sample}_1.fastq.gz",
        r2 = "raw_data/{sample}_2.fastq.gz"
    output:
        r1   = "trim/{sample}_1.trimmed.fq.gz",
        r2   = "trim/{sample}_2.trimmed.fq.gz",
        json = "trim/{sample}_fastp.json"
    threads: 4
    shell:
        "fastp -i {input.r1} -I {input.r2} "
        "-o {output.r1} -O {output.r2} "
        "--json {output.json} --thread {threads}"

rule star_index:
    input:
        fasta = "reference/TAIR10.fasta",
        gtf   = "reference/TAIR10.gtf"
    output:
        directory("reference/star_index")
    threads: 8
    shell:
        "STAR --runMode genomeGenerate --runThreadN {threads} "
        "--genomeDir {output} --genomeFastaFiles {input.fasta} "
        "--sjdbGTFfile {input.gtf} --genomeSAindexNbases 12"

rule align:
    input:
        r1    = "trim/{sample}_1.trimmed.fq.gz",
        r2    = "trim/{sample}_2.trimmed.fq.gz",
        index = "reference/star_index"
    output:
        bam = "bam/{sample}.bam",
        log = "bam/{sample}.Log.final.out"
    threads: 8
    shell:
        "STAR --runMode alignReads --runThreadN {threads} "
        "--readFilesCommand zcat --genomeDir {input.index} "
        "--readFilesIn {input.r1} {input.r2} "
        "--outSAMtype BAM SortedByCoordinate "
        "--outFileNamePrefix bam/{wildcards.sample}."

rule featurecounts:
    input:
        bam = "bam/{sample}.bam",
        gtf = "reference/TAIR10.gtf"
    output:
        "counts/{sample}.tsv"
    threads: 4
    shell:
        "featureCounts -p --countReadPairs -T {threads} "
        "-a {input.gtf} -o {output} {input.bam}"

rule multiqc:
    input:
        expand("trim/{sample}_fastp.json", sample=SAMPLES),
        expand("bam/{sample}.Log.final.out", sample=SAMPLES),
        expand("counts/{sample}.tsv", sample=SAMPLES)
    output:
        "qc/ATH_report.html"
    shell:
        "multiqc . -o qc/ -n ATH_report --force"
```

### Adding resources and conda envs to each rule

The Snakefile above runs, but it doesn't tell Slurm how much memory or time each rule needs, and it relies on whatever tools happen to be on `PATH`. In practice you want each rule to declare its own resources and its own conda env so the workflow is reproducible across machines. Add `resources:` and `conda:` blocks:

```python
rule align:
    input:
        r1    = "trim/{sample}_1.trimmed.fq.gz",
        r2    = "trim/{sample}_2.trimmed.fq.gz",
        index = "reference/star_index"
    output:
        bam = "bam/{sample}.bam",
        log = "bam/{sample}.Log.final.out"
    threads: 8
    resources:
        mem_mb         = 32000,
        runtime        = 240,           # minutes — translates to --time
        slurm_partition = "cpu-core-0",
        slurm_account  = "cpu-s5-bch709-6"
    conda: "envs/star.yaml"             # a per-rule env Snakemake materializes
    shell:
        "STAR --runMode alignReads --runThreadN {threads} "
        "--readFilesCommand zcat --genomeDir {input.index} "
        "--readFilesIn {input.r1} {input.r2} "
        "--outSAMtype BAM SortedByCoordinate "
        "--outFileNamePrefix bam/{wildcards.sample}."
```

The matching `envs/star.yaml`:

```yaml
name: star
channels: [bioconda, conda-forge]
dependencies:
  - star=2.7.10b
  - samtools=1.20
```

`--use-conda` makes Snakemake build (or reuse) that env automatically before running the rule on the compute node. You can ship the `envs/` directory with the Snakefile and the workflow becomes portable.

### Run it on Pronghorn

```bash
# from the project root that holds the Snakefile
snakemake \
    --executor slurm \
    --default-resources slurm_account=cpu-s5-bch709-6 slurm_partition=cpu-core-0 mem_mb=4000 runtime=60 \
    --cores 32 \
    --jobs 20 \
    --use-conda \
    --rerun-incomplete \
    --keep-going
```

Flag-by-flag:

| Flag | Why |
|------|-----|
| `--executor slurm` | use the Slurm plugin instead of running locally |
| `--default-resources` | values applied to every rule that doesn't override them |
| `--cores 32` | total CPU budget (Snakemake won't oversubscribe) |
| `--jobs 20` | how many *Slurm jobs* it can have outstanding at once |
| `--use-conda` | build per-rule conda envs from the `conda:` directives |
| `--rerun-incomplete` | finish anything that crashed mid-run rather than skipping it |
| `--keep-going` | one bad sample doesn't sink the entire cohort |

### Visualize what's about to run

```bash
snakemake --dag             | dot -Tpng > dag.png
snakemake --rulegraph       | dot -Tpng > rulegraph.png
snakemake --filegraph       | dot -Tpng > filegraph.png
snakemake --summary
```

`--dag` shows every concrete file/sample pair as a node; `--rulegraph` collapses to one node per rule (much smaller graph for large cohorts); `--summary` prints what's up-to-date and what would re-run.

### Where the logs live

Snakemake captures every job's stdout/stderr under `.snakemake/log/` and the per-job Slurm `slurm-<JID>.out` ends up in the directory you submitted from. To find why a single sample failed:

```bash
# 1. Snakemake's wrapper log
cat .snakemake/log/$(ls -t .snakemake/log/ | head -1)
# 2. The actual Slurm job log Snakemake submitted on your behalf
cat .snakemake/slurm_logs/<rule>/<wildcards>/<JID>.{out,err}
```

### Demonstrating partial reruns

Touch one input and only the dependent rules rerun:

```bash
# pretend the trim of sample6 changed:
touch trim/SRR1761511_1.trimmed.fq.gz
snakemake --dry-run            # shows: only align(SRR1761511), featurecounts(SRR1761511), multiqc will rerun
snakemake --executor slurm     # runs exactly those three jobs, skipping the other 5 samples
```

That's the killer Snakemake feature: failed sample 6 can be re-aligned without re-doing the other 5.

**Strengths**
- **File-level dependency tracking** — change one input, only the affected rules re-run. Snakemake compares timestamps + hashes.
- **Partial reruns are free** — `snakemake counts/SRR1761510.tsv` re-runs only what's needed to produce that file.
- **DAG visualization** — `snakemake --dag | dot -Tpng > dag.png` draws the pipeline.
- **Backend swap** — same Snakefile runs on Slurm, SGE, your laptop, AWS, Google Cloud, by changing one flag.
- **conda / Apptainer integration** — each rule can declare its own `conda:` env that Snakemake materializes.

**Weaknesses**
- **Python DSL learning curve** — wildcards, expand, scattergather, checkpoints have a real onboarding cost.
- **One process per file** — declaring "this rule produces 14 chromosome BAMs from 1 input BAM" needs `checkpoint` or scatter-gather logic that surprises beginners.
- **Pickier about exact filenames** — typo in a wildcard and Snakemake just says "MissingInputException".

---

## Approach 3: Nextflow

Nextflow is built on *dataflow channels*. Each *process* consumes a channel and emits a channel. The orchestrator wires them up in a `workflow {}` block. Nextflow runs natively on Slurm, AWS Batch, Kubernetes, Google Cloud, with no script change.

### Install Nextflow on Pronghorn

Nextflow needs Java 17+. Install both into a dedicated env:

```bash
micromamba create -n nf -c conda-forge openjdk=17 -y
micromamba activate nf

# Get the launcher
curl -fsSL https://get.nextflow.io | bash
mv nextflow ~/.local/bin/        # or any dir on your PATH
nextflow -v                       # → nextflow version 24.x.y
```

```groovy
// main.nf — orchestration via Nextflow
nextflow.enable.dsl=2

params.samples  = ['SRR1761506','SRR1761507','SRR1761508',
                   'SRR1761509','SRR1761510','SRR1761511']
params.outdir   = 'results'
params.gtf      = 'reference/TAIR10.gtf'
params.fasta    = 'reference/TAIR10.fasta'

process FASTQ_DUMP {
    tag "$sample"
    publishDir "${params.outdir}/raw_data", mode: 'symlink'

    input:
        val sample
    output:
        tuple val(sample), path("${sample}_1.fastq.gz"), path("${sample}_2.fastq.gz")
    script:
    """
    P=\${sample:0:6}
    wget -q https://ftp.sra.ebi.ac.uk/vol1/fastq/\$P/${sample}/${sample}_1.fastq.gz
    wget -q https://ftp.sra.ebi.ac.uk/vol1/fastq/\$P/${sample}/${sample}_2.fastq.gz
    """
}

process TRIM {
    tag "$sample"
    cpus 4
    publishDir "${params.outdir}/trim", mode: 'symlink'

    input:
        tuple val(sample), path(r1), path(r2)
    output:
        tuple val(sample), path("${sample}_1.trimmed.fq.gz"), path("${sample}_2.trimmed.fq.gz"), emit: reads
        path "${sample}_fastp.json", emit: json
    script:
    """
    fastp -i $r1 -I $r2 \
          -o ${sample}_1.trimmed.fq.gz -O ${sample}_2.trimmed.fq.gz \
          --json ${sample}_fastp.json --thread $task.cpus
    """
}

process STAR_INDEX {
    cpus 8

    input:
        path fasta
        path gtf
    output:
        path 'star_index'
    script:
    """
    mkdir -p star_index
    STAR --runMode genomeGenerate --runThreadN $task.cpus \
         --genomeDir star_index --genomeFastaFiles $fasta \
         --sjdbGTFfile $gtf --genomeSAindexNbases 12
    """
}

process ALIGN {
    tag "$sample"
    cpus 8
    publishDir "${params.outdir}/bam", mode: 'symlink'

    input:
        tuple val(sample), path(r1), path(r2)
        path  index
    output:
        tuple val(sample), path("${sample}.Aligned.sortedByCoord.out.bam"), emit: bam
        path "${sample}.Log.final.out", emit: log
    script:
    """
    STAR --runMode alignReads --runThreadN $task.cpus --readFilesCommand zcat \
         --genomeDir $index --readFilesIn $r1 $r2 \
         --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ${sample}.
    """
}

process FEATURECOUNTS {
    tag "$sample"
    cpus 4
    publishDir "${params.outdir}/counts", mode: 'copy'

    input:
        tuple val(sample), path(bam)
        path gtf
    output:
        path "${sample}.tsv"
    script:
    """
    featureCounts -p --countReadPairs -T $task.cpus \
                  -a $gtf -o ${sample}.tsv $bam
    """
}

process MULTIQC {
    publishDir "${params.outdir}/qc", mode: 'copy'

    input:
        path '*'
    output:
        path 'ATH_report.html'
    script:
    """
    multiqc . -o . -n ATH_report --force
    """
}

workflow {
    samples_ch = Channel.fromList(params.samples)

    FASTQ_DUMP(samples_ch)
    TRIM(FASTQ_DUMP.out)
    STAR_INDEX(file(params.fasta), file(params.gtf))
    ALIGN(TRIM.out.reads, STAR_INDEX.out)
    FEATURECOUNTS(ALIGN.out.bam, file(params.gtf))

    MULTIQC(
        TRIM.out.json
            .mix(ALIGN.out.log)
            .mix(FEATURECOUNTS.out)
            .collect()
    )
}
```

`nextflow.config`:

```groovy
process {
    executor = 'slurm'
    queue    = 'cpu-core-0'
    clusterOptions = '--account=cpu-s5-bch709-6'
}
```

Run:

```bash
nextflow run main.nf -profile slurm -resume
```

### Per-process resources and conda envs

The `nextflow.config` snippet above is minimal. A real pipeline declares per-process resources so each step asks Slurm for the right amount:

```groovy
process {
    executor = 'slurm'
    queue    = 'cpu-core-0'
    clusterOptions = '--account=cpu-s5-bch709-6'

    // defaults — every process inherits these unless overridden below
    cpus   = 2
    memory = 4.GB
    time   = '1h'

    withName: STAR_INDEX  { cpus = 8; memory = 32.GB; time = '4h';   conda = "${baseDir}/envs/star.yaml"  }
    withName: ALIGN       { cpus = 8; memory = 32.GB; time = '6h';   conda = "${baseDir}/envs/star.yaml"  }
    withName: TRIM        { cpus = 4; memory = 8.GB;  time = '2h';   conda = "${baseDir}/envs/fastp.yaml" }
    withName: FEATURECOUNTS { cpus = 4; memory = 8.GB; time = '2h';  conda = "${baseDir}/envs/subread.yaml" }
    withName: MULTIQC     { cpus = 2; memory = 4.GB;  time = '30m';  conda = "${baseDir}/envs/multiqc.yaml" }
}

// or, instead of conda:
// process.container = 'quay.io/biocontainers/star:2.7.10b--h9ee0642_0'
// singularity { enabled = true; autoMounts = true }
```

### Run it on Pronghorn (with resume + report)

```bash
micromamba activate nf
nextflow run main.nf \
    -profile slurm \
    -with-conda \
    -with-report  reports/exec_report.html \
    -with-trace   reports/trace.txt \
    -with-dag     reports/dag.svg \
    -resume
```

| Flag | Why |
|------|-----|
| `-profile slurm` | use the Slurm config block from `nextflow.config` |
| `-with-conda` | materialize the per-process conda envs |
| `-with-report` | one rich HTML execution summary (resource usage, per-process duration, retries) |
| `-with-trace` | per-task TSV — what ran where, how long, exit code |
| `-with-dag` | flowchart of channels and processes |
| `-resume` | reuse the `work/` cache; only re-run what's stale |

### Demonstrating `-resume`

```bash
# First run — every process executes
nextflow run main.nf -profile slurm -with-conda

# Pretend one input changed:
touch raw_data/SRR1761511_1.fastq.gz

# Second run — only the affected processes re-run, the others are cached
nextflow run main.nf -profile slurm -with-conda -resume
```

The cache lives in `work/` (one directory per process invocation, named by content hash). Don't `rm -rf work/` if you want resume to keep working — point Nextflow at a `-w` work dir on scratch instead of cleaning it.

### Where the logs live

Each process writes its own `.command.{sh,log,out,err}` inside `work/<hash>/`. The top-level `.nextflow.log` summarizes everything. To debug a single process:

```bash
# find the failing task's hash
nextflow log <run_name> -f hash,name,status,workdir | grep FAILED
# inspect its scratch dir
cd <workdir>
cat .command.sh        # the command Nextflow actually launched
cat .command.log       # tool stderr+stdout
```

### Use an existing nf-core pipeline (no code to write)

For RNA-Seq, ChIP-Seq, variant calling, etc., you almost certainly don't need to write `main.nf` yourself — there's already a peer-reviewed pipeline at [nf-co.re](https://nf-co.re).

```bash
# Example: nf-core/rnaseq on Pronghorn
mkdir -p ~/scratch/nfcore_rnaseq && cd ~/scratch/nfcore_rnaseq

# A samplesheet.csv with one row per sample (sample, fastq_1, fastq_2, strandedness)
cat > samplesheet.csv <<'EOF'
sample,fastq_1,fastq_2,strandedness
WT_rep1,/path/to/SRR1761506_1.fastq.gz,/path/to/SRR1761506_2.fastq.gz,reverse
WT_rep2,/path/to/SRR1761507_1.fastq.gz,/path/to/SRR1761507_2.fastq.gz,reverse
WT_rep3,/path/to/SRR1761508_1.fastq.gz,/path/to/SRR1761508_2.fastq.gz,reverse
ABA_rep1,/path/to/SRR1761509_1.fastq.gz,/path/to/SRR1761509_2.fastq.gz,reverse
ABA_rep2,/path/to/SRR1761510_1.fastq.gz,/path/to/SRR1761510_2.fastq.gz,reverse
ABA_rep3,/path/to/SRR1761511_1.fastq.gz,/path/to/SRR1761511_2.fastq.gz,reverse
EOF

# Pronghorn-friendly profile (slurm + singularity for reproducibility)
nextflow run nf-core/rnaseq \
    -r 3.14.0 \
    -profile singularity \
    --input samplesheet.csv \
    --outdir results \
    --genome 'TAIR10' \
    --aligner star_salmon \
    -c pronghorn.config       # contains the slurm executor + account/partition
```

What you get: trimmed FASTQs, STAR + Salmon alignment + quantification, featureCounts, MultiQC, Picard metrics, all packaged with versioned containers — all in one command, fully resumable, ready to publish.

**Strengths**
- **Resume by default** — `-resume` re-uses cached process outputs unless inputs changed. Built into the executor, no extra flags per process.
- **Backend portability** — same `main.nf` runs on Slurm, AWS Batch, Azure, Kubernetes, your laptop. Change the config, not the code.
- **Curated nf-core community** — [nf-core](https://nf-co.re) has 100+ peer-reviewed pipelines (rnaseq, sarek, chipseq, ...) you can run as-is.
- **Strong containers / conda support** — `process.container = 'biocontainers/star:2.7.10'` and you're done.

**Weaknesses**
- **Groovy DSL** — fewer biologists know Groovy than Python. Channel semantics (`mix`, `collect`, `combine`) take time to internalize.
- **Heavier runtime** — JVM startup is non-trivial for tiny pipelines.
- **Two languages to debug** — channel manipulation in Groovy + script bodies in bash.

---

## Common pitfalls when migrating from Slurm scripts

The mental model shift is bigger than the syntax shift. Things that bite people:

| What you used to do (Slurm chain) | What you do now (Snakemake / Nextflow) |
|-----------------------------------|----------------------------------------|
| `mkdir -p bam/` at the top of every script | The orchestrator creates output directories from each rule's `output:` paths. Don't `mkdir` — let it. |
| `samples.tsv` + `--array=1-N` array job | `expand("...{sample}...", sample=SAMPLES)` (Snakemake) or `Channel.fromList(params.samples)` (Nextflow). The fan-out is *declarative*. |
| Hard-code thread count `-t 8` in the command | Use `{threads}` (Snakemake) or `$task.cpus` (Nextflow) so the resource request and the tool flag stay in sync. |
| Re-submit job 4 manually when sample 5 dies | `snakemake counts/SRR....tsv` (just rebuild that target) or `nextflow run main.nf -resume` (cache replays the rest). |
| Activate env once with `micromamba activate ...` | Per-rule `conda:` (Snakemake) or per-process `conda` / `container` (Nextflow). The orchestrator activates the right env per job. |
| Read job logs from `slurm-<JID>.out` | `.snakemake/log/` and `.snakemake/slurm_logs/<rule>/<wildcards>/` (Snakemake), `work/<hash>/.command.{log,out,err}` (Nextflow). |
| One big `run_all.sh` for the whole cohort | One Snakefile / `main.nf` for the pipeline; cohorts are just data. The same workflow runs over 6 samples or 600. |
| Ad-hoc cleanup with `rm -rf bam/*.tmp` | Use `temp(...)` (Snakemake) or `publishDir mode: 'symlink'` (Nextflow) so intermediates auto-clean. |

---

## Side-by-side comparison

| Criterion | Slurm `--dependency` | Snakemake | Nextflow |
|-----------|---------------------|-----------|----------|
| Language to learn | bash + Slurm flags | Python DSL (rule, wildcards, expand) | Groovy DSL (process, channel, workflow) |
| Tracks file-level changes | ❌ no | ✅ yes (timestamp + hash) | ✅ yes (cache by content + params) |
| Partial reruns ("just redo sample X") | ❌ manual | ✅ specify the target file | ✅ `-resume` reuses cache |
| Visualize the DAG | ❌ no | ✅ `--dag` | ✅ `-with-dag` |
| Backend portability (Slurm → AWS → laptop) | ❌ no — Slurm only | ✅ yes — `--executor` flag | ✅ yes — config profile |
| Container / conda integration | manual `module load` or env activate | ✅ per-rule `conda:`, `container:` | ✅ per-process `container`, `conda` |
| Community pipelines you can grab | none | [snakemake-workflows](https://github.com/snakemake-workflows) | [nf-core](https://nf-co.re) |
| Best when | Cluster-only one-off, students learning Slurm, very small pipelines | Python-heavy lab, fine-grained file tracking, custom rules | Cross-platform pipelines, nf-core ready, multi-cloud |
| Worst when | Resuming a partial run, moving to a different cluster | DSL syntax fights you on dynamic scatter-gather | Tiny one-off scripts (JVM overhead), purely-Python team |
| Initial scaffolding effort | Lowest (just bash) | Medium (Snakefile + envs) | Medium-high (main.nf + config + container choices) |

---

## Decision tree

```
Is this a quick, one-off pipeline you'll run once or twice?
├─ Yes → Slurm --dependency chain (this course's HPC tutorials)
└─ No, will run many times / re-run on partial failure / share with others
   ├─ Will it ever leave Slurm (laptop, AWS, another cluster)?
   │  ├─ No   → Snakemake — least overhead, file-tracking, Python-friendly
   │  └─ Yes  → Nextflow — backend-agnostic, nf-core community
   └─ Is there already an nf-core pipeline for it (rnaseq, sarek, chipseq, atacseq, ...)?
      └─ Yes → Use that nf-core pipeline directly (don't write your own).
```

---

## Same problem, three commit sizes

For the 6-sample Arabidopsis RNA-Seq pipeline you ran in [HPC RNA-Seq](../HPC_RNA_SEQ/index.html), the orchestration cost is roughly:

| Approach | Lines of orchestration code | Tool installs |
|----------|----------------------------|---------------|
| Slurm `--dependency` (`run_all.sh` + 4 `.sh` scripts) | ~120 lines | Just the bioinformatics tools |
| Snakemake (`Snakefile`) | ~70 lines | + `snakemake>=8` |
| Nextflow (`main.nf` + `nextflow.config`) | ~110 lines | + `nextflow>=23` (Java/JVM) |

You don't *win* on lines-of-code with Snakemake/Nextflow. You win on **what the orchestrator handles for you**: failure recovery, partial reruns, portability, DAG visualization, community pipelines.

---

## Practical advice for this class

For the assignments and HPC tutorials in this course, **stick with the Slurm `--dependency` approach** — it's the right tool when you're learning Slurm itself. But once you start running the same pipeline on 50+ samples, on multiple species, or on more than one cluster, **graduate to Snakemake or Nextflow**. The decision is rarely "which is better in the abstract" — it's "which one matches my lab's existing skills + where I need to run it."

> ## Try it
> Take the Arabidopsis 6-sample example from [HPC RNA-Seq](../HPC_RNA_SEQ/index.html) and re-implement it as either a Snakefile or a Nextflow `main.nf`. Time how long the implementation takes. Then run it twice — note how long the *second* run takes vs the first (this is the resume-from-cache benefit you don't get with the Slurm chain).
{: .challenge}

---

## Resources

- [Snakemake docs](https://snakemake.readthedocs.io/) — start with the tutorial; the [executor plugin docs](https://snakemake.readthedocs.io/en/stable/executing/cluster.html) cover Slurm.
- [Nextflow docs](https://www.nextflow.io/docs/latest/) — start with [Hello world](https://www.nextflow.io/docs/latest/getstarted.html), then [nf-core/rnaseq](https://nf-co.re/rnaseq) for a full Slurm-ready RNA-Seq pipeline.
- [Awesome-pipeline](https://github.com/pditommaso/awesome-pipeline) — a curated list of orchestrators if you want to see the larger landscape (Cromwell, Toil, WDL, etc.).
- [nf-core](https://nf-co.re) — 100+ peer-reviewed Nextflow pipelines you can run on Pronghorn today.
