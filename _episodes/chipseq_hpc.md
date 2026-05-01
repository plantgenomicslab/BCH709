---
layout: page
title: ChIP-Seq on Pronghorn HPC (Parallelized with Slurm)
published: true
---

> ## 🗺️ HPC series — you are here
> **1.** [HPC Cluster basics](../HPC_cluster/index.html) — SSH, file transfer, Micromamba, Slurm, dependencies
> **2.** [Resequencing pipeline on HPC](../resequencing_hpc/index.html) — BWA-MEM2 + GATK, variant calling
> **3. 🔵 ChIP-Seq pipeline (this page)** — minimap2 + MACS3, peak calling
> **4.** [RNA-Seq pipeline on HPC](../HPC_RNA_SEQ/index.html) — STAR + featureCounts, DE analysis
{: .callout}

> ## ✅ Before you start — pre-class checklist
> Skim this once before running any code below:
>
> - [ ] You can `ssh <netid>@pronghorn.rc.unr.edu` and see the Pronghorn prompt
> - [ ] `sacctmgr show user $USER withassoc` shows account `cpu-s5-bch709-6` / partition `cpu-core-0`
> - [ ] `~/scratch` symlink exists and points under `/data/gpfs/assoc/bch709-6/<netid>`
> - [ ] `micromamba --version` works on the login node (shell hook in `~/.bashrc`)
> - [ ] You ran the laptop version of [ChIP-Seq Tutorial](../chipseq_tutorial/index.html) at least once so the biology makes sense
>
> Any box unchecked? → go back to the [HPC Cluster lesson](../HPC_cluster/index.html) first.
{: .prereq}

---

## Overview — Why Parallelize ChIP-Seq on HPC?

ChIP-Seq experiments usually have multiple conditions × multiple replicates × 1 input per condition. A 2-condition × 3-replicate study already means **8 samples** (6 ChIP + 2 Input), and each needs its own alignment, dedup, peak call, and coverage track.

Running this serially on a laptop takes 12–24 hours. With Slurm array jobs on Pronghorn, all 8 samples are processed in parallel and the whole pipeline finishes in roughly the time of a single sample — about 45 minutes for small references, 1–2 hours for the full human genome.

### Three HPC parallelization strategies used here

| Strategy | What it parallelizes | Slurm feature |
|----------|---------------------|---------------|
| **Array jobs** | Same command, different samples | `#SBATCH --array=0-N` |
| **Job dependencies** | Pipeline stages (align → dedup → peak call) | `sbatch --dependency=afterok:JOBID` |
| **Paired treatment/control** | MACS3 takes (ChIP_rep, matching Input) from one sample sheet | `awk` / `join` at submission time |

### Pipeline DAG

```
               [download array]                 ← parallel across all FASTQs
                       │
                [trim+align array]              ← parallel across all samples
                       │
               [filter+dedup array]             ← parallel across all samples
                       │
               [bamCoverage array]              ← per-sample BigWig
                       │
              [MACS3 peak call array]           ← ChIP vs Input pairs in parallel
                       │
          [QC + motif + IDR single job]         ← aggregate
```

---

## 0. Setup on Pronghorn

### Connect and check your Slurm account

```bash
ssh <YOUR_NET_ID>@pronghorn.rc.unr.edu

sacctmgr show user $USER format=user,account,partition,qos%-30 \
  withassoc where user=$USER
```

Use the values shown there (typically `cpu-s5-bch709-6` / `cpu-core-0` / `student`) in every `#SBATCH` line below.

### Create the ChIP-Seq environment

```bash
# Use python 3.10 — bioconda has prebuilt wheels of macs3 / deeptools /
# idr / pysam for py3.10 that are linked against the conda sysroot's
# glibc 2.17 and run on Pronghorn's CentOS 7 compute nodes.
# Building these from source under py3.11 hits gcc 15's C23 default
# (`__isoc23_strtol@GLIBC_2.38`) and the cluster cannot load them.
micromamba create -n chipseq_bch709 -c conda-forge -c bioconda python=3.10 -y
micromamba activate chipseq_bch709

# Everything from bioconda — single install, no compilers, no pip builds.
# `numpy<2.0` because deeptools / idr aren't NumPy-2 ready yet.
micromamba install -c conda-forge -c bioconda \
    fastqc 'fastp>=0.24' minimap2 \
    'samtools>=1.20' bedtools 'tabix>=1.11' \
    openjdk=17 'picard>=3' homer \
    'pysam>=0.22' macs3 deeptools idr \
    'numpy<2.0' -y

# Upgrade pip first — older pip can't find the prebuilt `tiktoken`
# manylinux wheel (a transitive multiqc dep), tries to build from Rust
# source, and fails on Pronghorn (no Rust compiler).
pip install --upgrade pip

# Pip-only tools: multiqc + tiktoken. `tiktoken<0.8` pin avoids the Rust build.
pip install --prefer-binary 'tiktoken<0.8' multiqc
```

**Patch `libcrypto` so `samtools` runs (do this now, not after it crashes):**

Bioconda's `samtools` is linked against `libcrypto.so.1.0.0`, but the current OpenSSL package ships `libcrypto.so.3` (or `.1.1`). Without a symlink, you get:
```
samtools: error while loading shared libraries: libcrypto.so.1.0.0: cannot open shared object file
```

Create the symlink once, right after activating:

```bash
# Must be run with the env ACTIVE — $CONDA_PREFIX points at the env folder
cd "$CONDA_PREFIX/lib"
if   [ -f libcrypto.so.1.1 ]; then ln -sf libcrypto.so.1.1 libcrypto.so.1.0.0
elif [ -f libcrypto.so.3   ]; then ln -sf libcrypto.so.3   libcrypto.so.1.0.0
fi
cd - > /dev/null

# Verify samtools works
samtools --version | head -1
# → should print:  samtools 1.xx   (no error about libcrypto)
```

If you see `samtools 1.xx` you're done. If `samtools --version` still errors, ask the instructor before moving on.

> ## 🔑 Activate once in your login shell — every `sbatch` inherits the environment
> **Do this ONCE in the login shell before running any `sbatch` command:**
>
> ```bash
> micromamba activate chipseq_bch709
> which fastp     # confirm: should print a path inside ~/micromamba/envs/chipseq_bch709/
> ```
>
> By default, `sbatch` submits jobs with `--export=ALL`, which means each Slurm job inherits your current shell's environment — including `PATH` pointing at the activated env. So you **don't** need to put `micromamba activate` inside every batch script.
>
> **Sanity check:** submit a tiny test job and confirm the tool is found on the compute node too:
>
> ```bash
> sbatch -A cpu-s5-bch709-6 -p cpu-core-0 --time=00:05:00 --wrap="which fastp && fastp --version"
> # check the slurm-<jobid>.out log — should show the same fastp path + version
> ```
>
> If you open a new SSH session, the activation is lost — **just run `micromamba activate chipseq_bch709` again** before submitting.
{: .callout}

### Create the project directory on scratch

```bash
mkdir -p ~/scratch/chipseq/{raw,trim,bam,bw,peaks,qc,logs,scripts}
cd ~/scratch/chipseq
```

| Folder | Contents |
|--------|----------|
| `raw/` | Downloaded FASTQ files |
| `trim/` | Adapter-trimmed FASTQ |
| `bam/` | Aligned / filtered / deduplicated BAMs |
| `bw/` | BigWig signal tracks |
| `peaks/` | MACS3 narrowPeak / broadPeak results |
| `qc/` | fingerprint, correlation, motif outputs |
| `logs/` | Slurm log files |
| `scripts/` | All batch scripts below |

### Define your sample sheet

`samples.tsv` lists every sample (ChIP and Input) with its ENCODE accession. A separate column names the matching input so MACS3 can look it up at submission time.

Each row is **TAB-separated** (NOT spaces). Downstream scripts use `cut -f1` and `awk -F'\t' '$3=="chip"'` — they break if you accidentally use spaces. Use `printf` so the `\t` is unambiguously a tab:

```bash
printf 'sample\taccession\ttype\tcontrol\n'                 >  ~/scratch/chipseq/samples.tsv
printf 'ctcf_rep1\tENCFF001HTO\tchip\tinput_k562\n'         >> ~/scratch/chipseq/samples.tsv
printf 'ctcf_rep2\tENCFF001HTP\tchip\tinput_k562\n'         >> ~/scratch/chipseq/samples.tsv
printf 'input_k562\tENCFF001HTT\tinput\t-\n'                >> ~/scratch/chipseq/samples.tsv
```

**Verify it's tab-separated** (every row should have exactly 4 fields):

```bash
awk -F'\t' '{print NF}' ~/scratch/chipseq/samples.tsv | sort -u
# Expected output:  4
# If you see 1, you pasted spaces instead of tabs — re-run the printf lines above.

column -t -s $'\t' ~/scratch/chipseq/samples.tsv          # pretty-print to verify content
```

| Column | Description |
|--------|-------------|
| `sample` | Short name used in output filenames |
| `accession` | ENCODE file accession for `wget` download |
| `type` | `chip` or `input` — MACS3 runs only on `chip` rows |
| `control` | Which `sample` row to use as the `-c` input for MACS3; `-` for input rows |

> **Scaling up:** add more rows for additional replicates or conditions. All the `--array` ranges below are derived from `wc -l samples.tsv` so everything scales automatically.
{: .callout}

---

## 1. Download FASTQ in parallel (Array Job)

### `scripts/01_download.sh`

```bash
#!/bin/bash
#SBATCH --job-name=chip_download
#SBATCH --account=cpu-s5-bch709-6
#SBATCH --partition=cpu-core-0
#SBATCH --cpus-per-task=2
#SBATCH --mem=4g
#SBATCH --time=12:00:00
#SBATCH --array=1-3                         # = number of samples (data rows in samples.tsv)
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=<YOUR_EMAIL>
#SBATCH -o logs/01_download_%A_%a.out

set -euo pipefail
cd ~/scratch/chipseq

# Row 1 is header → add 1 to task ID
ROW=$((SLURM_ARRAY_TASK_ID + 1))
LINE=$(sed -n "${ROW}p" samples.tsv)
SAMPLE=$(echo "$LINE" | cut -f1)
ACC=$(echo "$LINE" | cut -f2)

URL="https://www.encodeproject.org/files/${ACC}/@@download/${ACC}.fastq.gz"
echo "[task ${SLURM_ARRAY_TASK_ID}] Downloading ${SAMPLE} (${ACC})"

# --retry-all-errors + -C -: harden against SSL eof drops on large
# transfers; resume partials in place. See HPC_RNA_SEQ.md fastq-dump
# (#64) for the underlying issue.
curl -fsSL --retry 5 --retry-all-errors --retry-delay 10 --max-time 3600 \
    -C - -o raw/${SAMPLE}.fastq.gz "${URL}"

# Validate the gzip immediately — prevents silent truncation
gunzip -t raw/${SAMPLE}.fastq.gz
ls -lh raw/${SAMPLE}.fastq.gz
echo "${SAMPLE} OK"
```

**Submit (capture the job ID so downstream steps can depend on it):**

```bash
# Make sure the env is active in THIS shell (once per login session)
micromamba activate chipseq_bch709
cd ~/scratch/chipseq

DL_JID=$(sbatch --parsable scripts/01_download.sh)
echo "Download job: ${DL_JID}"
squeue -u $USER
```

---

## 2. Reference Preparation (Single Job)

```bash
#!/bin/bash
#SBATCH --job-name=chip_ref
#SBATCH --account=cpu-s5-bch709-6
#SBATCH --partition=cpu-core-0
#SBATCH --cpus-per-task=4
#SBATCH --mem=16g
#SBATCH --time=02:00:00
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=<YOUR_EMAIL>
#SBATCH -o logs/02_reference_%j.out

set -euo pipefail
# Environment is activated in the login shell before `sbatch` is called —
# the PATH (with fastp, minimap2, samtools, …) is inherited automatically.

cd ~/scratch/chipseq

# hg19 genome (example data is ENCODE-legacy aligned to hg19)
# Use HTTPS + curl with retries — UCSC redirects HTTP to HTTPS and a
# plain wget can fail on the redirect. Two mirror choices:
#   primary: hgdownload.soe.ucsc.edu (US west)
#   mirror : hgdownload-euro.soe.ucsc.edu (Europe — faster from EU)
UCSC_URL="https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz"
UCSC_MIRROR="https://hgdownload-euro.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz"
# --retry-all-errors + -C -: harden against SSL eof drops on large
# transfers; resume partials in place. See HPC_RNA_SEQ.md fastq-dump
# (#64) for the underlying issue.
curl -fsSL --retry 5 --retry-all-errors --retry-delay 10 --max-time 3600 \
    -C - -o hg19.fa.gz "${UCSC_URL}" \
    || curl -fsSL --retry 5 --retry-all-errors --max-time 3600 \
        -C - -o hg19.fa.gz "${UCSC_MIRROR}"
gunzip -f hg19.fa.gz
mv -f hg19.fa reference.fasta

samtools faidx reference.fasta
minimap2 -d reference.mmi reference.fasta

# TSS BED for later heatmaps
REFG_URL="https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/refGene.txt.gz"
REFG_MIRROR="https://hgdownload-euro.soe.ucsc.edu/goldenPath/hg19/database/refGene.txt.gz"
curl -fsSL --retry 3 --retry-delay 10 --max-time 600 -o refGene.txt.gz "${REFG_URL}" \
    || curl -fsSL --retry 3 --max-time 600 -o refGene.txt.gz "${REFG_MIRROR}"
gunzip -f refGene.txt.gz
awk 'BEGIN{OFS="\t"} {
    if ($4=="+") print $3, $5, $5+1, $2, "0", $4;
    else         print $3, $6-1, $6, $2, "0", $4
}' refGene.txt | sort -k1,1 -k2,2n | uniq > TSS.bed

echo "Reference prep done."
```

Save as `scripts/02_reference.sh` and submit (independent of download — runs in parallel with Step 1):

```bash
REF_JID=$(sbatch --parsable scripts/02_reference.sh)
echo "Reference job: ${REF_JID}"
```

---

## 3. Trim + Align (Array Job)

`scripts/03_align.sh` — trims adapters with fastp, aligns with minimap2 `-ax sr`, sorts+indexes in one pipe.

```bash
#!/bin/bash
#SBATCH --job-name=chip_align
#SBATCH --account=cpu-s5-bch709-6
#SBATCH --partition=cpu-core-0
#SBATCH --cpus-per-task=8
#SBATCH --mem=32g
#SBATCH --time=12:00:00
#SBATCH --array=1-3
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=<YOUR_EMAIL>
#SBATCH -o logs/03_align_%A_%a.out

set -euo pipefail
set -o pipefail
# Environment is activated in the login shell before `sbatch` is called —
# the PATH (with fastp, minimap2, samtools, …) is inherited automatically.

cd ~/scratch/chipseq

ROW=$((SLURM_ARRAY_TASK_ID + 1))
SAMPLE=$(sed -n "${ROW}p" samples.tsv | cut -f1)

# ---- fastp ----
fastp \
    --in1  raw/${SAMPLE}.fastq.gz \
    --out1 trim/${SAMPLE}.trimmed.fq.gz \
    --qualified_quality_phred 20 \
    --length_required 25 \
    --thread ${SLURM_CPUS_PER_TASK} \
    --html trim/${SAMPLE}_fastp.html \
    --json trim/${SAMPLE}_fastp.json

# ---- minimap2 align + sort ----
minimap2 -ax sr -t ${SLURM_CPUS_PER_TASK} reference.mmi trim/${SAMPLE}.trimmed.fq.gz \
    2> logs/${SAMPLE}.minimap2.log \
  | samtools sort -@ ${SLURM_CPUS_PER_TASK} -T /tmp/${SAMPLE}_sort \
      -o bam/${SAMPLE}.bam

samtools index -@ ${SLURM_CPUS_PER_TASK} bam/${SAMPLE}.bam
samtools quickcheck bam/${SAMPLE}.bam && echo "${SAMPLE}.bam OK"
samtools flagstat bam/${SAMPLE}.bam > bam/${SAMPLE}.flagstat
```

**Submit with dependencies** — wait for both download (raw FASTQ) and reference prep (index) to succeed:

```bash
ALIGN_JID=$(sbatch --parsable --dependency=afterok:${DL_JID}:${REF_JID} scripts/03_align.sh)
echo "Align: ${ALIGN_JID}"
```

---

## 4. Filter MAPQ + Remove Duplicates (Array Job)

```bash
#!/bin/bash
#SBATCH --job-name=chip_dedup
#SBATCH --account=cpu-s5-bch709-6
#SBATCH --partition=cpu-core-0
#SBATCH --cpus-per-task=4
#SBATCH --mem=16g
#SBATCH --time=06:00:00
#SBATCH --array=1-3
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=<YOUR_EMAIL>
#SBATCH -o logs/04_dedup_%A_%a.out

set -euo pipefail
# Environment is activated in the login shell before `sbatch` is called —
# the PATH (with fastp, minimap2, samtools, …) is inherited automatically.

cd ~/scratch/chipseq

ROW=$((SLURM_ARRAY_TASK_ID + 1))
SAMPLE=$(sed -n "${ROW}p" samples.tsv | cut -f1)

# ---- Filter by MAPQ>=30 and drop unmapped ----
samtools view -b -q 30 -F 4 -@ ${SLURM_CPUS_PER_TASK} \
    bam/${SAMPLE}.bam > bam/${SAMPLE}.filt.bam
samtools index -@ ${SLURM_CPUS_PER_TASK} bam/${SAMPLE}.filt.bam

# ---- Remove PCR duplicates ----
picard -Xmx12g MarkDuplicates \
    I=bam/${SAMPLE}.filt.bam \
    O=bam/${SAMPLE}.dedup.bam \
    M=bam/${SAMPLE}.dedup.metrics \
    REMOVE_DUPLICATES=true \
    VALIDATION_STRINGENCY=SILENT
samtools index -@ ${SLURM_CPUS_PER_TASK} bam/${SAMPLE}.dedup.bam

# ---- Cleanup ----
rm -f bam/${SAMPLE}.bam bam/${SAMPLE}.bam.bai
rm -f bam/${SAMPLE}.filt.bam bam/${SAMPLE}.filt.bam.bai
```

**Submit:**

```bash
DEDUP_JID=$(sbatch --parsable --dependency=afterok:${ALIGN_JID} scripts/04_dedup.sh)
```

---

## 5. BigWig Signal Tracks (Array Job)

Every sample — including inputs — gets its own RPKM-normalized BigWig.

```bash
#!/bin/bash
#SBATCH --job-name=chip_bw
#SBATCH --account=cpu-s5-bch709-6
#SBATCH --partition=cpu-core-0
#SBATCH --cpus-per-task=8
#SBATCH --mem=16g
#SBATCH --time=04:00:00
#SBATCH --array=1-3
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=<YOUR_EMAIL>
#SBATCH -o logs/05_bw_%A_%a.out

set -euo pipefail
# Environment is activated in the login shell before `sbatch` is called —
# the PATH (with fastp, minimap2, samtools, …) is inherited automatically.

cd ~/scratch/chipseq

ROW=$((SLURM_ARRAY_TASK_ID + 1))
SAMPLE=$(sed -n "${ROW}p" samples.tsv | cut -f1)

bamCoverage \
    --bam bam/${SAMPLE}.dedup.bam \
    --outFileName bw/${SAMPLE}.bw \
    --normalizeUsing RPKM \
    --binSize 10 \
    --numberOfProcessors ${SLURM_CPUS_PER_TASK}
```

Submit:

```bash
BW_JID=$(sbatch --parsable --dependency=afterok:${DEDUP_JID} scripts/05_bigwig.sh)
```

---

## 6. MACS3 Peak Calling (ChIP-only Array Job)

Only the `chip` rows in `samples.tsv` get peak-called. Each task looks up its `control` column to find the matching Input BAM.

```bash
#!/bin/bash
#SBATCH --job-name=chip_macs3
#SBATCH --account=cpu-s5-bch709-6
#SBATCH --partition=cpu-core-0
#SBATCH --cpus-per-task=4
#SBATCH --mem=16g
#SBATCH --time=04:00:00
#SBATCH --array=1-2                    # number of `chip`-type rows (ctcf_rep1, ctcf_rep2)
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=<YOUR_EMAIL>
#SBATCH -o logs/06_macs3_%A_%a.out

set -euo pipefail
# Environment is activated in the login shell before `sbatch` is called —
# the PATH (with fastp, minimap2, samtools, …) is inherited automatically.

cd ~/scratch/chipseq

# Pull only the `chip`-type rows (skip header and input rows)
CHIP_ROWS=($(awk -F'\t' 'NR>1 && $3=="chip" {print $1":"$4}' samples.tsv))
IDX=$((SLURM_ARRAY_TASK_ID - 1))
PAIR=${CHIP_ROWS[$IDX]}
SAMPLE=${PAIR%:*}
CONTROL=${PAIR#*:}

echo "Task ${SLURM_ARRAY_TASK_ID}: ChIP=${SAMPLE} Control=${CONTROL}"

macs3 callpeak \
    -t bam/${SAMPLE}.dedup.bam \
    -c bam/${CONTROL}.dedup.bam \
    --format BAM \
    --gsize hs \
    --name ${SAMPLE} \
    --outdir peaks/${SAMPLE} \
    --qvalue 0.05 \
    --nomodel --extsize 147    # skip model building; use for low-depth / small-reference data

echo "Peaks written to peaks/${SAMPLE}/"
wc -l peaks/${SAMPLE}/${SAMPLE}_peaks.narrowPeak
```

> **Dynamic array size tip:** if you don't want to manually count ChIP rows, set the array size at submission time:
> ```bash
> N_CHIP=$(awk -F'\t' 'NR>1 && $3=="chip"' samples.tsv | wc -l)
> sbatch --array=1-${N_CHIP} --dependency=afterok:${DEDUP_JID} scripts/06_macs3.sh
> ```
{: .callout}

Submit:

```bash
N_CHIP=$(awk -F'\t' 'NR>1 && $3=="chip"' samples.tsv | wc -l)
MACS_JID=$(sbatch --parsable --array=1-${N_CHIP} --dependency=afterok:${DEDUP_JID} scripts/06_macs3.sh)
```

---

## 7. QC Aggregation (Single Job — runs after everything else)

Combines fingerprint, correlation, IDR, and MultiQC into one post-processing job.

```bash
#!/bin/bash
#SBATCH --job-name=chip_qc
#SBATCH --account=cpu-s5-bch709-6
#SBATCH --partition=cpu-core-0
#SBATCH --cpus-per-task=8
#SBATCH --mem=16g
#SBATCH --time=04:00:00
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=<YOUR_EMAIL>
#SBATCH -o logs/07_qc_%j.out

set -euo pipefail
# Environment is activated in the login shell before `sbatch` is called —
# the PATH (with fastp, minimap2, samtools, …) is inherited automatically.

cd ~/scratch/chipseq

# ---- deepTools fingerprint (global IP quality) ----
BAMS=$(awk -F'\t' 'NR>1 {print "bam/"$1".dedup.bam"}' samples.tsv | tr '\n' ' ')
LABELS=$(awk -F'\t' 'NR>1 {print $1}' samples.tsv | tr '\n' ' ')

plotFingerprint \
    --bamfiles ${BAMS} \
    --labels ${LABELS} \
    --plotFile qc/fingerprint.png \
    --outRawCounts qc/fingerprint.tab \
    --numberOfProcessors ${SLURM_CPUS_PER_TASK}

# ---- Replicate correlation ----
multiBamSummary bins \
    --bamfiles ${BAMS} \
    --labels ${LABELS} \
    --outFileName qc/multibam.npz \
    --numberOfProcessors ${SLURM_CPUS_PER_TASK}

plotCorrelation \
    --corData qc/multibam.npz \
    --corMethod pearson \
    --skipZeros \
    --whatToPlot heatmap \
    --colorMap RdYlBu \
    --plotFile qc/correlation.png \
    --outFileCorMatrix qc/correlation.txt

# ---- TSS heatmap for every ChIP sample ----
CHIP_BWS=$(awk -F'\t' 'NR>1 && $3=="chip" {print "bw/"$1".bw"}' samples.tsv | tr '\n' ' ')
CHIP_LABELS=$(awk -F'\t' 'NR>1 && $3=="chip" {print $1}' samples.tsv | tr '\n' ' ')

computeMatrix reference-point \
    --scoreFileName ${CHIP_BWS} \
    --regionsFileName TSS.bed \
    --referencePoint center \
    --upstream 3000 --downstream 3000 \
    --binSize 50 \
    --numberOfProcessors ${SLURM_CPUS_PER_TASK} \
    --samplesLabel ${CHIP_LABELS} \
    -o qc/tss_matrix.gz

plotHeatmap \
    --matrixFile qc/tss_matrix.gz \
    --outFileName qc/tss_heatmap.png \
    --colorMap Blues

# ---- IDR (only if exactly 2 replicates per condition) ----
# Generalize this loop for your conditions.
if [ -f peaks/ctcf_rep1/ctcf_rep1_peaks.narrowPeak ] && \
   [ -f peaks/ctcf_rep2/ctcf_rep2_peaks.narrowPeak ]; then
    sort -k5,5rn peaks/ctcf_rep1/ctcf_rep1_peaks.narrowPeak > qc/rep1_sorted.np
    sort -k5,5rn peaks/ctcf_rep2/ctcf_rep2_peaks.narrowPeak > qc/rep2_sorted.np
    idr \
        --samples qc/rep1_sorted.np qc/rep2_sorted.np \
        --input-file-type narrowPeak \
        --rank p.value \
        --output-file qc/idr_peaks.txt \
        --plot \
        --log-output-file qc/idr.log
fi

# ---- MultiQC aggregates fastp + flagstat + MarkDup + MACS3 + deepTools fingerprint ----
# NOTE: Explicit --module whitelist matches what this pipeline actually emits
# (no rsem, no gatk). We *also* pass --exclude rsem --exclude gatk as
# belt-and-suspenders: MultiQC 1.34's rsem parser ValueErrors on stray
# '#'-prefixed files in sibling dirs (e.g. stale ~/scratch/rnaseq/ outputs
# from RNA-Seq lessons), and the gatk module hits a rich.panel
# AttributeError on BQSR scatters — neither is relevant here, but we don't
# want either to load even if the whitelist parser changes upstream.
# Re-evaluate once the upstream MultiQC bugs are fixed.
multiqc . -o qc/ -n chipseq_report --force \
    --module fastp \
    --module samtools \
    --module picard \
    --module deeptools \
    --module macs2 \
    --exclude rsem \
    --exclude gatk
```

> ⚠️ **Do NOT run a bare `multiqc .` from `~/scratch/` or any parent directory.**
> MultiQC 1.34 will scan everything below the cwd and try to parse files from sibling pipelines (e.g. `~/scratch/rnaseq/`), tripping known module bugs (`rsem.py` ValueError on stray `#` lines, the GATK BQSR `rich.panel` AttributeError). The safest path is just:
>
> ```bash
> sbatch scripts/07_qc.sh
> ```
>
> If you must run it interactively, mirror the script's flag set exactly — `cd` into the chipseq workspace and pass the same module whitelist + excludes:
>
> ```bash
> cd ~/scratch/chipseq
> multiqc . -o qc/ -n chipseq_report --force \
>     --module fastp --module samtools --module picard \
>     --module deeptools --module macs2 \
>     --exclude rsem --exclude gatk
> ```

Submit:

```bash
QC_JID=$(sbatch --parsable --dependency=afterok:${MACS_JID}:${BW_JID}:${REF_JID}:${ALIGN_JID} scripts/07_qc.sh)
```

---

## 8. Submit the Entire Pipeline with One Script

```bash
#!/bin/bash
# scripts/run_all.sh
set -euo pipefail
cd ~/scratch/chipseq

# Activate the env in THIS shell so every sbatch below inherits the PATH
# (sbatch --export=ALL is the default — the submitted jobs see the same tools)
export MAMBA_ROOT_PREFIX="${MAMBA_ROOT_PREFIX:-$HOME/micromamba}"
eval "$(micromamba shell hook --shell=bash)"
micromamba activate chipseq_bch709

N_SAMPLES=$(awk 'NR>1' samples.tsv | wc -l)
N_CHIP=$(awk -F'\t' 'NR>1 && $3=="chip"' samples.tsv | wc -l)

DL_JID=$(sbatch   --parsable --array=1-${N_SAMPLES}                                 scripts/01_download.sh)
REF_JID=$(sbatch  --parsable                                                         scripts/02_reference.sh)
ALIGN_JID=$(sbatch --parsable --array=1-${N_SAMPLES} --dependency=afterok:${DL_JID}:${REF_JID} scripts/03_align.sh)
DEDUP_JID=$(sbatch --parsable --array=1-${N_SAMPLES} --dependency=afterok:${ALIGN_JID}         scripts/04_dedup.sh)
BW_JID=$(sbatch    --parsable --array=1-${N_SAMPLES} --dependency=afterok:${DEDUP_JID}         scripts/05_bigwig.sh)
MACS_JID=$(sbatch  --parsable --array=1-${N_CHIP}    --dependency=afterok:${DEDUP_JID}         scripts/06_macs3.sh)
QC_JID=$(sbatch    --parsable --dependency=afterok:${MACS_JID}:${BW_JID}:${REF_JID}:${ALIGN_JID} scripts/07_qc.sh)

cat <<EOF
Submitted ChIP-Seq pipeline (${N_SAMPLES} samples, ${N_CHIP} ChIP):
  download        ${DL_JID}
  reference       ${REF_JID}
  trim+align      ${ALIGN_JID}
  filter+dedup    ${DEDUP_JID}
  bigwig          ${BW_JID}
  macs3           ${MACS_JID}
  qc              ${QC_JID}

Monitor:  squeue -u \$USER
Cancel all: scancel ${DL_JID} ${REF_JID} ${ALIGN_JID} ${DEDUP_JID} ${BW_JID} ${MACS_JID} ${QC_JID}
EOF
```

Run:

```bash
chmod +x scripts/*.sh
bash scripts/run_all.sh
```

### 🧑‍💻 Hands-on walkthrough — submit the pipeline step-by-step

If you want to see exactly what `run_all.sh` does (or debug one step), submit each stage manually. Every `sbatch` returns a **job ID** that the next step depends on.

**Do this first (login shell — one time):**

```bash
micromamba activate chipseq_bch709
cd ~/scratch/chipseq

N_SAMPLES=$(awk 'NR>1' samples.tsv | wc -l)
N_CHIP=$(awk -F'\t' 'NR>1 && $3=="chip"' samples.tsv | wc -l)
echo "Samples: $N_SAMPLES   |   ChIP: $N_CHIP"
```

**Then submit each step — each line is one command:**

```bash
# --- Step 1: download FASTQs (no prerequisites) ---
DL_JID=$(sbatch --parsable --array=1-${N_SAMPLES} scripts/01_download.sh)
echo "download     → $DL_JID"

# --- Step 2: reference prep (independent — runs in parallel with Step 1) ---
REF_JID=$(sbatch --parsable scripts/02_reference.sh)
echo "reference    → $REF_JID"

# --- Step 3: trim + align (waits for BOTH download and reference) ---
ALIGN_JID=$(sbatch --parsable --array=1-${N_SAMPLES} \
    --dependency=afterok:${DL_JID}:${REF_JID} scripts/03_align.sh)
echo "align        → $ALIGN_JID"

# --- Step 4: filter MAPQ + dedup (waits for align) ---
DEDUP_JID=$(sbatch --parsable --array=1-${N_SAMPLES} \
    --dependency=afterok:${ALIGN_JID} scripts/04_dedup.sh)
echo "dedup        → $DEDUP_JID"

# --- Step 5: bigwig (waits for dedup) ---
BW_JID=$(sbatch --parsable --array=1-${N_SAMPLES} \
    --dependency=afterok:${DEDUP_JID} scripts/05_bigwig.sh)
echo "bigwig       → $BW_JID"

# --- Step 6: MACS3 peak calling (ChIP only, waits for dedup) ---
MACS_JID=$(sbatch --parsable --array=1-${N_CHIP} \
    --dependency=afterok:${DEDUP_JID} scripts/06_macs3.sh)
echo "macs3        → $MACS_JID"

# --- Step 7: QC aggregation (waits for bigwig, MACS3, reference TSS.bed, and align outputs) ---
QC_JID=$(sbatch --parsable --dependency=afterok:${MACS_JID}:${BW_JID}:${REF_JID}:${ALIGN_JID} scripts/07_qc.sh)
echo "qc           → $QC_JID"

# Check that everything is queued
squeue -u $USER
# Steps 3-7 should show state PD with reason (Dependency)
```

After pasting the block, `squeue` shows all 7 jobs — some running, most pending. Slurm takes care of the ordering.

> ## Why copy the commands into your terminal, not a script?
> The hands-on walkthrough is literally what `run_all.sh` does — but by typing each line you *see* each job ID appear and can inspect things in between. Once you're comfortable, just run `bash scripts/run_all.sh` next time.
{: .callout}

---

## 9. Monitoring

```bash
squeue -u $USER                                               # live queue
sacct -u $USER --starttime=today --format=JobID,JobName,State,Elapsed,MaxRSS,ExitCode
tail -f logs/06_macs3_${MACS_JID}_1.out                       # follow one task's output
```

Cancel the whole pipeline if something's wrong:

```bash
scancel -u $USER
```

---

## 10. Resource Tuning

After the first successful run, review actual usage and shrink the requests:

```bash
sacct -u $USER --starttime=today \
      --format=JobID,JobName,State,Elapsed,MaxRSS,ReqMem,ReqCPUs
```

| If you see | Then | Why |
|------------|------|-----|
| `MaxRSS` ≪ `ReqMem` | Lower `--mem` | Job queues faster on smaller requests |
| `Elapsed` ≪ `ReqTime` | Lower `--time` | Same |
| `OUT_OF_MEMORY` | Raise `--mem` to ~2× `MaxRSS` | Leave headroom for Java GC |
| `TIMEOUT` on align | Raise `--time`, check fastq size | Very deep BAMs take longer |

> **Heap vs Slurm memory:** Picard and GATK obey `-Xmx` (Java heap), **not** `--mem`. Always set `-Xmx` 2–4 GB below `--mem` (e.g. `--mem=16g` → `-Xmx12g`) so the JVM has room for off-heap stack/GC.
{: .callout}

---

## 11. Expected Wall Time

For 3 samples (2 ChIP + 1 Input, ~35M reads each, full hg19):

| Step | Cores × array | Wall time |
|------|---------------|-----------|
| Download | 2 × 3 | ~10 min |
| Reference | 4 | ~10 min |
| Trim + align | 8 × 3 | ~25 min |
| Filter + dedup | 4 × 3 | ~10 min |
| BigWig | 8 × 3 | ~5 min |
| MACS3 | 4 × 2 | ~5 min |
| QC aggregation | 8 | ~5 min |
| **End-to-end** | | **~45 min** |

Versus **~10 hours** on a laptop for the same 3 samples.

---

## 12. HPC-Specific Errors

> **`Batch job submission failed: Invalid account or account/partition combination`**
>
> **Fix:** re-run `sacctmgr show user $USER ... withassoc` and copy values exactly into each script.
{: .callout}

> **MACS3 array task `sed: -e expression: unmatched ... samples.tsv`**
>
> **Cause:** `samples.tsv` has blank trailing lines or the task index doesn't match a real row.
> **Fix:**
> ```bash
> cat -A samples.tsv      # look for stray $ or tabs
> awk 'NF' samples.tsv > samples.tsv.clean && mv samples.tsv.clean samples.tsv
> ```
{: .callout}

> **MACS3 hangs or fails with `peaks: 0`**
>
> **Cause:** ChIP too shallow to build the fragment model; or input far deeper than ChIP.
> **Fix:** our script already passes `--nomodel --extsize 147`. If peaks are still zero, check FRiP and confirm the ChIP/Input pair is correct.
{: .callout}

> **`IDR AttributeError: module 'numpy' has no attribute 'int'`**
>
> **Cause:** IDR 2.0.3 incompatible with NumPy ≥ 1.24.
> **Fix:** `pip install "numpy<1.24"` (already in the setup above).
{: .callout}

> **Pipeline queued but downstream jobs stuck in `(Dependency)` state forever**
>
> **Cause:** A parent array task exited non-zero. `afterok` treats the whole array as failed if any single task fails.
> **Fix:** find the failing task with `sacct -u $USER --starttime=today | grep FAILED`, fix it, resubmit **only that step** and update downstream `--dependency` to the new job ID. See the next callout for the exact per-step commands when *every* upstream step is already `COMPLETED`.
{: .callout}

> **Re-running a single failed step (drop `--dependency=`)**
>
> When only one step failed and every upstream step is already in `COMPLETED` state, **submit just the failed script with no `--dependency=` flag**. The chained driver in Section 8 only needs `--dependency=...` because it submits everything at once with `sbatch --parsable` capturing job IDs that don't yet exist as completed jobs.
>
> Use these commands. Each is independent — pick the one for the step that failed, drop the rest:
>
> ```bash
> sbatch scripts/01_download.sh   # array — FASTQ download (one task per row in samples.tsv)
> sbatch scripts/02_reference.sh  # download + minimap2 index + TSS.bed
> sbatch scripts/03_align.sh      # array — fastp + minimap2 -ax sr
> sbatch scripts/04_dedup.sh      # array — MAPQ filter + samtools markdup
> sbatch scripts/05_bigwig.sh     # array — deepTools bamCoverage → bw/
> sbatch scripts/06_macs3.sh      # array (ChIP rows only) — MACS3 callpeak
> sbatch scripts/07_qc.sh         # aggregated QC — fingerprint + IDR + MultiQC
> ```
>
> **Re-run ONLY a subset of array tasks** (e.g. only the 2nd ChIP replicate failed in step 6):
>
> ```bash
> sbatch --array=2 scripts/06_macs3.sh
> ```
>
> overrides the `--array=` line in the script header and runs only the listed task IDs.
>
> **Special cases that need extra setup before re-running:**
>
> - **Step 6 (`06_macs3.sh`)** — needs both the ChIP and the matching Input dedup BAM (`bam/${SAMPLE}.dedup.bam` and `bam/${CONTROL}.dedup.bam`). If only one ChIP replicate failed, confirm the input BAM (e.g. `bam/input_k562.dedup.bam`) still exists before resubmitting.
> - **Step 7 (`07_qc.sh`)** — IDR step is gated on `peaks/ctcf_rep1/ctcf_rep1_peaks.narrowPeak` and `peaks/ctcf_rep2/...narrowPeak`. If MACS3 was rerun, those files exist — IDR will run. If you want the QC report even when MACS3 partly failed, resubmit with `--dependency=afterany` instead of `afterok`.
>
> **Sanity check before resubmitting** — verify upstream outputs exist:
>
> ```bash
> sacct -u $USER --format=JobID,JobName%-25,State,ExitCode --starttime today
> ls -lh ~/scratch/chipseq/{bam,bw,peaks}
> ```
{: .callout}

---

## 13. Next Steps

- **Differential binding:** with multiple conditions, run [DiffBind](https://bioconductor.org/packages/DiffBind/) in R on the dedup BAMs + narrowPeak files.
- **Motif discovery:** after the pipeline, run `findMotifsGenome.pl` (HOMER) on the PASS peaks.
- **Peak annotation:** use `ChIPseeker` in R to annotate peaks relative to genes/TSS.
- **Scaling to 100+ samples:** add more rows to `samples.tsv` — the `N_SAMPLES` / `N_CHIP` variables in `run_all.sh` pick up the new count automatically.
- **Workflow managers:** for production use, look at [Snakemake](https://snakemake.readthedocs.io/) or [Nextflow](https://www.nextflow.io/) — they handle this dependency bookkeeping natively.

---

## References

| Resource | Link |
|---------|------|
| Slurm Array Jobs | [slurm.schedmd.com/job_array.html](https://slurm.schedmd.com/job_array.html) |
| Slurm Job Dependencies | [slurm.schedmd.com/sbatch.html#OPT_dependency](https://slurm.schedmd.com/sbatch.html#OPT_dependency) |
| MACS3 | [macs3-project.github.io/MACS](https://macs3-project.github.io/MACS/) |
| deepTools | [deeptools.readthedocs.io](https://deeptools.readthedocs.io) |
| Basic ChIP-Seq (laptop) | [ChIP-Seq Tutorial]({{site.baseurl}}/episodes/chipseq_tutorial/) |
| HPC Cluster fundamentals | [HPC Cluster Lesson]({{site.baseurl}}/episodes/HPC_cluster/) |
| Basic Resequencing HPC | [Resequencing HPC]({{site.baseurl}}/episodes/resequencing_hpc/) |
