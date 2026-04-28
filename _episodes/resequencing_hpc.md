---
layout: page
title: Resequencing on Pronghorn HPC (Parallelized with Slurm)
published: true
---

> ## 🗺️ HPC series — you are here
> **1.** [HPC Cluster basics](../HPC_cluster/index.html) — SSH, file transfer, Micromamba, Slurm, dependencies
> **2. 🔵 Resequencing pipeline (this page)** — BWA-MEM2 + GATK, variant calling
> **3.** [ChIP-Seq pipeline on HPC](../chipseq_hpc/index.html) — minimap2 + MACS3, peak calling
> **4.** [RNA-Seq pipeline on HPC](../HPC_RNA_SEQ/index.html) — STAR + featureCounts, DE analysis
{: .callout}

> ## ✅ Before you start — pre-class checklist
> Skim this once before running any code below:
>
> - [ ] You can `ssh <netid>@pronghorn.rc.unr.edu` and see the Pronghorn prompt
> - [ ] `sacctmgr show user $USER withassoc` shows account `cpu-s5-bch709-6` / partition `cpu-core-0`
> - [ ] `~/scratch` symlink exists and points under `/data/gpfs/assoc/bch709-6/<netid>`
> - [ ] `micromamba --version` works on the login node (shell hook in `~/.bashrc`)
> - [ ] You ran the laptop version of [Resequencing Tutorial](../resequencing_tutorial/index.html) at least once so the biology makes sense
>
> Any box unchecked? → go back to the [HPC Cluster lesson](../HPC_cluster/index.html) first.
{: .prereq}

---

## Overview — Why Parallelize on HPC?

On your laptop, the [basic resequencing tutorial]({{site.baseurl}}/episodes/resequencing_tutorial/) runs **serially**, one step at a time, for each sample:

```
sample1 → align → markdup → BQSR → HaplotypeCaller → ...
sample2 → align → markdup → BQSR → HaplotypeCaller → ...
```

For 2 samples this takes 4-6 hours. For 100 samples on one laptop, it would take **weeks**. HPC lets you run all samples simultaneously, so 100 samples take about the same wall-clock time as 1 sample.

### Three HPC parallelization strategies

| Strategy | What it parallelizes | Slurm feature |
|----------|---------------------|---------------|
| **Array jobs** | Same command, different samples | `#SBATCH --array=0-N` |
| **Scatter-gather** | One sample, split by chromosome | Multiple jobs + one merge job |
| **Job dependencies** | Pipeline stages (align → markdup → ...) | `sbatch --dependency=afterok:JOBID` |

We use all three in this tutorial.

### Pipeline DAG (Directed Acyclic Graph)

```
                [align array]                      ← parallel across N samples
                       │
                [markdup array]                    ← parallel across N samples
                       │
                 [BQSR array]                      ← parallel across N samples
                       │
           [HaplotypeCaller scatter]               ← parallel across N × 7 chromosomes
                       │
                [gather per sample]                ← merge chr GVCFs → per-sample GVCF
                       │
          [CombineGVCFs + GenotypeGVCFs]           ← single job (serial)
                       │
           [Filter + Annotate + PLINK]             ← single job
```

---

## 0. Setup on Pronghorn

### Connect and check your account

```bash
ssh <YOUR_NET_ID>@pronghorn.rc.unr.edu

# Confirm your Slurm account/partition (required for sbatch)
sacctmgr show user $USER format=user,account,partition,qos%-30 \
  withassoc where user=$USER
```

You should see: `cpu-s5-bch709-6` / `cpu-core-0` / `student`. Use those values in every `#SBATCH` directive below.

### Create the resequencing environment

```bash
micromamba create -n reseq_bch709 -c conda-forge -c bioconda python=3.11 -y
micromamba activate reseq_bch709

micromamba install -c conda-forge -c bioconda \
    fastqc 'fastp>=0.24' bwa-mem2 \
    'samtools>=1.20' 'bcftools>=1.20' 'tabix>=1.11' \
    openjdk=17 'picard>=3' gatk4 'snpeff<=5.2' plink -y
# NOTE: snpeff is pinned to <=5.2 because 5.3.x+ requires Java 21,
# while we ship openjdk=17 (which Picard and GATK4 still target).

# Note: we deliberately do NOT install sra-tools. Bioconda's sra-tools 3.x is
# built against GLIBC 2.27+, which is newer than Pronghorn's system libc — the
# binary fails with "GLIBC_2.27 not found" on the compute nodes. The download
# script in Step 1 pulls FASTQ directly from ENA over HTTPS, so sra-tools is
# not needed.

# Upgrade pip first — older pip can't find the prebuilt `tiktoken`
# manylinux wheel (a transitive multiqc dep), tries to build it from
# Rust source, fails on Pronghorn (no Rust compiler).
pip install --upgrade pip

# MultiQC + pinned deps. `tiktoken<0.8` is the safety pin — older
# tiktoken has stable cp311 linux wheels.
pip install --prefer-binary \
    'numpy<2.0' 'pyarrow<17' 'tiktoken<0.8' multiqc
```

**Patch `libcrypto` so `samtools` / `bcftools` run (do this now, not after they crash):**

Bioconda's `samtools` / `bcftools` / `tabix` are linked against `libcrypto.so.1.0.0`, but the current OpenSSL package in the env ships `libcrypto.so.3` (or `.1.1`). Without a symlink, you'll see:
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

# Verify the three most-used HTS tools all run
samtools --version  | head -1    # → samtools 1.xx
bcftools --version  | head -1    # → bcftools 1.xx
tabix --version     | head -1    # → tabix (htslib) 1.xx
```

Expected output:

```
samtools 1.23.1
bcftools 1.23.1
tabix (htslib) 1.23.1
```

If all three print a version number (no `libcrypto` error), your environment is ready. If any still error, ask the instructor before moving on.

> ## 🔑 Activate once in your login shell — every `sbatch` inherits the environment
> **Do this ONCE in the login shell before running any `sbatch` command:**
>
> ```bash
> micromamba activate reseq_bch709
> which bwa-mem2   # confirm: should print a path inside ~/micromamba/envs/reseq_bch709/
> which samtools
> which gatk
> ```
>
> By default, `sbatch` submits jobs with `--export=ALL`, which means each Slurm job inherits your current shell's environment — including `PATH` pointing at the activated env. So you **don't** need to put `micromamba activate` inside every batch script.
>
> **Sanity check:** submit a tiny test job and confirm the tool is found on the compute node too:
>
> ```bash
> sbatch -A cpu-s5-bch709-6 -p cpu-core-0 --time=00:05:00 --wrap="which bwa-mem2 && bwa-mem2 version"
> # check the slurm-<jobid>.out log — should show the same path + version
> ```
>
> Expected `slurm-<jobid>.out`:
>
> ```
> ~/micromamba/envs/reseq_bch709/bin/bwa-mem2
> Looking to launch executable ".../bwa-mem2.avx2", simd = .avx2
> Launching executable ".../bwa-mem2.avx2"
> 2.2.1
> ```
>
> If you open a new SSH session, the activation is lost — **just run `micromamba activate reseq_bch709` again** before submitting.
{: .callout}

### Create the project directory on scratch

**Never put big BAM/VCF files in `$HOME`** — it has tight quotas and slow I/O. Use `~/scratch`.

```bash
mkdir -p ~/scratch/reseq/{raw,trim,bam,vcf,logs,scripts}
cd ~/scratch/reseq
```

| Folder | Contents |
|--------|----------|
| `raw/` | Downloaded FASTQ files |
| `trim/` | Adapter-trimmed FASTQ |
| `bam/` | Aligned / markdup / recalibrated BAMs |
| `vcf/` | GVCFs and final VCFs |
| `logs/` | Slurm `.out` / `.err` log files |
| `scripts/` | All the `.sh` batch scripts below |

### Define your sample sheet

Each row is **TAB-separated** (NOT spaces): `sample_id <TAB> SRR_accession`. Downstream scripts use `cut -f1` and `awk -F'\t'`, both of which expect real tabs — if you accidentally use spaces, every script breaks.

Use `printf` so the `\t` is unambiguously a tab character (a heredoc with literal tabs often gets converted to spaces during copy-paste from a browser):

```bash
printf 'sample1\tSRR519585\nsample2\tSRR519586\n' > ~/scratch/reseq/samples.tsv
```

**Verify it's actually tab-separated** (this should print just `2` — meaning every row has exactly 2 tab-separated fields):

```bash
awk -F'\t' '{print NF}' ~/scratch/reseq/samples.tsv | sort -u
# Expected output:  2
# If you see 1, you have spaces instead of tabs — re-run the printf above.

cat ~/scratch/reseq/samples.tsv          # quick visual check
```

> **Adding more samples later** is trivial — just append rows with `printf 'sampleN\tSRR_NNN\n' >> samples.tsv` (use `>>` to append, not `>` which overwrites!). The array size in every script automatically scales (you only change one number).
{: .callout}

---

## 1. Download FASTQ in parallel (Array Job)

One Slurm array task per sample. Each task downloads its own R1/R2 independently and in parallel.

### `scripts/01_download.sh`

```bash
#!/bin/bash
#SBATCH --job-name=download
#SBATCH --account=cpu-s5-bch709-6
#SBATCH --partition=cpu-core-0
#SBATCH --cpus-per-task=2
#SBATCH --mem=4g
#SBATCH --time=12:00:00
#SBATCH --array=1-2                          # one task per sample (rows in samples.tsv)
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=<YOUR_EMAIL>
#SBATCH -o logs/01_download_%A_%a.out        # %A=array job id, %a=task index

set -euo pipefail

cd ~/scratch/reseq

# Pick the N-th line from samples.tsv using the array task ID
LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" samples.tsv)
SAMPLE=$(echo "$LINE" | cut -f1)
SRR=$(echo "$LINE" | cut -f2)

echo "[task ${SLURM_ARRAY_TASK_ID}] Downloading ${SAMPLE} (${SRR}) from ENA"

# Download paired-end FASTQ from ENA (mirrors NCBI SRA, gzipped FASTQ ready-to-use).
# We deliberately avoid `prefetch` / `fasterq-dump` — bioconda's sra-tools 3.x is
# built against GLIBC 2.27+, which is newer than Pronghorn's system libc, so the
# binary fails with "GLIBC_2.27 not found" on the compute nodes.
mkdir -p raw

# Ask ENA for the exact fastq URLs (handles paired-end / single-end / multi-file).
META_URL="https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${SRR}&result=read_run&fields=fastq_ftp&format=tsv"
URLS=$(curl -fsSL --retry 3 --max-time 60 "${META_URL}" \
        | tail -n +2 \
        | awk -F'\t' '{print $NF}' \
        | tr ';' '\n' \
        | sed '/^$/d')
[ -n "${URLS}" ] || { echo "ERROR: ENA returned no fastq URLs for ${SRR}"; exit 1; }

i=1
for U in ${URLS}; do
    OUT="raw/${SAMPLE}_R${i}.fastq.gz"
    echo "[task ${SLURM_ARRAY_TASK_ID}] -> https://${U}"
    curl -fsSL --retry 3 --retry-delay 30 --max-time 3600 -o "${OUT}" "https://${U}"
    [ -s "${OUT}" ] || { echo "ERROR: download failed: ${U}"; exit 1; }
    i=$((i+1))
done

ls -lh raw/${SAMPLE}_R*.fastq.gz
```

**Submit (capture the job ID so downstream steps can depend on it):**

```bash
# Make sure the env is active in THIS shell (once per login session)
micromamba activate reseq_bch709
cd ~/scratch/reseq

DL_JID=$(sbatch --parsable scripts/01_download.sh)
echo "Download job: ${DL_JID}"
squeue -u $USER   # you should see 2 tasks (download_1, download_2) running in parallel
```

Expected output of `logs/01_download_<jobid>_1.out` (sample1 task — sample2's log is similar with SRR519586):

```
[task 1] Downloading sample1 (SRR519585) from ENA
[task 1] -> https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR519/SRR519585/SRR519585_1.fastq.gz
[task 1] -> https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR519/SRR519585/SRR519585_2.fastq.gz
-rw-r--r-- 1 <netid> rc-bch709-6 2.0G <date> raw/sample1_R1.fastq.gz
-rw-r--r-- 1 <netid> rc-bch709-6 2.0G <date> raw/sample1_R2.fastq.gz
```

Sample2 ends up around 2.6 G per mate (SRR519586 is a slightly deeper run). If either file shows `0` bytes, ENA hiccupped — re-submit just that array task with `sbatch --array=2 scripts/01_download.sh`.

| `#SBATCH` option | Purpose |
|-----------------|---------|
| `--array=1-2` | Launches 2 tasks, each with its own `$SLURM_ARRAY_TASK_ID` (1 and 2) |
| `%A_%a` in log name | `%A`=array job ID, `%a`=task index — keeps each task's log separate |

> **Scaling up:** for 100 samples, change `--array=1-2` to `--array=1-100`. Slurm will run as many in parallel as your QoS allows.
{: .callout}

---

## 2. Download and Index Reference (Single Job)

The reference is shared across all samples, so we only build it once.

### `scripts/02_reference.sh`

```bash
#!/bin/bash
#SBATCH --job-name=reference
#SBATCH --account=cpu-s5-bch709-6
#SBATCH --partition=cpu-core-0
#SBATCH --cpus-per-task=4
#SBATCH --mem=16g
#SBATCH --time=06:00:00
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=<YOUR_EMAIL>
#SBATCH -o logs/02_reference_%j.out

set -euo pipefail
# Environment is activated in the login shell before `sbatch` is called —
# the PATH (with bwa-mem2, samtools, gatk, …) is inherited automatically.

cd ~/scratch/reseq

# ---- Download TAIR10 from TAIR (www.arabidopsis.org) ----
# Primary: TAIR (canonical source) — uses self-signed cert, so curl needs -k
# Fallback: EBI Ensembl Plants mirror (most stable mirror; NCBI is last resort)
TAIR_URL="https://www.arabidopsis.org/api/download-files/download?filePath=Genes/TAIR10_genome_release/TAIR10_chromosome_files/TAIR10_chr_all.fas.gz"
EBI_URL="https://ftp.ebi.ac.uk/ensemblgenomes/pub/plants/release-60/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz"
NCBI_URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/735/GCF_000001735.4_TAIR10.1/GCF_000001735.4_TAIR10.1_genomic.fna.gz"

OUT=reference.fasta.gz
rm -f "$OUT"
for URL in "$TAIR_URL" "$EBI_URL" "$NCBI_URL"; do
    echo "[ref] trying $URL"
    if curl -kfsSL --retry 3 --max-time 900 -o "$OUT" "$URL" && [ -s "$OUT" ]; then
        echo "[ref] downloaded from $URL"
        break
    fi
    rm -f "$OUT"
done
[ -s "$OUT" ] || { echo "ERROR: TAIR10 download failed from all sources"; exit 1; }
gunzip -f "$OUT"

# Normalize chromosome names to match the 1001genomes VCF (1..5, Mt, Pt):
#   - TAIR fasta uses Chr1..Chr5, ChrM, ChrC
#   - NCBI RefSeq uses NC_003070.9 …
#   - EBI Ensembl already uses 1..5, Mt, Pt (no rename needed)
if grep -q '^>Chr' reference.fasta; then
    sed -i -E 's/^>Chr([0-9]+).*/>\1/; s/^>ChrM.*/>Mt/; s/^>ChrC.*/>Pt/' reference.fasta
elif grep -q '^>NC_' reference.fasta; then
    awk 'BEGIN{
        m["NC_003070.9"]="1"; m["NC_003071.7"]="2"; m["NC_003074.8"]="3"
        m["NC_003075.7"]="4"; m["NC_003076.8"]="5"
        m["NC_037304.1"]="Mt"; m["NC_000932.1"]="Pt"
    }
    /^>/{
        name=substr($1,2); print ">" (name in m ? m[name] : name); next
    }
    {print}' reference.fasta > reference.renamed.fa && mv -f reference.renamed.fa reference.fasta
fi

# Build all indices
bwa-mem2 index reference.fasta
samtools faidx reference.fasta
picard CreateSequenceDictionary R=reference.fasta O=reference.dict

# Download and prepare known variants (sites-only to avoid malformed VCF)
# 1001genomes.org is sometimes slow — give curl up to 30 min and retry 3×
KSITES_URL="https://1001genomes.org/data/GMI-MPI/releases/v3.1/1001genomes_snp-short-indel_only_ACGTN.vcf.gz"
# The 1001genomes VCF is ~19 GB; on Pronghorn outbound (~3 MB/s) the download
# takes ~3 hours. We use `-C -` to resume on retry instead of restarting at
# byte 0, and drop --max-time so a single attempt isn't capped (Slurm --time
# already bounds the whole job).
curl -fsSL --retry 5 --retry-delay 30 -C - \
    -o 1001genomes_snp-short-indel_only_ACGTN.vcf.gz "${KSITES_URL}"
# Atomic rename: write to a .tmp file, only promote to known_sites.vcf.gz once
# bgzip + tabix BOTH succeed and outputs are non-empty. If anything in this
# block fails, the original 1001genomes file is preserved so the student can
# re-run step 02 without re-downloading 19 GB.
zcat 1001genomes_snp-short-indel_only_ACGTN.vcf.gz \
    | awk 'BEGIN{OFS="\t"} /^##/{print; next} /^#CHROM/{print $1,$2,$3,$4,$5,$6,$7,$8; next} {print $1,$2,$3,$4,$5,$6,$7,$8}' \
    | bgzip > known_sites.vcf.gz.tmp
tabix -p vcf known_sites.vcf.gz.tmp
[ -s known_sites.vcf.gz.tmp ] && [ -s known_sites.vcf.gz.tmp.tbi ] || {
    echo "ERROR: known_sites.vcf.gz build failed (empty output). Re-run step 02." >&2
    rm -f known_sites.vcf.gz.tmp known_sites.vcf.gz.tmp.tbi
    exit 1
}
mv -f known_sites.vcf.gz.tmp known_sites.vcf.gz
mv -f known_sites.vcf.gz.tmp.tbi known_sites.vcf.gz.tbi
rm -f 1001genomes_snp-short-indel_only_ACGTN.vcf.gz

echo "Reference prep done."
```

**Submit** (independent of download, so it can run in parallel with Step 1):

```bash
REF_JID=$(sbatch --parsable scripts/02_reference.sh)
echo "Reference job: ${REF_JID}"
```

Expected output of `logs/02_reference_<jobid>.out` (truncated — first lines + tail):

```
[ref] trying https://www.arabidopsis.org/api/download-files/download?filePath=Genes/TAIR10_genome_release/TAIR10_chromosome_files/TAIR10_chr_all.fas.gz
[ref] downloaded from https://www.arabidopsis.org/api/download-files/download?filePath=...TAIR10_chr_all.fas.gz
Looking to launch executable ".../bwa-mem2.avx2", simd = .avx2
[bwa_index] Pack FASTA... 0.55 sec
* Entering FMI_search
ref seq len = 239337268
build suffix-array ticks = 62781215970
build fm-index ticks = 13953625410
Total time taken: 41.8652
INFO  <date>  CreateSequenceDictionary
...
[<date>] CreateSequenceDictionary OUTPUT=reference.dict REFERENCE=reference.fasta ...
[<date>] picard.sam.CreateSequenceDictionary done. Elapsed time: 0.02 minutes.
Reference prep done.
```

Quick sanity check after the job finishes — you should see the FASTA, the bwa-mem2 index files, the `.fai`, and the Picard `.dict`:

```bash
ls -lh ~/scratch/reseq/reference.*
```

Expected output:

```
-rw-r--r-- 1 <netid> rc-bch709-6 116M <date> reference.fasta
-rw-r--r-- 1 <netid> rc-bch709-6 229M <date> reference.fasta.0123
-rw-r--r-- 1 <netid> rc-bch709-6 7.4K <date> reference.fasta.amb
-rw-r--r-- 1 <netid> rc-bch709-6  230 <date> reference.fasta.ann
-rw-r--r-- 1 <netid> rc-bch709-6 371M <date> reference.fasta.bwt.2bit.64
-rw-r--r-- 1 <netid> rc-bch709-6  175 <date> reference.fasta.fai
-rw-r--r-- 1 <netid> rc-bch709-6  29M <date> reference.fasta.pac
-rw-r--r-- 1 <netid> rc-bch709-6 1.0K <date> reference.dict
```

After `known_sites.vcf.gz` finishes building (the long curl + bgzip + tabix step at the end of the script), confirm both the VCF and its tabix index exist:

```bash
ls -lh ~/scratch/reseq/known_sites.vcf.gz*
```

Expected output:

```
-rw-r--r-- 1 <netid> rc-bch709-6  52M <date> known_sites.vcf.gz
-rw-r--r-- 1 <netid> rc-bch709-6 220K <date> known_sites.vcf.gz.tbi
```

If the `.tbi` is missing, re-run step 02 — `BaseRecalibrator` in step 04 will fail without it.

---

## 3. Trim and Align in parallel (Array Job with Dependency)

`fastp` and `bwa-mem2` each accept multiple threads. We give each array task 8 cores to cut wall time.

### `scripts/03_align.sh`

```bash
#!/bin/bash
#SBATCH --job-name=align
#SBATCH --account=cpu-s5-bch709-6
#SBATCH --partition=cpu-core-0
#SBATCH --cpus-per-task=8
#SBATCH --mem=32g
#SBATCH --time=12:00:00
#SBATCH --array=1-2
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=<YOUR_EMAIL>
#SBATCH -o logs/03_align_%A_%a.out

set -euo pipefail
# Environment is activated in the login shell before `sbatch` is called —
# the PATH (with bwa-mem2, samtools, gatk, …) is inherited automatically.

cd ~/scratch/reseq

LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" samples.tsv)
SAMPLE=$(echo "$LINE" | cut -f1)

# ---- Trim ----
fastp \
    --in1  raw/${SAMPLE}_R1.fastq.gz \
    --in2  raw/${SAMPLE}_R2.fastq.gz \
    --out1 trim/${SAMPLE}_R1.trimmed.fq.gz \
    --out2 trim/${SAMPLE}_R2.trimmed.fq.gz \
    --detect_adapter_for_pe \
    --qualified_quality_phred 20 \
    --length_required 50 \
    --thread ${SLURM_CPUS_PER_TASK} \
    --html trim/${SAMPLE}_fastp.html \
    --json trim/${SAMPLE}_fastp.json

# ---- Align + sort ----
set -o pipefail
bwa-mem2 mem \
    -t ${SLURM_CPUS_PER_TASK} \
    -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:ILLUMINA\tLB:lib_${SAMPLE}\tPU:unit_${SAMPLE}" \
    reference.fasta \
    trim/${SAMPLE}_R1.trimmed.fq.gz \
    trim/${SAMPLE}_R2.trimmed.fq.gz \
  | samtools sort -@ ${SLURM_CPUS_PER_TASK} -T /tmp/${SAMPLE}_sort \
      -o bam/${SAMPLE}.bam

samtools index -@ ${SLURM_CPUS_PER_TASK} bam/${SAMPLE}.bam
samtools quickcheck bam/${SAMPLE}.bam && echo "BAM OK"
```

Expected output of `logs/03_align_<jobid>_1.out` (sample1 — fastp summary at the top, bwa-mem2 + samtools sort at the bottom; truncated):

```
Detecting adapter sequence for read1...
>TruSeq2_PE_r
AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAG

Read1 before filtering:
total reads: 29418122
total bases: 2235777272
Q20 bases: 2090582549(93.5059%)
Q30 bases: 1935916821(86.5881%)

Read1 after filtering:
total reads: 26525643
total bases: 2003802159
Q20 bases: 1942580727(96.9447%)
Q30 bases: 1815798338(90.6176%)

Filtering result:
reads passed filter: 53051286
reads failed due to low quality: 4575222
reads failed due to too many N: 64736
reads failed due to too short: 1145000
# ...
fastp v1.3.3, time used: 137 seconds
Looking to launch executable ".../bwa-mem2.avx2", simd = .avx2
	Time taken for main_mem function: 1797.12 sec
[bam_sort_core] merging from 2 files and 8 in-memory blocks...
BAM OK
```

The `trim/sample1_fastp.json` summary captures the same numbers in a parseable form (the HTML report is the human-friendly view of the same data):

```bash
ls -lh ~/scratch/reseq/trim/sample1_fastp.* ~/scratch/reseq/bam/sample1.bam
```

Expected output:

```
-rw-r--r-- 1 <netid> rc-bch709-6 444K <date> trim/sample1_fastp.html
-rw-r--r-- 1 <netid> rc-bch709-6 109K <date> trim/sample1_fastp.json
-rw-r--r-- 1 <netid> rc-bch709-6 2.9G <date> bam/sample1.bam
```

A quick `samtools flagstat` confirms the BAM looks reasonable (>98% of reads mapped, ~91% properly paired):

```bash
samtools flagstat ~/scratch/reseq/bam/sample1.bam
```

Expected output:

```
53076354 + 0 in total (QC-passed reads + QC-failed reads)
53051286 + 0 primary
0 + 0 secondary
25068 + 0 supplementary
0 + 0 duplicates
52354774 + 0 mapped (98.64% : N/A)
52329706 + 0 primary mapped (98.64% : N/A)
53051286 + 0 paired in sequencing
26525643 + 0 read1
26525643 + 0 read2
48411760 + 0 properly paired (91.25% : N/A)
52276894 + 0 with itself and mate mapped
52812 + 0 singletons (0.10% : N/A)
490848 + 0 with mate mapped to a different chr
129225 + 0 with mate mapped to a different chr (mapQ>=5)
```

**Submit with dependencies** — wait for **both** download (raw FASTQ) and reference prep (index) to succeed before aligning:

```bash
ALIGN_JID=$(sbatch --parsable --dependency=afterok:${DL_JID}:${REF_JID} scripts/03_align.sh)
echo "Align: ${ALIGN_JID}"
squeue -u $USER   # ALIGN_JID should show state PD with reason (Dependency)
```

If you skipped capturing `DL_JID` or `REF_JID` earlier, you can look them up with `squeue -u $USER` or just depend on whichever you know — Slurm only needs valid predecessor IDs.

> **`${SLURM_CPUS_PER_TASK}`** is automatically set by Slurm to match `--cpus-per-task`. Using it everywhere means changing one `#SBATCH` line automatically updates every thread count below.
{: .callout}

> **`-T /tmp/${SAMPLE}_sort`** puts `samtools sort` temp files on each node's local disk (fast) instead of shared scratch (slower, and can fill up).
{: .callout}

---

## 4. Mark Duplicates + BQSR (Array Job)

MarkDuplicates and BQSR aren't heavily multi-threaded, but we still get per-sample parallelism.

### `scripts/04_markdup_bqsr.sh`

```bash
#!/bin/bash
#SBATCH --job-name=markdup_bqsr
#SBATCH --account=cpu-s5-bch709-6
#SBATCH --partition=cpu-core-0
#SBATCH --cpus-per-task=4
#SBATCH --mem=16g
#SBATCH --time=08:00:00
#SBATCH --array=1-2
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=<YOUR_EMAIL>
#SBATCH -o logs/04_markdup_bqsr_%A_%a.out

set -euo pipefail
# Environment is activated in the login shell before `sbatch` is called —
# the PATH (with bwa-mem2, samtools, gatk, …) is inherited automatically.

cd ~/scratch/reseq

SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" samples.tsv | cut -f1)

# ---- MarkDuplicates ----
picard -Xmx12g MarkDuplicates \
    I=bam/${SAMPLE}.bam \
    O=bam/${SAMPLE}.markdup.bam \
    M=bam/${SAMPLE}.markdup.metrics \
    VALIDATION_STRINGENCY=SILENT
samtools index -@ ${SLURM_CPUS_PER_TASK} bam/${SAMPLE}.markdup.bam

# ---- Preflight: known_sites.vcf.gz must exist (built by step 02) ----
if [ ! -s known_sites.vcf.gz ] || [ ! -s known_sites.vcf.gz.tbi ]; then
    echo "ERROR: known_sites.vcf.gz (and .tbi) not found in $(pwd)." >&2
    echo "       Step 02 (02_reference.sh) did not finish — re-run it before step 04." >&2
    ls -lh 1001genomes_*.vcf.gz known_sites.vcf.gz* 2>/dev/null || true
    exit 1
fi

# ---- BaseRecalibrator ----
gatk --java-options "-Xmx12g" BaseRecalibrator \
    -I bam/${SAMPLE}.markdup.bam \
    -R reference.fasta \
    --known-sites known_sites.vcf.gz \
    -O bam/${SAMPLE}.recal.table

# ---- ApplyBQSR ----
gatk --java-options "-Xmx12g" ApplyBQSR \
    -I bam/${SAMPLE}.markdup.bam \
    -R reference.fasta \
    --bqsr-recal-file bam/${SAMPLE}.recal.table \
    -O bam/${SAMPLE}.recal.bam
samtools index -@ ${SLURM_CPUS_PER_TASK} bam/${SAMPLE}.recal.bam

# Clean up intermediate markdup BAM
rm -f bam/${SAMPLE}.bam bam/${SAMPLE}.bam.bai
rm -f bam/${SAMPLE}.markdup.bam bam/${SAMPLE}.markdup.bam.bai
```

Expected tail of `logs/04_markdup_bqsr_<jobid>_1.out` — Picard wraps up MarkDuplicates, then GATK starts BaseRecalibrator and ApplyBQSR (truncated):

```
INFO	<date>	MarkDuplicates	Marking 14734338 records as duplicates.
INFO	<date>	MarkDuplicates	Found 0 optical duplicate clusters.
INFO	<date>	MarkDuplicates	Written    50,000,000 records.  Elapsed time: 00:08:35s.
INFO	<date>	MarkDuplicates	Writing complete. Closing input iterator.
[<date>] picard.sam.markduplicates.MarkDuplicates done. Elapsed time: 13.21 minutes.
Using GATK jar .../gatk-package-4.6.2.0-local.jar
Running:
    java ... -Xmx12g ... BaseRecalibrator -I bam/sample1.markdup.bam -R reference.fasta --known-sites known_sites.vcf.gz -O bam/sample1.recal.table
INFO  BaseRecalibrator - The Genome Analysis Toolkit (GATK) v4.6.2.0
INFO  BaseRecalibrationEngine - The covariates being used here:
INFO  BaseRecalibrationEngine - 	ReadGroupCovariate
INFO  BaseRecalibrationEngine - 	QualityScoreCovariate
INFO  BaseRecalibrationEngine - 	ContextCovariate
INFO  BaseRecalibrationEngine - 	CycleCovariate
INFO  ProgressMeter - Starting traversal
INFO  ProgressMeter -        Current Locus  Elapsed Minutes       Reads Processed     Reads/Minute
INFO  ProgressMeter -            1:1523668              0.2                371000        2221335.2
# ...
INFO  ProgressMeter -            done.            <X.X>              <NNNNNNNN>          <rate>
INFO  BaseRecalibrator - Finished. Elapsed time: <NN.NN> minutes.
INFO  ApplyBQSR - Done. Elapsed time: <NN.NN> minutes.
```

Picard writes a metrics file — open `bam/sample1.markdup.metrics` to see the duplication rate per library (the `## METRICS CLASS` block is the part GATK and MultiQC parse):

```bash
head -10 ~/scratch/reseq/bam/sample1.markdup.metrics
```

Expected output:

```
## htsjdk.samtools.metrics.StringHeader
# MarkDuplicates INPUT=[bam/sample1.bam] OUTPUT=bam/sample1.markdup.bam METRICS_FILE=bam/sample1.markdup.metrics ...
## htsjdk.samtools.metrics.StringHeader
# Started on: <date>

## METRICS CLASS	picard.sam.DuplicationMetrics
LIBRARY	UNPAIRED_READS_EXAMINED	READ_PAIRS_EXAMINED	SECONDARY_OR_SUPPLEMENTARY_RDS	UNMAPPED_READS	UNPAIRED_READ_DUPLICATES	READ_PAIR_DUPLICATES	READ_PAIR_OPTICAL_DUPLICATES	PERCENT_DUPLICATION	ESTIMATED_LIBRARY_SIZE
lib_sample1	52812	26138447	25068	721580	16888	7358725	0	0.281567	37224722
```

`PERCENT_DUPLICATION` of ~0.28 (28%) is on the high side — these public SRA libraries were sequenced deeply on a small genome, so some duplication is expected. After BQSR finishes, the recalibrated BAM is what every downstream step uses:

```bash
ls -lh ~/scratch/reseq/bam/sample1.recal.bam*
```

Expected output:

```
# example output (your numbers will differ)
-rw-r--r-- 1 <netid> rc-bch709-6 2.7G <date> bam/sample1.recal.bam
-rw-r--r-- 1 <netid> rc-bch709-6 350K <date> bam/sample1.recal.bam.bai
```

**Submit:**

```bash
MARKDUP_JID=$(sbatch --parsable --dependency=afterok:${ALIGN_JID} scripts/04_markdup_bqsr.sh)
echo "MarkDup+BQSR: ${MARKDUP_JID}"
```

---

## 5. HaplotypeCaller — Scatter by Chromosome (2-D Array)

HaplotypeCaller is the slowest step. For Arabidopsis (7 chromosomes × N samples), we use a **2-dimensional array**: each task handles one `(sample, chromosome)` pair. This is the biggest wall-time win.

For N=2 samples and 7 chromosomes (1, 2, 3, 4, 5, Mt, Pt) → **14 parallel tasks** instead of 2 serial ones.

### `scripts/05a_haplotypecaller.sh`

```bash
#!/bin/bash
#SBATCH --job-name=hc_scatter
#SBATCH --account=cpu-s5-bch709-6
#SBATCH --partition=cpu-core-0
#SBATCH --cpus-per-task=4
#SBATCH --mem=16g
#SBATCH --time=12:00:00
#SBATCH --array=1-14                  # N_SAMPLES * N_CHROMS = 2 * 7
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=<YOUR_EMAIL>
#SBATCH -o logs/05a_hc_%A_%a.out

set -euo pipefail
# Environment is activated in the login shell before `sbatch` is called —
# the PATH (with bwa-mem2, samtools, gatk, …) is inherited automatically.

cd ~/scratch/reseq

# Map the linear task ID to (sample_index, chr_index)
CHROMS=(1 2 3 4 5 Mt Pt)
N_CHROMS=${#CHROMS[@]}

# 0-indexed within the array
IDX=$((SLURM_ARRAY_TASK_ID - 1))
SAMPLE_IDX=$((IDX / N_CHROMS + 1))     # 1-indexed line in samples.tsv
CHR_IDX=$((IDX % N_CHROMS))

SAMPLE=$(sed -n "${SAMPLE_IDX}p" samples.tsv | cut -f1)
CHR=${CHROMS[$CHR_IDX]}

echo "Task ${SLURM_ARRAY_TASK_ID}: sample=${SAMPLE}, chr=${CHR}"

mkdir -p vcf/scatter

gatk --java-options "-Xmx12g" HaplotypeCaller \
    -R reference.fasta \
    -I bam/${SAMPLE}.recal.bam \
    -O vcf/scatter/${SAMPLE}.${CHR}.g.vcf.gz \
    -L ${CHR} \
    -ERC GVCF \
    --native-pair-hmm-threads ${SLURM_CPUS_PER_TASK}
```

Expected tail of `logs/05a_hc_<jobid>_1.out` (sample1, chromosome 1 — the GVCF mode is verbose; the key lines are the ProgressMeter and the final summary):

```
# example output (your numbers will differ)
INFO  HaplotypeCaller - The Genome Analysis Toolkit (GATK) v4.6.2.0
INFO  HaplotypeCaller - HTSJDK Version: 4.2.0
INFO  HaplotypeCaller - Picard Version: 3.4.0
INFO  HaplotypeCallerEngine - Tool is in reference confidence mode and the annotation, the following changes will be made to any specified annotations: 'StrandBiasBySample' will be enabled. 'ChromosomeCounts', 'FisherStrand', 'StrandOddsRatio' and 'QualByDepth' annotations have been disabled
INFO  HaplotypeCaller - Defragmenting 100 events for 1 samples
INFO  ProgressMeter - Starting traversal
INFO  ProgressMeter -        Current Locus  Elapsed Minutes       Reads Processed     Reads/Minute
INFO  ProgressMeter -            1:1523668              0.5                 53000          106000.0
# ...
INFO  ProgressMeter -           1:30425192             14.8               2920000          197297.3
INFO  ProgressMeter -            traversal complete. Processed <N> total regions in <M.M> minutes.
INFO  HaplotypeCaller - Shutting down engine
INFO  HaplotypeCaller - Done. Elapsed time: <NN.NN> minutes.
```

After the array finishes, you should see 14 per-`(sample, chromosome)` GVCFs (2 samples × 7 chromosomes):

```bash
ls -lh ~/scratch/reseq/vcf/scatter/
```

Expected output:

```
# example output (your numbers will differ)
-rw-r--r-- 1 <netid> rc-bch709-6 12M <date> sample1.1.g.vcf.gz
-rw-r--r-- 1 <netid> rc-bch709-6 280K <date> sample1.1.g.vcf.gz.tbi
-rw-r--r-- 1 <netid> rc-bch709-6 8.0M <date> sample1.2.g.vcf.gz
# ... one pair per (sample, chr) — 14 GVCFs + 14 .tbi total
-rw-r--r-- 1 <netid> rc-bch709-6 1.1M <date> sample2.Pt.g.vcf.gz
-rw-r--r-- 1 <netid> rc-bch709-6  18K <date> sample2.Pt.g.vcf.gz.tbi
```

### `scripts/05b_gather_gvcf.sh` — Gather per-chromosome GVCFs into one per-sample GVCF

```bash
#!/bin/bash
#SBATCH --job-name=hc_gather
#SBATCH --account=cpu-s5-bch709-6
#SBATCH --partition=cpu-core-0
#SBATCH --cpus-per-task=2
#SBATCH --mem=8g
#SBATCH --time=02:00:00
#SBATCH --array=1-2                   # one task per sample
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=<YOUR_EMAIL>
#SBATCH -o logs/05b_gather_%A_%a.out

set -euo pipefail
# Environment is activated in the login shell before `sbatch` is called —
# the PATH (with bwa-mem2, samtools, gatk, …) is inherited automatically.

cd ~/scratch/reseq

SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" samples.tsv | cut -f1)

# Merge the 7 chromosome GVCFs for this sample into one
CHROMS=(1 2 3 4 5 Mt Pt)
INPUTS=""
for CHR in "${CHROMS[@]}"; do
    INPUTS+=" -I vcf/scatter/${SAMPLE}.${CHR}.g.vcf.gz"
done

gatk --java-options "-Xmx6g" GatherVcfs \
    ${INPUTS} \
    -O vcf/${SAMPLE}.g.vcf.gz

gatk IndexFeatureFile -I vcf/${SAMPLE}.g.vcf.gz

# Verify and clean up scatter files
ls -lh vcf/${SAMPLE}.g.vcf.gz
rm -rf vcf/scatter/${SAMPLE}.*.g.vcf.gz vcf/scatter/${SAMPLE}.*.g.vcf.gz.tbi
```

Expected tail of `logs/05b_gather_<jobid>_1.out` — `GatherVcfs` is fast (it just concatenates already-sorted VCFs) and `IndexFeatureFile` writes a `.tbi`:

```
# example output (your numbers will differ)
INFO  GatherVcfs - Checking inputs.
INFO  GatherVcfs - Gathering by copying gzip blocks. Will not be indexed.
INFO  GatherVcfs - Done.
[<date>] picard.vcf.GatherVcfs done. Elapsed time: 0.05 minutes.
INFO  IndexFeatureFile - Successfully wrote index to vcf/sample1.g.vcf.gz.tbi
[<date>] org.broadinstitute.hellbender.tools.IndexFeatureFile done. Elapsed time: 0.10 minutes.
-rw-r--r-- 1 <netid> rc-bch709-6  90M <date> vcf/sample1.g.vcf.gz
```

**Submit with chained dependencies:**

```bash
HC_JID=$(sbatch --parsable --dependency=afterok:${MARKDUP_JID} scripts/05a_haplotypecaller.sh)
GATHER_JID=$(sbatch --parsable --dependency=afterok:${HC_JID} scripts/05b_gather_gvcf.sh)
echo "HaplotypeCaller: ${HC_JID} | Gather: ${GATHER_JID}"
```

> **Array size formula:** `N_SAMPLES × N_CHROMS`. For 10 samples × 7 chroms, use `--array=1-70`. For human (25 chroms) × 10 samples, `--array=1-250`.
{: .callout}

> **Why scatter-gather beats per-sample:** HaplotypeCaller on a whole Arabidopsis genome takes ~60 min. Per chromosome it takes 5-15 min. 7 chromosomes in parallel finish in ~15 min (limited by the slowest chromosome). For humans the speedup is much larger.
{: .callout}

---

## 6. Joint Genotyping (Single Job)

Once all per-sample GVCFs exist, combine them and run joint genotyping. This step is serial but fast (minutes).

### `scripts/06_joint_genotype.sh`

```bash
#!/bin/bash
#SBATCH --job-name=joint_geno
#SBATCH --account=cpu-s5-bch709-6
#SBATCH --partition=cpu-core-0
#SBATCH --cpus-per-task=4
#SBATCH --mem=32g
#SBATCH --time=04:00:00
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=<YOUR_EMAIL>
#SBATCH -o logs/06_joint_%j.out

set -euo pipefail
# Environment is activated in the login shell before `sbatch` is called —
# the PATH (with bwa-mem2, samtools, gatk, …) is inherited automatically.

cd ~/scratch/reseq

# Build the -V argument list from samples.tsv (auto-scales with sample count)
INPUTS=""
while read SAMPLE SRR; do
    INPUTS+=" -V vcf/${SAMPLE}.g.vcf.gz"
done < samples.tsv

gatk --java-options "-Xmx24g" CombineGVCFs \
    -R reference.fasta \
    ${INPUTS} \
    -O vcf/cohort.g.vcf.gz

gatk --java-options "-Xmx24g" GenotypeGVCFs \
    -R reference.fasta \
    -V vcf/cohort.g.vcf.gz \
    -O vcf/cohort.vcf.gz

echo "Joint genotyping done. Variants:"
bcftools stats vcf/cohort.vcf.gz | grep "^SN"
```

Expected tail of `logs/06_joint_<jobid>.out` — `CombineGVCFs` then `GenotypeGVCFs`, ending with the `bcftools stats` summary:

```
# example output (your numbers will differ)
INFO  CombineGVCFs - The Genome Analysis Toolkit (GATK) v4.6.2.0
INFO  ProgressMeter - Starting traversal
INFO  ProgressMeter -            traversal complete. Processed <N> total variants in <M.M> minutes.
INFO  CombineGVCFs - Shutting down engine
INFO  CombineGVCFs - Done. Elapsed time: <NN.NN> minutes.
INFO  GenotypeGVCFs - The Genome Analysis Toolkit (GATK) v4.6.2.0
INFO  ProgressMeter - Starting traversal
INFO  ProgressMeter -            traversal complete. Processed <N> total variants in <M.M> minutes.
INFO  GenotypeGVCFs - Done. Elapsed time: <NN.NN> minutes.
Joint genotyping done. Variants:
SN	0	number of samples:	2
SN	0	number of records:	<NNNNNNN>
SN	0	number of no-ALTs:	0
SN	0	number of SNPs:	<NNNNNNN>
SN	0	number of MNPs:	0
SN	0	number of indels:	<NNNNNN>
SN	0	number of others:	0
SN	0	number of multiallelic sites:	<NNNNN>
```

For Arabidopsis with 2 deeply-sequenced samples, expect on the order of ~1–2 M raw cohort variants before filtering. You can confirm with a one-liner once the job completes:

```bash
bcftools view -H ~/scratch/reseq/vcf/cohort.vcf.gz | wc -l
```

**Submit:**

```bash
JOINT_JID=$(sbatch --parsable --dependency=afterok:${GATHER_JID} scripts/06_joint_genotype.sh)
```

> **For >100 samples** use `GenomicsDBImport` instead of `CombineGVCFs` — it's dramatically faster for large cohorts. See [GATK docs on joint calling](https://gatk.broadinstitute.org/hc/en-us/articles/360035890411).
{: .callout}

---

## 7. Filter, Annotate, and Population Analysis (Single Job)

Everything from here is fast and has no further parallelism to exploit.

### `scripts/07_filter_annotate.sh`

```bash
#!/bin/bash
#SBATCH --job-name=filter_ann
#SBATCH --account=cpu-s5-bch709-6
#SBATCH --partition=cpu-core-0
#SBATCH --cpus-per-task=4
#SBATCH --mem=16g
#SBATCH --time=04:00:00
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=<YOUR_EMAIL>
#SBATCH -o logs/07_filter_%j.out

set -euo pipefail
# Environment is activated in the login shell before `sbatch` is called —
# the PATH (with bwa-mem2, samtools, gatk, …) is inherited automatically.

cd ~/scratch/reseq

# ---- Split SNPs / Indels ----
gatk SelectVariants -R reference.fasta -V vcf/cohort.vcf.gz \
    --select-type-to-include SNP   -O vcf/cohort.snps.vcf.gz
gatk SelectVariants -R reference.fasta -V vcf/cohort.vcf.gz \
    --select-type-to-include INDEL -O vcf/cohort.indels.vcf.gz

# ---- Hard filter ----
gatk VariantFiltration -R reference.fasta -V vcf/cohort.snps.vcf.gz \
    --filter-expression "QD < 2.0"               --filter-name "QD2" \
    --filter-expression "FS > 60.0"              --filter-name "FS60" \
    --filter-expression "MQ < 40.0"              --filter-name "MQ40" \
    --filter-expression "MQRankSum < -12.5"      --filter-name "MQRankSum-12.5" \
    --filter-expression "ReadPosRankSum < -8.0"  --filter-name "ReadPosRankSum-8" \
    -O vcf/cohort.snps.filtered.vcf.gz

gatk VariantFiltration -R reference.fasta -V vcf/cohort.indels.vcf.gz \
    --filter-expression "QD < 2.0"                --filter-name "QD2" \
    --filter-expression "FS > 200.0"              --filter-name "FS200" \
    --filter-expression "ReadPosRankSum < -20.0"  --filter-name "ReadPosRankSum-20" \
    -O vcf/cohort.indels.filtered.vcf.gz

# ---- PASS-only VCF ----
bcftools view -f PASS -Oz -o vcf/cohort.snps.pass.vcf.gz vcf/cohort.snps.filtered.vcf.gz
tabix -p vcf vcf/cohort.snps.pass.vcf.gz

# ---- SnpEff annotation ----
snpEff -dataDir ${HOME}/snpeff_data Arabidopsis_thaliana \
    vcf/cohort.snps.pass.vcf.gz > vcf/cohort.snps.annotated.vcf 2> logs/snpeff.log

# ---- PLINK (requires 2+ samples) ----
cd vcf
plink --vcf cohort.snps.pass.vcf.gz --make-bed --out cohort       --allow-extra-chr
plink --bfile cohort --pca 10                    --out cohort.pca --allow-extra-chr
plink --bfile cohort --indep-pairwise 50 10 0.2  --out cohort.ld  --allow-extra-chr

echo "Pipeline complete."
```

Expected tail of `logs/07_filter_<jobid>.out` — `VariantFiltration` annotates the FILTER column, `bcftools view -f PASS` keeps only the survivors, then snpEff and PLINK run:

```
# example output (your numbers will differ)
INFO  VariantFiltration - The Genome Analysis Toolkit (GATK) v4.6.2.0
INFO  ProgressMeter - Starting traversal
INFO  ProgressMeter -            traversal complete. Processed <N> total variants in <M.M> minutes.
INFO  VariantFiltration - Done. Elapsed time: <N.NN> minutes.
[snpEff] Reading database for genome: 'Arabidopsis_thaliana'
[snpEff] Analysis done.
Pipeline complete.
```

snpEff also writes a summary HTML/CSV next to the working directory (`snpEff_summary.html`, `snpEff_genes.txt`) — open the HTML in a browser, or peek at the summary table that snpEff prints to stderr / `logs/snpeff.log`:

```bash
head -25 ~/scratch/reseq/logs/snpeff.log
```

Expected output:

```
# example output (your numbers will differ — these are the section headings snpEff always prints)
Number_of_variants_processed	<NNNNNN>
Number_of_effects	<NNNNNNN>
# Effects by impact
HIGH		<NNNN>
LOW		<NNNNNN>
MODERATE	<NNNNN>
MODIFIER	<NNNNNNN>
# Effects by functional class
MISSENSE	<NNNNN>
NONSENSE	<NNN>
SILENT		<NNNNN>
# Effects by type
intron_variant		<NNNNNN>
synonymous_variant	<NNNNN>
missense_variant	<NNNNN>
stop_gained		<NNN>
# ...
```

The PLINK PCA writes `cohort.pca.eigenvec` (per-sample PC scores) and `cohort.pca.eigenval` (variance per PC). With only 2 samples PCA is degenerate — for real cohorts (8+ samples) the first 2 PCs already separate populations:

```bash
cat ~/scratch/reseq/vcf/cohort.pca.eigenvec
```

Expected output:

```
# example output — columns are: FID  IID  PC1  PC2  ...  PC10
sample1 sample1  -0.7071  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
sample2 sample2   0.7071  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
```

Finally, count the PASS-filter SNPs in your cohort VCF — this is the "headline number" your downstream analysis (GWAS, popgen, annotation) starts from:

```bash
bcftools view -H ~/scratch/reseq/vcf/cohort.snps.pass.vcf.gz | wc -l
bcftools stats ~/scratch/reseq/vcf/cohort.snps.pass.vcf.gz | grep "^SN"
```

Expected output:

```
# example output (your numbers will differ; for 2 deep Arabidopsis samples expect ~0.5–1.5 M PASS SNPs)
<NNNNNNN>
SN	0	number of samples:	2
SN	0	number of records:	<NNNNNNN>
SN	0	number of SNPs:	<NNNNNNN>
SN	0	number of indels:	0
SN	0	number of multiallelic sites:	<NNNNN>
SN	0	number of multiallelic SNP sites:	<NNNNN>
```

**Submit:**

```bash
FILTER_JID=$(sbatch --parsable --dependency=afterok:${JOINT_JID} scripts/07_filter_annotate.sh)
echo "Filter+Annotate: ${FILTER_JID}"
```

---

## 8. Submit the Entire Pipeline with One Script

Instead of manually chaining each step, use a driver script that submits everything with correct dependencies in one go.

### `scripts/run_all.sh`

```bash
#!/bin/bash
set -euo pipefail
cd ~/scratch/reseq

# Activate the env in THIS shell so every sbatch below inherits the PATH
# (sbatch --export=ALL is the default — the submitted jobs see the same tools)
export MAMBA_ROOT_PREFIX="${MAMBA_ROOT_PREFIX:-$HOME/micromamba}"
eval "$(micromamba shell hook --shell=bash)"
micromamba activate reseq_bch709

# --parsable returns just the job ID, perfect for chaining
DL_JID=$(sbatch --parsable                                  scripts/01_download.sh)
REF_JID=$(sbatch --parsable                                 scripts/02_reference.sh)
ALIGN_JID=$(sbatch --parsable --dependency=afterok:${DL_JID}:${REF_JID} scripts/03_align.sh)
MARK_JID=$(sbatch --parsable --dependency=afterok:${ALIGN_JID}          scripts/04_markdup_bqsr.sh)
HC_JID=$(sbatch --parsable --dependency=afterok:${MARK_JID}             scripts/05a_haplotypecaller.sh)
GATHER_JID=$(sbatch --parsable --dependency=afterok:${HC_JID}           scripts/05b_gather_gvcf.sh)
JOINT_JID=$(sbatch --parsable --dependency=afterok:${GATHER_JID}        scripts/06_joint_genotype.sh)
FILT_JID=$(sbatch --parsable --dependency=afterok:${JOINT_JID}          scripts/07_filter_annotate.sh)

echo "Submitted pipeline:"
printf '  %-20s %s\n' "download"       "${DL_JID}"
printf '  %-20s %s\n' "reference"      "${REF_JID}"
printf '  %-20s %s\n' "align"          "${ALIGN_JID}"
printf '  %-20s %s\n' "markdup+bqsr"   "${MARK_JID}"
printf '  %-20s %s\n' "haplotypecaller" "${HC_JID}"
printf '  %-20s %s\n' "gather"         "${GATHER_JID}"
printf '  %-20s %s\n' "joint_genotype" "${JOINT_JID}"
printf '  %-20s %s\n' "filter+annotate" "${FILT_JID}"

echo ""
echo "Monitor with:  squeue -u \$USER"
echo "Cancel all:    scancel ${DL_JID} ${REF_JID} ${ALIGN_JID} ${MARK_JID} ${HC_JID} ${GATHER_JID} ${JOINT_JID} ${FILT_JID}"
```

**Run it:**

```bash
chmod +x scripts/*.sh
bash scripts/run_all.sh
```

You'll see 8 job IDs printed. Slurm queues them; each waits for its predecessor to succeed before starting.

| Dependency type | Meaning |
|----------------|---------|
| `afterok:JID` | Start only if `JID` finished successfully |
| `afterany:JID` | Start after `JID` finishes, regardless of exit status |
| `afternotok:JID` | Start only if `JID` failed (for error-handling scripts) |
| `singleton` | Run after all other jobs with the same name finish |

### 🧑‍💻 Hands-on walkthrough — submit the pipeline step-by-step

If you want to see exactly what `run_all.sh` does (or debug one step), you can submit each stage manually. Every `sbatch` returns a **job ID** that the next step depends on.

**Do this first (login shell — one time):**

```bash
micromamba activate reseq_bch709
cd ~/scratch/reseq
```

**Then submit each step — each line is one command:**

```bash
# --- Step 1: download FASTQs (no prerequisites) ---
DL_JID=$(sbatch --parsable scripts/01_download.sh)
echo "download     → $DL_JID"

# --- Step 2: reference prep (no prerequisites, runs in parallel with Step 1) ---
REF_JID=$(sbatch --parsable scripts/02_reference.sh)
echo "reference    → $REF_JID"

# --- Step 3: align (waits for BOTH download and reference) ---
ALIGN_JID=$(sbatch --parsable --dependency=afterok:${DL_JID}:${REF_JID} scripts/03_align.sh)
echo "align        → $ALIGN_JID"

# --- Step 4: markdup + BQSR (waits for align) ---
MARK_JID=$(sbatch --parsable --dependency=afterok:${ALIGN_JID} scripts/04_markdup_bqsr.sh)
echo "markdup+bqsr → $MARK_JID"

# --- Step 5a: HaplotypeCaller (waits for markdup) ---
HC_JID=$(sbatch --parsable --dependency=afterok:${MARK_JID} scripts/05a_haplotypecaller.sh)
echo "haplocaller  → $HC_JID"

# --- Step 5b: gather per-sample GVCFs (waits for HaplotypeCaller) ---
GATHER_JID=$(sbatch --parsable --dependency=afterok:${HC_JID} scripts/05b_gather_gvcf.sh)
echo "gather       → $GATHER_JID"

# --- Step 6: joint genotyping (waits for gather) ---
JOINT_JID=$(sbatch --parsable --dependency=afterok:${GATHER_JID} scripts/06_joint_genotype.sh)
echo "joint_geno   → $JOINT_JID"

# --- Step 7: filter + annotate (waits for joint genotyping) ---
FILT_JID=$(sbatch --parsable --dependency=afterok:${JOINT_JID} scripts/07_filter_annotate.sh)
echo "filter+ann   → $FILT_JID"

# Check that everything is queued
squeue -u $USER
# Steps 3-7 should show state PD with reason (Dependency)
```

After pasting the block, `squeue` shows all 8 jobs — some running, most pending. Slurm takes care of the ordering.

> ## Why copy the commands into your terminal, not a script?
> The hands-on walkthrough is literally what `run_all.sh` does — but by typing each line you *see* each job ID appear and can inspect things in between. Once you're comfortable, just run `bash scripts/run_all.sh` next time.
{: .callout}

---

## 9. Monitoring the Pipeline

### Live status

```bash
squeue -u $USER                                        # see all your jobs
squeue -u $USER --start                                # see expected start times for pending jobs
sacct -u $USER --starttime=today --format=JobID,JobName,State,Elapsed,MaxRSS,ExitCode
```

### Watch logs in real time

```bash
# Latest log from any step
tail -f logs/03_align_*.out | head -200

# Specifically task 2 of the align array job:
tail -f logs/03_align_${ALIGN_JID}_2.out
```

### Kill the whole pipeline if something is wrong

```bash
# Cancel all your jobs at once
scancel -u $USER

# Or cancel just the downstream ones (if e.g. alignment is still fine)
scancel ${MARK_JID} ${HC_JID} ${GATHER_JID} ${JOINT_JID} ${FILT_JID}
```

---

## 10. Resource Tuning

After your first successful run, check actual usage with `sacct` and right-size your requests:

```bash
sacct -u $USER --starttime=today \
      --format=JobID,JobName,State,Elapsed,MaxRSS,ReqMem,ReqCPUs
```

| Observed | Action |
|----------|--------|
| `MaxRSS` << `ReqMem` | Lower `--mem` — reserves less, starts sooner |
| `MaxRSS` ≈ `ReqMem` | Leave ~20% headroom to avoid `OUT_OF_MEMORY` |
| `Elapsed` << `ReqTime` | Lower `--time` |
| `TIMEOUT` | Raise `--time`, or check for stuck loops |
| `OUT_OF_MEMORY` | Raise `--mem` (suggest 2× `MaxRSS`) AND the `-Xmx` flag inside the script |

> **GATK / Picard heap tip:** Always set `-Xmx` 2-4 GB **below** your `#SBATCH --mem`. E.g., `--mem=16g` → `-Xmx12g`. Java needs headroom for off-heap stack/GC; if `-Xmx` equals `--mem`, jobs die with `Out of memory`.
{: .callout}

---

## 11. Expected Wall Time

For 2 Arabidopsis samples (~10-15x coverage each):

| Step | CPU cores × array size | Wall time |
|------|-----------------------|-----------|
| Download | 2 × 2 | ~10 min |
| Reference prep | 4 | ~15 min |
| Align + trim | 8 × 2 | ~30 min |
| MarkDup + BQSR | 4 × 2 | ~20 min |
| HaplotypeCaller (scatter) | 4 × 14 | ~15 min |
| GatherVcfs | 2 × 2 | ~2 min |
| Joint genotyping | 4 | ~5 min |
| Filter + annotate | 4 | ~5 min |
| **Total (pipelined)** | | **~90 min end-to-end** |

Compared to **~8 hours** on a laptop for the same 2 samples. For 100 samples, HPC finishes in about the same 90 minutes; the laptop would need ~1 week.

---

## 12. Common HPC-Specific Errors

> **`Batch job submission failed: Invalid account or account/partition combination`**
>
> **Cause:** Typo in `--account` or `--partition`.
> **Fix:** Re-check with `sacctmgr show user $USER format=user,account,partition,qos%-30 withassoc where user=$USER`.
{: .callout}

> **Array task fails with `sed: no input files` or empty SAMPLE variable**
>
> **Cause:** `samples.tsv` is missing, empty, or has trailing blank lines.
> **Fix:**
> ```bash
> cat -A samples.tsv        # show hidden chars, should have one $ per line, no blanks
> wc -l samples.tsv         # line count must match --array range
> ```
{: .callout}

> **`CANCELLED by 0 ExitCode 0:15` in a dependent job**
>
> **Cause:** A parent job in the dependency chain failed, so Slurm auto-cancels all downstream jobs.
> **Fix:** Find the first `FAILED` job with `sacct -u $USER`, fix it, then resubmit **from that step onward** (don't re-run successful steps).
{: .callout}

> **Jobs stuck in `(Dependency)` state after a parent job succeeded**
>
> **Cause:** Parent had a non-zero exit (even a warning) → `afterok` treats it as failed.
> **Fix:** If the "failure" was a harmless stderr warning, use `--dependency=afterany:${JID}` instead of `afterok`. But first check: was the job actually successful?
{: .callout}

> **`Disk quota exceeded` in `$HOME`**
>
> **Cause:** You accidentally wrote BAMs/VCFs to `~` instead of `~/scratch`.
> **Fix:** `du -sh ~/* | sort -h` to find the offender, then move it to scratch.
{: .callout}

> **`samtools sort` runs out of `/tmp` on compute nodes**
>
> **Cause:** Local `/tmp` is small (often <20 GB). Large BAM sorts can exceed it.
> **Fix:** Use Slurm's job-local scratch variable instead:
> ```bash
> samtools sort -T ${TMPDIR:-/tmp}/${SAMPLE}_sort ...
> ```
{: .callout}

---

## 13. What's Next

- **Scaling to 1000s of samples:** switch to `GenomicsDBImport` + VQSR instead of hard filters.
- **Split big BAMs further:** for human WGS, scatter HaplotypeCaller into 25-50 chunks per sample using interval lists.
- **GWAS:** once your cohort VCF is PASS-only and PLINK-converted, proceed to the [GWAS tutorial]({{site.baseurl}}/episodes/GWAS/).
- **Workflow managers:** for production pipelines, consider [Nextflow](https://www.nextflow.io/) or [Snakemake](https://snakemake.readthedocs.io/) — they handle all this dependency bookkeeping automatically.

---

## References

| Resource | Link |
|---------|------|
| Slurm Array Job guide | [slurm.schedmd.com/job_array.html](https://slurm.schedmd.com/job_array.html) |
| Slurm Job Dependencies | [slurm.schedmd.com/sbatch.html#OPT_dependency](https://slurm.schedmd.com/sbatch.html#OPT_dependency) |
| GATK parallelism docs | [gatk.broadinstitute.org/hc/en-us/articles/360035890411](https://gatk.broadinstitute.org/hc/en-us/articles/360035890411) |
| Pronghorn user guide | [UNR Research Computing](https://www.unr.edu/research-computing/hpc) |
| Basic resequencing (laptop) | [Resequencing Tutorial]({{site.baseurl}}/episodes/resequencing_tutorial/) |
| HPC Cluster fundamentals | [HPC Cluster Lesson]({{site.baseurl}}/episodes/HPC_cluster/) |
