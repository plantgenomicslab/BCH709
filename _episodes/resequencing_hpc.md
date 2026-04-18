---
layout: page
title: Resequencing on Pronghorn HPC (Parallelized with Slurm)
published: true
---

> ## Prerequisites
> You should already have:
> 1. Pronghorn account and SSH access — see the [HPC Cluster lesson]({{site.baseurl}}/episodes/HPC_cluster/)
> 2. Micromamba installed on Pronghorn
> 3. Scratch directory at `~/scratch` (symlinked to `/data/gpfs/assoc/bch709-6/<your_netid>`)
> 4. Familiarity with `sbatch`, `squeue`, `sacct`, and `#SBATCH` directives
>
> Run the basic [Resequencing Tutorial]({{site.baseurl}}/episodes/resequencing_tutorial/) **first** on your laptop to understand the workflow logically before parallelizing it.
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
micromamba create -n reseq_bch709 -c bioconda -c conda-forge python=3.11 -y
micromamba activate reseq_bch709

micromamba install -c bioconda -c conda-forge \
    fastqc fastp bwa-mem2 samtools bcftools tabix \
    openjdk=17 picard gatk4 snpeff plink multiqc -y
```

> **libcrypto fix (if samtools complains):**
> ```bash
> ENV_PREFIX=$(micromamba info --envs | awk '/reseq_bch709/ {print $NF}')
> ln -sf ${ENV_PREFIX}/lib/libcrypto.so.3 ${ENV_PREFIX}/lib/libcrypto.so.1.0.0
> ```
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

Make one file listing every sample you want to process. Each column is tab-separated: `sample_id  SRR_accession`.

```bash
cat > ~/scratch/reseq/samples.tsv <<'EOF'
sample1	SRR519585
sample2	SRR519586
EOF
```

> **Adding more samples later** is trivial — just append rows to `samples.tsv`. The array size in every script automatically scales (you only change one number).
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

echo "[task ${SLURM_ARRAY_TASK_ID}] Downloading ${SAMPLE} (${SRR})"

PREFIX=${SRR:0:6}  # first 6 chars for ENA directory structure
wget -q -O raw/${SAMPLE}_R1.fastq.gz \
    ftp://ftp.sra.ebi.ac.uk/vol1/fastq/${PREFIX}/${SRR}/${SRR}_1.fastq.gz
wget -q -O raw/${SAMPLE}_R2.fastq.gz \
    ftp://ftp.sra.ebi.ac.uk/vol1/fastq/${PREFIX}/${SRR}/${SRR}_2.fastq.gz

ls -lh raw/${SAMPLE}_R*.fastq.gz
```

**Submit:**

```bash
cd ~/scratch/reseq
sbatch scripts/01_download.sh
squeue -u $USER   # you should see 2 tasks (download_1, download_2) running in parallel
```

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
#SBATCH --time=02:00:00
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=<YOUR_EMAIL>
#SBATCH -o logs/02_reference_%j.out

set -euo pipefail
source activate reseq_bch709 2>/dev/null || micromamba activate reseq_bch709

cd ~/scratch/reseq

# Download TAIR10
wget -q https://ftp.ensemblgenomes.org/pub/plants/release-60/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz
gunzip -f Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz
mv -f Arabidopsis_thaliana.TAIR10.dna.toplevel.fa reference.fasta

# Build all indices
bwa-mem2 index reference.fasta
samtools faidx reference.fasta
picard CreateSequenceDictionary R=reference.fasta O=reference.dict

# Download and prepare known variants (sites-only to avoid malformed VCF)
wget -q https://1001genomes.org/data/GMI-MPI/releases/v3.1/1001genomes_snp-short-indel_only_ACGTN.vcf.gz
zcat 1001genomes_snp-short-indel_only_ACGTN.vcf.gz \
    | awk 'BEGIN{OFS="\t"} /^##/{print; next} /^#CHROM/{print $1,$2,$3,$4,$5,$6,$7,$8; next} {print $1,$2,$3,$4,$5,$6,$7,$8}' \
    | bgzip > known_sites.vcf.gz
tabix -p vcf known_sites.vcf.gz
rm -f 1001genomes_snp-short-indel_only_ACGTN.vcf.gz

echo "Reference prep done."
```

**Submit:**

```bash
sbatch scripts/02_reference.sh
```

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
source activate reseq_bch709 2>/dev/null || micromamba activate reseq_bch709

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

**Submit with dependency** (wait for reference prep to succeed first):

```bash
# Get the reference job ID from the previous step, or chain them
REF_JID=$(sbatch --parsable scripts/02_reference.sh)
ALIGN_JID=$(sbatch --parsable --dependency=afterok:${REF_JID} scripts/03_align.sh)
echo "Reference: ${REF_JID} | Align: ${ALIGN_JID}"
```

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
source activate reseq_bch709 2>/dev/null || micromamba activate reseq_bch709

cd ~/scratch/reseq

SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" samples.tsv | cut -f1)

# ---- MarkDuplicates ----
picard -Xmx12g MarkDuplicates \
    I=bam/${SAMPLE}.bam \
    O=bam/${SAMPLE}.markdup.bam \
    M=bam/${SAMPLE}.markdup.metrics \
    VALIDATION_STRINGENCY=SILENT
samtools index -@ ${SLURM_CPUS_PER_TASK} bam/${SAMPLE}.markdup.bam

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
source activate reseq_bch709 2>/dev/null || micromamba activate reseq_bch709

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
source activate reseq_bch709 2>/dev/null || micromamba activate reseq_bch709

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
source activate reseq_bch709 2>/dev/null || micromamba activate reseq_bch709

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
source activate reseq_bch709 2>/dev/null || micromamba activate reseq_bch709

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
