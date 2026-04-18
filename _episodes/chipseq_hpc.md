---
layout: page
title: ChIP-Seq on Pronghorn HPC (Parallelized with Slurm)
published: true
---

> ## Prerequisites
> You should already have:
> 1. Pronghorn account and SSH access — see the [HPC Cluster lesson]({{site.baseurl}}/episodes/HPC_cluster/)
> 2. Micromamba installed on Pronghorn
> 3. Scratch directory at `~/scratch` (symlinked to `/data/gpfs/assoc/bch709-6/<your_netid>`)
> 4. Familiarity with `sbatch`, `squeue`, `sacct`, and `#SBATCH` directives
>
> Run the basic [ChIP-Seq Tutorial]({{site.baseurl}}/episodes/chipseq_tutorial/) **first** on your laptop to understand the workflow logically before parallelizing it.
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
micromamba create -n chipseq_bch709 -c bioconda -c conda-forge python=3.11 -y
micromamba activate chipseq_bch709

# Conda installs for alignment/QC
micromamba install -c bioconda -c conda-forge \
    fastqc fastp minimap2 samtools bedtools tabix \
    openjdk=17 picard homer -y

# deepTools + MACS3 + MultiQC via pip (conda has conflicts)
pip install deeptools macs3 multiqc

# IDR (from GitHub — the pip "idr" is a different project)
pip install "numpy<1.24" git+https://github.com/nboley/idr.git
```

> **libcrypto fix (if samtools complains):**
> ```bash
> ENV_PREFIX=$(micromamba info --envs | awk '/chipseq_bch709/ {print $NF}')
> ln -sf ${ENV_PREFIX}/lib/libcrypto.so.3 ${ENV_PREFIX}/lib/libcrypto.so.1.0.0
> ```
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

```bash
cat > ~/scratch/chipseq/samples.tsv <<'EOF'
sample	accession	type	control
ctcf_rep1	ENCFF001HTO	chip	input_k562
ctcf_rep2	ENCFF001HTP	chip	input_k562
input_k562	ENCFF001HTT	input	-
EOF
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

wget -q -c "${URL}" -O raw/${SAMPLE}.fastq.gz

# Validate the gzip immediately — prevents silent truncation
gunzip -t raw/${SAMPLE}.fastq.gz
ls -lh raw/${SAMPLE}.fastq.gz
echo "${SAMPLE} OK"
```

**Submit:**

```bash
cd ~/scratch/chipseq
sbatch scripts/01_download.sh
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
#SBATCH -o logs/02_reference_%j.out

set -euo pipefail
source activate chipseq_bch709 2>/dev/null || micromamba activate chipseq_bch709

cd ~/scratch/chipseq

# hg19 genome (example data is ENCODE-legacy aligned to hg19)
wget -q http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz
gunzip -f hg19.fa.gz
mv -f hg19.fa reference.fasta

samtools faidx reference.fasta
minimap2 -d reference.mmi reference.fasta

# TSS BED for later heatmaps
wget -q http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/refGene.txt.gz
gunzip -f refGene.txt.gz
awk 'BEGIN{OFS="\t"} {
    if ($4=="+") print $3, $5, $5+1, $2, "0", $4;
    else         print $3, $6-1, $6, $2, "0", $4
}' refGene.txt | sort -k1,1 -k2,2n | uniq > TSS.bed

echo "Reference prep done."
```

Save as `scripts/02_reference.sh` and submit:

```bash
REF_JID=$(sbatch --parsable scripts/02_reference.sh)
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
source activate chipseq_bch709 2>/dev/null || micromamba activate chipseq_bch709

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

**Submit with dependency on reference prep:**

```bash
ALIGN_JID=$(sbatch --parsable --dependency=afterok:${REF_JID} scripts/03_align.sh)
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
source activate chipseq_bch709 2>/dev/null || micromamba activate chipseq_bch709

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
source activate chipseq_bch709 2>/dev/null || micromamba activate chipseq_bch709

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
source activate chipseq_bch709 2>/dev/null || micromamba activate chipseq_bch709

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
#SBATCH -o logs/07_qc_%j.out

set -euo pipefail
source activate chipseq_bch709 2>/dev/null || micromamba activate chipseq_bch709

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

# ---- MultiQC aggregates fastp + flagstat + MACS3 + fingerprint ----
multiqc . -o qc/ -n chipseq_report
```

Submit:

```bash
QC_JID=$(sbatch --parsable --dependency=afterok:${MACS_JID}:${BW_JID} scripts/07_qc.sh)
```

---

## 8. Submit the Entire Pipeline with One Script

```bash
#!/bin/bash
# scripts/run_all.sh
set -euo pipefail
cd ~/scratch/chipseq

N_SAMPLES=$(awk 'NR>1' samples.tsv | wc -l)
N_CHIP=$(awk -F'\t' 'NR>1 && $3=="chip"' samples.tsv | wc -l)

DL_JID=$(sbatch   --parsable --array=1-${N_SAMPLES}                                 scripts/01_download.sh)
REF_JID=$(sbatch  --parsable                                                         scripts/02_reference.sh)
ALIGN_JID=$(sbatch --parsable --array=1-${N_SAMPLES} --dependency=afterok:${DL_JID}:${REF_JID} scripts/03_align.sh)
DEDUP_JID=$(sbatch --parsable --array=1-${N_SAMPLES} --dependency=afterok:${ALIGN_JID}         scripts/04_dedup.sh)
BW_JID=$(sbatch    --parsable --array=1-${N_SAMPLES} --dependency=afterok:${DEDUP_JID}         scripts/05_bigwig.sh)
MACS_JID=$(sbatch  --parsable --array=1-${N_CHIP}    --dependency=afterok:${DEDUP_JID}         scripts/06_macs3.sh)
QC_JID=$(sbatch    --parsable --dependency=afterok:${MACS_JID}:${BW_JID}                       scripts/07_qc.sh)

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
> **Fix:** find the failing task with `sacct -u $USER --starttime=today | grep FAILED`, fix it, resubmit **only that step** and update downstream `--dependency` to the new job ID.
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
