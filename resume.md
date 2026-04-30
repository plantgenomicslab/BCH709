# BCH709 — session checkpoint (resume.md)

**Last updated:** 2026-04-29 (post-live-class triage, second pass)
**Branch:** `gh-pages` (deploys to https://plantgenomicslab.github.io/BCH709/)
**Repo / issues:** https://github.com/plantgenomicslab/BCH709

## Live-class triage on 2026-04-28/30 — issues #46–#69

Twenty-four new issues opened and closed across seven waves of live-class and post-audit triage:

| # | File | Fix |
|---|------|-----|
| #46 | `resequencing_hpc.md` step 5b | `GatherVcfs` → `MergeVcfs` (TAIR10 dict orders contigs `1,2,3,4,5,Pt,Mt`; `GatherVcfs` rejected the `Mt`-before-`Pt` array. `MergeVcfs` sorts against the dict.) |
| #47 | `resequencing_hpc.md` env recipe | `snpeff<=5.2` → `snpeff<5.2`. Bioconda repackaged `snpeff=5.2` to require Java 21 (class file 65); we ship `openjdk=17`. Pin now lands on `5.1d`. Existing envs need `micromamba install -n reseq_bch709 -c bioconda 'snpeff<5.2' -y`. |
| #48 | `index.md` | Linked `pipeline_orchestration.md` from the schedule (Week 14 supplement row). |
| #49 | `resequencing_hpc.md` sections 6/7 | Filled `<NNN>` placeholders with clegenbauer's real numbers: joint = 1,219,245 records / 1,005,193 SNPs / 215,941 indels / 15,281 multiallelic / 3,032 multiallelic SNP. PASS = 898,373. snpEff numbers left as placeholders (clegenbauer's snpEff failed on the Java mismatch above). |
| #50 | `HPC_RNA_SEQ.md` | VectorBase Anopheles URLs `Current_Release/VectorBase-54_…` returned 404. Pinned to `release-68/VectorBase-68_…` (verified 200). |
| #51 | `resequencing_hpc.md` sec 12 | Added "Re-running a single failed step (drop `--dependency=`)" callout listing per-script `sbatch` command + `--array=N` override + step-7 snpEff env-patch / step-5b script-recopy / step-2 script-recopy notes. |
| #52 | `HPC_cluster.md` Step 10 | Teaching-level callout explaining when to drop `--dependency=` (concept; pointers to the four HPC pipeline lessons). |
| #53 | `chipseq_hpc.md` sec 12 | Same per-script re-run callout adapted to the 7 chipseq scripts (`01_download`–`07_qc`) with macs3/IDR-specific gotchas. |
| #54 | `HPC_RNA_SEQ.md` | Same callout adapted to BOTH Arabidopsis and Drosophila pipelines (loop-driven, so `samples.txt` editing replaces `--array=N`). MultiQC `afterany` reminder included. |
| #55 | `resequencing_hpc.md` step 8 | Step 7 emits `cohort.snps.filtered.vcf.gz` + `cohort.indels.filtered.vcf.gz` separately; step 8 was reading the (never-produced) merged `cohort.filtered.vcf.gz`. Step 8 now does an idempotent `bcftools concat -a` to build the merged file before `bcftools stats`. |
| #56 | `chipseq_hpc.md` run_all | `07_qc.sh` dispatch lacked `${REF_JID}` and `${ALIGN_JID}` — partial re-runs of just steps 5/6/7 broke on missing `TSS.bed` (from step 2) or fastp JSONs (from step 3). Tightened to `afterok:${MACS_JID}:${BW_JID}:${REF_JID}:${ALIGN_JID}` in 3 dispatch sites. |
| #57 | `HPC_RNA_SEQ.md` ATH | Three fixes: (a) ATH project-setup `mkdir` now includes `logs/` (Slurm was rejecting jobs because `#SBATCH -o logs/...` had no parent dir); (b) added missing `featureCounts.sh` step + wired into `run_all.sh` (ATH was chaining `align→multiqc` directly while multiqc had `--module featureCounts` — silently empty report); (c) removed stale `--module fastqc` from both ATH and Drosophila multiqc (no FastQC step exists in either). |
| #58 | `resequencing_hpc.md` step 7 | Pre-download snpEff `Arabidopsis_thaliana` DB (guarded `[ ! -d ${HOME}/snpeff_data/Arabidopsis_thaliana ] && snpEff download …`). Was relying on snpEff's silent auto-download from compute nodes. Section 8's MultiQC source-file table also corrected — snpEff defaults emit `snpEff_summary.html` / `_genes.txt` to working dir, not `vcf/`. |
| #59 | `resequencing_hpc.md` step 7 | snpEff missing `-csvStats vcf/snpEff_summary.csv`. MultiQC's snpeff module ONLY parses the CSV (not HTML/genes.txt), so the final report silently dropped the variant-impact panel. Sec 8 source-file table updated to point at the CSV. Reported live by a student. |
| #60 | `BLAST.md` | `blastp -query your_protein.fasta` referenced a placeholder file the lesson never created. Replaced with `P01308.fasta` (insulin), already downloaded two lines earlier. |
| #61 | `Linux_Enviroment_and_command_line.md` | FASTA-handling section (line ~863) ran `awk/grep ... mrna.fa` without a preceding `cd ~/bch709_data`. The earlier `cd` was buried inside a collapsed solution block several sections back. Added explicit `cd` + one-line narrative. |
| #62 | `chipseq_tutorial.md` | Sec 9 (`multiBamSummary`) and Sec 15 (IDR) both referenced `chip_rep1.dedup.bam` / `chip_rep2.dedup.bam` and `chip_rep1_peaks.narrowPeak` / `chip_rep2_peaks.narrowPeak` — the laptop tutorial only processes a single ChIP + Input. Marked both as `⚠️ Demonstration only` and added a runnable pseudo-replicate variant for Sec 9 (split `chip.dedup.bam` in halves) so students can still see `multiBamSummary` + `plotCorrelation` work. Sec 15 kept as concept-only with pointer to ENCODE pipeline. |
| #63 | `resequencing_hpc.md` step 8 | Drop `--module gatk` from `08_multiqc.sh`. MultiQC 1.34's `gatk/base_recalibrator` parser hits a pydantic ValidationError ("points.0.name and points.1.name are both None") and a follow-on `rich.panel` AttributeError, killing the whole report. recal.table itself is fine — `unit_${SAMPLE}` ReadGroup is present. Sec 8 source-file table + report module list updated to direct students to `bam/*.recal.table` for manual BQSR inspection. Re-enable when upstream MultiQC bug is fixed. |
| #64 | `HPC_RNA_SEQ.md` Arabidopsis + Drosophila `fastq-dump.sh` | EBI HTTPS drops mid-stream on >1 GB SRR fastq transfers; `curl --retry 3` doesn't retry SSL eof (curl 7.71+ needs `--retry-all-errors`), no `-C -` resume, and `[ -s "${OUT}" ]` skip-check accepted partial downloads as complete. Hardened both fastq-dump blocks: `--retry 5 --retry-all-errors --retry-delay 30 -C -`; dropped the skip-check (curl `-C -` no-ops on complete files, resumes on partials). |
| #65 | `HPC_cluster.md` fastq-dump | Parity with #64. Hardened the EBI fastq curl + dropped `[ -s ]` skip-check. |
| #66 | `chipseq_hpc.md` downloads | Parity with #64 on TWO sites: ENCODE fastq loop in `01_download.sh` and UCSC hg19.fa.gz primary+mirror fallback in `02_reference.sh` (~948 MB). Both `||` mirror legs hardened. The smaller `refGene.txt.gz` (~5 MB) left unchanged. |
| #67 | `resequencing_hpc.md` downloads | Parity with #64 on THREE sites: EBI fastq loop, TAIR10 multi-mirror FASTA (preserved `-k` for self-signed cert and the legitimate `[ -s ]` mirror-fallback control flow), and the 19 GB 1001genomes known-sites VCF (added missing `--retry-all-errors` to existing `-C -`). |
| #68 | `HPC_RNA_SEQ.md` reference downloads | Parity with #64 on FOUR additional non-fastq download sites: Drosophila FlyBase+Ensembl fallback, Mouse GRCm39 (~900 MB), VectorBase Anopheles (already release-68 pinned in #50), and RefSeq plant.1.protein.faa.gz. Both `||` mirror legs hardened where present; `-O` (uppercase) preserved. |
| #69 | `resequencing_hpc.md` env recipe + step 7 | Final consolidation of the snpEff Java-version saga (#47, #58, #59, #63 lineage). Switched env to a single `reseq_bch709` with `openjdk=21` (no `snpeff<5.2` pin). Java is forward-compatible — Java 21 JVM runs Java 17 GATK/Picard jars unchanged, so a single 21 satisfies snpeff (any version) AND GATK 4.6 / Picard 3. The earlier inline patch line buried between `conda install` and `pip install … multiqc` was causing students to fix snpeff and stop, never reaching the multiqc install — moved the migration note into a standalone `.callout` AFTER the full setup block, with a 3-step verify→install→re-verify pattern (since `micromamba install openjdk=21` sometimes silently keeps Java 17 due to soft pins; included an aggressive `remove openjdk` form as fallback). Added `IMPORTANT: don't skip this — step 8 needs multiqc` comment above the multiqc pip line. |

### Student-side action items (not in lesson — verbal/Slack)

- **All resequencing students** must run `micromamba install -n reseq_bch709 -c bioconda 'snpeff<5.2' -y` once. Existing envs are NOT auto-updated by the lesson edit.
- Students who failed at step 2 with `curl --max-time 1800` timeout or "reference.dict already exists" had a stale local copy of `02_reference.sh`. Re-copy from current lesson — the lesson itself is correct.
- Students who failed at step 5b with the `GatherVcfs` Mt/Pt order error must re-copy `05b_gather_gvcf.sh` (now uses `MergeVcfs`).

## Current goal

Make every BCH709 lesson reproducible end-to-end on UNR Pronghorn so a student can complete the homework without instructor help. All fixes follow the lifecycle **open issue → edit → commit → push → close** so each change is traceable on GitHub.

## What's done in this session

### Pipeline-level work

- **MultiQC integrated into all four HPC pipelines** — runs as the final step, aggregates fastp / Picard / GATK / snpEff / bcftools / STAR / featureCounts artifacts into one HTML report per pipeline (#33 reseq, #35 RNA-Seq; chipseq already had it).
- **Shared workspace `~/scratch/rnaseq/`** — HPC_cluster.md and HPC_RNA_SEQ.md both work under that parent dir so the two lessons hand off cleanly without cluttering scratch root (#10, #14).
- **Atomic `known_sites.vcf.gz` build** — step 02 writes via `.tmp` files and verifies before promoting; step 04 has a preflight check that fails loud if the dict/index are missing (#27, live-class report).
- **Idempotent reference build** — `rm -f reference.dict` before `picard CreateSequenceDictionary` so re-running step 02 after a partial failure no longer aborts (#44, live-class report).
- **HPC_RNA_SEQ Drosophila section** expanded to a full student walkthrough — samples.txt, loop-driven trim/align scripts, featureCounts, MultiQC, run_all.sh, hands-on dependency-chain narrative, expected outputs at every step (#43, +553 lines).
- **`Expected output:` blocks** added across all four HPC pipelines so students can self-verify each step (resequencing #25, RNA-seq sections, chipseq sections, hpc_cluster sections).

### New pages

- **`_episodes/pipeline_orchestration.md`** (652 lines) — side-by-side comparison of three orchestrators using the same Arabidopsis 6-sample RNA-Seq pipeline:
  - Slurm `--dependency` chain (this course's HPC tutorials)
  - Snakemake (full Snakefile + per-rule `resources:` / `conda:` + `--executor slurm` runbook + DAG/rerun demos)
  - Nextflow (full `main.nf` + `nextflow.config` + per-process resources + `-resume` / `-with-report` demos + nf-core/rnaseq runbook)
  - Includes: comparison table, decision tree, "common pitfalls when migrating" mapping table.

### Lessons validated end-to-end

Each was walked top-to-bottom by a subagent that ran every executable example and committed fixes via the issue→push lifecycle.

| Lesson | File | Result |
|--------|------|--------|
| HPC Cluster basics | HPC_cluster.md | ✅ #5 hostname escape, #14 shared scratch path |
| HPC Resequencing | resequencing_hpc.md | ✅ #3 sra-tools, #4 snpeff Java, #24 19 GB curl, #27 atomic known_sites, #44 dict idempotent |
| HPC ChIP-Seq | chipseq_hpc.md | ✅ #6 gcc/cython, #16 pysam glibc, #26 macs3 cross-toolchain, #28 env simplification |
| HPC RNA-Seq | HPC_RNA_SEQ.md | ✅ #7,8,10–13,15,17–23 (path/sample/typo fixes), #35 MultiQC, #43 Drosophila expansion |
| Linux env | Linux_Enviroment_and_command_line.md | ✅ #29 grep MacOS, #30 FASTA case, #32 GFF3 outputs |
| Conda / Compile | compile.md | ✅ #37 versions, #38 BWA 0.7.17 GCC 10+ |
| GitHub & server | github_server.md | ✅ #31 xclip clipboard selection |
| Vibe coding | vibe_coding.md | ✅ #34 stale outputs, FASTA header, Gasch column layout |
| Sequence manipulation | seqkit_tutorial.md | ✅ #40 hairpin.fa.gz typos, glob, csvtk numeric sort |
| BLAST | BLAST.md | ✅ #36 ftp:// → https://, -parse_seqids, outfmt, stale stats |
| RNA-Seq tutorial (laptop) | RNA-seq_tutorial.md | ✅ #39 bc fallback, STAR --genomeSAindexNbases 10 |
| ChIP-Seq tutorial (laptop) | chipseq_tutorial.md | ✅ #41 chip.fastq.gz filename, #42 hg38→hg19 |
| Resequencing tutorial (laptop) | resequencing_tutorial.md | ✅ #45 Expected output blocks (no bug fixes — tutorial runs end-to-end as written) |

### All validation complete

Every lesson linked from `index.md` has been walked top-to-bottom by a subagent. Tool versions verified (fastp 1.3.3, bwa-mem2 2.2.1, samtools 1.23.1, gatk 4.6.2.0, snpEff 5.2, PLINK v1.9.0-b.8, etc.). All 45 GitHub issues opened during this session are CLOSED.

## What's still open

(All issues closed as of 2026-04-30. #46–#69 all CLOSED.)

## Live-class hand-off (verbal/Slack)

Existing students whose `reseq_bch709` env was built before #69 must run a 3-step Java upgrade — see the "Java versions and snpEff" callout in section 0 of `resequencing_hpc.md`. The first form (`micromamba install … openjdk=21 snpeff`) sometimes silently leaves Java 17 in place; if `java -version` still shows 17 after install, the aggressive form (`micromamba remove openjdk` then install) is required. Always verify with `micromamba run -n reseq_bch709 java -version` afterwards — don't assume the install worked.

## Outstanding follow-ups (no GitHub issues open — context for the next session)

1. **Resequencing snpEff `<NNN>` placeholders** — sections 7's snpEff `Number_of_variants_processed`, HIGH/MODERATE/LOW/MODIFIER, MISSENSE/NONSENSE/SILENT counts are still `<NNNNNN>`. clegenbauer's snpEff hit the Java 21 mismatch (now fixed in #47), so when a student re-runs step 7 with the patched env, capture the real numbers and fill them in. Steps 4 (ApplyBQSR) and 5 (HaplotypeCaller scatter) Expected output blocks also still hold placeholders.
2. **File-system invariants** — `/data/gpfs/assoc/bch709-6/Course_material` is a symlink to `/data/gpfs/assoc/pgl/Lecture/Course_material`. If that bch709-6 dir is ever reset, recreate the symlink and `test_mrna.fna` (gunzipped from `Athaliana_167_TAIR10.cds.fa.gz`) — several tutorials reference both paths.
3. **Snakemake / Nextflow homework** — `pipeline_orchestration.md` has a "Try it" challenge re-implementing the Arabidopsis pipeline; could become a graded assignment.
4. ~~`index.md` schedule~~ — done in #48 (Week 14 supplement row).
5. ~~VectorBase Anopheles + UCSC EU paths~~ — VectorBase done in #50 (release-68 pin). UCSC EU mirror was a planning note only — no committed lesson reference; chipseq lessons already use the working US mirror `hgdownload.soe.ucsc.edu` (verified 200, ~948 MB).

## Workflow rules in effect for this repo

These are saved in claude memory and enforced by every subagent:

- **Issue → fix → commit → push → close** for every markdown bug. Don't batch unrelated fixes.
- **Embed real outputs** as `Expected output:` blocks where useful for student self-verification. Sanitize netid → `<netid>`, real Slurm job IDs → `<jobid>`, full datetimes → `<date>`.
- **Group output additions** per file into one issue + one commit, separate from bug fixes.
- **Don't fragment** — one issue per logical bug, one commit per issue, push immediately.
- **gh CLI** is authenticated as `wyim-pgl`. Push protocol: `git@github.com:plantgenomicslab/BCH709.git`, branch `gh-pages`.

## Useful shortcuts when picking this back up

- Recent commits: `git log --oneline origin/gh-pages | head -20`
- Open issues: `gh issue list --repo plantgenomicslab/BCH709 --state open`
- Per-validation working dirs: `/data/gpfs/assoc/bch709-6/wyim/tutorial_runs/<lesson>/`
- Reseq sample data already on disk: `/data/gpfs/assoc/bch709-6/wyim/tutorial_runs/resequencing/reseq/`
- Pre-class smoke test (run all 3 HPC pipelines end-to-end): use the validated `tutorial_runs/<lesson>/` dirs as a sanity check.

## Quick environment health checks (one-liner per env)

```bash
samtools --version | head -1
bcftools --version | head -1     # reseq only
multiqc --version
which fastp idr macs3 deeptools  # chipseq specific
which STAR fasterq-dump featureCounts  # RNA-Seq specific
```
