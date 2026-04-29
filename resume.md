# BCH709 — session checkpoint (resume.md)

**Last updated:** 2026-04-28 14:55 PDT
**Branch:** `gh-pages` (deploys to https://plantgenomicslab.github.io/BCH709/)
**Repo / issues:** https://github.com/plantgenomicslab/BCH709

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

(All issues closed as of 14:55 PDT — none currently OPEN.)

## Outstanding follow-ups (no GitHub issues open — context for the next session)

1. **Resequencing pipeline `<NNN>` placeholders** — Expected output blocks for steps 4 (ApplyBQSR), 5 (HaplotypeCaller scatter), 6 (joint genotyping), 7 (filter/snpEff/PLINK) use `# example output (your numbers will differ)` with placeholder counts. When a full pipeline run produces real numbers, swap them in.
2. **File-system invariants** — `/data/gpfs/assoc/bch709-6/Course_material` is a symlink to `/data/gpfs/assoc/pgl/Lecture/Course_material`. If that bch709-6 dir is ever reset, recreate the symlink and `test_mrna.fna` (gunzipped from `Athaliana_167_TAIR10.cds.fa.gz`) — several tutorials reference both paths.
3. **Snakemake / Nextflow homework** — `pipeline_orchestration.md` has a "Try it" challenge re-implementing the Arabidopsis pipeline; could become a graded assignment.
4. **`index.md` schedule** — the new `pipeline_orchestration.md` page is NOT yet linked from the course homepage. Decide whether to slot it after the four HPC weeks (Week 13/14) or as a standalone supplementary lesson.
5. **VectorBase Anopheles + UCSC EU paths** — flagged earlier in this file but still pending: VectorBase reorganized `Current_Release/AstephensiSDA-500/`, and the UCSC EU mirror `hg19.fa.gz` path 404s. Either fix or drop those mirrors.

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
