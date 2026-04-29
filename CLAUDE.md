# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project overview

BCH709 is a bioinformatics course website (Jekyll + Carpentries lesson template) deployed at https://plantgenomicslab.github.io/BCH709/. **The deployed branch is `gh-pages`, not `main`** — pushes to `gh-pages` go live on the public site within a minute or two, so students see edits immediately. Treat every tutorial change as production-impacting.

## Build / preview / lint

```bash
bundle install        # one-time Ruby deps
make serve            # local Jekyll preview at http://localhost:4000
make site             # build only
make lesson-check     # validate lesson Markdown (Carpentries linter)
make lesson-md        # convert _episodes_rmd/*.Rmd → _episodes/*.md
make clean            # remove build artifacts
```

`./commit.sh` pulls, stages all, commits with message `"BCH709"`, and pushes — convenient for quick edits, but do not use it when the change deserves a real commit message (most tutorial fixes do).

## Lesson taxonomy

The 67 files in `_episodes/` are not interchangeable — they fall into three groups, each with different conventions:

- **HPC pipelines** (`HPC_cluster.md`, `resequencing_hpc.md`, `chipseq_hpc.md`, `HPC_RNA_SEQ.md`) — Slurm-driven tutorials that run on UNR Pronghorn. Students execute every code block; mistakes here cost classroom time. These share a `~/scratch/<topic>/` workspace convention so the lessons hand off to each other.
- **Laptop tutorials** (`Linux_Enviroment_and_command_line.md`, `compile.md`, `seqkit_tutorial.md`, `BLAST.md`, `RNA-seq_tutorial.md`, `resequencing_tutorial.md`, `chipseq_tutorial.md`, `github_server.md`, `vibe_coding.md`) — work on a single machine; same tools as the HPC pipelines but smaller datasets.
- **Conceptual / supplementary** (`pipeline_orchestration.md`, year-archived `*_2019.md`, etc.) — comparison pages, theory.

When picking which file to edit for a reported student bug, first identify the group; the same tool (e.g. `bwa-mem2`) appears in both an HPC and a laptop lesson and the fix usually applies to only one.

## HPC tutorial conventions

Every Slurm script block in the HPC pipeline lessons follows this header (do not deviate):

```bash
#SBATCH --account=cpu-s5-bch709-6
#SBATCH --partition=cpu-core-0
#SBATCH --time=...
#SBATCH --mail-user=<YOUR_EMAIL>          # placeholder — students replace
#SBATCH -o logs/<step>_%j.out
```

`<YOUR_EMAIL>` is the only placeholder students replace. Account / partition are stable for the semester — don't substitute them.

**Env manager:** the HPC tutorials use **micromamba**, not conda. Each pipeline has its own env, built once per student: `reseq_bch709` (resequencing), `chipseq_bch709` (chipseq), `RNASEQ_bch709` (RNA-Seq), `DEG_bch709` (DESeq2/EdgeR). Activate in the *login shell* before `sbatch` so jobs inherit the right `PATH` (`sbatch --export=ALL` is the default).

**Shared scratch layout** — HPC_cluster.md and HPC_RNA_SEQ.md both work under `~/scratch/rnaseq/` (raw_data/, trim/, star_index/, bam/, qc/) so the cluster basics lesson hands off cleanly to the RNA-Seq pipeline. Don't reintroduce per-lesson top-level scratch dirs.

**MultiQC is the final step** of every HPC pipeline — `08_multiqc.sh` (resequencing), `multiqc.sh` (RNA-Seq), the post-processing step (chipseq). It runs with `--dependency=afterany:` (not `afterok`) so the report still renders if a prior step partially fails.

**`Expected output:` blocks** are the project's self-verification convention. After a code block whose output is informative, add:

````markdown
```bash
some-command
```

Expected output:

```
sanitized output excerpt
```
````

Sanitize captured outputs: real netid → `<netid>`, real Slurm job IDs → `<jobid>`, full datetimes → `<date>`. Truncate >25 lines (first 10–15 + `# ...` + last 5).

**Idempotency** — every step's script should be safe to re-run. `bwa-mem2 index` and `samtools faidx` overwrite by default; Picard `CreateSequenceDictionary` does NOT — `rm -f reference.dict` precedes it. The `known_sites.vcf.gz` build writes via `.tmp` and verifies before promoting. Preserve these patterns when editing.

## Live-class fix workflow

Students run these tutorials in real time. When a bug is reported (or you find one validating):

1. `gh issue create --repo plantgenomicslab/BCH709 --title "[<file>.md] <one-line problem>" --body "..."` → record issue number.
2. Fix in `_episodes/<file>.md`.
3. `git add` + `git commit -m "<summary> (closes #<N>)"`.
4. `git push origin gh-pages` — verify with `git log --oneline -1 origin/gh-pages`.
5. The `closes #N` keyword auto-closes the issue.

Don't batch unrelated fixes. One issue per logical bug, one commit per issue, push immediately. Group `Expected output:` additions for one file into a single dedicated issue + commit.

## Cross-lesson coupling

- `HPC_cluster.md` "Download and preprocess public RNA-Seq data" produces `~/scratch/rnaseq/raw_data/` and `~/scratch/rnaseq/trim/` — `HPC_RNA_SEQ.md`'s "🔁 Already ran fastq-dump + trim in HPC_cluster?" callout consumes them. Editing path or filename suffix on either side requires updating the other.
- `index.md` schedule references most `_episodes/*.md` by `<basename>/index.html` (Jekyll permalink). New lessons are NOT auto-linked — add a row to `index.md` if it should appear in the schedule table.
- Course material (FASTAs, scripts) lives at `/data/gpfs/assoc/pgl/Lecture/Course_material/` on Pronghorn, with a symlink at `/data/gpfs/assoc/bch709-6/Course_material`. Several tutorials reference the latter path; if the bch709-6 dir is reset, recreate the symlink.

## Validation working dirs

When testing a tutorial change end-to-end on Pronghorn, work under `/data/gpfs/assoc/bch709-6/wyim/tutorial_runs/<lesson>/` (off the gh-pages tree). Each subdir corresponds to a lesson and may already hold artifacts (`reference.fasta`, BAMs, MultiQC reports, etc.) from previous validation runs that can be reused instead of re-downloaded.

## Repository state

`resume.md` at the repo root holds the current end-of-session checkpoint — open issues, lessons validated, what's still in flight. Read it first when picking up cross-session work; update it (and commit) before stopping.
