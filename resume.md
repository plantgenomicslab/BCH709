# BCH709 HPC pipelines — checkpoint (resume.md)

**Last updated:** 2026-04-28 (work session ongoing)
**Branch:** `gh-pages` (deploys to https://bch709.plantgenomicslab.org/)

## Current goal

Make the four HPC lessons (HPC_cluster, resequencing_hpc, chipseq_hpc,
HPC_RNA_SEQ) bulletproof so students can run the homework end-to-end
without instructor intervention.

## What's done in this session

1. **samples.tsv tab safety** — replaced heredoc with `printf '\t'` +
   `awk -F'\t' '{print NF}'` verification step. Browser copy-paste was
   silently converting tabs to spaces.
2. **NCBI SRA download for resequencing** — replaced `wget` from
   `ftp.sra.ebi.ac.uk` with `prefetch` + `fasterq-dump`. Added
   `sra-tools>=3.0` to the reseq env.
3. **conda → micromamba** — every `conda` invocation in HPC_RNA_SEQ.md
   converted to micromamba (DEG_bch709, blast, venn envs + tree).
4. **IDR install** — split the chipseq pip install into Step A
   (numpy/pyarrow/multiqc/etc.) + Step B (`pip install
   --no-build-isolation git+...idr.git`) so IDR's setup.py can
   `import numpy` from the active env.
5. **multiqc / tiktoken Rust build failure** — added `pip install
   --upgrade pip` first, `tiktoken<0.8` pin, and `--prefer-binary`
   in all 4 env recipes.
6. **TAIR10 Ensembl Plants flaky download** — multi-mirror fallback
   chain Ensembl → EBI → NCBI RefSeq with auto-rename of NC_*
   chromosome names to TAIR-style 1,2,…,Mt,Pt.
7. **All other download URLs hardened** — switched `wget -q` to
   `curl -fsSL --retry 3 --max-time NNN`, added HTTPS, added mirror
   fallbacks where appropriate (FlyBase → Ensembl, UCSC primary +
   EU mirror, NCBI mouse, VectorBase Anopheles, 1001genomes VCF).

## Broken URLs detected by HEAD-test (2026-04-28)

| Host | Status | What we did |
|---|---|---|
| `ftp.ensemblgenomes.org` | ❌ 000 from this network (entire host) | Reseq already has EBI mirror + NCBI fallback ✅ |
| `ftp.flybase.net` | ❌ 000 (entire host) | Drosophila section already has Ensembl r104 fallback ✅ |
| `vectorbase.org/.../AstephensiSDA-500/...` | ❌ 404 (path reorganized) | **TODO** — find current Anopheles stephensi release path on VectorBase |
| `hgdownload-euro.soe.ucsc.edu/.../hg19.fa.gz` | ❌ 404 | **TODO** — UCSC EU mirror uses a different layout for hg19; remove or replace |
| `ftp.ensemblgenomes.org/pub/plants/release-*/...` (any release) | ❌ 000 | (covered by EBI mirror — keep both URLs in case service returns) |

The remaining URLs all return HTTP 200 (or 206 for ENCODE which
blocks HEAD but allows GET).

## TODO (next session)

1. **VectorBase Anopheles** — find correct current path. The
   `Current_Release/AstephensiSDA-500/` folder has been removed.
   Likely candidates: VectorBase v60+ uses `ASTeI2/` or new genome
   release codes. Need to crawl `https://vectorbase.org/common/downloads/Current_Release/`
   to find the new directory.
2. **UCSC EU mirror** — verify the actual path on
   `hgdownload-euro.soe.ucsc.edu`. Either fix the path or drop the
   EU mirror fallback (the US primary is already reliable).
3. **Pre-class smoke test** — actually run all three pipelines
   end-to-end on Pronghorn before assigning as homework. Dry-runs
   only validate package resolution, not runtime correctness.
4. **Trinity 2.15.2 helper Perl modules** — confirm
   `align_and_estimate_abundance.pl` works (it needs `perl-dbi`
   + `perl-dbd-sqlite`, both already added).

## Open GitHub issue

A tracking issue for the broken URLs has been opened on
`plantgenomicslab/BCH709` so students/TAs can flag new breakage. See
issue title: "HPC homework — broken external download URLs (audit
2026-04-28)".

## Recent commit history

```
a47542e fix Arabidopsis TAIR10 download — multi-mirror fallback
d4b5329 fix multiqc tiktoken Rust build (pip too old to find wheel)
20734f3 fix IDR install: --no-build-isolation
5b04fbc samples.tsv tab safety, NCBI SRA download, conda → micromamba
e57ea32 fix env recipes (channel order, modern pins, dry-run verified)
784b07e HPC_RNA_SEQ.md homework cleanup
0177c05 fix missing mail directives in chipseq_hpc; reschedule RNA-Seq
f821bcd unify HPC_RNA_SEQ paths/account with HPC_cluster
```

(plus uncommitted work on this session's URL hardening — to be committed
together once VectorBase + UCSC EU paths are resolved or dropped).

## Quick environment health checks (one-liner per env)

```bash
# After install, students should run these to verify:
samtools --version | head -1
bcftools --version | head -1     # reseq only
multiqc --version
which fastp idr macs3 deeptools  # chipseq specific
which Trinity STAR fasterq-dump  # RNA-Seq specific
```
