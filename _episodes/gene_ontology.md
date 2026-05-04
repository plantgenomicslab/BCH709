---
layout: page
title: Gene Ontology & Hypergeometric Enrichment
published: true
---

> ## Paper Reading
> Recommended before class:
> - [Ashburner et al. (2000) Gene Ontology: tool for the unification of biology. *Nature Genetics* 25:25–29](https://www.nature.com/articles/ng0500_25)
> - [Subramanian et al. (2005) Gene set enrichment analysis. *PNAS* 102:15545](https://www.pnas.org/doi/10.1073/pnas.0506580102)
{: .prereq}

---

## 1. What is Gene Ontology (GO)?

The Gene Ontology project provides a **controlled vocabulary** of terms that describe gene-product attributes across all organisms. Every GO term is a node in a directed acyclic graph (DAG), where edges represent relationships such as `is_a`, `part_of`, and `regulates`.

GO is split into three orthogonal **namespaces**:

| Namespace                | Abbreviation | What it captures                                              |
|--------------------------|--------------|---------------------------------------------------------------|
| Biological Process       | BP           | The wider biological objective (e.g. *DNA repair*)            |
| Molecular Function       | MF           | The molecular activity (e.g. *ATP binding*)                   |
| Cellular Component       | CC           | Where in the cell the gene product acts (e.g. *mitochondrion*)|

A single gene can be annotated to many terms across all three namespaces. Annotations carry **evidence codes** (IEA, IDA, ISS, etc.) that record how the assignment was made — IEA is computational/electronic, while IDA is from a direct experiment. **Always check evidence codes before downstream analysis** — many enrichment surprises are driven by IEA-only annotations.

> ## DAG, not tree
> GO is *not* a tree. A child term can have multiple parents, and a gene annotated to a child is *implicitly* annotated to **all** ancestors via the **true-path rule**. This matters for enrichment: when you count how many of your DEGs are in "DNA repair", you must include genes annotated to its descendants too.
{: .callout}

---

## 2. The enrichment question

You have a list of "interesting" genes (e.g. 312 DEGs at padj < 0.05) and a much larger background of all expressed genes (~14,000). The question:

> Is GO term **X** represented in my DEG list **more than expected by chance**, given how often it appears in the background?

This is **over-representation analysis (ORA)**. The natural statistical test is the **hypergeometric test**, equivalent to Fisher's exact test on a 2×2 contingency table.

---

## 3. The hypergeometric model

Imagine the background genes as a bag of marbles:
- $N$ = total genes in background
- $K$ = number of background genes annotated to GO term X
- $n$ = number of DEGs (your draw, without replacement)
- $k$ = number of DEGs annotated to GO term X

The probability of drawing **exactly** $k$ "term-X" genes is:

$$
P(X = k) = \frac{\binom{K}{k}\binom{N-K}{n-k}}{\binom{N}{n}}
$$

The **enrichment p-value** is the probability of seeing **k or more** by chance:

$$
p = \sum_{i=k}^{\min(K, n)} \frac{\binom{K}{i}\binom{N-K}{n-i}}{\binom{N}{n}}
$$

This is a **one-sided** test (over-representation only).

### 3.1 Worked example

| Quantity                        | Value  |
|---------------------------------|--------|
| Background genes ($N$)          | 14,000 |
| Background annotated to "DNA repair" ($K$) | 250 |
| DEGs ($n$)                      | 312    |
| DEGs annotated to "DNA repair" ($k$)       | 18  |

Expected count under the null: $n \cdot K / N = 312 \cdot 250 / 14000 \approx 5.57$. Observed: 18 — over **3× the expected**.

In R:

```r
# phyper(k-1, K, N-K, n, lower.tail = FALSE)
phyper(18 - 1, 250, 14000 - 250, 312, lower.tail = FALSE)
# [1] 1.4e-06
```

So $p \approx 1.4 \times 10^{-6}$ for "DNA repair" enrichment.

> ## Why `k - 1`?
> `phyper(q, ...)` returns $P(X > q)$ when `lower.tail = FALSE`. To get $P(X \geq k)$ you must pass $k - 1$. Off-by-one here is one of the most common bugs in homemade GO tools.
{: .callout}

---

## 4. The background trap

The background set $N$ is **not** "all annotated genes in the genome" — it is the set of genes that **could have been called as DEGs** in your experiment. For RNA-Seq, this means genes expressed at detectable levels (e.g. TPM > 1 in at least one sample).

| Background choice                     | What happens                                              |
|---------------------------------------|-----------------------------------------------------------|
| Whole genome (e.g. 27,000 genes)      | Inflates significance — silent genes pull down term-X frequency, falsely making your DEGs look enriched |
| Expressed genes only (e.g. 14,000)    | **Correct** — matches the universe of "could have been picked" |
| Transcripts on the array/panel only   | Correct for targeted assays                               |

The infamous case: ribosome biogenesis is highly expressed in proliferating tissue. If half your "background" is silent in that tissue, ribosome-related GO terms will dominate every enrichment, regardless of DEG biology. **Always set `universe = expressed_genes` in `clusterProfiler::enrichGO`.**

```r
library(clusterProfiler)
library(org.At.tair.db)

deg_genes      <- read.table("deg_padj_005.txt")$V1
expressed_set  <- read.table("expressed_tpm_gt1.txt")$V1   # ← universe

ego <- enrichGO(
  gene          = deg_genes,
  universe      = expressed_set,                            # ← critical
  OrgDb         = org.At.tair.db,
  keyType       = "TAIR",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.10,
  readable      = TRUE
)
```

---

## 5. Multiple testing correction

GO BP alone has thousands of terms. Testing 5,000 terms at α = 0.05 yields ~250 false positives by chance. Two standard corrections:

| Method      | Controls                | Adjusted threshold (5,000 tests, α = 0.05) | Behavior          |
|-------------|-------------------------|--------------------------------------------|-------------------|
| Bonferroni  | FWER (any false positive) | $p < 0.05 / 5000 = 10^{-5}$              | Very conservative |
| Benjamini-Hochberg (BH) | FDR (expected false-positive proportion) | data-driven; ~5% of called hits expected to be FP | Standard for ORA |

In practice GO enrichment uses **BH-adjusted q-values < 0.05** — Bonferroni is too strict and erases real signal when terms are correlated (parents and children share genes).

### 5.1 Quick mental check

- Raw $p < 0.001$ → 5,000 × 0.001 = **5 expected false positives**.
- If you find 25 terms at raw $p < 0.001$, FDR ≈ 5/25 = **20%** — too noisy.
- Tighten to raw $p < 10^{-4}$ → 0.5 expected FP. If 10 terms pass, FDR ≈ 5%.

---

## 6. Other enrichment paradigms (when ORA isn't enough)

ORA throws away the *ranking* of genes — every DEG is treated equally above the cutoff. **GSEA** (Gene Set Enrichment Analysis) uses the full ranked list (e.g. by log₂ fold-change or Wald statistic), no threshold needed:

```r
geneList <- deg_results$log2FoldChange
names(geneList) <- deg_results$gene_id
geneList <- sort(geneList, decreasing = TRUE)

gsea <- gseGO(geneList = geneList,
              OrgDb    = org.At.tair.db,
              keyType  = "TAIR",
              ont      = "BP",
              pAdjustMethod = "BH")
```

Use ORA when you have a clean cutoff and few DEGs; use GSEA when the signal is broad and weak (subtle perturbations across many genes).

---

## 7. Hands-on assignment

1. Run DESeq2 on the toy A vs B counts in `~/scratch/rnaseq/counts/`.
2. Extract DEGs at `padj < 0.05`.
3. Run `enrichGO` **twice**: once with the whole-genome background and once with `universe = expressed_genes`. Compare top 10 terms — note shifts in q-value.
4. Visualize with `dotplot(ego)` and `cnetplot(ego)`.
5. Report: which background did the previous student probably use if their top hit was "ribosome biogenesis, q = 1e-30"?

---

## 8. Common pitfalls — checklist

- [ ] Background = expressed genes, not whole genome.
- [ ] Evidence codes filtered (drop IEA-only if you want experimental support).
- [ ] BH (q-value), not raw p, for cutoffs.
- [ ] `phyper(k - 1, ...)` — off-by-one.
- [ ] DAG semantics: redundant parent/child terms in the top hits often inflate apparent biology — collapse with `simplify(ego, cutoff = 0.7)`.
- [ ] Report **gene ratio** (k/n) and **bg ratio** (K/N) alongside p — large fold differences with tiny gene counts are noisy.

> ## Reading list
> - Huang et al. (2009) *Bioinformatics enrichment tools: paths toward the comprehensive functional analysis of large gene lists.* Nucleic Acids Res. 37:1.
> - Khatri et al. (2012) *Ten years of pathway analysis: current approaches and outstanding challenges.* PLoS Comput Biol 8:e1002375.
> - Wijesooriya et al. (2022) *Urgent need for consistent standards in functional enrichment analysis.* PLoS Comput Biol 18:e1009935.
{: .callout}
