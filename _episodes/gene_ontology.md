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

### 1.1 Relationship types in detail

Edges in the GO DAG are **typed** — and the type controls whether annotations propagate upward via the true-path rule. The five most common relations:

| Relation                | Reads as                        | Propagates annotations? | Example                                                              |
|-------------------------|---------------------------------|-------------------------|----------------------------------------------------------------------|
| `is_a`                  | "A is a kind of B"              | Yes                     | *mitotic cell cycle* `is_a` *cell cycle*                             |
| `part_of`               | "A is always part of B"         | Yes                     | *nucleolus* `part_of` *nucleus*                                      |
| `regulates`             | "A regulates B"                 | **No** (by default)     | *regulation of DNA repair* `regulates` *DNA repair*                  |
| `positively_regulates`  | "A activates B"                 | **No**                  | *positive regulation of apoptosis* `positively_regulates` *apoptosis*|
| `negatively_regulates`  | "A inhibits B"                  | **No**                  | *negative regulation of translation* `negatively_regulates` *translation*|
| `has_part`              | "B always has A as a part"      | Yes (downward)          | *ribosome* `has_part` *small ribosomal subunit*                      |
| `occurs_in`             | "process A happens in CC B"     | No                      | *protein folding* `occurs_in` *endoplasmic reticulum*                |

The crucial distinction is **regulators are not the regulated process**. If gene *X* is annotated to *positive regulation of DNA repair*, it is **not** counted as a "DNA repair" gene by default in enrichment — it's a regulator of DNA repair, not a repair enzyme. Most ORA implementations honor this: `clusterProfiler::enrichGO()` only follows `is_a` and `part_of` for propagation. If you do want regulators counted, you must explicitly include the *regulation of …* terms in your gene-set definition.

> ## Why this matters in practice
> A common student observation: "I have 10 known *DNA-repair* genes in my DEGs but the GO term *DNA repair* is not enriched." Often the genes are annotated to *regulation of DNA repair* (a sibling term), not *DNA repair* itself. The DAG edge between them is `regulates`, which doesn't propagate — so they do not count toward the parent term.
{: .callout}

### 1.2 Where do GO annotations come from?

A gene-to-term annotation is a curated statement. The pipeline:

1. **Gene Ontology Consortium (GOC)** — the umbrella project that defines the controlled vocabulary itself (the OBO file).
2. **Model-organism databases (MODs)** add experimental annotations for their species:
   - **TAIR** — *Arabidopsis thaliana*
   - **MGI** — mouse
   - **SGD** — *Saccharomyces cerevisiae*
   - **WormBase**, **FlyBase**, **ZFIN**, **RGD**, **PomBase**, **dictyBase**
3. **UniProt-GOA** centralizes annotations across species and adds two large automated streams:
   - **InterPro2GO** — if a protein has an InterPro domain (e.g. *Pfam:PF00069 protein kinase*), assign the linked GO term automatically (evidence code IEA).
   - **UniProtKB-KW2GO** — UniProt keywords (e.g. "Kinase") map to GO terms.
   - **Ensembl Compara** projection — propagate experimental annotations from a well-studied ortholog (e.g. *A. thaliana* TAIR-curated terms transferred to *Brassica rapa*).
4. **Phylogenetic annotation (PAINT / IBA)** — curators look at a gene family tree and assert that an ancestral function existed at a node, then propagate to all descendants with evidence code IBA.

Because each MOD curates differently, **the same gene can have a different GO annotation set in TAIR vs UniProt vs Ensembl Plants**. For *A. thaliana* tutorials in this course we use `org.At.tair.db`, which mirrors TAIR + UniProt-GOA. If a colleague reports different numbers for the same gene list, ask which database they queried.

### 1.3 Evidence codes — the full table

GO evidence codes record *how* an annotation was made. They are critical when filtering for high-confidence enrichment.

| Group                | Code | Name                                                | Trust    |
|----------------------|------|-----------------------------------------------------|----------|
| **Experimental**     | EXP  | Inferred from Experiment                            | High     |
|                      | IDA  | Inferred from Direct Assay                          | High     |
|                      | IPI  | Inferred from Physical Interaction                  | High     |
|                      | IMP  | Inferred from Mutant Phenotype                      | High     |
|                      | IGI  | Inferred from Genetic Interaction                   | High     |
|                      | IEP  | Inferred from Expression Pattern                    | Moderate |
| **Phylogenetic**     | IBA  | Inferred from Biological aspect of Ancestor (PAINT) | Moderate |
|                      | IBD  | Inferred from Biological aspect of Descendant       | Moderate |
|                      | IKR  | Inferred from Key Residues                          | Moderate |
| **Computational**    | ISS  | Inferred from Sequence/Structural Similarity        | Moderate |
|                      | ISO  | Inferred from Sequence Orthology                    | Moderate |
|                      | ISA  | Inferred from Sequence Alignment                    | Moderate |
|                      | ISM  | Inferred from Sequence Model                        | Moderate |
|                      | RCA  | Inferred from Reviewed Computational Analysis       | Moderate |
|                      | IEA  | Inferred from Electronic Annotation                 | **Low**  |
| **Author / curator** | TAS  | Traceable Author Statement                          | Moderate |
|                      | NAS  | Non-traceable Author Statement                      | Low      |
|                      | IC   | Inferred by Curator                                 | Moderate |
| **No data**          | ND   | No biological Data available                        | None     |

The vast majority of annotations in any organism — often **>80%** — are IEA. Filtering them out leaves a much smaller, high-confidence set.

```r
# Pull only experimentally-supported BP annotations for a gene list
library(org.At.tair.db)

ann <- AnnotationDbi::select(
  org.At.tair.db,
  keys     = deg_genes,
  columns  = c("GO", "EVIDENCE", "ONTOLOGY"),
  keytype  = "TAIR"
)

# Keep only "high-trust" experimental codes
exp_codes <- c("EXP", "IDA", "IPI", "IMP", "IGI", "IEP")
ann_exp   <- subset(ann, ONTOLOGY == "BP" & EVIDENCE %in% exp_codes)

# How much did filtering shrink the annotation set?
nrow(ann)        # all evidence
nrow(ann_exp)    # experimental only
```

> ## When to drop IEA
> - **Drop IEA** when you want a confirmatory, publication-ready statement ("genes are *experimentally* known to do X").
> - **Keep IEA** when you are exploring a poorly-studied genome (most plant species). Removing IEA can wipe out 90% of annotations and leave nothing to test.
{: .callout}

---

## 2. The enrichment question

You have a list of "interesting" genes (e.g. 312 DEGs at padj < 0.05) and a much larger background of all expressed genes (~14,000). The question:

> Is GO term **X** represented in my DEG list **more than expected by chance**, given how often it appears in the background?

This is **Over-Representation Analysis (ORA)**. The natural statistical test is the **hypergeometric test**, equivalent to Fisher's exact test on a 2×2 contingency table.

---

## 3. The hypergeometric model

Imagine the background genes as a bag of marbles:
- ***N*** = total genes in background
- ***K*** = number of background genes annotated to GO term X
- ***n*** = number of DEGs (your draw, without replacement)
- ***k*** = number of DEGs annotated to GO term X

> ## Combination notation (`C(n, k)`)
> `C(n, k)` — read "*n* choose *k*", also written `ⁿCₖ` or `(n k)` — is the number of ways to choose *k* items from *n* **without regard to order**:
>
> ```
>            n!
> C(n, k) = ───────────       (read: n-factorial divided by k! · (n−k)!)
>           k! · (n − k)!
> ```
>
> Example: `C(5, 2) = 5! / (2! · 3!) = 120 / (2 · 6) = 10`. The ten 2-element subsets of {1,2,3,4,5} are: {1,2}, {1,3}, {1,4}, {1,5}, {2,3}, {2,4}, {2,5}, {3,4}, {3,5}, {4,5}.
>
> In R: `choose(n, k)` (e.g. `choose(5, 2)` returns `10`).
{: .callout}

The probability of drawing **exactly** *k* "term-X" genes is:

```
              C(K, k) · C(N − K, n − k)
P(X = k) = ───────────────────────────────
                       C(N, n)
```

where `C(a, b)` denotes the binomial coefficient *a* choose *b*.

The **enrichment p-value** is the probability of seeing ***k* or more** by chance:

```
                    min(K, n)
                       ___    C(K, i) · C(N − K, n − i)
p  =  P(X ≥ k)  =      \      ─────────────────────────
                       /__              C(N, n)
                      i = k
```

This is a **one-sided** test (over-representation only).

### 3.1 Worked example

| Quantity                                    | Value  |
|---------------------------------------------|--------|
| Background genes (*N*)                      | 14,000 |
| Background annotated to "DNA repair" (*K*)  |    250 |
| DEGs (*n*)                                  |    312 |
| DEGs annotated to "DNA repair" (*k*)        |     18 |

Expected count under the null:
`E[k] = n · K / N = 312 × 250 / 14000 ≈ 5.57`. Observed: 18 — over **3× the expected**.

In R:

```r
# phyper(k-1, K, N-K, n, lower.tail = FALSE)
phyper(18 - 1, 250, 14000 - 250, 312, lower.tail = FALSE)
# [1] 1.29e-05
```

So `p ≈ 1.3 × 10⁻⁵` for "DNA repair" enrichment.

> ## Why `k - 1`?
> `phyper(q, ...)` returns `P(X > q)` when `lower.tail = FALSE`. To get `P(X ≥ k)` you must pass `q = k - 1`. Off-by-one here is one of the most common bugs in homemade GO tools.
{: .callout}

### 3.2 Computing the p-value — step by step

The hypergeometric distribution `Hypergeometric(N, K, n)` has the following standard properties (cf. [Wikipedia: Hypergeometric distribution](https://en.wikipedia.org/wiki/Hypergeometric_distribution)):

**Mean (expected count):**

```
E[X] = n · K / N
```

For our example: `E[X] = 312 × 250 / 14000 ≈ 5.57`.

**Variance:**

```
              n · K · (N − K) · (N − n)
Var(X) = ─────────────────────────────────
                  N² · (N − 1)
```

Plugging in:

```
Var(X) = 312 × 250 × 13,750 × 13,688
         ─────────────────────────────  ≈ 5.35
              14,000² × 13,999

SD(X)  = √5.35 ≈ 2.31
```

So the null distribution has mean ≈ 5.57 and SD ≈ 2.31. The observed 18 sits about **(18 − 5.57) / 2.31 ≈ 5.4 SD** above the mean — confidently in the upper tail.

> ## ⚠️ Why a *z*-test would mislead you here
> A naïve *p*-value from `pnorm(5.37, lower.tail = FALSE)` ≈ **3.9 × 10⁻⁸**. The true hypergeometric p ≈ **1.3 × 10⁻⁵** — about **330× larger**. When *E[X]* is small (here ≈ 5.6) the discrete count distribution has much heavier upper tails than the normal approximation, so z-tests *underestimate* the p-value and *overstate* significance. Always use `phyper`/Fisher exact for ORA.
{: .callout}

**The 2×2 contingency table (Fisher exact view).** ORA is equivalent to a one-sided Fisher exact test on this table:

|                  | In GO term     | Not in GO term            | Row total      |
|------------------|---------------:|--------------------------:|---------------:|
| **DEG**          |      *k* = 18  |        *n* − *k* = 294    |    *n* = 312   |
| **Not DEG**      | *K* − *k* = 232 | *N* − *K* − (*n* − *k*) = 13,456 | *N* − *n* = 13,688 |
| **Column total** |     *K* = 250  |        *N* − *K* = 13,750 |   *N* = 14,000 |

```r
mat <- matrix(c(18, 294, 232, 13456), nrow = 2)
fisher.test(mat, alternative = "greater")$p.value
# [1] 1.29e-05   ← matches phyper
```

**Why we don't compute by hand.** The p-value is a sum of 233 terms (`i = 18, 19, …, 250`), each a ratio of three combinations on huge populations:

```
                      C(250, 18) · C(13,750, 294)
P(X = 18)  =   ─────────────────────────────────────────
                          C(14,000, 312)
```

The numbers involved are *enormous* — orders of magnitude beyond what fits in a 64-bit double. From R's `choose()` and `lchoose()`:

| Combination | Approximate value | log₁₀ | Digits |
|-------------|------------------:|------:|-------:|
| `C(250, 18)`           | `≈ 1.21 × 10²⁷`     |   27.08 |    28 |
| `C(13,750, 294)`       | `≈ 4.0 × 10⁶¹⁵`     |  615.6  |   616 |
| `C(14,000, 312)`       | `≈ 5.1 × 10⁶⁴⁷`     |  647.7  |   648 |

```r
choose(250, 18)               # 1.214e+27 — fits in a double
choose(13750, 294)            # Inf — overflow
choose(14000, 312)            # Inf — overflow

# Fix: work in log-space with lchoose()
lchoose(14000, 312)           # 1491.51 (natural log → 648 decimal digits)
exp(lchoose(K, k) + lchoose(N-K, n-k) - lchoose(N, n))   # 9.49e-06  (= P(X=18))
```

`C(14,000, 312)` alone has **about 648 digits** — multiplying these three numbers by hand is hopeless, but their **logs** add and subtract cleanly. R/Python evaluate every probability through `lgamma` (log-Gamma) and `lchoose`, then exponentiate at the end.

In practice the **first term** `P(X = 18) ≈ 9.5 × 10⁻⁶` already accounts for **~74 %** of the total tail probability `1.29 × 10⁻⁵` — the sum converges quickly because each successive ratio `P(X = i+1) / P(X = i)` shrinks rapidly when *i* is far above the mean.

**Python equivalent (`scipy.stats.hypergeom`):**

```python
from scipy.stats import hypergeom
# survival function sf(x) = P(X > x); we want P(X >= 18) = sf(17)
hypergeom.sf(17, 14000, 250, 312)
# 1.29e-05
```

**Sanity check — three equivalent computations in R:**

```r
N <- 14000; K <- 250; n <- 312; k <- 18

# (a) phyper one-liner
p1 <- phyper(k - 1, K, N - K, n, lower.tail = FALSE)

# (b) explicit sum of dhyper across the upper tail
p2 <- sum(dhyper(k:min(K, n), K, N - K, n))

# (c) Fisher exact, one-sided "greater"
p3 <- fisher.test(matrix(c(k, n - k, K - k, N - K - (n - k)), nrow = 2),
                   alternative = "greater")$p.value

c(phyper = p1, dhyper_sum = p2, fisher = p3)
# phyper      dhyper_sum    fisher
# 1.285e-05   1.285e-05     1.285e-05
```

All three agree — this is the *correct* way to verify any homemade ORA implementation.

---

## 4. The background trap

The background set *N* is **not** "all annotated genes in the genome" — it is the set of genes that **could have been called as DEGs** in your experiment. For RNA-Seq, this means genes expressed at detectable levels (e.g. TPM > 1 in at least one sample).

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
| Bonferroni  | **FWER — Family-Wise Error Rate** (probability of any false positive) | `p < 0.05 / 5000 = 10⁻⁵`                 | Very conservative |
| Benjamini-Hochberg (BH) | **FDR — False Discovery Rate** (expected false-positive proportion among called hits) | data-driven; ~5% of called hits expected to be FP | Standard for ORA |

In practice GO enrichment uses **BH-adjusted q-values < 0.05** — Bonferroni is too strict and erases real signal when terms are correlated (parents and children share genes).

### 5.1 Quick mental check

- Raw `p < 0.001` → 5,000 × 0.001 = **5 expected false positives**.
- If you find 25 terms at raw `p < 0.001`, FDR ≈ 5/25 = **20%** — too noisy.
- Tighten to raw `p < 10⁻⁴` → 0.5 expected FP. If 10 terms pass, FDR ≈ 5%.

---

## 5b. GO Slim — high-level summarization

The full GO has tens of thousands of terms — far too granular for a single-figure summary. **GO Slim** is a curated **reduced subset** that retains only broad, high-level terms (e.g. *DNA metabolic process*, *response to stress*, *transport*). Use it when:

- You want a pie chart / bar chart of the **distribution** of cellular components or biological processes among your DEGs.
- You're presenting to a non-specialist audience and 1,500 leaf terms would be unreadable.
- You're comparing across very different gene lists where leaf-level overlaps are tiny.

There are several pre-built slims: **generic GO Slim**, **plant GO Slim** (Plant Ontology consortium — used for *A. thaliana*), **Yeast Slim**, **AGR Slim**, etc.

```r
library(clusterProfiler)
library(GSEABase)

# Load the plant-specific GO Slim (download goslim_plant.obo from GOC)
slim <- getOBOCollection("goslim_plant.obo")

# Map your DEGs' annotations onto the slim
slim_map <- groupGOTerms(deg_go_terms, slim, ontology = "BP")

# In clusterProfiler:
ego_slim <- enrichGO(
  gene     = deg_genes,
  universe = expressed_set,
  OrgDb    = org.At.tair.db,
  keyType  = "TAIR",
  ont      = "BP",
  level    = 3        # use only ancestor terms at depth 3 — a poor-man's slim
)
```

A simpler approximation: restrict the enrichment result to GO levels 2–3 of the DAG (close to the root). This collapses *response to chitin* and *response to bacterium* into the parent *response to biotic stimulus* — adequate for high-level summary plots.

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

### 6.1 How GSEA actually works

ORA asks: "of my **discrete** DEG list, are term-X genes over-represented?" GSEA asks a different question: "walking down a fully ranked gene list, do term-X genes pile up at the top (or bottom) more than chance?" No threshold is needed — every gene contributes.

**Algorithm (Subramanian 2005):**

1. **Rank** every gene in the experiment by a chosen statistic — typically log₂ fold-change, the Wald statistic, or `-log10(p) * sign(LFC)`.
2. Walk down the ranked list. For each gene, update a **running sum**:
   - `+hit_w` if the gene is in gene-set *S* (weighted by its rank statistic),
   - `−miss_w` if it isn't.
3. The **enrichment score (ES)** is the maximum deviation of the running sum from zero.
4. **Normalize** ES across gene-sets of different sizes → **NES** (normalized enrichment score).
5. **Permute** sample labels (preferred) or gene labels (`gseGO()` default) thousands of times; compare the observed NES to the null distribution to obtain a permutation p-value, then BH-adjust.
6. The **leading-edge subset** is the set of genes from *S* that appear *before* the running-sum reaches its peak — these are the "drivers" of the enrichment.

```r
library(clusterProfiler)
library(org.At.tair.db)

# Build a named, sorted vector — required input format
geneList <- deseq_results$log2FoldChange
names(geneList) <- deseq_results$gene_id
geneList <- sort(geneList, decreasing = TRUE)

set.seed(42)                     # GSEA is permutation-based; seed for reproducibility
gsea <- gseGO(
  geneList      = geneList,
  OrgDb         = org.At.tair.db,
  keyType       = "TAIR",
  ont           = "BP",
  minGSSize     = 10,
  maxGSSize     = 500,
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  verbose       = FALSE
)

head(as.data.frame(gsea))
```

The result table columns:

| Column            | Meaning                                                                 |
|-------------------|-------------------------------------------------------------------------|
| `ID`, `Description` | GO term ID and label                                                  |
| `setSize`         | Number of genes in the gene-set that were also in `geneList`            |
| `enrichmentScore` | Raw ES (signed, ranges roughly −1 to +1)                                |
| `NES`             | Normalized ES — comparable across gene-sets                             |
| `pvalue`          | Permutation p-value                                                     |
| `p.adjust`        | BH-adjusted q-value                                                     |
| `rank`            | Position of the leading-edge peak in the ranked list                    |
| `core_enrichment` | Slash-separated leading-edge gene IDs — the drivers of the enrichment   |

**Sign conventions:** positive NES = gene-set members enriched at the **top** of the ranked list (up-regulated in your contrast). Negative NES = enriched at the **bottom** (down-regulated).

> ## ORA vs GSEA at a glance
>
> |                          | ORA (`enrichGO`)                | GSEA (`gseGO`)                       |
> |--------------------------|---------------------------------|--------------------------------------|
> | Input                    | Discrete DEG list + universe    | Full ranked vector of all genes      |
> | Statistical model        | Hypergeometric / Fisher 2×2     | Kolmogorov–Smirnov-like running-sum  |
> | Threshold needed?        | Yes (e.g. padj < 0.05)          | No                                   |
> | Sensitive to subtle, broad shifts? | No                    | Yes                                  |
> | Easy to explain in a paper figure? | Very                  | Less so (requires the running-sum plot) |
{: .callout}

### 6.2 Tools beyond R

R + clusterProfiler is the dominant academic toolchain, but several alternatives are worth knowing.

| Tool             | Type         | Strengths                                                        | Weaknesses                                |
|------------------|--------------|------------------------------------------------------------------|-------------------------------------------|
| **DAVID**        | Web          | Old, well-known, friendly UI; functional clustering              | Annotations updated infrequently; not reproducible without recording the version |
| **PANTHER**      | Web          | GOC's official front-end; PANTHER Overrepresentation Test = current gold-standard ORA; also offers GSEA-style "Statistical enrichment test" | Web-only by default; harder to script    |
| **g:Profiler**   | Web + R (`gprofiler2`) | Multi-database (GO, KEGG, Reactome, TF, miRNA, HPA) in one shot; current annotation snapshot pinned per release | Default background is whole-genome — must override |
| **AgriGO v2**    | Web          | Plant-specific (incl. *A. thaliana*); built-in singular and parent-child enrichment | Web-only; smaller community now           |
| **topGO**        | R/Bioconductor | Implements `weight01` and `elim` algorithms that account for the GO DAG when computing p-values, reducing redundancy of parent/child hits | Steeper API than clusterProfiler; no direct GSEA |
| **GOATOOLS**     | Python       | Scriptable; produces publication-quality DAG figures with `go_plot.py`; Klopfenstein 2018 paper details algorithm | Pure ORA, no GSEA                         |
| **enrichR**      | R + web      | Wraps Enrichr's 200+ libraries (drug, disease, TF) — wider than just GO | Mouse/human-centric; sparse for *A. thaliana* |

For *A. thaliana* coursework, **clusterProfiler + org.At.tair.db** is sufficient and reproducible. Consider topGO when you want DAG-aware p-values (it down-weights a parent term whose signal is fully explained by a more specific child).

### 6.3 KEGG and Reactome — sister enrichment frameworks

GO is a vocabulary of *functions and locations*; **KEGG** and **Reactome** are databases of **pathways** — defined sequences of molecular events. They overlap with GO BP but are not the same:

|              | GO Biological Process                              | KEGG / Reactome pathway                              |
|--------------|----------------------------------------------------|------------------------------------------------------|
| Granularity  | Conceptual ("response to cold")                    | Mechanistic ("CBF-COR cold response pathway")        |
| Members      | Curated list of genes performing the function      | Reaction graph: enzymes, substrates, products, regulators |
| Visualization| DAG of related terms                               | Wiring diagram of the pathway                        |
| When best    | First-pass, broad biological themes                | Confirming a specific pathway hypothesis; therapeutic targeting |

Run them in addition to GO when the biology is pathway-shaped:

```r
library(clusterProfiler)
library(ReactomePA)

# KEGG (uses KEGG-internal IDs; clusterProfiler converts)
kk <- enrichKEGG(
  gene         = entrez_ids,
  organism     = "ath",                 # Arabidopsis thaliana
  pvalueCutoff = 0.05
)

# Reactome (organism = "arabidopsis")
rr <- enrichPathway(
  gene         = entrez_ids,
  organism     = "arabidopsis",
  pvalueCutoff = 0.05,
  readable     = TRUE
)
```

Reactome's plant coverage is thinner than KEGG's, so for *A. thaliana* prefer KEGG when your pathway of interest is metabolic (photosynthesis, glucosinolate biosynthesis, flavonoid biosynthesis). Reactome is richer for human/mouse signaling.

---

## 7. Semantic similarity & redundancy

A common shock: the top-20 hits of a GO enrichment look almost identical — *response to stress*, *response to abiotic stimulus*, *response to chemical*, *cellular response to stimulus*, … These are not independent findings. They are **parents and grandparents of the same handful of leaf terms**, and they share most of their genes.

**Semantic similarity** measures attempt to quantify how much two GO terms overlap *in meaning* (and therefore in gene membership). Three classic measures:

- **Resnik (1995)** — based on the *information content* (IC) of the lowest common ancestor (LCA): rare terms carry more information than common terms.
- **Lin (1998)** — Resnik's IC normalized by the average IC of the two terms.
- **Wang (2007)** — accounts for DAG topology: each term is described by a vector of contributions from all its ancestors, weighted by edge type.

You don't need to compute these by hand. `clusterProfiler::simplify()` wraps `GOSemSim::mgoSim()` (Wang by default) and collapses redundant terms:

```r
ego_simple <- simplify(
  ego,
  cutoff      = 0.7,        # similarity threshold
  by          = "p.adjust", # which column to break ties on
  select_fun  = min,        # keep the term with the smallest q-value
  measure     = "Wang"
)
```

For a **publication-quality figure** with genuinely non-redundant terms, the standard pipeline is to feed your enrichment result into **REVIGO** ([http://revigo.irb.hr/](http://revigo.irb.hr/)) — paste the GO IDs and p-values, choose the species, and REVIGO returns a clustered scatter plot and a treemap. The R-native equivalent is `rrvgo` (Sayols 2023):

```r
library(rrvgo)

simMatrix <- calculateSimMatrix(
  ego$ID,
  orgdb   = "org.At.tair.db",
  ont     = "BP",
  method  = "Rel"
)
scores <- setNames(-log10(ego$p.adjust), ego$ID)
reduced <- reduceSimMatrix(simMatrix, scores, threshold = 0.7, orgdb = "org.At.tair.db")

treemapPlot(reduced)
scatterPlot(simMatrix, reduced)
```

> ## Always reduce before publishing
> Reviewers will ask "why are *response to stress*, *response to abiotic stimulus*, and *response to chemical* all in your top hits?" Run `simplify()` or REVIGO first; report both the full and reduced tables in your supplement.
{: .callout}

---

## 8. Visualization gallery

`clusterProfiler` and the companion `enrichplot` package give you a uniform plotting API on any enrichment result (`enrichResult` for ORA, `gseaResult` for GSEA).

```r
library(enrichplot)
```

| Plot          | Call                                       | Best for                                                                 |
|---------------|--------------------------------------------|--------------------------------------------------------------------------|
| `barplot`     | `barplot(ego, showCategory = 15)`          | Quick overview of top terms by count or gene ratio                       |
| `dotplot`     | `dotplot(ego, showCategory = 20)`          | Same data but encodes gene ratio (x), q-value (color), set size (dot size) — **the standard figure** |
| `cnetplot`    | `cnetplot(ego, foldChange = geneList)`     | Gene-concept network: which **genes** drive multiple terms              |
| `emapplot`    | `emapplot(pairwise_termsim(ego))`          | Enrichment map: terms with shared genes cluster — visualizes redundancy |
| `goplot`      | `goplot(ego, showCategory = 10)`           | Embed top terms in the GO DAG itself, color by significance              |
| `treeplot`    | `treeplot(pairwise_termsim(ego))`          | Hierarchical clustering of terms with text-summary labels per cluster    |
| `heatplot`    | `heatplot(ego, foldChange = geneList)`     | Genes × terms heatmap colored by fold-change                             |
| `gseaplot2`   | `gseaplot2(gsea, geneSetID = 1:3)`         | The signature GSEA running-sum plot for the top 3 sets                   |
| `ridgeplot`   | `ridgeplot(gsea)`                          | Density of fold-changes per enriched set (GSEA only)                     |

Worked example — a "publication strip" of three plots from one ORA result:

```r
library(patchwork)

p1 <- dotplot(ego, showCategory = 15) + ggtitle("Top BP terms (ORA)")
p2 <- cnetplot(ego, foldChange = geneList, showCategory = 8)
p3 <- emapplot(pairwise_termsim(ego), showCategory = 30)

(p1 | p2) / p3
ggsave("figures/go_summary.png", width = 14, height = 10, dpi = 300)
```

> ## Which plot to lead with?
> - **Talk slide:** `dotplot` — it answers "what's enriched?" in one glance.
> - **Paper main figure:** `dotplot` + `cnetplot` side by side — the network shows reviewers that you have multiple genes per term, not a single-gene fluke.
> - **Supplement:** `emapplot` to demonstrate redundancy structure; `goplot` to show DAG context.
{: .callout}

---

## 9. Worked example: heat shock RNA-Seq in *A. thaliana*

A 4-hour 38 °C heat shock on 2-week-old seedlings, contrasted with 22 °C controls (3 reps each). DESeq2 gives 487 DEGs at padj < 0.05. Biologically we **expect** to see *response to heat*, *protein folding*, *unfolded protein response*, and *chaperone* terms — this is the classic heat-shock-protein (HSP) response.

```r
library(DESeq2)
library(clusterProfiler)
library(org.At.tair.db)
library(enrichplot)

# 1. DEGs and ranked list (assume `res` is a DESeq2 result)
sig <- subset(as.data.frame(res), padj < 0.05)
deg_genes     <- rownames(sig)
expressed_set <- rownames(subset(as.data.frame(res), baseMean > 1))

# Ranked vector for GSEA
ranked <- res$log2FoldChange
names(ranked) <- rownames(res)
ranked <- sort(ranked[!is.na(ranked)], decreasing = TRUE)

# 2. ORA on biological process
ego <- enrichGO(
  gene          = deg_genes,
  universe      = expressed_set,
  OrgDb         = org.At.tair.db,
  keyType       = "TAIR",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.10,
  readable      = TRUE          # convert TAIR IDs → gene symbols in output
)

# 3. Make IDs human-readable in stored result (idempotent — `readable=TRUE` already did it)
ego <- setReadable(ego, OrgDb = org.At.tair.db, keyType = "TAIR")

# 4. Collapse redundant parent/child terms
ego_s <- simplify(ego, cutoff = 0.7, by = "p.adjust", select_fun = min)

# 5. Inspect — expect HSP-flavored top hits
head(as.data.frame(ego_s)[, c("Description", "GeneRatio", "p.adjust", "Count")], 10)
```

Expected top hits (real *A. thaliana* heat-shock results look like this):

```
Description                            GeneRatio   p.adjust   Count
response to heat                       42/487      3.2e-29     42
protein folding                        31/487      8.1e-21     31
response to temperature stimulus       58/487      1.4e-19     58
chaperone-mediated protein folding     14/487      6.0e-13     14
unfolded protein binding (MF)          22/487      4.5e-12     22
response to misfolded protein          11/487      2.8e-09     11
heat acclimation                        9/487      1.1e-07      9
```

```r
# 6. Plots for the figure
dotplot(ego_s, showCategory = 12) + ggtitle("Heat-shock DEG GO BP")
cnetplot(ego_s, foldChange = ranked, showCategory = 6, circular = TRUE)

# 7. GSEA — captures the broader thermal response, not just the hard cutoff
set.seed(7)
gsea <- gseGO(
  geneList      = ranked,
  OrgDb         = org.At.tair.db,
  keyType       = "TAIR",
  ont           = "BP",
  minGSSize     = 10,
  pAdjustMethod = "BH"
)
gseaplot2(gsea, geneSetID = head(gsea$ID, 3),
          title = "GSEA: top heat-shock processes")

# 8. Save outputs
write.csv(as.data.frame(ego_s), "results/go_BP_heatshock.csv", row.names = FALSE)
write.csv(as.data.frame(gsea),  "results/gsea_BP_heatshock.csv", row.names = FALSE)
saveRDS(ego_s, "results/ego_simplified.rds")
```

**Biological narrative.** The leading edge of *response to heat* (column `core_enrichment` in the GSEA table) is dominated by **HSP70**, **HSP90**, **HSP101**, **HSP17.6**, **HSFA2**, and **BAG6** — the canonical heat-shock factors and chaperones. The fact that *cellular response to misfolded protein* and *unfolded protein binding* both surface confirms the protein-quality-control arm is engaged, not just transcriptional induction. If your data instead returned *photosynthesis* and *response to light intensity*, you would suspect a **temperature-coupled light artifact** — heat-shock chambers often raise irradiance — and the GO result would have caught the confounder before you wrote the discussion.

---

## 10. Example

A worked walk-through of the full ORA workflow on the course's toy dataset:

1. Run DESeq2 on the toy A vs B counts in `~/scratch/rnaseq/counts/`.
2. Extract DEGs at `padj < 0.05`.
3. Run `enrichGO` **twice**: once with the whole-genome background and once with `universe = expressed_genes`. Compare top 10 terms — note shifts in q-value.
4. Visualize with `dotplot(ego)` and `cnetplot(ego)`.
5. Reflection: which background did the previous student probably use if their top hit was "ribosome biogenesis, q = 1e-30"?

---

## 11. Common pitfalls — checklist

- [ ] Background = expressed genes, not whole genome.
- [ ] Evidence codes filtered (drop IEA-only if you want experimental support).
- [ ] BH (q-value), not raw p, for cutoffs.
- [ ] `phyper(k - 1, ...)` — off-by-one.
- [ ] DAG semantics: redundant parent/child terms in the top hits often inflate apparent biology — collapse with `simplify(ego, cutoff = 0.7)`.
- [ ] Report **gene ratio** (k/n) and **bg ratio** (K/N) alongside p — large fold differences with tiny gene counts are noisy.
- [ ] Remember `regulates` does **not** propagate — *regulation of X* genes do not count toward *X* by default.
- [ ] Confirm the species annotation snapshot (`packageVersion("org.At.tair.db")`) — different snapshots can give different top hits and reviewers will ask.
- [ ] For GSEA, set a random seed (`set.seed(...)`) — the permutation step is stochastic.
- [ ] Run REVIGO or `rrvgo` before publication-quality figures; submit both full and reduced tables in your supplement.
- [ ] If you are using KEGG/Reactome alongside GO, report each separately — don't union p-values across orthogonal databases.
- [ ] When a result looks "too clean" (all top hits are direct ancestors of one leaf term), you are looking at DAG correlation, not multiple independent signals.

> ## Reading list
> - Ashburner et al. (2000) *Gene Ontology: tool for the unification of biology.* Nature Genetics 25:25–29.
> - Rhee, Wood, Dolinski, Draghici (2008) *Use and misuse of the Gene Ontology annotations.* Nature Reviews Genetics 9:509–515. — the GOC's own how-to-not-misuse manifesto.
> - Huang et al. (2009) *Bioinformatics enrichment tools: paths toward the comprehensive functional analysis of large gene lists.* Nucleic Acids Res. 37:1.
> - Khatri et al. (2012) *Ten years of pathway analysis: current approaches and outstanding challenges.* PLoS Comput Biol 8:e1002375.
> - Subramanian et al. (2005) *Gene set enrichment analysis: a knowledge-based approach for interpreting genome-wide expression profiles.* PNAS 102:15545.
> - Klopfenstein et al. (2018) *GOATOOLS: A Python library for Gene Ontology analyses.* Scientific Reports 8:10872.
> - Sayols (2023) *rrvgo: a Bioconductor package for interpreting lists of Gene Ontology terms.* MicroPubl. Biol. — the R-native REVIGO replacement.
> - Wijesooriya et al. (2022) *Urgent need for consistent standards in functional enrichment analysis.* PLoS Comput Biol 18:e1009935.
{: .callout}
