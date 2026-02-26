---
layout: page
title: ESM3 Structure-Conditioned Loop Design
published: true
---

{% include gh_variables.html %}

# ESM3 Structure-Conditioned Loop Design: TIMP3 C-Loop

---

## Assignment Setup

### Create your GitHub repository

1. Create a **private** GitHub repository named `ESM3-TIMP3-StructCond`
2. Initialize with a README.md
3. Add Dr. Yim as a collaborator (GitHub username: `wyim-pgl`)
4. Your repository structure should look like this when complete:

```
ESM3-TIMP3-StructCond/
    README.md
    environment.yml
    01_seq_only_generation.py
    02_struct_conditioned_generation.py
    03_validation.py
    04_analysis.ipynb          # Figures + written analysis
    results/
        seq_only_loops.csv
        struct_cond_loops.csv
        validation_metrics.csv
    figures/
        fig1_plddt_distribution.png
        fig2_rmsd_vs_plddt.png
        fig3_aa_frequency_heatmap.png
```

### Environment setup

```yaml
# environment.yml
name: esm3-timp3
channels:
  - conda-forge
  - pytorch
  - nvidia
dependencies:
  - python=3.10
  - pytorch>=2.0
  - pytorch-cuda=12.1
  - pip
  - pip:
    - esm
    - biotite
    - py3Dmol
    - matplotlib
    - seaborn
    - pandas
    - scipy
```

```bash
conda env create -f environment.yml
conda activate esm3-timp3
```

---

## CRITICAL: The BOS/EOS Token Indexing Trap

Read this before writing any code.

When you call `model.encode()`, ESM3 wraps the sequence with special tokens:

```
Your sequence:     M   C   P   Q   ...   res_62   ...   res_187
Position:          0   1   2   3   ...   62       ...   187

After encode():
Token index:  0     1   2   3   4   ...   63       ...   188     189
Token:        BOS   M   C   P   Q   ...   res_62   ...   res_187  EOS
              4098                                                 4097
```

**Sequence position `i` = token index `i + 1`**. Tensor length = `len(sequence) + 2`.

When copying structure tokens between tensors, either:
- Work entirely in token space (index-to-index copy), OR
- Work in residue space and add +1 for every index

Mixing conventions shifts everything by one residue. This is the #1 source of silent bugs in multi-track prompts.

### Special token values

| Token | Value | Meaning |
|-------|-------|---------|
| BOS | 4098 | Beginning of sequence (index 0) |
| EOS | 4097 | End of sequence (last index) |
| MASK | 4096 | Masked position (model generates here) |
| 0-4095 | -- | Structure codebook (real structural info) |

---

## Part 1: Sequence-Only Generation [15 pts]

### Question 1.1: ESM3 Sequence-Only Generation

Write `01_seq_only_generation.py` that:

1. Loads ESM3 (`ESM3.from_pretrained("esm3_sm_open_v1")`)
2. Loads the TIMP3 mature sequence (UniProt P35625, signal peptide removed)
3. Masks the C-loop at positions 62-67 with `_`
4. Generates 50 sequences using `GenerationConfig(track="sequence")`
5. Extracts the 6-mer loop from each generated sequence
6. Saves results to `results/seq_only_loops.csv` with columns: `index`, `full_sequence`, `loop_sequence`
7. Prints summary statistics: number of unique loops, wildtype matches, per-position entropy

> ## Note
> This replicates your current pipeline. The point is to have a clean baseline for comparison.
{: .callout}

> ## Hint
> For reproducibility, set `temperature=0.7` and `num_steps=8`. Save the model loading and sequence setup in a shared utility if you want to avoid duplication with Part 2.
{: .callout}

> ## Reference
> Your existing code in `Generation/Training_with_Gen.ipynb`, Step 5b.
{: .callout}

### Question 1.2 (Optional Bonus, +5 pts)

Build a PSSM baseline from your experimental C-loop library (~47,000 sequences). Compute a position frequency matrix for the 6 loop positions, sample 50 sequences from that distribution, and save to `results/pssm_loops.csv`.

If you complete this, include the PSSM baseline in all downstream comparisons (Parts 3-4) for a 3-way analysis.

> ## Note
> This directly feeds into Q4.4. If you have the data, your answer to "when does a PSSM beat ESM3?" becomes empirical instead of speculative.
{: .callout}

> ## Hint
> With 85% positive binders and 47K sequences, the PSSM captures the experimental fitness landscape well. The question is whether ESM3 adds anything beyond what this simple model already knows.
{: .callout}

---

## Part 2: Structure-Conditioned Generation [25 pts]

### Question 2.1: Load the template structure

Write the structure-loading section of `02_struct_conditioned_generation.py`.

Download PDB 3CKI and extract the N-TIMP3 chain.

> ## Note
> 3CKI is the TACE(ADAM17)-N-TIMP3 complex (Wisniewska et al., 2008, *J. Mol. Biol.*). It has two chains:
> - **Chain A**: TACE catalytic domain
> - **Chain B**: N-TIMP3 (~125 residues, N-terminal domain only)
>
> This is NOT the full-length TIMP3 (188 residues). The C-loop (positions 62-67) falls within the N-terminal domain, so 3CKI covers your target region.
{: .callout}

> ## CAUTION
> PDB residue numbering may differ from your full-length sequence indexing. After loading, align the 3CKI chain B sequence against your TIMP3 sequence to find the correct position mapping. Do not assume they match.
{: .callout}

> ## Hint
> ```python
> timp3_chain = ProteinChain.from_rcsb("3CKI", chain_id="B")
> template = ESMProtein.from_protein_chain(timp3_chain)
> template_tokens = model.encode(template)
> ```
{: .callout}

### Question 2.2: Construct the multi-track prompt

In the same script, build a prompt with:
- **Sequence track**: C-loop masked, everything else is wildtype
- **Structure track**: C-loop masked (4096), flanking regions have real structure tokens copied from the 3CKI template

Your code must include a sanity check that prints:
- Number of structure-conditioned positions (should be all non-loop residues)
- Number of masked positions (should be exactly 6)
- A visual map showing which positions have structure conditioning

> ## CAUTION
> This is where the BOS/EOS offset matters. If your sanity check shows anything other than 6 masked structure positions, your indexing is wrong. Go back and re-read the BOS/EOS section above.
{: .callout}

> ## Hint
> Initialize an all-masked structure track, set BOS/EOS, then copy flanking tokens:
> ```python
> prompt_tokens.structure = torch.full_like(prompt_tokens.sequence, 4096)
> prompt_tokens.structure[0] = 4098   # BOS
> prompt_tokens.structure[-1] = 4097  # EOS
> # Copy flanking structure: residue i -> token index i + 1
> ```
{: .callout}

> ## Hint
> Sanity check pattern:
> ```python
> n_masked = (prompt_tokens.structure == 4096).sum().item()
> assert n_masked == 6, f"Expected 6 masked, got {n_masked}"
> ```
{: .callout}

### Question 2.3: Single-pass multi-track generation

Using the multi-track prompt from Q2.2, generate 50 sequences. Both sequence and structure tracks are masked at the loop positions, so ESM3 fills in both simultaneously in a single `generate()` call.

Save results to `results/struct_cond_loops.csv` with the same format as Part 1.

> ## Note
> The key difference from Part 1 is that the structure track now carries real structural context from 3CKI for the flanking regions. ESM3 sees the loop's structural neighborhood when generating the sequence, not just the amino acid context.
{: .callout}

> ## Hint
> Use `temperature=0.7` and `num_steps=8` (same as Part 1 for fair comparison).
{: .callout}

> ## Reference
> ESM3 GFP tutorial -- [https://virtualcellmodels.cziscience.com/tutorial/esm3-tutorial](https://virtualcellmodels.cziscience.com/tutorial/esm3-tutorial)
{: .callout}

---

## Part 3: Structural Validation [25 pts]

### Question 3.1

Write `03_validation.py` that:

1. Reads both CSV files from `results/`
2. Folds every generated sequence using ESM3 (`track="structure"`, `temperature=0.0` for deterministic prediction)
3. For each folded structure, calculates:
   - Mean pLDDT of the C-loop region (positions 62-67)
   - Loop backbone RMSD vs. the 3CKI template (superimpose on non-loop residues first, then measure loop RMSD)
   - Overall backbone RMSD
4. Saves all metrics to `results/validation_metrics.csv` with columns: `method`, `index`, `loop_sequence`, `loop_plddt`, `loop_rmsd`, `backbone_rmsd`

> ## Note
> Superimpose on the conserved (non-loop) backbone first, THEN measure loop RMSD without re-superimposing. This isolates loop geometry differences from global alignment noise.
{: .callout}

> ## Hint
> For pLDDT extraction:
> ```python
> folded_chain = folded_protein.to_protein_chain()
> loop_plddt = folded_chain.confidence[62:68].mean()
> ```
{: .callout}

> ## Hint
> For RMSD with biotite:
> ```python
> from biotite.structure import superimpose, rmsd
> # Select Ca atoms, exclude loop for superimposition, include loop for RMSD
> ```
{: .callout}

### Question 3.2

In `04_analysis.ipynb`, create three figures and save them to `figures/`:

**Figure 1** (`fig1_plddt_distribution.png`): Loop pLDDT distribution
- Box plot or violin plot, one per method
- Horizontal line at pLDDT = 70

**Figure 2** (`fig2_rmsd_vs_plddt.png`): Loop RMSD vs. pLDDT scatter
- Color by method
- Label the high-pLDDT / low-RMSD quadrant as the ideal region

**Figure 3** (`fig3_aa_frequency_heatmap.png`): Amino acid frequency comparison
- Side-by-side heatmaps (6 positions x 20 amino acids)
- Same color scale for both panels

> ## Hint
> Use `matplotlib.pyplot.subplots(1, 2, ...)` with `sharey=True` for the heatmaps so the scale is identical.
{: .callout}

---

## Part 4: Critical Analysis [30 pts]

Answer these in `04_analysis.ipynb` as markdown cells. Each answer: 3-5 sentences, supported by specific numbers from your results.

### Question 4.1: Sequence divergence [6 pts]

Does structure conditioning produce measurably different loop sequences compared to sequence-only generation? Show evidence from your amino acid frequency data. If the distributions are similar, what does that tell you about what ESM3's sequence track already knows about structure?

### Question 4.2: Structural quality [6 pts]

Compare pLDDT and RMSD distributions between methods. Define a threshold for "good" (justify it), report the pass rate for each method, and run a statistical test (Mann-Whitney U or t-test).

### Question 4.3: The 85% baseline problem [6 pts]

The experimental TIMP3 C-loop library has ~47,000 sequences with ~85% positive binders. That means almost any 6-mer works. How would you demonstrate that ESM3 generates sequences that are genuinely better than random? What metric beyond binary binder/non-binder would you use?

> ## Hint
> Think about what "better" could mean -- binding affinity (Ki), thermal stability, specificity for one MMP over another, resistance to proteolysis. The binary hit rate is not enough when the baseline is 85%.
{: .callout}

### Question 4.4: Model appropriateness [6 pts]

1.4 billion parameters for 6 amino acids. The experimental data already gives you a PSSM from 47,000 sequences. When would the PSSM win? When would ESM3 win? Consider data availability, in-distribution vs. out-of-distribution design, and compute cost.

### Question 4.5: Scaling beyond loop redesign [6 pts]

If you applied multi-track conditioning to a de novo protein design (not loop redesign), what changes? In this assignment you used single-pass generation where ESM3 fills sequence and structure simultaneously. Under what circumstances would a two-stage approach (generate structure first, then design sequence) be more appropriate? What fails when the target fold has no close homolog in ESM3's training data?

> ## Reference
> Hayes et al. (2025), particularly the discussion of how esmGFP was designed at 58% identity to known GFPs -- a case where structure conditioning was essential because sequence homology alone was insufficient.
{: .callout}

---

## Submission

1. Push all code, results, and figures to your GitHub repository
2. Ensure `README.md` includes:
   - Brief description of the project
   - How to set up the environment (`conda env create -f environment.yml`)
   - How to run each script in order
   - Summary of key findings (2-3 sentences)
3. Verify that all scripts run without errors on a clean checkout
4. Submit the repository URL

**Due date**: Check the course schedule

---

## Grading

| Component | Points |
|-----------|--------|
| Part 1: `01_seq_only_generation.py` + output CSV | 15 |
| Part 2: `02_struct_conditioned_generation.py` + output CSV + sanity checks | 25 |
| Part 3: `03_validation.py` + metrics CSV + 3 figures | 25 |
| Part 4: Written analysis (Q4.1-Q4.5, 6 pts each) | 30 |
| Code quality, README, repo organization | 5 |
| **Total** | **100** |

---

## References

- Hayes, T. et al. (2025). Simulating 500 million years of evolution with a language model. *Science*. DOI: 10.1126/science.ads0018
- Wisniewska, M. et al. (2008). Structural determinants of the ADAM inhibition by TIMP-3: crystal structure of the TACE-N-TIMP-3 complex. *J. Mol. Biol.*, 381(5), 1307-1319. PDB: 3CKI
- ESM3 GitHub: [https://github.com/evolutionaryscale/esm](https://github.com/evolutionaryscale/esm)
- ESM3 GFP Tutorial: [https://virtualcellmodels.cziscience.com/tutorial/esm3-tutorial](https://virtualcellmodels.cziscience.com/tutorial/esm3-tutorial)
- Lin, Z. et al. (2023). Evolutionary-scale prediction of atomic-level protein structure with a language model. *Science*, 379(6637), 1123-1130.

---

## Compute Notes

- 50 sequences per method on A100: ~30-60 min
- Structure-conditioned takes ~2x longer (two generation stages)
- If using T4/V100, reduce to 20 sequences and note in README
- Folding 100 sequences (validation): ~20-30 min additional
