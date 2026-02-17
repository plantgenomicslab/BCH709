---
layout: page
title: Vibe Coding in Life Sciences
published: true
---

{% include gh_variables.html %}

> ## Source Attribution
> This content is adapted from ["When Vibe Coding Meets Life Science"](https://www.linkedin.com/pulse/when-vibe-coding-meets-life-science-gozde-eskici-phd)
> by Gozde Eskici, Ph.D. (The Second Translation newsletter, April 14, 2025)
{: .callout}

---

## Table of Contents

1. [What is Vibe Coding?](#what-is-vibe-coding)
2. [Vibe Coding in Bioinformatics](#vibe-coding-in-bioinformatics)
3. [8 Tools Rewriting the Rules](#8-tools-rewriting-the-rules-of-life-sciences)
4. [Setting Up VS Code with AI Assistants](#setting-up-vs-code-with-ai-coding-assistants)
5. [BCH709 Lab Materials](#bch709-bioinformatics-vibe-coding-lab-materials)
   - [Step 0A: Brainstorming (Python)](#step-0a-brainstorming-prompt-analysis-1--python)
   - [Step 0A: Brainstorming (R)](#step-0a-brainstorming-prompt-analysis-2--r)
   - [Step 0B: How to Ask AI to Set Up Your Environment](#step-0b-how-to-ask-ai-to-set-up-your-environment)
     - [Environment Setup Prompt (Python)](#step-0b-environment-setup-prompt-analysis-1--python)
     - [Environment Setup Prompt (R)](#step-0b-environment-setup-prompt-analysis-2--r)
   - [Step 0C: Environment Creation Commands](#step-0c-environment-creation-commands)
   - [Step 0D: Research Project Design (Advanced)](#step-0d-research-project-design-prompt-advanced)
   - [Telling AI About Your Environment](#telling-ai-assistants-about-your-conda-environment)
     - [Persistent Configuration (Claude `/init`, Copilot, Gemini, ChatGPT)](#persistent-environment-configuration-per-ai-assistant)
     - [Unofficial Config Files (GEMINI.md, CODEX.md)](#unofficial-configuration-files-geminimd-codexmd)
     - [Single-Conversation Templates](#single-conversation-chat-message-templates)
   - [How to Write Effective Prompts](#how-to-write-effective-vibe-coding-prompts)
   - [Part 1: Vibe Coding Examples](#part-1-vibe-coding-examples)
     - [Example 1: Python GFF3 Analysis](#example-1-python-per-chromosome-feature-counts-from-mgi-gff3--qc)
     - [Example 2: R TPM Heatmap](#example-2-r-top-200-variable-genes-from-ice-plant-tpm--heatmap)
   - [Part 2: Homework Assignments](#part-2-homework-assignments)
     - [Homework 1: Python FASTA Analysis](#homework-1-python-mrna-fasta-analysis--gc-distribution-graph)
     - [Homework 2: R Clustering](#homework-2-r-z-score-clustering-of-cv-top-200-genes--pattern-visualization)
6. [Appendix: Prompt Templates](#appendix-effective-vibe-coding-prompt-template)
7. [Input/Output Prompt Checklist](#inputoutput-prompt-checklist)

---

## What is Vibe Coding?

On February 2nd, 2025, **Andrej Karpathy**, one of the most influential voices in AI, introduced a new term:

> "There's a new kind of coding I call 'vibe coding', where you fully give in to the vibes, embrace exponentials, and forget that the code even exists. It's possible because the LLMs (e.g. Cursor Composer w Sonnet) are getting too good. Also I just talk to Composer with SuperWhisper so I barely even touch the keyboard... I'm building a project or webapp, but it's not really coding – I just see stuff, say stuff, run stuff, and copy paste stuff, and it mostly works."

**IBM** followed with a formal definition:

> "Vibe coding is a fresh take in coding where users express their intention using plain speech and the AI transforms that thinking into executable code."

---

## Vibe Coding in Bioinformatics

In bioinformatics, coding isn't about building websites—it's about running genome pipelines, analyzing RNA-seq data, or scripting variant calling workflows. Historically, that's meant technical depth, time, and a dedicated computational team.

But what if a scientist could just say:

- "Compare these RNA-seq datasets."
- "Predict disease progression from this clinical data."
- "Simulate how this protein binds ligands."

...and the AI handles the rest?

Thanks to LLMs, Biopython, and Colab-powered interfaces, we're now close. The act of building has become more conversational, more iterative—more "vibey."

### Why This Matters

Bioinformatics has long been bottlenecked by translation—the gap between biological question and computational answer. Vibe coding changes that by:

| Benefit | Description |
|---------|-------------|
| **Faster iteration** | Rapid prototyping of experiments and product ideas |
| **Lower barriers** | Scientists can code without deep programming expertise |
| **Broader access** | More people can prototype, test, and scale ideas |
| **Leaner teams** | Smaller teams can accomplish more, especially at early stages |

---

## 8 Tools Rewriting the Rules of Life Sciences

### 1. Superbio.ai – No-Code AI Marketplace

Founded by Berke Buyukkucak and Ronjon Nag. Run cutting-edge AI tools for drug discovery, protein design, and literature review—no code needed.

**Link:** [superbio.ai](https://superbio.ai)

### 2. Recursion's LOWE – LLM-Orchestrated Wet Lab

Recursion's internal tool (unveiled by Chris Gibson): describe an assay, LOWE designs and executes it via robotics using their proprietary phenomics and chemistry stack.

**Link:** [recursion.com](https://recursion.com)

### 3. DrBioRight 2.0 – Cancer Proteomics Chatbot

Built at MD Anderson by the Han Liang Lab. Ask questions like "Which proteins in pathway X are altered in this tumor?" and get real answers with plots.

**Publication:** [Nature Communications (2025)](https://www.nature.com/articles/s41467-025-56650-0) | **Link:** [drbioright.org](https://drbioright.org)

### 4. BioChatter – Open Source Bio-AI Toolkit

From EMBL-EBI. Build custom AI assistants that connect to APIs, databases, and bio tools. Fully open-source and on-prem ready.

**Link:** [biochatter.org](https://biochatter.org)

### 5. OLAF – Conversational Bioinformatics OS

From Weill Cornell (Dylan Riffle et al.). Say "Analyze this RNA-seq file" and OLAF writes the code, runs it, and returns transparent, inspectable results.

**Publication:** [arXiv](https://arxiv.org/abs/2503.12465)

### 6. TinyBio – ChatGPT for Scientists

Acquired by Seqera. Started by Sasha Dagayev and Vishal Patel in 2022. Real-time code execution supporting 50+ bio libraries with self-healing error correction.

**Link:** [tinybio.cloud](https://tinybio.cloud)

### 7. Scispot (Scibot) – Lab AI Analyst

YC-backed. Their AI assistant Scibot makes lab data conversational: "Summarize this week's PCR results" produces instant dashboards.

**Link:** [scispot.com](https://scispot.com)

### 8. Synthace – Conversational Wet Lab Automation

Describe experiments in plain English; AI generates protocols and sends them directly to lab robots.

**Link:** [synthace.com](https://synthace.com)

> ## Key Takeaways
> - **Vibe coding** lets users build through intent, not syntax
> - In **biotech**, that means less friction, faster feedback, and broader access
> - These tools don't just "assist" scientists—they enable more with less code and more creativity
{: .callout}

---

## Setting Up VS Code with AI Coding Assistants

VS Code is the recommended editor for vibe coding. By installing AI extensions, you turn it into a conversational coding environment where you can write prompts, generate code, and iterate — all in one place.

```
 ┌──────────────────────────────────────────────────────────────────────────────┐
 │                         AI Coding Assistants                                 │
 │                                                                              │
 │   VS Code Extensions:                          Web-Based:                   │
 │   ┌────────────┐ ┌────────────┐ ┌────────────┐ ┌────────────┐              │
 │   │  Claude    │ │  GitHub   │ │  Gemini    │ │  ChatGPT / │              │
 │   │ (Anthropic)│ │  Copilot  │ │ Code Assist│ │  Codex     │              │
 │   │            │ │  (OpenAI) │ │  (Google)  │ │  (OpenAI)  │              │
 │   └─────┬──────┘ └─────┬─────┘ └─────┬──────┘ └─────┬──────┘              │
 │         │              │              │              │                      │
 │         └──────────────┼──────────────┼──────────────┘                      │
 │                        │              │                                      │
 │              ┌─────────▼──────────────▼──────────┐                          │
 │              │   Your prompt in plain English     │                          │
 │              └──────────────────┬─────────────────┘                          │
 │                                 │                                            │
 │              ┌──────────────────▼─────────────────┐                          │
 │              │   AI-generated Python / R code     │                          │
 │              └──────────────────┬─────────────────┘                          │
 │                                 │                                            │
 │              ┌──────────────────▼─────────────────┐                          │
 │              │   Execute & inspect results        │                          │
 │              └────────────────────────────────────┘                          │
 └──────────────────────────────────────────────────────────────────────────────┘
```

### Install Visual Studio Code

> ## Windows (WSL)
> 1. Download VS Code from [https://code.visualstudio.com/](https://code.visualstudio.com/)
> 2. Install on Windows (not inside WSL)
> 3. Install the **WSL** extension in VS Code
> 4. Open WSL terminal and type `code .` to launch VS Code connected to WSL
{: .solution}

> ## macOS
> 1. Download VS Code from [https://code.visualstudio.com/](https://code.visualstudio.com/)
> 2. Move to Applications folder
> 3. Open VS Code, press `Cmd+Shift+P`, type "Shell Command: Install 'code' command in PATH"
> 4. Now you can use `code .` from Terminal
{: .solution}

> For detailed VS Code configuration with conda environments, see the [Software Installation lesson](compile.html#using-conda-environments-in-vs-code).
{: .callout}

### Extension 1: Claude (Anthropic)

[![Claude VS Code Extension](https://img.shields.io/badge/VS_Code-Claude_(Anthropic)-7C3AED?style=for-the-badge&logo=anthropic&logoColor=white)](https://marketplace.visualstudio.com/items?itemName=anthropic.claude-code)

[Claude](https://marketplace.visualstudio.com/items?itemName=anthropic.claude-code) provides a chat panel and inline code generation powered by Anthropic's Claude models. It excels at understanding large code contexts and following detailed instructions.

**Install:**
```bash
$ code --install-extension anthropic.claude-code
```

**Setup:**
1. Open VS Code and click the Claude icon in the sidebar
2. Sign in with your Anthropic account or enter an API key from [console.anthropic.com](https://console.anthropic.com/)
3. Start a chat and paste your prompt

> ## Claude Code (CLI Alternative)
> Claude is also available as a command-line tool for terminal-based workflows:
> ```bash
> $ npm install -g @anthropic-ai/claude-code
> $ claude
> ```
> This is useful for working directly in the terminal without VS Code.
{: .solution}

### Extension 2: GitHub Copilot (OpenAI Codex)

[![GitHub Copilot Extension](https://img.shields.io/badge/VS_Code-GitHub_Copilot_(OpenAI)-000000?style=for-the-badge&logo=github&logoColor=white)](https://marketplace.visualstudio.com/items?itemName=GitHub.copilot)

[GitHub Copilot](https://marketplace.visualstudio.com/items?itemName=GitHub.copilot) provides real-time inline autocomplete suggestions as you type. [Copilot Chat](https://marketplace.visualstudio.com/items?itemName=GitHub.copilot-chat) adds a conversational chat panel for writing prompts.

**Install:**
```bash
$ code --install-extension GitHub.copilot
$ code --install-extension GitHub.copilot-chat
```

**Setup:**
1. You need a GitHub account
2. Open VS Code and sign in to GitHub when prompted
3. Copilot starts suggesting code automatically as you type; use `Tab` to accept

> ## Free for Students
> GitHub Copilot is **free** for verified students through [GitHub Education](https://education.github.com/).
> Apply with your university email (`.edu`) to get access.
{: .callout}

### Extension 3: Gemini Code Assist (Google)

[![Gemini Code Assist Extension](https://img.shields.io/badge/VS_Code-Gemini_Code_Assist_(Google)-4285F4?style=for-the-badge&logo=google&logoColor=white)](https://marketplace.visualstudio.com/items?itemName=google.gemini-code-assist)

[Gemini Code Assist](https://marketplace.visualstudio.com/items?itemName=google.gemini-code-assist) provides AI-powered code generation and a chat panel, backed by Google's Gemini models.

**Install:**
```bash
$ code --install-extension google.gemini-code-assist
```

**Setup:**
1. Open VS Code and click the Gemini icon in the sidebar
2. Sign in with your Google account
3. Start a chat — a free usage tier is available for individual developers

### ChatGPT and Codex (OpenAI) — Web-Based Alternative

[![ChatGPT](https://img.shields.io/badge/Web-ChatGPT_(OpenAI)-412991?style=for-the-badge&logo=openai&logoColor=white)](https://chatgpt.com/)
[![Codex](https://img.shields.io/badge/Web-Codex_(OpenAI)-412991?style=for-the-badge&logo=openai&logoColor=white)](https://chatgpt.com/codex)

You don't need VS Code to do vibe coding. [ChatGPT](https://chatgpt.com/) and [Codex](https://chatgpt.com/codex) are web-based tools by OpenAI that let you write prompts and generate code directly in the browser.

**ChatGPT:**
- Go to [chatgpt.com](https://chatgpt.com/) and sign in with an OpenAI account
- Paste your prompt; ChatGPT generates code you can copy into your editor or terminal
- Free tier available; Plus subscription unlocks GPT-4o and longer context

**Codex (OpenAI):**
- Available at [chatgpt.com/codex](https://chatgpt.com/codex)
- Specialized for code generation tasks
- Can execute code in a sandboxed environment and return results
- Requires ChatGPT Plus or Pro subscription

> ## When to Use Web-Based Tools vs. VS Code Extensions
>
> | Use Case | Recommended Tool |
> |----------|-----------------|
> | Quick one-off code generation | ChatGPT (web) |
> | Iterating on code in a project | VS Code + Claude / Copilot / Gemini |
> | Running code in a sandboxed cloud environment | Codex (web) |
> | Working on HPC cluster via terminal | Claude Code (CLI) |
{: .callout}

### How AI Assistants Fit into the Vibe Coding Workflow

```
 Step 1             Step 2              Step 3             Step 4
 ┌──────────┐      ┌──────────────┐    ┌──────────────┐   ┌──────────────┐
 │  Write    │      │  AI generates│    │  Run code in │   │  Check       │
 │  prompt   │─────▶│  code in     │───▶│  terminal or │──▶│  output and  │
 │  in chat  │      │  editor      │    │  notebook    │   │  revise      │
 │  panel    │      │              │    │              │   │  prompt      │
 └──────────┘      └──────────────┘    └──────────────┘   └───────┬──────┘
                                                                   │
      ◀────────────────────────────────────────────────────────────┘
                            Iterate until correct
```

> ## Typical Session (VS Code)
> 1. Open VS Code with your conda environment active
> 2. Open the AI chat panel (Claude, Copilot, or Gemini)
> 3. Paste your structured prompt (environment + input + task + output specs)
> 4. Review the generated code, click "Insert at Cursor" or copy to a `.py` / `.R` file
> 5. Run the script in the integrated terminal
> 6. Inspect results; refine the prompt if needed
{: .callout}

> ## Typical Session (Web-Based: ChatGPT / Codex)
> 1. Open [chatgpt.com](https://chatgpt.com/) or [chatgpt.com/codex](https://chatgpt.com/codex) in your browser
> 2. Paste your structured prompt
> 3. Copy the generated code into your local editor or terminal
> 4. Run the script in your conda environment: `conda activate bch709 && python script.py`
> 5. Inspect results; return to ChatGPT and refine the prompt if needed
{: .callout}

### Comparison: AI Coding Assistants

| Feature | Claude | GitHub Copilot | Gemini Code Assist | ChatGPT / Codex |
|---------|--------|----------------|-------------------|-----------------|
| **Provider** | Anthropic | GitHub / OpenAI | Google | OpenAI |
| **Type** | VS Code extension + CLI | VS Code extension | VS Code extension | Web-based |
| **Authentication** | API key or Anthropic account | GitHub account | Google account | OpenAI account |
| **Free for students** | Usage-based pricing | Free via GitHub Education | Free tier available | Free tier (GPT-4o mini) |
| **Inline autocomplete** | Yes | Yes | Yes | N/A (web) |
| **Chat panel** | Yes | Yes | Yes | Yes (browser) |
| **Code execution** | Via terminal | Via terminal | Via terminal | Codex sandbox |
| **Best for** | Detailed prompts; multi-file context | Real-time autocomplete | Google Cloud integration | Quick generation; no setup |

### Quick Install: All VS Code Extensions

Install all extensions in one command:
```bash
$ code --install-extension anthropic.claude-code && \
  code --install-extension GitHub.copilot && \
  code --install-extension GitHub.copilot-chat && \
  code --install-extension google.gemini-code-assist
```

> ## Which One Should I Use for BCH709?
> You can use **any** of these AI assistants for the lab exercises and homework. The prompts in this lesson are written in plain English and work with all AI coding tools — VS Code extensions and web-based tools alike.
>
> **Recommendation:** Try multiple tools during the semester and compare the results. Different AI models produce different code for the same prompt — that's part of the learning experience.
{: .callout}

---

# BCH709 Bioinformatics Vibe Coding Lab Materials

> ## Lab Overview
> - **Audience:** BCH709 Genome Informatics graduate students
> - **Goal:** Experience how prompt specificity transforms code quality and output
> - **Core Lesson:** The more specific your prompt, the closer the result to what you actually need
> - **Structure:** 2 Examples (Python, R) + 2 Homework Assignments (Python, R)
{: .prereq}

## The Vibe Coding Workflow

```
Natural-language prompt → AI generates code → Execute → Inspect results → Revise prompt → Repeat
```

> ## The Key Insight
> **"Saying exactly what you want"** is the core skill.
> A prompt controls not only the code but also the **execution environment**—without that, reproducibility breaks down.
{: .callout}

---

## Step 0: Brainstorming and Environment Setup

> ## Learning Objective
> Before writing any code, ask the AI about possible approaches and required tools first.
{: .objectives}

## Step 0A. Brainstorming Prompt (Analysis 1 — Python)

Copy and paste this prompt into the AI first. **This is a strategy question, not a code request.**

~~~
I am a beginner student in BCH709.
Before writing any code, brainstorm the approaches and libraries I need for the following analysis.

Analysis (Python, GFF3 analysis):
- Input: MGI.gff3.gz (GFF3, gzip), chrom.sizes (TSV: chrom, length_bp)
- Goal: Count genes, exons (preventing isoform overcounting), snRNAs, and lncRNAs per chromosome; compute density
- Output: TSV table + dropped_seqids.txt (QC artifact)

Requirements:
1) Break the analysis into functional units (input parsing, statistical computation, visualization, file output, QC).
2) For each functional unit, suggest 1–2 candidate libraries.
3) Pick one recommended combination for beginners and explain why.
4) List the exact conda-forge package names for that combination.
~~~

## Step 0A. Brainstorming Prompt (Analysis 2 — R)

Now do the same for the R analysis. Copy and paste this prompt into a **new conversation** (or continue the same one).

~~~
I am a beginner student in BCH709.
Before writing any code, brainstorm the approaches and libraries I need for the following analysis.

Analysis (R, RNA-Seq TPM analysis):
- Input: iceplant_TPM_DT_ZT.tab.gz (gzip TSV, gene_id + sample TPM values)
- Goal: Select top 200 genes by CV (with mean >= 1 filter), save log2(TPM+1) heatmap
- Output: TSV + heatmap PNG

Requirements:
1) Break the analysis into functional units (input parsing, statistical computation, visualization, file output, QC).
2) For each functional unit, suggest 1–2 candidate libraries.
3) Pick one recommended combination for beginners and explain why.
4) List the exact conda-forge package names for that combination.
~~~

> ## Why Brainstorming First?
> - You don't need to memorize package names — the AI suggests them
> - The AI produces a structured "function → package" mapping you can review
> - The conda install commands follow naturally in the next step (Step 0B)
{: .callout}

## Step 0B. How to Ask AI to Set Up Your Environment

> ## Learning Objective
> Instead of memorizing conda commands, learn to **describe what you need** and let the AI generate the installation plan for you.
{: .objectives}

Setting up a conda environment is a three-step process:

```
Step 1: Tell the AI what you want to do       (Step 0A — Brainstorming)
        ↓
Step 2: Ask the AI to generate install commands (Step 0B — Environment Prompt)
        ↓
Step 3: Copy-paste and run the commands         (Step 0C — Install)
```

> ## Key Idea
> You already told the AI what analysis you want to do in **Step 0A** (brainstorming). The AI knows which libraries it recommended — you don't need to list them again. Simply ask: **"Based on what you recommended, give me the install commands."**
{: .callout}

### Step 0B. Environment Setup Prompt (Analysis 1 — Python)

Once brainstorming is complete, use this prompt:

~~~
Using the library combination you just recommended, generate conda environment creation commands.

Conditions:
- Python environment name: bch709-python
- Pin Python 3.11
- Include import/library verification tests after installation
- Present the commands in copy-paste order so a beginner can just run them one by one
~~~

### Step 0B. Environment Setup Prompt (Analysis 2 — R)

Once brainstorming is complete, use this prompt:

~~~
Using the library combination you just recommended, generate conda environment creation commands.

Conditions:
- R environment name: bch709-R
- Pin R 4.3
- Include import/library verification tests after installation
- Present the commands in copy-paste order so a beginner can just run them one by one
~~~

## Step 0C. Environment Creation Commands

The AI will generate commands like the ones below. **Copy-paste and run them in your terminal.**

> ## Important
> The commands below are what the AI typically produces. Your results may vary slightly depending on which AI assistant you use — that's fine as long as the verification step passes.
{: .callout}

### Python Environment: `bch709-python`

```bash
# Create environment
conda create -n bch709-python -y python=3.11
conda activate bch709-python

# Install packages
conda install -c conda-forge -y \
  pandas numpy matplotlib seaborn biopython tqdm

# Verify installation
python -c "import pandas, numpy, matplotlib, Bio; print('bch709-python OK')"

# Create working directories
mkdir -p data results
```

### R Environment: `bch709-R`

```bash
# Create environment
conda create -n bch709-R -y -c conda-forge \
  r-base=4.3 r-data.table r-ggplot2 r-pheatmap r-viridislite r-scales
conda activate bch709-R

# Verify installation
R -q -e 'library(data.table); library(ggplot2); library(pheatmap); cat("bch709-R OK\n")'

# Create working directories
mkdir -p results
```

> ## Troubleshooting
> If the verification step fails:
> 1. Check that you activated the correct environment: `conda activate bch709-python`
> 2. Re-run the install command — sometimes packages fail to download on the first try
> 3. Ask the AI: *"I got this error when verifying: [paste error]. How do I fix it?"*
{: .callout}

### Data Downloads

```bash
# GFF3 (Example 1, Homework 1)
curl -L -o data/MGI.gff3.gz http://www.informatics.jax.org/downloads/mgigff3/MGI.gff3.gz

# Chromosome sizes (Example 1)
curl -L -o data/chrom.sizes https://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/mm39.chrom.sizes

# mRNA FASTA (Homework 1)
curl -L -O https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/mrna.fa.gz

# Ice plant TPM (Example 2, Homework 2)
git clone https://github.com/plantgenomicslab/Ice-plant-transcriptome-profiling
```

## Step 0D. Research Project Design Prompt (Advanced)

When starting a new bioinformatics research project, use this prompt to systematically explore analytical directions **before** writing any code. This is especially useful for thesis projects, grant proposals, or novel research questions.

> ## Learning Objective
> Design a comprehensive analysis framework by exploring multiple analytical directions grounded in existing literature.
{: .objectives}

### The Research Project Design Prompt

Copy and customize the following prompt. Replace `[Research Question]` with your specific question.

~~~
Design a bioinformatics analysis framework to address the following research question.

[Research Question]
[YOUR RESEARCH QUESTION HERE - e.g., "How do circadian-regulated genes in CAM plants differ from C3 plants at the regulatory level?"]

Your task is NOT to propose a finalized pipeline, but to explore and structure multiple analytical directions, explicitly grounded in existing literature.

Instructions:

A. Distinct Analytical Directions
1. Propose at least FIVE analysis directions that are clearly differentiated from commonly used or expected approaches in this field.
2. Each direction should represent a distinct analytical framing or inferential perspective, not a minor methodological variation.

B. Structured Evaluation of Each Direction
For EACH proposed analysis direction, provide the following in a clearly labeled structure:

1. Core idea
   - What is the central analytical concept?

2. Why it is interesting
   - What biological or conceptual insight could this reveal that standard analyses typically miss?

3. Relationship to prior work
   - Cite 1–3 representative references (author–year format is sufficient).
   - Explicitly state whether this direction:
     a) Extends existing approaches,
     b) Reinterprets prior findings, or
     c) Challenges an implicit assumption in the literature.
   - Avoid citing review articles unless they are used specifically to define or question a dominant paradigm.
   - If direct primary literature is sparse or absent, explicitly state this limitation and explain how the proposed analysis explores underexamined or emerging conceptual space rather than reiterating established findings.

4. Additional data needs
   - What new or orthogonal data, if any, would strengthen or enable this analysis?

5. Assumptions
   - What biological, evolutionary, or statistical assumptions does this analysis rely on?

6. Analysis difficulty
   - Rate as Low, Medium, or High, and briefly justify the rating.

C. Hypothesis Scope
- Include speculative or not-yet-validated hypotheses where appropriate.
- Do NOT exclude an analysis direction solely because it lacks direct experimental validation.
- Clearly distinguish between evidence-supported claims and conjectural interpretations.

D. Evidence Integration and Conflict Resolution
1. Identify at least three independent axes of evidence across the proposed analyses.
2. Describe how conclusions would be interpreted if these evidence axes yield conflicting or partially inconsistent results.
3. Specify how such inconsistencies would guide follow-up analyses, reframing of hypotheses, or narrowing of scope.

E. Critical Self-Assessment
- Identify where a skeptical reviewer is most likely to push back.
- Discuss risks related to reproducibility, overinterpretation, and literature bias.
- Explicitly distinguish what the data would demonstrate versus what would remain inferential or model-dependent.

Emphasize analytical reasoning, interpretive logic, and literature positioning over tool selection.
~~~

### Example Research Questions

Here are example research questions you can adapt:

| Domain | Example Research Question |
|--------|--------------------------|
| **Transcriptomics** | How do salt stress response genes in halophytes differ from glycophytes at the regulatory network level? |
| **Genomics** | What genomic signatures distinguish drought-tolerant crop varieties from susceptible ones? |
| **Metagenomics** | How does rhizosphere microbiome composition correlate with plant disease resistance? |
| **Comparative Genomics** | What is the evolutionary origin of C4 photosynthesis based on gene family expansion patterns? |
| **Single-cell** | How do cell-type-specific expression patterns change during plant development under stress? |

### When to Use This Prompt

| Situation | Use This Prompt? |
|-----------|-----------------|
| Starting a thesis project | Yes - explore directions before committing |
| Writing a grant proposal | Yes - identify novel angles |
| Class homework assignment | No - use simpler brainstorming prompts |
| Replicating a published analysis | No - follow the original methods |
| Exploring a new dataset | Yes - discover unexpected patterns |

> ## Key Insight
> This prompt forces you to think **beyond the obvious analysis**. Instead of jumping to "run DESeq2," you first ask: "What are five fundamentally different ways to approach this question?"
{: .callout}

> ## Warning: AI Limitations
> AI assistants may:
> - Cite papers that don't exist (hallucination) — always verify references
> - Miss recent publications (knowledge cutoff)
> - Oversimplify domain-specific nuances
>
> Use this prompt as a **starting point for exploration**, not as a definitive literature review.
{: .callout}

---

## Telling AI Assistants About Your Conda Environment

AI coding assistants (Claude, Copilot, Gemini, ChatGPT, Codex) don't know what packages you have installed. **You must tell them explicitly** — otherwise they'll assume arbitrary packages that may not be available in your environment.

### Why This Matters

```
 ┌─────────────────────────────────────────────────────────────────────────┐
 │  Without environment info:              With environment info:          │
 │  ┌─────────────────────────┐           ┌─────────────────────────┐     │
 │  │ AI assumes random       │           │ AI generates code that  │     │
 │  │ packages → code fails   │    vs     │ works in YOUR setup     │     │
 │  │ with ImportError        │           │                         │     │
 │  └─────────────────────────┘           └─────────────────────────┘     │
 └─────────────────────────────────────────────────────────────────────────┘
```

### The Environment Header Pattern

**Always start your prompt with this structure:**

~~~
Write [Python/R] code that runs in the [env-name] conda environment.
The following packages are installed: [list packages].
~~~

### Examples for Different AI Assistants

All AI assistants accept the same prompt format. Here are ready-to-use templates:

> ## For Claude (VS Code / CLI / Web)
> ~~~
> Write Python code that runs in the bch709-python conda environment.
>
> Installed packages:
> - pandas, numpy, matplotlib, seaborn
> - biopython, tqdm
>
> Do NOT use packages outside this list.
> ~~~
{: .solution}

> ## For GitHub Copilot Chat
> ~~~
> # Environment: bch709-python conda environment
> # Available: pandas, numpy, matplotlib, seaborn, biopython, tqdm
> # Task: [your task here]
> ~~~
> Copilot also reads comments in your code, so adding environment comments at the top of your file helps autocomplete suggestions.
{: .solution}

> ## For Gemini Code Assist
> ~~~
> Context: I'm working in a conda environment called bch709-python.
> Installed packages: pandas, numpy, matplotlib, seaborn, biopython, tqdm.
>
> Write code that [your task here].
> ~~~
{: .solution}

> ## For ChatGPT / Codex (Web)
> ~~~
> I'm using a conda environment with the following setup:
>
> Environment name: bch709-python
> Python version: 3.11
> Installed packages:
> - pandas
> - numpy
> - matplotlib
> - seaborn
> - biopython
> - tqdm
>
> Please write code that only uses these packages.
> [Your task here]
> ~~~
{: .solution}

### Persistent Environment Configuration (Per AI Assistant)

Instead of repeating environment info in every prompt, some AI assistants support **persistent configuration files** that automatically provide your environment context.

#### Claude Code: Create CLAUDE.md with `/init`

Claude Code can automatically generate a `CLAUDE.md` file that describes your project. This file is read automatically in every conversation — no need to paste environment info in prompts.

**How to use `/init`:**

```bash
# 1. Navigate to your project directory
$ cd ~/bch709

# 2. Activate your conda environment first
$ conda activate bch709-python

# 3. Launch Claude Code
$ claude

# 4. Inside Claude Code, run the /init command
> /init
```

`/init` will:
1. **Scan your project** — detect files, build systems, and installed packages
2. **Generate a `CLAUDE.md`** — a concise summary of your project environment
3. **Save it to your project root** — Claude reads it automatically from now on

> ## After Running /init
> Review the generated `CLAUDE.md` and edit if needed. A good `CLAUDE.md` is **short and focused** — it should contain only things Claude can't figure out by reading your code.
>
> For each line, ask: *"Would removing this cause Claude to make a mistake?"* If not, cut it.
{: .callout}

#### GitHub Copilot: copilot-instructions.md

Create `.github/copilot-instructions.md` in your repository. Copilot reads this file automatically.

#### Gemini Code Assist: VS Code Settings

Add custom instructions in VS Code: **Settings** → search "Gemini custom instructions" → enter your environment description.

#### ChatGPT / Codex: Custom Instructions

Go to [chatgpt.com](https://chatgpt.com/) → **Profile** → **Customize ChatGPT** → describe your environment in "What would you like ChatGPT to know about you?"

> ## Summary: One-Time Setup per AI
>
> | AI Assistant | Configuration Method | How to Set Up |
> |--------------|---------------------|---------------|
> | **Claude Code** | `CLAUDE.md` in project root | Run `/init` inside Claude Code |
> | **GitHub Copilot** | `.github/copilot-instructions.md` | Create file in repo |
> | **Gemini** | VS Code settings | Add `gemini.codeAssist.customInstructions` |
> | **ChatGPT/Codex** | Custom Instructions (web) | Profile → Customize ChatGPT |
{: .callout}

### Unofficial Configuration Files (Gemini.md, CODEX.md)

Unlike Claude (which officially supports `CLAUDE.md`), Gemini and ChatGPT/Codex don't automatically read configuration files from your project. However, you can still create **unofficial template files** in your project for reference and quick copy-paste.

#### GEMINI.md (Unofficial)

Create `GEMINI.md` in your project root:

```markdown
# Environment Configuration for Gemini

## Python Environment: bch709-python
- Python 3.11
- Installed packages: pandas, numpy, matplotlib, seaborn, biopython, tqdm
- Use ONLY these packages; do not assume other packages are available

## R Environment: bch709-R
- R 4.3
- Installed packages: data.table, ggplot2, pheatmap, viridisLite, scales

## Code Style
- Use pathlib for file paths
- Include error handling for file I/O
- Print progress messages to console

## Project Structure
- Input data: data/
- Output files: results/
- Scripts: scripts/
```

**How to use:** Open `GEMINI.md`, copy the content, and paste at the start of your Gemini conversation:

```
Here is my environment configuration:
[paste GEMINI.md content here]

Now, please write code that [your task]
```

#### CODEX.md (Unofficial)

Create `CODEX.md` in your project root:

```markdown
# Environment Configuration for ChatGPT/Codex

## Python Environment: bch709-python
- Python version: 3.11
- Conda environment name: bch709-python
- Installed packages:
  - pandas (data manipulation)
  - numpy (numerical computing)
  - matplotlib (plotting)
  - seaborn (statistical visualization)
  - biopython (bioinformatics)
  - tqdm (progress bars)

## R Environment: bch709-R
- R version: 4.3
- Conda environment name: bch709-R
- Installed packages:
  - data.table (fast data manipulation)
  - ggplot2 (visualization)
  - pheatmap (heatmaps)
  - viridisLite (color palettes)
  - scales (axis formatting)

## Important Constraints
- Generate code using ONLY the packages listed above
- Do NOT suggest installing additional packages
- Use gzip.open() for .gz files
- Save outputs to results/ directory
```

**How to use:** Copy and paste at the beginning of your ChatGPT/Codex conversation.

### Single-Conversation Chat Message Templates

When you don't want to use persistent Custom Instructions, use these **copy-paste templates** at the start of each conversation:

> ## Template for Gemini (Single Conversation)
> ```
> Before you generate any code, read my environment setup:
>
> I'm working in a conda environment with ONLY these packages installed:
>
> **Python (bch709-python):**
> - Python 3.11
> - pandas, numpy, matplotlib, seaborn, biopython, tqdm
>
> **R (bch709-R):**
> - R 4.3
> - data.table, ggplot2, pheatmap, viridisLite, scales
>
> IMPORTANT: Do NOT use any packages outside this list. If you need a package that isn't listed, tell me before generating code.
>
> Now, here's my task:
> [your task here]
> ```
{: .solution}

> ## Template for ChatGPT/Codex (Single Conversation)
> ```
> I'm a bioinformatics student. Please generate code using ONLY these installed packages:
>
> **Environment: bch709-python (conda)**
> - Python 3.11
> - pandas, numpy, matplotlib, seaborn, biopython, tqdm
>
> **Environment: bch709-R (conda)**
> - R 4.3
> - data.table, ggplot2, pheatmap, viridisLite, scales
>
> Rules:
> 1. Only use packages from the list above
> 2. If you need an unlisted package, ask me first
> 3. Use gzip.open() for .gz files
> 4. Save outputs to results/ directory
>
> Task:
> [your task here]
> ```
{: .solution}

> ## Template for GitHub Copilot Chat (Single Conversation)
> ```
> @workspace I'm using conda environment bch709-python with:
> - pandas, numpy, matplotlib, seaborn, biopython, tqdm
>
> Generate code using only these packages.
>
> Task: [your task here]
> ```
{: .solution}

### Quick Reference: Configuration Methods Summary

| AI Assistant | Persistent Config | Single-Conversation Method |
|--------------|-------------------|---------------------------|
| **Claude** | `CLAUDE.md` (official) | Paste env info at start of prompt |
| **Copilot** | `.github/copilot-instructions.md` | `@workspace` + env comment |
| **Gemini** | VS Code settings | Paste `GEMINI.md` content |
| **ChatGPT/Codex** | Custom Instructions (all conversations) | Paste `CODEX.md` content |

> ## Pro Tip: Create All Config Files at Once
> ```bash
> # Create all configuration files in your project
> mkdir -p .github
>
> # Claude (official)
> cat > CLAUDE.md << 'EOF'
> # Project Environment
> Python: bch709-python (pandas, numpy, matplotlib, seaborn, biopython, tqdm)
> R: bch709-R (data.table, ggplot2, pheatmap, viridisLite, scales)
> EOF
>
> # Copilot (official)
> cat > .github/copilot-instructions.md << 'EOF'
> Python 3.11 environment: pandas, numpy, matplotlib, seaborn, biopython, tqdm
> R 4.3 environment: data.table, ggplot2, pheatmap
> EOF
>
> # Gemini (unofficial - for copy-paste)
> cat > GEMINI.md << 'EOF'
> Context: bch709-python conda env with pandas, numpy, matplotlib, seaborn, biopython, tqdm
> EOF
>
> # ChatGPT/Codex (unofficial - for copy-paste)
> cat > CODEX.md << 'EOF'
> Environment: bch709-python (Python 3.11)
> Packages: pandas, numpy, matplotlib, seaborn, biopython, tqdm
> EOF
>
> echo "Created: CLAUDE.md, .github/copilot-instructions.md, GEMINI.md, CODEX.md"
> ```
{: .callout}

> ## Pro Tip: Save Your Environment Prompt
> Create a text file with your environment description that you can quickly copy-paste:
> ```bash
> $ cat > ~/env_prompt.txt << 'EOF'
> Write Python code that runs in the bch709-python conda environment.
> Installed packages: pandas, numpy, matplotlib, seaborn, biopython, tqdm.
> EOF
> ```
> Then just `cat ~/env_prompt.txt` and paste before each prompt.
{: .callout}

> ## Warning
> **If you don't specify the environment, the AI will assume an arbitrary one** and may generate code that:
> - Uses packages you don't have installed
> - Assumes different package versions
> - Imports modules with different names (e.g., `sklearn` vs `scikit-learn`)
{: .callout}

---

## How to Write Effective Vibe Coding Prompts

Writing a good prompt is like writing a recipe: the more specific your instructions, the better the result. Here's a step-by-step guide to crafting prompts that produce working code on the first try.

### The 5-Part Prompt Structure

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                    Effective Prompt = 5 Essential Parts                     │
│                                                                             │
│  1. Environment  →  "Write Python code in bch709-python conda env"         │
│  2. Input        →  "Read data/file.gz (gzip TSV, columns: a, b, c)"       │
│  3. Task         →  "Compute X using formula Y, filter by Z"               │
│  4. Output       →  "Save to results/out.tsv (cols, decimals, sorting)"    │
│  5. QC/Console   →  "Print top 10 rows, save dropped items to log.txt"     │
└─────────────────────────────────────────────────────────────────────────────┘
```

### Step-by-Step Prompt Construction

#### Step 1: Environment (Who Are You?)

Tell the AI what tools you have.

**Bad:**
```
Write Python code to analyze my data.
```

**Good:**
```
Write Python code that runs in the bch709-python conda environment.
Installed packages: pandas, numpy, matplotlib, seaborn, biopython, tqdm.
```

#### Step 2: Input (What Are You Reading?)

Describe the input file precisely.

**Bad:**
```
Read the GFF file.
```

**Good:**
```
Input: data/MGI.gff3.gz
- Format: GFF3 (9 tab-separated columns), gzip compressed
- Columns: seqid, source, type, start, end, score, strand, phase, attributes
- seqid = chromosome name (e.g., chr1, chr2, chrX)
- type = feature type (gene, exon, mRNA, snRNA, lnc_RNA, etc.)
```

#### Step 3: Task (What Should You Do?)

Define computations with explicit formulas.

**Bad:**
```
Find the most variable genes.
```

**Good:**
```
Task:
1. Compute mean_tpm = row-wise mean of all TPM columns
2. Compute sd_tpm = row-wise standard deviation
3. Compute CV = sd_tpm / mean_tpm
4. Filter: keep only genes where mean_tpm >= 1
5. Select: top 200 genes by CV (descending)
```

#### Step 4: Output (What Files Should You Create?)

Specify exact filenames, formats, columns, and formatting.

**Bad:**
```
Save the results.
```

**Good:**
```
Output: results/cv_top200.tsv
- Format: TSV with header
- Columns: gene_id, mean_tpm, sd_tpm, cv
- Round numeric values to 4 decimal places
- Sort by cv descending
```

#### Step 5: QC/Console (What Should You Print?)

Tell the AI what to display for verification.

**Bad:**
```
Print something.
```

**Good:**
```
Console output:
- Print number of genes that passed the mean >= 1 filter
- Print number of genes that were filtered out
- Print top 10 rows of the result table
- Print "Saved: [filename]" for each output file
```

### Complete Prompt Examples

> ## Example 1: Python GFF3 Analysis (Complete Prompt)
> ```
> Write Python code that runs in the bch709-python conda environment.
> Installed packages: pandas, numpy, matplotlib, seaborn, biopython, tqdm.
>
> **Input:**
> - GFF3 file: data/MGI.gff3.gz (gzip, 9 tab-separated columns)
>   - seqid = chromosome (chr1, chr2, ..., chrX, chrY)
>   - type = feature type (gene, exon, mRNA, snRNA, lnc_RNA, etc.)
> - Chromosome sizes: data/chrom.sizes (TSV: chrom, length_bp)
>
> **Task:**
> 1. Only include seqids that exist in chrom.sizes
> 2. Log seqids NOT in chrom.sizes to a QC file
> 3. Count genes per chromosome (type == "gene")
> 4. Count unique exons per chromosome (unique start, end, strand tuples to prevent isoform overcounting)
> 5. Count snRNA (type == "snRNA")
> 6. Count lncRNA (type == "lnc_RNA" OR "lncRNA")
> 7. Compute density: gene_per_Mb = n_gene / (chrom_length_bp / 1e6)
>
> **Output 1:** results/chr_feature_counts.tsv
> - Columns: chrom, chrom_length_bp, n_gene, n_exon_unique, n_snRNA, n_lncRNA, gene_per_Mb
> - Round densities to 4 decimal places
> - Sort by gene_per_Mb descending
>
> **Output 2:** results/dropped_seqids.txt
> - One seqid per line, sorted alphabetically
>
> **Console:**
> - Print number of dropped seqids and number of dropped feature lines
> - Print top 5 rows of the result table
> ```
{: .solution}

> ## Example 2: R Heatmap Analysis (Complete Prompt)
> ```
> Write R code that runs in the bch709-R conda environment.
> Installed packages: data.table, ggplot2, pheatmap, viridisLite, scales.
>
> **Input:**
> - TPM file: Ice-plant-transcriptome-profiling/iceplant_TPM_DT_ZT.tab.gz
>   - First column: gene_id
>   - Remaining columns: DT_ZT{time}_rep{1-3} (18 sample columns)
>
> **Task:**
> 1. Compute mean_tpm = row-wise mean of all 18 sample columns
> 2. Compute sd_tpm = row-wise standard deviation
> 3. Compute CV = sd_tpm / mean_tpm
> 4. Filter: keep only genes where mean_tpm >= 1
> 5. Select: top 200 genes by CV descending
> 6. Transform for heatmap: log2(TPM + 1)
>
> **Output 1:** results/iceplant_cv_top200.tsv
> - Columns: gene_id, mean_tpm, sd_tpm, cv
> - Round to 4 decimal places
>
> **Output 2:** results/iceplant_cv_top200_heatmap.png
> - Size: 1800 × 1200 pixels, dpi 200
> - Data: log2(TPM + 1) values
> - Rows: gene_id (maintain CV descending order, cluster_rows = FALSE)
> - Columns: original sample order (cluster_cols = FALSE)
> - X-axis labels: rotated 90 degrees
> - Title: "Ice plant log2(TPM+1), CV top200 (mean>=1)"
>
> **Console:**
> - Print top 10 rows of the CV table
> - Print "Saved: [filename]" for each output
> ```
{: .solution}

### Prompt Writing Checklist

Use this checklist before sending your prompt:

> ## Before You Send Your Prompt
> **Environment:**
> - [ ] Specified language (Python/R)
> - [ ] Specified conda environment name
> - [ ] Listed installed packages
>
> **Input:**
> - [ ] Specified filename and path
> - [ ] Specified format (TSV, CSV, GFF3, FASTA, etc.)
> - [ ] Specified if gzip compressed
> - [ ] Described column structure
>
> **Task:**
> - [ ] Defined formulas (CV = sd/mean, etc.)
> - [ ] Specified filter criteria (mean >= 1, etc.)
> - [ ] Explained any deduplication logic
>
> **Output:**
> - [ ] Specified filename and path
> - [ ] Listed column names
> - [ ] Specified decimal places
> - [ ] Specified sorting order
> - [ ] Specified plot dimensions and format (if applicable)
>
> **QC:**
> - [ ] Specified what to print to console
> - [ ] Specified any QC files to save
{: .checklist}

### Common Prompt Mistakes and Fixes

| Mistake | Problem | Fix |
|---------|---------|-----|
| "Analyze the data" | AI doesn't know what analysis | Specify exact computation: "Compute CV = sd/mean" |
| "Save the results" | AI chooses random filename | Specify: "Save to results/output.tsv" |
| "Make a nice plot" | AI chooses arbitrary colors/size | Specify: "1800x1200 px, dpi 200, blue-white-red colors" |
| "Filter low genes" | AI doesn't know threshold | Specify: "Filter: keep genes where mean_tpm >= 1" |
| "Count exons" | AI may double-count isoforms | Specify: "Count unique (start, end, strand) tuples" |

### Iteration Strategy

If the first prompt doesn't work perfectly, follow this pattern:

```
First prompt → AI generates code → Run code → Check output

If errors:
  "The code produced an error: [paste error message]
   Please fix: [describe the issue]"

If wrong output:
  "The output is incorrect. Expected: [describe expected]
   Actual: [describe actual]
   Please modify: [specific change needed]"

If missing feature:
  "The code works but needs these additions:
   1. [addition 1]
   2. [addition 2]"
```

> ## Pro Tip: Start with a Brainstorming Prompt
> Before writing code, ask the AI to help you plan:
> ```
> I need to analyze [data type]. Before writing code:
> 1. What libraries do you recommend?
> 2. What are the key steps in this analysis?
> 3. What edge cases should I handle?
> ```
> Use the AI's response to write a more specific code-generation prompt.
{: .callout}

---

## Part 1: Vibe Coding Examples

---

## Example 1 (Python): Per-Chromosome Feature Counts from MGI GFF3 + QC

### Background

Extract chromosome-level feature counts (genes, exons, snRNAs, lncRNAs) from the MGI GFF3 file and cross-reference against an external chrom.sizes file to verify data integrity.

**GFF3 file structure (9 tab-separated columns):**

```
chr1  MGI  gene  3214482  3671498  .  -  .  ID=MGI:MGI:1918911;Name=Xkr4;biotype=protein_coding
```

| Column | Content |
|--------|---------|
| 1 | chromosome (seqid) |
| 2 | source |
| 3 | feature type (gene, mRNA, exon, snRNA, lnc_RNA, etc.) |
| 4 | start position |
| 5 | end position |
| 6 | score |
| 7 | strand (+/−) |
| 8 | phase |
| 9 | attributes (key=value pairs, semicolon-delimited) |

**chrom.sizes file structure:**

```
chr1    195465000
chr2    182105000
chrX    171031299
```

### Critical Design Decisions

> ## Definitions That Must Appear in Your Prompt
> 1. **Exon overcounting prevention:** Count unique (start, end, strand) intervals per chromosome
> 2. **Chromosome length source:** Use external chrom.sizes file; log excluded seqids to `dropped_seqids.txt`
> 3. **snRNA/lncRNA definition:** Count lines where type is `snRNA`, `lnc_RNA`, or `lncRNA`
{: .checklist}

---

### Stage 1: Vague Prompt

**Prompt:**
> "Extract the gene, exon, snRNA, and lncRNA counts per chromosome from the GFF3 file."

**AI-generated code:**
```python
import gzip, re
from collections import defaultdict

counts = defaultdict(lambda: defaultdict(int))

with gzip.open("data/MGI.gff3.gz", "rt") as f:
    for line in f:
        if line.startswith("#"):
            continue
        fields = line.strip().split("\t")
        if len(fields) < 9:
            continue
        chrom = fields[0]
        ftype = fields[2]
        if ftype in ("gene", "exon", "snRNA", "lncRNA"):
            counts[chrom][ftype] += 1

for chrom in sorted(counts):
    print(chrom, dict(counts[chrom]))
```

> ## Problems with Stage 1
> - No chromosome length reference—scaffolds included with primary chromosomes
> - Exons double-counted across isoforms
> - Misses `lnc_RNA` (standard GFF3 spelling)
> - No QC tracking of dropped seqids
> - No file output
{: .callout}

---

### Stage 2: Specifying Output Format

**Prompt:**
~~~
Write Python code that runs in the bch709-python conda environment.
pandas and numpy are installed.

Input:
- GFF3 gzip: data/MGI.gff3.gz (9 columns, tab-separated)
- Chrom sizes: data/chrom.sizes (TSV: chrom, length_bp)

Task:
- Only include chromosomes present in chrom.sizes
- Save GFF seqids not found in chrom.sizes to results/dropped_seqids.txt
- Count genes per chromosome (type == "gene")
- Count exons per chromosome: to prevent isoform overcounting, define count as number of unique (start, end, strand) intervals
- Count snRNA/lncRNA: lines where type is "snRNA", "lnc_RNA", or "lncRNA"

Output:
- results/chr_feature_counts.tsv
- Columns: chrom, chrom_length_bp, n_gene, n_exon_unique, n_snRNA, n_lncRNA
- Include header
- Print top 5 rows to console
~~~

**AI-generated code:**
```python
from pathlib import Path
import gzip
import pandas as pd
from collections import defaultdict

DATA_GFF  = Path("data/MGI.gff3.gz")
DATA_SIZES = Path("data/chrom.sizes")
OUT_TSV   = Path("results/chr_feature_counts.tsv")
OUT_DROP  = Path("results/dropped_seqids.txt")
OUT_TSV.parent.mkdir(parents=True, exist_ok=True)

# Load chrom sizes
sizes = pd.read_csv(DATA_SIZES, sep="\t", header=None, names=["chrom","length_bp"])
chrom_set = set(sizes["chrom"])
chrom_len = dict(zip(sizes["chrom"], sizes["length_bp"]))

n_gene = defaultdict(int)
exon_intervals = defaultdict(set)
n_snRNA = defaultdict(int)
n_lncRNA = defaultdict(int)
dropped_seqids = set()

SNRNA_TYPES  = {"snRNA"}
LNCRNA_TYPES = {"lnc_RNA", "lncRNA"}

with gzip.open(DATA_GFF, "rt") as f:
    for line in f:
        if not line or line.startswith("#"):
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 9:
            continue
        seqid, _, ftype, start, end, _, strand, _, _ = fields

        if seqid not in chrom_set:
            dropped_seqids.add(seqid)
            continue

        if ftype == "gene":
            n_gene[seqid] += 1
        elif ftype == "exon":
            exon_intervals[seqid].add((int(start), int(end), strand))
        elif ftype in SNRNA_TYPES:
            n_snRNA[seqid] += 1
        elif ftype in LNCRNA_TYPES:
            n_lncRNA[seqid] += 1

rows = []
for chrom in sizes["chrom"]:
    rows.append({
        "chrom": chrom,
        "chrom_length_bp": int(chrom_len[chrom]),
        "n_gene": n_gene.get(chrom, 0),
        "n_exon_unique": len(exon_intervals.get(chrom, set())),
        "n_snRNA": n_snRNA.get(chrom, 0),
        "n_lncRNA": n_lncRNA.get(chrom, 0),
    })

df = pd.DataFrame(rows)
df.to_csv(OUT_TSV, sep="\t", index=False)
OUT_DROP.write_text("\n".join(sorted(dropped_seqids)) + "\n")

print("Saved:", OUT_TSV)
print("Dropped seqids:", len(dropped_seqids))
print(df.head(5).to_string(index=False))
```

**Improvements:** TSV output, unique exon definition, dropped seqid QC

**Still missing:** Density not computed, no sorting, unclear zero-feature handling

---

### Stage 3: Adding Detailed Constraints

**Prompt:**
~~~
Extend the Stage 2 results.

Additional output requirements:
- Add density columns:
  - gene_per_Mb = n_gene / (chrom_length_bp / 1e6)
  - exon_unique_per_Mb, snRNA_per_Mb, lncRNA_per_Mb (same formula)
- Sort by gene_per_Mb descending
- Round to 4 decimal places
- Include all chromosomes from chrom.sizes, filling counts with 0 where no features exist
- Save unique dropped seqids sorted alphabetically to dropped_seqids.txt
- Print to console: number of dropped seqids, number of excluded feature lines, top 5 rows of the result table
~~~

**AI-generated code (final):**
```python
from pathlib import Path
import gzip
import pandas as pd
from collections import defaultdict

DATA_GFF  = Path("data/MGI.gff3.gz")
DATA_SIZES = Path("data/chrom.sizes")
OUT_TSV   = Path("results/chr_feature_counts.tsv")
OUT_DROP  = Path("results/dropped_seqids.txt")
OUT_TSV.parent.mkdir(parents=True, exist_ok=True)

# 1) Chrom sizes
sizes = pd.read_csv(DATA_SIZES, sep="\t", header=None, names=["chrom","length_bp"])
chrom_set = set(sizes["chrom"])
chrom_len = dict(zip(sizes["chrom"], sizes["length_bp"]))

# 2) Counters
n_gene = defaultdict(int)
exon_intervals = defaultdict(set)
n_snRNA = defaultdict(int)
n_lncRNA = defaultdict(int)
dropped_seqids = set()
dropped_lines = 0

SNRNA_TYPES  = {"snRNA"}
LNCRNA_TYPES = {"lnc_RNA", "lncRNA"}

with gzip.open(DATA_GFF, "rt") as f:
    for line in f:
        if not line or line.startswith("#"):
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 9:
            continue
        seqid, _, ftype, start, end, _, strand, _, _ = fields

        if seqid not in chrom_set:
            dropped_seqids.add(seqid)
            dropped_lines += 1
            continue

        if ftype == "gene":
            n_gene[seqid] += 1
        elif ftype == "exon":
            exon_intervals[seqid].add((int(start), int(end), strand))
        elif ftype in SNRNA_TYPES:
            n_snRNA[seqid] += 1
        elif ftype in LNCRNA_TYPES:
            n_lncRNA[seqid] += 1

# 3) Build result table (include all chroms from chrom.sizes; fill 0 where no features)
rows = []
for chrom in sizes["chrom"]:
    L = float(chrom_len[chrom])
    g  = n_gene.get(chrom, 0)
    ex = len(exon_intervals.get(chrom, set()))
    sn = n_snRNA.get(chrom, 0)
    ln = n_lncRNA.get(chrom, 0)
    Mb = L / 1e6 if L > 0 else 1

    rows.append({
        "chrom": chrom,
        "chrom_length_bp": int(L),
        "n_gene": g,
        "n_exon_unique": ex,
        "n_snRNA": sn,
        "n_lncRNA": ln,
        "gene_per_Mb":        round(g  / Mb, 4),
        "exon_unique_per_Mb": round(ex / Mb, 4),
        "snRNA_per_Mb":       round(sn / Mb, 4),
        "lncRNA_per_Mb":      round(ln / Mb, 4),
    })

df = pd.DataFrame(rows).sort_values("gene_per_Mb", ascending=False)
df.to_csv(OUT_TSV, sep="\t", index=False)
OUT_DROP.write_text("\n".join(sorted(dropped_seqids)) + "\n")

print(f"Saved: {OUT_TSV}")
print(f"Saved: {OUT_DROP}")
print(f"Dropped seqids: {len(dropped_seqids)}")
print(f"Dropped feature lines: {dropped_lines}")
print(df.head(5).to_string(index=False))
```

**Sample output:**
```
Saved: results/chr_feature_counts.tsv
Saved: results/dropped_seqids.txt
Dropped seqids: 47
Dropped feature lines: 1823

 chrom  chrom_length_bp  n_gene  n_exon_unique  n_snRNA  n_lncRNA  gene_per_Mb  ...
 chr11      122082543     2847        42156        12        45      23.3224  ...
 chr19       61431566     1892        31245         8        32      30.8024  ...
 chr17       94987271     2456        38912        10        38      25.8563  ...
  chr1      195465000     3215        52341        18        67      16.4488  ...
  chr2      182105000     2987        48562        15        58      16.4027  ...
```

---

### Example 1: Comparison Summary

| Aspect | Stage 1 (Vague) | Stage 2 (Format) | Stage 3 (Detailed) |
|--------|-----------------|------------------|-------------------|
| Chromosome scope | Everything | chrom.sizes only | chrom.sizes + zero-fill |
| Exon definition | Duplicate-counted | Unique interval | Unique interval |
| QC artifact | None | dropped_seqids.txt | Count + line count + file |
| Density | None | None | 4 per_Mb columns |
| Sorting | None | None | gene_per_Mb descending |
| **Reusability** | **Low** | **Medium** | **High (publication-ready)** |

### QC Interpretation Questions

> ## Questions Students Must Answer
> 1. What seqids ended up in `dropped_seqids.txt`? (Alternative contigs? Unplaced scaffolds? Mitochondrial?)
> 2. What fraction of total genes were dropped? Could this affect conclusions?
> 3. If the prompt had NOT specified using chrom.sizes, what errors could have occurred **silently**?
{: .challenge}

> ## Key Message
> `dropped_seqids.txt` is not a failure log—it is a **declaration of what this analysis accepts as data**.
{: .callout}

---

## Example 2 (R): Top 200 Variable Genes from Ice Plant TPM + Heatmap

### Background

Analyze time-course (Zeitgeber Time) TPM expression data from ice plant (*Mesembryanthemum crystallinum*). Extract the **top 200 genes by coefficient of variation (CV)** and visualize their expression patterns.

**Data source:** [Ice-plant-transcriptome-profiling](https://github.com/plantgenomicslab/Ice-plant-transcriptome-profiling)

**Data structure:**

| Column | Description |
|--------|-------------|
| gene_id | Gene identifier (e.g., Mc01g0010) |
| DT_ZT2_rep1 through DT_ZT22_rep3 | 6 time points × 3 replicates = 18 samples |

- DT = Day Treatment, ZT = Zeitgeber Time
- Time points: ZT2, ZT6, ZT10, ZT14, ZT18, ZT22 (4-hour intervals)

### Critical Design Decisions

| Decision | Rationale |
|----------|-----------|
| **mean >= 1 filter** | Genes with near-zero TPM produce infinite/extreme CV values |
| **CV = sd / mean** | Measures relative variability independent of expression level |
| **log2(TPM + 1)** | Reduces influence of extreme values in heatmap |

---

### Stage 1: Vague Prompt

**Prompt:**
> "Pick the 200 most variable genes from the TPM data and draw a heatmap."

**AI-generated code:**
```r
data <- read.table("iceplant_TPM_DT_ZT.tab", header=TRUE, sep="\t", row.names=1)
vars <- apply(data, 1, var)
top200 <- head(sort(vars, decreasing=TRUE), 200)
heatmap(as.matrix(data[names(top200), ]))
```

> ## Problems with Stage 1
> - Uses variance instead of CV → biased toward highly expressed genes
> - No mean >= 1 filter → noise genes in top ranks
> - No file output, no image size/resolution control
> - No log2 transformation → color scale dominated by extremes
{: .callout}

---

### Stage 2: Specifying Output Format

**Prompt:**
~~~
Write R code that runs in the bch709-R conda environment.
data.table and pheatmap are installed.

Input:
- Ice-plant-transcriptome-profiling/iceplant_TPM_DT_ZT.tab.gz (gzip TSV)
- First column: gene_id; remaining columns: numeric TPM

Task:
- Compute mean_tpm and sd_tpm for each gene
- CV = sd_tpm / mean_tpm
- Keep only genes with mean_tpm >= 1
- Select top 200 by CV descending

Output:
- Save to results/iceplant_cv_top200.tsv
- Columns: gene_id, mean_tpm, sd_tpm, cv
- Round to 4 decimal places
- Print top 10 to console
~~~

**AI-generated code:**
```r
library(data.table)

dir.create("results", showWarnings = FALSE, recursive = TRUE)

f <- "Ice-plant-transcriptome-profiling/iceplant_TPM_DT_ZT.tab.gz"
dt <- fread(f)
stopifnot("gene_id" %in% names(dt))

sample_cols <- setdiff(names(dt), "gene_id")
mat <- as.matrix(dt[, ..sample_cols])
mode(mat) <- "numeric"

mean_tpm <- rowMeans(mat, na.rm = TRUE)
sd_tpm   <- apply(mat, 1, sd, na.rm = TRUE)

keep <- mean_tpm >= 1
res <- data.table(
  gene_id  = dt$gene_id[keep],
  mean_tpm = mean_tpm[keep],
  sd_tpm   = sd_tpm[keep]
)
res[, cv := sd_tpm / mean_tpm]
setorder(res, -cv)

top200 <- res[1:min(200, .N)]
top200[, `:=`(
  mean_tpm = round(mean_tpm, 4),
  sd_tpm   = round(sd_tpm, 4),
  cv       = round(cv, 4)
)]

fwrite(top200, "results/iceplant_cv_top200.tsv", sep = "\t")
print(top200[1:min(10, .N)])
cat("Saved: results/iceplant_cv_top200.tsv\n")
```

**Improvements:** Explicit CV definition, mean >= 1 filter, TSV output

**Still missing:** No heatmap generated, no plot specifications

---

### Stage 3: Adding Detailed Constraints

**Prompt:**
~~~
Using the top 200 genes from Stage 2, generate a heatmap.

Heatmap data:
- Values: log2(TPM + 1)
- Rows: gene_id (maintain CV descending order, cluster_rows = FALSE)
- Columns: original sample column order (cluster_cols = FALSE)

Plot specifications:
- Save to: results/iceplant_cv_top200_heatmap.png
- Size: 1800 × 1200 px, dpi 200
- Rotate x-axis labels 90 degrees
- Title: "Ice plant log2(TPM+1), CV top200 (mean>=1)"
- Also print top 10 summary table to console
~~~

**AI-generated code (final):**
```r
library(data.table)
library(pheatmap)

dir.create("results", showWarnings = FALSE, recursive = TRUE)

f <- "Ice-plant-transcriptome-profiling/iceplant_TPM_DT_ZT.tab.gz"
dt <- fread(f)
stopifnot("gene_id" %in% names(dt))

sample_cols <- setdiff(names(dt), "gene_id")
mat <- as.matrix(dt[, ..sample_cols])
mode(mat) <- "numeric"

mean_tpm <- rowMeans(mat, na.rm = TRUE)
sd_tpm   <- apply(mat, 1, sd, na.rm = TRUE)

keep <- mean_tpm >= 1
res <- data.table(
  gene_id  = dt$gene_id[keep],
  mean_tpm = mean_tpm[keep],
  sd_tpm   = sd_tpm[keep]
)
res[, cv := sd_tpm / mean_tpm]
setorder(res, -cv)
top200 <- res[1:min(200, .N)]

# Summary TSV
top200_out <- copy(top200)
top200_out[, `:=`(
  mean_tpm = round(mean_tpm, 4),
  sd_tpm   = round(sd_tpm, 4),
  cv       = round(cv, 4)
)]
fwrite(top200_out, "results/iceplant_cv_top200.tsv", sep = "\t")

# Heatmap matrix
idx <- match(top200$gene_id, dt$gene_id)
submat <- mat[idx, , drop = FALSE]
submat <- log2(submat + 1)
rownames(submat) <- top200$gene_id

# Save (1800×1200 px, dpi 200)
png("results/iceplant_cv_top200_heatmap.png",
    width = 1800, height = 1200, res = 200)
pheatmap(
  submat,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  fontsize_col = 6,
  fontsize_row = 3,
  angle_col = 90,
  main = "Ice plant log2(TPM+1), CV top200 (mean>=1)"
)
dev.off()

print(top200_out[1:min(10, .N)])
cat("Saved: results/iceplant_cv_top200.tsv\n")
cat("Saved: results/iceplant_cv_top200_heatmap.png\n")
```

---

### Example 2: Comparison Summary

| Aspect | Stage 1 (Vague) | Stage 2 (Format) | Stage 3 (Detailed) |
|--------|-----------------|------------------|-------------------|
| Variability metric | Variance | CV (sd/mean) | CV (sd/mean) |
| Filtering | None | mean >= 1 | mean >= 1 |
| Data transformation | None | None | log2(TPM+1) |
| File output | None | TSV | TSV + PNG (size/dpi specified) |
| **Reusability** | **Low** | **Medium** | **High** |

### Interpretation Points

> ## Key Insights
> - Without `mean >= 1` filter, genes with near-zero expression but a single spike dominate top ranks
> - CV captures **relative variability independent of absolute expression** → fair comparison across expression levels
> - `log2(TPM+1)` transformation equalizes heatmap color distribution, making patterns visible
{: .callout}

---

## Part 2: Homework Assignments

---

## Homework 1 (Python): mRNA FASTA Analysis + GC Distribution Graph

### Problem Description

Extract sequence information from the UCSC human mRNA FASTA file (mrna.fa.gz), analyze GC content distribution, and produce a graph and an HTML report page.

> ## Objective
> **Write a single prompt using vibe coding that produces the desired result in one shot.**
> The goal is to get all outputs correct from a single, well-crafted prompt.
{: .objectives}

### Input Data

```bash
curl -L -O https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/mrna.fa.gz
```

**FASTA file structure:**
```
>BC032353 /gb=BC032353 /gi=34783011 /ug=Hs.444613 /len=1873
CGGCAGCGGCTGCGGGGAGGATGGCGGCGACGGCGACTTTTAAATATATTGG...
```

### Expected Output

| Output | Specification |
|--------|---------------|
| `results/mrna_metrics.tsv` | accession, length, gc_content (4 decimals, sorted by gc_content desc, top 20) |
| `results/gc_content_distribution.png` | Histogram + density curve (1600×900 px, dpi 200) |
| `results/gc_content_distribution.svg` | Same graph in SVG format |
| `docs/gc_content_distribution.html` | Title + summary statistics + embedded PNG |

**Plot specifications:**
- Bins: 0.00 to 1.00, step 0.01
- Bar color: light gray; border: dark gray
- Density curve: blue, line width 2.0
- Mean: red dashed vertical line; Median: green dashed vertical line
- Caption showing n, mean, median, sd

### Chart Selection Discussion (In-Class Activity)

~~~
I want to show the distribution of gc_content values.
Compare histogram, density curve, and ECDF — which is most appropriate?

My situation:
- Continuous values between 0 and 1, very large sample size (tens of thousands)
- I want to convey distribution shape, central tendency (mean, median), and tails
- Must be reproducible at publication-quality level

For each option, give pros and cons, then recommend your top pick with justification.
~~~

Reference: [R Graph Gallery — Distribution section](https://r-graph-gallery.com/)

### Grading Criteria

| Criterion | Weight | Description |
|-----------|--------|-------------|
| Prompt quality | 40% | Input format, accession parsing, gzip handling, output specs (filename, columns, decimals, sorting, graph size/color/font) |
| Code correctness | 40% | Correct gzip streaming parse, accurate computation, all 4 files generated |
| Result interpretation | 20% | Explain 3 possible biological/technical reasons for high top-20 GC content values |

---

## Homework 2 (R): Z-Score Clustering of CV Top 200 Genes + Pattern Visualization

### Problem Description

Using the 200 genes from `results/iceplant_cv_top200.tsv` (from Example 2):

1. **Z-score normalize** (row-wise)
2. **Hierarchical clustering** (ward.D2 method, euclidean distance)
3. **Cut tree at k=4** to assign clusters
4. **Save cluster mean expression pattern** line plots
5. **Save cluster assignment** table as TSV

> ## Objective
> Practice writing prompts that **precisely control R visualization output**.
{: .objectives}

### Input Data

```bash
# TSV from Example 2 + original TPM data
# results/iceplant_cv_top200.tsv (gene_id list)
# Ice-plant-transcriptome-profiling/iceplant_TPM_DT_ZT.tab.gz (original TPM)
```

### Expected Output

**1. Clustered Heatmap** — `results/cv_top200_cluster_heatmap.pdf`

| Specification | Value |
|---------------|-------|
| Data | Mean TPM per ZT → row-wise Z-score |
| Z-score | (value − row_mean) / row_sd |
| Rows | gene_id (hierarchical clustering, ward.D2, euclidean) |
| Columns | ZT2, ZT6, ZT10, ZT14, ZT18, ZT22 (chronological, cluster_cols = FALSE) |
| Colors | blue (low) → white (0) → red (high) |
| Annotation | k=4 cutree as color bar |
| Size | 8 × 12 inches |

**2. Cluster Mean Pattern Line Plot** — `results/cluster_patterns.pdf`

| Specification | Value |
|---------------|-------|
| Layout | 2×2 panel (facet_wrap) |
| Y-axis | Mean TPM ± SD |
| X-axis | ZT (2, 6, 10, 14, 18, 22) |
| Panel titles | "Cluster 1 (n=XX genes)" |
| Size | 10 × 8 inches |

**3. Assignment Table** — `results/cluster_assignment.tsv`

| Column | Description |
|--------|-------------|
| gene_id | Gene identifier |
| cluster | 1–4 |
| peak_ZT | ZT with highest mean TPM |
| trough_ZT | ZT with lowest mean TPM |
| amplitude | max − min (2 decimal places) |

Sort by cluster ascending, then amplitude descending within cluster.

### Prompt-Writing Hints

~~~
Analysis procedure:
1. From iceplant_TPM_DT_ZT.tab.gz, extract only the top 200 genes (by gene_id list)
2. Parse ZT time point from sample names (regex: digits after "ZT")
3. Average the 3 replicates per ZT → gene × ZT matrix (200 rows × 6 columns)
4. Z-score normalize: for each row (gene), compute (value - mean) / sd
5. Hierarchical clustering: dist(euclidean) → hclust(ward.D2)
6. cutree(k=4) to assign 4 clusters

For the heatmap:
- Use pheatmap or ComplexHeatmap
- Show cluster assignment as annotation_row color bar
- Columns (ZT) in chronological order (cluster_cols = FALSE)

For the line plot:
- Use original TPM values (no log transformation) for mean ± SD
- ggplot2 facet_wrap(~cluster, ncol=2)
~~~

### Grading Criteria

| Criterion | Weight | Description |
|-----------|--------|-------------|
| Prompt quality | 40% | Z-score definition, clustering method (ward.D2, euclidean), k=4, replicate-to-mean procedure, output specs |
| Code correctness | 40% | Accurate normalization, clustering, cutree, summary statistics, all 3 files generated |
| Result interpretation | 20% | Interpret 4 cluster patterns in context of CAM photosynthesis (2 sentences per cluster) |

---

## Appendix: Effective Vibe Coding Prompt Template

~~~
[Environment]: Write [language] code that runs in the [env name] conda environment.
[Installed packages] are available.

**Input specification:**
- File: [filename] ([format], [delimiter], [gzip?])
- Structure: [column descriptions, special parsing rules]
- Additional inputs: [chrom.sizes, gene lists, or other reference files]

**Analysis conditions:**
- [Filter criteria (e.g., mean >= 1)]
- [Computation method (e.g., CV = sd/mean)]
- [Definitions (e.g., exon count = unique intervals only)]

**Output 1 — Table:**
- Filename: [filename]
- Columns: [list column names]
- Decimal places: [number]
- Sorting: [criterion, direction]
- Filter: [top N]

**Output 2 — Plot:**
- Filename: [filename]
- Size: [px or inches], dpi: [value]
- Colors: [specify explicitly]
- Axes/labels/legend: [specify explicitly]

**Output 3 — QC:**
- [Dropped items filename]
- [Summary information to print to console]
~~~

---

## Input/Output Prompt Checklist

> ## When Specifying Input
> - Filename and path
> - File format (GFF3, TSV, FASTA, etc.)
> - Compression (gzip or not)
> - Delimiter (tab, comma, space)
> - Header presence
> - Data structure (column names, what rows represent)
> - Special structures (e.g., GFF3 attribute parsing rules)
> - External reference files (chrom.sizes, etc.)
{: .checklist}

> ## When Specifying Output
> - File format (TSV, CSV, PDF, PNG, SVG, HTML)
> - Filename
> - Column names and order
> - Decimal places
> - Sorting criterion (ascending/descending)
> - Filter conditions (top N, minimum threshold, etc.)
> - Plot: size, resolution, colors, font, legend position, axis range
> - QC artifacts (dropped items, summary statistics)
{: .checklist}

> ## When Specifying Analysis Definitions
> - Metric definitions (CV = sd/mean, Z-score = (x−mean)/sd)
> - Filter rules (mean >= 1)
> - Deduplication handling (unique intervals, etc.)
> - Transformation methods (log2(TPM+1))
> - Clustering parameters (method, distance metric, k)
{: .checklist}

---

> ## Data Sources
> - **MGI GFF3:** [http://www.informatics.jax.org/downloads/mgigff3/MGI.gff3.gz](http://www.informatics.jax.org/downloads/mgigff3/MGI.gff3.gz)
> - **Human mRNA:** [https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/mrna.fa.gz](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/mrna.fa.gz)
> - **Ice Plant TPM:** [https://github.com/plantgenomicslab/Ice-plant-transcriptome-profiling](https://github.com/plantgenomicslab/Ice-plant-transcriptome-profiling)
{: .callout}
