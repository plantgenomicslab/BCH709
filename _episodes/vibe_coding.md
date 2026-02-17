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
   - [Step 0: Project Setup (GitHub Repo)](#step-0-project-setup)
   - [Step 0A: Brainstorming (Python)](#step-0a-brainstorming-prompt-analysis-1--python)
   - [Step 0A: Brainstorming (R)](#step-0a-brainstorming-prompt-analysis-2--r)
   - [Step 0B: How to Ask AI to Set Up Your Environment](#step-0b-how-to-ask-ai-to-set-up-your-environment)
     - [Environment Setup Prompt (Python)](#step-0b-environment-setup-prompt-analysis-1--python)
     - [Environment Setup Prompt (R)](#step-0b-environment-setup-prompt-analysis-2--r)
   - [Step 0C: Environment Creation Commands](#step-0c-environment-creation-commands)
     - [Set Up AI Configuration Files](#set-up-ai-configuration-files)
   - [Step 0D: Research Project Design (Advanced)](#step-0d-research-project-design-prompt-advanced)
   - [Telling AI About Your Environment](#telling-ai-assistants-about-your-conda-environment)
     - [Persistent Configuration Reference](#persistent-configuration-reference)
   - [How to Write Effective Prompts](#how-to-write-effective-vibe-coding-prompts)
   - [Part 1: Vibe Coding Examples](#part-1-vibe-coding-examples)
     - [Example 1: Python Yeast GFF3 Analysis](#example-1-python-per-chromosome-feature-counts-from-yeast-gff3--qc)
     - [Example 2: R Yeast Stress Heatmap](#example-2-r-top-200-variable-genes-from-yeast-stress-data--heatmap)
   - [Part 2: Homework Assignments](#part-2-homework-assignments)
     - [Homework 1: Python Yeast FASTA Analysis](#homework-1-python-yeast-mrna-fasta-analysis--gc-distribution-graph)
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

## Step 0: Project Setup

> ## Learning Objective
> Before writing any code, create a GitHub repository for your project, then ask the AI about possible approaches and required tools.
{: .objectives}

### Create a GitHub Repository

Start by creating a new repository on GitHub to keep your work organized and version-controlled.

#### 1. Log in to GitHub CLI

First, make sure you are logged in to GitHub from the command line. You only need to do this once.

```bash
$ gh auth login
```

Follow the prompts:
- Select **GitHub.com**
- Select **HTTPS** as the protocol
- Select **Login with a web browser**
- Copy the one-time code, press Enter, and paste it in the browser window that opens
- Authorize the GitHub CLI

To verify you are logged in:

```bash
$ gh auth status
```

```output
github.com
  ✓ Logged in to github.com account your-username
```

#### 2. Create a Project Directory

```bash
# Create a new directory for your vibe coding project
$ mkdir ~/bch709_vibe_coding

# Move into the directory
$ cd ~/bch709_vibe_coding
```

#### 3. Initialize Git and Create the GitHub Repo

```bash
# Initialize a git repository
$ git init

# Create a README file so the repo is not empty
$ echo "# BCH709 Vibe Coding" > README.md

# Stage and make the first commit
$ git add README.md
$ git commit -m "Initial commit"

# Create the repo on GitHub and push
$ gh repo create bch709_vibe_coding --public --source=. --remote=origin --push
```

```output
✓ Created repository your-username/bch709_vibe_coding on GitHub
✓ Added remote origin
✓ Pushed commits to origin/main
```

#### 4. Verify Everything Worked

```bash
# Check that the remote is set up
$ git remote -v
```

```output
origin  https://github.com/your-username/bch709_vibe_coding.git (fetch)
origin  https://github.com/your-username/bch709_vibe_coding.git (push)
```

You can also visit `https://github.com/your-username/bch709_vibe_coding` in your browser to see the repo.

#### 5. Create Project Folders

Set up a directory structure for your data, scripts, and results:

```bash
$ mkdir -p data results scripts
```

```
~/bch709_vibe_coding/
├── README.md
├── data/          ← input files (GFF3, chrom.sizes, expression data, etc.)
├── results/       ← output files (TSV, PNG, PDF)
└── scripts/       ← your Python and R scripts
```

#### 6. Save the Project Structure to GitHub

```bash
# Stage the new folders
$ git add -A

# Commit
$ git commit -m "Add project folder structure"

# Push to GitHub
$ git push
```

> ## Why Start with a GitHub Repo?
> - All your code, data, and results stay in one place
> - You can track changes and revert mistakes with `git log` and `git diff`
> - AI assistants like Claude Code can read your project structure via `CLAUDE.md`
> - You can submit your homework by sharing the repo link
{: .callout}

> ## Saving Your Work
> After making changes, save them to GitHub:
> ```bash
> $ git add -A
> $ git commit -m "Describe what you changed"
> $ git push
> ```
> Do this regularly — after finishing each analysis step or before closing your terminal.
{: .callout}

From here, all commands assume you are working inside `~/bch709_vibe_coding`.

---

## Step 0A. Brainstorming Prompt (Analysis 1 — Python)

Copy and paste this prompt into the AI first. **This is a strategy question, not a code request.**

~~~
I am a beginner student in BCH709.
Before writing any code, brainstorm the approaches and libraries I need for the following analysis.

Analysis (Python, GFF3 analysis):
- Input: saccharomyces_cerevisiae.gff.gz (GFF3, gzip), chrom.sizes (TSV: chrom, length_bp)
- Goal: Count genes, exons (preventing isoform overcounting), tRNAs, and snoRNAs per chromosome; compute density
- Output: TSV table + dropped_seqids.txt (QC artifact)

Requirements:
1) Break the analysis into functional units (input parsing, statistical computation, visualization, file output, QC).
2) For each functional unit, suggest 1–2 candidate libraries.
3) Pick one recommended combination for beginners and explain why.
4) List the exact conda-forge package names for that combination.
5) Provide the conda environment creation command (environment name, Python version, and all packages in one command).
6) (Optional) Provide import verification commands to confirm all packages installed correctly.
~~~

## Step 0A. Brainstorming Prompt (Analysis 2 — R)

Now do the same for the R analysis. Copy and paste this prompt into a **new conversation** (or continue the same one).

~~~
I am a beginner student in BCH709.
Before writing any code, brainstorm the approaches and libraries I need for the following analysis.

Analysis (R, Yeast stress response expression analysis):
- Input: gasch2000.txt (TSV, gene_id + log2 expression ratios across ~170 stress conditions)
- Goal: Select top 200 genes by CV, generate a heatmap of stress response patterns
- Output: TSV + heatmap PNG

Requirements:
1) Break the analysis into functional units (input parsing, statistical computation, visualization, file output, QC).
2) For each functional unit, suggest 1–2 candidate libraries.
3) Pick one recommended combination for beginners and explain why.
4) List the exact conda-forge package names for that combination.
5) Provide the conda environment creation command (environment name, R version, and all packages in one command).
6) (Optional) Provide library() verification commands to confirm all packages installed correctly.
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
- Conda environment name: bch709_vibe_coding
- Pin Python 3.11
- Include import/library verification tests after installation
- Present the commands in copy-paste order so a beginner can just run them one by one
~~~

### Step 0B. Environment Setup Prompt (Analysis 2 — R)

Once brainstorming is complete, use this prompt:

~~~
Using the library combination you just recommended, generate conda environment creation commands.

Conditions:
- Conda environment name: bch709_vibe_coding
- Pin R 4.3
- Install R packages into the SAME environment (bch709_vibe_coding) that already has Python
- Include library() verification tests after installation
- Present the commands in copy-paste order so a beginner can just run them one by one
~~~

> ## One Environment for Everything
> We use a **single conda environment** (`bch709_vibe_coding`) that matches the GitHub repo name. This keeps things simple — one project, one environment, one name.
{: .callout}

## Step 0C. Environment Creation Commands

The AI will generate commands like the ones below. **Copy-paste and run them in your terminal.**

> ## Important
> The commands below are what the AI typically produces. Your results may vary slightly depending on which AI assistant you use — that's fine as long as the verification step passes.
{: .callout}

### Create the Environment (Python + R)

```bash
# Create environment with Python and R together
conda create -n bch709_vibe_coding -y -c conda-forge \
  python=3.11 r-base=4.3 \
  pandas numpy matplotlib seaborn biopython tqdm \
  r-data.table r-ggplot2 r-pheatmap r-viridislite r-scales

# Activate
conda activate bch709_vibe_coding

# Verify Python packages
python -c "import pandas, numpy, matplotlib, Bio; print('Python OK')"

# Verify R packages
R -q -e 'library(data.table); library(ggplot2); library(pheatmap); cat("R OK\n")'
```

> ## Troubleshooting
> If the verification step fails:
> 1. Check that you activated the correct environment: `conda activate bch709_vibe_coding`
> 2. Re-run the install command — sometimes packages fail to download on the first try
> 3. Ask the AI: *"I got this error when verifying: [paste error]. How do I fix it?"*
{: .callout}

### Data Downloads

```bash
# Yeast GFF3 from SGD (Example 1, Homework 1)
curl -L -o data/saccharomyces_cerevisiae.gff.gz http://sgd-archive.yeastgenome.org/curation/chromosomal_feature/saccharomyces_cerevisiae.gff.gz

# Yeast chromosome sizes from UCSC sacCer3 (Example 1)
curl -L -o data/chrom.sizes https://hgdownload.soe.ucsc.edu/goldenPath/sacCer3/bigZips/sacCer3.chrom.sizes

# Yeast mRNA FASTA (Homework 1)
curl -L -o data/mrna.fa.gz https://hgdownload.soe.ucsc.edu/goldenPath/sacCer3/bigZips/mrna.fa.gz

# Yeast genome FASTA (Homework 1, reference)
curl -L -o data/sacCer3.fa.gz https://hgdownload.soe.ucsc.edu/goldenPath/sacCer3/bigZips/sacCer3.fa.gz

# Yeast stress response expression data - Gasch et al. (2000) (Example 2, Homework 2)
curl -L -o data/gasch2000.txt https://www.shackett.org/files/gasch2000.txt
```

### Set Up AI Configuration Files

Now that your environment is installed, let AI create configuration files so every AI assistant knows your setup **before** you ask it to write code.

```bash
# Make sure you're in your project directory with the environment active
cd ~/bch709_vibe_coding
conda activate bch709_vibe_coding

# Launch Claude Code
claude
```

Inside Claude Code, run:

```
> /init
```

Then ask:

~~~
Look at my project and conda environment (bch709_vibe_coding).
Create the following configuration files for my project:

1. .github/copilot-instructions.md — for GitHub Copilot
2. GEMINI.md — for Gemini Code Assist (unofficial, for copy-paste)
3. CODEX.md — for ChatGPT/Codex (unofficial, for copy-paste)

Each file should describe my conda environment, installed packages,
project structure, and coding constraints.
~~~

> ## Why Do This Now?
> From this point on, every time you ask an AI to write code, it will **read these files first** and know exactly what packages you have. No more "ModuleNotFoundError" surprises.
{: .callout}

---

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

### Persistent Configuration Reference

If you followed [Step 0C](#step-0c-environment-creation-commands), your AI configuration files are already set up. Here's what each file does:

> ## Configuration Files Summary
>
> | AI Assistant | File / Setting | What It Does |
> |--------------|----------------|--------------|
> | **Claude Code** | `CLAUDE.md` (created by `/init`) | Auto-read every session — no prompting needed |
> | **GitHub Copilot** | `.github/copilot-instructions.md` | Auto-read for autocomplete and chat |
> | **Gemini** | `GEMINI.md` (copy-paste) or VS Code settings | Paste at start of conversation |
> | **ChatGPT/Codex** | `CODEX.md` (copy-paste) or Custom Instructions | Paste at start of conversation |
{: .callout}

> ## If You Don't Have Config Files
> If you skipped [Set Up AI Configuration Files](#set-up-ai-configuration-files), you must include environment info at the **start of every prompt**:
>
> ~~~
> Write [Python/R] code that runs in the bch709_vibe_coding conda environment.
> Installed packages: [list your packages].
> Do NOT use packages outside this list.
> ~~~
>
> Without this, the AI will assume arbitrary packages and your code will fail with `ImportError`.
{: .callout}

---

## How to Write Effective Vibe Coding Prompts

Writing a good prompt is like writing a recipe: the more specific your instructions, the better the result. Here's a step-by-step guide to crafting prompts that produce working code on the first try.

### The 5-Part Prompt Structure

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                    Effective Prompt = 5 Essential Parts                     │
│                                                                             │
│  1. Tool         →  "Write Python code" (AI reads your config files)       │
│  2. Input        →  "Read data/file.gz (gzip TSV, columns: a, b, c)"       │
│  3. Task         →  "Compute X using formula Y, filter by Z"               │
│  4. Output       →  "Save to results/out.tsv (cols, decimals, sorting)"    │
│  5. QC/Console   →  "Print top 10 rows, save dropped items to log.txt"     │
└─────────────────────────────────────────────────────────────────────────────┘
```

### Step-by-Step Prompt Construction

#### Step 1: Environment (Which Tool Do We Need?)

If you set up [AI configuration files](#set-up-ai-configuration-files), your AI already knows your environment. Just tell it which language to use:

**Bad:**
```
Write Python code to analyze my data.
```

**Good (with config files):**
```
Write Python code for Analysis 1.
```

**Good (without config files):**
```
Write Python code that runs in the bch709_vibe_coding conda environment.
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
Input: data/saccharomyces_cerevisiae.gff.gz
- Format: GFF3 (9 tab-separated columns), gzip compressed
- Columns: seqid, source, type, start, end, score, strand, phase, attributes
- seqid = chromosome name (e.g., chrI, chrII, chrXVI)
- type = feature type (gene, exon, mRNA, tRNA, snoRNA, etc.)
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
> Write Python code that runs in the bch709_vibe_coding conda environment.
> Installed packages: pandas, numpy, matplotlib, seaborn, biopython, tqdm.
>
> **Input:**
> - GFF3 file: data/saccharomyces_cerevisiae.gff.gz (gzip, 9 tab-separated columns)
>   - seqid = chromosome (chrI, chrII, ..., chrXVI, chrM)
>   - type = feature type (gene, exon, mRNA, tRNA, snoRNA, etc.)
> - Chromosome sizes: data/chrom.sizes (TSV: chrom, length_bp)
>
> **Task:**
> 1. Only include seqids that exist in chrom.sizes
> 2. Log seqids NOT in chrom.sizes to a QC file
> 3. Count genes per chromosome (type == "gene")
> 4. Count unique exons per chromosome (unique start, end, strand tuples to prevent isoform overcounting)
> 5. Count tRNA (type == "tRNA")
> 6. Count snoRNA (type == "snoRNA")
> 7. Compute density: gene_per_Mb = n_gene / (chrom_length_bp / 1e6)
>
> **Output 1:** results/chr_feature_counts.tsv
> - Columns: chrom, chrom_length_bp, n_gene, n_exon_unique, n_tRNA, n_snoRNA, gene_per_Mb
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
> Write R code that runs in the bch709_vibe_coding conda environment.
> Installed packages: data.table, ggplot2, pheatmap, viridisLite, scales.
>
> **Input:**
> - Expression file: data/gasch2000.txt (TSV, log2 ratios)
>   - First column: UID (systematic gene name, e.g., YAL001C)
>   - Skip columns: NAME, description, GWEIGHT
>   - Remaining columns: ~170 stress condition columns (log2 expression ratios)
>
> **Task:**
> 1. Parse gene_id from UID column, skip metadata columns
> 2. Compute row-wise sd across all condition columns
> 3. Compute CV = sd / abs(mean) for each gene
> 4. Filter: remove genes with all NA or zero variance
> 5. Select: top 200 genes by CV descending
>
> **Output 1:** results/yeast_stress_cv_top200.tsv
> - Columns: gene_id, mean_expr, sd_expr, cv
> - Round to 4 decimal places
>
> **Output 2:** results/yeast_stress_cv_top200_heatmap.png
> - Size: 1800 × 1200 pixels, dpi 200
> - Data: log2 expression ratios (already log-transformed)
> - Rows: gene_id (maintain CV descending order, cluster_rows = FALSE)
> - Columns: original condition order (cluster_cols = FALSE)
> - X-axis labels: rotated 90 degrees
> - Title: "Yeast stress response, CV top200 (Gasch et al. 2000)"
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

## Example 1 (Python): Per-Chromosome Feature Counts from Yeast GFF3 + QC

> ## Research Question
> **How are genomic features (genes, exons, tRNAs, snoRNAs) distributed across yeast chromosomes, and does feature density correlate with chromosome size?**
{: .objectives}

### Background

Extract chromosome-level feature counts (genes, exons, tRNAs, snoRNAs) from the *Saccharomyces cerevisiae* (yeast) GFF3 file and cross-reference against an external chrom.sizes file to verify data integrity.

**GFF3 file structure (9 tab-separated columns):**

```
chrI  SGD  gene  335  649  .  +  .  ID=YAL069W;Name=YAL069W;gene=YAL069W
```

| Column | Content |
|--------|---------|
| 1 | chromosome (seqid) |
| 2 | source |
| 3 | feature type (gene, mRNA, exon, tRNA, snoRNA, etc.) |
| 4 | start position |
| 5 | end position |
| 6 | score |
| 7 | strand (+/−) |
| 8 | phase |
| 9 | attributes (key=value pairs, semicolon-delimited) |

**chrom.sizes file structure:**

```
chrI     230218
chrII    813184
chrIV    1531933
```

### Critical Design Decisions

> ## Definitions That Must Appear in Your Prompt
> 1. **Exon overcounting prevention:** Count unique (start, end, strand) intervals per chromosome
> 2. **Chromosome length source:** Use external chrom.sizes file; log excluded seqids to `dropped_seqids.txt`
> 3. **tRNA/snoRNA definition:** Count lines where type is `tRNA` or `snoRNA`
{: .checklist}

---

### Stage 1: Vague Prompt

**Prompt:**
> "Extract the gene, exon, tRNA, and snoRNA counts per chromosome from the GFF3 file."

**AI-generated code:**
```python
import gzip, re
from collections import defaultdict

counts = defaultdict(lambda: defaultdict(int))

with gzip.open("data/saccharomyces_cerevisiae.gff.gz", "rt") as f:
    for line in f:
        if line.startswith("#"):
            continue
        fields = line.strip().split("\t")
        if len(fields) < 9:
            continue
        chrom = fields[0]
        ftype = fields[2]
        if ftype in ("gene", "exon", "tRNA", "snoRNA"):
            counts[chrom][ftype] += 1

for chrom in sorted(counts):
    print(chrom, dict(counts[chrom]))
```

> ## Problems with Stage 1
> - No chromosome length reference—all seqids included without filtering
> - Exons double-counted across isoforms
> - No QC tracking of dropped seqids
> - No file output
{: .callout}

---

### Stage 2: Specifying Output Format

**Prompt:**
~~~
Write Python code that runs in the bch709_vibe_coding conda environment.
pandas and numpy are installed.

Input:
- GFF3 gzip: data/saccharomyces_cerevisiae.gff.gz (9 columns, tab-separated)
- Chrom sizes: data/chrom.sizes (TSV: chrom, length_bp)

Task:
- Only include chromosomes present in chrom.sizes
- Save GFF seqids not found in chrom.sizes to results/dropped_seqids.txt
- Count genes per chromosome (type == "gene")
- Count exons per chromosome: to prevent isoform overcounting, define count as number of unique (start, end, strand) intervals
- Count tRNA and snoRNA: lines where type is "tRNA" or "snoRNA"

Output:
- results/chr_feature_counts.tsv
- Columns: chrom, chrom_length_bp, n_gene, n_exon_unique, n_tRNA, n_snoRNA
- Include header
- Print top 5 rows to console
~~~

**AI-generated code:**
```python
from pathlib import Path
import gzip
import pandas as pd
from collections import defaultdict

DATA_GFF  = Path("data/saccharomyces_cerevisiae.gff.gz")
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
n_tRNA = defaultdict(int)
n_snoRNA = defaultdict(int)
dropped_seqids = set()

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
        elif ftype == "tRNA":
            n_tRNA[seqid] += 1
        elif ftype == "snoRNA":
            n_snoRNA[seqid] += 1

rows = []
for chrom in sizes["chrom"]:
    rows.append({
        "chrom": chrom,
        "chrom_length_bp": int(chrom_len[chrom]),
        "n_gene": n_gene.get(chrom, 0),
        "n_exon_unique": len(exon_intervals.get(chrom, set())),
        "n_tRNA": n_tRNA.get(chrom, 0),
        "n_snoRNA": n_snoRNA.get(chrom, 0),
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
  - exon_unique_per_Mb, tRNA_per_Mb, snoRNA_per_Mb (same formula)
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

DATA_GFF  = Path("data/saccharomyces_cerevisiae.gff.gz")
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
n_tRNA = defaultdict(int)
n_snoRNA = defaultdict(int)
dropped_seqids = set()
dropped_lines = 0

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
        elif ftype == "tRNA":
            n_tRNA[seqid] += 1
        elif ftype == "snoRNA":
            n_snoRNA[seqid] += 1

# 3) Build result table (include all chroms from chrom.sizes; fill 0 where no features)
rows = []
for chrom in sizes["chrom"]:
    L = float(chrom_len[chrom])
    g  = n_gene.get(chrom, 0)
    ex = len(exon_intervals.get(chrom, set()))
    tr = n_tRNA.get(chrom, 0)
    sn = n_snoRNA.get(chrom, 0)
    Mb = L / 1e6 if L > 0 else 1

    rows.append({
        "chrom": chrom,
        "chrom_length_bp": int(L),
        "n_gene": g,
        "n_exon_unique": ex,
        "n_tRNA": tr,
        "n_snoRNA": sn,
        "gene_per_Mb":        round(g  / Mb, 4),
        "exon_unique_per_Mb": round(ex / Mb, 4),
        "tRNA_per_Mb":        round(tr / Mb, 4),
        "snoRNA_per_Mb":      round(sn / Mb, 4),
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
Dropped seqids: 1
Dropped feature lines: 42

  chrom  chrom_length_bp  n_gene  n_exon_unique  n_tRNA  n_snoRNA  gene_per_Mb  ...
  chrIII       316620      174          210         10         5     549.5146  ...
  chrI         230218      117          136          4         3     508.2190  ...
  chrVI        270161      136          170          7         4     503.2285  ...
  chrIX        439888      218          275         10         4     495.5724  ...
  chrV         576874      282          356         17         6     488.8408  ...
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
> 1. What seqids ended up in `dropped_seqids.txt`? (Mitochondrial genome? 2-micron plasmid?)
> 2. What fraction of total genes were dropped? Could this affect conclusions?
> 3. If the prompt had NOT specified using chrom.sizes, what errors could have occurred **silently**?
{: .challenge}

> ## Key Message
> `dropped_seqids.txt` is not a failure log—it is a **declaration of what this analysis accepts as data**.
{: .callout}

---

## Example 2 (R): Top 200 Variable Genes from Yeast Stress Data + Heatmap

> ## Research Question
> **Which yeast genes show the most variable expression across environmental stress conditions, and do they reveal a coordinated stress response program?**
{: .objectives}

### Background

Analyze the classic yeast stress response microarray dataset from Gasch et al. (2000). Extract the **top 200 genes by coefficient of variation (CV)** across ~170 environmental stress conditions and visualize their expression patterns.

**Data source:** [Gasch et al. (2000)](https://pubmed.ncbi.nlm.nih.gov/11102521/) — "Genomic expression programs in the response of yeast cells to environmental changes." *Mol Biol Cell* 11(12):4241-4257.

**Data structure:**

| Column | Description |
|--------|-------------|
| UID | Systematic gene name (e.g., YAL001C) |
| NAME | Gene common name (e.g., TFC3) |
| (description) | Functional description |
| GWEIGHT | Gene weight (skip) |
| Remaining columns | ~170 stress conditions (log2 expression ratios) |

- Values are **log2 ratios** (positive = upregulated, negative = downregulated)
- Conditions include: heat shock, oxidative stress, osmotic shock, amino acid starvation, nitrogen depletion, stationary phase

### Critical Design Decisions

| Decision | Rationale |
|----------|-----------|
| **Remove all-NA rows** | Some genes have missing data across conditions |
| **CV = sd / abs(mean)** | Measures relative variability; use absolute mean since values can be negative (log2 ratios) |
| **Data already log-transformed** | No additional log transformation needed for heatmap |

---

### Stage 1: Vague Prompt

**Prompt:**
> "Pick the 200 most variable genes from the stress data and draw a heatmap."

**AI-generated code:**
```r
data <- read.table("data/gasch2000.txt", header=TRUE, sep="\t", row.names=1)
vars <- apply(data, 1, var, na.rm=TRUE)
top200 <- head(sort(vars, decreasing=TRUE), 200)
heatmap(as.matrix(data[names(top200), ]))
```

> ## Problems with Stage 1
> - Uses variance instead of CV → biased toward genes with large absolute changes
> - No NA handling → may crash or produce misleading results
> - No file output, no image size/resolution control
> - Metadata columns (NAME, GWEIGHT) not removed → treated as data
{: .callout}

---

### Stage 2: Specifying Output Format

**Prompt:**
~~~
Write R code that runs in the bch709_vibe_coding conda environment.
data.table and pheatmap are installed.

Input:
- data/gasch2000.txt (TSV, log2 expression ratios)
- First column: UID (gene_id); skip NAME, description, GWEIGHT columns
- Remaining columns: ~170 stress conditions (numeric)

Task:
- Compute mean_expr and sd_expr for each gene across all conditions
- CV = sd_expr / abs(mean_expr)
- Remove genes with all NA or zero variance
- Select top 200 by CV descending

Output:
- Save to results/yeast_stress_cv_top200.tsv
- Columns: gene_id, mean_expr, sd_expr, cv
- Round to 4 decimal places
- Print top 10 to console
~~~

**AI-generated code:**
```r
library(data.table)

dir.create("results", showWarnings = FALSE, recursive = TRUE)

f <- "data/gasch2000.txt"
dt <- fread(f, header = TRUE)

# Extract gene_id and numeric condition columns (skip NAME, description, GWEIGHT)
gene_ids <- dt[[1]]
meta_cols <- c("UID", "NAME", "GWEIGHT")
desc_col <- names(dt)[3]  # description column
skip_cols <- c(meta_cols, desc_col)
cond_cols <- setdiff(names(dt), skip_cols)
mat <- as.matrix(dt[, ..cond_cols])
mode(mat) <- "numeric"

mean_expr <- rowMeans(mat, na.rm = TRUE)
sd_expr   <- apply(mat, 1, sd, na.rm = TRUE)

keep <- !is.na(sd_expr) & sd_expr > 0
res <- data.table(
  gene_id   = gene_ids[keep],
  mean_expr = mean_expr[keep],
  sd_expr   = sd_expr[keep]
)
res[, cv := sd_expr / abs(mean_expr)]
res <- res[is.finite(cv)]
setorder(res, -cv)

top200 <- res[1:min(200, .N)]
top200[, `:=`(
  mean_expr = round(mean_expr, 4),
  sd_expr   = round(sd_expr, 4),
  cv        = round(cv, 4)
)]

fwrite(top200, "results/yeast_stress_cv_top200.tsv", sep = "\t")
print(top200[1:min(10, .N)])
cat("Saved: results/yeast_stress_cv_top200.tsv\n")
```

**Improvements:** Explicit CV definition, NA handling, TSV output

**Still missing:** No heatmap generated, no plot specifications

---

### Stage 3: Adding Detailed Constraints

**Prompt:**
~~~
Using the top 200 genes from Stage 2, generate a heatmap.

Heatmap data:
- Values: log2 expression ratios (already log-transformed in the data)
- Rows: gene_id (maintain CV descending order, cluster_rows = FALSE)
- Columns: original condition order (cluster_cols = FALSE)

Plot specifications:
- Save to: results/yeast_stress_cv_top200_heatmap.png
- Size: 1800 × 1200 px, dpi 200
- Rotate x-axis labels 90 degrees
- Title: "Yeast stress response, CV top200 (Gasch et al. 2000)"
- Also print top 10 summary table to console
~~~

**AI-generated code (final):**
```r
library(data.table)
library(pheatmap)

dir.create("results", showWarnings = FALSE, recursive = TRUE)

f <- "data/gasch2000.txt"
dt <- fread(f, header = TRUE)

# Extract gene_id and numeric condition columns
gene_ids <- dt[[1]]
meta_cols <- c("UID", "NAME", "GWEIGHT")
desc_col <- names(dt)[3]
skip_cols <- c(meta_cols, desc_col)
cond_cols <- setdiff(names(dt), skip_cols)
mat <- as.matrix(dt[, ..cond_cols])
mode(mat) <- "numeric"

mean_expr <- rowMeans(mat, na.rm = TRUE)
sd_expr   <- apply(mat, 1, sd, na.rm = TRUE)

keep <- !is.na(sd_expr) & sd_expr > 0
res <- data.table(
  gene_id   = gene_ids[keep],
  mean_expr = mean_expr[keep],
  sd_expr   = sd_expr[keep]
)
res[, cv := sd_expr / abs(mean_expr)]
res <- res[is.finite(cv)]
setorder(res, -cv)
top200 <- res[1:min(200, .N)]

# Summary TSV
top200_out <- copy(top200)
top200_out[, `:=`(
  mean_expr = round(mean_expr, 4),
  sd_expr   = round(sd_expr, 4),
  cv        = round(cv, 4)
)]
fwrite(top200_out, "results/yeast_stress_cv_top200.tsv", sep = "\t")

# Heatmap matrix
idx <- match(top200$gene_id, gene_ids)
submat <- mat[idx, , drop = FALSE]
rownames(submat) <- top200$gene_id

# Save (1800×1200 px, dpi 200)
png("results/yeast_stress_cv_top200_heatmap.png",
    width = 1800, height = 1200, res = 200)
pheatmap(
  submat,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  fontsize_col = 4,
  fontsize_row = 3,
  angle_col = 90,
  main = "Yeast stress response, CV top200 (Gasch et al. 2000)"
)
dev.off()

print(top200_out[1:min(10, .N)])
cat("Saved: results/yeast_stress_cv_top200.tsv\n")
cat("Saved: results/yeast_stress_cv_top200_heatmap.png\n")
```

---

### Example 2: Comparison Summary

| Aspect | Stage 1 (Vague) | Stage 2 (Format) | Stage 3 (Detailed) |
|--------|-----------------|------------------|-------------------|
| Variability metric | Variance | CV (sd/abs(mean)) | CV (sd/abs(mean)) |
| Filtering | None | Remove NA/zero-variance | Remove NA/zero-variance |
| Data transformation | None | None | Already log2 (no extra transform) |
| File output | None | TSV | TSV + PNG (size/dpi specified) |
| **Reusability** | **Low** | **Medium** | **High** |

### Interpretation Points

> ## Key Insights
> - Removing NA and zero-variance genes prevents infinite/undefined CV values from dominating results
> - CV captures **relative variability independent of absolute expression** → fair comparison across expression levels
> - The Gasch 2000 data is already log2-transformed, so no additional transformation is needed for the heatmap
> - Genes with high CV across stress conditions are likely part of the **Environmental Stress Response (ESR)** — a conserved transcriptional program in yeast
{: .callout}

---

## Part 2: Homework Assignments

---

## Homework 1 (Python): Yeast mRNA FASTA Analysis + GC Distribution Graph

> ## Research Question
> **What is the GC content distribution of yeast mRNA sequences, and are there distinct GC-content subpopulations?**
{: .objectives}

### Problem Description

Extract sequence information from the UCSC yeast (*Saccharomyces cerevisiae*) mRNA FASTA file (mrna.fa.gz), analyze GC content distribution, and produce a graph and an HTML report page.

> ## Objective
> **Write a single prompt using vibe coding that produces the desired result in one shot.**
> The goal is to get all outputs correct from a single, well-crafted prompt.
{: .objectives}

### Input Data

```bash
# Yeast mRNA FASTA
curl -L -o data/mrna.fa.gz https://hgdownload.soe.ucsc.edu/goldenPath/sacCer3/bigZips/mrna.fa.gz

# Yeast genome FASTA (for reference)
curl -L -o data/sacCer3.fa.gz https://hgdownload.soe.ucsc.edu/goldenPath/sacCer3/bigZips/sacCer3.fa.gz
```

**FASTA file structure:**
```
>BC001547 /gb=BC001547 /gi=12654078 /ug=Sc.3456 /len=1254
ATGTCTGCTCCAGCTAGCAGTGAAACTTTATTCAGAAACTGCTTAG...
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
| Result interpretation | 20% | Explain 3 possible biological/technical reasons for high top-20 GC content values in yeast |

---

## Homework 2 (R): Z-Score Clustering of CV Top 200 Genes + Pattern Visualization

> ## Research Question
> **Do the top 200 most variable yeast stress-response genes cluster into distinct expression patterns, and what biological processes characterize each cluster?**
{: .objectives}

### Problem Description

Using the 200 genes from `results/yeast_stress_cv_top200.tsv` (from Example 2):

1. **Z-score normalize** (row-wise)
2. **Hierarchical clustering** (ward.D2 method, euclidean distance)
3. **Cut tree at k=4** to assign clusters
4. **Save cluster mean expression pattern** box plots per stress category
5. **Save cluster assignment** table as TSV

> ## Objective
> Practice writing prompts that **precisely control R visualization output**.
{: .objectives}

### Input Data

```bash
# TSV from Example 2 + original expression data
# results/yeast_stress_cv_top200.tsv (gene_id list)
# data/gasch2000.txt (original log2 expression ratios)
```

### Expected Output

**1. Clustered Heatmap** — `results/cv_top200_cluster_heatmap.pdf`

| Specification | Value |
|---------------|-------|
| Data | Log2 expression ratios → row-wise Z-score |
| Z-score | (value − row_mean) / row_sd |
| Rows | gene_id (hierarchical clustering, ward.D2, euclidean) |
| Columns | Stress conditions (original order, cluster_cols = FALSE) |
| Colors | blue (low) → white (0) → red (high) |
| Annotation | k=4 cutree as color bar |
| Size | 8 × 12 inches |

**2. Cluster Mean Pattern Box Plot** — `results/cluster_patterns.pdf`

| Specification | Value |
|---------------|-------|
| Layout | 2×2 panel (facet_wrap) |
| Y-axis | Mean log2 expression ratio ± SD |
| X-axis | Stress categories (heat, oxidative, osmotic, nutrient) |
| Panel titles | "Cluster 1 (n=XX genes)" |
| Size | 10 × 8 inches |

**3. Assignment Table** — `results/cluster_assignment.tsv`

| Column | Description |
|--------|-------------|
| gene_id | Gene identifier (e.g., YAL001C) |
| cluster | 1–4 |
| peak_condition | Condition with highest expression |
| trough_condition | Condition with lowest expression |
| amplitude | max − min (2 decimal places) |

Sort by cluster ascending, then amplitude descending within cluster.

### Prompt-Writing Hints

~~~
Analysis procedure:
1. From gasch2000.txt, extract only the top 200 genes (by gene_id list from results/yeast_stress_cv_top200.tsv)
2. Skip metadata columns (NAME, description, GWEIGHT), keep only numeric condition columns
3. Z-score normalize: for each row (gene), compute (value - mean) / sd
4. Hierarchical clustering: dist(euclidean) → hclust(ward.D2)
5. cutree(k=4) to assign 4 clusters

For the heatmap:
- Use pheatmap or ComplexHeatmap
- Show cluster assignment as annotation_row color bar
- Columns in original condition order (cluster_cols = FALSE)

For the box/line plot:
- Group conditions by stress type (heat, oxidative, osmotic, nutrient depletion)
- Use original log2 ratios for mean ± SD
- ggplot2 facet_wrap(~cluster, ncol=2)
~~~

### Grading Criteria

| Criterion | Weight | Description |
|-----------|--------|-------------|
| Prompt quality | 40% | Z-score definition, clustering method (ward.D2, euclidean), k=4, metadata column handling, output specs |
| Code correctness | 40% | Accurate normalization, clustering, cutree, summary statistics, all 3 files generated |
| Result interpretation | 20% | Interpret 4 cluster patterns in context of yeast Environmental Stress Response (ESR) (2 sentences per cluster) |

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
> - **Yeast GFF3 (SGD):** [http://sgd-archive.yeastgenome.org/curation/chromosomal_feature/saccharomyces_cerevisiae.gff.gz](http://sgd-archive.yeastgenome.org/curation/chromosomal_feature/saccharomyces_cerevisiae.gff.gz)
> - **Yeast chrom.sizes (UCSC sacCer3):** [https://hgdownload.soe.ucsc.edu/goldenPath/sacCer3/bigZips/sacCer3.chrom.sizes](https://hgdownload.soe.ucsc.edu/goldenPath/sacCer3/bigZips/sacCer3.chrom.sizes)
> - **Yeast mRNA FASTA (UCSC sacCer3):** [https://hgdownload.soe.ucsc.edu/goldenPath/sacCer3/bigZips/mrna.fa.gz](https://hgdownload.soe.ucsc.edu/goldenPath/sacCer3/bigZips/mrna.fa.gz)
> - **Yeast genome FASTA (UCSC sacCer3):** [https://hgdownload.soe.ucsc.edu/goldenPath/sacCer3/bigZips/sacCer3.fa.gz](https://hgdownload.soe.ucsc.edu/goldenPath/sacCer3/bigZips/sacCer3.fa.gz)
> - **Yeast stress response (Gasch 2000):** [https://www.shackett.org/files/gasch2000.txt](https://www.shackett.org/files/gasch2000.txt) — Original paper: [Gasch et al. (2000) Mol Biol Cell](https://pubmed.ncbi.nlm.nih.gov/11102521/)
{: .callout}
