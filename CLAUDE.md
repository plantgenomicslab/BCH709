# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

BCH709 is a bioinformatics course website built with Jekyll using the Carpentries lesson template. The site is hosted on GitHub Pages at https://plantgenomicslab.github.io/BCH709/.

## Build Commands

```bash
# Install dependencies (requires Ruby and Bundler)
bundle install

# Run local development server
make serve

# Build site without serving
make site

# Run using Docker (alternative)
make docker-serve

# Clean build artifacts
make clean

# Validate lesson Markdown
make lesson-check

# Run unit tests on checking tools
make unittest

# Convert RMarkdown to Markdown (if using .Rmd files in _episodes_rmd/)
make lesson-md
```

## Content Structure

- `_episodes/` - Main lesson content as Markdown files (e.g., Linux commands, RNA-Seq, BLAST, genome assembly)
- `_episodes_rmd/` - RMarkdown source files that get converted to `_episodes/` via `make lesson-md`
- `_episodes_2019/` - Archived lessons from previous years
- `_extras/` - Supplementary materials
- `_includes/` - Jekyll includes (HTML partials)
- `_layouts/` - Jekyll page layouts
- `fig/` - Images for lessons
- `files/` - Downloadable files
- `data/` - Data files for lessons
- `syllabus/` - Course syllabus PDFs

## Key Files

- `index.md` - Course homepage with syllabus, schedule, policies
- `setup.md` - Student setup instructions (Windows/Mac/Linux)
- `_config.yml` - Jekyll configuration (title, URLs, collections)
- `Gemfile` - Ruby dependencies (github-pages gem)

## Quick Commit Workflow

```bash
./commit.sh  # Pull, add all, commit with "BCH709", push
./push.sh    # Set permissions, add all, commit, push
```

## Episode Front Matter

Lesson files in `_episodes/` use YAML front matter. Example structure visible in existing files - follow the same pattern when adding new lessons.
