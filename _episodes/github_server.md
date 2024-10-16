## GitHub

### Introduction to Git and GitHub
Git is a distributed version control system that allows you to track changes in your code and collaborate with others. GitHub is a cloud-based platform that hosts Git repositories, making it easier to share and collaborate on projects. Version control using Git is the most reasonable way to keep track of changes in code, manuscripts, presentations, and data analysis projects.

### Why Use GitHub in Bioinformatics?
- **Version Control:** Track changes in scripts, pipelines, and documentation.
- **Collaboration:** Work with peers and share your work with the community.
- **Reproducibility:** Maintain a history of your analyses for reproducibility.
- **Integration:** Connect with other tools and platforms used in bioinformatics.

### Why Version Control?
Version control is essential in creating any project that takes longer than 5 minutes to complete. Even if your memory is longer than 5 minutes, next month you are not likely to be able to retrace your steps.  
![github-workflow](./fig/git_overview.png)

![github](./fig/github_dri.png)

## 1. Setting Up Git and GitHub
### Mac
```bash
brew install git
git --version
```

### Linux (WSL)
```bash
sudo apt-get update
sudo apt-get install git
git --version
```

## 2. Configure Git
Open your terminal or command prompt and set your Git username and email:
```bash
git config --global user.name "Your Name"
git config --global user.email "youremail@example.com"
```

### Generate SSH Keys 
SSH keys allow secure communication with GitHub without entering your password each time.

#### Generate an SSH key:
```bash
ssh-keygen -t ed25519 -C "youremail@example.com"
```
Press Enter to accept the default file location and enter a passphrase if desired.

#### Start the SSH agent and add your key:
```bash
eval "$(ssh-agent -s)"
ssh-add ~/.ssh/id_ed25519
```

#### Add the SSH key to your GitHub account:
Copy the SSH key to your clipboard:
```bash
cat ~/.ssh/id_ed25519.pub
```
Go to GitHub > Settings > SSH and GPG keys > New SSH key. Paste the key and save.

## 3. Basic Git Commands
### Understanding basic Git commands is essential for managing your projects.

#### Initialize a Repository
```bash
git init
```

#### Clone a Repository
```bash
git clone git@github.com:wyim-pgl/test_files.git
```

#### Check Repository Status
```bash
git status
```

#### Add Changes
```bash
git add filename
```

#### Commit Change
```bash
git commit -m "Your commit message"
```

#### Push Changes to GitHub
```bash
git push origin main
```

### Pull Updates from GitHub
```bash
git pull origin main
```

## 4. Working with Repositories
Repositories (repos) are where your project files and version history are stored.

### 4.1 Creating a New Repository on GitHub
- Log in to GitHub and click the **+** icon in the top-right corner.
- Select **New repository**.
- Enter a repository name, description (optional), and choose to make it **Public** or **Private**.
- (Optional) Initialize with a README, .gitignore, or license.
- Click **Create repository**.

### 4.2 Forking a Repository
If you want to contribute to someone else's project, you can create a fork:
- Navigate to the repository you want to fork.
- Click the **Fork** button in the top-right corner.
- This creates a copy of the repository in your GitHub account.

### 4.3 Cloning Your Repository Locally
Clone your repository (either your own or a forked one) to work on it locally:
```bash
git clone git@github.com:your-username/repository-name.git
```

### 4.4 Creating Branches
Branches allow you to work on features or experiments without affecting the main codebase.
```bash
git checkout -b feature-branch
```

### 4.5 Committing Changes
After making changes, you need to commit them:
```bash
git add .
git commit -m "Describe your changes here"
```

### 4.6 Pushing Changes to GitHub
To upload your changes to the remote repository:
```bash
git push origin feature-branch
```

### 4.7 Pulling Changes from GitHub
To update your local repository with changes from the remote repository:
```bash
git pull origin main
```
This ensures that your local branch is up-to-date with the latest changes from the main branch.

### 4.8 Merging Branches
After completing work on a branch, merge it back into the main branch.

Switch to the main branch:
```bash
git checkout main
```

Merge the feature branch:
```bash
git merge feature-branch
```

### 4.9 Deleting Branches
After merging, you can delete the branch:
```bash
git branch -d feature-branch
```

## 5. Collaborating with Others
GitHub facilitates collaboration through features like pull requests, issues, and code reviews.

### 5.1 Forking a Repository
- Navigate to the repository you want to fork.
- Click the **Fork** button in the top-right corner.
- This creates a copy of the repository under your GitHub account.

### 5.2 Creating a Pull Request
After making changes in your forked repository:
- Push your changes to a branch in your fork.
- Go to the original repository and click **Pull requests** > **New pull request**.
- Select your fork and branch, then submit the pull request.

### 5.3 Managing Issues
Issues are used to track tasks, enhancements, and bugs.
- Go to the **Issues** tab in your repository.
- Click **New issue**.
- Fill in the title and description, then submit.

### 5.4 Code Reviews
Collaborators can review pull requests, suggest changes, and approve merges.
- Use comments to provide feedback.
- Request changes if necessary.
- Approve and merge once the code meets standards.

## 6. Best Practices for Bioinformatics Projects
### 6.1 Organize Your Repository
Structure your repository to make it easy to navigate. A typical bioinformatics repo might include:
- `README.md`: Project overview and instructions.
- `data/`: Raw and processed datasets.
- `scripts/`: Analysis scripts (e.g., Python, R).
- `results/`: Output files and visualizations.
- `docs/`: Additional documentation.

### 6.2 Write Clear Commit Messages
Use descriptive commit messages to explain what changes were made and why.

**Good Example:**
```bash
Add script for RNA-seq data normalization
```

**Bad Example:**
```bash
Update stuff
```

### 6.3 Use .gitignore Files
Exclude unnecessary files (e.g., large datasets, temporary files) from your repository by creating a `.gitignore` file.

**Example `.gitignore`:**
```bash
# Ignore data files
/data/raw/
/data/processed/

# Ignore temporary files
*.tmp
*.log
```

### 6.4 Document Your Work
Maintain clear documentation to help others understand and reproduce your analyses.
- Update the `README.md` with project details.
- Comment your scripts thoroughly.
- Use Markdown for well-formatted documentation.

### 6.5 Version Control for Data
While Git handles code effectively, managing large datasets can be challenging. Consider using [Git LFS (Large File Storage)](https://git-lfs.github.com/) for large files.

**Installing Git LFS:**
```bash
git lfs install
```

**Tracking a Large File:**
```bash
git lfs track "*.csv"
```

## HISAT2 installation from GitHub
```bash
https://github.com/DaehwanKimLab/hisat2
```
![github-workflow](./fig/hisat2_git.png)
![github-workflow](./fig/hisat2_git2.png)

### Download
```bash
cd ~/bch709/bin
git clone https://github.com/DaehwanKimLab/hisat2.git
cd hisat2
ls -algh
less MANUAL
```
![github-workflow](./fig/hisat2_git3.png)

### Compile
```bash
make -j <YOUR CPU>
```

## HISAT2 installation from binary
![github-workflow](./fig/hisat2_binary.png)

### BWA
> ## BWA Source Code Installation
> If you prefer to install from source, follow the instructions below:
> ```bash
> cd ~/bch709/bin
> curl -OL http://sourceforge.net/projects/bio-bwa/files/bwa-0.7.17.tar.bz2
> tar xvf bwa-0.7.17.tar.bz2
> cd bwa-0.7.17
> make
> ```

### Please check README.md
```bash
less <FILENAME>
```
#### Tell me your error

> ## How to solve it?
> ```bash
> sudo apt install zlib1g-dev
> brew install zlib
> ```

**Test your installation by running:**
```bash
./bwa
```

### BWA GitHub Installation
***Search BWA on Google***

### How to run BWA?
```bash
bwa index <YOUR_GENOME_SEQUENCE>
bwa mem  <YOUR_GENOME_SEQUENCE> <SEQUENCING_READS>
```

### Reference:
- Conda documentation: https://docs.conda.io/en/latest/
- Conda-forge: https://conda-forge.github.io/
- BioConda: https://bioconda.github.io/


