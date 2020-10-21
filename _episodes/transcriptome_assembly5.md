---
layout: page
title: 9_Transcriptome Assembly
published: true
---

### Conda env
```
conda activate rnaseq
```


#### File prepare
```
mkdir -p /data/gpfs/assoc/bch709-1/<YOURID>/rnaseq_assembly/trim

cd /data/gpfs/assoc/bch709-1/<YOURID>/rnaseq_assembly/trim

cp /data/gpfs/assoc/bch709-1/Course_material/2020/RNASeq_trimmed_fastq/*.gz .

cd ../

ls


```


### Youtube Video for Slurm
https://www.youtube.com/watch?v=U42qlYkzP9k

https://www.youtube.com/watch?v=LRJMQO7Ercw

### More Slurm command
```
# Show the overall status of each partition
sinfo

# Submit a job
sbatch <JOBFILE>

# See the entire job queue
squeue

# See only jobs for a given user
squeue -u <username>

# Count number of running / in queue jobs
squeue -u <username> | wc -l

# Get estimated start times for your jobs (when Sherlock is busy)
squeue --start -u <username>

# Show the status of a currently running job
sstat -j <jobID>

# Show the final status of a finished job
sacct -j <jobID>

# Kill a job with ID $PID
scancel <jobID>

# Kill ALL jobs for a user
scancel -u <username>
```


***PLEASE CHECK YOUR MULTIQC after job running***


## Running Trinity

Trinity is run via the script: 'Trinity' found in the base installation directory.

Usage info is as follows:

     ###############################################################################
     #
     #     ______  ____   ____  ____   ____  ______  __ __
     #    |      ||    \ |    ||    \ |    ||      ||  |  |
     #    |      ||  D  ) |  | |  _  | |  | |      ||  |  |
     #    |_|  |_||    /  |  | |  |  | |  | |_|  |_||  ~  |
     #      |  |  |    \  |  | |  |  | |  |   |  |  |___, |
     #      |  |  |  .  \ |  | |  |  | |  |   |  |  |     |
     #      |__|  |__|\_||____||__|__||____|  |__|  |____/
     #
     ###############################################################################
     #
     # Required:
     #
     #  --seqType <string>      :type of reads: ( fa, or fq )
     #
     #  --max_memory <string>      :suggested max memory to use by Trinity where limiting can be enabled. (jellyfish, sorting, etc)
     #                            provided in Gb of RAM, ie.  '--max_memory 10G'
     #
     #  If paired reads:
     #      --left  <string>    :left reads, one or more file names (separated by commas, not spaces)
     #      --right <string>    :right reads, one or more file names (separated by commas, not spaces)
     #
     #  Or, if unpaired reads:
     #      --single <string>   :single reads, one or more file names, comma-delimited (note, if single file contains pairs, can use flag: --run_as_paired )
     #
     #  Or,
     #      --samples_file <string>         tab-delimited text file indicating biological replicate relationships.
     #                                   ex.
     #                                        cond_A    cond_A_rep1    A_rep1_left.fq    A_rep1_right.fq
     #                                        cond_A    cond_A_rep2    A_rep2_left.fq    A_rep2_right.fq
     #                                        cond_B    cond_B_rep1    B_rep1_left.fq    B_rep1_right.fq
     #                                        cond_B    cond_B_rep2    B_rep2_left.fq    B_rep2_right.fq
     #
     #                      # if single-end instead of paired-end, then leave the 4th column above empty.
     #
     ####################################
     ##  Misc:  #########################
     #
     #  --SS_lib_type <string>          :Strand-specific RNA-Seq read orientation.
     #                                   if paired: RF or FR,
     #                                   if single: F or R.   (dUTP method = RF)
     #                                   See web documentation.
     #
     #  --CPU <int>                     :number of CPUs to use, default: 2
     #  --min_contig_length <int>       :minimum assembled contig length to report
     #                                   (def=200)
     #
     #  --long_reads <string>           :fasta file containing error-corrected or circular consensus (CCS) pac bio reads
     #
     #  --genome_guided_bam <string>    :genome guided mode, provide path to coordinate-sorted bam file.
     #                                   (see genome-guided param section under --show_full_usage_info)
     #
     #  --jaccard_clip                  :option, set if you have paired reads and
     #                                   you expect high gene density with UTR
     #                                   overlap (use FASTQ input file format
     #                                   for reads).
     #                                   (note: jaccard_clip is an expensive
     #                                   operation, so avoid using it unless
     #                                   necessary due to finding excessive fusion
     #                                   transcripts w/o it.)
     #
     #  --trimmomatic                   :run Trimmomatic to quality trim reads
     #                                        see '--quality_trimming_params' under full usage info for tailored settings.
     #
     #
     #  --no_normalize_reads            :Do *not* run in silico normalization of reads. Defaults to max. read coverage of 50.
     #                                       see '--normalize_max_read_cov' under full usage info for tailored settings.
     #                                       (note, as of Sept 21, 2016, normalization is on by default)
     #
     #
     #
     #  --output <string>               :name of directory for output (will be
     #                                   created if it doesn't already exist)
     #                                   default( your current working directory: "/Users/bhaas/GITHUB/trinityrnaseq/trinity_out_dir"
     #                                    note: must include 'trinity' in the name as a safety precaution! )
     #
     #  --full_cleanup                  :only retain the Trinity fasta file, rename as ${output_dir}.Trinity.fasta
     #
     #  --cite                          :show the Trinity literature citation
     #
     #  --version                       :reports Trinity version (BLEEDING_EDGE) and exits.
     #
     #  --show_full_usage_info          :show the many many more options available for running Trinity (expert usage).
     #
     #
     ###############################################################################
     #
     #  *Note, a typical Trinity command might be:
     #
     #        Trinity --seqType fq --max_memory 50G --left reads_1.fq  --right reads_2.fq --CPU 6
     #
     #
     #    and for Genome-guided Trinity:
     #
     #        Trinity --genome_guided_bam rnaseq_alignments.csorted.bam --max_memory 50G
     #                --genome_guided_max_intron 10000 --CPU 6
     #
     #     see: /Users/bhaas/GITHUB/trinityrnaseq/sample_data/test_Trinity_Assembly/
     #          for sample data and 'runMe.sh' for example Trinity execution
     #
     #     For more details, visit: http://trinityrnaseq.github.io
     #
     ###############################################################################


<a name='strand_specific_assembly'>
Trinity performs best with strand-specific data, in which case sense and antisense transcripts can be resolved.  For protocols on strand-specific RNA-Seq, see: [Borodina T, Adjaye J, Sultan M. A strand-specific library preparation protocol for RNA sequencing. Methods Enzymol. 2011;500:79-98. PubMed PMID: 21943893](http://www.ncbi.nlm.nih.gov/pubmed/21943893).


If you have strand-specific data, specify the library type.  There are four library types:

- Paired reads:
    * *RF*: first read (/1) of fragment pair is sequenced as anti-sense (reverse(*R*)), and second read (/2) is in the sense strand (forward(*F*)); typical of the dUTP/UDG sequencing method.
    * *FR*: first read (/1) of fragment pair is sequenced as sense (forward), and second read (/2) is in the antisense strand (reverse)

- Unpaired (single) reads:
    * *F*: the single read is in the sense (forward) orientation
    * *R*: the single read is in the antisense (reverse) orientation

By setting the *--SS_lib_type* parameter to one of the above, you are indicating that the reads are strand-specific.  By default, reads are treated as not strand-specific.

![strand specific specification](https://raw.githubusercontent.com/wiki/trinityrnaseq/trinityrnaseq/images/strand_specificity.jpg)

Other important considerations:

- Whether you use Fastq or Fasta formatted input files, be sure to keep the reads oriented as they are reported by Illumina, if the data are strand-specific. This is because, Trinity will properly orient the sequences according to the specified library type.  If the data are not strand-specific, no worries because the reads will be parsed in both orientations.

- If you have both paired and unpaired data, and the data are NOT strand-specific, you can combine the unpaired data with the left reads of the paired fragments.  Be sure that the unpaired reads have a /1 as a suffix to the accession value similarly to the left fragment reads.  The right fragment reads should all have /2 as the accession suffix.  Then, run Trinity using the --left and --right parameters as if all the data were paired.

- If you have multiple paired-end library fragment sizes, set the '--group_pairs_distance' according to the larger insert library.  Pairings that exceed that distance will be treated as if they were unpaired by the Butterfly process.

- by setting the '--CPU option', you are indicating the maximum number of threads to be used by processes within Trinity. Note that Inchworm alone will be internally capped at 6 threads, since performance will not improve for this step beyond that setting)


<a name="typical_trinity_command_line"></a>
## Typical Trinity Command Line

A typical Trinity command for assembling non-strand-specific RNA-seq data would be like so, running the entire process on a single high-memory server (aim for \~1G RAM per \~1M \~76 base Illumina paired reads, but often *much* less memory is required):

Run Trinity like so:

     Trinity --seqType fq --max_memory 50G \
             --left reads_1.fq.gz  --right reads_2.fq.gz --CPU 6

If you have multiple sets of fastq files, such as corresponding to multiple tissue types or conditions, etc., you can indicate them to Trinity like so:

     Trinity --seqType fq --max_memory 50G  \
             --left condA_1.fq.gz,condB_1.fq.gz,condC_1.fq.gz \
             --right condA_2.fq.gz,condB_2.fq.gz,condC_2.fq.gz \
             --CPU 6

or better yet, create a 'samples.txt' file that describes the data like so:


     #      --samples_file <string>         tab-delimited text file indicating biological replicate relationships.
     #                                   ex.
     #                                        cond_A    cond_A_rep1    A_rep1_left.fq    A_rep1_right.fq
     #                                        cond_A    cond_A_rep2    A_rep2_left.fq    A_rep2_right.fq
     #                                        cond_B    cond_B_rep1    B_rep1_left.fq    B_rep1_right.fq
     #                                        cond_B    cond_B_rep2    B_rep2_left.fq    B_rep2_right.fq

This same samples file can then be used later on with other downstream analysis steps, including expression quantification and differential expression analysis.


>note that fastq files can be gzip-compressed as shown above, in which case they should require a '.gz' extension.



<a name="typical_options"></a>

### Nano setting
```
set regexp
set nowrap
set softwrap
##set const
## Nanorc files
include "/usr/share/nano/nanorc.nanorc"

## C/C++
include "/usr/share/nano/c.nanorc"

## HTML
include "/usr/share/nano/html.nanorc"

## TeX
include "/usr/share/nano/tex.nanorc"

## Quoted emails (under e.g. mutt)
include "/usr/share/nano/mutt.nanorc"

## Patch files
include "/usr/share/nano/patch.nanorc"

## Manpages
include "/usr/share/nano/man.nanorc"

## Groff
include "/usr/share/nano/groff.nanorc"

## Perl
include "/usr/share/nano/perl.nanorc"

## Python
include "/usr/share/nano/python.nanorc"

## Ruby
include "/usr/share/nano/ruby.nanorc"

## Java
include "/usr/share/nano/java.nanorc"

## Assembler
include "/usr/share/nano/asm.nanorc"

## Bourne shell scripts
include "/usr/share/nano/sh.nanorc"

## POV-Ray
include "/usr/share/nano/pov.nanorc"

```


### Trinity run
```

nano trinity.sh

```
###^ means CTRL key
###M- means ALT key



```bash
#!/bin/bash
#SBATCH --job-name=<TRINITY>
#SBATCH --time=15:00
#SBATCH --account=cpu-s2-bch709-0
#SBATCH --partition=cpu-s2-core-0
#SBATCH --cpus-per-task=8
#SBATCH --mem=80g
#SBATCH --mail-type=fail
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-user=<YOUR ID>@unr.edu
#SBATCH --output=<TRINITY>.out


Trinity --seqType fq  --CPU 24 --max_memory 80G --left trim/DT1_R1_val_1.fq.gz,trim/DT2_R1_val_1.fq.gz,trim/DT3_R1_val_1.fq.gz,trim/WT1_R1_val_1.fq.gz,trim/WT2_R1_val_1.fq.gz,trim/WT3_R1_val_1.fq.gz --right trim/DT1_R2_val_2.fq.gz,trim/DT2_R2_val_2.fq.gz,trim/DT3_R2_val_2.fq.gz,trim/WT1_R2_val_2.fq.gz,trim/WT2_R2_val_2.fq.gz,trim/WT3_R2_val_2.fq.gz

```

***Please consider relative path and absolute path***

### Submit job
```
chmod 775 trinity.sh

sbatch trinity.sh
```


>## de Bruijn
>de Bruijn graph construction 
> - Draw (by hand) the de Bruijn graph for the following reads using k=3 (assume all reads are from the forward strand, no sequencing errors)
>AGT   
>ATG  
>CAT  
>GTA  
>GTT  
>TAC    
>TAG  
>TGT  
>TTA  
{: .prereq}


```
[('AG', 'GT'),
 ('AT', 'TG'),
 ('CA', 'AT'),
 ('GT', 'TA'),
 ('GT', 'TT'),
 ('TA', 'AC'),
 ('TA', 'AG'),
 ('TG', 'GT'),
 ('TT', 'TA')]

```
![debj_graph]({{site.baseurl}}/fig/debruijn_graph.png)


### python code
https://colab.research.google.com/drive/12AIJ21eGQ2npxeHcU3h4fT7OxEzxxzin#scrollTo=FN4AvrS0g3uN


### Please check the result
```bash

cd /data/gpfs/assoc/bch709-1/<YOURID>/rnaseq_assembly/trinity_out_dir/

egrep -c ">" Trinity.fasta 

TrinityStats.pl /data/gpfs/assoc/bch709-1/<YOURID>/rnaseq_assembly/trinity_out_dir/Trinity.fasta  > <YOURID>.trinity.stat

cat <YOURID>.trinity.stat
```



### RNA-Seq reads Count Analysis
```bash
pwd
### Your current location is 
## /data/gpfs/assoc/bch709-1/wyim/rnaseq_slurm
align_and_estimate_abundance.pl


nano reads_count.sh
```

### RNA-Seq reads Count Analysis job script
```
#!/bin/bash
#SBATCH --job-name=<JOBNAME>
#SBATCH --cpus-per-task=16
#SBATCH --time=15:00
#SBATCH --mem=16g
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-type=fail
#SBATCH --mail-user=<YOUR ID>@unr.edu
#SBATCH -o <JOBNAME>.out # STDOUT
#SBATCH -e <JOBNAME>.err # STDERR
#SBATCH --account=cpu-s2-bch709-1 
#SBATCH --partition=cpu-s2-core-0


align_and_estimate_abundance.pl --transcripts trinity_out_dir/Trinity.fasta --seqType fq --left trim/DT1_R1_val_1.fq.gz --right trim/DT1_R2_val_2.fq.gz --est_method RSEM --aln_method bowtie2 --trinity_mode --prep_reference --output_dir rsem_outdir_test  --thread_count  16
```
![RSEM]({{site.baseurl}}/fig/RSEM4.png)


### Assignment

WT1,WT2,WT3,DT1,DT2,DT3


```bash
align_and_estimate_abundance.pl --transcripts <TRINITY_FASTA> --seqType fq --left <INPUT_R1> --right <INPUT_R2> --est_method RSEM --aln_method bowtie2 --trinity_mode --prep_reference --output_dir <OUTPUT_FOLDER> --thread_count  16
```
#### Upload six (WT1,WT2,WT3,DT1,DT2,DT3) RSEM.genes.results file
```bash
<OUTPUT_FOLDER>/RSEM.genes.results
```


### Job submission
```bash
sbatch reads_count.sh  
```

### Job check
```
squeue
```
### Job running check
```bash
## do ```ls``` first
ls
cat slurm_<JOBID>.out 

```


### RSEM results check
```bash
less rsem_outdir_test/RSEM.genes.results
```

### Expression values and Normalization

CPM, RPKM, FPKM, TPM, RLE, MRN, Q, UQ, TMM, VST, RLOG, VOOM ... Too many...  

CPM: Controls for sequencing depth when dividing by total count. Not for within-sample comparison or DE.  

Counts per million (CPM) mapped reads are counts scaled by the number of fragments you sequenced (N) times one million. This unit is related to the FPKM without length normalization and a factor of 10^3:  
![CPM]({{site.baseurl}}/fig/CPM.png)

RPKM/FPKM: Controls for sequencing depth and gene length. Good for technical replicates, not good for sample-sample due to compositional bias. Assumes total RNA output is same in all samples. Not for DE.  

TPM: Similar to RPKM/FPKM. Corrects for sequencing depth and gene length. Also comparable between samples but no correction for compositional bias.  

TMM/RLE/MRN: Improved assumption: The output between samples for a core set only of genes is similar. Corrects for compositional bias. Used for DE. RLE and MRN are very similar and correlates well with sequencing depth. edgeR::calcNormFactors() implements TMM, TMMwzp, RLE & UQ.   DESeq2::estimateSizeFactors implements median ratio method (RLE). Does not correct for gene length.  

VST/RLOG/VOOM: Variance is stabilised across the range of mean values. For use in exploratory analyses. Not for DE. vst() and rlog() functions from DESeq2. voom() function from Limma converts data to normal distribution.  

geTMM: Gene length corrected TMM.  

**For DEG using DEG R packages (DESeq2, edgeR, Limma etc), use raw counts  
For visualisation (PCA, clustering, heatmaps etc), use TPM or TMM  
For own analysis with gene length correction, use TPM (maybe geTMM?)  
Other solutions: spike-ins/house-keeping genes**


#### FPKM
![FPKM]({{site.baseurl}}/fig/FPKM.png)

X = mapped reads count
N = number of reads
L = Length of transcripts

```bash
head -n 2 rsem_outdir_test/RSEM.genes.results
```
'length' is this transcript's sequence length (poly(A) tail is not counted). 'effective_length' counts only the positions that can generate a valid fragment.



#### Reads count
```bash
samtools flagstat rsem_outdir_test/bowtie2.bam

cat rsem_outdir_test/RSEM.genes.results | egrep -v FPKM | awk '{ sum+=$5} END {print sum}'
```
### Call Python
```bash
python
```

### FPKM 
Fragments per Kilobase of transcript per million mapped reads


```python
X = 2012.00
Number_Reads_mapped = 559779
Length = 650.92
fpkm= X*(1000/Length)*(1000000/Number_Reads_mapped)
fpkm
```

#### ten to the ninth power = 10\*\*9


```python
fpkm=X/(Number_Reads_mapped*Length)*10**9
fpkm
```


### TPM
 Transcripts Per Million

![TPM]({{site.baseurl}}/fig/TPM.png)

![TPM2]({{site.baseurl}}/fig/TPM2.png)


### Sum of FPKM
```bash
cat rsem_outdir_test/RSEM.genes.results | egrep -v FPKM | awk '{ sum+=$7} END {print sum}'
```
### TPM calculation from FPKM

```python
FPKM = 5521.839239919676
SUM_FPKM = 931616
TPM=(FPKM/SUM_FPKM)*10**6
TPM
```

### TPM calculation from reads count
```bash
cat rsem_outdir_test/RSEM.genes.results | egrep -v FPKM | awk '{ sum+=$5/$4} END {print sum}'
```

```python
sum_count_per_length = 521.494
X = 2012.00
Length = 650.92
TPM = (X/Length)*(1/sum_count_per_length )*10**6
TPM
```

### Paper read
[Li et al., 2010, RSEM](http://bioinformatics.oxfordjournals.org/content/26/4/493.long)  

[Dillies et al., 2013](http://bib.oxfordjournals.org/content/14/6/671.full)


## DEG calculation

### Conda env
```
conda activate rnaseq
```


#### File prepare
```

cd /data/gpfs/assoc/bch709-1/<YOURID>/rnaseq_slurm/
```


### Type sample file
```bash
nano sample.txt
```

###^ means CTRL key
###M- means ALT key


```bash
WT<TAB>WT_REP1<TAB>trim/WT1_R1_val_1.fq.gz<TAB>trim/WT1_R2_val_2.fq.gz
WT<TAB>WT_REP2<TAB>trim/WT2_R1_val_1.fq.gz<TAB>trim/WT2_R2_val_2.fq.gz
WT<TAB>WT_REP3<TAB>trim/WT3_R1_val_1.fq.gz<TAB>trim/WT3_R2_val_2.fq.gz
DT<TAB>DT_REP1<TAB>trim/DT1_R1_val_1.fq.gz<TAB>trim/DT1_R2_val_2.fq.gz
DT<TAB>DT_REP2<TAB>trim/DT2_R1_val_1.fq.gz<TAB>trim/DT2_R2_val_2.fq.gz
DT<TAB>DT_REP3<TAB>trim/DT3_R1_val_1.fq.gz<TAB>trim/DT3_R2_val_2.fq.gz
```

### Change <TAB> charater to Tab
```bash
sed 's/<TAB>/\t/g' sample.txt

sed -i 's/<TAB>/\t/g' sample.txt
```
### Linux command explaination
https://explainshell.com/   

https://explainshell.com/explain?cmd=cat+rsem_outdir_test%2FRSEM.genes.results+%7C+egrep+-v+FPKM+%7C+awk+%27%7B+sum%2B%3D%245%7D+END+%7Bprint+sum%7D%27#  

### SED AWK explaination
https://emb.carnegiescience.edu/sites/default/files/140602-sedawkbash.key_.pdf


### Job file create
```bash
nano alignment.sh
```

### Run alignment 
```bash
#!/bin/bash
#SBATCH --job-name=<JOB_NAME>
#SBATCH --cpus-per-task=32
#SBATCH --time=15:00
#SBATCH --mem=64g
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-type=fail
#SBATCH --mail-user=<YOUR ID>@nevada.unr.edu
#SBATCH -o <JOB_NAME>.out # STDOUT
#SBATCH -e <JOB_NAME>.err # STDERR
#SBATCH --account=cpu-s2-bch709-1 
#SBATCH --partition=cpu-s2-core-0

align_and_estimate_abundance.pl --thread_count 64 --transcripts trinity_out_dir/Trinity.fasta --seqType fq  --est_method RSEM --aln_method bowtie2  --trinity_mode --prep_reference --samples_file sample.txt
```



### Install R-packages
```bash
conda install -c r -c conda-forge -c anaconda -c bioconda  bioconductor-ctc bioconductor-deseq2 bioconductor-edger bioconductor-biobase  bioconductor-qvalue  r-ape  r-gplots  r-fastcluster

```


### abundance_estimates_to_matrix

```
mkdir DEG && cd DEG
nano abundance.sh
```

```bash
#!/bin/bash
#SBATCH --job-name=<JOB_NAME>
#SBATCH --cpus-per-task=1
#SBATCH --time=15:00
#SBATCH --mem=12g
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-type=fail
#SBATCH --mail-user=<YOUR ID>@nevada.unr.edu
#SBATCH -o <JOB_NAME>.out # STDOUT
#SBATCH -e <JOB_NAME>.err # STDERR
#SBATCH --account=cpu-s2-bch709-1 
#SBATCH --partition=cpu-s2-core-0

abundance_estimates_to_matrix.pl  --est_method RSEM --gene_trans_map none --name_sample_by_basedir  --cross_sample_norm TMM ../WT_REP1/RSEM.isoforms.results ../WT_REP2/RSEM.isoforms.results ../WT_REP3/RSEM.isoforms.results   ../DT_REP1/RSEM.isoforms.results ../DT_REP2/RSEM.isoforms.results ../DT_REP3/RSEM.isoforms.results
```
```
sbatch abundance.sh
```
![RSEM]({{site.baseurl}}/fig/RSEM_result.png)



### PtR (Quality Check Your Samples and Biological Replicates)

Once you've performed transcript quantification for each of your biological replicates, it's good to examine the data to ensure that your biological replicates are well correlated, and also to investigate relationships among your samples. If there are any obvious discrepancies among your sample and replicate relationships such as due to accidental mis-labeling of sample replicates, or strong outliers or batch effects, you'll want to identify them before proceeding to subsequent data analyses (such as differential expression).  


```bash

cut -f 1,2 ../sample.txt >> samples_ptr.txt

nano ptr.sh
```

```bash

#!/bin/bash
#SBATCH --job-name=<JOB_NAME>
#SBATCH --cpus-per-task=1
#SBATCH --time=15:00
#SBATCH --mem=12g
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-type=fail
#SBATCH --mail-user=<YOUR ID>@nevada.unr.edu
#SBATCH -o <JOB_NAME>.out # STDOUT
#SBATCH -e <JOB_NAME>.err # STDERR
#SBATCH --account=cpu-s2-bch709-1 
#SBATCH --partition=cpu-s2-core-0

PtR  --matrix RSEM.isoform.counts.matrix --samples samples_ptr.txt --CPM --log2 --min_rowSums 10  --compare_replicates

PtR  --matrix RSEM.isoform.counts.matrix --samples samples_ptr.txt --CPM  --log2 --min_rowSums 10   --sample_cor_matrix

PtR  --matrix RSEM.isoform.counts.matrix --samples samples_ptr.txt --CPM  --log2 --min_rowSums 10   --center_rows --prin_comp 3
```

***Please transfer results to your local computer***

### DEG calculation
```bash
nano deseq.sh
```
```bash

#!/bin/bash
#SBATCH --job-name=<JOB_NAME>
#SBATCH --cpus-per-task=1
#SBATCH --time=15:00
#SBATCH --mem=12g
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-type=fail
#SBATCH --mail-user=<YOUR ID>@nevada.unr.edu
#SBATCH -o <JOB_NAME>.out # STDOUT
#SBATCH -e <JOB_NAME>.err # STDERR
#SBATCH --account=cpu-s2-bch709-1 
#SBATCH --partition=cpu-s2-core-0


run_DE_analysis.pl --matrix RSEM.isoform.counts.matrix --samples_file samples_ptr.txt --method DESeq2 
```
```bash
nano edgeR.sh
```
```bash

#!/bin/bash
#SBATCH --job-name=<JOB_NAME>
#SBATCH --cpus-per-task=1
#SBATCH --time=15:00
#SBATCH --mem=12g
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-type=fail
#SBATCH --mail-user=<YOUR ID>@nevada.unr.edu
#SBATCH -o <JOB_NAME>.out # STDOUT
#SBATCH -e <JOB_NAME>.err # STDERR
#SBATCH --account=cpu-s2-bch709-1 
#SBATCH --partition=cpu-s2-core-0

run_DE_analysis.pl --matrix RSEM.isoform.counts.matrix --samples_file samples_ptr.txt --method edgeR
```

```bash
cd DESeq2.XXXXX.dir
analyze_diff_expr.pl --matrix ../RSEM.isoform.TMM.EXPR.matrix  -P 0.001 -C 1  --samples ../samples_ptr.txt
wc -l RSEM.isoform.counts.matrix.DT_vs_WT.DESeq2.DE_results.P0.001_C1.DE.subset
cd ../

cd edgeR.XXXXX.dir
analyze_diff_expr.pl --matrix ../RSEM.isoform.TMM.EXPR.matrix  -P 0.001 -C 1  --samples ../samples_ptr.txt
wc -l RSEM.isoform.counts.matrix.DT_vs_WT.edgeR.DE_results.P0.001_C1.DE.subset
```

## DESeq2 vs EdgeR Normalization method
DESeq and EdgeR are very similar and both assume that no genes are differentially expressed. DEseq uses a "geometric" normalisation strategy, whereas EdgeR is a weighted mean of log ratios-based method. Both normalise data initially via the calculation of size / normalisation factors.


### DESeq
DESeq: This normalization method is included in the DESeq Bioconductor package (version 1.6.0) and is based on the hypothesis that most genes are not DE. A DESeq scaling factor for a given lane is computed as the median of the ratio, for each gene, of its read count over its geometric mean across all lanes. The underlying idea is that non-DE genes should have similar read counts across samples, leading to a ratio of 1. **Assuming most genes are not DE, the median of this ratio for the lane provides an estimate of the correction factor that should be applied to all read counts of this lane to fulfill the hypothesis.** By calling the estimateSizeFactors() and sizeFactors() functions in the DESeq Bioconductor package, this factor is computed for each lane, and raw read counts are divided by the factor associated with their sequencing lane.  
[DESeq2](https://www.ncbi.nlm.nih.gov/pubmed/22988256)

### EdgeR
Trimmed Mean of M-values (TMM): This normalization method is implemented in the edgeR Bioconductor package (version 2.4.0). It is also based on the hypothesis that most genes are not DE. The TMM factor is computed for each lane, with one lane being considered as a reference sample and the others as test samples. For each test sample, TMM is computed as the weighted mean of log ratios between this test and the reference, after exclusion of the most expressed genes and the genes with the largest log ratios. **According to the hypothesis of low DE, this TMM should be close to 1. If it is not, its value provides an estimate of the correction factor that must be applied to the library sizes (and not the raw counts) in order to fulfill the hypothesis.** The calcNormFactors() function in the edgeR Bioconductor package provides these scaling factors. To obtain normalized read counts, these normalization factors are re-scaled by the mean of the normalized library sizes. Normalized read counts are obtained by dividing raw read counts by these re-scaled normalization factors.  
[EdgeR](https://www.ncbi.nlm.nih.gov/pubmed/22988256)

## DESeq2 vs EdgeR Statistical tests for differential expression
### DESeq2
DESeq2 uses raw counts, rather than normalized count data, and models the normalization to fit the counts within a Generalized Linear Model (GLM) of the negative binomial family with a logarithmic link. Statistical tests are then performed to assess differential expression, if any.  

### EdgeR
Data are normalized to account for sample size differences and variance among samples. The normalized count data are used to estimate per-gene fold changes and to perform statistical tests of whether each gene is likely to be differentially expressed.  
EdgeR uses an exact test under a negative binomial distribution (Robinson and Smyth, 2008). The statistical test is related to Fisher's exact test, though Fisher uses a different distribution.  


## Negative binormal
### DESeq2 
ϕ was assumed to be a function of μ determined by nonparametric regression. The recent version used in this paper follows a more versatile procedure. Firstly, for each transcript, an estimate of the dispersion is made, presumably using maximum likelihood. Secondly, the estimated dispersions for all transcripts are fitted to the functional form:  
ϕ=a+bμ(DESeq parametric fit), using a gamma-family generalised linear model  (Using regression)

### EdgeR
edgeR recommends a “tagwise dispersion” function, which estimates the dispersion on a gene-by-gene basis, and implements an empirical Bayes strategy for squeezing the estimated dispersions towards the common dispersion. Under the default setting, the degree of squeezing is adjusted to suit the number of biological replicates within each condition: more biological replicates will need to borrow less information from the complete set of transcripts and require less squeezing.  

***[DEG software comparison paper](https://bmcgenomics.biomedcentral.com/articles/10.1186/1471-2164-13-484)***

*** DEseq uses a "geometric" normalisation strategy, whereas EdgeR is a weighted mean of log ratios-based method. Both normalise data initially via the calculation of size / normalisation factors ***




### Draw Venn Diagram
```bash
conda create -n venn python=2.7  
conda activate venn  
conda install -c bioconda  bedtools intervene   
conda install -c r r-UpSetR r-corrplot r-Cairo  

``` 

```bash
cd ../
pwd
# /data/gpfs/assoc/bch709-1/wyim/rnaseq_slurm/DEG
mkdir Venn 
cd Venn

###DESeq2
cut -f 1 ../DEG/DESeq2.#####.dir/RSEM.isoform.counts.matrix.DT_vs_WT.DESeq2.DE_results.P0.001_C1.DT-UP.subset  | grep -v sample > DESeq.UP.subset
cut -f 1 ../DEG/DESeq2.####.dir/RSEM.isoform.counts.matrix.DT_vs_WT.DESeq2.DE_results.P0.001_C1.WT-UP.subset  | grep -v sample > DESeq.DOWN.subset

###edgeR
cut -f 1 ../DEG/edgeR.####.dir/RSEM.isoform.counts.matrix.DT_vs_WT.edgeR.DE_results.P0.001_C1.DT-UP.subset   | grep -v sample > edgeR.UP.subset
cut -f 1 ../DEG/edgeR.####.dir/RSEM.isoform.counts.matrix.DT_vs_WT.edgeR.DE_results.P0.001_C1.WT-UP.subset   | grep -v sample > edgeR.DOWN.subset


### Drawing
intervene venn -i DESeq.DOWN.subset DESeq.UP.subset edgeR.DOWN.subset edgeR.UP.subset  --type list --save-overlaps
intervene upset -i DESeq.DOWN.subset DESeq.UP.subset edgeR.DOWN.subset edgeR.UP.subset  --type list --save-overlaps
intervene pairwise  -i DESeq.DOWN.subset DESeq.UP.subset edgeR.DOWN.subset edgeR.UP.subset  --type list
```
