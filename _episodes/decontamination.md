---
layout: page
title: metagenome
published: true
---

## making the working directory
```bash
echo $USER

cd /data/gpfs/assoc/bch709-3/${USER}

mkdir -p /data/gpfs/assoc/bch709-3/${USER}/metagenome/fastq
mkdir  /data/gpfs/assoc/bch709-3/${USER}/metagenome/database 

cd /data/gpfs/assoc/bch709-3/${USER}/metagenome

```
# Environment creation
```bash
conda create -n BCH709_metagenome -c bioconda -c conda-forge -c anaconda mamba -y 

conda activate BCH709_metagenome

mamba install -c bioconda -c conda-forge -c anaconda multiqc=1.9 kraken2 -y

git clone https://github.com/jenniferlu717/bracken
cd bracken/
bash install_bracken.sh
cd /data/gpfs/assoc/bch709-3/${USER}/metagenome
```



## Activate environment
```bash
conda activate BCH709_metagenome
```


## Link fastq file
```bash
ln -s /data/gpfs/assoc/bch709-3/Course_materials/metagenome/*.fastq /data/gpfs/assoc/bch709-3/${USER}/metagenome/fastq
```

## Link database file
```bash
ln -s /data/gpfs/assoc/bch709-3/Course_materials/metagenome/database/* /data/gpfs/assoc/bch709-3/${USER}/metagenome/database
```

# Investigate taxa

Here we introduce a software called Kraken2. This tool uses k-mers to assign a taxonomic labels in form of NCBI Taxonomy to the sequence (if possible). The taxonomic label is assigned based on similar k-mer content of the sequence in question to the k-mer content of reference genome sequence. The result is a classification of the sequence in question to the most likely taxonomic label. If the k-mer content is not similar to any genomic sequence in the database used, it will not assign any taxonomic label.

## Kraken2 
Kraken 2, KrakenUniq and Bracken indexes
Kraken 2 is a fast and memory efficient tool for taxonomic assignment of metagenomics sequencing reads. Bracken is a related tool that additionally estimates relative abundances of species or genera. See the Kraken 2 manual for more information about the individual libraries and their relationship to public repositories like Refseq. See also the Kraken protocol for advice on how to use it.
https://www.nature.com/articles/s41596-022-00738-y

https://benlangmead.github.io/aws-indexes/k2


![](https://i.imgur.com/6Ck2722.png)

## The most recent database
Kraken 2 database (Standard 46G) includes archaea, bacteria, viral, plasmid, human1, UniVec_Core has been downloaded in below. We already decompress and linked by `ln` command above.
https://benlangmead.github.io/aws-indexes/k2

```bash
/data/gpfs/assoc/bch709-3/Course_materials/metagenome/k2_standard_20220926.tar.gz
```


### Create file list
```bash
cd /data/gpfs/assoc/bch709-3/${USER}/metagenome/fastq


ls -1 *.fastq | sed 's/\.fastq//g' | sort -u > /data/gpfs/assoc/bch709-3/${USER}/metagenome/filelist

cat /data/gpfs/assoc/bch709-3/${USER}/metagenome/filelist
```

#### Copy templet
```bash
cp /data/gpfs/assoc/bch709-3/Course_materials/metagenome/run.sh /data/gpfs/assoc/bch709-3/${USER}/metagenome/metagenome.sh

sed -i "s/16g/128g/g; s/\-\-cpus\-per\-task\=2/\-\-cpus\-per\-task\=2/g; s/\[NAME\]/Metagenome/g; s/\[youremail\]/${USER}\@unr.edu\,${USER}\@nevada.unr.edu/g" /data/gpfs/assoc/bch709-3/${USER}/metagenome/metagenome.sh

```


i=SRR19419488

sbatch -c 2 --mem=128g --wrap="kraken2  --threads 2 --report ${i}.kreport --db /data/gpfs/assoc/bch709-3/${USER}/metagenome/database ${i}.fastq"

KRAKEN_DB=/data/gpfs/assoc/bch709-3/${USER}/metagenome/database 
READ_LEN=75
CLASSIFICATION_LVL=S
THRESHOLD=2 
BRACKEN_OUTPUT_FILE=wastewater
SAMPLE=SRR19419488


python /data/gpfs/assoc/bch709-3/${USER}/metagenome/bracken/src/est_abundance.py  -i ${SAMPLE}.kreport -k ${KRAKEN_DB}/database${READ_LEN}mers.kmer_distrib -l ${CLASSIFICATION_LVL} -t ${THRESHOLD} -o ${BRACKEN_OUTPUT_FILE}_${SAMPLE}.bracken

python /data/gpfs/assoc/bch709-3/${USER}/metagenome/bracken/analysis_scripts/combine_bracken_outputs.py


### Donwload most recent database
```bash
https://benlangmead.github.io/aws-indexes/k2
wget https://genome-idx.s3.amazonaws.com/kraken/k2_standard_16gb_20200919.tar.gz
```

```bash
cd /data/gpfs/assoc/bch709-1/<YOURID>/RNA-Seq_example

mkdir Kraken2 && cd Kraken2
kraken2-inspect --db EXAMPLE_DB | head -5

kraken2  --threads 24 --report SRR11968954 --db /data/gpfs/assoc/bch709-1/Course_material/database/ /data/gpfs/assoc/bch709-1/<YOURID>/RNA-Seq_example/Drosophila/raw_data/SRR11968954.fastq.gz
```
```
"C"/"U": a one letter code indicating that the sequence was either classified or unclassified.
```
```bash
https://fbreitwieser.shinyapps.io/pavian/
```






