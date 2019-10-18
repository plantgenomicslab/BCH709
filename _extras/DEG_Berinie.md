---
layout: page
title: Mastering-markdown
published: true
---
## Before you start
Please prepare following files and software. We cannot help the installation and download.
When you excute example code under code part, please change < > to your system.

> ## Prerequisite software
> We recommend to use Conda. For installation please refer [this link](https://plantgenomicslab.github.io/BCH709/conda/index.html)
> ### List of software
>\_libgcc_mutex=0.1=main
>\_r-mutex=1.0.1=anacondar_1
>alsa-lib=1.1.5=h516909a_1001
>asn1crypto=0.24.0=py37_1003
>binutils_impl_linux-64=2.31.1=h6176602_1
>binutils_linux-64=2.31.1=h6176602_9
>bioconductor-biobase=2.42.0=r351h14c3975_1
>bioconductor-biocgenerics=0.28.0=r351_1
>bioconductor-ebseq=1.24.0=r351_0
>bowtie=1.2.3=py37hc9558a2_0
>bowtie2=2.3.5=py37he860b03_0
>bwidget=1.9.11=0
>bzip2=1.0.8=h516909a_1
>ca-certificates=2019.9.11=hecc5488_0
>cairo=1.16.0=hfb77d84_1002
>certifi=2019.9.11=py37_0
>cffi=1.12.3=py37h8022711_0
>chardet=3.0.4=py37_1003
>click=7.0=py_0
>colormath=3.0.0=py_2
>cryptography=2.7=py37h72c5cf5_0
>curl=7.65.3=hf8cf82a_0
>cutadapt=2.5=py37h516909a_0
>cycler=0.10.0=py_1
>dbus=1.13.6=he372182_0
>decorator=4.4.0=py_0
>dnaio=0.3=py37h14c3975_1
>expat=2.2.5=he1b5a44_1003
>fastqc=0.11.8=1
>font-ttf-dejavu-sans-mono=2.37=hab24e00_0
>fontconfig=2.13.1=h86ecdb6_1001
>freetype=2.10.0=he983fc9_1
>fribidi=1.0.5=h516909a_1002
>future=0.17.1=py37_1000
>gcc_impl_linux-64=7.3.0=habb00fd_1
>gcc_linux-64=7.3.0=h553295d_9
>gettext=0.19.8.1=hc5be6a0_1002
>gfortran_impl_linux-64=7.3.0=hdf63c60_1
>gfortran_linux-64=7.3.0=h553295d_9
>gffread=0.11.4
>giflib=5.1.7=h516909a_1
>glib=2.58.3=h6f030ca_1002
>graphite2=1.3.13=hf484d3e_1000
>gsl=2.5=h294904e_1
>gst-plugins-base=1.14.5=h0935bb2_0
>gstreamer=1.14.5=h36ae1b5_0
>gxx_impl_linux-64=7.3.0=hdf63c60_1
>gxx_linux-64=7.3.0=h553295d_9
>harfbuzz=2.4.0=h9f30f68_3
>htslib=1.9=h4da6232_3
>icu=64.2=he1b5a44_1
>idna=2.8=py37_1000
>jellyfish=2.2.10=h6bb024c_1
>jemalloc=5.2.1=he1b5a44_0
>jinja2=2.10.1=py_0
>jpeg=9c=h14c3975_1001
>kiwisolver=1.1.0=py37hc9558a2_0
>krb5=1.16.3=h05b26f9_1001
>lcms2=2.9=h2e4bb80_0
>libblas=3.8.0=12_openblas
>libcblas=3.8.0=12_openblas
>libcurl=7.65.3=hda55be3_0
>libdeflate=1.3=h516909a_0
>libedit=3.1.20170329=hf8c457e_1001
>libffi=3.2.1=he1b5a44_1006
>libgcc=7.2.0=h69d50b8_2
>libgcc-ng=9.1.0=hdf63c60_0
>libgfortran-ng=7.3.0=hdf63c60_0
>libiconv=1.15=h516909a_1005
>liblapack=3.8.0=12_openblas
>libopenblas=0.3.7=h6e990d7_1
>libpng=1.6.37=hed695b0_0
>libssh2=1.8.2=h22169c7_2
>libstdcxx-ng=9.1.0=hdf63c60_0
>libtiff=4.0.10=h57b8799_1003
>libuuid=2.32.1=h14c3975_1000
>libxcb=1.13=h14c3975_1002
>libxml2=2.9.9=hee79883_5
>lz4-c=1.8.3=he1b5a44_1001
>lzstring=1.0.4=py_1001
>make=4.2.1=h14c3975_2004
>markdown=3.1.1=py_0
>markupsafe=1.1.1=py37h14c3975_0
>matplotlib=3.1.1=py37_1
>matplotlib-base=3.1.1=py37he7580a8_1
>multiqc=1.7=py_4
>mysql-connector-c=6.1.11=hd2bbab6_1003
>ncurses=6.1=hf484d3e_1002
>networkx=2.3=py_0
>numpy=1.17.2=py37h95a1406_0
>openjdk=11.0.1=h46a85a0_1017
>openssl=1.1.1c=h516909a_0
>pango=1.42.4=ha030887_1
>pcre=8.41=hf484d3e_1003
>perl=5.26.2=h516909a_1006
>perl-app-cpanminus=1.7044=pl526_1
>perl-carp=1.38=pl526_3
>perl-constant=1.33=pl526_1
>perl-cpan-meta=2.150010=pl526_0
>perl-cpan-meta-requirements=2.140=pl526_0
>perl-cpan-meta-yaml=0.018=pl526_0
>perl-data-dumper=2.173=pl526_0
>perl-encode=2.88=pl526_1
>perl-exporter=5.72=pl526_1
>perl-extutils-cbuilder=0.280230=pl526_1
>perl-extutils-makemaker=7.36=pl526_1
>perl-extutils-manifest=1.72=pl526_0
>perl-extutils-parsexs=3.35=pl526_0
>perl-file-path=2.16=pl526_0
>perl-file-temp=0.2304=pl526_2
>perl-getopt-long=2.50=pl526_1
>perl-ipc-cmd=1.02=pl526_0
>perl-json-pp=4.04=pl526_0
>perl-locale-maketext-simple=0.21=pl526_2
>perl-module-build=0.4224=pl526_3
>perl-module-corelist=5.20190524=pl526_0
>perl-module-load=0.32=pl526_1
>perl-module-load-conditional=0.68=pl526_2
>perl-module-metadata=1.000036=pl526_0
>perl-params-check=0.38=pl526_1
>perl-parent=0.236=pl526_1
>perl-perl-ostype=1.010=pl526_1
>perl-scalar-list-utils=1.52=pl526h516909a_0
>perl-text-abbrev=1.02=pl526_0
>perl-text-parsewords=3.30=pl526_0
>perl-version=0.9924=pl526_0
>pigz=2.3.4=0
>pip=19.2.3=py37_0
>pixman=0.38.0=h516909a_1003
>pthread-stubs=0.4=h14c3975_1001
>pycparser=2.19=py37_1
>pyopenssl=19.0.0=py37_0
>pyparsing=2.4.2=py_0
>pyqt=5.9.2=py37hcca6a23_4
>pysocks=1.7.1=py37_0
>python=3.7.3=h33d41f4_1
>python-dateutil=2.8.0=py_0
>pyyaml=5.1.2=py37h516909a_0
>qt=5.9.7=h0c104cb_3
>r-assertthat=0.2.1=r35h6115d3f_1
>r-base=3.5.1=hbc19587_1011
>r-bibtex=0.4.2=r35hcdcec82_1003
>r-bitops=1.0_6=r35hcdcec82_1003
>r-blockmodeling=0.3.4=r35h9bbef5b_2
>r-catools=1.17.1.2=r35h0357c0b_1
>r-cli=1.1.0=r35h6115d3f_2
>r-codetools=0.2_16=r35h6115d3f_1001
>r-crayon=1.3.4=r35h6115d3f_1002
>r-digest=0.6.21=r35h0357c0b_0
>r-doparallel=1.0.15=r35h6115d3f_0
>r-dorng=1.7.1=r35h6115d3f_1001
>r-evaluate=0.14=r35h6115d3f_1
>r-foreach=1.4.7=r35h6115d3f_0
>r-gdata=2.18.0=r35h6115d3f_1002
>r-glue=1.3.1=r35hcdcec82_1
>r-gplots=3.0.1.1=r35h6115d3f_1
>r-gtools=3.8.1=r35hcdcec82_1003
>r-iterators=1.0.12=r35h6115d3f_0
>r-kernsmooth=2.23_15=r35h9bbef5b_1004
>r-lattice=0.20_38=r35hcdcec82_1002
>r-magrittr=1.5=r35h6115d3f_1002
>r-matrix=1.2_17=r35hcdcec82_1
>r-pkgmaker=0.27=r35h6115d3f_1001
>r-praise=1.0.0=r35h6115d3f_1003
>r-r6=2.4.0=r35h6115d3f_2
>r-registry=0.5_1=r35h6115d3f_1
>r-rlang=0.4.0=r35hcdcec82_1
>r-rngtools=1.4=r35h6115d3f_1
>r-stringi=1.4.3=r35h0e574ca_3
>r-stringr=1.4.0=r35h6115d3f_1
>r-testthat=2.2.1=r35h0357c0b_0
>r-withr=2.1.2=r35h6115d3f_1001
>r-xtable=1.8_4=r35h6115d3f_2
>readline=8.0=hf8c457e_0
>requests=2.22.0=py37_1
>rsem=1.3.2=pl526r351hc0aa232_0
>salmon=0.14.1=ha0cc327_2
>samtools=1.9=h10a08f8_12
>setuptools=41.2.0=py37_0
>simplejson=3.16.0=py37h14c3975_1002
>sip=4.19.8=py37hf484d3e_1000
>six=1.12.0=py37_1000
>spectra=0.0.11=py_1
>sqlite=3.29.0=hcee41ef_1
>star=2.7.2b=0
>subread=1.6.4=h84994c4_1
>tbb=2019.8=hc9558a2_0
>tk=8.6.9=hed695b0_1003
>tktable=2.10=h555a92e_2
>tornado=6.0.3=py37h516909a_0
>trim-galore=0.6.4=0
>trimmomatic=0.39=1
>trinity=2.8.5=h8b12597_2
>ucsc-bigwigsummary=357=1
>urllib3=1.25.5=py37_0
>wheel=0.33.6=py37_0
>xopen=0.8.2=py37_0
>xorg-fixesproto=5.0=h14c3975_1002
>xorg-inputproto=2.3.2=h14c3975_1002
>xorg-kbproto=1.0.7=h14c3975_1002
>xorg-libice=1.0.10=h516909a_0
>xorg-libsm=1.2.3=h84519dc_1000
>xorg-libx11=1.6.8=h516909a_0
>xorg-libxau=1.0.9=h14c3975_0
>xorg-libxdmcp=1.1.3=h516909a_0
>xorg-libxext=1.3.4=h516909a_0
>xorg-libxfixes=5.0.3=h516909a_1004
>xorg-libxi=1.7.10=h516909a_0
>xorg-libxrender=0.9.10=h516909a_1002
>xorg-libxtst=1.2.3=h14c3975_1002
>xorg-recordproto=1.14.2=h14c3975_1002
>xorg-renderproto=0.11.1=h14c3975_1002
>xorg-xextproto=7.3.0=h14c3975_1002
>xorg-xproto=7.0.31=h14c3975_1007
>xz=5.2.4=h14c3975_1001
>yaml=0.1.7=h14c3975_1001
>zlib=1.2.11=h516909a_1006
>zstd=1.4.0=h3b9ef0a_0
>
{: .prereq}

## Prerequisite files
Please download Arabidopsisa genome and GFF3 from [Phytozome](https://phytozome.jgi.doe.gov/pz/portal.html)
- Athaliana_167_TAIR10.gene.gff3
- Athaliana_447_TAIR10.fa


### GFF to GTF conversion
```bash
gffread Athaliana_167_TAIR10.gene.gff3 -T -o Athaliana_167_TAIR10.gene.gtf
```

### Reference indexing by STAR
We recommned to put your `Athaliana_167_TAIR10.gene.gff3` `Athaliana_167_TAIR10.gene.gtf` and `Athaliana_447_TAIR10.fa` to reference folder.
https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3530905/

```bash
STAR  --runThreadN <YOUR THREAD> --runMode genomeGenerate --genomeDir reference/ --genomeFastaFiles  reference/Athaliana_447_TAIR10.fa --sjdbGTFfile  reference/Athaliana_167_TAIR10.gene.gtf  --sjdbOverhang 99 
```

### Sequencing Reads Trimming
We recommend to put your raw reads under ```raw``` folder and trim reads under `trim` folder.
We used ```trim_galore``` to remove adaptor sequence and low quality base.
https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/

Below `EXAMPLE` command shows execution for `KR24D1_1.fq.gz` and `KR24D1_2.fq.gz`
```bash
trim_galore --paired --three_prime_clip_R1 20 --three_prime_clip_R2 20 --cores 8 --max_n 40 --gzip -o trim raw/KR24D1_1.fq.gz raw/KR24D1_2.fq.gz
```

Below`EXAMPLE` command shows the example for [Slurm HPC enviroment](https://en.wikipedia.org/wiki/Slurm_Workload_Manager)
```bash
sbatch -N 1  -c 8 --mem=16g -A <YOURACCOUNT> -p <YOURPARTITION> --wrap="trim_galore --paired --three_prime_clip_R1 20 --three_prime_clip_R2 20 --cores 8 --max_n 40 --gzip -o trim raw/KR24D1_1.fq.gz raw/KR24D1_2.fq.gz"
sbatch -N 1  -c 8 --mem=16g -A <YOURACCOUNT> -p <YOURPARTITION> --wrap="trim_galore --paired --three_prime_clip_R1 20 --three_prime_clip_R2 20 --cores 8 --max_n 40 --gzip -o trim raw/KR24D2_1.fq.gz raw/KR24D2_2.fq.gz"
sbatch -N 1  -c 8 --mem=16g -A <YOURACCOUNT> -p <YOURPARTITION> --wrap="trim_galore --paired --three_prime_clip_R1 20 --three_prime_clip_R2 20 --cores 8 --max_n 40 --gzip -o trim raw/KR24D3_1.fq.gz raw/KR24D3_2.fq.gz"
sbatch -N 1  -c 8 --mem=16g -A <YOURACCOUNT> -p <YOURPARTITION> --wrap="trim_galore --paired --three_prime_clip_R1 20 --three_prime_clip_R2 20 --cores 8 --max_n 40 --gzip -o trim raw/KR24D4_1.fq.gz raw/KR24D4_2.fq.gz"
sbatch -N 1  -c 8 --mem=16g -A <YOURACCOUNT> -p <YOURPARTITION> --wrap="trim_galore --paired --three_prime_clip_R1 20 --three_prime_clip_R2 20 --cores 8 --max_n 40 --gzip -o trim raw/KR24W1_1.fq.gz raw/KR24W1_2.fq.gz"
sbatch -N 1  -c 8 --mem=16g -A <YOURACCOUNT> -p <YOURPARTITION> --wrap="trim_galore --paired --three_prime_clip_R1 20 --three_prime_clip_R2 20 --cores 8 --max_n 40 --gzip -o trim raw/KR24W2_1.fq.gz raw/KR24W2_2.fq.gz"
sbatch -N 1  -c 8 --mem=16g -A <YOURACCOUNT> -p <YOURPARTITION> --wrap="trim_galore --paired --three_prime_clip_R1 20 --three_prime_clip_R2 20 --cores 8 --max_n 40 --gzip -o trim raw/KR24W3_1.fq.gz raw/KR24W3_2.fq.gz"
sbatch -N 1  -c 8 --mem=16g -A <YOURACCOUNT> -p <YOURPARTITION> --wrap="trim_galore --paired --three_prime_clip_R1 20 --three_prime_clip_R2 20 --cores 8 --max_n 40 --gzip -o trim raw/KR24W4_1.fq.gz raw/KR24W4_2.fq.gz"
sbatch -N 1  -c 8 --mem=16g -A <YOURACCOUNT> -p <YOURPARTITION> --wrap="trim_galore --paired --three_prime_clip_R1 20 --three_prime_clip_R2 20 --cores 8 --max_n 40 --gzip -o trim raw/KRWTD1_1.fq.gz raw/KRWTD1_2.fq.gz"
sbatch -N 1  -c 8 --mem=16g -A <YOURACCOUNT> -p <YOURPARTITION> --wrap="trim_galore --paired --three_prime_clip_R1 20 --three_prime_clip_R2 20 --cores 8 --max_n 40 --gzip -o trim raw/KRWTD2_1.fq.gz raw/KRWTD2_2.fq.gz"
sbatch -N 1  -c 8 --mem=16g -A <YOURACCOUNT> -p <YOURPARTITION> --wrap="trim_galore --paired --three_prime_clip_R1 20 --three_prime_clip_R2 20 --cores 8 --max_n 40 --gzip -o trim raw/KRWTD3_1.fq.gz raw/KRWTD3_2.fq.gz"
sbatch -N 1  -c 8 --mem=16g -A <YOURACCOUNT> -p <YOURPARTITION> --wrap="trim_galore --paired --three_prime_clip_R1 20 --three_prime_clip_R2 20 --cores 8 --max_n 40 --gzip -o trim raw/KRWTD4_1.fq.gz raw/KRWTD4_2.fq.gz"
sbatch -N 1  -c 8 --mem=16g -A <YOURACCOUNT> -p <YOURPARTITION> --wrap="trim_galore --paired --three_prime_clip_R1 20 --three_prime_clip_R2 20 --cores 8 --max_n 40 --gzip -o trim raw/KRWTW1_1.fq.gz raw/KRWTW1_2.fq.gz"
sbatch -N 1  -c 8 --mem=16g -A <YOURACCOUNT> -p <YOURPARTITION> --wrap="trim_galore --paired --three_prime_clip_R1 20 --three_prime_clip_R2 20 --cores 8 --max_n 40 --gzip -o trim raw/KRWTW2_1.fq.gz raw/KRWTW2_2.fq.gz"
sbatch -N 1  -c 8 --mem=16g -A <YOURACCOUNT> -p <YOURPARTITION> --wrap="trim_galore --paired --three_prime_clip_R1 20 --three_prime_clip_R2 20 --cores 8 --max_n 40 --gzip -o trim raw/KRWTW3_1.fq.gz raw/KRWTW3_2.fq.gz"
sbatch -N 1  -c 8 --mem=16g -A <YOURACCOUNT> -p <YOURPARTITION> --wrap="trim_galore --paired --three_prime_clip_R1 20 --three_prime_clip_R2 20 --cores 8 --max_n 40 --gzip -o trim raw/KRWTW4_1.fq.gz raw/KRWTW4_2.fq.gz"
```




### Sequencing Reads Alignment
Below`EXAMPLE` command shows the example for [Slurm HPC enviroment](https://en.wikipedia.org/wiki/Slurm_Workload_Manager)
Following command will map to genome index `reference` folder with trimmed reads under `trim` folder.
For options, please check the manual.
http://labshare.cshl.edu/shares/gingeraslab/www-data/dobin/STAR/STAR.posix/doc/STARmanual.pdf


```bash
sbatch -A <YOURACCOUNT> -p <YOURPARTITION> -N 1 -c 8 --mem=16g --wrap="STAR --runMode alignReads --runThreadN 8 --readFilesCommand zcat --outFilterMultimapNmax 10 --alignIntronMin 61 --alignIntronMax 10000 --genomeDir reference/ --readFilesIn trim/KR24D1_1_val_1.fq.gz trim/KR24D1_2_val_2.fq.gz --outSAMtype BAM SortedByCoordinate --outFileNamePrefix bam/KR24D1.bam"

sbatch -A <YOURACCOUNT> -p <YOURPARTITION> -N 1 -c 8 --mem=16g --wrap="STAR --runMode alignReads --runThreadN 8 --readFilesCommand zcat --outFilterMultimapNmax 10 --alignIntronMin 61 --alignIntronMax 10000 --genomeDir reference/ --readFilesIn trim/KR24D2_1_val_1.fq.gz trim/KR24D2_2_val_2.fq.gz --outSAMtype BAM SortedByCoordinate --outFileNamePrefix bam/KR24D2.bam"

sbatch -A <YOURACCOUNT> -p <YOURPARTITION> -N 1 -c 8 --mem=16g --wrap="STAR --runMode alignReads --runThreadN 8 --readFilesCommand zcat --outFilterMultimapNmax 10 --alignIntronMin 61 --alignIntronMax 10000 --genomeDir reference/ --readFilesIn trim/KR24D3_1_val_1.fq.gz trim/KR24D3_2_val_2.fq.gz --outSAMtype BAM SortedByCoordinate --outFileNamePrefix bam/KR24D3.bam"

sbatch -A <YOURACCOUNT> -p <YOURPARTITION> -N 1 -c 8 --mem=16g --wrap="STAR --runMode alignReads --runThreadN 8 --readFilesCommand zcat --outFilterMultimapNmax 10 --alignIntronMin 61 --alignIntronMax 10000 --genomeDir reference/ --readFilesIn trim/KR24D4_1_val_1.fq.gz trim/KR24D4_2_val_2.fq.gz --outSAMtype BAM SortedByCoordinate --outFileNamePrefix bam/KR24D4.bam"

sbatch -A <YOURACCOUNT> -p <YOURPARTITION> -N 1 -c 8 --mem=16g --wrap="STAR --runMode alignReads --runThreadN 8 --readFilesCommand zcat --outFilterMultimapNmax 10 --alignIntronMin 61 --alignIntronMax 10000 --genomeDir reference/ --readFilesIn trim/KR24W1_1_val_1.fq.gz trim/KR24W1_2_val_2.fq.gz --outSAMtype BAM SortedByCoordinate --outFileNamePrefix bam/KR24W1.bam"

sbatch -A <YOURACCOUNT> -p <YOURPARTITION> -N 1 -c 8 --mem=16g --wrap="STAR --runMode alignReads --runThreadN 8 --readFilesCommand zcat --outFilterMultimapNmax 10 --alignIntronMin 61 --alignIntronMax 10000 --genomeDir reference/ --readFilesIn trim/KR24W2_1_val_1.fq.gz trim/KR24W2_2_val_2.fq.gz --outSAMtype BAM SortedByCoordinate --outFileNamePrefix bam/KR24W2.bam"

sbatch -A <YOURACCOUNT> -p <YOURPARTITION> -N 1 -c 8 --mem=16g --wrap="STAR --runMode alignReads --runThreadN 8 --readFilesCommand zcat --outFilterMultimapNmax 10 --alignIntronMin 61 --alignIntronMax 10000 --genomeDir reference/ --readFilesIn trim/KR24W3_1_val_1.fq.gz trim/KR24W3_2_val_2.fq.gz --outSAMtype BAM SortedByCoordinate --outFileNamePrefix bam/KR24W3.bam"

sbatch -A <YOURACCOUNT> -p <YOURPARTITION> -N 1 -c 8 --mem=16g --wrap="STAR --runMode alignReads --runThreadN 8 --readFilesCommand zcat --outFilterMultimapNmax 10 --alignIntronMin 61 --alignIntronMax 10000 --genomeDir reference/ --readFilesIn trim/KR24W4_1_val_1.fq.gz trim/KR24W4_2_val_2.fq.gz --outSAMtype BAM SortedByCoordinate --outFileNamePrefix bam/KR24W4.bam"

sbatch -A <YOURACCOUNT> -p <YOURPARTITION> -N 1 -c 8 --mem=16g --wrap="STAR --runMode alignReads --runThreadN 8 --readFilesCommand zcat --outFilterMultimapNmax 10 --alignIntronMin 61 --alignIntronMax 10000 --genomeDir reference/ --readFilesIn trim/KRWTD1_1_val_1.fq.gz trim/KRWTD1_2_val_2.fq.gz --outSAMtype BAM SortedByCoordinate --outFileNamePrefix bam/KRWTD1.bam"

sbatch -A <YOURACCOUNT> -p <YOURPARTITION> -N 1 -c 8 --mem=16g --wrap="STAR --runMode alignReads --runThreadN 8 --readFilesCommand zcat --outFilterMultimapNmax 10 --alignIntronMin 61 --alignIntronMax 10000 --genomeDir reference/ --readFilesIn trim/KRWTD2_1_val_1.fq.gz trim/KRWTD2_2_val_2.fq.gz --outSAMtype BAM SortedByCoordinate --outFileNamePrefix bam/KRWTD2.bam"

sbatch -A <YOURACCOUNT> -p <YOURPARTITION> -N 1 -c 8 --mem=16g --wrap="STAR --runMode alignReads --runThreadN 8 --readFilesCommand zcat --outFilterMultimapNmax 10 --alignIntronMin 61 --alignIntronMax 10000 --genomeDir reference/ --readFilesIn trim/KRWTD3_1_val_1.fq.gz trim/KRWTD3_2_val_2.fq.gz --outSAMtype BAM SortedByCoordinate --outFileNamePrefix bam/KRWTD3.bam"

sbatch -A <YOURACCOUNT> -p <YOURPARTITION> -N 1 -c 8 --mem=16g --wrap="STAR --runMode alignReads --runThreadN 8 --readFilesCommand zcat --outFilterMultimapNmax 10 --alignIntronMin 61 --alignIntronMax 10000 --genomeDir reference/ --readFilesIn trim/KRWTD4_1_val_1.fq.gz trim/KRWTD4_2_val_2.fq.gz --outSAMtype BAM SortedByCoordinate --outFileNamePrefix bam/KRWTD4.bam"

sbatch -A <YOURACCOUNT> -p <YOURPARTITION> -N 1 -c 8 --mem=16g --wrap="STAR --runMode alignReads --runThreadN 8 --readFilesCommand zcat --outFilterMultimapNmax 10 --alignIntronMin 61 --alignIntronMax 10000 --genomeDir reference/ --readFilesIn trim/KRWTW1_1_val_1.fq.gz trim/KRWTW1_2_val_2.fq.gz --outSAMtype BAM SortedByCoordinate --outFileNamePrefix bam/KRWTW1.bam"

sbatch -A <YOURACCOUNT> -p <YOURPARTITION> -N 1 -c 8 --mem=16g --wrap="STAR --runMode alignReads --runThreadN 8 --readFilesCommand zcat --outFilterMultimapNmax 10 --alignIntronMin 61 --alignIntronMax 10000 --genomeDir reference/ --readFilesIn trim/KRWTW2_1_val_1.fq.gz trim/KRWTW2_2_val_2.fq.gz --outSAMtype BAM SortedByCoordinate --outFileNamePrefix bam/KRWTW2.bam"

sbatch -A <YOURACCOUNT> -p <YOURPARTITION> -N 1 -c 8 --mem=16g --wrap="STAR --runMode alignReads --runThreadN 8 --readFilesCommand zcat --outFilterMultimapNmax 10 --alignIntronMin 61 --alignIntronMax 10000 --genomeDir reference/ --readFilesIn trim/KRWTW3_1_val_1.fq.gz trim/KRWTW3_2_val_2.fq.gz --outSAMtype BAM SortedByCoordinate --outFileNamePrefix bam/KRWTW3.bam"

sbatch -A <YOURACCOUNT> -p <YOURPARTITION> -N 1 -c 8 --mem=16g --wrap="STAR --runMode alignReads --runThreadN 8 --readFilesCommand zcat --outFilterMultimapNmax 10 --alignIntronMin 61 --alignIntronMax 10000 --genomeDir reference/ --readFilesIn trim/KRWTW4_1_val_1.fq.gz trim/KRWTW4_2_val_2.fq.gz --outSAMtype BAM SortedByCoordinate --outFileNamePrefix bam/KRWTW4.bam"
```

### RNA-Seq Quantification
We recommend to use `featurecounts`
http://subread.sourceforge.net/
```bash
featureCounts -p -Q 10 -M -s 0 -T 24 -a reference/Athaliana_167_TAIR10.gene.gtf bam/KR24D1.bamAligned.sortedByCoord.out.bam bam/KR24D2.bamAligned.sortedByCoord.out.bam bam/KR24D3.bamAligned.sortedByCoord.out.bam bam/KR24D4.bamAligned.sortedByCoord.out.bam bam/KR24W1.bamAligned.sortedByCoord.out.bam bam/KR24W2.bamAligned.sortedByCoord.out.bam bam/KR24W3.bamAligned.sortedByCoord.out.bam bam/KR24W4.bamAligned.sortedByCoord.out.bam bam/KRWTD1.bamAligned.sortedByCoord.out.bam bam/KRWTD2.bamAligned.sortedByCoord.out.bam bam/KRWTD3.bamAligned.sortedByCoord.out.bam bam/KRWTD4.bamAligned.sortedByCoord.out.bam bam/KRWTW1.bamAligned.sortedByCoord.out.bam bam/KRWTW2.bamAligned.sortedByCoord.out.bam bam/KRWTW3.bamAligned.sortedByCoord.out.bam bam/KRWTW4.bamAligned.sortedByCoord.out.bam -o featureCount.cnt
```

### FeatureCounts output formating
```bash
cut -f1,6-  featureCount.cnt  | egrep -v "#" >> featureCount.cnt_for_tpm
cut -f1,7-  featureCount.cnt  | egrep -v "#" | sed 's/\.bamAligned\.sortedByCoord\.out\.bam//g; s/\.TAIR10//g' >> featureCount.cnt_for_DEG
```

### TPM calculation
Please download following script
```bash
wget https://pastebin.com/raw/3iv87VJe  -O tpm_raw_exp_calculator.py
```
Calculate raw count to TPM
```bash
chmod 775 tpm_raw_exp_calculator.py
tpm_raw_exp_calculator.py -count featureCount.cnt_for_tpm
```

### DEG calculation
We recommand to use R and DESeq2 Packages.
http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html

Please download below code

```bash
wget https://pastebin.com/LBvFdWFr -O rnaseq_plot_funcs.R
chmod 775 rnaseq_plot_funcs.R
```
Pleas use `R`

```R

if (! require(DESeq2)) {
   source("https://bioconductor.org/biocLite.R")
   biocLite("DESeq2")
   library(DESeq2)
}

data = read.table("featureCount.cnt_for_DEG", header=T, row.names=1, com='')
col_ordering = c(1,2,3,4,5,6,7,8)
rnaseqMatrix = data[,col_ordering]
rnaseqMatrix = round(rnaseqMatrix)
rnaseqMatrix = rnaseqMatrix[rowSums(cpm(rnaseqMatrix) > 1) >= 2,]
conditions = data.frame(conditions=factor(c(rep("KR24D", 4), rep("KR24W", 4))))
rownames(conditions) = colnames(rnaseqMatrix)
ddsFullCountTable <- DESeqDataSetFromMatrix(
    countData = rnaseqMatrix,
    colData = conditions,
    design = ~ conditions)
dds = DESeq(ddsFullCountTable)
contrast=c("conditions","KR24D","KR24W")
res = results(dds, contrast)
baseMeanA <- rowMeans(counts(dds, normalized=TRUE)[,colData(dds)$conditions == "KR24D"])
baseMeanB <- rowMeans(counts(dds, normalized=TRUE)[,colData(dds)$conditions == "KR24W"])
res = cbind(baseMeanA, baseMeanB, as.data.frame(res))
res = cbind(sampleA="KR24D", sampleB="KR24W", as.data.frame(res))
res$padj[is.na(res$padj)]  <- 1
res = as.data.frame(res[order(res$pvalue),])
write.table(res, file='featureCount.cnt_for_DEG.KR24D_vs_KR24W.DESeq2.DE_results', sep='        ', quote=FALSE)
write.table(rnaseqMatrix, file='featureCount.cnt_for_DEG.KR24D_vs_KR24W.DESeq2.count_matrix', sep='     ', quote=FALSE)
source("rnaseq_plot_funcs.R")
pdf("featureCount.cnt_for_DEG.KR24D_vs_KR24W.DESeq2.DE_results.MA_n_Volcano.pdf")
plot_MA_and_Volcano(rownames(res), log2(res$baseMean+1), res$log2FoldChange, res$padj)
dev.off()



data = read.table("featureCount.cnt_for_DEG", header=T, row.names=1, com='')
col_ordering = c(1,2,3,4,9,10,11,12)
rnaseqMatrix = data[,col_ordering]
rnaseqMatrix = round(rnaseqMatrix)
rnaseqMatrix = rnaseqMatrix[rowSums(cpm(rnaseqMatrix) > 1) >= 2,]
conditions = data.frame(conditions=factor(c(rep("KR24D", 4), rep("KRWTD", 4))))
rownames(conditions) = colnames(rnaseqMatrix)
ddsFullCountTable <- DESeqDataSetFromMatrix(
    countData = rnaseqMatrix,
    colData = conditions,
    design = ~ conditions)
dds = DESeq(ddsFullCountTable)
contrast=c("conditions","KR24D","KRWTD")
res = results(dds, contrast)
baseMeanA <- rowMeans(counts(dds, normalized=TRUE)[,colData(dds)$conditions == "KR24D"])
baseMeanB <- rowMeans(counts(dds, normalized=TRUE)[,colData(dds)$conditions == "KRWTD"])
res = cbind(baseMeanA, baseMeanB, as.data.frame(res))
res = cbind(sampleA="KR24D", sampleB="KRWTD", as.data.frame(res))
res$padj[is.na(res$padj)]  <- 1
res = as.data.frame(res[order(res$pvalue),])
write.table(res, file='featureCount.cnt_for_DEG.KR24D_vs_KRWTD.DESeq2.DE_results', sep='        ', quote=FALSE)
write.table(rnaseqMatrix, file='featureCount.cnt_for_DEG.KR24D_vs_KRWTD.DESeq2.count_matrix', sep='     ', quote=FALSE)
source("rnaseq_plot_funcs.R")
pdf("featureCount.cnt_for_DEG.KR24D_vs_KRWTD.DESeq2.DE_results.MA_n_Volcano.pdf")
plot_MA_and_Volcano(rownames(res), log2(res$baseMean+1), res$log2FoldChange, res$padj)
dev.off()

if (! require(DESeq2)) {
   source("https://bioconductor.org/biocLite.R")
   biocLite("DESeq2")
   library(DESeq2)
}

data = read.table("featureCount.cnt_for_DEG", header=T, row.names=1, com='')
col_ordering = c(5,6,7,8,13,14,15,16)
rnaseqMatrix = data[,col_ordering]
rnaseqMatrix = round(rnaseqMatrix)
rnaseqMatrix = rnaseqMatrix[rowSums(cpm(rnaseqMatrix) > 1) >= 2,]
conditions = data.frame(conditions=factor(c(rep("KR24W", 4), rep("KRWTW", 4))))
rownames(conditions) = colnames(rnaseqMatrix)
ddsFullCountTable <- DESeqDataSetFromMatrix(
    countData = rnaseqMatrix,
    colData = conditions,
    design = ~ conditions)
dds = DESeq(ddsFullCountTable)
contrast=c("conditions","KR24W","KRWTW")
res = results(dds, contrast)
baseMeanA <- rowMeans(counts(dds, normalized=TRUE)[,colData(dds)$conditions == "KR24W"])
baseMeanB <- rowMeans(counts(dds, normalized=TRUE)[,colData(dds)$conditions == "KRWTW"])
res = cbind(baseMeanA, baseMeanB, as.data.frame(res))
res = cbind(sampleA="KR24W", sampleB="KRWTW", as.data.frame(res))
res$padj[is.na(res$padj)]  <- 1
res = as.data.frame(res[order(res$pvalue),])
write.table(res, file='featureCount.cnt_for_DEG.KR24W_vs_KRWTW.DESeq2.DE_results', sep='        ', quote=FALSE)
write.table(rnaseqMatrix, file='featureCount.cnt_for_DEG.KR24W_vs_KRWTW.DESeq2.count_matrix', sep='     ', quote=FALSE)
source("rnaseq_plot_funcs.R")
pdf("featureCount.cnt_for_DEG.KR24W_vs_KRWTW.DESeq2.DE_results.MA_n_Volcano.pdf")
plot_MA_and_Volcano(rownames(res), log2(res$baseMean+1), res$log2FoldChange, res$padj)
dev.off()


if (! require(DESeq2)) {
   source("https://bioconductor.org/biocLite.R")
   biocLite("DESeq2")
   library(DESeq2)
}

data = read.table("featureCount.cnt_for_DEG", header=T, row.names=1, com='')
col_ordering = c(9,10,11,12,13,14,15,16)
rnaseqMatrix = data[,col_ordering]
rnaseqMatrix = round(rnaseqMatrix)
rnaseqMatrix = rnaseqMatrix[rowSums(cpm(rnaseqMatrix) > 1) >= 2,]
conditions = data.frame(conditions=factor(c(rep("KRWTD", 4), rep("KRWTW", 4))))
rownames(conditions) = colnames(rnaseqMatrix)
ddsFullCountTable <- DESeqDataSetFromMatrix(
    countData = rnaseqMatrix,
    colData = conditions,
    design = ~ conditions)
dds = DESeq(ddsFullCountTable)
contrast=c("conditions","KRWTD","KRWTW")
res = results(dds, contrast)
baseMeanA <- rowMeans(counts(dds, normalized=TRUE)[,colData(dds)$conditions == "KRWTD"])
baseMeanB <- rowMeans(counts(dds, normalized=TRUE)[,colData(dds)$conditions == "KRWTW"])
res = cbind(baseMeanA, baseMeanB, as.data.frame(res))
res = cbind(sampleA="KRWTD", sampleB="KRWTW", as.data.frame(res))
res$padj[is.na(res$padj)]  <- 1
res = as.data.frame(res[order(res$pvalue),])
write.table(res, file='featureCount.cnt_for_DEG.KRWTD_vs_KRWTW.DESeq2.DE_results', sep='        ', quote=FALSE)
write.table(rnaseqMatrix, file='featureCount.cnt_for_DEG.KRWTD_vs_KRWTW.DESeq2.count_matrix', sep='     ', quote=FALSE)
source("rnaseq_plot_funcs.R")
pdf("featureCount.cnt_for_DEG.KRWTD_vs_KRWTW.DESeq2.DE_results.MA_n_Volcano.pdf")
plot_MA_and_Volcano(rownames(res), log2(res$baseMean+1), res$log2FoldChange, res$padj)
dev.off()
```


###DEG filteration
Please download `featureCount.cnt_for_*.DESeq2.DE_results` files and open at the Excel to filter it out, based on your creteria.



