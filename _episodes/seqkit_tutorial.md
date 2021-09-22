---
layout: page
title: seqkit_tutorial
published: true
---


# Tutorial

## Some manipulations on big genomes

A script [memusg](https://github.com/shenwei356/memusg) is
used to check the peek memory usage of seqkit. Usage: `memusg [-t] command`.

1. Human genome

        $ seqkit stat hsa.fa
        file    format  type  num_seqs        sum_len  min_len       avg_len      max_len
        hsa.fa  FASTA   DNA        194  3,099,750,718      970  15,978,096.5  248,956,422

1. Build FASTA index (***optional***, when using flag `-2` (`--two-pass`),
   some commands will automaticlly build it).
   For some commands, including `subseq`, `split`, `sort` and `shuffle`,
   when input files are (plain or gzipped) FASTA files or stdin,
   FASTA index would be optional used for
   rapid acccess of sequences and reducing memory occupation.
   ***ATTENTION***: the `.seqkit.fai` file created by SeqKit is a little different from .fai file
   created by samtools. SeqKit uses full sequence head instead of just ID as key.

        $ memusg -t seqkit faidx --id-regexp "^(.+)$"  hsa.fa -o hsa.fa.seqkit.fai

        elapsed time: 10.011s
        peak rss: 177.21 MB

    Create common .fai file:

        $ memusg -t seqkit faidx hsa.fa -o hsa.fa.fai2

        elapsed time: 10.454s
        peak rss: 172.82 MB


    Performance of samtools:

        $ memusg -t samtools faidx hsa.fa

        elapsed time: 9.574s
        peak rss: 1.45 MB

    Exactly same content:

        $ md5sum hsa.fa.fai*
        21e0c25b4d817d1c19ee8107191b9b31  hsa.fa.fai
        21e0c25b4d817d1c19ee8107191b9b31  hsa.fa.fai2

1. Sorting by sequence length

        $ memusg -t seqkit sort --by-length --reverse --two-pass hsa.fa > hsa.sorted.fa
        [INFO] create and read FASTA index ...
        [INFO] read sequence IDs and lengths from FASTA index ...
        [INFO] 194 sequences loaded
        [INFO] sorting ...
        [INFO] output ...

        elapsed time: 4.892s
        peak rss: 500.15 MB

    Detail:

        $ seqkit fx2tab --length hsa.sorted.fa --name --only-id | cut -f 1,4 | more
        1       248956422
        2       242193529
        3       198295559
        4       190214555
        5       181538259
        6       170805979
        7       159345973
        X       156040895
        8       145138636
        9       138394717
        11      135086622
        10      133797422
        12      133275309
        13      114364328
        14      107043718
        15      101991189
        16      90338345
        17      83257441
        18      80373285
        20      64444167
        19      58617616
        Y       57227415
        22      50818468
        21      46709983
        KI270728.1      1872759
        KI270727.1      448248
        ...

        real    0m10.697s
        user    0m11.153s
        sys     0m0.917s

1. Shuffling sequences

        $ memusg -t seqkit shuffle hsa.fa --two-pass > hsa.shuffled.fa
        [INFO] create and read FASTA index ...
        [INFO] read sequence IDs from FASTA index ...
        [INFO] 194 sequences loaded
        [INFO] shuffle ...
        [INFO] output ...

        elapsed time: 6.632s
        peak rss: 528.3 MB


1. Spliting into files with single sequence

        $ memusg -t seqkit split --by-id hsa.fa --two-pass
        [INFO] split by ID. idRegexp: ^([^\s]+)\s?
        [INFO] create and read FASTA index ...
        [INFO] read sequence IDs from FASTA index ...
        [INFO] 194 sequences loaded
        [INFO] write 1 sequences to file: hsa.id_KI270743.1.fa
        [INFO] write 1 sequences to file: hsa.id_KI270706.1.fa
        [INFO] write 1 sequences to file: hsa.id_KI270717.1.fa
        [INFO] write 1 sequences to file: hsa.id_KI270718.1.fa
        [INFO] write 1 sequences to file: hsa.id_KI270468.1.fa
        ...

        elapsed time: 18.807s
        peak rss: 1.36 GB

1. Geting subsequence of some chromesomes

        $ memusg -t seqkit subseq -r 1:10 --chr X --chr Y  hsa.fa
        >X_1-10 X dna_sm:chromosome chromosome:GRCh38:X:1:156040895:1 REF
        nnnnnnnnnn
        >Y_1-10 Y dna_sm:chromosome chromosome:GRCh38:Y:2781480:56887902:1 REF
        NNNNNNNNNN

        elapsed time: 1.276s
        peak rss: 640.92 MB


1. Geting CDS sequence of chr 1 by GTF files

        $ memusg -t seqkit subseq --gtf Homo_sapiens.GRCh38.84.gtf.gz --chr X --feature cds  hsa.fa > chrX.gtf.cds.fa
        [INFO] read GTF file ...
        [INFO] 22420 GTF features loaded

        elapsed time: 8.643s
        peak rss: 846.14 MB


## Remove contaminated reads

1. Mapping with reads on some potential contaminate genomes, and get the reads IDs list.

        $ wc -l contaminate.list
        244 contaminate.list

        $ head -n 2 contaminate.list
        HWI-D00523:240:HF3WGBCXX:1:1101:2574:2226
        HWI-D00523:240:HF3WGBCXX:1:1101:12616:2205

1. Remove contaminated reads

        $ seqkit grep -f contaminate.list -v reads_1.fq.gz -o reads_1.clean.fq.gz
        $ seqkit grep -f contaminate.list -v reads_2.fq.gz -o reads_2.clean.fq.gz

        $ seqkit stat *.fq.gz
        file                  seq_format   seq_type   num_seqs   min_len   avg_len   max_len
        reads_1.clean.fq.gz   FASTQ        DNA           2,256       226       227       229
        reads_1.fq.gz         FASTQ        DNA           2,500       226       227       229
        reads_2.clean.fq.gz   FASTQ        DNA           2,256       223       224       225
        reads_2.fq.gz         FASTQ        DNA           2,500       223       224       225



## Handling of aligned sequences

1. Some mock sequences (usually they will be much longer)

        $ cat seqs.fa
        >seq1
        ACAACGTCTACTTACGTTGCATCGTCATGCTGCATTACGTAGTCTGATGATG
        >seq2
        ACACCGTCTACTTTCATGCTGCATTACGTAGTCTGATGATG
        >seq3
        ACAACGTCTACTTACGTTGCATCGTCATGCTGCACTGATGATG
        >seq4
        ACAACGTCTACTTACGTTGCATCTTCGGTCATGCTGCATTACGTAGTCTGATGATG

1. Run multiple sequence alignment (clustalo)

        clustalo -i seqs.fa -o seqs.msa.fa --force --outfmt fasta --threads=4

1. Convert FASTA format to tabular format.

        $ seqkit fx2tab seqs.msa.fa
        seq1    ACAACGTCTACTTACGTTGCAT----CGTCATGCTGCATTACGTAGTCTGATGATG
        seq2    ---------------ACACCGTCTACTTTCATGCTGCATTACGTAGTCTGATGATG
        seq3    ACAACGTCTACTTACGTTGCATCGTCATGCTGCACTGATGATG-------------
        seq4    ACAACGTCTACTTACGTTGCATCTTCGGTCATGCTGCATTACGTAGTCTGATGATG

    or

        $ seqkit fx2tab seqs.msa.fa | cut -f 2
        ACAACGTCTACTTACGTTGCAT----CGTCATGCTGCATTACGTAGTCTGATGATG
        ---------------ACACCGTCTACTTTCATGCTGCATTACGTAGTCTGATGATG
        ACAACGTCTACTTACGTTGCATCGTCATGCTGCACTGATGATG-------------
        ACAACGTCTACTTACGTTGCATCTTCGGTCATGCTGCATTACGTAGTCTGATGATG

    For me, it's useful when 1) manually assembling Sanger sequencing result,
    2) designing site specific PCR primers.


1. Remove gaps

        $ seqkit seq seqs.msa.fa -g
        >seq1
        ACAACGTCTACTTACGTTGCATCGTCATGCTGCATTACGTAGTCTGATGATG
        >seq2
        ACACCGTCTACTTTCATGCTGCATTACGTAGTCTGATGATG
        >seq3
        ACAACGTCTACTTACGTTGCATCGTCATGCTGCACTGATGATG
        >seq4
        ACAACGTCTACTTACGTTGCATCTTCGGTCATGCTGCATTACGTAGTCTGATGATG



## Play with miRNA hairpins

### Dataset

[`hairpin.fa.gz`](ftp://mirbase.org/pub/mirbase/21/hairpin.fa.gz)
from [The miRBase Sequence Database -- Release 21](ftp://mirbase.org/pub/mirbase/21/)


### Quick glance

1. Sequence number

        $ seqkit stat hairpin.fa.gz
        file           format  type  num_seqs    sum_len  min_len  avg_len  max_len
        hairpin.fa.gz  FASTA   RNA     28,645  2,949,871       39      103    2,354

1. First 10 bases

        $ zcat hairpin.fa.gz | seqkit subseq -r 1:10 | seqkit sort -s | seqkit seq -s | head -n 10
        AAAAAAAAAA
        AAAAAAAAAA
        AAAAAAAAAG
        AAAAAAAAAG
        AAAAAAAAAG
        AAAAAAAAAU
        AAAAAAAAGG
        AAAAAAACAU
        AAAAAAACGA
        AAAAAAAUUA

    hmm, nothing special, non-coding RNA~


### Repeated hairpin sequences

We may want to check how may identical hairpins among different species there are.
`seqkit rmdup` could remove duplicated sequences by sequence content,
and save the replicates to another file (here is `duplicated.fa.gz`),
as well as replicating details (`duplicated.detail.txt`,
1th column is the repeated number,
2nd column contains sequence IDs seperated by comma).

    $ seqkit rmdup -s -i hairpin.fa.gz -o clean.fa.gz -d duplicated.fa.gz -D duplicated.detail.txt

    $ head -n 5 duplicated.detail.txt
    18      dre-mir-430c-1, dre-mir-430c-2, dre-mir-430c-3, dre-mir-430c-4, dre-mir-430c-5, dre-mir-430c-6, dre-mir-430c-7, dre-mir-430c-8, dre-mir-430c-9, dre-mir-430c-10, dre-mir-430c-11, dre-mir-430c-12, dre-mir-430c-13, dre-mir-430c-14, dre-mir-430c-15, dre-mir-430c-16, dre-mir-430c-17, dre-mir-430c-18
    16      hsa-mir-29b-2, mmu-mir-29b-2, rno-mir-29b-2, ptr-mir-29b-2, ggo-mir-29b-2, ppy-mir-29b-2, sla-mir-29b, mne-mir-29b, ppa-mir-29b-2, bta-mir-29b-2, mml-mir-29b-2, eca-mir-29b-2, aja-mir-29b, oar-mir-29b-1, oar-mir-29b-2, rno-mir-29b-3
    15      dme-mir-125, dps-mir-125, dan-mir-125, der-mir-125, dgr-mir-125-1, dgr-mir-125-2, dmo-mir-125, dpe-mir-125-2, dpe-mir-125-1, dpe-mir-125-3, dse-mir-125, dsi-mir-125, dvi-mir-125, dwi-mir-125, dya-mir-125
    13      hsa-mir-19b-1, ggo-mir-19b-1, age-mir-19b-1, ppa-mir-19b-1, ppy-mir-19b-1, ptr-mir-19b-1, mml-mir-19b-1, sla-mir-19b-1, lla-mir-19b-1, mne-mir-19b-1, bta-mir-19b, oar-mir-19b, chi-mir-19b
    13      hsa-mir-20a, ssc-mir-20a, ggo-mir-20a, age-mir-20, ppa-mir-20, ppy-mir-20a, ptr-mir-20a, mml-mir-20a, sla-mir-20, lla-mir-20, mne-mir-20, bta-mir-20a, eca-mir-20a

The result shows the most conserved miRNAs among different species,
`mir-29b`, `mir-125`, `mir-19b-1` and `mir-20a`.
And the `dre-miR-430c` has the most multicopies in *Danio rerio*.

### Hairpins in different species

1. Before spliting by species, let's take a look at the sequence names.

        $ seqkit seq hairpin.fa.gz -n | head -n 3
        cel-let-7 MI0000001 Caenorhabditis elegans let-7 stem-loop
        cel-lin-4 MI0000002 Caenorhabditis elegans lin-4 stem-loop
        cel-mir-1 MI0000003 Caenorhabditis elegans miR-1 stem-loop

    The first three letters (e.g. `cel`) are the abbreviation of species names.
    So we could split hairpins by the first letters by defining custom
    sequence ID parsing regular expression `^([\w]+)\-`.

    By default, `seqkit` takes the first non-space letters as sequence ID.
    For example,

    |   FASTA head                                                  |     ID                                            |
    |:--------------------------------------------------------------|:--------------------------------------------------|
    | >123456 gene name                                             | 123456                                            |
    | >longname                                                     | longname                                          |
    | >gi&#124;110645304&#124;ref&#124;NC_002516.2&#124; Pseudomona | gi&#124;110645304&#124;ref&#124;NC_002516.2&#124; |

    But for some sequences from NCBI,
    e.g. `>gi|110645304|ref|NC_002516.2| Pseudomona`, the ID is `NC_002516.2`.
    In this case, we could set sequence ID parsing regular expression by flag
    `--id-regexp "\|([^\|]+)\| "` or just use flag `--id-ncbi`. If you want
    the `gi` number, then use `--id-regexp "^gi\|([^\|]+)\|"`.

1. Split sequences by species.
A custom ID parsing regular expression is used, `^([\w]+)\-`.

        $ seqkit split hairpin.fa.gz -i --id-regexp "^([\w]+)\-" --two-pass

    ***To reduce memory usage when splitting big file, we should always use flag `--two-pass`***

2. Species with most miRNA hairpins. Third column is the sequences number.

        $ cd hairpin.fa.gz.split/;
        $ seqkit stat hairpin.id_* | csvtk space2tab | csvtk -t sort -k num_seqs:nr | csvtk -t pretty| more
        file                     format   type   num_seqs   sum_len   min_len   avg_len   max_len
        hairpin.id_hsa.fasta     FASTA    RNA    1,881      154,242   82        82        82
        hairpin.id_mmu.fasta     FASTA    RNA    1,193      107,370   90        90        90
        hairpin.id_bta.fasta     FASTA    RNA    808        61,408    76        76        76
        hairpin.id_gga.fasta     FASTA    RNA    740        42,180    57        57        57
        hairpin.id_eca.fasta     FASTA    RNA    715        89,375    125       125       125
        hairpin.id_mtr.fasta     FASTA    RNA    672        231,840   345       345       345

    Here, a CSV/TSV tool [csvtk](https://github.com/shenwei356/csvtk)
    is used to sort and view the result.

For human miRNA hairpins

1. Length distribution.
 `seqkit fx2tab` could show extra information like sequence length, GC content.
 A distribution ploting script is used, (
 [plot_distribution.py](https://github.com/shenwei356/bio_scripts/blob/master/plot/plot_distribution.py) )

        $ seqkit fx2tab hairpin.id_hsa.fa.gz -l | cut -f 3  | plot_distribution.py -o hairpin.id_hsa.fa.gz.lendist.png

    ![hairpin.id_hsa.fa.gz.lendist.png](/files/hairpin/hairpin.id_hsa.fa.gz.lendist.png)


## Bacteria genome

### Dataset

[Pseudomonas aeruginosa PAO1](http://www.ncbi.nlm.nih.gov/nuccore/110645304),
files:

-  Genbank file [`PAO1.gb`](/files/PAO1/PAO1.gb)
-  Genome FASTA file [`PAO1.fasta`](/files/PAO1/PAO1.fasta)
-  GTF file [`PAO1.gtf`](/files/PAO1/PAO1.gtf) was created with [`extract_features_from_genbank_file.py`](https://github.com/shenwei356/bio_scripts/blob/master/file_formats/extract_features_from_genbank_file.py), by

        extract_features_from_genbank_file.py  PAO1.gb -t . -f gtf > PAO1.gtf


### Motif distribution

Motifs

    $ cat motifs.fa
    >GTAGCGS
    GTAGCGS
    >GGWGKTCG
    GGWGKTCG

1. Sliding. Remember flag `--id-ncbi`, do you?
  By the way, do not be scared by the long flag `--circle-genome`, `--step`
  and so on. They have short ones, `-c`, `-s`

        $ seqkit sliding --id-ncbi --circular-genome --step 20000 --window 200000 PAO1.fasta -o PAO1.fasta.sliding.fa

        $ seqkit stat PAO1.fasta.sliding.fa
        file                   format  type  num_seqs     sum_len  min_len  avg_len  max_len
        PAO1.fasta.sliding.fa  FASTA   DNA        314  62,800,000  200,000  200,000  200,000

1. Locating motifs

        $ seqkit locate --id-ncbi --ignore-case --degenerate --pattern-file motifs.fa  PAO1.fasta.sliding.fa -o  PAO1.fasta.sliding.fa.motifs.tsv

1. Ploting distribution ([plot_motif_distribution.R](/files/PAO1/plot_motif_distribution.R))

        # preproccess
        $ perl -ne 'if (/_sliding:(\d+)-(\d+)\t(.+)/) {$loc= $1 + 100000; print "$loc\t$3\n";} else {print}' PAO1.fasta.sliding.fa.motifs.tsv  > PAO1.fasta.sliding.fa.motifs.tsv2

        # plot
        $ ./plot_motif_distribution.R

    Result

    ![motif_distribution.png](files/PAO1/motif_distribution.png)


### Find multicopy genes

1. Get all CDS sequences

        $ seqkit subseq --id-ncbi --gtf PAO1.gtf --feature cds PAO1.fasta -o PAO1.cds.fasta

        $ seqkit stat *.fasta
        file            format  type  num_seqs    sum_len    min_len    avg_len    max_len
        PAO1.cds.fasta  FASTA   DNA      5,572  5,593,306         72    1,003.8     16,884
        PAO1.fasta      FASTA   DNA          1  6,264,404  6,264,404  6,264,404  6,264,404

1. Get duplicated sequences

        $ seqkit rmdup --by-seq --ignore-case PAO1.cds.fasta -o PAO1.cds.uniq.fasta --dup-seqs-file PAO1.cds.dup.fasta --dup-num-file PAO1.cds.dup.text

        $ cat PAO1.cds.dup.text
        6       NC_002516.2_500104:501120:-, NC_002516.2_2556948:2557964:+, NC_002516.2_3043750:3044766:-, NC_002516.2_3842274:3843290:-, NC_002516.2_4473623:4474639:+, NC_002516.2_5382796:5383812:-
        2       NC_002516.2_2073555:2075438:+, NC_002516.2_4716660:4718543:+
        2       NC_002516.2_2072935:2073558:+, NC_002516.2_4716040:4716663:+
        2       NC_002516.2_2075452:2076288:+, NC_002516.2_4718557:4719393:+

### Flanking sequences

1. Get CDS and 1000 bp upstream sequence

        $ seqkit subseq --id-ncbi --gtf PAO1.gtf --feature cds PAO1.fasta --up-stream 1000

1. Get 1000 bp upstream sequence of CDS, *NOT* including CDS.

        $ seqkit subseq --id-ncbi --gtf PAO1.gtf --feature cds PAO1.fasta --up-stream 1000 --only-flank

<div id="disqus_thread"></div>
<script>
/**
* RECOMMENDED CONFIGURATION VARIABLES: EDIT AND UNCOMMENT THE SECTION BELOW TO INSERT DYNAMIC VALUES FROM YOUR PLATFORM OR CMS.
* LEARN WHY DEFINING THESE VARIABLES IS IMPORTANT: https://disqus.com/admin/universalcode/#configuration-variables
*/
/*
var disqus_config = function () {
this.page.url = PAGE_URL; // Replace PAGE_URL with your page's canonical URL variable
this.page.identifier = PAGE_IDENTIFIER; // Replace PAGE_IDENTIFIER with your page's unique identifier variable
};
*/
(function() { // DON'T EDIT BELOW THIS LINE
var d = document, s = d.createElement('script');

s.src = '//fastakit.disqus.com/embed.js';

s.setAttribute('data-timestamp', +new Date());
(d.head || d.body).appendChild(s);
})();
</script>
<noscript>Please enable JavaScript to view the <a href="https://disqus.com/?ref_noscript" rel="nofollow">comments powered by Disqus.</a></noscript>



    # Benchmark


## Softwares

1. [seqkit](https://github.com/shenwei356/seqkit). (Go).
   Version [v0.3.1.1](https://github.com/shenwei356/seqkit/releases/tag/v0.3.1.1).
   Compiled with Go 1.7rc5.
1. [fasta_utilities](https://github.com/jimhester/fasta_utilities). (Perl).
   Version [3dcc0bc](https://github.com/jimhester/fasta_utilities/tree/3dcc0bc6bf1e97839476221c26984b1789482579).
   Lots of dependencies to install.
1. [fastx_toolkit](http://hannonlab.cshl.edu/fastx_toolkit/). (Perl).
   Version [0.0.13](http://hannonlab.cshl.edu/fastx_toolkit/fastx_toolkit_0.0.13_binaries_Linux_2.6_amd64.tar.bz2).
   Can't handle multi-line FASTA files.
1. [seqmagick](http://seqmagick.readthedocs.io/en/latest/index.html#installation). (Python).
   Version 0.6.1
1. [seqtk](https://github.com/lh3/seqtk). (C).
   Version [1.1-r92-dirty](https://github.com/lh3/seqtk/tree/fb85aad4ce1fc7b3d4543623418a1ae88fe1cea6).

Not used:

1. [pyfaidx](https://github.com/mdshw5/pyfaidx). (Python).
   Version [0.4.7.1](https://pypi.python.org/packages/source/p/pyfaidx/pyfaidx-0.4.7.1.tar.gz#md5=f33604a3550c2fa115ac7d33b952127d). *Not used, because it exhausted my memory (10G) when computing reverse-complement on a 5GB fasta file of 250 bp.*

A Python script [memusg](https://github.com/shenwei356/memusg) was used
to compute running time and peak memory usage of a process.

## Features

Categories          |Features               |seqkit  |fasta_utilities|fastx_toolkit|pyfaidx|seqmagick|seqtk
:-------------------|:----------------------|:------:|:-------------:|:-----------:|:-----:|:-------:|:---:
**Formats support** |Multi-line FASTA       |Yes     |Yes            |--           |Yes    |Yes      |Yes
                    |FASTQ                  |Yes     |Yes            |Yes          |--     |Yes      |Yes
                    |Multi-line  FASTQ      |Yes     |Yes            |--           |--     |Yes      |Yes
                    |Validating sequences   |Yes     |--             |Yes          |Yes    |--       |--
                    |Supporting RNA         |Yes     |Yes            |--           |--     |Yes      |Yes
**Functions**       |Searching by motifs    |Yes     |Yes            |--           |--     |Yes      |--
                    |Sampling               |Yes     |--             |--           |--     |Yes      |Yes
                    |Extracting sub-sequence|Yes     |Yes            |--           |Yes    |Yes      |Yes
                    |Removing duplicates    |Yes     |--             |--           |--     |Partly   |--
                    |Splitting              |Yes     |Yes            |--           |Partly |--       |--
                    |Splitting by seq       |Yes     |--             |Yes          |Yes    |--       |--
                    |Shuffling              |Yes     |--             |--           |--     |--       |--
                    |Sorting                |Yes     |Yes            |--           |--     |Yes      |--
                    |Locating motifs        |Yes     |--             |--           |--     |--       |--
                    |Common sequences       |Yes     |--             |--           |--     |--       |--
                    |Cleaning bases         |Yes     |Yes            |Yes          |Yes    |--       |--
                    |Transcription          |Yes     |Yes            |Yes          |Yes    |Yes      |Yes
                    |Translation            |--      |Yes            |Yes          |Yes    |Yes      |--
                    |Filtering by size      |Indirect|Yes            |--           |Yes    |Yes      |--
                    |Renaming header        |Yes     |Yes            |--           |--     |Yes      |Yes
**Other features**  |Cross-platform         |Yes     |Partly         |Partly       |Yes    |Yes      |Yes
                    |Reading STDIN          |Yes     |Yes            |Yes          |--     |Yes      |Yes
                    |Reading gzipped file   |Yes     |Yes            |--           |--     |Yes      |Yes
                    |Writing gzip file      |Yes     |--             |--           |--     |Yes      |--

**Note 2**: See [usage](http://shenwei356.github.io/seqkit/usage/) for detailed options of seqkit.

## Datasets

All test data is available here: [seqkit-benchmark-data.tar.gz](http://app.shenwei.me/data/seqkit/seqkit-benchmark-data.tar.gz)  (2.2G)

### dataset_A.fa - large number of short sequences

Dataset A is reference genomes DNA sequences of gastrointestinal tract from
[NIH Human Microbiome Project](http://hmpdacc.org/):
[`Gastrointestinal_tract.nuc.fsa`](http://downloads.hmpdacc.org/data/reference_genomes/body_sites/Gastrointestinal_tract.nuc.fsa) (FASTA format, ~2.7G).

### dataset_B.fa - small number of large sequences

Dataset B is Human genome from [ensembl](http://uswest.ensembl.org/info/data/ftp/index.html).

- Genome DNA:  [`Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz`](ftp://ftp.ensembl.org/pub/release-84/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz) (Gzipped FASTA file, ~900M)
. Decompress it and rename to dataset_B.fa (~2.9G).
- GTF file:  [`Homo_sapiens.GRCh38.84.gtf.gz`](ftp://ftp.ensembl.org/pub/release-84/gtf/homo_sapiens/Homo_sapiens.GRCh38.84.gtf.gz) (~44M)
- BED file: `Homo_sapiens.GRCh38.84.bed.gz` was converted from `Homo_sapiens.GRCh38.84.gtf.gz` by  [`gtf2bed`](http://bedops.readthedocs.org/en/latest/content/reference/file-management/conversion/gtf2bed.html?highlight=gtf2bed)  with command

        $ zcat Homo_sapiens.GRCh38.84.gtf.gz | gtf2bed --do-not-sort | gzip -c > Homo_sapiens.GRCh38.84.bed.gz


### dataset_C.fq – Illumina single end reads (SE100)

Dataset C is Illumina single end (SE 100bp) reads file (~2.2G).


Summary

    $ seqkit stat *.fa
    file          format  type   num_seqs        sum_len  min_len       avg_len      max_len
    dataset_A.fa  FASTA   DNA      67,748  2,807,643,808       56      41,442.5    5,976,145
    dataset_B.fa  FASTA   DNA         194  3,099,750,718      970  15,978,096.5  248,956,422
    dataset_C.fq  FASTQ   DNA   9,186,045    918,604,500      100           100          100
    
### Sequence ID list

Parts of sequences IDs was sampled and shuffled from original data.
They were used in test of extracting sequences by ID list.

Commands:

    $ seqkit sample -p 0.3  dataset_A.fa | seqkit seq --name --only-id | shuf > ids_A.txt
    $ seqkit sample -p 0.3  dataset_B.fa | seqkit seq --name --only-id | shuf > ids_B.txt    
    $ seqkit sample -p 0.03 dataset_C.fq | seqkit seq --name --only-id | shuf > ids_C.txt

Numbers:

    $ wc -l ids*.txt
        20138 ids_A.txt
        58 ids_B.txt
    2754516 ids_C.txt

### BED file

Only BED data of chromosome 19 was used in test of subsequence with BED file:

    $ zcat Homo_sapiens.GRCh38.84.bed.gz | grep -E "^19" | gzip -c > chr19.bed.gz


## Platform

PC:

- CPU: Intel Core i5-3320M @ 2.60GHz, two cores/4 threads
- RAM: DDR3 1600MHz, 12GB
- SSD: SAMSUNG 850 EVO 250G, SATA-3
- OS: Fedora 24 (Scientific KDE spin),  Kernal: 4.6.4-301.fc24.x86_64

Softwares:

- Perl: perl 5, version 22, subversion 2 (v5.22.2) built for x86_64-linux-thread-multi
- Python: Python 2.7.11 (default, Jul 10 2016, 20:58:20) [GCC 6.1.1 20160621 (Red Hat 6.1.1-3)] on linux2

## Tests

Automatic benchmark and plotting scripts are available at:  [https://github.com/shenwei356/seqkit/tree/master/benchmark](https://github.com/shenwei356/seqkit/tree/master/benchmark).

All tests were repeated 3 times,
and average time and peak memory ware used for plotting.

All data were readed once before tests began to minimize the influence of page cache.

Output sequences of all softwares were not wrapped to fixed length.

### Test 1. Reverse Complement

`revcom_biogo` ([source](https://github.com/shenwei356/seqkit/blob/master/benchmark/revcom_biogo.go),
 [binary](https://github.com/shenwei356/seqkit/blob/master/benchmark/revcom_biogo?raw=true) ),
 a tool written in Golang (compiled with Go 1.6.3) using [biogo](https://github.com/biogo/biogo)
 (Version [7ebd71b](https://github.com/biogo/biogo/commit/7ebd71bd9afc52cdab7a7128467ae1a936b68958))
 package,
 was also used for comparison of FASTA file parsing performance.

*Note that some softwares (fasta_utilities and biogo) have different converting rules of computing complement sequence on ambiguous bases, there fore the results are different from others.*

[Commands](https://github.com/shenwei356/seqkit/blob/master/benchmark/run_benchmark_01_revcom.sh)

### Test 2. Extracting sequences by ID list

[Commands](https://github.com/shenwei356/seqkit/blob/master/benchmark/run_benchmark_02_exctact_by_id_list.sh)

### Test 3. Sampling by number

*Note that different softwares have different sampling strategies, the peak memory  depends on size of sampled sequences and the results may not be the same.*

[Commands](https://github.com/shenwei356/seqkit/blob/master/benchmark/run_benchmark_03_sampling.sh)

### Test 4. Removing duplicates by sequence content

[Commands](https://github.com/shenwei356/seqkit/blob/master/benchmark/run_benchmark_04_remove_duplicated_seqs_by_seq.sh)

### Test 5. Subsequence with BED file

[Commands](https://github.com/shenwei356/seqkit/blob/master/benchmark/run_benchmark_05_subseq_with_bed.sh)

## Results

seqkit version: v0.3.1.1

FASTA:

![benchmark-5tests.tsv.png](benchmark/benchmark.5tests.tsv.png)

FASTQ:

![benchmark-5tests.tsv.png](benchmark/benchmark.5tests.tsv.C.png)

### Test of multiple threads:

From the results, 2 threads/CPU is enough, so the default threads of seqkit is 2.

![benchmark-5tests.tsv.png](benchmark/seqkit_multi_threads/benchmark.5tests.tsv.png)

### Tests on different file sizes

Files are generated by replicating Human genome chr1 for N times.

![benchmark.seqkit.files_size.tsv.png](benchmark/seqkit_file_size/benchmark.seqkit.files_size.tsv.png)

<div id="disqus_thread"></div>
<script>
/**
* RECOMMENDED CONFIGURATION VARIABLES: EDIT AND UNCOMMENT THE SECTION BELOW TO INSERT DYNAMIC VALUES FROM YOUR PLATFORM OR CMS.
* LEARN WHY DEFINING THESE VARIABLES IS IMPORTANT: https://disqus.com/admin/universalcode/#configuration-variables
*/
/*
var disqus_config = function () {
this.page.url = PAGE_URL; // Replace PAGE_URL with your page's canonical URL variable
this.page.identifier = PAGE_IDENTIFIER; // Replace PAGE_IDENTIFIER with your page's unique identifier variable
};
*/
(function() { // DON'T EDIT BELOW THIS LINE
var d = document, s = d.createElement('script');

s.src = '//fastakit.disqus.com/embed.js';

s.setAttribute('data-timestamp', +new Date());
(d.head || d.body).appendChild(s);
})();
</script>
<noscript>Please enable JavaScript to view the <a href="https://disqus.com/?ref_noscript" rel="nofollow">comments powered by Disqus.</a></noscript>


# Usage and Examples

## Technical details and guides for use

### FASTA/Q format parsing

SeqKit uses author's lightweight and high-performance bioinformatics packages
[bio](https://github.com/shenwei356/bio) for FASTA/Q parsing,
which has [high performance](https://github.com/shenwei356/bio#fastaq-parsing)
close to the
famous C lib [klib](https://github.com/attractivechaos/klib/) ([kseq.h](https://github.com/attractivechaos/klib/blob/master/kseq.h)).

![](https://github.com/shenwei356/bio/raw/master/benchmark/benchmark.tsv.png)

### Sequence formats and types

SeqKit seamlessly support FASTA and FASTQ format.
Sequence format is automatically detected.
All subcommands except for `faidx` can handle both formats.
And only when some commands (`subseq`, `split`, `sort` and `shuffle`)
which utilise FASTA index to improve perfrmance for large files in two pass mode
(by flag `--two-pass`), only FASTA format is supported.


Sequence type (DNA/RNA/Protein) is automatically detected by leading subsequences
of the first sequences in file or STDIN. The length of the leading subsequences
is configurable by global flag `--alphabet-guess-seq-length` with default value
of 10000. If length of the sequences is less than that, whole sequences will
be checked.

### Sequence ID

By default, most softwares, including `seqkit`, take the leading non-space
letters as sequence identifier (ID). For example,

|   FASTA header                                                  |     ID                                            |
|:--------------------------------------------------------------|:--------------------------------------------------|
| >123456 gene name                                             | 123456                                            |
| >longname                                                     | longname                                          |
| >gi&#124;110645304&#124;ref&#124;NC_002516.2&#124; Pseudomona | gi&#124;110645304&#124;ref&#124;NC_002516.2&#124; |

But for some sequences from NCBI,
e.g. `>gi|110645304|ref|NC_002516.2| Pseudomona`, the ID is `NC_002516.2`.
In this case, we could set sequence ID parsing regular expression by global flag
`--id-regexp "\|([^\|]+)\| "` or just use flag `--id-ncbi`. If you want
the `gi` number, then use `--id-regexp "^gi\|([^\|]+)\|"`.

### FASTA index

For some commands, including `subseq`, `split`, `sort` and `shuffle`,
when input files are (plain or gzipped) FASTA files,
FASTA index would be optional used for
rapid access of sequences and reducing memory occupation.

ATTENTION: the `.seqkit.fai` file created by SeqKit is a little different from `.fai` file
created by `samtools`. SeqKit uses full sequence head instead of just ID as key.

### Parallelization of CPU intensive jobs

The validation of sequences bases and complement process of sequences
are parallelized for large sequences.

Parsing of line-based files, including BED/GFF file and ID list file are also parallelized.

The Parallelization is implemented by multiple goroutines in golang
 which are similar to but much
lighter weight than threads. The concurrency number is configurable with global
flag `-j` or `--threads` (default value: 1 for single-CPU PC, 2 for others).

### Memory occupation

Most of the subcommands do not read whole FASTA/Q records in to memory,
including `stat`, `fq2fa`, `fx2tab`, `tab2fx`, `grep`, `locate`, `replace`,
 `seq`, `sliding`, `subseq`.

Note that when using `subseq --gtf | --bed`, if the GTF/BED files are too
big, the memory usage will increase.
You could use `--chr` to specify chromesomes and `--feature` to limit features.

Some subcommands need to store sequences or heads in memory, but there are
strategy to reduce memory occupation, including `rmdup` and `common`.
When comparing with sequences, MD5 digest could be used to replace sequence by
flag `-m` (`--md5`).

Some subcommands could either read all records or read the files twice by flag
`-2` (`--two-pass`), including `sample`, `split`, `shuffle` and `sort`.
They use FASTA index for rapid acccess of sequences and reducing memory occupation.

### Reproducibility

Subcommands `sample` and `shuffle` use random function, random seed could be
given by flag `-s` (`--rand-seed`). This makes sure that sampling result could be
reproduced in different environments with same random seed.

## seqkit

Usage

```
SeqKit -- a cross-platform and ultrafast toolkit for FASTA/Q file manipulation

Version: 0.3.1.1

Author: Wei Shen <shenwei356@gmail.com>

Documents  : http://shenwei356.github.io/seqkit
Source code: https://github.com/shenwei356/seqkit

Usage:
  seqkit [command]

Available Commands:
  common      find common sequences of multiple files by id/name/sequence
  faidx       create FASTA index file
  fq2fa       covert FASTQ to FASTA
  fx2tab      covert FASTA/Q to tabular format (with length/GC content/GC skew)
  grep        search sequences by pattern(s) of name or sequence motifs
  head        print first N FASTA/Q records
  locate      locate subsequences/motifs
  rename      rename duplicated IDs
  replace     replace name/sequence by regular expression
  rmdup       remove duplicated sequences by id/name/sequence
  sample      sample sequences by number or proportion
  seq         transform sequences (revserse, complement, extract ID...)
  shuffle     shuffle sequences
  sliding     sliding sequences, circular genome supported
  sort        sort sequences by id/name/sequence/length
  split       split sequences into files by id/seq region/size/parts
  stat        simple statistics of FASTA files
  subseq      get subsequences by region/gtf/bed, including flanking sequences
  tab2fx      covert tabular format to FASTA/Q format
  version     print version information and check for update

Flags:
      --alphabet-guess-seq-length int   length of sequence prefix of the first FASTA record based on which seqkit guesses the sequence type (0 for whole seq) (default 10000)
      --id-ncbi                         FASTA head is NCBI-style, e.g. >gi|110645304|ref|NC_002516.2| Pseud...
      --id-regexp string                regular expression for parsing ID (default "^([^\\s]+)\\s?")
  -w, --line-width int                  line width when outputing FASTA format (0 for no wrap) (default 60)
  -o, --out-file string                 out file ("-" for stdout, suffix .gz for gzipped out) (default "-")
      --quiet                           be quiet and do not show extra information
  -t, --seq-type string                 sequence type (dna|rna|protein|unlimit|auto) (for auto, it automatically detect by the first sequence) (default "auto")
  -j, --threads int                     number of CPUs. (default value: 1 for single-CPU PC, 2 for others) (default 2)

Use "seqkit [command] --help" for more information about a command.

```

### Datasets

Datasets from [The miRBase Sequence Database -- Release 21](ftp://mirbase.org/pub/mirbase/21/)

- [`hairpin.fa.gz`](ftp://mirbase.org/pub/mirbase/21/hairpin.fa.gz)
- [`mature.fa.gz`](ftp://mirbase.org/pub/mirbase/21/mature.fa.gz)
- [`miRNA.diff.gz`](ftp://mirbase.org/pub/mirbase/21/miRNA.diff.gz)

Human genome from [ensembl](http://uswest.ensembl.org/info/data/ftp/index.html)
(For `seqkit subseq`)

- [`Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz`](ftp://ftp.ensembl.org/pub/release-84/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz)
- [`Homo_sapiens.GRCh38.84.gtf.gz`](ftp://ftp.ensembl.org/pub/release-84/gtf/homo_sapiens/Homo_sapiens.GRCh38.84.gtf.gz)
- `Homo_sapiens.GRCh38.84.bed.gz` is converted from `Homo_sapiens.GRCh38.84.gtf.gz`
by [`gtf2bed`](http://bedops.readthedocs.org/en/latest/content/reference/file-management/conversion/gtf2bed.html?highlight=gtf2bed)
with command

        zcat Homo_sapiens.GRCh38.84.gtf.gz | gtf2bed --do-not-sort | gzip -c > Homo_sapiens.GRCh38.84.bed.gz

Only DNA and gtf/bed data of Chr1 were used:

- `chr1.fa.gz`

            seqkit grep -p 1 Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz -o chr1.fa.gz

- `chr1.gtf.gz`

            zcat Homo_sapiens.GRCh38.84.gtf.gz | grep -w '^1' | gzip -c > chr1.gtf.gz

- `chr1.bed.gz`

            zcat Homo_sapiens.GRCh38.84.bed.gz | grep -w '^1' | gzip -c > chr1.bed.gz


## seq

Usage

```
transform sequences (revserse, complement, extract ID...)

Usage:
  seqkit seq [flags]

Flags:
  -p, --complement                complement sequence (blank for Protein sequence)
      --dna2rna                   DNA to RNA
  -G, --gap-letter string         gap letters (default "- ")
  -l, --lower-case                print sequences in lower case
  -n, --name                      only print names
  -i, --only-id                   print ID instead of full head
  -q, --qual                      only print qualities
  -g, --remove-gaps               remove gaps
  -r, --reverse                   reverse sequence)
      --rna2dna                   RNA to DNA
  -s, --seq                       only print sequences
  -u, --upper-case                print sequences in upper case
  -v, --validate-seq              validate bases according to the alphabet
  -V, --validate-seq-length int   length of sequence to validate (0 for whole seq) (default 10000)

```

Examples

1. Read and print

    - From file:

            $ seqkit seq hairpin.fa.gz
            >cel-let-7 MI0000001 Caenorhabditis elegans let-7 stem-loop
            UACACUGUGGAUCCGGUGAGGUAGUAGGUUGUAUAGUUUGGAAUAUUACCACCGGUGAAC
            UAUGCAAUUUUCUACCUUACCGGAGACAGAACUCUUCGA

            $ seqkit seq read_1.fq.gz
            @HWI-D00523:240:HF3WGBCXX:1:1101:2574:2226 1:N:0:CTGTAG
            TGAGGAATATTGGTCAATGGGCGCGAGCCTGAACCAGCCAAGTAGCGTGAAGGATGACTG
            CCCTACGGGTTGTAAACTTCTTTTATAAAGGAATAAAGTGAGGCACGTGTGCCTTTTTGT
            ATGTACTTTATGAATAAGGATCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGA
            TCCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGCAGGCGGT
            +
            HIHIIIIIHIIHGHHIHHIIIIIIIIIIIIIIIHHIIIIIHHIHIIIIIGIHIIIIHHHH
            HHGHIHIIIIIIIIIIIGHIIIIIGHIIIIHIIHIHHIIIIHIHHIIIIIIIGIIIIIII
            HIIIIIGHIIIIHIIIH?DGHEEGHIIIIIIIIIIIHIIHIIIHHIIHIHHIHCHHIIHG
            IHHHHHHH<GG?B@EHDE-BEHHHII5B@GHHF?CGEHHHDHIHIIH

    - From stdin:

            zcat hairpin.fa.gz | seqkit seq


1. Sequence types

    - By default, `seqkit seq` automatically detect the sequence type

            $ echo -e ">seq\nacgtryswkmbdhvACGTRYSWKMBDHV" | seqkit stat
            file  format  type  num_seqs  sum_len  min_len  avg_len  max_len
            -     FASTA   DNA          1       28       28       28       28

            $ echo -e ">seq\nACGUN ACGUN" | seqkit stat
            file  format  type  num_seqs  sum_len  min_len  avg_len  max_len
            -     FASTA   RNA          1       11       11       11       11

            $ echo -e ">seq\nabcdefghijklmnpqrstvwyz" | seqkit stat
            file  format  type     num_seqs  sum_len  min_len  avg_len  max_len
            -     FASTA   Protein         1       23       23       23       23

            $ echo -e "@read\nACTGCN\n+\n@IICCG" | seqkit stat
            file  format  type  num_seqs  sum_len  min_len  avg_len  max_len
            -     FASTQ   DNA          1        6        6        6        6

    - You can also set sequence type by flag `-t` (`--seq-type`).
      But this only take effect on subcommands `seq` and `locate`.

            $ echo -e ">seq\nabcdefghijklmnpqrstvwyz" | seqkit seq -t dna
            [INFO] when flag -t (--seq-type) given, flag -v (--validate-seq) is automatically switched on
            [ERRO] error when parsing seq: seq (invalid DNAredundant letter: e)


1. Only print names

    - Full name:

            $ seqkit seq hairpin.fa.gz -n
            cel-let-7 MI0000001 Caenorhabditis elegans let-7 stem-loop
            cel-lin-4 MI0000002 Caenorhabditis elegans lin-4 stem-loop
            cel-mir-1 MI0000003 Caenorhabditis elegans miR-1 stem-loop

    - Only ID:

            $ seqkit seq hairpin.fa.gz -n -i
            cel-let-7
            cel-lin-4
            cel-mir-1

    - Custom ID region by regular expression (this could be applied to all subcommands):

            $ seqkit seq hairpin.fa.gz -n -i --id-regexp "^[^\s]+\s([^\s]+)\s"
            MI0000001
            MI0000002
            MI0000003

1. Only print seq (global flag `-w` defines the output line width, 0 for no wrap)

        $ seqkit seq hairpin.fa.gz -s -w 0
        UACACUGUGGAUCCGGUGAGGUAGUAGGUUGUAUAGUUUGGAAUAUUACCACCGGUGAACUAUGCAAUUUUCUACCUUACCGGAGACAGAACUCUUCGA
        AUGCUUCCGGCCUGUUCCCUGAGACCUCAAGUGUGAGUGUACUAUUGAUGCUUCACACCUGGGCUCUCCGGGUACCAGGACGGUUUGAGCAGAU
        AAAGUGACCGUACCGAGCUGCAUACUUCCUUACAUGCCCAUACUAUAUCAUAAAUGGAUAUGGAAUGUAAAGAAGUAUGUAGAACGGGGUGGUAGU

1. Convert multi-line FASTQ to 4-line FASTQ

        $ seqkit seq reads_1.fq.gz -w 0

1. Reverse comlement sequence

        $ seqkit seq hairpin.fa.gz -r -p
        >cel-let-7 MI0000001 Caenorhabditis elegans let-7 stem-loop
        UCGAAGAGUUCUGUCUCCGGUAAGGUAGAAAAUUGCAUAGUUCACCGGUGGUAAUAUUCC
        AAACUAUACAACCUACUACCUCACCGGAUCCACAGUGUA

1. Remove gaps and to lower/upper case

        $ echo -e ">seq\nACGT-ACTGC-ACC" | seqkit seq -g -u
        >seq
        ACGTACTGCACC

1. RNA to DNA

        $ echo -e ">seq\nUCAUAUGCUUGUCUCAAAGAUUA" | seqkit seq --rna2dna
        >seq
        TCATATGCTTGTCTCAAAGATTA


## subseq

Usage

```
get subsequences by region/gtf/bed, including flanking sequences.

Recommendation: use plain FASTA file, so seqkit could utilize FASTA index.

The definition of region is 1-based and with some custom design.

Examples:

 1-based index    1 2 3 4 5 6 7 8 9 10
negative index    0-9-8-7-6-5-4-3-2-1
           seq    A C G T N a c g t n
           1:1    A
           2:4      C G T
         -4:-2                c g t
         -4:-1                c g t n
         -1:-1                      n
          2:-2      C G T N a c g t
          1:-1    A C G T N a c g t n

Usage:
  seqkit subseq [flags]

Flags:
      --bed string        by BED file
      --chr value         select limited sequence with sequence IDs (multiple value supported, case ignored) (default [])
  -d, --down-stream int   down stream length
      --feature value     select limited feature types (multiple value supported, case ignored, only works with GTF) (default [])
      --gtf string        by GTF (version 2.2) file
  -f, --only-flank        only return up/down stream sequence
  -r, --region string     by region. e.g 1:12 for first 12 bases, -12:-1 for last 12 bases, 13:-1 for cutting first 12 bases. type "seqkit subseq -h" for more examples
  -u, --up-stream int     up stream length

```

Examples

***Recommendation: use plain FASTA file, so seqkit could utilize FASTA index.***

1. First 12 bases

        $ zcat hairpin.fa.gz | seqkit subseq -r 1:12

1. Last 12 bases

        $ zcat hairpin.fa.gz | seqkit subseq -r -12:-1

1. Subsequences without first and last 12 bases

        $ zcat hairpin.fa.gz | seqkit subseq -r 13:-13

1. Get subsequence by GTF file

        $ cat t.fa
        >seq
        actgACTGactgn
        $ cat t.gtf
        seq     test    CDS     5       8       .       .       .       gene_id "A"; transcript_id "";
        seq     test    CDS     5       8       .       -       .       gene_id "B"; transcript_id "";
        $ seqkit

        $ seqkit subseq --gtf t.gtf t.fa
        >seq_5:8:. A
        ACTG
        >seq_5:8:- B
        CAGT

    Human genome example:

    ***AVOID loading all data from Homo_sapiens.GRCh38.84.gtf.gz,
    the uncompressed data are so big and may exhaust your RAM.***

    We could specify chromesomes and features.

        $ seqkit subseq --gtf Homo_sapiens.GRCh38.84.gtf.gz --chr 1 --feature cds  hsa.fa > chr1.gtf.cds.fa

        $ seqkit stat chr1.gtf.cds.fa
        file             format  type  num_seqs    sum_len  min_len  avg_len  max_len
        chr1.gtf.cds.fa  FASTA   DNA     65,012  9,842,274        1    151.4   12,045

1. Get CDS and 3bp up-stream sequences

        $ seqkit subseq --gtf t.gtf t.fa -u 3
        >seq_5:8:._us:3 A
        ctgACTG
        >seq_5:8:-_us:3 B
        agtCAGT

1. Get 3bp up-stream sequences of CDS, not including CDS

        $ seqkit subseq --gtf t.gtf t.fa -u 3 -f
        >seq_5:8:._usf:3 A
        ctg
        >seq_5:8:-_usf:3 B
        agt

1. Get subsequences by BED file.

    ***AVOID loading all data from Homo_sapiens.GRCh38.84.gtf.gz,
    the uncompressed data are so big and may exhaust your RAM.***

        $  seqkit subseq --bed Homo_sapiens.GRCh38.84.bed.gz --chr 1 hsa.fa >  chr1.bed.gz.fa

    We may need to remove duplicated sequences

        $ seqkit subseq --bed Homo_sapiens.GRCh38.84.bed.gz --chr 1 hsa.fa | seqkit rmdup > chr1.bed.rmdup.fa
        [INFO] 141060 duplicated records removed

    Summary:

        $ seqkit stat chr1.gz.*.gz
        file               seq_format   seq_type   num_seqs   min_len   avg_len     max_len
        chr1.gz.fa         FASTA        DNA         231,974         1   3,089.5   1,551,957
        chr1.gz.rmdup.fa   FASTA        DNA          90,914         1   6,455.8   1,551,957


## sliding

Usage

```
sliding sequences, circular genome supported

Usage:
  seqkit sliding [flags]

Flags:
  -C, --circular-genome   circular genome
  -s, --step int        step size
  -W, --window int      window size

```

Examples

1. General use

        $ echo -e ">seq\nACGTacgtNN" | seqkit sliding -s 3 -W 6
        >seq_sliding:1-6
        ACGTac
        >seq_sliding:4-9
        TacgtN

2. Circular genome

        $ echo -e ">seq\nACGTacgtNN" | seqkit sliding -s 3 -W 6 -C
        >seq_sliding:1-6
        ACGTac
        >seq_sliding:4-9
        TacgtN
        >seq_sliding:7-2
        gtNNAC
        >seq_sliding:10-5
        NACGTa

3. Generate GC content for ploting

        $ zcat hairpin.fa.gz | seqkit fx2tab | head -n 1 | seqkit tab2fx | seqkit sliding -s 5 -W 30 | seqkit fx2tab -n -g
        cel-let-7_sliding:1-30          50.00
        cel-let-7_sliding:6-35          46.67
        cel-let-7_sliding:11-40         43.33
        cel-let-7_sliding:16-45         36.67
        cel-let-7_sliding:21-50         33.33
        cel-let-7_sliding:26-55         40.00
        ...

## stat

Usage

```
simple statistics of FASTA files

Usage:
  seqkit stat [flags]

```

Eexamples

1. General use

        $ seqkit stat *.f{a,q}.gz
        file           format  type  num_seqs    sum_len  min_len  avg_len  max_len
        hairpin.fa.gz  FASTA   RNA     28,645  2,949,871       39      103    2,354
        mature.fa.gz   FASTA   RNA     35,828    781,222       15     21.8       34
        reads_1.fq.gz  FASTQ   DNA      2,500    567,516      226      227      229
        reads_2.fq.gz  FASTQ   DNA      2,500    560,002      223      224      225

## fq2fa

Usage

```
covert FASTQ to FASTA

Usage:
  seqkit fq2fa [flags]

```

Examples

    seqkit fq2fa reads_1.fq.gz -o reads1_.fa.gz


## fx2tab & tab2fx

Usage (fx2tab)

```
covert FASTA/Q to tabular format, and provide various information,
like sequence length, GC content/GC skew.

Usage:
  seqkit fx2tab [flags]

Flags:
  -B, --base-content value   print base content. (case ignored, multiple values supported) e.g. -B AT -B N (default [])
  -g, --gc                   print GC content
  -G, --gc-skew              print GC-Skew
  -H, --header-line          print header line
  -l, --length               print sequence length
  -n, --name                 only print names (no sequences and qualities)
  -i, --only-id              print ID instead of full head

```

Usage (tab2fx)

```
covert tabular format (first two/three columns) to FASTA/Q format

Usage:
  seqkit tab2fx [flags]

Flags:
  -p, --comment-line-prefix value   comment line prefix (default [#,//])


```

Examples

1. Default output

        $ seqkit fx2tab hairpin.fa.gz | head -n 2
        cel-let-7 MI0000001 Caenorhabditis elegans let-7 stem-loop      UACACUGUGGAUCCGGUGAGGUAGUAGGUUGUAUAGUUUGGAAUAUUACCACCGGUGAACUAUGCAAUUUUCUACCUUACCGGAGACAGAACUCUUCGA
        cel-lin-4 MI0000002 Caenorhabditis elegans lin-4 stem-loop      AUGCUUCCGGCCUGUUCCCUGAGACCUCAAGUGUGAGUGUACUAUUGAUGCUUCACACCUGGGCUCUCCGGGUACCAGGACGGUUUGAGCAGAU


1. Print sequence length, GC content, and only print names (no sequences),
we could also print title line by flag `-T`.

        $ seqkit fx2tab hairpin.fa.gz -l -g -n -i -H | head -n 4 | csvtk -t -C '&' pretty
        #name       seq   qual   length   GC
        cel-let-7                99       43.43
        cel-lin-4                94       54.26
        cel-mir-1                96       40.62

1. Use fx2tab and tab2fx in pipe

        $ zcat hairpin.fa.gz | seqkit fx2tab | seqkit tab2fx

        $ zcat reads_1.fq.gz | seqkit fx2tab | seqkit tab2fx

1. Sort sequences by length (use `seqkit sort -l`)

        $ zcat hairpin.fa.gz | seqkit fx2tab -l | sort -t"`echo -e '\t'`" -n -k4,4 | seqkit tab2fx
        >cin-mir-4129 MI0015684 Ciona intestinalis miR-4129 stem-loop
        UUCGUUAUUGGAAGACCUUAGUCCGUUAAUAAAGGCAUC
        >mmu-mir-7228 MI0023723 Mus musculus miR-7228 stem-loop
        UGGCGACCUGAACAGAUGUCGCAGUGUUCGGUCUCCAGU
        >cin-mir-4103 MI0015657 Ciona intestinalis miR-4103 stem-loop
        ACCACGGGUCUGUGACGUAGCAGCGCUGCGGGUCCGCUGU

        $ seqkit sort -l hairpin.fa.gz

    Sorting or filtering by GC (or other base by -flag `-B`) content could also achieved in similar way.

1. Get first 1000 sequences

        $ seqkit fx2tab hairpin.fa.gz | head -n 1000 | seqkit tab2fx

        $ seqkit fx2tab reads_1.fq.gz | head -n 1000 | seqkit tab2fx

**Extension**

After converting FASTA to tabular format with `seqkit fx2tab`,
it could be handled with CSV/TSV tools,
 e.g. [csvtk](https://github.com/shenwei356/csvtkt), a cross-platform, efficient and practical CSV/TSV toolkit

- `csvtk grep` could be used to filter sequences (similar with `seqkit grep`)
- `csvtk inter`
computates intersection of multiple files. It could achieve similar function
as `seqkit common -n` along with shell.
- `csvtk join` joins multiple CSV/TSV files by multiple IDs.
- [csv_melt](https://github.com/shenwei356/datakit/blob/master/csv_melt)
provides melt function, could be used in preparation of data for ploting.


## grep

Usage

```
search sequences by pattern(s) of name or sequence motifs

Usage:
  seqkit grep [flags]

Flags:
  -n, --by-name               match by full name instead of just id
  -s, --by-seq                match by seq
  -d, --degenerate            pattern/motif contains degenerate base
      --delete-matched        delete matched pattern to speedup
  -i, --ignore-case           ignore case
  -v, --invert-match          invert the sense of matching, to select non-matching records
  -p, --pattern value         search pattern (multiple values supported) (default [])
  -f, --pattern-file string   pattern file
  -r, --use-regexp            patterns are regular expression

```

Examples

1. Extract human hairpins (i.e. sequences with name starting with `hsa`)

        $ zcat hairpin.fa.gz | seqkit grep -r -p ^hsa
        >hsa-let-7a-1 MI0000060 Homo sapiens let-7a-1 stem-loop
        UGGGAUGAGGUAGUAGGUUGUAUAGUUUUAGGGUCACACCCACCACUGGGAGAUAACUAU
        ACAAUCUACUGUCUUUCCUA
        >hsa-let-7a-2 MI0000061 Homo sapiens let-7a-2 stem-loop
        AGGUUGAGGUAGUAGGUUGUAUAGUUUAGAAUUACAUCAAGGGAGAUAACUGUACAGCCU
        CCUAGCUUUCCU

1. Remove human and mice hairpins.

        $ zcat hairpin.fa.gz | seqkit grep -r -p ^hsa -p ^mmu -v

1. Extract new entries by information from miRNA.diff.gz

    1. Get IDs of new entries.

            $ zcat miRNA.diff.gz | grep ^# -v | grep NEW | cut -f 2 > list
            $ more list
            cfa-mir-486
            cfa-mir-339-1
            pmi-let-7


    2. Extract by ID list file

            $ zcat hairpin.fa.gz | seqkit grep -f list > new.fa

1. Extract sequences starting with AGGCG

        $ zcat hairpin.fa.gz | seqkit grep -s -r -i -p ^aggcg

1. Extract sequences with TTSAA (AgsI digest site) in SEQUENCE. Base S stands for C or G.

        $ zcat hairpin.fa.gz | seqkit grep -s -d -i -p TTSAA

    It's equal to but simpler than:

        $ zcat hairpin.fa.gz | seqkit grep -s -r -i -p TT[CG]AA


## locate

Usage

```
locate subsequences/motifs

Motifs could be EITHER plain sequence containing "ACTGN" OR regular
expression like "A[TU]G(?:.{3})+?[TU](?:AG|AA|GA)" for ORFs.
Degenerate bases like "RYMM.." are also supported by flag -d.

By default, motifs are treated as regular expression.
When flag -d given, regular expression may be wrong.
For example: "\w" will be wrongly converted to "\[AT]".

Usage:
  seqkit locate [flags]

Flags:
  -d, --degenerate                pattern/motif contains degenerate base
  -i, --ignore-case               ignore case
  -P, --only-positive-strand      only search at positive strand
  -p, --pattern value             search pattern/motif (multiple values supported) (default [])
  -f, --pattern-file string       pattern/motif file (FASTA format)
  -V, --validate-seq-length int   length of sequence to validate (0 for whole seq) (default 10000)

```

Examples

1. Locate ORFs.

        $ zcat hairpin.fa.gz | seqkit locate -i -p "A[TU]G(?:.{3})+?[TU](?:AG|AA|GA)"
        seqID   patternName     pattern strand  start   end     matched
        cel-lin-4       A[TU]G(?:.{3})+?[TU](?:AG|AA|GA)        A[TU]G(?:.{3})+?[TU](?:AG|AA|GA)        +  136      AUGCUUCCGGCCUGUUCCCUGAGACCUCAAGUGUGA
        cel-mir-1       A[TU]G(?:.{3})+?[TU](?:AG|AA|GA)        A[TU]G(?:.{3})+?[TU](?:AG|AA|GA)        +  54       95      AUGGAUAUGGAAUGUAAAGAAGUAUGUAGAACGGGGUGGUAG
        cel-mir-1       A[TU]G(?:.{3})+?[TU](?:AG|AA|GA)        A[TU]G(?:.{3})+?[TU](?:AG|AA|GA)        -  43       51      AUGAUAUAG

1. Locate Motif.

        $ zcat hairpin.fa.gz | seqkit locate -i -d -p AUGGACUN
        seqID         patternName   pattern    strand   start   end   matched
        cel-mir-58a   AUGGACUN      AUGGACUN   +        81      88    AUGGACUG
        ath-MIR163    AUGGACUN      AUGGACUN   -        122     129   AUGGACUC

    Notice that `seqkit grep` only searches in positive strand, but `seqkit loate` could recognize both strand


## rmdup

Usage

```
remove duplicated sequences by id/name/sequence

Usage:
  seqkit rmdup [flags]

Flags:
    -n, --by-name                by full name instead of just id
    -s, --by-seq                 by seq
    -D, --dup-num-file string    file to save number and list of duplicated seqs
    -d, --dup-seqs-file string   file to save duplicated seqs
    -i, --ignore-case            ignore case
    -m, --md5                    use MD5 instead of original seqs to reduce memory usage when comparing by seqs

```

Examples

Similar to `common`.

1. General use

        $ zcat hairpin.fa.gz | seqkit rmdup -s -o clean.fa.gz
        [INFO] 2226 duplicated records removed

        $ zcat reads_1.fq.gz | seqkit rmdup -s -o clean.fa.gz
        [INFO] 1086 duplicated records removed

1. Save duplicated sequences to file

        $ zcat hairpin.fa.gz | seqkit rmdup -s -i -m -o clean.fa.gz -d duplicated.fa.gz -D duplicated.detail.txt

        $ cat duplicated.detail.txt   # here is not the entire list
        3       hsa-mir-424, mml-mir-424, ppy-mir-424
        3       hsa-mir-342, mml-mir-342, ppy-mir-342
        2       ngi-mir-932, nlo-mir-932
        2       ssc-mir-9784-1, ssc-mir-9784-2

## common

Usage

```
find common sequences of multiple files by id/name/sequence

Usage:
  seqkit common [flags]

Flags:
    -n, --by-name       match by full name instead of just id
    -s, --by-seq        match by sequence
    -i, --ignore-case   ignore case
    -m, --md5           use MD5 instead of original seqs to reduce memory usage when comparing by seqs

```

Examples

1. By ID (default)

        seqkit common file*.fa -o common.fasta

1. By full name

        seqkit common file*.fa -n -o common.fasta

1. By sequence

        seqkit common file*.fa -s -i -o common.fasta

1. By sequence (***for large sequences***)

        seqkit common file*.fa -s -i -o common.fasta --md5


## split

Usage

```
split sequences into files by name ID, subsequence of given region,
part size or number of parts.

The definition of region is 1-based and with some custom design.

Examples:

 1-based index    1 2 3 4 5 6 7 8 9 10
negative index    0-9-8-7-6-5-4-3-2-1
           seq    A C G T N a c g t n
           1:1    A
           2:4      C G T
         -4:-2                c g t
         -4:-1                c g t n
         -1:-1                      n
          2:-2      C G T N a c g t
          1:-1    A C G T N a c g t n

Usage:
  seqkit split [flags]

Flags:
Flags:
  -i, --by-id              split squences according to sequence ID
  -p, --by-part int        split squences into N parts
  -r, --by-region string   split squences according to subsequence of given region. e.g 1:12 for first 12 bases, -12:-1 for last 12 bases. type "seqkit split -h" for more examples
  -s, --by-size int        split squences into multi parts with N sequences
  -d, --dry-run            dry run, just print message and no files will be created.
  -f, --force              overwrite output directory
  -k, --keep-temp          keep tempory FASTA and .fai file when using 2-pass mode
  -m, --md5                use MD5 instead of region sequence in output file when using flag -r (--by-region)
  -O, --out-dir string     output directory (default value is infile.split)
  -2, --two-pass           two-pass mode read files twice to lower memory usage. (only for FASTA format)

```

Examples

1. Split sequences into parts with at most 10000 sequences

        $ seqkit split hairpin.fa.gz -s 10000
        [INFO] split into 10000 seqs per file
        [INFO] write 10000 sequences to file: hairpin.fa.part_001.gz
        [INFO] write 10000 sequences to file: hairpin.fa.part_002.gz
        [INFO] write 8645 sequences to file: hairpin.fa.part_003.gz

1. Split sequences into 4 parts

        $ seqkit split hairpin.fa.gz -p 4
        [INFO] split into 4 parts
        [INFO] read sequences ...
        [INFO] read 28645 sequences
        [INFO] write 7162 sequences to file: hairpin.fa.part_001.gz
        [INFO] write 7162 sequences to file: hairpin.fa.part_002.gz
        [INFO] write 7162 sequences to file: hairpin.fa.part_003.gz
        [INFO] write 7159 sequences to file: hairpin.fa.part_004.gz


    ***To reduce memory usage when spliting big file, we should alwasy use flag `--two-pass`***

        $ seqkit split hairpin.fa.gz -p 4 -2
        [INFO] split into 4 parts
        [INFO] read and write sequences to tempory file: hairpin.fa.gz.fa ...
        [INFO] create and read FASTA index ...
        [INFO] read sequence IDs from FASTA index ...
        [INFO] 28645 sequences loaded
        [INFO] write 7162 sequences to file: hairpin.part_001.fa.gz
        [INFO] write 7162 sequences to file: hairpin.part_002.fa.gz
        [INFO] write 7162 sequences to file: hairpin.part_003.fa.gz
        [INFO] write 7159 sequences to file: hairpin.part_004.fa.gz

1. Split sequences by species. i.e. by custom IDs (first three letters)

        $ seqkit split hairpin.fa.gz -i --id-regexp "^([\w]+)\-" -2
        [INFO] split by ID. idRegexp: ^([\w]+)\-
        [INFO] read and write sequences to tempory file: hairpin.fa.gz.fa ...
        [INFO] create and read FASTA index ...
        [INFO] create FASTA index for hairpin.fa.gz.fa
        [INFO] read sequence IDs from FASTA index ...
        [INFO] 28645 sequences loaded
        [INFO] write 48 sequences to file: hairpin.id_cca.fa.gz
        [INFO] write 3 sequences to file: hairpin.id_hci.fa.gz
        [INFO] write 106 sequences to file: hairpin.id_str.fa.gz
        [INFO] write 1 sequences to file: hairpin.id_bkv.fa.gz
        ...

1. Split sequences by sequence region (for example, sequence barcode)

        $ seqkit split hairpin.fa.gz -r 1:3 -2
        [INFO] split by region: 1:3
        [INFO] read and write sequences to tempory file: hairpin.fa.gz.fa ...
        [INFO] read sequence IDs and sequence region from FASTA file ...
        [INFO] create and read FASTA index ...
        [INFO] write 463 sequences to file: hairpin.region_1:3_AUG.fa.gz
        [INFO] write 349 sequences to file: hairpin.region_1:3_ACU.fa.gz
        [INFO] write 311 sequences to file: hairpin.region_1:3_CGG.fa.gz

    **If region is too long, we could use falg `--md5`**,
    i.e. use MD5 instead of region sequence in output file.

    Sequence suffix could be defined as `-r -12:-1`

## sample

Usage

```
sample sequences by number or proportion.

Usage:
  seqkit sample [flags]

Flags:
  -n, --number int         sample by number (result may not exactly match)
  -p, --proportion float   sample by proportion
  -s, --rand-seed int      rand seed for shuffle (default 11)
  -2, --two-pass           2-pass mode read files twice to lower memory usage. Not allowed when reading from stdin

```

Examples

1. Sample by proportion

        $ zcat hairpin.fa.gz | seqkit sample -p 0.1 -o sample.fa.gz
        [INFO] sample by proportion
        [INFO] 2814 sequences outputed

1. Sample by number

        $ zcat hairpin.fa.gz | seqkit sample -n 1000 -o sample.fa.gz
        [INFO] sample by number
        [INFO] 949 sequences outputed

    ***To reduce memory usage when spliting big file, we could use flag `--two-pass`***

    ***We can also use `seqkit sample -p` followed with `seqkit head -n`:***

        $ zcat hairpin.fa.gz | seqkit sample -p 0.1 | seqkit head -n 1000 -o sample.fa.gz

1. Set rand seed to reproduce the result

        $ zcat hairpin.fa.gz | seqkit sample -p 0.1 -s 11

1. Most of the time, we could shuffle after sampling

        $ zcat hairpin.fa.gz | seqkit sample -p 0.1 | seqkit shuffle -o sample.fa.gz

Note that when sampling on FASTQ files, make sure using same random seed by
flag `-s` (`--rand-seed`)

## head

Usage

```
print first N FASTA/Q records

Usage:
  seqkit head [flags]

Flags:
  -n, --number int   print first N FASTA/Q records (default 10)

```

Examples

1. FASTA

        $ seqkit head -n 1 hairpin.fa.gz
        >cel-let-7 MI0000001 Caenorhabditis elegans let-7 stem-loop
        UACACUGUGGAUCCGGUGAGGUAGUAGGUUGUAUAGUUUGGAAUAUUACCACCGGUGAAC
        UAUGCAAUUUUCUACCUUACCGGAGACAGAACUCUUCGA

1. FASTQ

        $ seqkit head -n 1 reads_1.fq.gz
        @HWI-D00523:240:HF3WGBCXX:1:1101:2574:2226 1:N:0:CTGTAG
        TGAGGAATATTGGTCAATGGGCGCGAGCCTGAACCAGCCAAGTAGCGTGAAGGATGACTG
        CCCTACGGGTTGTAAACTTCTTTTATAAAGGAATAAAGTGAGGCACGTGTGCCTTTTTGT
        ATGTACTTTATGAATAAGGATCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGA
        TCCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGCAGGCGGT
        +
        HIHIIIIIHIIHGHHIHHIIIIIIIIIIIIIIIHHIIIIIHHIHIIIIIGIHIIIIHHHH
        HHGHIHIIIIIIIIIIIGHIIIIIGHIIIIHIIHIHHIIIIHIHHIIIIIIIGIIIIIII
        HIIIIIGHIIIIHIIIH?DGHEEGHIIIIIIIIIIIHIIHIIIHHIIHIHHIHCHHIIHG
        IHHHHHHH<GG?B@EHDE-BEHHHII5B@GHHF?CGEHHHDHIHIIH



## replace

Usage

```
replace name/sequence/by regular expression.

Note that the replacement supports capture variables.
e.g. $1 represents the text of the first submatch.
ATTENTION: use SINGLE quote NOT double quotes in *nix OS.

Examples: Adding space to all bases.

    seqkit replace -p "(.)" -r '$1 ' -s

Or use the \ escape character.

    seqkit replace -p "(.)" -r "\$1 " -s

more on: http://shenwei356.github.io/seqkit/usage/#replace

Usage:
  seqkit replace [flags]

Flags:
  -s, --by-seq               replace seq
  -i, --ignore-case          ignore case
  -p, --pattern string       search regular expression
  -r, --replacement string   replacement. supporting capture variables.  e.g. $1 represents the text of the first submatch. ATTENTION: use SINGLE quote NOT double quotes in *nix OS or use the \ escape character. record number is also supported by "{NR}"

```

Examples

1. Remove descriptions

        $ echo -e ">seq1 abc-123\nACGT-ACGT" | seqkit replace -p " .+"
        >seq1
        ACGT-ACGT

1. Replace "-" with "="

        $ echo -e ">seq1 abc-123\nACGT-ACGT" | seqkit replace -p "\-" -r '='
        >seq1 abc=123
        ACGT-ACGT

1. Remove gaps in sequences.

        $ echo -e ">seq1 abc-123\nACGT-ACGT" | seqkit replace -p " |-" -s
        >seq1 abc-123
        ACGTACGT

1. Add space to every base. **ATTENTION: use SINGLE quote NOT double quotes in *nix OS**

        $ echo -e ">seq1 abc-123\nACGT-ACGT" | seqkit replace -p "(.)" -r '$1 ' -s
        >seq1 abc-123
        A C G T - A C G T

1. Transpose sequence with [csvtk](https://github.com/shenwei356/csvtk)

        $ echo -e ">seq1\nACTGACGT\n>seq2\nactgccgt" | seqkit replace -p "(.)" -r     "\$1 " -s | seqkit seq -s -u | csvtk space2tab | csvtk -t transpose
        A       A
        C       C
        T       T
        G       G
        A       C
        C       C
        G       G
        T       T

1. Rename with number of record

        echo -e ">abc\nACTG\n>123\nATTT" |  seqkit replace -p .+ -r "seq_{NR}"
        >seq_1
        ACTG
        >seq_2
        ATTT


## shuffle

Usage

```
shuffle sequences.

By default, all records will be readed into memory.
For FASTA format, use flag -2 (--two-pass) to reduce memory usage. FASTQ not
supported.

Firstly, seqkit reads the sequence IDs. If the file is not plain FASTA file,
seqkit will write the sequences to tempory files, and create FASTA index.

Secondly, seqkit shuffles sequence IDs and extract sequences by FASTA index.

Usage:
  seqkit shuffle [flags]

Flags:
  -k, --keep-temp       keep tempory FASTA and .fai file when using 2-pass mode
  -s, --rand-seed int   rand seed for shuffle (default 23)
  -2, --two-pass        two-pass mode read files twice to lower memory usage. (only for FASTA format)

```

Examples

1. General use.

        $ seqkit shuffle hairpin.fa.gz > shuffled.fa
        [INFO] read sequences ...
        [INFO] 28645 sequences loaded
        [INFO] shuffle ...
        [INFO] output ...

1. ***For big genome, you'd better use two-pass mode*** so seqkit could use
   FASTA index to reduce memory usage

        $ time seqkit shuffle -2 hsa.fa > shuffle.fa
        [INFO] create and read FASTA index ...
        [INFO] create FASTA index for hsa.fa
        [INFO] read sequence IDs from FASTA index ...
        [INFO] 194 sequences loaded
        [INFO] shuffle ...
        [INFO] output ...

        real    0m35.080s
        user    0m45.521s
        sys     0m3.411s

Note that when sampling on FASTQ files, make sure using same random seed by
flag `-s` (`--rand-seed`) for read 1 and 2 files.

## sort

Usage

```
sort sequences by id/name/sequence/length.

By default, all records will be readed into memory.
For FASTA format, use flag -2 (--two-pass) to reduce memory usage. FASTQ not
supported.

Firstly, seqkit reads the sequence head and length information.
If the file is not plain FASTA file,
seqkit will write the sequences to tempory files, and create FASTA index.

Secondly, seqkit sort sequence by head and length information
and extract sequences by FASTA index.

Usage:
  seqkit sort [flags]

Flags:
  -l, --by-length               by sequence length
  -n, --by-name                 by full name instead of just id
  -s, --by-seq                  by sequence
  -i, --ignore-case             ignore case
  -k, --keep-temp               keep tempory FASTA and .fai file when using 2-pass mode
  -r, --reverse                 reverse the result
  -L, --seq-prefix-length int   length of sequence prefix on which seqkit sorts by sequences (0 for whole sequence) (default 10000)
  -2, --two-pass                two-pass mode read files twice to lower memory usage. (only for FASTA format)

```

Examples

***For FASTA format, use flag -2 (--two-pass) to reduce memory usage***

1. sort by ID

        $ echo -e ">seq1\nACGTNcccc\n>SEQ2\nacgtnAAAA" | seqkit sort --quiet
        >SEQ2
        acgtnAAAA
        >seq1
        ACGTNcccc

1. sort by ID, ignoring case.

        $ echo -e ">seq1\nACGTNcccc\n>SEQ2\nacgtnAAAA" | seqkit sort --quiet -i
        >seq1
        ACGTNcccc
        >SEQ2
        acgtnAAAA

1. sort by seq, ignoring case.

        $ echo -e ">seq1\nACGTNcccc\n>SEQ2\nacgtnAAAA" | seqkit sort --quiet -s -i
        >SEQ2
        acgtnAAAA
        >seq1
        ACGTNcccc

1. sort by sequence length

        $ echo -e ">seq1\nACGTNcccc\n>SEQ2\nacgtnAAAAnnn\n>seq3\nacgt" | seqkit sort --quiet -l
        >seq3
        acgt
        >seq1
        ACGTNcccc
        >SEQ2
        acgtnAAAAnnn

<div id="disqus_thread"></div>
<script>
/**
* RECOMMENDED CONFIGURATION VARIABLES: EDIT AND UNCOMMENT THE SECTION BELOW TO INSERT DYNAMIC VALUES FROM YOUR PLATFORM OR CMS.
* LEARN WHY DEFINING THESE VARIABLES IS IMPORTANT: https://disqus.com/admin/universalcode/#configuration-variables
*/
/*
var disqus_config = function () {
this.page.url = PAGE_URL; // Replace PAGE_URL with your page's canonical URL variable
this.page.identifier = PAGE_IDENTIFIER; // Replace PAGE_IDENTIFIER with your page's unique identifier variable
};
*/
(function() { // DON'T EDIT BELOW THIS LINE
var d = document, s = d.createElement('script');

s.src = '//fastakit.disqus.com/embed.js';

s.setAttribute('data-timestamp', +new Date());
(d.head || d.body).appendChild(s);
})();
</script>
<noscript>Please enable JavaScript to view the <a href="https://disqus.com/?ref_noscript" rel="nofollow">comments powered by Disqus.</a></noscript>
