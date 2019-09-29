# ASCII-style alignment pileup

## Description

Generates an ASCII-style pileup of read alignments in one or more BAM files
against one or more regions specified in a BED file.

## Usage

```sh
ascii_alignment_pileup.R [-hv] [OPTIONS] bed bam [bam2 ...]
```

## Requirements

* R 3.6.0
* GenomicAlignments 1.20.0
* optparse 1.6.2
* rtracklayer 1.44.0
* tools 3.6.0

## Input files

* [BED](https://www.ensembl.org/info/website/upload/bed.html) file
* [BAM](https://samtools.github.io/hts-specs/) file(s)
* Optional: [FASTA](https://en.wikipedia.org/wiki/FASTA_format) file compressed with [`bgzip`](http://www.htslib.org/doc/bgzip.html)
* Optional: [GFF/GTF/GFF3](https://en.wikipedia.org/wiki/General_feature_format) file

## Output files

* Custom file format:

```console
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>        hsa-let-7a-1
.....>>>>>>>>>>>>>>>>>>>>>>.....................................................        hsa-let-7a-5p
........................................................>>>>>>>>>>>>>>>>>>>>>...        hsa-let-7a-3p
TGGGATGAGGTAGTAGGTTGTATAGTTTTAGGGTCACACCCACCACTGGGAGATAACTATACAATCTACTGTCTTTCCTA        9:94175957-94176036:+
TACCATGAGGTAGTAGGTTGTATAGTT.....................................................        1
...CATGAGGTAGTAGGTTGTATAGTT.....................................................        10
...GA-GAGGTAGTAGGTTGTATAGTT.....................................................        2
....A-GAGGTAGTAGGTTGTATAGTT.....................................................        19
.....TGAGGTAGTAGGTTGTATAGTT.....................................................        17
......GAGGTAGTAGGTTGTATAGTT.....................................................        33
.......AGGTAGTAGGTTGTATAGTT.....................................................        9
.......AGGTAGTAGGTTGTATAGTTT....................................................        2
........GGTAGTAGGTTGTATAGTT.....................................................        7
...................................................GATAACTATACAATCTACTGTCTT.....        1
......................................................AACTATACAATCTACT..........        1
........................................................CTATACAATCTACTGTCTTTCT..        28
........................................................CTATACAATCTACTGTCTTTC-T.        22
........................................................CTATACAATCTACTGTCTTTCC..        19
........................................................CTATACAATCTACTGTCTTTC...        12
........................................................CTATACAATCTACTGTCTTTCTT.        2
........................................................CTATACAATCTACTGTC.......        1
........................................................CTATACAATCTACTGTCTT.....        1
........................................................CTATACAATCTACTGTCTTTCG..        1
.........................................................TATACAATCTACTGTCTTTCT..        4
.........................................................TATACAATCTACTGTCTTTC-T.        4
.........................................................TATACAATCTACTGTCTTTC...        2
.........................................................TATACAATCTACTGTCTTTCC..        1
.........................................................TATACAATCTACTGTCTTTCCT.        1
```

## Options

```console
--reference=FILE
        Reference genome sequence in FASTA format. The file *MUST* be compressed
    with BGZIP. If supplied, the reference sequence for the query region(s) will
    be added to the output. Note that on the first run with a specific reference
    genome file, an FAI index is generated which will take some time.

--annotations=FILE
        Annotation file in GFF/GTF format used to annotate sequences. If
    supplied, features overlapping the query region(s) will be visualized in the
    output. Ensure that the argument to option `annotation-name-field`
    corresponds to a field in the annotations, otherwise the script will fail.

--output-directory=DIR
        Output directory. One output file will be created for each region in
    `--bed` and the filenames will be generated from the basenames of the
    supplied BAM file(s) and the name field (4th column) of the BED file.
    [default "."]

--maximum-region-width=INT
        Maximum input region width. Use with care as wide regions will use
    excessive resources. [default 200]

--do-not-collapse-alignments
        Show alignments of reads with identical sequences individually.

--minimum-count=INT
        Alignments of reads with less copies than the specified number will not
    be printed. Option is not considered if `do-not-collapse-alignments` is
    set. [default 1]

--annotation-name-field=STR
        Annotation field used to populate the `name` column in the output.
    [default "Name"]

--padding-character=CHAR
        Character used for padding alignments. [default "."]

--indel-character=CHAR
        Character to denote insertions and deletions in alignments.
    [default "-"]

-h, --help
        Show this information and die.

-v, --verbose
        Print log messages to STDOUT.
```

## Tags

* BED
* BAM
* GFF
* FASTA
* alignments
* pileup
* ASCII
* microRNA
* miRNA
* miR
* isomiR

## Version

29-SEP-2019: 1.0.0

## Contact

* Author: Alexander Kanitz
* Affiliation: Biozentrum, University of Basel
* Email: alexander.kanitz@alumni.ethz.ch

