# Mapping and post processing

In this section we will cover mapping and post-processing of alignment files
Bellow are links to sections we will cover



## 3.1 Mapping aka read alignment

The next step is to align the reads to genome assembly.
This is done using Bowtie2 tool. 
The resulting .sam files are next transformed to .bam files and filtered for best aligned reads using samtools. 
PCR duplicates are removed.
The resulting BAM files and your primary output and can be then used for the next steps.

> **NOTE**: Note that bowtie2 does not support compressed files thus we need to decompress.

```bash
gunzip -k -d GFP_Mouse_ES_merged.fastq.gz

gunzip -k -d Nanog_Mouse_ES_merged.fastq.gz
```



## Alignment

The NGS reads are aligned with Alignment tool against the reference genome sequence. 

In ChIP-Seq experiments, it is usually more appropriate to eliminate reads mapping to multiple locations and to perform soft-clipping to eliminate any low quality sequence from the reads.

![](https://alexpmagalhaes.github.io/ChIPseq_course/img/Alignment_errors.png)


In theory, this sounds like a very simple case of string matching. We take the sequence read and figure out where it originated from in the reference genome. However, in practice, this is **actually quite difficult!** This is because:

* The reference genome we are searching is large and complex (e.g. the human genome is ~3,200,000,000bp). 
* By contrast, the reads we are searching for are much smaller (50-150bp), and they are on the range of millions for a given sample.
* We have to consider non-exact matching of the read to the reference due to natural variation and sequencing errors.
* We have to consider non-unique alignment due to the short length of reads and high percentage of repetitive regions in the genome (e.g. repetitive regions = >50% of the human genome).

There are many different tools that have been developed for alignment of next-generation sequencing data, and some that are more suitable to different technologies. A popular tool commonly used with ChIP-seq data, and the one that we will be using in this workshop is [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml).
Bowtie2 is a fast and accurate alignment tool that supports gaped, local and paired-end alignment modes and works best for reads that are **at least 50 bp** (shorter read lengths should use [Bowtie1](http://bowtie-bio.sourceforge.net/index.shtml)). 

![](https://alexpmagalhaes.github.io/ChIPseq_course/img/softclipping.png)

By default, Bowtie2 will perform a global *end-to-end read alignment*, which aligns from the first to the last base of the read. This alignment is best for reads that have already been trimmed for quality and adapters (e.g. reads where nucleotide bases of poor quality or matching adapter sequences have been removed from the ends of the reads prior to alignment). However, Bowtie2 also has a _local alignment mode_, which, in contrast to end-to-end alignment, ignores portions at the ends of the reads that do match well to the reference genome. This is referred to as **soft-clipping** and allows for a more accurate alignment. The procedure can carry a small penalty for each soft-clipped base, but amounts to a significantly smaller penalty than mismatching bases. In contrast to trimming, which removes the unwanted sequence (hard-clipping), soft-clipping retains the soft-clipped base in the sequence and simply marks it. _We will use this option since we did not trim our reads._

The command to run the alignment is simply bowtie2. Some additional arguments that we will need for aligning reads to the genome using Bowtie2 are described below:

* `-p`: number of processors/cores
* `-q`: reads are in FASTQ format
* `--local`: local alignment feature to perform soft-clipping
* `-x`: /path/to/genome_indices_directory
* `-U`: /path/to/FASTQ_file
* `-S`: /path/to/output/SAM_file

Bowtie2 does not generate log summary files. Rather this information gets printed to screen. If we want to capture that and save it in a file we can access later we can use the 2> operator. To redirect the standard error from the bowtie2 command we could do the following:


```bash
mkdir sam #folder where we will save the resulting SAM file

#Align GFP chip control
bowtie2 -p 6 -q --local \
 -x /project/pcpool_data/molmed/Genomes/mm10/mm10 \
 -U  ./GFP_Mouse_ES_merged.fastq \ 
 -S ./sam/GFP_Mouse_ES_chip.sam 2> GFP_Mouse_ES_chip_bowtie2.log
 
 #Align the Nanog chip sample
 bowtie2 -p 6 -q --local \
 -x /project/pcpool_data/molmed/Genomes/mm10/mm10 \
 -U ./Nanog_Mouse_ES_merged.fastq \ 
 -S ./sam/Nanog_Mouse_ES_chip.sam 2> Nanog_Mouse_ES_chip_bowtie2.log

```

## Alignment output: SAM/BAM file format

The output from the Bowtie2 aligner is an unsorted SAM file, also known as **Sequence Alignment/Map format**. The SAM file is a **tab-delimited text file** that contains information for each individual read and its alignment to the genome. While we will go into some features of the SAM format, the paper by [Heng Li et al.](http://bioinformatics.oxfordjournals.org/content/25/16/2078.full) provides a lot more detail on the specification.

The file begins with a **header**, which is optional. The header is used to describe source of data, reference sequence, method of alignment, etc., this will change depending on the aligner being used. Each section begins with character ‘@’ followed by [**a two-letter record type code**](https://www.samformat.info/sam-format-header). These are followed by two-letter tags and values.

Following the header is the **alignment section**. Each line corresponds to the alignment information for a single read. Each alignment line has **11 mandatory fields for essential mapping information** and a variable number of other fields for aligner-specific information. 

![](https://alexpmagalhaes.github.io/ChIPseq_course/img/sam_bam.png)

An example read mapping is displayed above. *Note that the example above spans two lines, but in the actual file it is a single line.* Let's go through the fields one at a time. 

- **`QNAME`:** Query name or read name - this is the same read name present in the header of the FASTQ file
- **`FLAG`:** numerical value providing information about read mapping and whether the read is part of a pair. Specifics regarding the FLAG values are [available](https://www.samformat.info/sam-format-flag-single).
- **`RNAME`:** is the reference sequence name, giving the chromosome to which the read maps. The example read is from chromosome 1, which explains why we see 'chr1'. 
- **`POS`:** refers to the 1-based leftmost position of the alignment. 
- **`MAPQ`:** is giving us the alignment quality, the scale of which will depend on the aligner being used. 
- **`CIGAR`:** is a sequence of letters and numbers that represent the *edits or operations* required to match the read to the reference. The letters are operations that are used to indicate which bases align to the reference (i.e. match, mismatch, deletion, insertion), and the numbers indicate the associated base lengths for each 'operation'. For example, in the SAM image below, the `100M` represents all 100 bp match the genome (no insertions, deletions, or gaps).

Now to the remaining fields in our SAM file:

![](https://alexpmagalhaes.github.io/ChIPseq_course/img/sam_bam3.png)

The next three fields are more pertinent to paired-end data. 

- **`MRNM`:** is the mate reference name. 
- **`MPOS`:** is the mate position (1-based, leftmost). 
- **`ISIZE`:** is the inferred insert size.

Finally, you have the raw sequence data from the original FASTQ file stored for each read:

- **`SEQ`:** is the raw sequence
- **`QUAL`:** is the associated quality values for each position in the read.


### Converting file format from SAM to BAM

While the SAM alignment file from Bowtie2 is human-readable, we need a BAM alignment file for downstream analysis. A BAM file is a binary equivalent version of the SAM file, in other words, the same file in a compressed format. Therefore, BAM file is not human-readable, and it is much smaller. BAM file is the typical format used in bioinformatics tools. We will use [Samtools](http://samtools.github.io) to convert the file format from SAM to BAM. Samtools is a program that consists of many utilities for working with the Sequence Alignment/Map (SAM) format. Here, we will use the `samtools view` command to convert our SAM file into its binary compressed version (BAM) and save it to file.

> NOTE: Once we generate the BAM file, we don't need to retain the SAM file anymore - we can delete it to save space.


We outline below the parameters to use with the command `samtools view`, and what each does:

* `-h`: include header in output
* `-S`: input is in SAM format
* `-b`: output BAM format
* `-o`: /path/to/output/file

> **NOTE**: You can find detailed instructions for different samtools functions and additional parameter options in this [manual](http://www.htslib.org/doc/samtools-1.2.html). 

```bash
mkdir bam

samtools view -@ 4 -h -S -b \
-o ./bam/GFP_Mouse_ES_chip.bam \
./sam/GFP_Mouse_ES_chip.sam

samtools view -@ 4 -h -S -b \
-o ./bam/Nanog_Mouse_ES_chip.bam \
./sam/Nanog_Mouse_ES_chip.sam

```


> _**NOTE:** After performing read alignment, it's useful to evaluate the mapping rate for each sample by taking look at the log files. Additionally, it is common to aggregate QC metrics and visualize them with plots using tools such as [MultiQC](http://multiqc.info). This is important to do prior to moving on to the next steps of the analysis._


To continue with the tutorial please go to [Mapping and post-processing](https://alexpmagalhaes.github.io/ChIPseq_course/mapping.md)

To go back to the home page follow this [Link](https://alexpmagalhaes.github.io/ChIPseq_course/index.md)