# Mapping and post processing

In this section we will cover mapping and post-processing of alignment files
Bellow are links to sections we will cover

<!-- TOC -->
* [Mapping and post processing](#mapping-and-post-processing)
  * [3.1 Mapping aka read alignment](#31-mapping-aka-read-alignment)
  * [Alignment](#alignment)
  * [Alignment output: SAM/BAM file format](#alignment-output-sambam-file-format)
    * [Converting file format from SAM to BAM](#converting-file-format-from-sam-to-bam)
  * [3.2 Filtering reads aka postprocessing](#32-filtering-reads-aka-postprocessing)
    * [Multi-mapping reads](#multi-mapping-reads)
    * [Duplicate reads](#duplicate-reads)
    * [Filtering workflow](#filtering-workflow)
    * [1. Sort BAM files by genomic coordinates](#1-sort-bam-files-by-genomic-coordinates)
    * [2. Filter the reads to keep only uniquely mapping reads](#2-filter-the-reads-to-keep-only-uniquely-mapping-reads)
<!-- TOC -->


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

<p align="center">
	<img src="https://alexpmagalhaes.github.io/ChIPseq_course/img/Alignment_errors.png" width="600" alt="">
</p>



In theory, this sounds like a very simple case of string matching. We take the sequence read and figure out where it originated from in the reference genome. However, in practice, this is **actually quite difficult!** This is because:

* The reference genome we are searching is large and complex (e.g. the human genome is ~3,200,000,000bp). 
* By contrast, the reads we are searching for are much smaller (50-150bp), and they are on the range of millions for a given sample.
* We have to consider non-exact matching of the read to the reference due to natural variation and sequencing errors.
* We have to consider non-unique alignment due to the short length of reads and high percentage of repetitive regions in the genome (e.g. repetitive regions = >50% of the human genome).

There are many different tools that have been developed for alignment of next-generation sequencing data, and some that are more suitable to different technologies. A popular tool commonly used with ChIP-seq data, and the one that we will be using in this workshop is [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml).
Bowtie2 is a fast and accurate alignment tool that supports gaped, local and paired-end alignment modes and works best for reads that are **at least 50 bp** (shorter read lengths should use [Bowtie1](http://bowtie-bio.sourceforge.net/index.shtml)). 

<p align="center">
	<img src="https://alexpmagalhaes.github.io/ChIPseq_course/img/softclipping.png" width="700" alt="">
</p>



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

![SAM1](https://alexpmagalhaes.github.io/ChIPseq_course/img/sam_bam.png)

An example read mapping is displayed above. *Note that the example above spans two lines, but in the actual file it is a single line.* Let's go through the fields one at a time. 

- **`QNAME`:** Query name or read name - this is the same read name present in the header of the FASTQ file
- **`FLAG`:** numerical value providing information about read mapping and whether the read is part of a pair. Specifics regarding the FLAG values are [available](https://www.samformat.info/sam-format-flag-single).
- **`RNAME`:** is the reference sequence name, giving the chromosome to which the read maps. The example read is from chromosome 1, which explains why we see 'chr1'. 
- **`POS`:** refers to the 1-based leftmost position of the alignment. 
- **`MAPQ`:** is giving us the alignment quality, the scale of which will depend on the aligner being used. 
- **`CIGAR`:** is a sequence of letters and numbers that represent the *edits or operations* required to match the read to the reference. The letters are operations that are used to indicate which bases align to the reference (i.e. match, mismatch, deletion, insertion), and the numbers indicate the associated base lengths for each 'operation'. For example, in the SAM image below, the `100M` represents all 100 bp match the genome (no insertions, deletions, or gaps).

Now to the remaining fields in our SAM file:

![SAM2](https://alexpmagalhaes.github.io/ChIPseq_course/img/sam_bam3.png)

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

## 3.2 Filtering reads aka postprocessing


A key issue when working with a ChIP-seq data is to **move forward with only the uniquely mapping reads**.  Allowing for multi-mapped reads increases the number of usable reads and the sensitivity of peak detection; however, the number of false positives may also increase [[1]](https://www.ncbi.nlm.nih.gov/pubmed/21779159/). To increase our confidence in peak calling and improve data reproducibility, we need to **filter out both multi-mapping reads and duplicate reads**.

### Multi-mapping reads

Multi-mapping reads are reads that are mapping to multiple loci on the reference genome.

<p align="center">
 <img src="https://alexpmagalhaes.github.io/ChIPseq_course/img/Multimapping_reads.png" width="500" alt="">
</p>


### Duplicate reads

Duplicate reads are reads that map at the exact same location, with the same coordinates and the same strand. These duplicates can arise from experimental artifacts, but can also contribute to genuine ChIP-signal.

* **The bad kind of duplicates:** If initial starting material is low, this can lead to over-amplification of this material before sequencing. Any biases in PCR will compound this problem and can lead to artificially enriched regions. 
* **The good kind of duplicates:** You can expect some biological duplicates with ChIP-seq since you are only sequencing a small part of the genome. This number can increase if your depth of coverage is excessive or if your protein only binds to few sites. If there are a good proportion of biological duplicates, removal can lead to an underestimation of the ChIP signal. 
    
To get an idea on **what to expect in terms of duplication rate**, we encourage you to take a look at the [ENCODE quality metrics for complexity](https://www.encodeproject.org/data-standards/terms/#library). Different metrics are described and there is also a table which describes how to classify a sample as good, moderate or bad, based on these values.

<p align="center">
 <img src="https://alexpmagalhaes.github.io/ChIPseq_course/img/Duplicate_reads.png" width="500" alt="">
</p>

> #### Some additional notes on duplicates
> Most peak calling algorithms also implement methods to deal with duplicate reads. While they are commonly removed prior to peak calling, another option is to leave them now and deal with them later. **Skip the duplicate filtering at this step if**:
> * You are planning on performing a differential binding analysis.
> * You are expecting binding in repetitive regions (also, use paired-end sequencing) 
> * You have included [UMIs](https://www.illumina.com/techniques/sequencing/ngs-library-prep/multiplexing/unique-molecular-identifiers.html) into your experimental setup.


### Filtering workflow

The older version of Bowtie2 had an argument that allowed us to easily perform filtering during the alignment process. But the latest Bowtie2 does not have this option. As a result, the filtering will be done with [samtools](https://www.htslib.org/).
This **lesson will consist of two steps**:

1. Sort BAM files by genomic coordinates (using `samtools`).
2. Filter the reads to keep only uniquely mapping reads (using `gatk`). This will also remove any unmapped reads.


### 1. Sort BAM files by genomic coordinates

Before we can do the filtering, we need to sort our BAM alignment files by genomic coordinates (instead of by name). To perform the sorting, we could use [Samtools](http://www.htslib.org/), a tool we previously used when converting our SAM file to a BAM file. 

```bash
samtools sort -@ 4 -O BAM -o ./bam/Nanog_Mouse_ES_chip.sorted.bam ./bam/Nanog_Mouse_ES_chip.bam

samtools sort -@ 4 -O BAM -o ./bam/GFP_Mouse_ES_chip.sorted.bam ./bam/GFP_Mouse_ES_chip.bam
```

### 2. Filter the reads to keep only uniquely mapping reads

Next, we can filter the sorted BAM files to keep only uniquely mapping reads. We will use the `gatk` command with the following parameters:



```bash
gatk MarkDuplicates \
  --java-options "-XX:ParallelGCThreads=6 -Xmx8g" \
  --INPUT ./bam/GFP_Mouse_ES_chip.sorted.bam \
  --OUTPUT ./bam/GFP_Mouse_ES_chip.filtered.bam \
  --METRICS_FILE GFP_Mouse_ES_chip.dedup-metrics.txt \
  --REMOVE_DUPLICATES true \
  --VALIDATION_STRINGENCY LENIENT
  
gatk MarkDuplicates \
  --java-options "-XX:ParallelGCThreads=6 -Xmx8g" \
  --INPUT ./bam/Nanog_Mouse_ES_chip.sorted.bam \
  --OUTPUT ./bam/Nanog_Mouse_ES_chip.filtered.bam \
  --METRICS_FILE Nanog_Mouse_ES_chip.dedup-metrics.txt \
  --REMOVE_DUPLICATES true \
  --VALIDATION_STRINGENCY LENIENT

```


> ### Filtering out Blacklisted Regions
> Although we do not perform this step, it is common practice to apply an additional level of filtering to our BAM files. That is, we remove alignments that occur with defined Blacklisted Regions. **We will filter out blacklist regions post-peak calling.**
> 
> Blacklisted regions represent artifact regions that tend to show artificially high signal (excessive unstructured anomalous reads mapping). These regions are often found at specific types of repeats such as centromeres, telomeres and satellite repeats and typically appear uniquely mappable so simple mappability filters applied above do not remove them. The ENCODE and modENCODE consortia have compiled blacklists for various species and genome versions including human, mouse, worm and fly. These blacklisted regions (coordinate files) can be filtered out from our alignment files before proceeding to peak calling.
> 
> If we wanted to filter blacklist regions at this point in our workflow, we would use the following code:
> 
> ``` 
> # DO NOT RUN
> bedtools intersect -v Nanog_Mouse_ES_chip.filtered.bam -b mm10-blacklist.v2.bed > Nanog_Mouse_ES_chip.filtered.blacklist_filtered.bam
> ```
> 
> _bedtools is a suite of tools that we will discuss in more detail in a later lesson when blacklist filtering is applied._




To continue with the tutorial please go to [Peak Calling](https://alexpmagalhaes.github.io/ChIPseq_course/peak_calling)

To go back to the home page follow this [Link](https://alexpmagalhaes.github.io/ChIPseq_course/index)