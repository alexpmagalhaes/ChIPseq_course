# Sequence QC and Trimming
In this section we will cover sequence quality control and trimming
Bellow are links to sections we will cover

<!-- TOC -->
* [Sequence QC and Trimming](#sequence-qc-and-trimming)
  * [2.1 Sequence QC](#21-sequence-qc)
  * [Quality control of sequence reads](#quality-control-of-sequence-reads)
    * [Unmapped read data: FASTQ file format](#unmapped-read-data-fastq-file-format)
  * [Assessing sequence read quality with FastQC](#assessing-sequence-read-quality-with-fastqc)
    * [Running FastQC](#running-fastqc-)
    * [Interpreting the FastQC HTML report](#interpreting-the-fastqc-html-report)
      * [Per base sequence quality](#per-base-sequence-quality)
      * [Sequence length distribution](#sequence-length-distribution)
      * [Sequence duplication levels](#sequence-duplication-levels)
      * [Over-represented sequences](#over-represented-sequences)
      * [K-mer content](#k-mer-content)
    * [Conclusion](#conclusion)
  * [2.2 Trimming](#22-trimming)
    * [To trim or not to trim?](#to-trim-or-not-to-trim)
<!-- TOC -->

## 2.1 Sequence QC

Read quality is the first step in all the analyses of sequenced reads. 

Fastq files can be downloaded from Sequence Read Archive (SRA).

In the folder `/project/pcpool_data/molmed/fastq/ChIP-seq/` you can find fastq files from the Chen et al. paper.


* SRR001996: GFP
* SRR001997: GFP
* SRR001998: GFP
* SRR001999: GFP


* SRR002004: nanog
* SRR002005: nanog
* SRR002006: nanog
* SRR002007: nanog
* SRR002008: nanog
* SRR002009: nanog
* SRR002010: nanog
* SRR002011: nanog

We already downloaded the fastq for you so fell free copy the relevant files to the work folder

To do so on the Terminal window make a copy of the fastq to your results folder

```bash
cp /project/pcpool_data/molmed/fastq/ChIP-seq/*.fastq.gz ./
```

If you need to download a different dataset you can download via `wget` or `curl`

One resource to help you get the right code to do it is [SRA Explorer](https://sra-explorer.info/)

The exemple bellow is for **Nanog** and **GFP** samples that we will use for this tutorial.

```bash
#DO NOT RUN ME

# curl allows a user to download the fastq files
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR002/SRR002006/SRR002006.fastq.gz -o SRR002006_Illumina_sequencing_of_Mouse_ES_nanog_genomic_fragment_library.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR001/SRR001998/SRR001998.fastq.gz -o SRR001998_Illumina_sequencing_of_Mouse_ES_GFP_genomic_fragment_library.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR002/SRR002007/SRR002007.fastq.gz -o SRR002007_Illumina_sequencing_of_Mouse_ES_nanog_genomic_fragment_library.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR002/SRR002008/SRR002008.fastq.gz -o SRR002008_Illumina_sequencing_of_Mouse_ES_nanog_genomic_fragment_library.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR001/SRR001996/SRR001996.fastq.gz -o SRR001996_Illumina_sequencing_of_Mouse_ES_GFP_genomic_fragment_library.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR002/SRR002011/SRR002011.fastq.gz -o SRR002011_Illumina_sequencing_of_Mouse_ES_nanog_genomic_fragment_library.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR002/SRR002010/SRR002010.fastq.gz -o SRR002010_Illumina_sequencing_of_Mouse_ES_nanog_genomic_fragment_library.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR002/SRR002009/SRR002009.fastq.gz -o SRR002009_Illumina_sequencing_of_Mouse_ES_nanog_genomic_fragment_library.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR001/SRR001997/SRR001997.fastq.gz -o SRR001997_Illumina_sequencing_of_Mouse_ES_GFP_genomic_fragment_library.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR002/SRR002004/SRR002004.fastq.gz -o SRR002004_Illumina_sequencing_of_Mouse_ES_nanog_genomic_fragment_library.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR002/SRR002005/SRR002005.fastq.gz -o SRR002005_Illumina_sequencing_of_Mouse_ES_nanog_genomic_fragment_library.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR001/SRR001999/SRR001999.fastq.gz -o SRR001999_Illumina_sequencing_of_Mouse_ES_GFP_genomic_fragment_library.fastq.gz

# since the sequencing was done in fragmented library you will need to merge the fastq files
# cat if used to merge the fastq files into one single fastq file
cat *_GFP_*.fastq.gz > GFP_Mouse_ES_merged.fastq.gz
cat *_nanog_*.fastq.gz > Nanog_Mouse_ES_merged.fastq.gz
```

In the project folders there is already fastq files for Nanog and GFP control fastq files, we will work with those from now on.

## Quality control of sequence reads

### Unmapped read data: FASTQ file format

The [FASTQ](https://en.wikipedia.org/wiki/FASTQ_format) file format is the default file format for sequence reads generated from next-generation sequencing technologies. This file format evolved from FASTA in that it contains sequence data, but also contains quality information. Similar to FASTA, the FASTQ file begins with a header line. The difference is that the FASTQ header is denoted by a `@` character. For a single record (sequence read) there are four lines, each of which are described below:

| Line | Description                                                                                                  |
|------|--------------------------------------------------------------------------------------------------------------|
| 1    | Always begins with '@' and then information about the read                                                   |
| 2    | The actual DNA sequence                                                                                      |
| 3    | Always begins with a '+' and sometimes the same info in line 1                                               |
| 4    | Has a string of characters which represent the quality scores; must have same number of characters as line 2 |

Let's use the following read as an example:

```
@HWI-ST330:304:H045HADXX:1:1101:1111:61397
CACTTGTAAGGGCAGGCCCCCTTCACCCTCCCGCTCCTGGGGGANNNNNNNNNNANNNCGAGGCCCTGGGGTAGAGGGNNNNNNNNNNNNNNGATCTTGG
+
@?@DDDDDDHHH?GH:?FCBGGB@C?DBEGIIIIAEF;FCGGI#########################################################
```

As mentioned previously, line 4 has characters encoding the quality of each nucleotide in the read. The legend below provides the mapping of quality scores (Phred-33) to the quality encoding characters. *Different quality encoding scales exist (differing by offset in the ASCII table), but note the most commonly used one is fastqsanger.*

 ```
 Quality encoding: !"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHI
                   |         |         |         |         |
    Quality score: 0........10........20........30........40                                
```
 
Using the quality encoding character legend, the first nucleotide in the read (C) is called with a quality score of 31, and our Ns are called with a score of 2. **As you can tell by now, this is a bad read.** 

Each quality score represents the probability that the corresponding nucleotide call is incorrect. This quality score is logarithmically based and is calculated as:

	Q = -10 x log10(P), where P is the probability that a base call is erroneous

These probability values are the results from the base calling algorithm and dependent on how much signal was captured for the base incorporation. The score values can be interpreted as follows:

| Phred Quality Score | Probability of incorrect base call | Base call accuracy |
|:--------------------|:----------------------------------:|-------------------:|
| 10	                 |              1 in 10               |               	90% |
| 20	                 |              1 in 100              |               	99% |
| 30	                 |             1 in 1000              |             	99.9% |
| 40	                 |            1 in 10,000             |            	99.99% |
| 50	                 |            1 in 100,000            |           	99.999% |
| 60	                 |           1 in 1,000,000           |          	99.9999% |

Therefore, for the first nucleotide in the read (C), there is less than a 1 in 1000 chance that the base was called incorrectly. However, for the end of the read, there is greater than 50% probability that the base is called incorrectly.

## Assessing sequence read quality with FastQC

Now we understand what information is stored in a FASTQ file, the next step is to **generate quality metrics for our sequence data**.

[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) provides a simple way to do some quality control checks on raw sequence data coming from high throughput sequencing pipelines. It contains a modular set of analyses, allowing us to have a quick glimpse of whether the data is problematic before doing any further analysis.

The main features of FastQC are:

* Imports data in FASTQ files (or BAM files)
* Evaluates the reads automatically and identifies potential issues in the data
* Generates an HTML-based quality report with graphs and tables


### Running FastQC  

To start quality checking we will use `fastqc` and make html report qc report


```bash
mkdir qc_report

fastqc --outdir ./qc_report --threads 6 ./GFP_Mouse_ES_merged.fastq.gz
fastqc --outdir ./qc_report --threads 6 ./Nanog_Mouse_ES_merged.fastq.gz
```

> _NOTE_: FastQC can also accept multiple files as input. In order to do this, you just need to separate each file name with a space. FastQC also recognizes multiple files with the use of wildcard characters.
**You may have noticed that it takes a few minutes to run this single sample through FastQC. How can we speed it up?**


**Let's take a closer look at the files generated by FastQC:**


You should see **two output files** generated:

1. The first is an **HTML file** which is a self-contained document with various graphs embedded into it. Each of the graphs evaluates different quality aspects of our data, we will discuss in more detail in this lesson.

2. Alongside the HTML file is a **zip file** (with the same name as the HTML file, but with .zip added to the end). This file contains the different plots from the report as separate image files but also contains data files which are designed to be easily parsed to allow for a more detailed and automated evaluation of the raw data on which the QC report is built.

### Interpreting the FastQC HTML report

Now we can take a look at the various metrics generated by FastQC and assess the quality of our sequencing data!

Generally it is a good idea to take a look at **basic statistics** for the sample. Keep track of the total number of reads sequenced for each sample, and make sure the read length and %GC content is as expected.

![](https://alexpmagalhaes.github.io/ChIPseq_course/img/03_fastqc_statistics.png)
	
On the left side of the report, a summary of all the modules is given. **Don't take the yellow "WARNING"s and red "FAIL"s too seriously - they should be interpreted as flags for modules to check out.** Below we will discuss some of the relevant plots to keep an eye out for when using FastQC to assess ChIP-seq data.

![](https://alexpmagalhaes.github.io/ChIPseq_course/img/03_fastqc_summary.png)

> **Please note that the evaluation of FASTQC metrics will remain the same most NGS data, however there are minor differences when working with some types like Hi-C or ATAC-seq data.**

#### Per base sequence quality
The **[Per base sequence quality](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/2%20Per%20Base%20Sequence%20Quality.html)** plot is an important analysis module in FastQC for ChIP-seq; it provides the distribution of quality scores across all bases at each position in the reads. This information helps determine whether there are any problems during the sequencing of your data. Generally, we might observe a decrease in quality towards the ends of the reads, but we shouldn't see any quality drops at the beginning or in the middle of the reads.

![](https://alexpmagalhaes.github.io/ChIPseq_course/img/03_fastqc_sequence_quality.png)

Based on the sequence quality plot of our sample, the majority of the reads have high quality. Particularly, the quality remains high towards the end of the read and we do not observe any unexpected quality drop in the middle of the sequence. If that happens, we will need to contact the sequencing facility for further investigation. 

#### Sequence length distribution
The **[Sequence length distribution](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/7%20Sequence%20Length%20Distribution.html)** module generates a graph showing the distribution of fragment sizes in the file which was analysed. Most high-throughput sequencers generate reads of uniform length, and this will produce a simple graph showing a peak only at one size, but for variable length FastQ files this will show the relative amounts of each different size of sequence fragment.

![](https://alexpmagalhaes.github.io/ChIPseq_course/img/03_fastqc_seq_length.png)

We observe a unique sequence length (36 bp) in our data, which is as expected. In some cases of ChIP-seq data, we have observed variable sequence lengths that is not due to the sequencer. Rather, this can be a result of the adapter removal in the previous steps suggesting that there was adapter contamination in your sample.

#### Sequence duplication levels

The **[Sequence duplication level](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/8%20Duplicate%20Sequences.html)** plot indicates what proportion of your library corresponds to duplicates. Duplicates are typically removed in the ChIP-seq workflow (even if there is a chance that they are biological duplicates). Therefore, if we observe a large amount of duplication, the number of reads available for mapping and peak calling will be reduced. We don't observe concerning duplication levels in our data.

![](https://alexpmagalhaes.github.io/ChIPseq_course/img/03_fastqc_duplication_level.png)

#### Over-represented sequences

**[Over-represented sequences](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/9%20Overrepresented%20Sequences.html)** could come from actual biological significance, or biases introduced during the sequencing. With ChIP-seq, you expect to see over-represented sequences in the immuno-precipitation sample, because that's exactly what you're doing - enriching for particular sequences based on binding affinity. However, lack of over-represented sequences in FastQC report doesn't mean you have a bad experiment. If you observe over-represented sequences in the input sample, that usually suggests some bias in the protocol to specific regions. 


#### K-mer content

The **[Kmer content](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/11%20Kmer%20Content.html)** module examines whether there is positional bias for any small fragment of sequence. It measures the number of each 7-mer at each position in your library and then uses a binomial test to look for significant deviations from an even coverage at all positions. Any K-mers with positionally biased enrichment are reported and the top 6 most biased K-mers are additionally plotted to show their distribution.

![](https://alexpmagalhaes.github.io/ChIPseq_course/img/03_fastqc_kmer_content.png)


### Conclusion

Based on the metrics and plots generated by FastQC, we have a poor quality sample to move forward with. However, these libraries are very old, and they are a reflexion of sequencing technology from 2010s. Typically, you would generate these reports for all samples in your dataset and keep track of any samples that did not fare as well. We focused on the interpretation of selected plots that are particularly relevant for ChIP-seq, but if you would like to go through the remaining plots and metrics, FastQC has a well-documented [manual page](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) with [more details](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/) (under the Analysis Modules directory) about all the plots in the report.



## 2.2 Trimming

### To trim or not to trim?

Trimming is the process of removing unwanted sequence prior to sequence alignment. In general, trimming is an optional step. The 3' ends of the sequence may contain part of the Illumina sequencing adapter. This adapter contamination may prevent the reads from aligning to the reference genome correctly, thus adversely impacting the downstream analysis. You could evaluate potential adapter contamination either from the FastQC report (in "Overrepresented sequences" or "Adapter Contamination" sections), or from the size distribution of your sequencing libraries. If you suspect that your reads are contaminated with adapters, you should run an adapter removal tool. We list some additional considerations below whether trimming is needed:

  * If the read length is 25bp, there is no need to trim - adapter sequences will not be included in reads of inserts >25 bp.
  * If you perform trimming, then there is no need to use soft-clipping during the alignment.
  * After trimming, a minimum read length of 25bp should be imposed, as reads smaller than this are hard to align accurately.
  * Should you choose not to trim reads, you will need to use `--local` when running Bowtie2 - this will perform "soft-clipping" to ignore parts of the reads that are of low quality.

There are many programs to do QC, and many specific tools for each one. For now, we are going to focus on Cutadapt  

```bash
#DO NOT RUN ME
cutadapt 
    -q 30 #defines the quality
    --minimum-length 25 #minimum length of the read
    --cores 6 #number of cores
    --adapter $adaptors #adapter sequence
    --output $fq1_trimmed #output file (usually is good to attach a reference to the name like "_trimmed")
    $FQ #input file

```

> _NOTE_:  It is important to note that the libraries covered in this tutorial are 25bp and thus no not contain any adapters, so we will skip this step.

To continue with the tutorial please go to [Mapping and post-processing](https://alexpmagalhaes.github.io/ChIPseq_course/coverage.md)

To go back to the home page follow this [Link](https://alexpmagalhaes.github.io/ChIPseq_course/index.md)