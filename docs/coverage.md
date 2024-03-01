## Peak visualization

In this section we will cover bigwig conversion for visualization
Bellow are links to sections we will cover

<!-- TOC -->
  * [Peak visualization](#peak-visualization)
  * [4.1 File formats for peak visualization](#41-file-formats-for-peak-visualization)
    * [BedGraph format](#bedgraph-format)
    * [Wiggle and bigWig formats](#wiggle-and-bigwig-formats)
  * [4.2 Creating bigWig files](#42-creating-bigwig-files)
    * [Normalization](#normalization)
    * [bamCoverage from deepTools](#bamcoverage-from-deeptools)
    * [bamCompare from deepTools](#bamcompare-from-deeptools)
  * [4.3 Profile plots](#43-profile-plots)
  * [Evaluating signal in Nanog binding sites](#evaluating-signal-in-nanog-binding-sites)
    * [1. Create the matrix](#1-create-the-matrix)
    * [2. Drawing the profile plot](#2-drawing-the-profile-plot)
<!-- TOC -->

Now that we have identified regions in the genome that are enriched through some interaction with Nanog, we can **take those regions and visually assess the amount of signal observed**. This can be done in one of two ways:

1. Uploading the data to a genome viewer such as the Broad's [Integrative Genome Viewer (IGV)](https://software.broadinstitute.org/software/igv/) or the [UCSC Genome Browser](https://genome.ucsc.edu/cgi-bin/hgGateway), to explore the genome and specific loci for pileups.
2. Create profile plots and heatmaps to look at the signal aggregated over all binding sites. 

**In order to perform any of assessments described above, you will need a file in the appropriate format.** 

The goal of this lesson is to introduce you to **different file formats used for ChIP-seq data visualization and how to generate these files** using [`deepTools`](https://deeptools.readthedocs.io/en/develop/index.html).


## 4.1 File formats for peak visualization

There are several different types of file formats that can hold data associated with high-throughput sequencing data, these file formats have a distinct structure and hold specific types of data. We have already encountered:

* the sequence data format - FASTQ
* the alignment file formats - SAM and BAM
* the peak call format - BED, narrowPeak 

In this section we want to introduce you to a few additional formats that can be used for visualizing peaks.

The commonality among these file formats is that they represent the peak location in a manner similar to the BED format (shown below). 

<p align="center">
  <img src="https://alexpmagalhaes.github.io/ChIPseq_course/img/bed.png" width="500" alt="">
</p>


### BedGraph format

In addition to accommodating peak calls (discrete data), the BedGraph format also allows display of continuous data as a track on a genome browser. This display type is useful for including and plotting some quantitative information, e.g. intensity. For the purposes of visualization, bedGraph and bigWig are practically interchangeable, with the bigWig file being a lot smaller for a given dataset.

<p align="center">
  <img src="https://alexpmagalhaes.github.io/ChIPseq_course/img/bedgraph.png" width="650" alt="">
</p>

### Wiggle and bigWig formats

The Wiggle format (wig) also allows the display of continuous data. This format is "line" oriented, and has declaration lines and data lines. It also requires a separate wiggle track definition line as a header. There are two options for how data in wiggle files are represented: variableStep and fixedStep. These formats were developed to allow the file to be written as compactly as possible.

<p align="center">
  <img src="https://alexpmagalhaes.github.io/ChIPseq_course/img/wiggle.png" width="650" alt="">
</p>

The bigWig format is an indexed *binary* form of the wiggle file format, and is useful for large amounts of dense and continuous data to be displayed in a genome browser as a graphical track. As mentioned above, the visual representation of this format is very similar to bedGraph.


## 4.2 Creating bigWig files

For this workshop, we will focus on creating bigWig files, as we will be using them in the next lesson for qualitative assessment. BigWig files have a much smaller data footprint compared to BAM files, especially as your bin size (a parameter described below) increases. The general procedure is to take our **alignment files (BAM) and convert them into bigWig files**, and we will do this using [`deepTools`](http://deeptools.readthedocs.org/en/latest/content/list_of_tools.html). The coverage is calculated as the number of reads per bin, where bins are short consecutive counting windows of a defined size. It is possible to extend the length of the reads to better reflect the actual fragment length. 

<p align="center">
<img src="https://alexpmagalhaes.github.io/ChIPseq_course/img/bam_to_bigwig.png" width="700" alt="">
</p>

*Image acquired from the [deepTools documentation](http://deeptools.readthedocs.io/en/latest/content/tools/bamCoverage.html?highlight=bigwig)*

**`deepTools` is a suite of Python tools developed for the efficient analysis of high-throughput sequencing data**, such as ChIP-seq, RNA-seq, or MNase-seq. `deepTools` has a wide variety of commands that go beyond what we will cover in this workshop. We encourage you to look through the documentation and explore more on your own time. 

First lets make a folder for the files

```bash
cd ..
mkdir -p visualization/bigWig
```

We then need to **create an index file for the BAM file**. Often, when working with BAM files you will find that many tools require an index (an associated `.bai` file). You can think of an index similar to that which is located at the back of a textbook - when you are interested in a particular subject, you look for the keyword in the index and identify the pages that contain the relevant information. Similarly, indexing the BAM file aims to achieve fast retrieval of alignments overlapping a specified region without going through the whole alignment file. Essentially, a `bai` file along with the `bam` ensures that downstream applications are able to use the information with the `bam` file much more speedily.

We will use [SAMtools](http://samtools.sourceforge.net/) again, specifically the **`samtools index`** command, to index the BAM file.

Create an index for the `Nanog_Mouse_ES_chip.filtered.bam` and `GFP_Mouse_ES_chip.filtered.bam` file that we created in earlier

```bash
samtools index ./bam/GFP_Mouse_ES_chip.filtered.bam
samtools index ./bam/Nanog_Mouse_ES_chip.filtered.bam
```

Finally, let's make sure we have the required modules loaded to use `deepTools`:


### Normalization

The methods for bigWig creation (`bamCoverage` and `bamCompare`) allows for normalization, which is great if we want **to compare different samples to each other, and they vary in terms of sequencing depth**. DeepTools offers different **methods of normalization** as listed below, each is performed per bin. The default is no normalization.

> NOTE: We will not normalize the data we are working with because we are following the methods described in [Baizabal, 2018](https://doi.org/10.1016/j.neuron.2018.04.033). However, it is highly recommended to choose one of the methods described.

* Reads Per Kilobase per Million mapped reads (RPKM)
  * number of reads per bin / (number of mapped reads (in millions) * bin length (kb))
* Counts per million (CPM); this is similar to CPM in RNA-seq
  * number of reads per bin / number of mapped reads (in millions)
* Bins Per Million mapped reads (BPM); same as TPM in RNA-seq
  * number of reads per bin / sum of all reads per bin (in millions)
* Reads per genomic content (RPGC)
  * number of reads per bin / scaling factor for 1x average coverage 
  * scaling factor is determined from the sequencing depth: (total number of mapped reads * fragment length) / effective genome size
  * this option requires an effectiveGenomeSize
  

**Spike-in normalization**

Another option for normalization is to **normalize each sample using a scale factor/ normalization factor**. By default, the scale factor is set to 1 in deepTools. In an earlier lesson, we described the spike-in strategy and how to compute a normalization factor for each individual sample. Those values can be used when creating bigWig files. The `--scaleFactor` parameter takes in the user provided value and each bin is multiplied by this value. **If normalizing with a scale factor, be sure that none of the other normalization methods are applied.**


### bamCoverage from deepTools

This command takes a **BAM file as input** and evaluates which areas of the genome have reads associated with them, i.e. how much of the genome is "covered" with reads. The coverage is calculated as the number of reads per bin, where bins are short consecutive sections of the genome (bins) that can be defined by the user. The **output of this command can be either a bedGraph or a bigWig file**. We will be generating a bigWig file, since that format has a much smaller data footprint, especially as the bin size increases.

These are some parameters of bamCoverage that are worth considering:
* `normalizeUsing`: Possible choices: RPKM, CPM, BPM, RPGC. By default, no normalization is applied.
* `binSize`: size of bins in bases (default is 50)
* `--effectiveGenomeSize`: the portion of the genome that is mappable. It is useful to consider this when computing your scaling factor.
* `smoothLength`: defines a window, larger than the `binSize`, to average the number of reads over. This helps produce a more continuous plot.
* `centerReads`: reads are centered with respect to the fragment length as specified by `extendReads`. This option is useful to get a sharper signal around enriched regions.


We will be using the bare minimum of parameters as shown in the code below. We decrease the bin size to increase the resolution of the track (this also means larger file size). If you are interested, feel free to test out some of the other parameters to create different bigWig files. You can load them into a genome viewer like IGV and observe the differences.

Let's **create a bigWig file with a `binSize` of 20.

```bash
bamCoverage -p 6 -b ./bam/Nanog_Mouse_ES_chip.filtered.bam \
-o ./visualization/bigWig/Nanog_Mouse_ES_chip.bw \
--binSize 20
```


### bamCompare from deepTools

As an alternate to calculating genome coverage with `bamCoverage`, we could use `bamCompare`. `bamCompare` will **create a bigWig file in which we compare the ChIP against the input**. The command is quite similar to `bamCoverage`, except that it requires two files as input (`b1` and `b2`). Below, we show you an example of how you would run `bamCompare`. The default `--operation` used to compare the two samples is the **log2 ratio**, however you also have the option to add, subtract and average. Any of the parameters described above for `bamCoverage` can also be used. 

```bash
## DO NOT RUN

bamCompare -p 6 -b1 ./bam/Nanog_Mouse_ES_chip.filtered.bam \
-b2 ./bam/GFP_Mouse_ES_chip.filtered.bam \
-o ./visualization/bigWig/Nanog_Mouse_ES_chip_bc.bw \
--binSize 20
```

You can find more details about the difference between `bamCompare` and `bamCoverage` [linked here](https://deeptools.readthedocs.io/en/develop/content/help_faq.html#when-should-i-use-bamcoverage-or-bamcompare).

## 4.3 Profile plots
The profile plot allows us to **evaluate read density over sets of genomic regions**. Typically, these regions are genes (start and end coordinates), but any other regions defined in the BED file will work. 

We will use the `plotProfile` command that is part of the [deepTools suite](https://deeptools.readthedocs.io/en/develop/content/tools/plotProfile.html), to create our profile plots. To use this, a matrix generated by `computeMatrix` is required. You will need to **create a matrix for every new plot you generate**, so this step of the workflow can often take some time. 

> *NOTE*: In this lesson, we teach you how to create a matrix for the first figure, and then we provide the pre-computed matrices for you to use to save time.

<p align="center">
<img src="https://alexpmagalhaes.github.io/ChIPseq_course/img/computeMatrix_overview.png" width="700" alt="">
</p>

_Image source: [deepTools documentation](https://deeptools.readthedocs.io/en/develop/content/tools/computeMatrix.html?highlight=computeMatrix)_

> *NOTE:* The matrix we generate can also be used as input to [`plotHeatmap`](https://deeptools.readthedocs.io/en/develop/content/tools/plotHeatmap.html?highlight=plotheatmap). Heatmaps are also a very popular data visualization option for peak call data. We will not be plotting heatmaps in this lesson, but encourage you to explore the command and the numerous parameters available for optimization.


TThe `computeMatrix` command can take some time, so we want to take advantage of the multi-threading.

```bash
mkdir -p ./visualization/profile
```

## Evaluating signal in Nanog binding sites

### 1. Create the matrix
The first step in generating the profile plot is to create the matrix. The `computeMatrix` command accepts multiple bigWig files and multiple region files (BED format) to create a count matrix. The command can also filter and sort regions according to their scores. For each window, `computeMatrix` will calculate scores based on the read density values in the bigWig files.

Below we describe the **parameters** we will be using:

* `reference-point`: The reference point for plotting. Here, we use the center of the consensus peaks (default is TSS).
* `-b`, `a`: Specify a window around the reference point (before and after). We have used +/- 4000 bp. 
   > **How do I choose the window size?** For narrow peaks, we can use a smaller window as we expect a more punctuate binding profile. Broader peaks require larger windows to capture the whole profile shape. For either profile, you may need to play around with the window size and see what works best with your data. A good starting point is +/- 2kb.
* `-R`: The region file will be the BED file we generated for WT replicate overlap.
* `-S`: The list of bigWig files (WT replicates), that we have generated for you.
* `--skipZeros`: Do not include regions with only scores of zero
* `-o`: output file name
* `-p`: number of cores


Let's create a matrix for the WT replicates:

```bash
computeMatrix reference-point --referencePoint center \
-b 3000 -a 3000 \
-R ./macs2/nanog_summits.bed \
-S ./visualization/bigWig/Nanog_Mouse_ES_chip.bw \
--skipZeros \
-o ./visualization/profile/Nanog_Mouse_ES_chip.gz \
-p 6
```
> _Runtime estimate: 8-10 minutes_

### 2. Drawing the profile plot
Once you have computed the matrix, you can create the **profile plot**. First, make a directory designated for the figures we will be creating, and then we will run `plotProfile`. _The `plotProfile` command will take a shorter amount of time to run._ 

> **NOTE:** `plotProfile` has many options to optimize your figure, including the ability to change the type of lines plotted, and to plot by group rather than sample. We encourage you to explore the [documentation](https://deeptools.readthedocs.io/en/develop/content/tools/plotProfile.html?highlight=plotProfile) to find out more detail.

```bash

# Create figures directory under visualization
mkdir ./visualization/figures

# Plot the profiles
plotProfile -m ./visualization/profile/Nanog_Mouse_ES_chip.gz \
-out ./visualization/figures/plot1_Nanog_Mouse_ES.png \
--regionsLabel "" \
--perGroup \
--colors red blue \
--samplesLabel "Nanog_Mouse_ES" \
--refPointLabel "Nanog binding sites"
```

The figure should look like the one displayed below. We observe that the **replicate 2 has a much higher signal** present in these regions. This is not uncommon in ChIP-seq data. There will likely be one replicate that exhibits stronger signal. What is encouraging to see is that there is a **decent amount of signal in both replicates**, so we have some confidence in the regions we identified.


<p align="center">
<img src="https://alexpmagalhaes.github.io/ChIPseq_course/img/plot1_Nanog_Mouse_ES.png" width="500" alt="">
</p>


To continue with the tutorial please go to [ChIP quality assessment](https://alexpmagalhaes.github.io/ChIPseq_course/chipseeker)

To go back to the home page follow this [Link](https://alexpmagalhaes.github.io/ChIPseq_course/index.md)