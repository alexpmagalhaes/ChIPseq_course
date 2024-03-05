## Peak annotation

In this section we will cover peak annotation, comparison and visualization in the context of R

Bellow are links to sections we will cover

<!-- TOC -->
  * [Peak annotation](#peak-annotation)
  * [Important considerations](#important-considerations)
    * [Bioconductor Resources for ChIP-Seq](#bioconductor-resources-for-chip-seq)
    * [General Purpose Resources](#general-purpose-resources)
  * [5.1 Peak annotation](#51-peak-annotation)
    * [Read the peak files from MACS2](#read-the-peak-files-from-macs2)
    * [Annotation with TxDb file](#annotation-with-txdb-file)
    * [Write resulting table](#write-resulting-table)
  * [5.2 GO term enrichment analysis](#52-go-term-enrichment-analysis)
    * [Make lists of gene Ids identified in peak annotation step](#make-lists-of-gene-ids-identified-in-peak-annotation-step)
    * [Run enrichment analysis and plot](#run-enrichment-analysis-and-plot)
  * [5.3 Motif analysis](#53-motif-analysis)
    * [Extracting sequences under peaks](#extracting-sequences-under-peaks)
    * [Writing to FASTA file](#writing-to-fasta-file)
  * [5.4 MEME-ChIP analysis](#54-meme-chip-analysis)
    * [Parsing back FIMO results](#parsing-back-fimo-results)
    * [FIMO to R](#fimo-to-r)
    * [FIMO to valid GFF3](#fimo-to-valid-gff3)
    * [Scanning for known motifs](#scanning-for-known-motifs-)
    * [Get motifs from JASPAR with TFBStools](#get-motifs-from-jaspar-with-tfbstools)
    * [Motif scanning with motifmatchr](#motif-scanning-with-motifmatchr)
    * [Exporting motif matches](#exporting-motif-matches)

<!-- TOC -->

## Important considerations

Effective analysis of ChIP-seq data in R can be achieved with many different methods
The following list of methods and resources are available:

### Bioconductor Resources for ChIP-Seq
#### General Purpose Resources

- `GenomicRanges` ([Link](http://bioconductor.org/packages/release/bioc/html/GenomicRanges.html)): High-level infrastructure for range data
- `Rsamtools` ([Link](http://bioconductor.org/packages/release/bioc/html/Rsamtools.html)): BAM support
- `DiffBind` ([Link](http://www.bioconductor.org/packages/release/bioc/html/DiffBind.html)): Differential binding analysis of ChIP-Seq peak data
- `rtracklayer` ([Link](http://bioconductor.org/packages/release/bioc/html/rtracklayer.html)): Annotation imports, interface to online genome browsers
- `DESeq` ([Link](http://bioconductor.org/packages/release/bioc/html/DESeq.html)): RNA-Seq analysis
- `edgeR` ([Link](http://bioconductor.org/packages/release/bioc/html/edgeR.html)): RNA-Seq analysis
- `chipseq` ([Link](http://bioconductor.org/packages/release/bioc/html/chipseq.html)): Utilities for ChIP-Seq analysis
- `ChIPpeakAnno` ([Link](http://bioconductor.org/packages/release/bioc/html/ChIPpeakAnno.html)): Annotating peaks with genome context information
- `MotifDb` ([Link](http://www.bioconductor.org/packages/release/bioc/html/MotifDb.html)): Collection of motif databases
- `motifStack` ([Link](http://www.bioconductor.org/packages/release/bioc/html/motifStack.html)): Stacked logo plots
- `PWMEnrich` ([Link](http://www.bioconductor.org/packages/release/bioc/html/PWMEnrich.html)): Position Weight Matrix (PWM) enrichment analysis
- `Bioc Workflow` ([Link](http://www.bioconductor.org/help/workflows/generegulation/)): Overview of resources
- ...

Key reference: Pepke *et al*. (2009) Computation for ChIP-seq and RNA-seq studies. [_Nature Methods_](http://www.ncbi.nlm.nih.gov/pubmed/19844228) 6(11s): S22-S32.

For this we will need to load the following libraries

```{r}

library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(clusterProfiler)
library(org.Mm.eg.db)
library(BSgenome)
library(BSgenome.Mmusculus.UCSC.mm10)
library(GenomicRanges)
library(rtracklayer)
library(JASPAR2020)
library(TFBSTools)
library(motifmatchr)

```
Then we will set the working directory to the results folder
and set the Mouse mm10 TxDB annotation as an object
and get a list of genes

```{r}
#set the working directory
setwd("Whaever folder you choose") 

#set Mouse mm10 TxDB annotation as an object
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGenet

#list of genes
allGeneGR <- genes(TxDb.Mmusculus.UCSC.mm10.knownGene)
```

## 5.1 Peak annotation

Now we will be using `ChIPseeker` to annotate the peaks 

### Read the peak files from MACS2

```{r, eval=FALSE}
#read macsPeak file from narrowPeak (BED file)
macs_peak <- readPeakFile("macs2/nanog_peaks.narrowPeak")

#read macsPeak file from excel file

macsPeaks <- "macs2/nanog_peaks.xls"
macsPeaks_DF <- read.delim(macsPeaks,comment.char="#")
macsPeaks_GR <- GRanges(seqnames=macsPeaks_DF[,"chr"],
                        IRanges(macsPeaks_DF[,"start"],macsPeaks_DF[,"end"]))
mcols(macsPeaks_GR) <- macsPeaks_DF[,c("abs_summit", "fold_enrichment")]
```

### Annotation with TxDb file

```{r, eval=FALSE}

#annotate the peaks
peakAnno <- annotatePeak("macs2/nanog_peaks.narrowPeak", #read the macs2 peak file
 TxDb=txdb, #set annotation file 
 annoDb="org.Mm.eg.db", #set annotation database so we can have gene names
 tssRegion=c(-1000, 500)) #define the TSS region 
```

### Write resulting table
```{r}
#define peak annotation file as data frame
df <- as.data.frame(peakAnno)

#write it as CSV file
write.csv(df, "ChIPseeker/nanog_peaks.peaks.annotated.csv")
```

## 5.2 GO term enrichment analysis
The following will perform Gene Ontology (GO) term enrichment analysis for each annotated peak set.


### Make lists of gene Ids identified in peak annotation step
```{r, eval=FALSE}
#store gene Ids in list
gene_ids <- df$geneId

#store mouse mm10 gene Ids in list
allGeneIDs <- allGeneGR$gene_id
```

### Run enrichment analysis and plot

```{r, eval=FALSE}
#get GO enrichment for *Biological process*
GO_result <- enrichGO(gene = gene_ids, 
                      universe = allGeneIDs,
                      OrgDb = org.Mm.eg.db,
                      ont = "BP")

dotplot(GO_result)
```

Plot should look something like this:

<p align="center">
<img src="https://alexpmagalhaes.github.io/ChIPseq_course/img/GoResult.png" width="500" alt="">
</p>

## 5.3 Motif analysis

Motif analysis is useful for identifying the DNA-binding motifs in the ChIP-seq peaks. 
When the motif of the ChIPed protein is already known, motif analysis provides validation of the success of the experiment. When the motif is not known beforehand, identifying a centrally located motif in a large fraction of the peaks is an indication of a successful experiment.


### Extracting sequences under peaks
Enrichment analysis of known DNA binding motifs or *de novo* discovery of novel motifs requires the DNA sequences of the identifed peak regions. To parse the corresponding sequences from the reference genome, the `getSeq` function from the *Biostrings* package can be used.  The following example parses the sequences for each peak set and saves the results to  separate  FASTA files, one  for  each  peak  set. In  addition, the sequences in  the FASTA files are  ranked  (sorted) by increasing p-values as expected by some motif discovery tools, such as *BCRANK*.

First we need to load the BSgenome object for the genome we are working on, UCSC’s mm10 build for the mouse genome, BSgenome.Mmusculus.UCSC.mm10.

```{r, eval=FALSE}
macsSummits_GR <- GRanges(seqnames(macsPeaks_GR), IRanges(macsPeaks_GR$abs_summit,
    macsPeaks_GR$abs_summit), score = macsPeaks_GR$fold_enrichment)
macsSummits_GR <- resize(macsSummits_GR, 100, fix = "center")
```

Once we have re-centered our peaks we can use the getSeq function with our GRanges of resized common peaks and the BSgenome object for mm10.

The getSeq function returns a DNAStringSet object containing sequences under peaks.
```{r, eval=FALSE}
peaksSequences <- getSeq(BSgenome.Mmusculus.UCSC.mm10, macsSummits_GR)
names(peaksSequences) <- paste0(seqnames(macsSummits_GR), ":", start(macsSummits_GR),
    "-", end(macsSummits_GR))
```

### Writing to FASTA file
The writeXStringSet function allows the user to write DNA/RNA/AA(aminoacid)StringSet objects out to file.

By default the writeXStringSet function writes the sequence information in FASTA format (as required for MEME-ChIP).

```{r, eval=FALSE}
writeXStringSet(peaksSequences, file = "Nanog_peaks.fa")
```

## 5.4 MEME-ChIP analysis

Now the file "mycMel_rep1.fa" contains sequences around the geometric center of peaks suitable for Motif analysis in MEME-ChIP. 

In your own work you will typically run this from your own computer with MEME installed locally but today we will upload our generated FASTA file to their [web portal](http://meme-suite.org/tools/meme-chip). 

Results files from MEME-ChIP can be found [here](https://alexpmagalhaes.github.io/ChIPseq_course/resources/meme.zip)

---
### Parsing back FIMO results

We can retrieve the locations of Myc motifs identified in MEME-ChIP from the FIMO output.
FIMO reports Myc motif locations as a GFF3 file which we should be able to vizualise in IGV. 

Sadly, this GFF file's naming conventions cause only a fraction of motifs to be reported.

### FIMO to R

Fortunately we can parse our motif's GFF file into R and address this using the **import** function in  the **rtracklayer** package.

```{r, echo=TRUE,collapse=F,eval=FALSE}

motifGFF <- import("~/Downloads/fimo.gff")
```

---
### FIMO to valid GFF3

We can give the sequences some more sensible names and export the GFF to file to visualise in IGV.

```{r, echo=TRUE,collapse=F,eval=FALSE}
motifGFF$Name <- paste0(seqnames(motifGFF),":",
                        start(motifGFF),"-",end(motifGFF))
motifGFF$ID <- paste0(seqnames(motifGFF),":",
                      start(motifGFF),"-",end(motifGFF))
export.gff3(motifGFF,con="~/Downloads/fimoUpdated.gff")
```
---
### Scanning for known motifs 

We saw previously we can scan sequences using some of the Biostrings functionality **matchPattern**.
Often with ChIPseq we may know the motif we are looking for or we can use a set of known motifs from a database such as a [JASPAR](http://jaspar.genereg.net).


---
### Get motifs from JASPAR with TFBStools

We can access the model for the our motif of interest using the **TFBSTools** package and its **getMatrixByName** function.

```{r, echo=TRUE,collapse=F,eval=TRUE}
pfm <- getMatrixByName(JASPAR2020, 
                       name="MYC")
pfm
```

---
### Motif scanning with motifmatchr

With this PWM we can use the **motifmatchr** package to scan our summits for the Myc motif and return the positions of the motifs.
We will need to provide our PWM, GRanges to scan within and BSGenome object to extract sequence from. 
We also set the **out** paramter to positions for this instance.

```{r, echo=TRUE,collapse=F,eval=TRUE}

NanogMotifs <- matchMotifs(pfm,
                         macsSummits_GR,BSgenome.Mmusculus.UCSC.mm10, 
                         out = "positions")
NanogMotifs
```

---
### Exporting motif matches

We can export the Myc motif positions within peaks for use later in IGV or for metaplot vizualisation.


```{r}
export.bed(MycMotifs[[1]],con = "NanogMotifsMotifs.bed")
```

# You are done!

Notebook with all the commands can be found [here](https://github.com/alexpmagalhaes/ChIPseq_course/blob/main/docs/resources/ChipExercises.Rmd)

To go back to the home page follow this [Link](https://alexpmagalhaes.github.io/ChIPseq_course/index)


<details>
<summary>References</summary>
<br>

Furey, T.S. 2012. ChIP-seq and beyond: new and improved methodologies to detect and characterize protein-DNA interactions. [_Nat Rev Genet_](http://www.ncbi.nlm.nih.gov/pubmed/23090257) 13: 840-852

Ku, C.S., Naidoo, N., Wu, M. and Soong, R. 2011. Studying the epigenome using next generation sequencing. [_J Med Genet_](http://www.ncbi.nlm.nih.gov/pubmed/21825079) 48: 721-730

Girke, T. 2015. systemPipeR: NGS Workflow and Report Generation Environment. UC Riverside. https://github.com/tgirke/systemPipeR

Kharchenko, P.V., Tolstorukov, M.Y. and Park, P.J. 2008. Design and analysis of ChIP-seq experiments for DNA-binding proteins. [_Nature Biotechnology_](http://www.ncbi.nlm.nih.gov/pubmed/19029915) 26(12): 1351-9

Bailey, T., Krajewski, P., Ladunga, I., Lefebvre, C., Li, Q., Liu, T., Madrigal, P., Taslim, C. and Zhang, J. 2013. Practical Guidelines for the Comprehensive Analysis of ChIP-seq Data. [_PLoS Computational Biology_](http://www.ncbi.nlm.nih.gov/pubmed/24244136) 9(11): e1003326 

Lun, A.T. and Smyth, G.K. 2014. _De novo_ detection of differentially bound regions for ChIP-seq data using peaks and windows: controlling error rates correctly. [_Nucleic Acids Research_](http://www.ncbi.nlm.nih.gov/pubmed/24852250) 42(11): e95

Pepke, S., Wold, B. and Mortazavi, A. 2009. Computation for ChIP-seq and RNA-seq studies. [_Nature Methods_](http://www.ncbi.nlm.nih.gov/pubmed/19844228) 6(11s): S22-S32

Langmead, B. and Salzberg, S.L. 2012. Fast gapped-read alignment with Bowtie 2. [_Nature Methods_](http://www.ncbi.nlm.nih.gov/pubmed/22388286) 9(4): 357-9

Zhang, Y., Liu, T., Meyer, C.A., Eeckhoute, J., Johnson, D.S., Bernstein, B.E., Nussbaum, C., Myers, R.M., Brown, M., Li, W. and Liu, X.S. 2008. Model-based analysis of ChIP-Seq (MACS). [_Genome Biology_](http://www.ncbi.nlm.nih.gov/pubmed/18798982) 9(9): R137

Zhu, L., Gazin, C., Lawson, N., Pages, H., Lin, S., Lapointe, D. and Green, M. 2010. ChIPpeakAnno: a Bioconductor package to annotate ChIP-seq and ChIP-chip data. [_BMC Bioinformatics_](http://www.ncbi.nlm.nih.gov/pubmed/20459804) 11(1): 237

Yu, G., Wang, L. and He, Q. 2015. ChIPseeker: an R/Bioconductor package for ChIP peak annotation, comparison and visualization. [_Bioinformatics_](http://www.ncbi.nlm.nih.gov/pubmed/25765347) 31(14): 2382-2383

Robinson, M.D., McCarthy, D.J. and Smyth, G.K. 2010. edgeR: a Bioconductor package for differential expression analysis of digital gene expression data. [_Bioinformatics_](http://www.ncbi.nlm.nih.gov/pubmed/19910308) 26(1): 139-40

Love, M.I., Huber, W. and Anders, S. 2014. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. [_Genome Biology_](http://www.ncbi.nlm.nih.gov/pubmed/25516281) 15(12): 550

Ameur, A. 2010. BCRANK: Predicting binding site consensus from ranked DNA sequences. R package version 1.32.0 [[Link](https://www.bioconductor.org/packages/release/bioc/html/BCRANK.html)]

</details>

