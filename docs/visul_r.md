
## GO term enrichment analysis
The following will perform Gene Ontology (GO) term enrichment analysis for each annotated peak set.
```{r, eval=FALSE}
args <- systemArgs(sysma="param/macs2.param", mytargets="targets_bam_ref.txt")
args_anno <- systemArgs(sysma="param/annotate_peaks.param", mytargets="targets_macs.txt")
annofiles <- outpaths(args_anno)
gene_ids <- sapply(names(annofiles), function(x) unique(as.character(read.delim(annofiles[x])[,"geneId"])))
load("data/GO/catdb.RData")
BatchResult <- GOCluster_Report(catdb=catdb, setlist=gene_ids, method="all", id_type="gene", CLSZ=2, cutoff=0.9, gocats=c("MF", "BP", "CC"), recordSpecGO=NULL)
```

## Motif analysis
### Parse DNA sequences of peak regions from genome
Motif analysis is useful for identifying the DNA-binding motifs in the ChIP-seq peaks. When the motif of the ChIPed protein is already known, motif analysis provides validation of the success of the experiment. When the motif is not known beforehand, identifying a centrally located motif in a large fraction of the peaks is an indication of a successful experiment.

Enrichment analysis of known DNA binding motifs or *de novo* discovery of novel motifs requires the DNA sequences of the identifed peak regions. To parse the corresponding sequences from the reference genome, the `getSeq` function from the *Biostrings* package can be used.  The following example parses the sequences for each peak set and saves the results to  separate  FASTA files, one  for  each  peak  set. In  addition, the sequences in  the FASTA files are  ranked  (sorted) by increasing p-values as expected by some motif discovery tools, such as *BCRANK*.
```{r, eval=FALSE}
library(Biostrings); library(seqLogo); library(BCRANK)
args <- systemArgs(sysma="param/annotate_peaks.param", mytargets="targets_macs.txt")
rangefiles <- infile1(args)
for(i in seq(along=rangefiles)) {
    df <- read.delim(rangefiles[i], comment="#")
    peaks <- as(df, "GRanges")
    names(peaks) <- paste0(as.character(seqnames(peaks)), "_", start(peaks), "-", end(peaks))
    peaks <- peaks[order(values(peaks)$X.log10.pvalue, decreasing=TRUE)]
    pseq <- getSeq(FaFile("./data/tair10.fasta"), peaks)
    names(pseq) <- names(peaks)
    writeXStringSet(pseq, paste0(rangefiles[i], ".fasta")) 
}
```

### Motif discovery with BCRANK
The Bioconductor package *BCRANK* (Ameur A, [2010](https://www.bioconductor.org/packages/3.3/bioc/vignettes/BCRANK/inst/doc/BCRANK.pdf)) is one of the many tools available for *de novo* discovery of DNA binding motifs in the peak regions of ChIP-Seq  experiments. The given example applies this method on the first peak sample set and plots the sequence logo of the highest ranking motif.
```{r, eval=FALSE}
set.seed(0)
BCRANKout <- bcrank(paste0(rangefiles[1], ".fasta"), restarts=25, use.P1=TRUE, use.P2=TRUE)
toptable(BCRANKout)
topMotif <- toptable(BCRANKout, 1)
weightMatrix <- pwm(topMotif, normalize = FALSE)
weightMatrixNormalized <- pwm(topMotif, normalize = TRUE)
pdf("results/seqlogo.pdf")
seqLogo(weightMatrixNormalized)
dev.off()
```
![](ChIPseq_files/figure-html/seqlogo.png)

<div align="center">**Figure 6:** One of the motifs identifed by BCRANK. </div>

<div align="right">[Back to Table of Contents]()</div>

# Version Information
```{r}
sessionInfo()
```

<div align="right">[Back to Table of Contents]()</div>


# References

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