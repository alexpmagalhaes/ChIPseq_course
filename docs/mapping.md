# Mapping and post processing

## 3.1 Mapping aka read alignment

The next step is to align the reads to genome assembly.
This is done using Bowtie2 tool. 
The resulting .sam files are next transformed to .bam files and filtered for best aligned reads using samtools. 
PCR duplicates are removed.
The resulting BAM files and your primary output and can be then used for the next steps.

Note that bowtie2 does not support compressed files thus we need to decompress.

```
cd /project/pcpool_data/molmed/fastq/ChIP-seq/GFP/
gunzip -k -d GFP_Mouse_ES_merged.fastq.gz

cd /project/pcpool_data/molmed/fastq/ChIP-seq/nanog/
gunzip -k -d Nanog_Mouse_ES_merged.fastq.gz
```

