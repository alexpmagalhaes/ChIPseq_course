# Sequence QC and Trimming

## 2.1 Sequence QC


Read quality is the first step in all the analyses of sequenced reads. 

Fastq files can be downloaded from Sequence Read Archive (SRA).

In the folder `/project/pcpool_data/molmed/fastq/ChIP-seq/` you can find fastq files from the Chen et al paper.


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

If you need to download a different dataset you can download via `wget`

Exemple bellow is for nanog and GFP samples you will use.

```
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
```
To start quality trimming we will use `fastqc`

First we need to go to your project folder and make a folder for results and activate conda environment:

```
mkdir results
cd results
conda activate /project/pcpool_data/molmed/molmod
```

Second we will make qc report with `fastqc`

```
mkdir qc_report

fastqc --outdir ./qc_report --threads 10 /project/pcpool_data/molmed/fastq/ChIP-seq/GFP/*
fastqc --outdir ./qc_report --threads 10 /project/pcpool_data/molmed/fastq/ChIP-seq/nanog/*
```
Please check the html reports, we will discuss the results.

Before trimming we should merge the runs into a single fastq file.
Due to the nature of sequencing run and age of the data sets this library is fragmented.

```
cat /project/pcpool_data/molmed/fastq/ChIP-seq/nanog/*.fastq.gz > /project/pcpool_data/molmed/fastq/ChIP-seq/nanog/Nanog_Mouse_ES_merged.fastq.gz
cat /project/pcpool_data/molmed/fastq/ChIP-seq/GFP/*.fastq.gz > /project/pcpool_data/molmed/fastq/ChIP-seq/GFP/GFP_Mouse_ES_merged.fastq.gz
```

If you want you can run `fastqc` again on the merged fastq files


## 2.2 Trimming

There are many programs to do QC, and many specific tools for each one. For now we are going to focus on Cutadapt  

```
cutadapt 
    -q 30 #defines the quality
    --minimum-length 25 #minimum length of the read
    --cores 10 #number of cores
    --adapter $adaptors #adapter sequence
    --output $fq1_trimmed #output file (usualy is good to attach a reference to the name like "_trimmed")
    $FQ #input file

```

### It is important to note that the libraries covered in this tutorial are 25bp and thus no not contain any adapters so we will skip this step.

To continue with the tutorial please go to [Mapping and post processing](https://alexpmagalhaes.github.io/ChIPseq_course/mapping.md)

To go back to the home page follow this [Link](https://alexpmagalhaes.github.io/ChIPseq_course/index.md)