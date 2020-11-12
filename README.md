# scRNA_ReproducibilityStudy

A reproducibility study of a single cell dataset, published in this paper:

Tattikota, Sudhir Gopal, et al. "A single-cell survey of Drosophila blood." Elife 9 (2020): e54818.

## 1. Data retrieval from NCBI
10X genomics single cell RNAseq dataset from SRA [GSE146596](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE146596):

Here we focus on 4 10X genomics datasets: GSM4396377, GSM4396378, GSM4396379, GSM4396380. Each data has 4 runs (ran in 4 lanes), and the SRR id: SRR11263248 - SRR11263263

Use SRA Toolkit to download fastq files, with --split-files to separate index (I1), barcode file (R1), and sequences (R2):

```
#Download
for i in {11263248..11263263};do d="SRR"$i;echo $d;fastq-dump --split-files $d;done

#Combine. The _3.fastq file is about 30G
cat SRR11263248_1.fastq SRR11263249_1.fastq SRR11263250_1.fastq SRR11263251_1.fastq > GSM4396377_1.fastq
cat SRR11263248_2.fastq SRR11263249_2.fastq SRR11263250_2.fastq SRR11263251_2.fastq > GSM4396377_2.fastq
cat SRR11263248_3.fastq SRR11263249_3.fastq SRR11263250_3.fastq SRR11263251_3.fastq > GSM4396377_3.fastq
cat SRR11263252_1.fastq SRR11263253_1.fastq SRR11263254_1.fastq SRR11263255_1.fastq > GSM4396378_1.fastq
cat SRR11263252_2.fastq SRR11263253_2.fastq SRR11263254_2.fastq SRR11263255_2.fastq > GSM4396378_2.fastq
cat SRR11263252_3.fastq SRR11263253_3.fastq SRR11263254_3.fastq SRR11263255_3.fastq > GSM4396378_3.fastq
cat SRR11263256_1.fastq SRR11263257_1.fastq SRR11263258_1.fastq SRR11263259_1.fastq > GSM4396379_1.fastq
cat SRR11263256_2.fastq SRR11263257_2.fastq SRR11263258_2.fastq SRR11263259_2.fastq > GSM4396379_2.fastq
cat SRR11263256_3.fastq SRR11263257_3.fastq SRR11263258_3.fastq SRR11263259_3.fastq > GSM4396379_3.fastq
cat SRR11263260_1.fastq SRR11263261_1.fastq SRR11263262_1.fastq SRR11263263_1.fastq > GSM4396380_1.fastq
cat SRR11263260_2.fastq SRR11263261_2.fastq SRR11263262_2.fastq SRR11263263_2.fastq > GSM4396380_2.fastq
cat SRR11263260_3.fastq SRR11263261_3.fastq SRR11263262_3.fastq SRR11263263_3.fastq > GSM4396380_3.fastq
```

Other platforms such as Dropseq, InDrops may not be compatible with Cell Ranger pipeline:

https://www.reddit.com/r/bioinformatics/comments/8ht6kc/converting_dropseq_and_fluidigm_fastqs_into/

## 2, Prepare genome and annotation files for Cell Ranger

As claimed in the original paper, they used fruit fly genome r6.27:

Download the same genome as the eLife paper used: FlyBase r6.27:
ftp://ftp.flybase.net/genomes/Drosophila_melanogaster/dmel_r6.27_FB2019_02/fasta/

Then download gtf annotation file from:
ftp://ftp.flybase.net/genomes/Drosophila_melanogaster/dmel_r6.27_FB2019_02/gtf/



## 3, Build STAR reference for genome mapping:

## 4, Run Cell Ranger count:

