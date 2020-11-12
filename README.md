# scRNA_ReproducibilityStudy

A reproducibility study of a single cell dataset, published in this paper:

Tattikota, Sudhir Gopal, et al. "A single-cell survey of Drosophila blood." *Elife* 9 (2020): e54818.

## 1. Data retrieval from NCBI
10X genomics single cell RNAseq dataset from SRA [GSE146596](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE146596):

Here we focus on 4 10X genomics datasets: GSM4396377, GSM4396378, GSM4396379, GSM4396380. Each data has 4 runs (ran in 4 lanes), and the SRR id: SRR11263248 - SRR11263263

Use [SRA Toolkit](https://ncbi.github.io/sra-tools/) to download fastq files, with **--split-files** to separate index (I1), barcode file (R1), and sequences (R2): _1, _2, _3:

The R1 file contains 16bp barcode + 10bp UMI sequences. The R2 file contains single-end RNAseq read for each barcode+UNI.

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

Other platforms such as Drop-seq, inDrops may not be compatible with Cell Ranger pipeline:

https://www.reddit.com/r/bioinformatics/comments/8ht6kc/converting_dropseq_and_fluidigm_fastqs_into/

## 2, Prepare Drosophila genome and annotation files for Cell Ranger

As claimed in the original paper, they used fruit fly genome r6.27:

Download the same genome as the eLife paper used: FlyBase r6.27: [dmel-all-chromosome-r6.27.fasta](ftp://ftp.flybase.net/genomes/Drosophila_melanogaster/dmel_r6.27_FB2019_02/fasta/)

Then download gtf annotation file: [dmel-all-r6.27.gtf](ftp://ftp.flybase.net/genomes/Drosophila_melanogaster/dmel_r6.27_FB2019_02/gtf/)

After comparing the GTF I downloaded with the published single cell counts (gene names), I realized they removed mirRNA genes and snoRNA genes.

```
grep -v snoRNA dmel-all-r6.27.gtf|grep -v mir- > dmel-all-r6.27.filtered.gtf.t

```

Optional: An optional step to use `cellranger mkgtf` to keep the genes (usually by a group tag, such as 'protein_coding') we need, based on the features in column 9:
```
#For example, I have a GTF like this:
X  FlyBase  mRNA  19961689  19968479  .  +  .  gene_id "FBgn0031081"; gene_symbol "Nep3"; transcript_id "FBtr0070000"; transcript_symbol "Nep3-RA"; group "protein_coding"

#If I would like to keep GTF lines with gene_symbol == "Nep3":
cellranger mkgtf dmel-all-r6.27.filtered.gtf dmel-all-r6.27.filtered.Test.gtf --attribute=gene_symbol:Nep3
#If I would like to add group feature:
cellranger mkgtf dmel-all-r6.27.filtered.gtf dmel-all-r6.27.filtered.Test.gtf --attribute=gene_symbol:Nep3 --attribute=group:protein_coding

#I got an error:
Invalid strand in GTF line 116680: 3R FlyBase mRNA 21360390 21377399...gene_id "FBgn0002781"; gene_symbol "mod(mdg4)"; transcript_id "FBtr0084079"; transcript_symbol "mod(mdg4)-RT"; group "protein_coding"

#This was because mod(mdg4) has transcripts from both strands, so column 7 in GTF doesn't have strad information. Remove it:
awk '$12!="\"mod\(mdg4\)\"\;"'  dmel-all-r6.27.filtered.gtf.t > dmel-all-r6.27.filtered.gtf

#Final files:
-rw-r-----+ 1 ysun43 user 145942246 Nov  7 16:43 dmel-all-chromosome-r6.27.fasta
-rw-rw----+ 1 ysun43 user  77446320 Nov 11 20:45 dmel-all-r6.27.filtered.gtf
```

## 3, Build STAR reference (index) for genome mapping:

`Cellranger mkref` is a wrapper program to build SRAT index from genome fasta and GTF annotation.

For mouse and human, 10X genomices provides pre-built [genome files](https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/2.0)

```
#Check manual
cellranger mkref --help

#Run the following command in the path containing both .fasta and .gtf files:
#--genome: speficy the folder name of the genome index
#--fasta: speficy the genome file
#--genes: speficy the filtered gtf file
#--nthreads: number of CPUs
#--memgb: Maximum memory (GB) used when aligning reads with STAR. Default 16
#--ref-version=<str> Optional reference version string to include with reference.

cellranger mkref --genome=dm6 --fasta=dmel-all-chromosome-r6.27.fasta --genes=dmel-all-r6.27.filtered.gtf --nthreads 20 --memgb=96 --ref-version=r6.27
```

For fly genome, it will take ~5 mins to finish.
```
#An overview of the dm6 folder:
dm6/
    ├── mm10/
       ├── genome.fa
       └── genome.fa.fai
    ├── genes/
       └── genes.gtf
    ├── reference.json
    ├── genes/
       └── genes.pickle
    └── star/
       ...
```

## 4, Run Cell Ranger count:
Since we downloaded data from NCBI, rather than generating FASTQ files from bcl2fastq (see more details [here](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/2.0/using/mkfastq)), all data need to be renamed to meet the software criteria.

```
#If you don't follow the naming criteria, Cell Ranger will generate the following error when running cellranger count:

error: No input FASTQs were found for the requested parameters.
If your files came from bcl2fastq or mkfastq:
 - Make sure you are specifying the correct --sample(s), i.e. matching the sample sheet
 - Make sure your files follow the correct naming convention, e.g. SampleName_S1_L001_R1_001.fastq.gz (and the R2 version)
 - Make sure your --fastqs points to the correct location.
```

Rename downloaded data:
```
mv GSM4396377_1.fastq GSM4396377_S1_L001_I1_001.fastq
mv GSM4396377_2.fastq GSM4396377_S1_L001_R1_001.fastq
mv GSM4396377_3.fastq GSM4396377_S1_L001_R2_001.fastq
mv GSM4396378_1.fastq GSM4396378_S1_L001_I1_001.fastq
mv GSM4396378_2.fastq GSM4396378_S1_L001_R1_001.fastq
mv GSM4396378_3.fastq GSM4396378_S1_L001_R2_001.fastq
mv GSM4396379_1.fastq GSM4396379_S1_L001_I1_001.fastq
mv GSM4396379_2.fastq GSM4396379_S1_L001_R1_001.fastq
mv GSM4396379_3.fastq GSM4396379_S1_L001_R2_001.fastq
mv GSM4396380_1.fastq GSM4396380_S1_L001_I1_001.fastq
mv GSM4396380_2.fastq GSM4396380_S1_L001_R1_001.fastq
mv GSM4396380_3.fastq GSM4396380_S1_L001_R2_001.fastq

#Compress all files:
for i in *fastq;do gzip $i;done
```

Check the directory hierarchy:

```
scRNA/
    ├── Annotation
        ├── dm6/
            └── ...
        ├── dmel-all-chromosome-r6.27.fasta
        ├── dmel-all-r6.27.filtered.gtf
        └── Log.out
    ├── Data/
        ├── GSM4396377_S1_L001_I1_001.fastq.gz
        ├── GSM4396377_S1_L001_R1_001.fastq.gz
        ├── GSM4396377_S1_L001_R2_001.fastq.gz
        ├── GSM4396378_S1_L001_I1_001.fastq.gz
        ├── GSM4396378_S1_L001_R1_001.fastq.gz
        ├── GSM4396378_S1_L001_R2_001.fastq.gz
        ├── GSM4396379_S1_L001_I1_001.fastq.gz
        ├── GSM4396379_S1_L001_R1_001.fastq.gz
        ├── GSM4396379_S1_L001_R2_001.fastq.gz
        ├── GSM4396380_S1_L001_I1_001.fastq.gz
        ├── GSM4396380_S1_L001_R1_001.fastq.gz
        └── GSM4396380_S1_L001_R2_001.fastq.gz
    ├── Counting/
    ├── Aggr/
```

Go to the Counting folder, and submit cellranger count scripts:

```
#Check manual
cellranger count --help

#Run cellranger count
#--id: output folder name
#--transcriptome: index path of the organism, built by cellranger mkref
#--fastqs: path to all fastq files
#--sample: same as the prefix of all fastq files
#--expect-cells: optional parameter, default 3000
#--localcores: CPU
#--localmem: memory in GB
#--nosecondary: only calculate the counts. Downstream analysis will be done in R.

cellranger count --id GSM4396377 --transcriptome=../Annotation/dm6 --fastqs=../Data --sample=GSM4396377 --localcores=20 --localmem=80 --nosecondary
cellranger count --id GSM4396378 --transcriptome=../Annotation/dm6 --fastqs=../Data --sample=GSM4396378 --localcores=20 --localmem=80 --nosecondary
cellranger count --id GSM4396379 --transcriptome=../Annotation/dm6 --fastqs=../Data --sample=GSM4396379 --localcores=20 --localmem=80 --nosecondary
cellranger count --id GSM4396380 --transcriptome=../Annotation/dm6 --fastqs=../Data --sample=GSM4396380 --localcores=20 --localmem=80 --nosecondary
```

Check output files:
```
GSM4396377/
    └── outs/
       ├── web_summary.html                                     #Summary file
       ├── metrics_summary.csv                                  #Data metrics
       ├── possorted_genome_bam.bam, bai                        #Genome mapping file
       ├── filtered_gene_bc_matrices/
           └── barcodes.tsv.gz、features.tsv.gz、matrix.mtx.gz  #Results for downstream analysis in R (e.g. Seurat、Scater、Monocle)
       ├── filtered_feature_bc_matrix.h5                        #Filtered barcodes in h5 format
       ├── raw_feature_bc_matrix                                #Raw barcode info
       ├── raw_feature_bc_matrix.h5                             #Raw barcode info in HDF5 format
       ├── analysis/
           ...                                                  #If you choose to run analysis, this path contains results
       ├── molecule_info.h5                                     #Input file for cellranger aggr program
       └── cloupe.cloupe                                        #For 10X visualization tool Loupe Cell Browser
```

## 5, Run Cell Ranger aggr to merge multiple datasets

Go to the Aggr folder, and prepare the .csv file storing data information:


