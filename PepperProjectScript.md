# Scripts of the Pepper Resequencing Project

## Author: Jiantao Zhao (Email: jz426@cornell.edu; jiantaozhao426@gmail.com)

### [Google scholar link](https://scholar.google.com/citations?user=ZCdj_VcAAAAJ&hl=en&oi=ao)

## Supervisor: Zhangjun Fei (zf25@cornell.edu)

### [Lab website](http://bioinfo.bti.cornell.edu/lab/homepage/member.shtml)

### [Google scholar link](https://scholar.google.com/citations?user=QPD9ztMAAAAJ&hl=en&oi=ao)

### Date: 12/5/2020

**Note**：This file contains the majority of the scripts that I used in the pepper resequencing project, which included detailed pipelines from raw data processing to intense population genetic and genomic analyses. Some of the codes might be contributed from our lab members, especially Honghe Sun. Please be aware that the majority of the codes should work fine from my personal side and some of the codes might not be updated accordingly in different levels of analyses. For instance, in our first round of sumbitting the manuscript to Nature Genetics, the SNPs were called both in the single reference genome approach and the graph-pangenome approach. If someone wants to testify the scripts for their own analyses and meet specific issues, either check whether the software is successfully installed by checking the example data. Alternatively, you are always welcomed to send me directly via email or wechat (ID:JiantaoZhao). Please also be noted that I am not responsible for any issues that you may meet when testing the pipeline, as I have not checked the scripts in details for each line. I am more than happy for a direct discussions.

# Table of Contents

- [Scripts of the Pepper Project](#scripts-of-the-pepper-project)
  - [Author: Jiantao Zhao (Email: jz426@cornell.edu; jiantaozhao426@gmail.com)](#author-jiantao-zhao-email-jz426cornelledu-jiantaozhao426gmailcom)
    - [Google scholar link](#google-scholar-link)
  - [Supervisor: Zhangjun Fei (zf25@cornell.edu)](#supervisor-zhangjun-fei-zf25cornelledu)
    - [Lab website](#lab-website)
    - [Google scholar link](#google-scholar-link-1)
    - [Date: 12/5/2020](#date-1252020)
- [Table of Contents](#table-of-contents)
  - [Reference genome evaluation](#reference-genome-evaluation)
    - [Remove contamination](#remove-contamination)
    - [BUSCO evaluation](#busco-evaluation)
    - [Qualiy evaluation using Merqury](#qualiy-evaluation-using-merqury)
      - [Zhangshugang](#zhangshugang)
      - [Calculate the genome size using R](#calculate-the-genome-size-using-r)
      - [Grif\_1614 Merqury](#grif_1614-merqury)
      - [PI632928 Merqury](#pi632928-merqury)
    - [Genome size estimation using Jellyfish](#genome-size-estimation-using-jellyfish)
      - [Zhangshugang](#zhangshugang-1)
      - [Grif\_1614](#grif_1614)
      - [Calculate the genome size using R](#calculate-the-genome-size-using-r-1)
      - [PI632928](#pi632928)
    - [LTR Assembly Index (LAI)](#ltr-assembly-index-lai)
    - [NCBI sumbission](#ncbi-sumbission)
      - [Remove contamination for PI\_632928](#remove-contamination-for-pi_632928)
      - [Remove contamination for Grif\_1614](#remove-contamination-for-grif_1614)
  - [Reference Genome analysis](#reference-genome-analysis)
    - [Circus plot of S8/Zhangshugang](#circus-plot-of-s8zhangshugang)
    - [Genome synteny using mummer](#genome-synteny-using-mummer)
    - [Genome synteny using MCScanX](#genome-synteny-using-mcscanx)
    - [Map the Chiltepin reads to the Chiltepin genome to check the read count](#map-the-chiltepin-reads-to-the-chiltepin-genome-to-check-the-read-count)
    - [Align sequence to S8](#align-sequence-to-s8)
  - [vg pangenome](#vg-pangenome)
    - [Trimmomatic](#trimmomatic)
    - [SV estimation](#sv-estimation)
    - [vg mapping](#vg-mapping)
    - [vg mapping on cornell server](#vg-mapping-on-cornell-server)
    - [SNP calling using sentieon](#snp-calling-using-sentieon)
    - [SNP calling using sentieon for tomato as outgroup](#snp-calling-using-sentieon-for-tomato-as-outgroup)
    - [SNP calling using sentieon on cornell server](#snp-calling-using-sentieon-on-cornell-server)
    - [Join the vcf file](#join-the-vcf-file)
    - [Extract the high quality SNPs and Indels](#extract-the-high-quality-snps-and-indels)
    - [SNP and Indel snpeff](#snp-and-indel-snpeff)
      - [SNP missing rate cutoff](#snp-missing-rate-cutoff)
      - [SNP check](#snp-check)
      - [SNP check to old S8 references](#snp-check-to-old-s8-references)
      - [SNP check to vg reference for bwa mapping](#snp-check-to-vg-reference-for-bwa-mapping)
      - [SNP check to vg reference for vg mapping](#snp-check-to-vg-reference-for-vg-mapping)
    - [SNP LD prune](#snp-ld-prune)
    - [Extract 4DTV SNPs and IQ-tree](#extract-4dtv-snps-and-iq-tree)
    - [Newhyvbrid](#newhyvbrid)
    - [HyDe](#hyde)
      - [CHA and BAC Miss0.1](#cha-and-bac-miss01)
      - [Annuum Miss0.1](#annuum-miss01)
    - [Population structure](#population-structure)
    - [Fst](#fst)
      - [Fst permutation for BAB.BAP](#fst-permutation-for-babbap)
      - [Fst permutation for PUB with the other groups](#fst-permutation-for-pub-with-the-other-groups)
      - [Top 1% Fst](#top-1-fst)
    - [LD estimation](#ld-estimation)
    - [LD for Annuum clade](#ld-for-annuum-clade)
    - [Pi](#pi)
    - [Pi with only window of 100000](#pi-with-only-window-of-100000)
      - [Permutation test using 4DTV SNPs](#permutation-test-using-4dtv-snps)
      - [T-test in R](#t-test-in-r)
    - [Effective population size](#effective-population-size)
      - [ANN](#ann)
      - [GLA](#gla)
      - [BAB](#bab)
      - [BAP](#bap)
      - [CHA](#cha)
      - [CHN](#chn)
      - [FRU](#fru)
      - [PUB](#pub)
      - [Final combined plot](#final-combined-plot)
    - [GWAS](#gwas)
      - [GLA.ANN](#glaann)
        - [Shape variations](#shape-variations)
      - [BAB.BAP](#babbap)
        - [Shape variations](#shape-variations-1)
      - [Seed size](#seed-size)
        - [Extract the phenotype](#extract-the-phenotype)
        - [GLA.ANN](#glaann-1)
          - [Manhattan plot](#manhattan-plot)
        - [FRU](#fru-1)
          - [Manhattan plot](#manhattan-plot-1)
        - [CHN](#chn-1)
          - [Manhattan plot](#manhattan-plot-2)
        - [Annuum](#annuum)
          - [Manhattan plot](#manhattan-plot-3)
        - [BAB.BAP](#babbap-1)
          - [Manhattan plot](#manhattan-plot-4)
      - [CHN fruit shape](#chn-fruit-shape)
          - [Manhattan plot](#manhattan-plot-5)
    - [SweeD](#sweed)
      - [Installation](#installation)
      - [BAP with 10k grid](#bap-with-10k-grid)
      - [BAP with 100k size](#bap-with-100k-size)
      - [ANN with 10k grid](#ann-with-10k-grid)
      - [ANN with 100k grid](#ann-with-100k-grid)
      - [Overlaped regions](#overlaped-regions)
      - [Comparison with Cao et al. 2022](#comparison-with-cao-et-al-2022)
    - [Introgression](#introgression)
    - [Introgression with LD prune](#introgression-with-ld-prune)
  - [Mapping quality evaluation](#mapping-quality-evaluation)
    - [Mapping quality of the 1314 collection](#mapping-quality-of-the-1314-collection)
    - [Mapping quality of the 500 core collection](#mapping-quality-of-the-500-core-collection)
    - [Genotype check of C. baccatum](#genotype-check-of-c-baccatum)
  - [SNP calling and evaluation](#snp-calling-and-evaluation)
    - [SNPs calling and quality control of the 1300 pepper collection](#snps-calling-and-quality-control-of-the-1300-pepper-collection)
    - [SNP calling quality check](#snp-calling-quality-check)
    - [SNP count for each group](#snp-count-for-each-group)
    - [SNP prunning](#snp-prunning)
    - [SNP effects and distribution of different subgroups.](#snp-effects-and-distribution-of-different-subgroups)
      - [Upset plot in R](#upset-plot-in-r)
  - [4DTV SNPs](#4dtv-snps)
  - [IQ-tree](#iq-tree)
    - [IQ-Tree for the ~1300 accessions](#iq-tree-for-the-1300-accessions)
    - [IQ-tree for 504 accessions](#iq-tree-for-504-accessions)
    - [IQ-tree for 504 accessions with prunned SNPs](#iq-tree-for-504-accessions-with-prunned-snps)
  - [PCA analysis](#pca-analysis)
  - [Population structure analysis](#population-structure-analysis)
    - [Population structure use pruned 4DTV SNPs](#population-structure-use-pruned-4dtv-snps)
  - [Genome-wide LD](#genome-wide-ld)
  - [Fst differentiation](#fst-differentiation)
    - [Fst with accessions no taxonomy reassignment](#fst-with-accessions-no-taxonomy-reassignment)
    - [Fst with accessions of pure background](#fst-with-accessions-of-pure-background)
    - [Fst permutation test](#fst-permutation-test)
      - [BAB\_CHA](#bab_cha)
  - [Nucleotide diversity](#nucleotide-diversity)
    - [Nucleotide diversity with accessions no taxonomy reassignment](#nucleotide-diversity-with-accessions-no-taxonomy-reassignment)
    - [Nucleotide diversity with accessions of pure background](#nucleotide-diversity-with-accessions-of-pure-background)
  - [Treemix analysis](#treemix-analysis)
  - [ABBA-BABA tests](#abba-baba-tests)
  - [Introgression](#introgression-1)
    - [Introgression Chr07:4,900,001-10,000,000](#introgression-chr074900001-10000000)
    - [Introgression Chr07:258,150,001-259,300,000](#introgression-chr07258150001-259300000)
    - [Introgression Chr12:250550001-255900000](#introgression-chr12250550001-255900000)
    - [PUB.BAC.FRUCHN.GLA](#pubbacfruchngla)
  - [SweeD](#sweed-1)
    - [Installation](#installation-1)
    - [BAC with 10k grid](#bac-with-10k-grid)
    - [BAC with pruned SNPs](#bac-with-pruned-snps)
    - [GLA with 10k grid](#gla-with-10k-grid)
  - [XP-CLR for C. annuum](#xp-clr-for-c-annuum)
    - [XP-CLR for C. annuum with all C. annuum var. grabriusculum](#xp-clr-for-c-annuum-with-all-c-annuum-var-grabriusculum)
    - [XP-CLR for C. bac. baccatum versus C. bac. pendulum](#xp-clr-for-c-bac-baccatum-versus-c-bac-pendulum)
    - [Enrichment analysis](#enrichment-analysis)
  - [Effective population size](#effective-population-size-1)
    - [SMC for GLA with different geographical background](#smc-for-gla-with-different-geographical-background)
  - [GWAS analysis of fruit shape](#gwas-analysis-of-fruit-shape)
    - [GWAS for _C. annuum_](#gwas-for-c-annuum)
    - [Extract the genotypes for CaOvate](#extract-the-genotypes-for-caovate)
    - [Manhattan plot for GLA.ANN](#manhattan-plot-for-glaann)
    - [Manhattan plot for BAB.BAP](#manhattan-plot-for-babbap)
  - [Transcriptome data](#transcriptome-data)
  - [Personal scripts](#personal-scripts)
    - [4DTV\_scan.pl](#4dtv_scanpl)
    - [assign\_genetic.pl](#assign_geneticpl)
    - [PI\_slide\_win.pl](#pi_slide_winpl)
    - [merge.pl](#mergepl)


## Reference genome evaluation 

### Remove contamination

```bash
# working directory
cd /data/zhaojiantao/pepper/References/S8

# Rename the chromosomal names of the reference genome
cat pepper_S8.fasta | \
    sed 's/Superscaffold1_54_332676983/Chr01/g' | \
    sed 's/Superscaffold2_27_177319215/Chr02/g' | \
    sed 's/Superscaffold3_40_289790774/Chr03/g' | \
    sed 's/Superscaffold4_26_248932513/Chr04/g' | \
    sed 's/Superscaffold5_64_254874144/Chr05/g' | \
    sed 's/Superscaffold6_82_253233553/Chr06/g' | \
    sed 's/Superscaffold7_54_266382521/Chr07/g' | \
    sed 's/Superscaffold8_30_174326481/Chr08/g' | \
    sed 's/Superscaffold9_43_278428952/Chr09/g' | \
    sed 's/Superscaffold10_31_210332287/Chr10/g' | \
    sed 's/Superscaffold11_37_275185330/Chr11/g' | \
    sed 's/Superscaffold12_60_260676116/Chr12/g' \
    > pepper_S8.fasta

# Index the reference genome
samtools faidx pepper_S8.fasta

# Remove mitochondrion contaminations
awk '{print $1}' pepper_S8.fasta.fai | grep -v -f mitochondrion.id | tr -s "\r\n" " " > clean.id

# Rename chromosome names for genome gff file: pepper_s8.gff3
zcat pepper_S8.gff3.gz | \
    sed 's/Superscaffold10/Chr10/g' | \
    sed 's/Superscaffold11/Chr11/g' | \
    sed 's/Superscaffold12/Chr12/g' | \
    sed 's/Superscaffold1/Chr01/g' | \
    sed 's/Superscaffold2/Chr02/g' | \
    sed 's/Superscaffold3/Chr03/g' | \
    sed 's/Superscaffold4/Chr04/g' | \
    sed 's/Superscaffold5/Chr05/g' | \
    sed 's/Superscaffold6/Chr06/g' | \
    sed 's/Superscaffold7/Chr07/g' | \
    sed 's/Superscaffold8/Chr08/g' | \
    sed 's/Superscaffold9/Chr09/g' \
    > Ca.S8.gff

# Remove contaminations
samtools faidx Zhangshugang.clean.fasta

# Create a bed file containing the chr, start and end+500 position and name as contamination.bed
for i in Contig00220 Contig00221 Contig00309 Contig00415 Contig00416 Contig00417 Contig00721 ; do 
    sed -i "/$i/,+1d" Neorickii.genome.fa
done

bedtools maskfasta -mc M -fi Zhangshugang.clean.fasta -bed contamination.bed -fo masked.fasta
sed 's/M//g' masked.fasta > Zhangshugang.updated.new.fasta
samtools faidx masked.fasta
samtools faidx Zhangshugang.updated.new.fasta

# Update the gff files
grep Chr01 Ca.S8.updated.gff  | awk '$4<332426690 {print $0}' > Ca.S8.updated.Chr01.before.gff
grep Chr01 Ca.S8.updated.gff  | awk '$5>332536023 {print $1,$2,$3,$4-61608,$5-61608,$6,$7,$8,$9}' | \
    sed 's/ /\t/g' > Ca.S8.updated.Chr01.after.gff
grep -v Chr01 Ca.S8.updated.gff | grep -v Chr09 | grep -v Chr10 | \
    grep -v Chr11 | grep -v Chr12 | grep -v Contig > Ca.S8.updated.Chr02to08.gff
grep Chr09 Ca.S8.updated.gff  | awk '$5<278410513 {print $0}' > Ca.S8.updated.Chr09.before.gff
grep Chr10 Ca.S8.updated.gff > Ca.S8.updated.Chr10.gff
grep Chr11 Ca.S8.updated.gff > Ca.S8.updated.Chr11.gff
grep Chr12 Ca.S8.updated.gff  | awk '$4>939741 {print $1,$2,$3,$4-939741,$5-939741,$6,$7,$8,$9}' | \
    sed 's/ /\t/g' > Ca.S8.updated.Chr12.after.gff
grep Contig Ca.S8.updated.gff > Ca.S8.updated.Contig.gff

cat Ca.S8.updated.Chr01.before.gff Ca.S8.updated.Chr01.after.gff Ca.S8.updated.Chr02to08.gff \
    Ca.S8.updated.Chr09.before.gff Ca.S8.updated.Chr10.gff Ca.S8.updated.Chr11.gff \
    Ca.S8.updated.Chr12.after.gff Ca.S8.updated.Contig.gff | \
    sed 's/gene-//g' | sed 's/rna-//g' \
    > Ca.S8.updated.new.gff

# Extract cds
for i in Ca.S8.updated.Chr01.before Ca.S8.updated.Chr01.after Ca.S8.updated.Chr02to08 \
    Ca.S8.updated.Chr09.before Ca.S8.updated.Chr10 Ca.S8.updated.Chr11 Ca.S8.updated.Chr12.after ; do
    gffread -x $i.cds.fasta -g Zhangshugang.updated.new.fasta $i.gff
done

cat Ca.S8.updated.Chr01.before.cds.fasta Ca.S8.updated.Chr01.after.cds.fasta \
    Ca.S8.updated.Chr02to08.cds.fasta Ca.S8.updated.Chr09.before.cds.fasta \
    Ca.S8.updated.Chr10.cds.fasta Ca.S8.updated.Chr11.cds.fasta \
    Ca.S8.updated.Chr12.after.cds.fasta Ca.S8.updated.Contig.cds.fasta \
    > Zhangshugang.updated.new.cds.fasta

# Alternatively
# cds
gffread -x Zhangshugang.updated.new.cds.fasta -g Zhangshugang.updated.new.fasta Ca.S8.updated.new.gff

# protein
gffread -y Zhangshugang.updated.new.protein.fasta -g Zhangshugang.updated.new.fasta Ca.S8.updated.new.gff

# Upload to NCBI

# Update the vcf file
# For Chr01
zcat Chr01.clean.vcf.gz | grep "#" > header
zcat Chr01.clean.vcf.gz | grep -v "#" | \
    awk '$2<332426690 {print $0}' > Chr01.clean.vcf-1
zcat Chr01.clean.vcf.gz | grep -v "#" | \
    awk '$2>332536023 {print $1,$2-61608,$0}' | \
    sed 's/ /\t/g' | cut -f1,2,5- \
    > Chr01.clean.vcf-2
cat header Chr01.clean.vcf-1 Chr01.clean.vcf-2 | \
    gzip > Chr01.clean.new.vcf.gz

# Chr09
zcat Chr09.clean.vcf.gz | grep "#" > header
zcat Chr09.clean.vcf.gz | grep -v "#" | awk '$2<278410513 {print $0}' | \
    cat header - | gzip > Chr09.clean.new.vcf.gz

# Chr12
zcat Chr12.clean.vcf.gz | grep "#" > header
zcat Chr12.clean.vcf.gz | grep -v "#" | \
    awk '$2>939741 {print $1,$2-939741,$0}' | \
    sed 's/ /\t/g' | cut -f1,2,5- | cat header - | gzip > Chr12.clean.new.vcf.gz

# Cross-check of the genotype
## Chr01
# Extract the SNP position infos
grep -v "#" Chr01.clean.vcf | awk '{print $1,$2-1,$2}' | sed 's/ /\t/g' > Chr01.bed
# Extract the geno from references
bedtools getfasta -fi \
    /data/zhaojiantao/pepper/References/S8/Zhangshugang.fasta \
    -bed Chr01.bed -name -tab \
    > Chr01.bed.out

# Extract the reference genotypes from vcf
bcftools query -f '%CHROM %POS %REF\n' Chr01.clean.vcf.gz > Chr01.vcf.bed

# Check the consistency of genotypes
paste Chr01.bed.out  Chr01.vcf.bed | awk '{if ($2==$5) print "TRUE"; else print "FALSE"}' | grep FALSE | wc 
# Number of inconsistency: 0

## Chr12 
# Extract the SNP position infos
grep -v "#" Chr12.clean.vcf | awk '{print $1,$2-1,$2}' | sed 's/ /\t/g' > Chr12.bed
# Extract the geno from references
bedtools getfasta -fi \
    /data/zhaojiantao/pepper/References/S8/Zhangshugang.fasta \
    -bed Chr12.bed -name -tab \
    > Chr12.bed.out

# Extract the reference genotypes from vcf
bcftools query -f '%CHROM %POS %REF\n' Chr12.clean.vcf.gz > Chr12.vcf.bed

# Check the consistency of genotypes
paste Chr12.bed.out  Chr12.vcf.bed | awk '{if ($2==$5) print "TRUE"; else print "FALSE"}' | grep FALSE | wc 
# Number of inconsistency: 0

```

### BUSCO evaluation

```bash
# Working directory
cd /data/zhaojiantao/pepper/References/S8

# Genome completeness
busco -i /data/zhaojiantao/pepper/References/S8/Zhangshugang.updated.new.fasta \
    -l /data/zhaojiantao/database/embryophyta_odb10 \
    -c 20 -m genome -f --out Zhangshugang.embryophyta.genome

# Gene completeness
busco -i /data/zhaojiantao/pepper/References/S8/Zhangshugang.updated.new.protein.fasta \
    -l /data/zhaojiantao/database/embryophyta_odb10 \
    -c 20 -m proteins -f --out Zhangshugang.embryophyta.gene

```

### Qualiy evaluation using Merqury

#### Zhangshugang

```bash
# [Merqury](https://github.com/marbl/merqury)
cd /data/zhaojiantao/pepper/References/S8/Merqury

# Get the right k size
best_k.sh 3020000000
# Best K-mer size = 20

# Build meryl dbs
meryl k=20 count output S8_R1.meryl S8_R1.fastq.gz
meryl k=20 count output S8_R2.meryl S8_R2.fastq.gz

# Merge the files
meryl union-sum output S8.genome.meryl S8_R*.meryl

# Calculate QV score
merqury.sh S8.genome.meryl ../pepper_S8.fasta S8.QV

```

#### Calculate the genome size using R

```r
histo <- read.table("S8.reads.histo") 
dim(histo)
# Estimate genome size
sum(as.numeric(histo[2:100896,1]*histo[2:100896,2]))/45

```

#### Grif_1614 Merqury

```bash
# [Merqury](https://github.com/marbl/merqury)
cd /data/zhaojiantao/pepper/pepper_raw_data/Grif_1614.HiFI/Merqury

# Get the right k size
best_k.sh 3927487386
# Best K-mer size = 21

# Build meryl dbs
meryl k=21 count output threads=80 Grif_1614.genome.meryl ../Grif1614.fastq.gz

# Calculate QV score
merqury.sh Grif_1614.genome.meryl /data/zhaojiantao/pepper/References/Grif_1614/Grif_1614.fa Grif_1614

# Check the result summary
Grif_1614.QV.qv
Grif_1614.QV.completeness.stats

```

#### PI632928 Merqury

```bash
# [Merqury](https://github.com/marbl/merqury)
cd /data/zhaojiantao/pepper/pepper_raw_data/PI632928_Illumina/Merqury

# Build meryl dbs
meryl k=21 count output B554.R1.meryl B554.R1.fastq.gz
meryl k=21 count output B554.R2.meryl B554.R2.fastq.gz

# Merge the files
meryl union-sum output B554.genome.meryl B554.R*.meryl

# Calculate QV score
merqury.sh B554.genome.meryl /data/zhaojiantao/pepper/References/PI_632928/PI_632928.fa PI_632928.QV

```

### Genome size estimation using Jellyfish

#### Zhangshugang

```bash
# k-mer evaluation
# Obtain the k-mer frequencies of your raw data (fastq) using JellyFish
jellyfish count -C -m 21 -s 1000000000 -t 30 <(zcat S8_R1.fastq.gz) <(zcat S8_R2.fastq.gz) -o S8.reads.jf

# Export the kmer count histogram
jellyfish histo -t 10 -h 1000000 S8.reads.jf > S8.reads.histo

# Run Genomescope
genomescope.R S8.reads.histo 21 150 S8.Kmer

```

#### Grif_1614

```bash
# Working directory
cd /data/zhaojiantao/pepper/pepper_raw_data/Grif_1614.HiFI/Jellyfish

# k-mer evaluation
# Obtain the k-mer frequencies of your raw data (fastq) using JellyFish
jellyfish bc -m 21 -s 100000000000 -t 80 -C -o Grif1614.bc ../Grif1614.fastq
# Jellyfish count
jellyfish count -m 21 -s 4000000000 -t 80 -C --bc Grif1614.bc ../Grif1614.fastq -o Grif_1614.reads.jf

# Export the kmer count histogram
jellyfish histo -t 100 -h 1000000 Grif_1614.reads.jf > Grif_1614.reads.histo

# Run Genomescope
genomescope.R Grif_1614.reads.histo 21 14210 Grif_1614.GenomeScope

# If you met an error showing: Failed to converge, which could possible due to some of the reads count x coverage exceeded the limited of maximum number. You can simply select a portion (eg. 3/4) of your data and run it again. 

```

#### Calculate the genome size using R

```r
histo <- read.table("Grif_1614.reads.histo") 
dim(histo)
# Estimate genome size
sum(as.numeric(histo[2:100896,1]*histo[2:100896,2]))/20

```

#### PI632928

```bash
# Working directory
cd /data/zhaojiantao/pepper/pepper_raw_data/PI632928_Illumina/Jellyfish

jellyfish count -C -m 21 -s 3000000000 -t 80 PI632928.fastq -o PI632928.reads.jf

# Export the kmer count histogram
jellyfish histo -t 100 -h 1000000 PI632928.reads.jf > PI632928.reads.histo

# Run Genomescope
genomescope.R PI632928.reads.histo 21 150 PI632928.Genomescope

```


### [LTR Assembly Index (LAI)](https://github.com/oushujun/LTR_retriever)

```bash
# Working directory
/data/zhaojiantao/pepper/References/S8

# Install [RepeatMasker](http://www.repeatmasker.org/RepeatMasker/RepeatMasker-4.1.2-p1.tar.gz)
wget http://www.repeatmasker.org/RepeatMasker/RepeatMasker-4.1.2-p1.tar.gz
tar -xvzf RepeatMasker-4.1.2-p1.tar.gz
# Install and update RepeatMasker Libraries
cd RepeatMasker/Libraries
rm Dfam.h5

# Install LTR_finder from https://github.com/xzhub/LTR_Finder
wget https://github.com/xzhub/LTR_Finder/archive/refs/heads/master.zip
unzip master.zip
cd LTR_Finder-master/source

# Index the genome
gt suffixerator -db ../Zhangshugang.fasta -indexname ../Zhangshugang.fasta -tis -suf -lcp -des -ssp -sds -dna

gt ltrharvest \
    -index ../Zhangshugang.fasta \
    -similar 85 -vic 10 -seed 20 -seqids yes \
    -minlenltr 100 -maxlenltr 7000 -mintsd 4 -maxtsd 6 \
    -motif TGCA -motifmis 1 \
    > Zhangshugang.harvest.scn

# Use [LTR_FINDER_parallel](https://github.com/oushujun/LTR_FINDER_parallel) to accelerate the speed
LTR_FINDER_parallel \
    -seq ../Zhangshugang.fasta \
    -harvest_out \
    -size 1000000 \
    -time 300 \
    -threads 50

cat Zhangshugang.harvest.scn Zhangshugang.fasta.finder.combine.scn > Zhangshugang.rawLTR.scn

LTR_retriever \
    -genome ../Zhangshugang.fasta \
    -inharvest Zhangshugang.rawLTR.scn \
    -threads 60

mv Zhangshugang.fasta.mod.out.LAI Zhangshugang.fasta.mod.out.default.LAI

# Add -unlock
perl /data/zhaojiantao/tools/LTR_retriever-master/LAI \
    -genome Zhangshugang.fasta \
    -intact Zhangshugang.fasta.pass.list \
    -all Zhangshugang.fasta.out \
    -t 10 -q -unlock

```

### NCBI sumbission

```bash
# For Grif_1614 genome assembly
cd /data/zhaojiantao/pepper/References/NCBI.submit

# Submit via ascp for Grif_1614
ascp -i /data/zhaojiantao/pepper/References/NCBI.submit/aspera.openssh \
    -QT -l1000m -k1 -d Grif_1614 subasp@upload.ncbi.nlm.nih.gov:uploads/zf25_cornell.edu_cFD55kCx

# Submit via ascp for PI_632928
ascp -i /data/zhaojiantao/pepper/References/NCBI.submit/aspera.openssh \
    -QT -l1000m -k1 -d PI_632928 subasp@upload.ncbi.nlm.nih.gov:uploads/zf25_cornell.edu_cFD55kCx

```

#### Remove contamination for PI_632928

```bash
# Working directory
cd /data/zhaojiantao/pepper/References/PI_632928

# Replace the contig ids
seqkit replace -k chr.ids --pattern "^(\w+)" --replacement "{kv}" -w 0 PI_632928.fa > ../../NCBI.submit/PI_632928/PI_632928.fa

# Extract the contig id
grep HiC_scaffold Genome.without.contig.merge/PI_632928.fa.fai | cut -f1 > Chr00.ids
grep -v HiC_scaffold Genome.without.contig.merge/PI_632928.fa.fai | cut -f1 > Chr.ids

# Extract the Chr00 sequences
seqkit grep -n -f Chr00.ids Genome.without.contig.merge/PI_632928.fa -w 0 -o Chr00.old.fa
seqkit grep -n -f Chr.ids Genome.without.contig.merge/PI_632928.fa -w 0 -o Chr.fa

# Rename the seq ids
seqkit replace -k Genome.without.contig.merge/contig.ids --pattern "^(\w+)" --replacement "{kv}" -w 0 Chr00.old.fa 

# Remove contamination
bedtools maskfasta -mc M -fi PI_639682_genome_submit.fa -bed contamination.bed2 -fo PI_639682_genome_submit.fa-1
sed 's/M//g' PI_639682_genome_submit.fa-1 > PI_639682.withctg.updated.fa

# If there are any contaminations, replace with N
bedtools maskfasta -mc N -fi PI_632928.fa -bed contamination.bed -fo PI_632928.fa2

# Alternatively, replace it with others and remove it
bedtools maskfasta -mc M -fi PI_632928.old.fa -bed contamination.bed2 -fo PI_632928.old.fa2
sed 's/M//g' PI_632928.old.fa2 > PI_632928.updated.fa

```

#### Remove contamination for Grif_1614

```bash
# Working directory
cd /data/zhaojiantao/pepper/References/Grif_1614/Genome.without.contig.merge

# Extract the contig id
grep -v scaffold62 Grif_1614.fa.fai | cut -f1 > Chr.ids
seqkit grep -n -f Chr.ids Grif_1614.fa -w 0 -o Chr.fa

# Prepare the bed file for masking contaminations
awk '{print $1,$2-1,$3}' contamination.ctg | sed 's/ /\t/g' > contamination.ctg2

# Replace contaminations with N
bedtools maskfasta -mc N -fi Chr.fa -bed contamination.ctg2 -fo Grif_1614_genome_submit.fa

# Re-submit to NCBI
cp Grif_1614_genome_submit.fa /data/zhaojiantao/pepper/References/NCBI.submit/Grif_1614.update1

ascp -i /data/zhaojiantao/pepper/References/NCBI.submit/aspera.openssh \
    -QT -l1000m -k1 -d /data/zhaojiantao/pepper/References/NCBI.submit/Grif_1614.update1 subasp@upload.ncbi.nlm.nih.gov:uploads/zf25_cornell.edu_cFD55kCx

ascp -i /data/zhaojiantao/pepper/References/NCBI.submit/aspera.openssh \
    -QT -l1000m -k1 -d /data/zhaojiantao/pepper/References/NCBI.submit/Grif_1614.update3 subasp@upload.ncbi.nlm.nih.gov:uploads/zf25_cornell.edu_cFD55kCx


# Working directory
cd /data/zhaojiantao/pepper/References/Grif_1614

# Update the annotation
# Remove the following genes from Chr10 due to contaminations named as Contamination.genes
Cpu10g10630
Cpu10g10640
Cpu10g10650
Cpu10g10750
Cpu10g10760
Cpu10g10770
Cpu10g10780
Cpu00g04710
Cpu00g04720
Cpu00g04730

# Remove the genes from gff3
grep -vf Contamination.genes Grif_1614.old.gff3 > Grif_1614.gff3

# Clean genome file
bedtools maskfasta -mc N -fi Grif_1614.old.fa -bed Genome.without.contig.merge/contamination.ctg2 -fo Grif_1614-1.fa
bedtools maskfasta -mc M -fi Grif_1614-1.fa -bed Contamination.bed2 -fo Grif_1614-2.fa
sed 's/M//g' Grif_1614-2.fa > Grif_1614.fa

# Extract the CDS and pep sequences
gffread -x Grif_1614.CDS.fa -g Grif_1614.fa Grif_1614.gff3
gffread -y Grif_1614.pep.fa -g Grif_1614.fa Grif_1614.gff3
gffread -w Grif_1614.mRNA.fa -g Grif_1614.fa Grif_1614.gff3

```


## Reference Genome analysis
### Circus plot of S8/Zhangshugang

```bash
# Working directory
cd /data/zhaojiantao/pepper/References/S8/Circos

# Prepare genome bed file
grep Super ../pepper_S8.fasta.ann | sed 's/_/ /g' | awk '{print $2,$4}' | \
    sed 's/Superscaffold/chr/g' | sed 's/ /\t/g' \
    > pepper_S8.genome.bed

# Split the genome into 1 Mb window
bedtools makewindows -g pepper_S8.genome.bed -w 500000 | awk '{print $1,$2+1,$3}' | sed 's/ /\t/g' > S8.win500kb.bed

# Calculate the gene numbers within 1 Mb window
awk '{if ($3=="gene") print $1,$4-1,$5}' ../gff | sed 's/Chr0/chr/g' | \
    sed 's/Chr/chr/g' | grep -v chr0. | sed 's/ /\t/g' | \
    bedtools coverage -a S8.win500kb.bed -b - -counts | \
    sed 's/ /\t/g' \
    > S8.gene.Win500kb.bed.txt

# Calculate the transposon percentage within 1 Mb window
zcat ../repeat/all.gff.gff3.gz | grep -v unanchor | grep Transposon | grep -v "#" | awk '{print $1,$4-1,$5}' | \
    sed 's/Superscaffold/chr/g' | sed 's/ /\t/g' | \
    bedtools coverage -a S8.win500kb.bed -b - | \
    awk '{print $1,$2,$3,$7}' | \
    sed 's/ /\t/g' \
    > S8.TE.Win500kb.bed.txt

# Run Circos
circos -conf circos.conf -debug_group legend

# Combine the plots together
smc++ plot combined.pdf \
    ANN/smc.out/model.final.json \
    BAB/smc.out/model.final.json \
    BAP/smc.out/model.final.json \
    CHA/smc.out/model.final.json \
    CHN/smc.out/model.final.json \
    FRU/smc.out/model.final.json \
    GLA1/smc.out/model.final.json \
    GLA2/smc.out/model.final.json \
    PUB/smc.out/model.final.json
```

### Genome synteny using mummer

```bash
# Working directory
# Remove unanchorred sequences
head -24 ../References/S8/pepper_S8.fasta | sed 's/Superscaffold1_54_332676983/Chr01/g' | \
    sed 's/Superscaffold2_27_177319215/Chr02/g' | sed 's/Superscaffold3_40_289790774/Chr03/g' | \
    sed 's/Superscaffold4_26_248932513/Chr04/g' | sed 's/Superscaffold5_64_254874144/Chr05/g' | \
    sed 's/Superscaffold6_82_253233553/Chr06/g' | sed 's/Superscaffold7_54_266382521/Chr07/g' | \
    sed 's/Superscaffold8_30_174326481/Chr08/g' | sed 's/Superscaffold9_43_278428952/Chr09/g' | \
    sed 's/Superscaffold10_31_210332287/Chr10/g' | sed 's/Superscaffold11_37_275185330/Chr11/g' | \
    sed 's/Superscaffold12_60_260676116/Chr12/g' > S8.fa

head -24 ../References/CM334/v2.0/CM334.fasta | sed 's/chr/Chr/g' > CM334.fa
sed -n '1,44153421p' ../References/Zunla-1/v2.0/Zunla-1.fasta > Zunla-1.fa
sed -n '1,40892390p' ../References/Chiltepin/v2.0/Chiltepin.fasta | cut -f1 > Chiltepin.fa
head -24 ../References/PI159236/PI159236.fasta > PI159236.fa
sed -n '1,44535272p' ../References/CM334F1/CM334F1.fasta | sed 's/chr/Chr0/g' | \
    sed 's/Chr010/Chr10/g' | sed 's/Chr011/Chr11/g' | sed 's/Chr012/Chr12/g' > CM334F1.fa

# Use (Mummer)[https://github.com/mummer4/mummer]
wget https://github.com/mummer4/mummer/releases/download/v4.0.0rc1/mummer-4.0.0rc1.tar.gz
tar -xvzf mummer-4.0.0rc1.tar.gz
cd mummer-4.0.0rc1
./configure --prefix=/data/zhaojiantao/tools/mummer-4.0.0rc1/
make
make install

# Genome-wide running takes a long time to complete, split the fasta by chromosomes
for i in $(seq -w 1 12); do 
    nucmer --minmatch=1000 --mincluster=1000 --minalign=1000 -t 10 S8.Chr$i.fa Chiltepin.Chr$i.fa -p S8.Chiltepin.Chr$i.nucmer
    delta-filter -r S8.Chiltepin.Chr$i.nucmer.delta -q > S8.Chiltepin.Chr$i.nucmer.filter.delta
done

# Combine the files
head -2 S8.Chiltepin.Chr01.nucmer.filter.delta | sed 's/Chr01.//g' | \
    cat - <(cat S8.Chiltepin.Chr*.nucmer.filter.delta | grep -v zhaojiantao | grep -v NUCMER) \
    > S8.Chiltepin.filter.delta

zcat S8.Chiltepin.Chr01.nucmer.filter.delta.gz | head -2 | sed 's/Chr01.//g' | \
    cat - <(zcat S8.Chiltepin.Chr*.nucmer.filter.delta.gz | grep -v zhaojiantao | grep -v NUCMER) \
    > S8.Chiltepin.filter.delta

# Generate the plot
mummerplot --filter S8.Chiltepin.filter.delta --postscript --medium -p S8.Chiltepin.filter.delta

# Convert ps to pdf
ps2pdf S8.Chiltepin.filter.delta.ps S8.Chiltepin.filter.delta.pdf
convert -density 300 S8.Chiltepin.filter.delta.pdf S8.Chiltepin.filter.delta.png

# Alternatively, generate the plot using (DotPlotly)[https://github.com/tpoorten/dotPlotly/tree/master/example]
for i in $(seq -w 1 12); do 
    delta-filter -r S8.Chiltepin.Chr$i.nucmer.delta > S8.Chiltepin.Chr$i.nucmer.filter.delta
    show-coords -c S8.Chiltepin.Chr$i.nucmer.filter.delta > S8.Chiltepin.Chr$i.nucmer.filter.delta.coords
done

for i in $(seq -w 2 12); do  
    sed -n '6,$p' S8.Chiltepin.Chr$i.nucmer.filter.delta.coords > S8.Chiltepin.Chr$i.nucmer.filter.delta2.coords
done

cat S8.Chiltepin.Chr01.nucmer.filter.delta.coords S8.Chiltepin.Chr*.nucmer.filter.delta2.coords | \
    sed 's/S8.Chr01/S8/g' | sed 's/Chiltepin.Chr01/Chiltepin/g' | \
    sed 's/Chr0//g' | sed 's/Chr//g' > S8.Chiltepin.coords

mummerCoordsDotPlotly.R -i S8.Chiltepin.coords -o S8.Chiltepin.nucmer.plot \
    -m 1000 -q 300000 -k 12 -s -t -p 5 -r 1,2,3,4,5,6,7,8,9,10,11,12

```

### Genome synteny using [MCScanX](https://github.com/wyp1125/MCScanX)

```bash
# Working directory
cd /data/zhaojiantao/pepper/MCScanX

# Use S8 as the reference
# Extract the gene id
cat /data/zhaojiantao/pepper/References/S8/gene/pepper_S8.gff3 | grep -v unanchor | \
    sed 's/;/\t/g' | grep mRNA | awk '{print "zs"$1,$9,$4,$5}' | sed 's/ID=//g' | sed 's/Superscaffold//g' > S8.gene.gff
# Extract the sequences
seqkit grep -j 20 -f <(cut -f2 -d ' ' S8.gene.gff) /data/zhaojiantao/pepper/References/S8/gene/pepper_S8.pep.fa > S8.fa

# Create index as the reference
makeblastdb -in S8.fa -dbtype prot -out index/S8 -parse_seqids

# Create index as the reference
makeblastdb -in Chiltepin.fa -dbtype prot -out index/Chiltepin -parse_seqids

# CA59
cat /data/zhaojiantao/pepper/References/CA59/Ca_59.v1a.gff3 | grep -v scaffold | \
    sed 's/;/\t/g' | grep mRNA | awk '{print "ca"$1,$9,$4,$5}' | \
    sed 's/ID=//g' | sed 's/Ca_59Chr0//g' | sed 's/Ca_59Chr//g' > CA59.gene.gff
# Extract the sequences
seqkit grep -j 20 -f <(cut -f2 -d ' ' CA59.gene.gff) \
    /data/zhaojiantao/pepper/References/CA59/Ca_59.protein.fa > CA59.fa

# PBC81
zcat /data/zhaojiantao/pepper/References/PBC81/PBC81.chromosome.gff3.gz | \
    grep gene | awk '{print "pb"$1,$9,$4,$5}' | \
    grep Chr | sed 's/Chr0//g' | sed 's/Chr//g' | \
    sed 's/ID=//g' | sed 's/ /\t/g' > PBC81.gene.gff
# Extract the sequences
seqkit grep -j 20 -f <(cut -f2 PBC81.gene.gff) \
    /data/zhaojiantao/pepper/References/PBC81/PBC81.PEP.fa | \
    cut -f1 | sed 's/*//g' > PBC81.fa

# Chiltepin
zcat /data/zhaojiantao/pepper/References/Chiltepin/v2.0/Chiltepin.genes.gff.gz | \
    grep mRNA | awk '{print "ch"$1,$9,$4,$5}' | \
    grep Chr | grep -v Chr00 | sed 's/Chr0//g' | sed 's/Chr//g' | \
    sed 's/ID=//g' | sed 's/;//g' | sed 's/ /\t/g' \
    > Chiltepin.gene.gff
# Extract the sequences
seqkit grep -j 20 -f <(cut -f2 Chiltepin.gene.gff) \
    /data/zhaojiantao/pepper/References/Chiltepin/v2.0/Chiltepin.PEP.fasta | \
    cut -f1 | sed 's/*//g' > Chiltepin.fa

# CM334
cat /data/zhaojiantao/pepper/References/CM334/annotation/Pepper.v.1.55.chr.annotated.gff | \
    grep mRNA | sed 's/;/\t/g' | \
    awk '{print "cm"$1,$9,$4,$5}' | grep -v chr00 | \
    sed 's/Pepper.v.1.55.chr0//g' | sed 's/Pepper.v.1.55.chr//g' | \
    sed 's/ID=TC.//g' | sed 's/ /\t/g' > CM334.gene.gff
# Extract the sequences
seqkit grep -j 20 -f <(cut -f2 CM334.gene.gff) \
    /data/zhaojiantao/pepper/References/CM334/annotation/Pepper.v.1.55.PEP.fa | \
    cut -f1 | sed 's/*//g' > CM334.fa

# PI159236
zcat /data/zhaojiantao/pepper/References/PI159236/PI159236.chromosome.gff3.gz | \
    grep mRNA | sed 's/;/\t/g' | grep Chr | \
    awk '{print "pi"$1,$9,$4,$5}' | sed 's/ID=TC.//g' | \
    sed 's/piChr0/pi/g' | sed 's/piChr/pi/g' | sed 's/ /\t/g' \
    > PI159236.gene.gff
# Extract the sequences
seqkit grep -j 20 -f <(cut -f2 PI159236.gene.gff) \
    /data/zhaojiantao/pepper/References/PI159236/PI159236.PEP.fa | \
    cut -f1 | sed 's/*//g' > PI159236.fa

# Zunla-1
cat /data/zhaojiantao/pepper/References/Zunla-1/v2.0/Zunla-1.genes.gff | \
    grep mRNA | sed 's/;/\t/g' | grep -v Chr00 | \
    awk '{print "zu"$1,$9,$4,$5}' | sed 's/ID=//g' | sed 's/zuChr0/zu/g' | \
    sed 's/zuChr/zu/g' | sed 's/ /\t/g' \
    > Zunla-1.gene.gff
# Extract the sequences
seqkit grep -j 20 -f <(cut -f2 Zunla-1.gene.gff) \
    /data/zhaojiantao/pepper/References/Zunla-1/v2.0/Zunla-1.PEP.fa | \
    cut -f1 | sed 's/*//g' > Zunla-1.fa

# C. bac. pendulum PI_632928
cat /data/zhaojiantao/pepper/References/PI_632928/PI_632928.gff3 | \
    grep mRNA | sed 's/;/\t/g' | grep -v Chr00 | \
    awk '{print "cb"$1,$9,$4,$5}' | sed 's/ID=//g' | sed 's/cbChr0/cb/g' | \
    sed 's/cbChr/cb/g' | sed 's/ /\t/g' \
    > PI_632928.gene.gff
# Extract the sequences
seqkit grep -j 20 -f <(cut -f2 PI_632928.gene.gff) \
    /data/zhaojiantao/pepper/References/PI_632928/PI_632928.pep.fa | \
    cut -f1 | sed 's/*//g' > PI_632928.fa
# Remove "." from the sequences
perl fasta2rawseq.cgi PI_632928.fa
perl rawseq2fasta_protein.cgi new.fasta > AA
mv AA PI_632928.fa

# C.pubescens
cat /data/zhaojiantao/pepper/References/Grif_1614/Grif_1614.gff3 | \
    grep mRNA | sed 's/;/\t/g' | grep -v Chr00 | \
    awk '{print "cp"$1,$9,$4,$5}' | sed 's/ID=//g' | sed 's/cpChr0/cp/g' | \
    sed 's/cpChr/cp/g' | sed 's/ /\t/g' \
    > Grif_1614.gene.gff
# Extract the sequences
seqkit grep -j 20 -f <(cut -f2 Grif_1614.gene.gff) \
    /data/zhaojiantao/pepper/References/Grif_1614/Grif_1614.pep.fa | \
    cut -f1 | sed 's/*//g' > Grif_1614.fa
perl fasta2rawseq.cgi Grif_1614.fa
perl rawseq2fasta_protein.cgi new.fasta > AA
mv AA Grif_1614.fa

# XLL
cat /data/zhaojiantao/pepper/References/XLL/XLL.gff3 | \
    grep mRNA | sed 's/;/\t/g' | grep -v Chr00 | \
    awk '{print "xl"$1,$9,$4,$5}' | sed 's/ID=//g' | sed 's/xlChr0/xl/g' | \
    sed 's/xlChr/xl/g' | sed 's/ /\t/g' \
    > XLL.gene.gff
# Extract the sequences
seqkit grep -j 20 -f <(cut -f2 XLL.gene.gff) \
    /data/zhaojiantao/pepper/References/XLL/XLL.pep.fa | \
    cut -f1 | sed 's/*//g' > XLL.fa
perl fasta2rawseq.cgi XLL.fa
perl rawseq2fasta_protein.cgi new.fasta > AA
mv AA XLL.fa

# Create index as the reference
makeblastdb -in PI_632928.fa -dbtype prot -out index/PI_632928 -parse_seqids
makeblastdb -in Grif_1614.fa -dbtype prot -out index/Grif_1614 -parse_seqids

# Chiltepin versus Zunla-1

# Align protein sequences to S8
blastp -query PBC81.fa -db index/S8 -out S8_PBC81.blast -evalue 1e-10 -num_threads 60 -outfmt 6 -max_target_seqs 1
blastp -query XLL.fa -db index/S8 -out S8_XLL.blast -evalue 1e-10 -num_threads 60 -outfmt 6 -max_target_seqs 1

blastp -query S8.fa -db index/PI_632928 -out PI_632928_S8.blast -evalue 1e-10 -num_threads 60 -outfmt 6 -max_target_seqs 1
blastp -query S8.fa -db index/Grif_1614 -out Grif_1614_S8.blast -evalue 1e-10 -num_threads 60 -outfmt 6 -max_target_seqs 1
blastp -query Zhangshugang.fa -db index/PI_632928 -out PI_632928_Zhangshugang.blast -evalue 1e-10 -num_threads 60 -outfmt 6 -max_target_seqs 1
blastp -query Zhangshugang.fa -db index/Grif_1614 -out Grif_1614_Zhangshugang.blast -evalue 1e-10 -num_threads 60 -outfmt 6 -max_target_seqs 1

blastp -query CA59.fa -db index/S8 -out S8_CA59.blast -evalue 1e-10 -num_threads 60 -outfmt 6 -max_target_seqs 1

blastp -query Chiltepin.fa -db index/S8 -out S8_Chiltepin.blast -evalue 1e-10 -num_threads 60 -outfmt 6 -max_target_seqs 1
blastp -query CM334.fa -db index/S8 -out S8_CM334.blast -evalue 1e-10 -num_threads 60 -outfmt 6 -max_target_seqs 1
blastp -query PI159236.fa -db index/S8 -out S8_PI159236.blast -evalue 1e-10 -num_threads 60 -outfmt 6 -max_target_seqs 1
blastp -query Zunla-1.fa -db index/S8 -out S8_Zunla-1.blast -evalue 1e-10 -num_threads 60 -outfmt 6 -max_target_seqs 1

blastp -query Chiltepin.fa -db index/Zunla-1 -out Zunla-1_Chiltepin.blast -evalue 1e-10 -num_threads 60 -outfmt 6 -max_target_seqs 5
blastp -query Chiltepin.fa -db index/CM334 -out CM334_Chiltepin.blast -evalue 1e-10 -num_threads 60 -outfmt 6 -max_target_seqs 5
blastp -query PBC81.fa -db index/CM334 -out CM334_PBC81.blast -evalue 1e-10 -num_threads 30 -outfmt 6 -max_target_seqs 5
blastp -query Chiltepin.fa -db index/PI159236 -out PI159236_Chiltepin.blast -evalue 1e-10 -num_threads 60 -outfmt 6 -max_target_seqs 5

blastp -query Chiltepin.fa -db index/PBC81 -out PBC81_Chiltepin.blast -evalue 1e-10 -num_threads 60 -outfmt 6 -max_target_seqs 5
blastp -query CA59.fa -db index/PBC81 -out PBC81_CA59.blast -evalue 1e-10 -num_threads 60 -outfmt 6 -max_target_seqs 5
blastp -query CM334.fa -db index/PBC81 -out PBC81_CM334.blast -evalue 1e-10 -num_threads 60 -outfmt 6 -max_target_seqs 5
blastp -query Zunla-1.fa -db index/PBC81 -out PBC81_Zunla-1.blast -evalue 1e-10 -num_threads 60 -outfmt 6 -max_target_seqs 5
blastp -query PI159236.fa -db index/PBC81 -out PBC81_PI159236.blast -evalue 1e-10 -num_threads 60 -outfmt 6 -max_target_seqs 5

blastp -query CA59.fa -db index/Chiltepin -out Chiltepin_CA59.blast -evalue 1e-10 -num_threads 60 -outfmt 6 -max_target_seqs 5
blastp -query CM334.fa -db index/Chiltepin -out Chiltepin_CM334.blast -evalue 1e-10 -num_threads 60 -outfmt 6 -max_target_seqs 5
blastp -query Zunla-1.fa -db index/Chiltepin -out Chiltepin_Zunla-1.blast -evalue 1e-10 -num_threads 60 -outfmt 6 -max_target_seqs 5
blastp -query PI159236.fa -db index/Chiltepin -out Chiltepin_PI159236.blast -evalue 1e-10 -num_threads 60 -outfmt 6 -max_target_seqs 5
blastp -query PI_632928.fa -db index/Chiltepin -out Chiltepin_PI_632928.blast -evalue 1e-10 -num_threads 60 -outfmt 6 -max_target_seqs 5
blastp -query PI_632928.fa -db index/Grif_1614 -out Grif_1614_PI_632928.blast -evalue 1e-10 -num_threads 60 -outfmt 6 -max_target_seqs 1
blastp -query Grif_1614.fa -db index/Chiltepin -out Chiltepin_Grif_1614.blast -evalue 1e-10 -num_threads 60 -outfmt 6 -max_target_seqs 5

# Prepare the gene infos
cat S8.gene.gff PBC81.gene.gff | sed 's/ /\t/g' > S8_PBC81.gff
cat S8.gene.gff XLL.gene.gff | sed 's/ /\t/g' > S8_XLL.gff

cat Grif_1614.gene.gff PI_632928.gene.gff | sed 's/ /\t/g' > Grif_1614_PI_632928.gff
cat Chiltepin.gene.gff PBC81.gene.gff | sed 's/ /\t/g' > Chiltepin_PBC81.gff
cat Chiltepin.gene.gff PI_632928.gene.gff | sed 's/ /\t/g' > Chiltepin_PI_632928.gff
cat Chiltepin.gene.gff PBC81.gene.gff | sed 's/ /\t/g' > Chiltepin_PBC81.gff

cat CA59.gene.gff PBC81.gene.gff | sed 's/ /\t/g' > PBC81_CA59.gff
cat CM334.gene.gff PBC81.gene.gff | sed 's/ /\t/g' > PBC81_CM334.gff
cat Zunla-1.gene.gff PBC81.gene.gff | sed 's/ /\t/g' > PBC81_Zunla-1.gff
cat PI159236.gene.gff PBC81.gene.gff | sed 's/ /\t/g' > PBC81_PI159236.gff

cat CA59.gene.gff Chiltepin.gene.gff | sed 's/ /\t/g' > Chiltepin_CA59.gff
cat CM334.gene.gff Chiltepin.gene.gff | sed 's/ /\t/g' > Chiltepin_CM334.gff
cat Zunla-1.gene.gff Chiltepin.gene.gff | sed 's/ /\t/g' > Chiltepin_Zunla-1.gff
cat PI159236.gene.gff Chiltepin.gene.gff | sed 's/ /\t/g' > Chiltepin_PI159236.gff
cat PI159236.gene.gff Chiltepin.gene.gff | sed 's/ /\t/g' > Chiltepin_PI159236.gff
cat PI_632928.gene.gff Chiltepin.gene.gff | sed 's/ /\t/g' > Chiltepin_PI_632928.gff
cat Grif_1614.gene.gff Chiltepin.gene.gff | sed 's/ /\t/g' > Chiltepin_Grif_1614.gff

cat S8.gene.gff CA59.gene.gff | sed 's/ /\t/g' > S8_CA59.gff
cat S8.gene.gff Chiltepin.gene.gff | sed 's/ /\t/g' > S8_Chiltepin.gff
cat S8.gene.gff CM334.gene.gff | sed 's/ /\t/g' > S8_CM334.gff
cat S8.gene.gff PI159236.gene.gff | sed 's/ /\t/g' > S8_PI159236.gff
cat S8.gene.gff Zunla-1.gene.gff | sed 's/ /\t/g' > S8_Zunla-1.gff
cat Zunla-1.gene.gff Chiltepin.gene.gff | sed 's/ /\t/g' > Zunla-1_Chiltepin.gff
cat CM334.gene.gff Chiltepin.gene.gff | sed 's/ /\t/g' > CM334_Chiltepin.gff
cat PI159236.gene.gff Chiltepin.gene.gff | sed 's/ /\t/g' > PI159236_Chiltepin.gff
cat PBC81.gene.gff Chiltepin.gene.gff | sed 's/ /\t/g' > PBC81_Chiltepin.gff

# MCScan
MCScanX S8_PBC81 -a
MCScanX S8_XLL -a

MCScanX S8_PI_632928 -a
MCScanX S8_Grif_1614 -a
MCScanX PI_632928_S8 -a
MCScanX Grif_1614_S8 -a
MCScanX PI_632928_Zhangshugang -a
MCScanX Grif_1614_Zhangshugang -a
MCScanX Grif_1614_PI_632928 -a

MCScanX S8_CA59 -a
MCScanX S8_Chiltepin -a
MCScanX S8_CM334 -a
MCScanX S8_PI159236 -a
MCScanX S8_Zunla-1 -a
MCScanX Zunla-1_Chiltepin -a
MCScanX CM334_Chiltepin -a
MCScanX PI159236_Chiltepin -a

MCScanX PBC81_Chiltepin -a
MCScanX PBC81_CA59 -a
MCScanX PBC81_CM334 -a
MCScanX PBC81_Zunla-1 -a
MCScanX PBC81_PI159236 -a

MCScanX Chiltepin_CA59 -a
MCScanX Chiltepin_CM334 -a
MCScanX Chiltepin_Zunla-1 -a
MCScanX Chiltepin_PI159236 -a
MCScanX Chiltepin_PI_632928 -a
MCScanX Chiltepin_Grif_1614 -a

# Generate the plot
# Important: modify the control file from the program folder instead of creating one using vim
# S8 versus PBC81
java /data/zhaojiantao/tools/MCScanX-master/downstream_analyses/dot_plotter.java \
    -g S8_PBC81.gff \
    -s S8_PBC81.collinearity \
    -c S8_PBC81.dot.ctl \
    -o S8_PBC81.png

# S8 versus XLL
java /data/zhaojiantao/tools/MCScanX-master/downstream_analyses/dot_plotter.java \
    -g S8_XLL.gff \
    -s S8_XLL.collinearity \
    -c S8_XLL.dot.ctl \
    -o S8_XLL.png

# Zhangshugang versus PI_632928
java /data/zhaojiantao/tools/MCScanX-master/downstream_analyses/dot_plotter.java \
    -g PI_632928_Zhangshugang.gff \
    -s PI_632928_Zhangshugang.collinearity \
    -c PI_632928_Zhangshugang.dot.ctl \
    -o PI_632928_Zhangshugang.png

# Zhangshugang vs Grif_1614
java /data/zhaojiantao/tools/MCScanX-master/downstream_analyses/dot_plotter.java \
    -g Grif_1614_Zhangshugang.gff \
    -s Grif_1614_Zhangshugang.collinearity \
    -c Grif_1614_Zhangshugang.dot.ctl \
    -o Grif_1614_Zhangshugang.png
    
# Grif_1614 versus PI_632928
java /data/zhaojiantao/tools/MCScanX-master/downstream_analyses/dot_plotter.java \
    -g Grif_1614_PI_632928.gff \
    -s Grif_1614_PI_632928.collinearity \
    -c Grif_1614_PI_632928.dot.ctl \
    -o Grif_1614_PI_632928.png

# PI_632928 vs S8
java /data/zhaojiantao/tools/MCScanX-master/downstream_analyses/dot_plotter.java \
    -g PI_632928_S8.gff \
    -s PI_632928_S8.collinearity \
    -c PI_632928_S8.dot.ctl \
    -o PI_632928_S8.png

# S8 versus Grif_1614
java /data/zhaojiantao/tools/MCScanX-master/downstream_analyses/dot_plotter.java \
    -g S8_Grif_1614.gff \
    -s S8_Grif_1614.collinearity \
    -c S8_Grif_1614.dot.ctl \
    -o S8_Grif_1614.png

# Grif_1614 vs S8
java /data/zhaojiantao/tools/MCScanX-master/downstream_analyses/dot_plotter.java \
    -g Grif_1614_S8.gff \
    -s Grif_1614_S8.collinearity \
    -c Grif_1614_S8.dot.ctl \
    -o Grif_1614_S8.png

# S8 versus CA59
java /data/zhaojiantao/tools/MCScanX-master/downstream_analyses/dot_plotter.java \
    -g S8_CA59.gff \
    -s S8_CA59.collinearity \
    -c S8_CA59.dot.ctl \
    -o S8_CA59.png

# S8 versus Chiltepin
java /data/zhaojiantao/tools/MCScanX-master/downstream_analyses/dot_plotter.java \
    -g S8_Chiltepin.gff \
    -s S8_Chiltepin.collinearity \
    -c S8_Chiltepin.dot.ctl \
    -o S8_Chiltepin.png

# S8 versus CM334
java /data/zhaojiantao/tools/MCScanX-master/downstream_analyses/dot_plotter.java \
    -g S8_CM334.gff \
    -s S8_CM334.collinearity \
    -c S8_CM334.dot.ctl \
    -o S8_CM334.png

# S8 versus PI159236
java /data/zhaojiantao/tools/MCScanX-master/downstream_analyses/dot_plotter.java \
    -g S8_PI159236.gff \
    -s S8_PI159236.collinearity \
    -c S8_PI159236.dot.ctl \
    -o S8_PI159236.png

# S8 versus Zunla-1
java /data/zhaojiantao/tools/MCScanX-master/downstream_analyses/dot_plotter.java \
    -g S8_Zunla-1.gff \
    -s S8_Zunla-1.collinearity \
    -c S8_Zunla-1.dot.ctl \
    -o S8_Zunla-1.png

# PBC81 versus CA59
java /data/zhaojiantao/tools/MCScanX-master/downstream_analyses/dot_plotter.java \
    -g PBC81_CA59.gff \
    -s PBC81_CA59.collinearity \
    -c PBC81_CA59.dot.ctl \
    -o PBC81_CA59.png

# PBC81 versus CM334
java /data/zhaojiantao/tools/MCScanX-master/downstream_analyses/dot_plotter.java \
    -g PBC81_CM334.gff \
    -s PBC81_CM334.collinearity \
    -c PBC81_CM334.dot.ctl \
    -o PBC81_CM334.png

# PBC81 versus Zunla-1
java /data/zhaojiantao/tools/MCScanX-master/downstream_analyses/dot_plotter.java \
    -g PBC81_Zunla-1.gff \
    -s PBC81_Zunla-1.collinearity \
    -c PBC81_Zunla-1.dot.ctl \
    -o PBC81_Zunla-1.png

# PBC81 versus PI159236
java /data/zhaojiantao/tools/MCScanX-master/downstream_analyses/dot_plotter.java \
    -g PBC81_PI159236.gff \
    -s PBC81_PI159236.collinearity \
    -c PBC81_PI159236.dot.ctl \
    -o PBC81_PI159236.png

# PBC81 versus Chiltepin
java /data/zhaojiantao/tools/MCScanX-master/downstream_analyses/dot_plotter.java \
    -g PBC81_Chiltepin.gff \
    -s PBC81_Chiltepin.collinearity \
    -c PBC81_Chiltepin.dot.ctl \
    -o PBC81_Chiltepin.png

# Chiltepin versus CA59
java /data/zhaojiantao/tools/MCScanX-master/downstream_analyses/dot_plotter.java \
    -g Chiltepin_CA59.gff \
    -s Chiltepin_CA59.collinearity \
    -c Chiltepin_CA59.dot.ctl \
    -o Chiltepin_CA59.png

# Chiltepin versus Zunla-1
java /data/zhaojiantao/tools/MCScanX-master/downstream_analyses/dot_plotter.java \
    -g Chiltepin_Zunla-1.gff \
    -s Chiltepin_Zunla-1.collinearity \
    -c Chiltepin_Zunla-1.dot.ctl \
    -o Chiltepin_Zunla-1.png

# Chiltepin versus PI159236
java /data/zhaojiantao/tools/MCScanX-master/downstream_analyses/dot_plotter.java \
    -g Chiltepin_PI159236.gff \
    -s Chiltepin_PI159236.collinearity \
    -c Chiltepin_PI159236.dot.ctl \
    -o Chiltepin_PI159236.png

# Chiltepin versus PI_632928
java /data/zhaojiantao/tools/MCScanX-master/downstream_analyses/dot_plotter.java \
    -g Chiltepin_PI_632928.gff \
    -s Chiltepin_PI_632928.collinearity \
    -c Chiltepin_PI_632928.dot.ctl \
    -o Chiltepin_PI_632928.png

# Chiltepin versus Grif_1614
java /data/zhaojiantao/tools/MCScanX-master/downstream_analyses/dot_plotter.java \
    -g Chiltepin_Grif_1614.gff \
    -s Chiltepin_Grif_1614.collinearity \
    -c Chiltepin_Grif_1614.dot.ctl \
    -o Chiltepin_Grif_1614.png

# Use MCsan python version
# Prepare the gene bed file
# S8
sed 's/;/\t/g' /data/zhaojiantao/pepper/References/S8/Zhangshugang.gff | \
    awk '$3=="mRNA" {print $1,$4,$5,$9,0,$7}' | sed 's/ID=//g' | sed 's/ /\t/g' \
    > Zhangshugang.bed
# CA59
sed 's/;/\t/g' /data/zhaojiantao/pepper/References/CA59/Ca_59.v1a.gff3 | \
    awk '$3=="mRNA" {print $1,$4,$5,$9,0,$7}' | sed 's/ID=//g' | sed 's/ /\t/g' | \
    sed 's/Ca_59Chr/Chr/g' > CA59.bed

# Prepare the cds file
# NCEBR
sed 's/mRNA://g' /data/zhaojiantao/pepper/References/S8/Zhangshugang.cds.fasta > Zhangshugang.cds

# Pimpi
sed 's/mRNA://g' /data/zhaojiantao/pepper/References/CA59/Ca_59.cds.fa > CA59.cds 

# Run MCScan
python -m jcvi.compara.catalog ortholog Zhangshugang CA59 --no_strip_names

```

### Map the Chiltepin reads to the Chiltepin genome to check the read count

```bash
# Download the mate-paire raw Reads from [NCBI](https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SAMN01906552&o=acc_s%3Aa)
# SRR653457 has the largest insert size
fasterq-dump -p SRR653457

# Use NxTrim to remove mate-pair 
git clone https://github.com/sequencing/NxTrim.git
cd NxTrim
make

cd /data/zhaojiantao/pepper/References/Chiltepin/v2.0/SRA
# Extract the region
samtools faidx ../Chiltepin.Chr01.fasta
samtools faidx ../Chiltepin.Chr01.fasta Chr01:198941948-199029039 > Chr01.region.fasta
bwa index Chr01.region.fasta

# Trimm the adaptors
nxtrim -1 SRR653457_1.fastq -2 SRR653457_2.fastq -O SRR653457
# Reverse complementary sequences using [fastx_tookit](https://github.com/Debian/fastx-toolkit)
fastx_reverse_complement -z -i <(zcat SRR653457_R2.unknown.fastq.gz) -o SRR653457_R2.unknown.reverse.fastq.gz

# Combine the mate-pare and read pairs had unknown orientation
zcat SRR653457_R1.mp.fastq.gz SRR653457_R1.unknown.fastq.gz | gzip > SRR653457_R1.mp.unknown.fastq.gz
zcat SRR653457_R2.mp.fastq.gz SRR653457_R2.unknown.reverse.fastq.gz | gzip > SRR653457_R2.mp.unknown.fastq.gz
bwa aln -e 1 -t 80 ../Chiltepin.Chr01.fasta test_R1.unknown.fastq.gz > test_R1.sai
bwa aln -e 1 -t 80 ../Chiltepin.Chr01.fasta test_R2.unknown.reverse.fastq.gz > test_R2.sai


# Map the read to the references
bwa aln -e -1 -t 80 -n 1 ../Chiltepin.Chr01.fasta SRR653457.mp.uk.fastq.gz -f SRR653457.sai
bwa sampe -f SRR653457.sam ../Chiltepin.Chr01.fasta SRR653457.sai SRR653457.mp.uk.fastq.gz
samtools view -Shu SRR653457.sam | samtools sort -O BAM -@ 80 > SRR653457.sort.bam

# Extract the region
samtools view -b SRR653457.sort.bam "Chr01:198942247-199028668" > Chr01.region.bam
samtools index Chr01.region.bam

head -100000 SRR653457_1.fastq > test_1.fastq
head -100000 SRR653457_2.fastq > test_2.fastq
nxtrim -1 test_1.fastq -2 test_2.fastq --separate -O test
zcat SRR653457_R1.unknown.fastq.gz | head -1000 | gzip > test_R1.unknown.fastq.gz
zcat SRR653457_R2.unknown.fastq.gz | head -1000 > test_R2.unknown.fastq
fastx_reverse_complement -z -i test_R2.unknown.fastq -o test_R2.unknown.reverse.fastq.gz

bwa aln -e 1 -t 80 ../Chiltepin.Chr01.fasta test_R1.unknown.fastq.gz > test_R1.sai
bwa aln -e 1 -t 80 ../Chiltepin.Chr01.fasta test_R2.unknown.reverse.fastq.gz > test_R2.sai
bwa sampe ../Chiltepin.Chr01.fasta test_R1.sai test_R2.sai test_R1.unknown.fastq.gz test_R2.unknown.reverse.fastq.gz > test.sam
samtools view -Shu test.sam | samtools sort -O BAM -@ 80 > test.sort.bam

```

### Align sequence to S8

```bash
# Working directory
cd /data/zhaojiantao/pepper/References/S8

# Install blast
wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.11.0/ncbi-blast-2.11.0+-x64-linux.tar.gz
tar -xvzf ncbi-blast-2.11.0+-x64-linux.tar.gz
mv ncbi-blast-2.11.0+ blast-2.11.0+
export PATH="/data/zhaojiantao/tools/blast-2.11.0+/bin:$PATH"

# Build S8 database
# Nuclei database
makeblastdb -in /data/zhaojiantao/pepper/References/S8/pepper_S8.fasta \
    -dbtype nucl -parse_seqids \
    -out /data/zhaojiantao/pepper/References/S8/blastdb/index

makeblastdb -in /data/zhaojiantao/pepper/References/S8/Zhangshugang.fasta \
    -dbtype nucl -parse_seqids \
    -out /data/zhaojiantao/pepper/References/S8/blastdb/Zhangshugang.index

makeblastdb -in /data/zhaojiantao/pepper/References/S8/Zhangshugang.cds.fa \
    -dbtype nucl -parse_seqids \
    -out /data/zhaojiantao/pepper/References/S8/blastdb/S8.cds.index

makeblastdb -in /data/zhaojiantao/pepper/References/CM334/assemblies/Pepper_1.55.fasta \
    -dbtype nucl -parse_seqids \
    -out /data/zhaojiantao/pepper/References/CM334/assemblies/blastdb/index
makeblastdb -in /data/zhaojiantao/pepper/References/Zunla-1/v1.0/Zunla-1.fasta \
    -dbtype nucl -parse_seqids \
    -out /data/zhaojiantao/pepper/References/Zunla-1/v1.0/blastdb/index

# Protein database
makeblastdb -in /data/zhaojiantao/pepper/References/CM334/annotation/Pepper.v.1.55.PEP.fa \
    -dbtype prot -parse_seqids -out /data/zhaojiantao/pepper/References/CM334/annotation/blastdb/index

makeblastdb -in /data/zhaojiantao/pepper/References/CM334/v2.0/CM334.PEP.fa \
    -dbtype prot -parse_seqids -out /data/zhaojiantao/pepper/References/CM334/v2.0/blastdb/index

makeblastdb -in /data/zhaojiantao/pepper/References/S8/gene/pepper_S8.pep.fa \
    -dbtype prot -parse_seqids -out /data/zhaojiantao/pepper/References/S8/gene/blastdb/index  
makeblastdb -in /data/zhaojiantao/pepper/References/Zunla-1/v1.0/Zunla-1.pep.fa \
    -dbtype prot -parse_seqids -out /data/zhaojiantao/pepper/References/Zunla-1/v1.0/blastdb_protein/index

# Blast the position on CM334
blastn -query AT3-2.fa \
    -db /data/zhaojiantao/pepper/References/CM334/assemblies/blastdb/index \
    -evalue 1e-20 -num_threads 6 \
    -out AT3-2.CM334.blast

blastn -query P1.306622565.fa \
    -db /data/zhaojiantao/pepper/References/S8/blastdb/index \
    -evalue 1e-20 -num_threads 6 \
    -out P1.306622565.start.blast

blastn -query pATM.fa \
    -db /data/zhaojiantao/pepper/References/Zunla-1/v1.0/blastdb/index \
    -evalue 1e-20 -num_threads 6 \
    -out ZL1_1826.Zunla-1.blast

blastp -query XP_016580728.1.fa \
    -db /data/zhaojiantao/pepper/References/CM334/annotation/blastdb/index \
    -evalue 1e-20 -num_threads 6 \
    -out XP_016580728.CM334.pep.blast

blastp -query Ca9g33797.pro.fa \
    -db /data/zhaojiantao/pepper/References/CM334/v2.0/blastdb/index \
    -evalue 1e-20 -num_threads 6 \
    -out Ca9g33797.pro.blast

# Find the gene ID
grep Ca2g12681 /data/zhaojiantao/pepper/References/S8/gene/pepper_S8.cds.fa | sed 's/>//g' > gene.id
grep CA11g14610 /data/zhaojiantao/pepper/References/CM334/annotation/Pepper.v.1.55.PEP.fa | sed 's/>//g' > gene.id
seqtk subseq /data/zhaojiantao/pepper/References/CM334/annotation/Pepper.v.1.55.PEP.fa gene.id > CA11g14610.CM334.fa
grep CA04g00860 /data/zhaojiantao/pepper/References/CM334/annotation/Pepper_1.55.gene_models.gff3

grep CA04g00860 /data/zhaojiantao/pepper/References/S8/gene/pepper_S8.gff3

seqtk subseq \
    /data/zhaojiantao/pepper/References/Zunla-1/v1.0/Zunla-1.cds.fasta \
    Zunla.pungency.genes.list \
    > Zunla.pungency.genes.cds.fa

# Split the fasta
awk '/^>/{s=++d".fasta"} {print > s}' Zunla.pungency.genes.cds.fa

# Blast for each gene
for i in $(seq 1 154); do
    # to S8
    blastn -query $i.fasta \
        -db /data/zhaojiantao/pepper/References/S8/blastdb/index \
        -evalue 1e-20 -num_threads 6 \
        -out $i.blast
    # to cds
    blastn -query $i.fasta \
        -db /data/zhaojiantao/pepper/References/S8/blastdb/S8.cds.index \
        -evalue 1e-20 -num_threads 60 -outfmt 6 \
        -out $i.blast3    
done

#----------------------------------------------------------------------------------------------------
# For the comparison of Cao et al. 2022



# Find genome synteny using [mcscan](https://github.com/tanghaibao/jcvi/wiki/MCscan-(Python-version))
# Convert the GFF to BED file and rename them.
cd /data/zhaojiantao/pepper/References/MCScanX

cut -f1 ../CM334/v2.0/Annuum.v.2.0.gff3 | sed 's/PGAv.1.6.scaffold/chr/g' | \
    paste - ../CM334/v2.0/Annuum.v.2.0.gff3 | grep gene | \
    awk '{print $1,$5,$6,$10,$7,$8}' | sed 's/ID=//g' | sed 's/ /\t/g' \
    > CM334v2.bed
    
zcat ../S8/pepper_S8.gff3.gz | sed 's/;/\t/g' | grep ID=gene | \
    awk '{print $1,$4,$5,$9,$6,$7}' | sed 's/Superscaffold/chr/g' | \
    sed 's/ID=gene-//g' | sed 's/ /\t/g' \
    > S8.bed

# Reformat fasta files
python -m jcvi.formats.fasta format ../CM334/v2.0/Annuum.v.2.0.CDS.fa.gz CM334v2.cds
python -m jcvi.formats.fasta format ../S8/gene/pepper_S8.cds.fa S8.cds
sed 's/\.1//g' S8.cds | sed 's/rna-//g' > S8.cds2
rm S8.cds
mv S8.cds2 S8.cds

# Pairwise synteny search
python -m jcvi.compara.catalog ortholog CM334v2 S8 --no_strip_names
python -m jcvi.compara.catalog ortholog CM334v2 S8 --cscore=.99 --no_strip_names
python -m jcvi.compara.synteny screen --minspan=30 --simple CM334v2.S8.anchors CM334v2.S8.anchors.new 

```
## vg pangenome

### Trimmomatic

```bash
# Working directory
cd /data/zhaojiantao/pepper/vg/00fastq

# Trim adaptors and low quality reads
for i in $(cat sample.list) ; do
    trimmomatic PE $i.R1.fastq.gz $i.R2.fastq.gz \
        $i.R1P.fq $i.R1U.fq \
        $i.R2P.fq $i.R2U.fq \
        -threads 60 \
        ILLUMINACLIP:/data/zhaojiantao/tools/miniconda3/share/trimmomatic-0.39-2/adapters/TruSeq3-PE-2.fa:2:30:10:1:TRUE \
        SLIDINGWINDOW:4:20 LEADING:3 TRAILING:3 MINLEN:40 ;
    # Remove temporary files
    rm $i.R1.fastq.gz $i.R2.fastq.gz $i.R1U.fq $i.R2U.fq
done

# From the cornell server
# Trim adaptors and low quality reads
for i in bailajiao ; do
    trimmomatic PE $i.R1.fastq.gz $i.R2.fastq.gz \
        /workdir/zf25/01fastq.clean/$i.R1P.fq /workdir/zf25/01fastq.clean/$i.R1U.fq \
        /workdir/zf25/01fastq.clean/$i.R2P.fq /workdir/zf25/01fastq.clean/$i.R2U.fq \
        -threads 60 \
        ILLUMINACLIP:/home/zf25/miniconda3/share/trimmomatic-0.39-2/adapters/TruSeq3-PE-2.fa:2:30:10:1:TRUE \
        SLIDINGWINDOW:4:20 LEADING:3 TRAILING:3 MINLEN:40 ;
    # Remove temporary files
    rm /workdir/zf25/01fastq.clean/$i.R1U.fq /workdir/zf25/01fastq.clean/$i.R2U.fq
done

nohup sh trim.sh 1>trim.log1 2>trim.log &
done

```

### SV estimation

```bash
# Working directory
cd /data/zhaojiantao/pepper/vg/PanReference

# PI_632928
grep -v "#" m6a.4gir.vcf | grep PASS | cut -f4,5,10 | grep -v "\." | awk '{print length($1),length($2)}' > PI_632928.alt.length

# Grif_1614
grep -v "#" m6a.4gir.vcf | grep PASS | cut -f4,5,11 | grep -v "\." | awk '{print length($1),length($2)}'  > Grif_1614.alt.length

# Calculate SV <=10k
awk '$1<=10000 {sum += $7 } END { print sum }' PI_632928.alt.length

```


### vg mapping

```bash
# Working directory
cd /data/zhaojiantao/pepper/vg/01bam

# vg mapping
for i in PI_355813 PI_355815 PI_355817 PI_355819 PI_355820 PI_357467 PI_358811 PI_358968; do
    # Trim
    trimmomatic PE ../00fastq/$i.R1.fastq.gz ../00fastq/$i.R2.fastq.gz \
    ../00fastq/$i.R1P.fq ../00fastq/$i.R1U.fq \
    ../00fastq/$i.R2P.fq ../00fastq/$i.R2U.fq \
    -threads 60 \
    ILLUMINACLIP:/data/zhaojiantao/tools/miniconda3/share/trimmomatic-0.39-2/adapters/TruSeq3-PE-2.fa:2:30:10:1:TRUE \
    SLIDINGWINDOW:4:20 LEADING:3 TRAILING:3 MINLEN:40 2>../00fastq/$i.trim.log
    # Remove temporary files
    rm ../00fastq/$i.R1.fastq.gz ../00fastq/$i.R2.fastq.gz ../00fastq/$i.R1U.fq ../00fastq/$i.R2U.fq
    # Generate gam file
    vg giraffe -t 60 \
        -Z m6a.giraffe.gbz \
        -m m6a.min -d m6a.dist \
        -f ../00fastq/$i.R1P.fq \
        -f ../00fastq/$i.R2P.fq \
        > $i.gam
    # Generate bam file
    vg surject -t 60 \
        -N $i \
        -R $i \
        -x m6a.giraffe.gbz \
        -b -i $i.gam \
        > $i.gam.bam
    # Sort the bam
    perl sort_bam_header_chrID.pl $i $i.gam.bam ${i}_gam.bam
    # Index the bam
 	samtools index -@ 60 ${i}_gam.bam
done

```


### vg mapping on cornell server

```bash
# Working directory
cd /workdir/zf25

# vg mapping
for i in $(cat list) ; do echo "
date
trimmomatic PE /workdir/zf25/00fastq/$i.R1.fastq.gz /workdir/zf25/00fastq/$i.R2.fastq.gz \
    /workdir/zf25/00fastq/$i.R1P.fq /workdir/zf25/00fastq/$i.R1U.fq \
    /workdir/zf25/00fastq/$i.R2P.fq /workdir/zf25/00fastq/$i.R2U.fq \
    -threads 110 \
    ILLUMINACLIP:/home/zf25/miniconda3/share/trimmomatic-0.39-2/adapters/TruSeq3-PE-2.fa:2:30:10:1:TRUE \
    SLIDINGWINDOW:4:20 LEADING:3 TRAILING:3 MINLEN:40 2>/workdir/zf25/00fastq/$i.trim.log

rm /workdir/zf25/00fastq/$i.R1U.fq /workdir/zf25/00fastq/$i.R2U.fq
rm /workdir/zf25/00fastq/$i.R1.fastq.gz /workdir/zf25/00fastq/$i.R2.fastq.gz

vg giraffe -t 110 \
    -Z /workdir/zf25/db/m6a.giraffe.gbz \
    -m /workdir/zf25/db/m6a.min \
    -d /workdir/zf25/db/m6a.dist \
    -f 00fastq/$i.R1P.fq \
    -f 00fastq/$i.R2P.fq \
    > $i.gam

vg surject -t 110 \
    -N $i \
    -R $i \
    -x /workdir/zf25/db/m6a.giraffe.gbz \
    -b -i $i.gam \
    > $i.gam.bam

perl /home/zf25/sort_bam_header_chrID.pl $i $i.gam.bam ${i}_gam.bam

samtools index -@ 110 ${i}_gam.bam
date " 
done | sed 's/         / /g' | sed 's/     / /g' > batch6.sh

```


### SNP calling using sentieon

```bash
# Working directory
cd /data/zhaojiantao/pepper/vg/02gvcf

# Run sentieon
# Call variants
for i in $(cat list) ; do
    # Collect read information to remove or mark duplicates
    sentieon driver -t 100 -i ../01bam/${i}_gam.bam \
        --algo LocusCollector --fun score_info $i.score.gz
    # Remove or mark duplicates
    sentieon driver -t 100 -i ../01bam/${i}_gam.bam \
        --algo Dedup --score_info $i.score.gz \
        --metrics $i.dedup_metrix.txt $i.deduped.bam
    # Call gvcf
    sentieon driver -t 100 \
        -r ../PanReference/Cann_Zhangshugang.chr.fa \
        -i $i.deduped.bam \
        --algo Haplotyper \
        --emit_mode gvcf \
        $i.gvcf.gz
    # Remove unnecessary files
    rm $i.score.gz* $i.dedup_metrix.txt $i.deduped.bam*
done

```

### SNP calling using sentieon for tomato as outgroup

```bash
# Working directory
cd /data/zhaojiantao/pepper/vg/02gvcf

# Run sentieon
# Download tomato from NCBI as outgroup
for i in ERR418102 ERR418103 ERR418107 ERR418106; do
    # Download SRA
    prefetch $i
    # Convert SRA to fastq
    fasterq-dump --split-files  $i/$i.sra
    # Remove adaptor and low-quality sequences
    trimmomatic PE ${i}_1.fastq ${i}_2.fastq \
            $i/$i.R1P.fq $i/$i.R1U.fq \
            $i/$i.R2P.fq $i/$i.R2U.fq \
            -threads 100 \
            ILLUMINACLIP:/data/zhaojiantao/tools/miniconda3/share/trimmomatic-0.39-2/adapters/TruSeq3-PE-2.fa:2:30:10:1:TRUE \
            SLIDINGWINDOW:4:20 LEADING:3 TRAILING:3 MINLEN:40
    # Remove temporary file
    rm $i/$i.sra $i/$i.R1U.fq $i/$i.R2U.fq ${i}_1.fastq ${i}_2.fastq 
    # vg mapping
    vg giraffe -t 100 \
        -Z /workdir/zf25/db/m6a.giraffe.gbz \
        -m /workdir/zf25/db/m6a.min \
        -d /workdir/zf25/db/m6a.dist \
        -f 00fastq/$i.R1P.fq \
        -f 00fastq/$i.R2P.fq \
        > $i.gam
    # Generate bam
    vg surject -t 100 \
        -N $i \
        -R $i \
        -x /workdir/zf25/db/m6a.giraffe.gbz \
        -b -i $i.gam \
        > $i.gam.bam
    # Sort the bam file
    perl /home/zf25/sort_bam_header_chrID.pl $i $i.gam.bam ${i}_gam.bam
    # Index the bam
    samtools index -@ 110 ${i}_gam.bam        
    # Call variants
    # Collect read information to remove or mark duplicates
    sentieon driver -t 100 -i ../01bam/${i}_gam.bam \
        --algo LocusCollector --fun score_info $i.score.gz
    # Remove or mark duplicates
    sentieon driver -t 100 -i ../01bam/${i}_gam.bam \
        --algo Dedup --score_info $i.score.gz \
        --metrics $i.dedup_metrix.txt $i.deduped.bam
    # Call gvcf
    sentieon driver -t 100 \
        -r ../PanReference/Cann_Zhangshugang.chr.fa \
        -i $i.deduped.bam \
        --algo Haplotyper \
        --emit_mode gvcf \
        $i.gvcf.gz
    # Remove unnecessary files
    rm $i.score.gz* $i.dedup_metrix.txt $i.deduped.bam*
done

```

### SNP calling using sentieon on cornell server
```bash
# work directory
cd /workdir/zf25/02gvcf

for i in bailajiao ; do echo "
    # Collect read information to remove or mark duplicates
    /programs/sentieon-genomics-202112/bin/sentieon driver -t 110 -i ../01bam/${i}_gam.bam \
        --algo LocusCollector --fun score_info $i.score.gz
    # Remove or mark duplicates
    /programs/sentieon-genomics-202112/bin/sentieon driver -t 110 -i ../01bam/${i}_gam.bam \
        --algo Dedup --score_info $i.score.gz \
        --metrics $i.dedup_metrix.txt $i.deduped.bam
    # Call gvcf
    /programs/sentieon-genomics-202112/bin/sentieon driver -t 110 \
        -r /home/zf25/PanReference/Cann_Zhangshugang.chr.fa \
        -i $i.deduped.bam \
        --algo Haplotyper \
        --emit_mode gvcf \
        $i.gvcf.gz
    # Remove unnecessary files
    rm $i.score.gz* $i.dedup_metrix.txt $i.deduped.bam* "
done | sed 's/         / /g' | sed 's/     / /g' | sed 's/    //g' > batch1.sh

```

### Join the vcf file

```bash
# Working directory
cd /data/zhaojiantao/pepper/vg/02gvcf

# Joint variant call
sentieon driver -t 104 \
    -r /data/zhaojiantao/pepper/vg/PanReference/Cann_Zhangshugang.chr.fa --algo GVCFtyper \
    ../03vcf/500.vcf.gz *gvcf.gz

```

### Extract the high quality SNPs and Indels

```bash
# Working directory
cd /data/zhaojiantao/pepper/vg/03vcf

# Create the FASTA sequence directory file of the reference to be used in gatk
java -jar /data/zhaojiantao/tools/picard/build/libs/picard.jar CreateSequenceDictionary \
    R=/data/zhaojiantao/pepper/vg/PanReference/Cann_Zhangshugang.chr.fa \
    O=/data/zhaojiantao/pepper/vg/PanReference/Cann_Zhangshugang.chr.dict

# Split the file into chromosomes
for i in $(seq -w 1 12); do
    bcftools view 500.vcf.gz --regions CannZSG_Chr$i | bgzip -@ 100 > Chr$i.vcf.gz;
    tabix -p vcf Chr$i.vcf.gz
done

# For contigs
zcat 500.vcf.gz | grep "#" | grep CannZSG_Contig | \
    sed 's/,/ /g' | sed 's/##contig=<ID=//g' | sed 's/length=//g' | awk '{print $1,1,$2}' | sed 's/ /\t/g' \
    > Contig.list
bcftools view 500.vcf.gz --regions-file Contig.list | bgzip -@ 20 > Chr00.vcf.gz
tabix -p vcf Chr00.vcf.gz

# Working directory
cd /data/zhaojiantao/pepper/vg/04SNP

# SNPs
for i in $(seq -w 0 12); do
    # Select SNP
    gatk SelectVariants \
        -R /data/zhaojiantao/pepper/vg/PanReference/Cann_Zhangshugang.chr.fa \
        -V ../03vcf/Chr$i.vcf.gz \
        -select-type SNP \
        -O Chr$i.SNP.all.vcf.gz
    # Hard filter for SNPs
    gatk VariantFiltration \
        -R /data/zhaojiantao/pepper/vg/PanReference/Cann_Zhangshugang.chr.fa \
        -V Chr$i.SNP.all.vcf.gz \
        --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || SOR >3.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
        --filter-name "my_snp_filter" \
        -O Chr$i.filter.vcf.gz
    # Extract the hard filtered
    bcftools view Chr$i.filter.vcf.gz \
        --apply-filters PASS | \
        bgzip -@ 20 \
        > Chr$i.filtered.vcf.gz
    # Index
    tabix -p vcf Chr$i.filtered.vcf.gz
    # Remove loci with missing <0.1 and maf < 0.05
    vcftools --gzvcf Chr$i.filtered.vcf.gz --maf 0.05 --max-missing 0.1 --recode --stdout | bgzip -@ 20 > Chr$i.filtered2.vcf.gz
    # Index
    tabix -p vcf Chr$i.filtered2.vcf.gz
    # Extract the high quality SNPs
    bcftools view Chr$i.filtered2.vcf.gz --types snps -m 2 -M 2 | bgzip -@ 20 > Chr$i.biSNP.vcf.gz
    tabix -p vcf Chr$i.biSNP.vcf.gz
    # Final SNP set with missing cut of 0.3
    bcftools view -i 'F_MISSING<0.3' Chr$i.biSNP.vcf.gz | bgzip -@ 20 > Chr$i.SNP.vcf.gz
done

# Indels
for i in $(seq -w 0 12); do
    # Select Indel
    gatk SelectVariants \
        -R /data/zhaojiantao/pepper/vg/PanReference/Cann_Zhangshugang.chr.fa \
        -V ../03vcf/Chr$i.vcf.gz \
        -select-type INDEL \
        -O Chr$i.INDEL.all.vcf.gz
    # GATK hard filter
    gatk VariantFiltration \
        -R /data/zhaojiantao/pepper/vg/PanReference/Cann_Zhangshugang.chr.fa \
        -V Chr$i.INDEL.all.vcf.gz \
        --filter-expression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0" \
        --filter-name "my_snp_filter" \
        -O Chr$i.INDEL.filter.vcf.gz
    # Extract the Indels passing the hard filter and further filter
    bcftools view Chr$i.INDEL.filtered.vcf.gz --apply-filters PASS \
        --types indels -m 2 -M 2 -i 'F_MISSING<0.3' --min-af 0.05 | \
        bgzip -@ 20 \
        > Chr$i.INDEL.filtered2.vcf.gz
    # Rename the chromosome ids
    bcftools annotate --rename-chrs chr.update.txt Chr$i.INDEL.filtered2.vcf.gz | bgzip -@ 20 > Chr$i.INDEL.final.vcf.gz &
done

## Missing<0.3
for i in $(seq -w 0 12); do 
bcftools view -i 'F_MISSING<0.3' ../Chr$i.SNP.vcf.gz | bgzip -@ 20 > Chr$i.SNP.vcf.gz
bcftools query -f 'x' Chr$i.SNP.vcf.gz | wc -c > Chr$i.SNP.count
bcftools query -f '%CHROM %POS\n' Chr$i.SNP.vcf.gz | sed 's/CannZSG_//g' > Chr$i.SNP.info
grep -f Chr$i.SNP.info ../Chr$i.SNP.info.old | wc -l > Chr$i.SNP.consis
done

# Add SNP id to the vcf
for i in $(seq -w 0 12); do
bcftools annotate --set-id '%CHROM\_%POS' Chr$i.SNP.vcf.gz | bgzip -@ 20 > Chr$i.SNP.new.vcf.gz
mv Chr$i.SNP.new.vcf.gz Chr$i.SNP.vcf.gz
done

# Count the number of SNPs and Indels
for i in $(seq -w 0 12); do
bcftools query -f 'x' Chr$i.SNP.vcf.gz | wc -c > Chr$i.SNP.count
bcftools query -f 'x' Chr$i.INDEL.final.vcf.gz | wc -c > Chr$i.INDEL.final.count
done

# Check the consistency of SNPs
for i in $(seq -w 1 12); do
# SNPs from pan
bcftools query -f '%CHROM %POS %REF %ALT\n' Chr$i.SNP.vcf.gz | sed 's/CannZSG_//g' | cut -f1,2 -d ' ' > Chr$i.SNP.info
# SNP from old
bcftools query -f '%CHROM %POS %REF %ALT\n' ../../vcf.clean/Chr$i.clean.vcf.gz  | cut -f1,2 -d ' ' > Chr$i.SNP.info.old
# Rotal consistent
grep -f Chr$i.SNP.info Chr$i.SNP.info.old | wc -l > Chr$i.SNP.consis
grep -vf Chr$i.SNP.info Chr$i.SNP.info.old | wc -l > Chr$i.SNP.inconsis
done

# Check the consistency of Indels
for i in $(seq -w 0 12); do
# Indels from pan
bcftools query -f '%CHROM %POS %REF %ALT\n' Chr$i.INDEL.vcf.gz | sed 's/CannZSG_//g' | cut -f1,2 -d ' ' > Chr$i.Indels.info
# Indels from old
bcftools query -f '%CHROM %POS %REF %ALT\n' ../../Indels/Chr$i.fltr5.indel.vcf.gz  | cut -f1,2 -d ' ' > Chr$i.Indels.info.old
# Total consistent
grep -f Chr$i.Indels.info Chr$i.Indels.info.old | wc -l > Chr$i.Indels.consis
grep -vf Chr$i.Indels.info Chr$i.Indels.info.old | wc -l > Chr$i.Indels.inconsis
done

# Save the file to another folder
mv *SNP* ../04SNP
mv *INDEL* ../05INDEL

```

### SNP and Indel snpeff

```bash
# Working directory
cd /data/zhaojiantao/tools/snpEff

# Prepare and bulid the genome index
cp /data/zhaojiantao/pepper/References/S8/Zhangshugang.fa data/genomes/Zhangshugang.fa
cp /data/zhaojiantao/pepper/References/S8/Zhangshugang.fa data/Zhangshugang/Zhangshugang.fa
cp /data/zhaojiantao/pepper/References/S8/Zhangshugang.fa data/Zhangshugang/sequences.fa
cp /data/zhaojiantao/pepper/References/S8/Zhangshugang.gff3 data/Zhangshugang/genes.gff

echo "Zhangshugang.genome : Zhangshugang" >> snpEff.config

# Build the database
java -jar snpEff.jar build -gff3 -v Zhangshugang

# Working directory for SNPs
cd /data/zhaojiantao/pepper/vg/04SNP/Miss_0.3

# Rename the chromosome ids
for i in $(seq -w 0 12); do
    bcftools annotate --rename-chrs Old/chr.update.txt Old/Chr$i.SNP.vcf.gz | bgzip -@ 20 > Chr$i.SNP.vcf.gz &
done

# Calculate the SNP effects
for i in $(seq -w 0 12); do
    cd /data/zhaojiantao/pepper/vg/04SNP/Miss_0.3/SnpEff.Chr$i
    java -Xmx8g -jar /data/zhaojiantao/tools/snpEff/snpEff.jar \
        -c /data/zhaojiantao/tools/snpEff/snpEff.config \
        Zhangshugang \
        ../Chr$i.SNP.vcf.gz | bgzip -@ 20 \
        > Chr$i.ann.vcf.gz &
done

# Working directory for Indels
cd /data/zhaojiantao/pepper/vg/05INDEL

# Rename the chromosome ids
for i in $(seq -w 0 12); do
    bcftools annotate --rename-chrs Old/chr.update.txt Old/Chr$i.INDEL.vcf.gz | bgzip -@ 20 > Chr$i.INDEL.vcf.gz &
done

# Calculate indel effects
for i in $(seq -w 0 12); do
    cd /data/zhaojiantao/pepper/vg/05INDEL/SnpEff.indel.Chr$i
    java -Xmx8g -jar /data/zhaojiantao/tools/snpEff/snpEff.jar \
        -c /data/zhaojiantao/tools/snpEff/snpEff.config \
        Zhangshugang \
        ../Chr$i.INDEL.final.vcf.gz | bgzip -@ 20 \
        > Chr$i.indel.ann.vcf.gz &
done

# Calculate the uniq variants by taking the first effect
cd /data/zhaojiantao/pepper/vg/04SNP/Miss_0.3
for i in $(seq -w 0 12); do
    zcat SnpEff.Chr$i/Chr$i.ann.vcf.gz | grep -v "#" | cut -f8 | sed 's/|/\t/g' | cut -f2 | sort | uniq -c > SnpEff.Chr$i/Chr$i.ann.count &
done

# Calculate the uniq variants by taking the first effect
cd /data/zhaojiantao/pepper/vg/05INDEL
for i in $(seq -w 0 12); do
    zcat SnpEff.indel.Chr$i/Chr$i.indel.ann.vcf.gz | grep -v "#" | cut -f8 | sed 's/;/\t/g' | \
    awk '{for(i=1;i<=NF;i++){if ($i ~ /ANN/){print $i}}}' | \
    sed 's/|/\t/g' | cut -f2 | sort | uniq -c > SnpEff.indel.Chr$i/Chr$i.indel.ann.count &
done

```

#### SNP missing rate cutoff

```bash
# Working directory
cd /data/zhaojiantao/pepper/vg/04SNP

bcftools view -i 'F_MISSING<0.1' Chr01.SNP.vcf.gz | bgzip > Chr01.Miss_0.1.SNP.vcf.gz &
bcftools view -i 'F_MISSING<0.2' Chr01.SNP.vcf.gz | bgzip > Chr01.Miss_0.2.SNP.vcf.gz &
bcftools view -i 'F_MISSING<0.3' Chr01.SNP.vcf.gz | bgzip > Chr01.Miss_0.3.SNP.vcf.gz &
bcftools view -i 'F_MISSING<0.4' Chr01.SNP.vcf.gz | bgzip > Chr01.Miss_0.4.SNP.vcf.gz &
bcftools view -i 'F_MISSING<0.5' Chr01.SNP.vcf.gz | bgzip > Chr01.Miss_0.5.SNP.vcf.gz &
bcftools view -i 'F_MISSING<0.6' Chr01.SNP.vcf.gz | bgzip > Chr01.Miss_0.6.SNP.vcf.gz &
bcftools view -i 'F_MISSING<0.7' Chr01.SNP.vcf.gz | bgzip > Chr01.Miss_0.7.SNP.vcf.gz &
bcftools view -i 'F_MISSING<0.8' Chr01.SNP.vcf.gz | bgzip > Chr01.Miss_0.8.SNP.vcf.gz &
bcftools view -i 'F_MISSING<0.9' Chr01.SNP.vcf.gz | bgzip > Chr01.Miss_0.9.SNP.vcf.gz &

for i in $(seq 1 9); do 
tabix -p vcf Chr01.Miss_0.$i.SNP.vcf.gz &
bcftools query -f 'x' Chr01.Miss_0.$i.SNP.vcf.gz  | wc -c > Chr01.Miss_0.$i.SNP.count
done

# Compare with vcftools
for i in $(seq 1 9); do 
vcftools --gzvcf Chr01.SNP.vcf.gz --max-missing 0.$i --maf 0.05 --out Chr01.Miss_0.$i.SNP.vcftools.count &
done


```

#### SNP check

```bash
# Working directory
cd /data/zhaojiantao/pepper/vg/04SNP/SNP.check

# Rename the sample id
for i in chaotiaojiao Grif_1615 Grif_9355 PI_424732; do
    # Rename
    bcftools view ${i}_bwa.gvcf.gz | \
        sed "s/${i}_1/${i}_bwa/g" | bgzip -@ 20 \
        > ${i}_bwa.rename.gvcf.gz
    # Index
    tabix -p vcf ${i}_bwa.rename.gvcf.gz
done

for i in chaotiaojiao Grif_1615 Grif_9355 PI_424732; do
    # Rename
    bcftools view ${i}_gam.gvcf.gz | \
        sed "s/$i/${i}_gam/g" | bgzip -@ 20 \
        > ${i}_gam.rename.gvcf.gz
    # Index
    tabix -p vcf ${i}_gam.rename.gvcf.gz
done

for i in chaotiaojiao Grif_1615 Grif_9355 PI_424732; do
    # Rename
    bcftools view $i.sppdup.g.vcf.gz | \
        sed "s/$i/$i.sppdup/g" | bgzip -@ 20 \
        > $i.sppdup.g.rename.vcf.gz
    # Index
    tabix -p vcf $i.sppdup.g.rename.vcf.gz
done

for i in chaotiaojiao Grif_1615 Grif_9355 PI_424732; do
    # Rename
    bcftools view $i.gvcf.gz | \
        sed "s/$i/$i.vg/g" | bgzip -@ 20 \
        > $i.rename.gvcf.gz
    # Index
    tabix -p vcf $i.rename.gvcf.gz
done

# Join the gvcf files
sentieon driver -t 100 \
    -r /data/zhaojiantao/pepper/vg/PanReference/Cann_Zhangshugang.chr.fa --algo GVCFtyper \
    Check.vcf.gz *rename*gz
# Hard filter
gatk VariantFiltration \
    -R /data/zhaojiantao/pepper/vg/PanReference/Cann_Zhangshugang.chr.fa \
    -V Check.vcf.gz \
    --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || SOR >3.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
    --filter-name "my_snp_filter" \
    -O Check.filter.vcf.gz
# Extract the hard filtered
bcftools view Check.filter.vcf.gz \
    --apply-filters PASS | \
    bgzip -@ 20 \
    > Check.filtered.vcf.gz
# Index
tabix -p vcf Check.filtered.vcf.gz
# Extract biallelic    
bcftools view Check.filtered.vcf.gz --types snps -m 2 -M 2 | bgzip -@ 20 > Check.filtered.biSNP.vcf.gz
tabix -p vcf Check.filtered.biSNP.vcf.gz

# Convert to plink biallelic
plink --vcf Check.filtered.biSNP.vcf.gz \
    --recode 12 --transpose \
    --double-id \
    --allow-extra-chr \
    --out Check.filtered.biSNP

# Count the number of missing
awk '$5=="0" && $6=="0" {print}' Check.filtered.biSNP.tped | wc -l > A01.missing &
awk '$7=="0" && $8=="0" {print}' Check.filtered.biSNP.tped | wc -l > A02.missing &
awk '$9=="0" && $10=="0" {print}' Check.filtered.biSNP.tped | wc -l > A03.missing &
awk '$11=="0" && $12=="0" {print}' Check.filtered.biSNP.tped | wc -l > A04.missing &
awk '$13=="0" && $14=="0" {print}' Check.filtered.biSNP.tped | wc -l > A05.missing &
awk '$15=="0" && $16=="0" {print}' Check.filtered.biSNP.tped | wc -l > A06.missing &
awk '$17=="0" && $18=="0" {print}' Check.filtered.biSNP.tped | wc -l > A07.missing &
awk '$19=="0" && $20=="0" {print}' Check.filtered.biSNP.tped | wc -l > A08.missing &
awk '$21=="0" && $22=="0" {print}' Check.filtered.biSNP.tped | wc -l > A09.missing &
awk '$23=="0" && $24=="0" {print}' Check.filtered.biSNP.tped | wc -l > A10.missing &
awk '$25=="0" && $26=="0" {print}' Check.filtered.biSNP.tped | wc -l > A11.missing &
awk '$27=="0" && $28=="0" {print}' Check.filtered.biSNP.tped | wc -l > A12.missing &

# Extract the results from previously mapping
for i in $(seq -w 0 12); do
bcftools view -S sample.list ../../../vcf.clean/Chr$i.clean.vcf.gz | bgzip -@ 20 > Pepper4.Chr$i.vcf.gz &
done

# Convert to plink biallelic
for i in $(seq -w 0 12); do
plink --vcf Pepper4.Chr$i.vcf.gz \
    --recode 12 --transpose \
    --double-id \
    --out Pepper4.Chr$i &
done

# Missing count
cat Pepper4.Chr*.tped | awk '$5=="0" && $6=="0" {print}' | wc -l > Pepper4.A01.missing &
cat Pepper4.Chr*.tped | awk '$7=="0" && $8=="0" {print}' | wc -l > Pepper4.A02.missing &
cat Pepper4.Chr*.tped | awk '$9=="0" && $10=="0" {print}' | wc -l > Pepper4.A03.missing &
cat Pepper4.Chr*.tped | awk '$11=="0" && $12=="0" {print}' | wc -l > Pepper4.A04.missing &

```

#### SNP check to old S8 references

```bash
# Working directory
cd /data/zhaojiantao/pepper/vg/04SNP/SNP.check/notused

# Join the gvcf files
sentieon driver -t 100 \
    -r /data/zhaojiantao/pepper/References/S8/pepper_S8.fasta --algo GVCFtyper \
    Check.vcf.gz *rename*gz
# Hard filter
gatk VariantFiltration \
    -R /data/zhaojiantao/pepper/References/S8/pepper_S8.fasta \
    -V Check.vcf.gz \
    --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || SOR >3.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
    --filter-name "my_snp_filter" \
    -O Check.filter.vcf.gz
# Extract the hard filtered
bcftools view Check.filter.vcf.gz \
    --apply-filters PASS | \
    bgzip -@ 20 \
    > Check.filtered.vcf.gz
# Index
tabix -p vcf Check.filtered.vcf.gz
# Extract biallelic    
bcftools view Check.filtered.vcf.gz --types snps -m 2 -M 2 | bgzip -@ 20 > Check.filtered.biSNP.vcf.gz
tabix -p vcf Check.filtered.biSNP.vcf.gz


```

#### SNP check to vg reference for bwa mapping

```bash
# Working directory
cd /data/zhaojiantao/pepper/vg/04SNP/SNP.check

# Join the gvcf files
sentieon driver -t 100 \
    -r /data/zhaojiantao/pepper/vg/PanReference/Cann_Zhangshugang.chr.fa --algo GVCFtyper \
    Check_bwa.vcf.gz *bwa.rename.gvcf.gz
# Hard filter
gatk VariantFiltration \
    -R /data/zhaojiantao/pepper/vg/PanReference/Cann_Zhangshugang.chr.fa \
    -V Check_bwa.vcf.gz \
    --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || SOR >3.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
    --filter-name "my_snp_filter" \
    -O Check_bwa.filter.vcf.gz
# Extract the hard filtered
bcftools view Check_bwa.filter.vcf.gz \
    --apply-filters PASS | \
    bgzip -@ 20 \
    > Check_bwa.filtered.vcf.gz
# Index
tabix -p vcf Check_bwa.filtered.vcf.gz
# Extract biallelic    
bcftools view Check_bwa.filtered.vcf.gz --types snps -m 2 -M 2 | bgzip -@ 20 > Check_bwa.filtered.biSNP.vcf.gz
tabix -p vcf Check_bwa.filtered.biSNP.vcf.gz

```

#### SNP check to vg reference for vg mapping

```bash
# Working directory
cd /data/zhaojiantao/pepper/vg/04SNP/SNP.check

# Join the gvcf files
sentieon driver -t 100 \
    -r /data/zhaojiantao/pepper/vg/PanReference/Cann_Zhangshugang.chr.fa --algo GVCFtyper \
    Check_gam.vcf.gz *gam.rename.gvcf.gz
# Hard filter
gatk VariantFiltration \
    -R /data/zhaojiantao/pepper/vg/PanReference/Cann_Zhangshugang.chr.fa \
    -V Check_gam.vcf.gz \
    --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || SOR >3.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
    --filter-name "my_snp_filter" \
    -O Check_gam.filter.vcf.gz
# Extract the hard filtered
bcftools view Check_gam.filter.vcf.gz \
    --apply-filters PASS | \
    bgzip -@ 20 \
    > Check_gam.filtered.vcf.gz
# Index
tabix -p vcf Check_gam.filtered.vcf.gz
# Extract biallelic    
bcftools view Check_gam.filtered.vcf.gz --types snps -m 2 -M 2 | bgzip -@ 20 > Check_gam.filtered.biSNP.vcf.gz
tabix -p vcf Check_gam.filtered.biSNP.vcf.gz

```

### SNP LD prune

```bash
# Working directory
cd /data/zhaojiantao/pepper/vg/04SNP

# Add the SNP positions in the vcf file
for i in $(seq -w 0 12); do
zcat Chr00.SNP.vcf.gz | grep "#" > header
bgzip -d Chr$i.SNP.vcf.gz -@ 20 -c | grep -v "#" > Chr$i.geno
cut -f1,2 Chr$i.geno | awk '{print $1,$2,$1"_"$2}' | sed 's/ /\t/g' | \
    paste - Chr$i.geno | cut -f1-3,7- | cat header - | bgzip -@ 20 \
    > Chr$i.SNP.withid.vcf.gz
done

# Working directory
cd /data/zhaojiantao/pepper/vg/08SNPprune

# LD-prune
for i in $(seq -w 1 12); do
# LD prune
plink --vcf ../04SNP/Miss_0.3/Chr$i.SNP.vcf.gz \
    --indep-pairwise 50 10 0.2 \
    --allow-extra-chr --double-id \
    --out Chr$i
# Extract the prunned SNPs
plink --vcf ../04SNP/Miss_0.3/Chr$i.SNP.vcf.gz \
    --extract Chr$i.prune.in \
    --recode vcf --double-id --allow-extra-chr \
    --out Chr$i.prune
done

# Use single sample id and bgzip the file
grep "#" Chr00.prune.vcf | grep -v "#CHROM" | cat - <(zcat ../04SNP/Miss_0.3/Chr00.SNP.withid.vcf.gz | grep "#CHROM" ) > header

for i in $(seq -w 1 12); do
grep -v "#" Chr$i.prune.vcf | cat header - | bgzip -@ 20 > Chr$i.prune.vcf.gz &
done

# Remove temporary files
rm *log *nosex *prune.in Chr*.prune.vcf header

```


### Extract 4DTV SNPs and IQ-tree

```bash
# Working directory
cd /data/zhaojiantao/pepper/vg/064DTV

# Call 4dtv SNPs
4DTV_scan.pl -g ../PanReference/Cann_Zhangshugang.chr.fa \
    -f ../PanReference/Cann_Zhangshugang.gff3 \
    -o ffds.list

# Select the ids
for i in $(seq -w 1 12); do
grep Chr$i ffds.list | cut -f1,2 > Chr$i.ffds.list
done
grep -v Chr ffds.list | cut -f1,2 > Chr00.ffds.list

# Extract the 4DTV SNPs
for i in $(seq -w 1 12); do
bcftools view ../04SNP/Chr$i.SNP.vcf.gz \
    --regions-file Chr$i.ffds.list | \
    bgzip -@ 20 \
    > Chr$i.ffds.SNP.vcf.gz &
done

# Merge the file
bcftools concat Chr*.ffds.SNP.vcf.gz | bgzip -@ 20 > 500.ffds.SNP.vcf.gz
tabix -p vcf 500.ffds.SNP.vcf.gz

# Use tomato as outgroups
for i in ERR418102 ERR418107 ERR418103 ERR418106; do
    # Download SRA
    prefetch $i
    # Convert SRA to fastq
    fasterq-dump --split-files  $i/$i.sra
    # Remove adaptor and low-quality sequences
    trimmomatic PE ${i}_1.fastq ${i}_2.fastq \
            $i/$i.R1P.fq $i/$i.R1U.fq \
            $i/$i.R2P.fq $i/$i.R2U.fq \
            -threads 100 \
            ILLUMINACLIP:/data/zhaojiantao/tools/miniconda3/share/trimmomatic-0.39-2/adapters/TruSeq3-PE-2.fa:2:30:10:1:TRUE \
            SLIDINGWINDOW:4:20 LEADING:3 TRAILING:3 MINLEN:40
    # Remove temporary file
    rm $i/$i.sra $i/$i.R1U.fq $i/$i.R2U.fq ${i}_1.fastq ${i}_2.fastq 
    # vg mapping
    vg giraffe -t 100 \
        -Z m6a.giraffe.gbz \
        -m m6a.min \
        -d m6a.dist \
        -f $i/$i.R1P.fq \
        -f $i/$i.R2P.fq \
        > $i.gam
    # Generate bam
    vg surject -t 100 \
        -N $i \
        -R $i \
        -x m6a.giraffe.gbz \
        -b -i $i.gam \
        > $i.gam.bam
    # Sort the bam file
    perl sort_bam_header_chrID.pl $i $i.gam.bam ${i}_gam.bam
    # Index the bam
    samtools index -@ 100 ${i}_gam.bam        
    # Call variants
    # Collect read information to remove or mark duplicates
    sentieon driver -t 100 -i ${i}_gam.bam \
        --algo LocusCollector --fun score_info $i.score.gz
    # Remove or mark duplicates
    sentieon driver -t 100 -i ${i}_gam.bam \
        --algo Dedup --score_info $i.score.gz \
        --metrics $i.dedup_metrix.txt $i.deduped.bam
    # Call gvcf
    sentieon driver -t 100 \
        -r ../PanReference/Cann_Zhangshugang.chr.fa \
        -i $i.deduped.bam \
        --algo Haplotyper \
        --emit_mode gvcf \
        $i.gvcf.gz
    # Remove unnecessary files
    rm $i.score.gz* $i.dedup_metrix.txt $i.deduped.bam*
done

# Join the gvcf files
sentieon driver -t 100 \
    -r /data/zhaojiantao/pepper/vg/PanReference/Cann_Zhangshugang.chr.fa --algo GVCFtyper --emit_mode all \
    ERR418102.vcf.gz ERR418102.gvcf.gz

# Hard filter
gatk VariantFiltration \
    -R /data/zhaojiantao/pepper/vg/PanReference/Cann_Zhangshugang.chr.fa \
    -V Tomato.vcf.gz \
    --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || SOR >3.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
    --filter-name "my_snp_filter" \
    -O Tomato.filter.vcf.gz
# Extract the hard filtered
bcftools view Tomato.filter.vcf.gz \
    --apply-filters PASS | \
    bgzip -@ 20 \
    > Tomato.filtered.vcf.gz
# Index
tabix -p vcf Tomato.filtered.vcf.gz
# Extract biallelic    
bcftools view Tomato.filtered.vcf.gz --types snps -m 2 -M 2 | bgzip -@ 20 > Tomato.filtered.biSNP.vcf.gz
bcftools view Tomato.filter.vcf.gz --types snps -m 2 -M 2 | bgzip -@ 20 > Tomato.filter.biSNP.vcf.gz
bcftools view 504.4DTV.vcf.gz --types snps -m 2 -M 2 | bgzip -@ 20 > 504.4DTV.biSNP.vcf.gz

tabix -p vcf Tomato.filtered.biSNP.vcf.gz


# Filter the 4DTV with missing > 0.3 and maf < 0.05

## Missing<0.3
bcftools view -i 'F_MISSING<0.3' 500.ffds.SNP.vcf.gz | bgzip -@ 20 > 500.ffds.filter.SNP.vcf.gz
tabix -p vcf 500.ffds.filter.SNP.vcf.gz
bcftools query -f '%CHROM %POS\n' 500.ffds.filter.SNP.vcf.gz | sed 's/ /\t/g' > 500.ffds.filter.SNP.info

# Extract the common SNPs from tomato
bcftools view Tomato.vcf.gz \
    --regions-file 500.ffds.filter.SNP.info | \
    bgzip -@ 20 \
    > Tomato.ffds.filter.SNP.vcf.gz
tabix -p vcf Tomato.ffds.filtered.SNP.vcf.gz

bcftools view Tomato.ffds.filtered.SNP.vcf.gz \
    --regions-file comm.info | \
    bgzip -@ 20 \
    > Tomato.ffds.filtered2.SNP.vcf.gz
tabix -p vcf 500.ffds.filtered2.SNP.vcf.gz

# Extract the tomato info
bcftools query -f '%CHROM %POS\n' Tomato.ffds.filtered.SNP.vcf.gz | sed 's/ /\t/g' > Tomato.ffds.filtered.SNP.info
bcftools query -f '%CHROM %POS\n' Tomato.filtered.biSNP.vcf.gz | sed 's/ /\t/g' > Tomato.filtered.biSNP.SNP.info
bcftools query -f '%CHROM %POS\n' 500.ffds.filtered.SNP.vcf.gz | sed 's/ /\t/g' > 500.ffds.filtered.SNP.info

Tomato.filtered.biSNP.vcf.gz

# Find the overlapped SNPs
comm -12 <(sort 500.ffds.filter.SNP.info) <(sort Tomato.ffds.filter.SNP.info) > comm.SNP.info
comm -12 <(sort 500.ffds.filtered.SNP.info) <(sort Tomato.filtered.biSNP.SNP.info) > comm.SNP.info

# Extract the final results
bcftools view Tomato.ffds.filter.SNP.vcf.gz \
    --regions-file comm.SNP.info | \
    bgzip -@ 20 \
    > Tomato.ffds.filtered.SNP.vcf.gz

bcftools view Tomato.filtered.biSNP.vcf.gz \
    --regions-file comm.SNP.info | \
    bgzip -@ 20 \
    > Tomato.ffds.filtered.biSNP.vcf.gz
tabix -p vcf Tomato.ffds.filtered.biSNP.vcf.gz

bcftools view 500.ffds.filter.SNP.vcf.gz \
    --regions-file comm.SNP.info | \
    bgzip -@ 20 \
    > 500.ffds.filtered.SNP.vcf.gz

bcftools view 500.ffds.filtered.SNP.vcf.gz \
    --regions-file comm.SNP.info | \
    bgzip -@ 20 \
    > 500.ffds.filtered.biSNP.vcf.gz
tabix -p vcf 500.ffds.filtered.biSNP.vcf.gz

# Merge tomato and pepper vcfs
bcftools merge Tomato.ffds.filtered.SNP.vcf.gz 500.ffds.filtered.SNP.vcf.gz -O z -o 504.4DTV.filtered.vcf.gz
bcftools merge Tomato.ffds.filtered.biSNP.vcf.gz 500.ffds.filtered.biSNP.vcf.gz -O z -o 504.4DTV.filtered.biSNP.vcf.gz


# Number of SNPs with F_MISSING<0.3
bcftools view -i 'F_MISSING<0.2' 504.4DTV.filtered.vcf.gz | bgzip -@ 20 > 504.4DTV.filtered2.vcf.gz
bcftools view -i 'F_MISSING<0.1' 504.4DTV.filtered.vcf.gz | bgzip -@ 20 > 504.4DTV.filtered1.vcf.gz

# LD-prune
plink --vcf 504.4DTV.Miss0.1.vcf.gz \
    --indep-pairwise 50 10 0.4 \
    --allow-extra-chr --double-id \
    --out 504.4DTV.Miss0.1.LD0.4

# Extract the prunned SNPs
plink --vcf 504.4DTV.Miss0.1.vcf.gz \
    --extract 504.4DTV.Miss0.1.LD0.4.prune.in \
    --recode vcf --double-id --allow-extra-chr \
    --out 504.4DTV.Miss0.1.LD0.4

zcat 504.4DTV.Miss0.1.vcf.gz | grep "#CHROM" > AA
grep -v "#" 504.4DTV.Miss0.1.LD0.4.vcf > BB
grep "#" 504.4DTV.Miss0.1.LD0.4.vcf |grep -v "#CHROM" | cat - AA BB | bgzip -@ 20 > 504.4DTV.Miss0.1.LD0.4.gz

# Convert vcf to phylip
python /data/zhaojiantao/tools/PythonScripts/vcf2phylip.py -i 504.4DTV.Miss0.1.LD0.4.gz
python /data/zhaojiantao/tools/PythonScripts/vcf2phylip.py -i 504.4DTV.filtered.biSNP.vcf.gz

# Run IQ-tree
iqtree -s 504.4DTV.Miss0.1.LD0.4.gz.phy -nt 100 -m MFP
iqtree -s 504.4DTV.filtered.biSNP.min4.phy -nt 100 -m MFP

# Save the tree for phylogenetic tree plot
cp 504.4DTV.Miss0.1.LD0.4.gz.phy.treefile 504.4DTV.vg.phy.txt

# LD-prune
plink --vcf 504.4DTV.biSNP.vcf.gz \
    --indep-pairwise 50 10 0.4 \
    --allow-extra-chr --double-id \
    --out 504.4DTV.biSNP.Miss0.3.LD0.2

# Extract the prunned SNPs
plink --vcf 504.4DTV.biSNP.vcf.gz \
    --extract 504.4DTV.biSNP.Miss0.3.LD0.2.prune.in \
    --recode vcf --double-id --allow-extra-chr \
    --out 504.4DTV.biSNP.Miss0.3.LD0.2

zcat 504.4DTV.biSNP.vcf.gz | grep "#CHROM" > AA
grep -v "#" 504.4DTV.biSNP.Miss0.3.LD0.2.vcf > BB
grep "#" 504.4DTV.biSNP.Miss0.3.LD0.2.vcf |grep -v "#CHROM" | cat - AA BB | bgzip -@ 20 > 504.4DTV.biSNP.Miss0.3.LD0.2.vcf.gz

# Convert vcf to phylip
python /data/zhaojiantao/tools/PythonScripts/vcf2phylip.py -i 504.4DTV.biSNP.Miss0.3.LD0.2.vcf.gz

# Run IQ-tree
iqtree -s 504.4DTV.biSNP.Miss0.3.LD0.2.phy -nt 100 -m MFP


# LD-prune
plink --vcf 504.4DTV.Miss0.1.vcf.gz \
    --indep-pairwise 50 10 0.5 \
    --allow-extra-chr --double-id \
    --out 504.4DTV.Miss0.1.LD0.5

# Extract the prunned SNPs
plink --vcf 504.4DTV.Miss0.1.vcf.gz \
    --extract 504.4DTV.Miss0.1.LD0.5.prune.in \
    --recode vcf --double-id --allow-extra-chr \
    --out 504.4DTV.Miss0.1.LD0.5

zcat 504.4DTV.Miss0.1.vcf.gz | grep "#CHROM" > AA
grep -v "#" 504.4DTV.Miss0.1.LD0.5.vcf> BB
grep "#" 504.4DTV.Miss0.1.LD0.5.vcf |grep -v "#CHROM" | cat - AA BB | bgzip -@ 20 > 504.4DTV.Miss0.1.LD0.5.vcf.gz

# Convert vcf to phylip
python /data/zhaojiantao/tools/PythonScripts/vcf2phylip.py -i 504.4DTV.Miss0.1.LD0.5.vcf.gz

# Run IQ-tree
iqtree -s 504.4DTV.Miss0.1.LD0.5.min4.phy -nt 100 -m MFP
# Save the tree for phylogenetic tree plot
cp  504.4DTV.LD.min4.phy.treefile 504.4DTV.vg.phy.txt

iqtree -s 504.4DTV.Miss1.LD4.min4.phy -nt 70 -m MFP
iqtree -s 504.4DTV.LD0.min4.phy -nt AUTO -m MFP


# LD-prune
plink --vcf 501.4DTV.filtered.vcf.gz \
    --indep-pairwise 50 10 0.2 \
    --allow-extra-chr --double-id \
    --out 501.4DTV.filtered

# Extract the prunned SNPs
plink --vcf 501.4DTV.filtered.vcf.gz \
    --extract 501.4DTV.filtered.prune.in \
    --recode vcf --double-id --allow-extra-chr \
    --out 501.4DTV.LD

zcat 501.4DTV.filtered.vcf.gz | grep "#CHROM" > AA
grep -v "#" 501.4DTV.LD.vcf > BB
grep "#" 501.4DTV.LD.vcf | grep -v "#CHROM" | cat - AA BB | bgzip -@ 20 > 501.4DTV.LD.vcf.gz

# Convert vcf to phylip
python /data/zhaojiantao/tools/PythonScripts/vcf2phylip.py -i 501.4DTV.LD.vcf.gz

# Run IQ-tree
iqtree -s 501.4DTV.LD.min4.phy -nt 60 -m MFP
# Save the tree for phylogenetic tree plot
cp  501.4DTV.LD.min4.phy.treefile 501.4DTV.vg.phy.txt

iqtree -s 501.4DTV.Miss1.LD4.min4.phy -nt 70 -m MFP
iqtree -s 501.4DTV.LD0.min4.phy -nt AUTO -m MFP


```

### Newhyvbrid

```bash
# Working directory
cd /data/zhaojiantao/pepper/vg/064DTV/newhybrid

# Test for potential hybrids between C. chacoense and C. baccatum var. baccatum
cat ../../Pop.list/CHA.list ../../Pop.list/BAB.list > CHA.BAC.list
# Add additional 17 BAP accessions

# Extract the non-missing SNPs
bcftools view -i 'F_MISSING=0' ../504.4DTV.LD1.vcf.gz | bgzip -@ 20 > 504.4DTV.LD1.Nomiss.vcf.gz
# Extract the CHA.BAC list
bcftools view -S CHA.BAC.list 504.4DTV.LD1.Nomiss.vcf.gz | bgzip -@ 20 > CHA.BAC.vcf.gz

# Convert to plink
plink --vcf CHA.BAC.vcf.gz \
    --recode 12 --allow-extra-chr \
    --output-missing-genotype 9 \
    --out CHA.BAC

# Combine every two columns
cut -f7- -d ' ' CHA.BAC.ped | sed 's! \([^ ]\+\)\( \|$\)!\1 !g' > CHA.BAC.ped-1

# Add the sample number
seq 43 | paste - CHA.BAC.ped-1 > CHA.BAC.ped-2

# Add the SNP information
sed 's/Ca.S8.//g' CHA.BAC.map | awk '{print $2}' | sed '1i\LocusNames' > CHA.BAC.snps

echo "NumIndivs 43
NumLoci 1046
Digits 1
Format Lumped" \
> header

deal_table.pl -transpose CHA.BAC.snps | \
    cat header - CHA.BAC.ped-2 | sed 's/ /\t/g' \
    > CHA.BAC.NewHybrids.txt

# Check the data format
head -10 CHA.BAC.NewHybrids.txt | cut -f1-5

# Test with regular SNPs


# Newhybrid has to be run on either window or mac
E:\Project\Pepper\Pepper\newhybrids\newhybrids\bin\PC\newhybrids.exe -d E:\Project\Pepper\Pepper\newhybrids\newhybrids\test_data\test\TestAFLP.txt -c E:\Project\Pepper\Pepper\newhybrids\newhybrids\test_data\test\TwoGensGtypFreq.txt --burn-in 100 --num-sweeps 100 --no-gui

```

### HyDe

#### CHA and BAC Miss0.1

```bash
# Working directory
cd /data/zhaojiantao/pepper/vg/064DTV/HyDe/CHA.BAC

# Extract the CHA + BAC accessions
bcftools view -S CHA.BAC.list ../../504.4DTV.Miss0.1.vcf.gz | bgzip -@ 20 > CHA.BAC.vcf.gz

# Convert to phy
python /data/zhaojiantao/tools/PythonScripts/vcf2phylip.py -i CHA.BAC.vcf.gz
sed '1d' CHA.BAC.min4.phy > CHA.BAC.min4.txt

# Run a full hybridization detection analysis
run_hyde_mp.py -i CHA.BAC.min4.txt -m map.txt -o C.pubescens -n 136 -t 4 -s 145492 -j 100

# Test all individuals for the hybrid populations specified
individual_hyde_mp.py -i CHA.BAC.min4.txt -m map.txt -tr hyde-out.txt -o C.pubescens -n 136 -t 4 -s 145492 -j 100

```


#### Annuum Miss0.1

```bash
# Working directory
cd /data/zhaojiantao/pepper/vg/064DTV/HyDe/Annuum

# Extract the Annuum accessions
bcftools view -S Annuum.list ../../504.4DTV.Miss0.1.vcf.gz | bgzip -@ 20 > Annuum.vcf.gz

# Convert to phy
python /data/zhaojiantao/tools/PythonScripts/vcf2phylip.py -i Annuum.vcf.gz
sed '1d' Annuum.min4.phy > Annuum.min4.txt

# Run a full hybridization detection analysis
run_hyde_mp.py -i Annuum.min4.txt -m map.txt -o CHA -n 327 -t 4 -s 145492 -j 100 -tr triples.txt

# Test all individuals for the hybrid populations specified
individual_hyde_mp.py -i Annuum.min4.txt -m map.txt -tr hyde-out.txt -o CHA -n 327 -t 4 -s 145492 -j 100
mv mv hyde-ind.txt hyde-ind.all.txt
# Only for the significant
individual_hyde_mp.py -i Annuum.min4.txt -m map.txt -tr hyde-out-filtered.txt -o CHA -n 327 -t 4 -s 145492 -j 100

# Bootstrap
bootstrap_hyde_mp.py -i Annuum.min4.txt -m map.txt -tr hyde-out-filtered.txt -o CHA -n 327 -t 4 -s 145492 -j 100

```

### Population structure

```bash
# Working directory
cd /data/zhaojiantao/pepper/vg/07Structure

# Requies python 2.7
conda activate fastStructure
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/data/zhaojiantao/tools/miniconda3/pkgs/gsl-2.6-he838d99_2/lib

# Make bed for the ped/map plink file
plink --vcf ../064DTV/504.4DTV.LD1.vcf.gz --make-bed --double-id --allow-extra-chr --out 504.4DTV

# fastStructure results at K = 7
structure.py -K 7 --input=504.4DTV --output=K7 --full --seed=100

# fastStructure
for i in $(seq 1 20); do 
    python2 /data/zhaojiantao/tools/fastStructure/structure.py -K $i --input=504.4DTV --output=output --full --seed=100
done

# Choose the optimal K (K = 7)
python2 /data/zhaojiantao/tools/fastStructure/chooseK.py --input=output

```

### Fst

```bash
# Working directory
cd /data/zhaojiantao/pepper/vg/09Fst

# Calculate Fst between PUB and the other subgroups
for i in PUB; do
    for j in CHA BAB BAP FRU CHN GLA ANN; do
        for chr in $(seq -w 1 12); do
        vcftools --gzvcf ../08SNPprune/Chr$chr.prune.vcf.gz \
            --weir-fst-pop ../Pop.list/$i.list \
            --weir-fst-pop ../Pop.list/$j.list \
            --fst-window-size 100000 \
            --fst-window-step 10000 \
            --out $i.$j.Chr$chr &
        done
    done
done

# Calculate Fst between CHA and the other subgroups
for i in CHA; do
    for j in BAB BAP FRU CHN GLA ANN; do
        for chr in $(seq -w 1 12); do
        vcftools --gzvcf ../08SNPprune/Chr$chr.prune.vcf.gz \
            --weir-fst-pop ../Pop.list/$i.list \
            --weir-fst-pop ../Pop.list/$j.list \
            --fst-window-size 100000 \
            --fst-window-step 10000 \
            --out $i.$j.Chr$chr &
        done
    done
done

# Calculate Fst between BAB and the other subgroups
for i in BAB; do
    for j in BAP FRU CHN GLA ANN; do
        for chr in $(seq -w 1 12); do
        vcftools --gzvcf ../08SNPprune/Chr$chr.prune.vcf.gz \
            --weir-fst-pop ../Pop.list/$i.list \
            --weir-fst-pop ../Pop.list/$j.list \
            --fst-window-size 100000 \
            --fst-window-step 10000 \
            --out $i.$j.Chr$chr &
        done
    done
done

# Calculate Fst between BAP and the other subgroups
for i in BAP; do
    for j in FRU CHN GLA ANN; do
        for chr in $(seq -w 1 12); do
        vcftools --gzvcf ../08SNPprune/Chr$chr.prune.vcf.gz \
            --weir-fst-pop ../Pop.list/$i.list \
            --weir-fst-pop ../Pop.list/$j.list \
            --fst-window-size 100000 \
            --fst-window-step 10000 \
            --out $i.$j.Chr$chr &
        done
    done
done

# Calculate Fst between FRU and the other subgroups
for i in FRU; do
    for j in CHN GLA ANN; do
        for chr in $(seq -w 1 12); do
        vcftools --gzvcf ../08SNPprune/Chr$chr.prune.vcf.gz \
            --weir-fst-pop ../Pop.list/$i.list \
            --weir-fst-pop ../Pop.list/$j.list \
            --fst-window-size 100000 \
            --fst-window-step 10000 \
            --out $i.$j.Chr$chr &
        done
    done
done

# Calculate Fst between CHN and the other subgroups
for i in CHN; do
    for j in GLA ANN; do
        for chr in $(seq -w 1 12); do
        vcftools --gzvcf ../08SNPprune/Chr$chr.prune.vcf.gz \
            --weir-fst-pop ../Pop.list/$i.list \
            --weir-fst-pop ../Pop.list/$j.list \
            --fst-window-size 100000 \
            --fst-window-step 10000 \
            --out $i.$j.Chr$chr &
        done
    done
done

# Calculate Fst between CHN and the other subgroups
for i in GLA; do
    for j in ANN; do
        for chr in $(seq -w 1 12); do
        vcftools --gzvcf ../08SNPprune/Chr$chr.prune.vcf.gz \
            --weir-fst-pop ../Pop.list/$i.list \
            --weir-fst-pop ../Pop.list/$j.list \
            --fst-window-size 100000 \
            --fst-window-step 10000 \
            --out $i.$j.Chr$chr &
        done
    done
done

# Combine the results
for i in $(ll *weir.fst | awk '{print $9}' | sed 's/.Chr/\t/g' | cut -f1 | sort | uniq); do 
    cat $i.Chr*.windowed.weir.fst | \
        grep -v CHROM | \
        awk '{print $7=$1"_"$2"_"$3,$1,$2,$3,$8=($2+$3)/2,$4,$5,$6}' | \
        awk '$8 > 0 {print $0}' | \
        awk '{gsub(/Chr0/, "", $2)} 1' | \
        awk '{gsub(/Chr/, "", $2)} 1' | \
        sed '1i\ID CHROM BIN_START BIN_END BIN_Middle N_VARIANTS WEIGHTED_FST MEAN_FST' | \
        sed 's/ /\t/g' \
        > $i.fst.txt
done

# Remove chromosomal level files
rm *.windowed.weir.fst *log 

# Calculate the mean Fst value
for i in $(ll *fst.txt.gz | awk '{print $9}' | sed 's/.fst.txt.gz//g'); do
    echo $i; zcat $i.fst.txt.gz | grep -v MEAN | awk '{sum+=$8} END {print sum/NR}'
done

```

#### Fst permutation for BAB.BAP

```bash
# Working directory
cd /data/zhaojiantao/pepper/vg/09Fst/permutations

# Prepare the list for permutation test
cat ../../Pop.list/BAB.list ../../Pop.list/BAP.list > BAB.BAP.list

# Randomize the samples 
for i in $(seq -w 1 100); do
    sort --random-sort BAB.BAP.list | head -9 > BAB.BAP.Random$i.BAB.list
    grep -vf BAB.BAP.Random$i.BAB.list BAB.BAP.list > BAB.BAP.Random$i.BAP.list
done

# Run the permutation test with 100 replicates
for i in $(seq -w 1 100); do
    for j in $(seq -w 1 12); do
    vcftools --gzvcf ../../08SNPprune/Chr$j.prune.vcf.gz \
        --weir-fst-pop BAB.BAP.Random$i.BAB.list \
        --weir-fst-pop BAB.BAP.Random$i.BAP.list \
        --fst-window-size 100000 \
        --fst-window-step 10000 \
        --out BAB.BAP.Random$i.Chr$j &
    done
done

# Calculate the mean value
for i in $(seq -w 1 100); do
    cat BAB.BAP.Random$i.Chr*.windowed.weir.fst | \
        grep -v CHROM | awk '$6>=0 {print}' | awk '{sum+=$6} END {print sum/NR}' \
        > BAB.BAP.Random$i.meanFst
done

# Combine the permutation result
cat *Random*.meanFst > BAB.BAP.Random.meanFst

# Remove tempory file
rm *log *weir.fst BAB.BAP.Random*.list
rm BAB.BAP.Random*.meanFst

# Generate scripts for significant t test in R
echo "
data <- read.table("BAB.BAP.Random.meanFst", header=FALSE)
AA <- t.test(x=data, mu=0.31, conf.level=0.95)
BAPrs <- capture.output(print(AA))
writeLines(BAPrs, con = file("t.test.txt")) " | \
sed 's/BAB.BAP.Random.meanFst/"BAB.BAP.Random.meanFst"/g' |\
sed 's/t.test.txt/"BAB.BAP.t.test.txt"/g' \
> BAB.BAP.t.test.R

# Significant t test in R
Rscript BAB.BAP.t.test.R

```

#### Fst permutation for PUB with the other groups

```bash
# Working directory
cd /data/zhaojiantao/pepper/vg/09Fst/permutations

# Run Fst permutation
for group in CHA BAB BAP FRU CHN GLA ANN; do
    for i in $(seq -w 1 100); do
        for j in $(seq -w 1 12); do
    # Prepare the list for permutation test
    cat ../../Pop.list/PUB.list ../../Pop.list/$group.list > PUB.$group.list
    # Randomize the samples 
    sort --random-sort PUB.$group.list | head -38 > PUB.$group.Random$i.PUB.list
    grep -vf PUB.$group.Random$i.PUB.list PUB.$group.list > PUB.$group.Random$i.$group.list
    
    # Run the permutation test with 100 replicates
    vcftools --gzvcf ../../08SNPprune/Chr$j.prune.vcf.gz \
        --weir-fst-pop PUB.$group.Random$i.PUB.list \
        --weir-fst-pop PUB.$group.Random$i.$group.list \
        --fst-window-size 100000 \
        --fst-window-step 10000 \
        --out PUB.$group.Random$i.Chr$j
    # Calculate the mean value
    cat PUB.$group.Random$i.Chr*.windowed.weir.fst | \
        grep -v CHROM | awk '$6>=0 {print}' | awk '{sum+=$6} END {print sum/NR}' \
        > PUB.$group.Random$i.meanFst
        done
    done
done

# Combine the permutation result
for group in CHA BAB BAP FRU CHN GLA ANN; do
cat PUB.*Random*.meanFst > PUB.$group.Random.meanFst
done

# Remove tempory file
rm *log *weir.fst PUB.$group.Random*.list
rm PUB.$group.Random*.meanFst

# Generate scripts for significant t test in R
echo "
data <- read.table("PUB.$group.Random.meanFst", header=FALSE)
AA <- t.test(x=data, mu=0.31, conf.level=0.95, alternative = "two.sided", )
$grouprs <- capture.output(print(AA))
writeLines($grouprs, con = file("t.test.txt")) " | \
sed 's/PUB.$group.Random.meanFst/"PUB.$group.Random.meanFst"/g' |\
sed 's/t.test.txt/"PUB.$group.t.test.txt"/g' \
> PUB.$group.t.test.R

# Significant t test in R
Rscript PUB.$group.t.test.R

```

#### Top 1% Fst

```bash
# Working directory
cd /data/zhaojiantao/pepper/vg/09Fst

for i in $(ls *fst.txt | sed 's/.txt//g' | grep -v CHA | grep -v BAB | grep -v GLA); do 
cut -f1,8 $i.txt | grep -v ID | sed 's/_/\t/g' | awk '{print $1"_"$2,$1,$2,$3,$4}' | sed 's/ /\t/g' > $i-2.txt
done

# Top 1% signal
wc -l *fst-2.txt
sort -k5gr BAP.ANN.fst-2.txt | head -2854 > BAP.ANN.fst.top1.txt
sort -k5gr BAP.CHN.fst-2.txt | head -2896 > BAP.CHN.fst.top1.txt
sort -k5gr BAP.FRU.fst-2.txt | head -2896 > BAP.FRU.fst.top1.txt
sort -k5gr CHN.ANN.fst-2.txt | head -2893 > CHN.ANN.fst.top1.txt
sort -k5gr FRU.ANN.fst-2.txt | head -2894 > FRU.ANN.fst.top1.txt
sort -k5gr FRU.CHN.fst-2.txt | head -2890 > FRU.CHN.fst.top1.txt
sort -k5gr PUB.ANN.fst-2.txt | head -2805 > PUB.ANN.fst.top1.txt
sort -k5gr PUB.BAP.fst-2.txt | head -2865 > PUB.BAP.fst.top1.txt
sort -k5gr PUB.CHN.fst-2.txt | head -2883 > PUB.CHN.fst.top1.txt
sort -k5gr PUB.FRU.fst-2.txt | head -2882 > PUB.FRU.fst.top1.txt

# Sort the file and add header in local excel
for i in $(ls *top1.txt | sed 's/.txt//g'); do perl merge.pl $i.txt > $i.merged.txt ; done

# Remove temporary files
rm rm *fst-2.txt *top1

# Calculate the total differentiated size for the five domesticated species
for i in $(ls *top1.merged.txt | grep -v BAB | grep -v CHA | grep -v GLA); do
    sed '1d' $i | awk '{sum += $3-$2+1 } END { print sum }'
done

# Zhangshugang bed
awk '$3=="gene" {print $1,$4,$5,$9}' ../../References/S8/Zhangshugang.gff | \
    sed 's/;/ /g' | cut -f1-4 -d ' ' | sed 's/ID=//g' | sed 's/ /\t/g' \
    > Zhangshugang.bed

# Find the overlapped genes
for i in $(ls *top1.merged.txt | sed 's/.merged.txt//g' | grep -v BAB | grep -v CHA | grep -v GLA); do
    bedtools intersect -a Zhangshugang.bed \
        -b <(sed '1d' $i.merged.txt) | cut -f4 \
        > $i.overlapped.genes
done

```

### LD estimation

```bash
# Working directory
cd /data/zhaojiantao/pepper/vg/10LD

# Use all SNPs
# Use PopLDdecay
for i in $(seq -w 1 12); do
    for j in PUB CHA BAB BAP FRU CHN GLA ANN; do
        PopLDdecay -InVCF ../04SNP/Miss_0.3/Chr$i.SNP.vcf.gz -SubPop ../Pop.list/$j.list \
            -MaxDist 2000 -MAF 0.05 -OutStat $j.Chr$i &
    done &
done

# Plot the LD
for i in PUB CHA BAB BAP FRU CHN GLA ANN; do
    for j in $(seq -w 1 12); do
        /data/zhaojiantao/tools/PopLDdecay-master/bin/Plot_OnePop.pl -inFile $i.Chr$j.stat.gz -output $i.Chr$j &
    done &
done

# Add chromosome information to each stat files
for i in PUB CHA BAB BAP FRU CHN GLA ANN; do
    zcat $i.Chr01.bin.gz | awk '{print $1,$2,1}' > $i.Chr01.bin.gz-2 ;
    zcat $i.Chr02.bin.gz | awk '{print $1,$2,2}' > $i.Chr02.bin.gz-2 ;
    zcat $i.Chr03.bin.gz | awk '{print $1,$2,3}' > $i.Chr03.bin.gz-2 ;
    zcat $i.Chr04.bin.gz | awk '{print $1,$2,4}' > $i.Chr04.bin.gz-2 ;
    zcat $i.Chr05.bin.gz | awk '{print $1,$2,5}' > $i.Chr05.bin.gz-2 ;
    zcat $i.Chr06.bin.gz | awk '{print $1,$2,6}' > $i.Chr06.bin.gz-2 ;
    zcat $i.Chr07.bin.gz | awk '{print $1,$2,7}' > $i.Chr07.bin.gz-2 ;
    zcat $i.Chr08.bin.gz | awk '{print $1,$2,8}' > $i.Chr08.bin.gz-2 ;
    zcat $i.Chr09.bin.gz | awk '{print $1,$2,9}' > $i.Chr09.bin.gz-2 ;
    zcat $i.Chr10.bin.gz | awk '{print $1,$2,10}' > $i.Chr10.bin.gz-2 ;
    zcat $i.Chr11.bin.gz | awk '{print $1,$2,11}' > $i.Chr11.bin.gz-2 ;
    zcat $i.Chr12.bin.gz | awk '{print $1,$2,12}' > $i.Chr12.bin.gz-2 ;
done

# Combine the results for each subgroup
for i in PUB CHA BAB BAP FRU CHN GLA ANN; do
cat $i.Chr*.bin.gz-2 | grep -v "#" > $i.ld
done

# Calculate the average LD across 12 chromosomes
for i in PUB CHA BAB BAP FRU CHN GLA ANN; do
    paste $i.Chr*.bin.gz-2 | sed '1d' | \
        awk '{print $1,($2+$5+$8+$11+$14+$17+$20+$23+$26+$29+$32+$35)/12}' \
        > $i.average.ld
done

# Combine the results
paste *average.ld | awk '{print $1,$2,$4,$6,$8,$10,$12,$14,$16,$18}' | \
    sed '1i\Dist A06 A03 A04 A02 A08 A07 A05 A01' | sed 's/ /\t/g' > Group.average.ld.txt


# Remove unnecessary files
rm *png *stat.gz *pdf *bin.gz *bin.gz-2

# Find the ld at LD =0.3
for i in PUB CHA BAB BAP FRU CHN GLA ANN; do
    echo $i ;
    cat $i.average.ld | awk '$2>0.2999 && $2<0.3001 {print $0}' 
done

# Find the minimum ld
for i in Ca.annuum Ca.bac.baccatum Ca.bac.pendulum Ca.chacoense Ca.chinense Ca.frutescens Ca.glabriusculum Ca.pubescens; do 
    cat $i.ld | awk '$3==1 {print $0}' | sort -nk2 | head -1 | cut -f2 -d ' ' ; 
done

```

### LD for Annuum clade

```bash
# Working directory
cd /data/zhaojiantao/pepper/vg/10LD

# Use PopLDdecay
for i in $(seq -w 1 12); do
    for j in FRU.high FRU.low FRU.mid CHN.high CHN.low CHN.mid ANN.high ANN.low ANN.mid ANN.30 ANN.60 CHN.30 CHN.60 FRU.30 FRU.60; do
        PopLDdecay -InVCF ../04SNP/Miss_0.3/Chr$i.SNP.vcf.gz -SubPop $j.list \
            -MaxDist 2000 -MAF 0.05 -OutStat $j.Chr$i &
    done &
done

# Plot the LD
for i in FRU.high FRU.low FRU.mid CHN.high CHN.low CHN.mid ANN.high ANN.low ANN.mid ANN.30 ANN.60 CHN.30 CHN.60 FRU.30 FRU.60; do
    for j in $(seq -w 1 12); do
        /data/zhaojiantao/tools/PopLDdecay-master/bin/Plot_OnePop.pl -inFile $i.Chr$j.stat.gz -output $i.Chr$j &
    done &
done

# Add chromosome information to each stat files
for i in FRU.high FRU.low FRU.mid CHN.high CHN.low CHN.mid ANN.high ANN.low ANN.mid ANN.30 ANN.60 CHN.30 CHN.60 FRU.30 FRU.60; do
    zcat $i.Chr01.bin.gz | awk '{print $1,$2,1}' > $i.Chr01.bin.gz-2 ;
    zcat $i.Chr02.bin.gz | awk '{print $1,$2,2}' > $i.Chr02.bin.gz-2 ;
    zcat $i.Chr03.bin.gz | awk '{print $1,$2,3}' > $i.Chr03.bin.gz-2 ;
    zcat $i.Chr04.bin.gz | awk '{print $1,$2,4}' > $i.Chr04.bin.gz-2 ;
    zcat $i.Chr05.bin.gz | awk '{print $1,$2,5}' > $i.Chr05.bin.gz-2 ;
    zcat $i.Chr06.bin.gz | awk '{print $1,$2,6}' > $i.Chr06.bin.gz-2 ;
    zcat $i.Chr07.bin.gz | awk '{print $1,$2,7}' > $i.Chr07.bin.gz-2 ;
    zcat $i.Chr08.bin.gz | awk '{print $1,$2,8}' > $i.Chr08.bin.gz-2 ;
    zcat $i.Chr09.bin.gz | awk '{print $1,$2,9}' > $i.Chr09.bin.gz-2 ;
    zcat $i.Chr10.bin.gz | awk '{print $1,$2,10}' > $i.Chr10.bin.gz-2 ;
    zcat $i.Chr11.bin.gz | awk '{print $1,$2,11}' > $i.Chr11.bin.gz-2 ;
    zcat $i.Chr12.bin.gz | awk '{print $1,$2,12}' > $i.Chr12.bin.gz-2 ;
done

# Combine the results for each subgroup
for i in FRU.high FRU.low FRU.mid CHN.high CHN.low CHN.mid ANN.high ANN.low ANN.mid ANN.30 ANN.60 CHN.30 CHN.60 FRU.30 FRU.60; do
cat $i.Chr*.bin.gz-2 | grep -v "#" > $i.ld
done

# Calculate the average LD across 12 chromosomes
for i in FRU.high FRU.low FRU.mid CHN.high CHN.low CHN.mid ANN.high ANN.low ANN.mid ANN.30 ANN.60 CHN.30 CHN.60 FRU.30 FRU.60; do
    paste $i.Chr*.bin.gz-2 | sed '1d' | \
        awk '{print $1,($2+$5+$8+$11+$14+$17+$20+$23+$26+$29+$32+$35)/12,"Group"}' | sed "s/Group/$i/g" \
        > $i.average.ld
done

cat *average.ld | awk '{gsub(/\./, " ", $3); print}' | cat - *average.ld-1 | sed '1i\Bin R2 Group1 Group2' | sed 's/ /\t/g' > ANN.ld.txt

```


### Pi

```bash
# Working directory
cd /data/zhaojiantao/pepper/vg/11Pi

# Estimate Pi for BAB accessions estimation
# Calculate per-site nucleotide diversity
for i in PUB CHA BAB BAP FRU CHN GLA ANN; do
    for j in $(seq -w 1 12); do
        vcftools --gzvcf ../04SNP/Miss_0.3/Chr$j.SNP.vcf.gz --keep ../Pop.list/$i.list --maf 0.05 --site-pi --out $i.Chr$j &
    done &
done

# Calculate windowed pi 
for i in PUB CHA BAB BAP FRU CHN GLA ANN; do
    for j in $(seq -w 1 12); do
        perl PI_slide_win.pl -w 100000 -s 10000 $i.Chr$j.sites.pi &
    done &
done

# Combine the results of each chromosome
for i in PUB CHA BAB BAP FRU CHN GLA ANN; do
    cat $i.Chr*.w100000.s10000.txt | \
        grep -v "#Chr" | \
        awk '{print $1"_"$3"_"$4,$1,$3,$4,($3+$4)/2,$6,$7}' | \
        awk '{gsub(/Chr0/, "", $2)} 1' | \
        awk '{gsub(/Chr/, "", $2)} 1' | \
        sed '1i\ID CHROM BIN_START BIN_END BIN_Middle N_VARIANTS PI' | \
        sed 's/ /\t/g' \
        > $i.windowed.pi.txt &
done

# Find the common ids
awk '{print $1,$2,$3,$7,"A"}' PUB.windowed.pi.txt > PUB.windowed.pi.txt-2
awk '{print $1,$2,$3,$7,"B"}' CHA.windowed.pi.txt > CHA.windowed.pi.txt-2
awk '{print $1,$2,$3,$7,"C"}' BAB.windowed.pi.txt > BAB.windowed.pi.txt-2
awk '{print $1,$2,$3,$7,"D"}' BAP.windowed.pi.txt > BAP.windowed.pi.txt-2
awk '{print $1,$2,$3,$7,"E"}' GLA.windowed.pi.txt > GLA.windowed.pi.txt-2
awk '{print $1,$2,$3,$7,"F"}' ANN.windowed.pi.txt > ANN.windowed.pi.txt-2
awk '{print $1,$2,$3,$7,"G"}' FRU.windowed.pi.txt > FRU.windowed.pi.txt-2
awk '{print $1,$2,$3,$7,"H"}' CHN.windowed.pi.txt > CHN.windowed.pi.txt-2

cat *windowed.pi.txt-2 | grep -v ID | sed '1i\SNP Chromosome Position PI Taxonomy' | sed 's/ /\t/g' > Group.pi.vg.txt

# Calculate the mean value of Pi
for i in *windowed.pi.txt ; do
    echo $i ;
    grep -v ID $i | awk '{ total += $7 } END { print total/NR }'
done

# Calculate the mean value of pi for all 
grep -v SNP Group.pi.vg.txt | awk '{ total += $4 } END { print total/NR }' 

```

### Pi with only window of 100000

```bash
# Working directory
cd /data/zhaojiantao/pepper/vg/11Pi/window-pi

# Calculate windowed nucleotide diversity
for i in PUB CHA BAB BAP FRU CHN GLA ANN; do
    for j in $(seq -w 1 12); do
        vcftools --gzvcf ../../04SNP/Miss_0.3/Chr$j.SNP.vcf.gz --keep ../../Pop.list/$i.list --maf 0.05 --window-pi 100000 --out $i.Chr$j &
    done &
done

# Find the common ids
cat PUB.Chr*.windowed.pi | grep -v CHROM | awk '{print $5,"A"}' > PUB.windowed.pi-1
cat CHA.Chr*.windowed.pi | grep -v CHROM | awk '{print $5,"B"}' > CHA.windowed.pi-1
cat BAB.Chr*.windowed.pi | grep -v CHROM | awk '{print $5,"C"}' > BAB.windowed.pi-1
cat BAP.Chr*.windowed.pi | grep -v CHROM | awk '{print $5,"D"}' > BAP.windowed.pi-1
cat GLA.Chr*.windowed.pi | grep -v CHROM | awk '{print $5,"E"}' > GLA.windowed.pi-1
cat ANN.Chr*.windowed.pi | grep -v CHROM | awk '{print $5,"F"}' > ANN.windowed.pi-1
cat FRU.Chr*.windowed.pi | grep -v CHROM | awk '{print $5,"G"}' > FRU.windowed.pi-1
cat CHN.Chr*.windowed.pi | grep -v CHROM | awk '{print $5,"H"}' > CHN.windowed.pi-1

# Combine the results
cat *pi-1 | sed '1i\PI Group' | sed 's/ /\t/g' > Combined.windowed.pi.txt

# Pair-wise significant t-test in R
R
group.pi <- read.table("Combined.windowed.pi.txt", header = TRUE)
pairwise_ttests <- pairwise.t.test(group.pi$PI, group.pi$Group, method = "t.test", p.adjust.method = "bonferroni")
pairwise_ttests

```

#### Permutation test using 4DTV SNPs

```bash
# Working directory
cd /data/zhaojiantao/pepper/vg/11Pi/permutations

# Run the permutation test
for group in PUB CHA BAB BAP FRU CHN GLA ANN; do
    # Generate the random 10 accession list
    for i in $(seq -w 1 100); do
    sort --random-sort ../../Pop.list/$group.list | head -10 > $group.Random$i.list
    done
    # Calculate per-site nucleotide diversity
    for i in $(seq -w 1 100); do
        for j in $(seq -w 1 12); do
            vcftools --gzvcf ../../064DTV/Chr$j.ffds.SNP.vcf.gz \
            --keep $group.Random$i.list --maf 0.05 --site-pi --out Random$i.Chr$j &
        done
    done
    # Calculate the average site pi
    for i in $(seq -w 1 100); do
    cat Random$i.Chr*.sites.pi | grep -v PI | grep -v nan | awk '{sum += $3} END {print sum/NR}' > Random$i.average.pi
    done
    # Combine the summary
    cat Random*.average.pi > $group.averaged.pi
    # Remove temporary files
    rm *log *sites.pi *list *average.pi
done

```

#### T-test in R
```r
data <- read.table("PUB.averaged.pi", header=FALSE)
mu <- mean(data$V1)
t.test(x=data, mu=mu, conf.level=0.95) 

data <- read.table("CHA.averaged.pi", header=FALSE)
mu <- mean(data$V1)
t.test(x=data, mu=mu, conf.level=0.95) 

data <- read.table("BAB.averaged.pi", header=FALSE)
mu <- mean(data$V1)
t.test(x=data, mu=mu, conf.level=0.95) 

data <- read.table("BAP.averaged.pi", header=FALSE)
mu <- mean(data$V1)
t.test(x=data, mu=mu, conf.level=0.95) 

data <- read.table("FRU.averaged.pi", header=FALSE)
mu <- mean(data$V1)
t.test(x=data, mu=mu, conf.level=0.95) 

data <- read.table("CHN.averaged.pi", header=FALSE)
mu <- mean(data$V1)
t.test(x=data, mu=mu, conf.level=0.95) 

data <- read.table("GLA.averaged.pi", header=FALSE)
mu <- mean(data$V1)
t.test(x=data, mu=mu, conf.level=0.95) 

data <- read.table("ANN.averaged.pi", header=FALSE)
mu <- mean(data$V1)
t.test(x=data, mu=mu, conf.level=0.95) 

```



### Effective population size

Software link: [SMC++](https://github.com/popgenmethods/smcpp#tips-for-using-smc)

#### ANN

```bash
# Working directory
cd /data/zhaojiantao/pepper/vg/12SMC++/ANN

# Convert vcf to smc format
for i in $(seq -w 1 12); do
    smc++ vcf2smc --cores 10 --missing-cutoff 10000 \
        /data/zhaojiantao/pepper/vg/04SNP/Miss_0.3/Chr$i.SNP.vcf.gz Chr$i.smc.gz Chr$i \
        ANN:PI_555642,PI_410407,chaotiaojiao,zhudachang,changanyingjiao,Grif_9122,PI_653659,PI_439378,PI_593500,Grif_9381,PI_645493,Grif_12454,PI_640503,PI_438538,PI_645508,PI_260611,PI_631152,PI_631129,qixingjiao,CXJ38,sanweijiao &
done

# Run smc++
smc++ estimate -o smc.out 6.96e-9 --cores 30 *smc.gz

# Generate the plots
smc++ plot ANN.smc20.jpeg -g 1 -c smc.out/model.final.json

```

#### GLA

```bash
# Working directory
cd /data/zhaojiantao/pepper/vg/12SMC++/GLA

# Convert vcf to smc 
for i in $(seq -w 1 12); do
    for j in $(cat GLA.list) ; do
        smc++ vcf2smc --cores 10 --missing-cutoff 10000 -d $j $j \
            /data/zhaojiantao/pepper/vg/04SNP/Miss_0.3/Chr$i.SNP.vcf.gz Chr$i.$j.smc.gz Chr$i \
            GLA:Grif_9142,PI_661082,PI_631136,PI_632930,PI_639662,PI_555616,PI_593527,PI_674459,PI_406948,PI_438567
    done &
done

# Run smc++
smc++ estimate -o smc.out 6.96e-9 --cores 60 *smc.gz

# Generate the plots
smc++ plot GLA.smc.jpeg -g 1 -c smc.out/model.final.json

```

#### BAB

```bash
# Working directory
cd /data/zhaojiantao/pepper/vg/12SMC++/BAB

# Convert vcf to smc 
for i in $(seq -w 1 12); do
    for j in $(cat BAB.list) ; do
        smc++ vcf2smc --cores 10 --missing-cutoff 10000 -d $j $j \
            /data/zhaojiantao/pepper/vg/04SNP/Miss_0.3/Chr$i.SNP.vcf.gz Chr$i.$j.smc.gz Chr$i \
            BAB:PI_260595,PI_446909,PI_439528,PI_441656,PI_639129,PI_633751,PI_337524,PI_441699,PI_441685
    done &
done

# Run smc++
smc++ estimate -o smc.out 6.96e-9 --cores 60 *smc.gz

# Generate the plots
smc++ plot BAB.smc.jpeg -g 1 -c smc.out/model.final.json

```

#### BAP

```bash
# Working directory
cd /data/zhaojiantao/pepper/vg/12SMC++/BAP

# Convert vcf to smc 
for i in $(seq -w 1 12); do
    for j in $(cat BAP.list) ; do
        smc++ vcf2smc --cores 10 --missing-cutoff 10000 -d $j $j \
            /data/zhaojiantao/pepper/vg/04SNP/Miss_0.3/Chr$i.SNP.vcf.gz Chr$i.$j.smc.gz Chr$i \
            BAP:PI_260538,PI_260567,PI_281310,PI_241679,Grif_9201,PI_585250,PI_632926,PI_257161,PI_439398,PI_159242
    done &
done

# Run smc++
smc++ estimate -o smc.out 6.96e-9 --cores 10 *smc.gz

# Generate the plots
smc++ plot BAP.smc.jpeg -g 1 -c smc.out/model.final.json

```

#### CHA

```bash
# Working directory
cd /data/zhaojiantao/pepper/vg/12SMC++/CHA

# Convert vcf to smc 
for i in $(seq -w 1 12); do
    for j in $(cat CHA.list) ; do
        smc++ vcf2smc --cores 10 --missing-cutoff 10000 -d $j $j \
            /data/zhaojiantao/pepper/vg/04SNP/Miss_0.3/Chr$i.SNP.vcf.gz Chr$i.$j.smc.gz Chr$i \
            CHA:PI_260429,PI_439414,PI_260427,PI_260430,PI_260431,PI_439415,PI_639659,PI_659106,PI_639652,PI_555612
    done &
done

# Run smc++
smc++ estimate -o smc.out 6.96e-9 --cores 20 *smc.gz

# Generate the plots
smc++ plot CHA.smc.jpeg -g 1 -c smc.out/model.final.json

```

#### CHN

```bash
# Working directory
cd /data/zhaojiantao/pepper/vg/12SMC++/CHN

# Convert vcf to smc 
for i in $(seq -w 1 12); do
    for j in $(cat CHN.list) ; do
        smc++ vcf2smc --cores 10 --missing-cutoff 10000 -d $j $j \
            /data/zhaojiantao/pepper/vg/04SNP/Miss_0.3/Chr$i.SNP.vcf.gz Chr$i.$j.smc.gz Chr$i \
            CHN:PI_260486,PI_441605,PI_257124,PI_439427,PI_355815,PI_666556,PI_438635,PI_238051,PI_281424,PI_159261
    done &
done

# Run smc++
smc++ estimate -o smc.out 6.96e-9 --cores 20 *smc.gz

# Generate the plots
smc++ plot CHN.smc.jpeg -g 1 -c smc.out/model.final.json

```

#### FRU

```bash
# Working directory
cd /data/zhaojiantao/pepper/vg/12SMC++/FRU

# Convert vcf to smc 
for i in $(seq -w 1 12); do
    for j in $(cat FRU.list) ; do
        smc++ vcf2smc --cores 10 --missing-cutoff 10000 -d $j $j \
            /data/zhaojiantao/pepper/vg/04SNP/Miss_0.3/Chr$i.SNP.vcf.gz Chr$i.$j.smc.gz Chr$i \
            FRU:PI_645681,PI_441690,xiaomila,PI_257104,Grif_9226,PI_360732,PI_631144,PI_368080,PI_439512,PI_281421
    done &
done

# Run smc++
smc++ estimate -o smc.out 6.96e-9 --cores 20 *smc.gz

# Generate the plots
smc++ plot FRU.smc.jpeg -g 1 -c smc.out/model.final.json

```

#### PUB

```bash
# Working directory
cd /data/zhaojiantao/pepper/vg/12SMC++/PUB

# Convert vcf to smc 
for i in $(seq -w 1 12); do
    for j in $(cat PUB.list) ; do
        smc++ vcf2smc --cores 10 --missing-cutoff 10000 -d $j $j \
            /data/zhaojiantao/pepper/vg/04SNP/Miss_0.3/Chr$i.SNP.vcf.gz Chr$i.$j.smc.gz Chr$i \
            PUB:Grif_1614,PI_585273,PI_593633,PI_585266,PI_585268,PI_585275,PI_593618,PI_593626,PI_593639,PI_593632
    done &
done

# Run smc++
smc++ estimate -o smc.out 6.96e-9 --cores 20 *smc.gz

# Generate the plots
smc++ plot PUB.smc.jpeg -g 1 -c smc.out/model.final.json 

```


#### Final combined plot

```bash
# Working directory
cd /data/zhaojiantao/pepper/vg/12SMC++

# Combined final plot
smc++ plot Final.smc.jpeg -g 1 */smc.out/model.final.json -x 100 14000 -y 600 20000

```

### GWAS

#### GLA.ANN

```bash
# Working directory
cd /data/zhaojiantao/pepper/vg/13GWAS/GLA.ANN

# Convert vcf to tped
for i in $(seq -w 1 12); do
    plink --vcf ../../04SNP/Miss_0.3/Chr$i.SNP.vcf.gz \
        --keep GLA.ANN.list \
        --allow-extra-chr \
        --recode \
        --double-id \
        --geno 0.7 \
        --maf 0.05 \
        --out GLA.ANN.Chr$i
    # Rename the chr id
    sed 's/CannZSG_//g' GLA.ANN.Chr$i.map > GLA.ANN.Chr$i.map2
    mv GLA.ANN.Chr$i.map2 GLA.ANN.Chr$i.map
    # Make bed
    plink --file GLA.ANN.Chr$i --make-bed --out GLA.ANN.Chr$i
done

# Calculate the effective number of SNPs
for i in $(seq -w 1 12); do
    java -jar /data/zhaojiantao/tools/GeneticTypeIError/gec.jar \
        --effect-number \
        --plink-binary GLA.ANN.Chr$i
    mv gec.sum GLA.ANN.Chr$i.gec.sum
done

# Calculate the effective number of SNPs of all chromosomes
cat *gec.sum | grep -v Effective_Number | awk '{ sum += $2 } END { print sum }' 

# Population structure
plink --vcf <(zcat ../../064DTV/504.4DTV.Miss0.1.vcf.gz | grep -v Chr00) \
    --keep GLA.ANN.list \
    --allow-extra-chr \
    --recode \
    --double-id \
    --geno 0.7 \
    --maf 0.05 \
    --out GLA.ANN.ffds

# Calculate PCA
plink --file GLA.ANN.ffds --pca --out GLA.ANN.ffds

# Remove unnecessary files
rm *log *nosex

# Prepare the PCA file as cofactor
awk '{print $1,$2,1,$3,$4,$5,$6,$7}' GLA.ANN.ffds.eigenvec > GLA.ANN.ffds.PCA.txt

# Transpose the ped file
for i in $(seq -w 1 12); do
    plink --file GLA.ANN.Chr$i --recode12 --output-missing-genotype 0 --transpose --out GLA.ANN.Chr$i
done

# Calculate kinship
plink --file GLA.ANN.ffds --recode12 --output-missing-genotype 0 --transpose --out GLA.ANN.ffds --allow-extra-chr
# Rename the chr id
cut -f1 -d ' ' GLA.ANN.ffds.tped | sed 's/CannZSG_Chr0//g' | sed 's/CannZSG_Chr//g' | \
    paste - GLA.ANN.ffds.tped | sed 's/\t/ /g' | cut -f1,3- -d ' ' > GLA.ANN.ffds.tped2
mv GLA.ANN.ffds.tped2 GLA.ANN.ffds.tped

# Estimate kinship
emmax-kin -v -d 10 GLA.ANN.ffds

# EMMAX
for i in $(seq -w 1 12); do
for j in FruitLength FruitWidth FruitRatio ; do
    emmax -v -d 10 -t GLA.ANN.Chr$i -p $j.txt -k GLA.ANN.ffds.BN.kinf -c GLA.ANN.ffds.PCA.txt -o GLA.ANN.$j.Chr$i &
done
done

# Prepare the output for manhattan plot
for i in FruitLength FruitWidth FruitRatio ; do
    cat GLA.ANN.$i.Chr*.ps | cut -f1 | sed 's/_/\t/g' | sed 's/Chr0//g' | sed 's/Chr//g' | \
        paste - <(cat GLA.ANN.$i.Chr*.ps) | \
        awk '{print $1,$2,$3,$5}' | sed '1i\CHR BP SNP P' | sed 's/ /\t/g' \
        > GLA.ANN.$i.GWAS.txt
done

# Calculate the LD for the peak SNP
plink --file GLA.ANN.Chr02 \
    --r2 \
    --ld-snp 2_160951499 \
    --ld-window-kb 1000 \
    --ld-window 1050600 \
    --ld-window-r2 0.351 \
    --out 2_160951499

for i in BAP.long BAP.bell ; do
    vcftools --gzvcf ../vcf.clean/Chr10.clean.vcf.gz --keep $i.list --site-pi --out $i.Chr10
    perl PI_slide_win.pl -w 100000 -s 10000 $i.Chr10.sites.pi
done

# Step 15.2 Calculate windowed pi 
for i in Ca.annuum Ca.bac.baccatum Ca.bac.pendulum Ca.chacoense Ca.chinense Ca.frutescens Ca.glabriusculum Ca.glabriusculum_II Ca.pubescens ; do
    for j in $(seq -w 1 12); do
    perl PI_slide_win.pl -w 100000 -s 10000 $i.Chr$j.sites.pi
done
done

```

##### Shape variations

```bash
# Working directory
cd /data/zhaojiantao/pepper/vg/13GWAS/GLA.ANN/Shape

# Extrac the samples
for i in $(seq -w 1 12); do
plink --file ../GLA.ANN.Chr$i --keep sample.list --make-bed --recode --out GLA.ANN.Chr$i
done

# Calculate the effective number of SNPs
for i in $(seq -w 1 12); do
    java -jar /data/zhaojiantao/tools/GeneticTypeIError/gec.jar \
        --effect-number \
        --plink-binary GLA.ANN.Chr$i
    mv gec.sum GLA.ANN.Chr$i.gec.sum
done

# Merge the chromosomal-level map/ped files
ls *map | sed 's/.map//g' > mergefiles.txt
plink --merge-list mergefiles.txt --make-bed --recode --out GLA.ANN

# LD prune
plink --file GLA.ANN \
    --indep-pairwise 50 10 0.2 \
    --allow-extra-chr --double-id \
    --out GLA.ANN
# Extract the prunned SNPs
plink --file GLA.ANN \
    --extract GLA.ANN.prune.in \
    --recode --double-id \
    --out GLA.ANN.prune

# Calculate PCA
plink --file GLA.ANN.prune --pca --out GLA.ANN.prune

# Prepare the PCA file as cofactor
awk '{print $1,$2,1,$3,$4,$5,$6,$7}' GLA.ANN.prune.eigenvec > GLA.ANN.prune.PCA.txt

# Transpose the ped file
for i in $(seq -w 1 12); do
    plink --file GLA.ANN.Chr$i --recode12 --output-missing-genotype 0 --transpose --out GLA.ANN.Chr$i &
done

# Calculate kinship
plink --file GLA.ANN.prune --recode12 --output-missing-genotype 0 --transpose --out GLA.ANN.prune --allow-extra-chr

# Estimate kinship
emmax-kin -v -d 10 GLA.ANN.prune

# EMMAX
for i in $(seq -w 1 12); do
    emmax -v -d 10 -t GLA.ANN.Chr$i -p Shape.normal.txt -k GLA.ANN.prune.BN.kinf -c GLA.ANN.prune.PCA.txt -o GLA.ANN.Chr$i &
done

# Prepare Manhattan plot
cat GLA.ANN.Chr*.ps | cut -f1 | sed 's/_/\t/g' | sed 's/Chr0//g' | sed 's/Chr//g' | \
    paste - <(cat GLA.ANN.Chr*.ps) | \
    awk '{print $1,$2,$3,$5}' | sed '1i\CHR BP SNP P' | sed 's/ /\t/g' \
    > GLA.ANN.Shape.GWAS.txt

```


#### BAB.BAP

```bash
# Working directory
cd /data/zhaojiantao/pepper/vg/13GWAS/BAB.BAP

# Convert vcf to tped
for i in $(seq -w 1 12); do
    plink --vcf ../../04SNP/Miss_0.3/Chr$i.SNP.vcf.gz \
        --keep BAB.BAP.list \
        --allow-extra-chr \
        --recode \
        --double-id \
        --geno 0.7 \
        --maf 0.05 \
        --out BAB.BAP.Chr$i
    # Rename the chr id
    sed 's/CannZSG_//g' BAB.BAP.Chr$i.map > BAB.BAP.Chr$i.map2
    mv BAB.BAP.Chr$i.map2 BAB.BAP.Chr$i.map
    # Make bed
    plink --file BAB.BAP.Chr$i --make-bed --out BAB.BAP.Chr$i
done

# Calculate the effective number of SNPs
for i in $(seq -w 1 12); do
    java -jar /data/zhaojiantao/tools/GeneticTypeIError/gec.jar \
        --effect-number \
        --plink-binary BAB.BAP.Chr$i
    mv gec.sum BAB.BAP.Chr$i.gec.sum
done

# Calculate the effective number of SNPs of all chromosomes
cat *gec.sum | grep -v Effective_Number | awk '{ sum += $2 } END { print sum }' 

# Population structure
plink --vcf <(zcat ../../064DTV/504.4DTV.Miss0.1.vcf.gz | grep -v Chr00) \
    --keep BAB.BAP.list \
    --allow-extra-chr \
    --recode \
    --double-id \
    --geno 0.7 \
    --maf 0.05 \
    --out BAB.BAP.ffds

# Calculate PCA
plink --file BAB.BAP.ffds --pca --out BAB.BAP.ffds --allow-extra-chr

# Remove unnecessary files
rm *log *nosex

# Prepare the PCA file as cofactor
awk '{print $1,$2,1,$3,$4,$5,$6,$7}' BAB.BAP.ffds.eigenvec > BAB.BAP.ffds.PCA.txt

# Transpose the ped file
for i in $(seq -w 1 12); do
    plink --file BAB.BAP.Chr$i --recode12 --output-missing-genotype 0 --transpose --out BAB.BAP.Chr$i
done

# Calculate kinship
plink --file BAB.BAP.ffds --recode12 --output-missing-genotype 0 --transpose --allow-extra-chr --out BAB.BAP.ffds
# Rename the chr id
cut -f1 -d ' ' BAB.BAP.ffds.tped | sed 's/CannZSG_Chr0//g' | sed 's/CannZSG_Chr//g' | \
    paste - BAB.BAP.ffds.tped | sed 's/\t/ /g' | cut -f1,3- -d ' ' > BAB.BAP.ffds.tped2
mv BAB.BAP.ffds.tped2 BAB.BAP.ffds.tped

# Estimate kinship
emmax-kin -v -d 10 BAB.BAP.ffds

# EMMAX
for i in $(seq -w 1 12); do
for j in FruitLength FruitWidth FruitRatio ; do
    emmax -v -d 10 -t BAB.BAP.Chr$i -p $j.txt -k BAB.BAP.ffds.BN.kinf -c BAB.BAP.ffds.PCA.txt -o BAB.BAP.$j.Chr$i &
done
done

# Prepare the output for manhattan plot
for i in FruitLength FruitWidth FruitRatio ; do
    cat BAB.BAP.$i.Chr*.ps | cut -f1 | sed 's/_/\t/g' | sed 's/Chr0//g' | sed 's/Chr//g' | \
        paste - <(cat BAB.BAP.$i.Chr*.ps) | \
        awk '{print $1,$2,$3,$5}' | sed '1i\CHR BP SNP P' | sed 's/ /\t/g' \
        > BAB.BAP.$i.GWAS.txt
done

# Calculate the LD for the peak SNP
plink --file BAB.BAP.Chr02 \
    --r2 \
    --ld-snp Chr02_172705522 \
    --ld-window-kb 1000 \
    --ld-window 1050600 \
    --ld-window-r2 0.351 \
    --out Chr02_172705522

plink --file BAB.BAP.Chr12 \
    --r2 \
    --ld-snp Chr12_236447804 \
    --ld-window-kb 1000 \
    --ld-window 1050600 \
    --ld-window-r2 0.351 \
    --out Chr12_236447804

bedtools intersect -a ../../../References/S8/Zhangshugang.gff -b AA.bed | \
    grep gene | sed 's/;/\t/g' | cut -f9 | sed 's/ID=//g' > Candidate.genes

for i in BAP.long BAP.bell ; do
    vcftools --gzvcf ../../04SNP/Miss_0.3/Chr10.SNP.vcf.gz --keep $i.list --site-pi --out $i.Chr10 &
    perl PI_slide_win.pl -w 100000 -s 10000 $i.Chr10.sites.pi
done

# Prepare the pi plot for the candidate region
grep -v "#" BAP.bell.Chr10.sites.pi.w100000.s10000.txt | \
    awk '$3>120500000 && $4<122600000 {print $4,$7,"bell"}' > BAP.bell.Chr10.candidate.pi
grep -v "#" BAP.long.Chr10.sites.pi.w100000.s10000.txt | \
    awk '$3>120500000 && $4<122600000 {print $4,$7,"long"}' > BAP.long.Chr10.candidate.pi
grep -v "CHROM" ../../11Pi/BAB.windowed.pi.txt | grep Chr10 | \
    awk '$3>120500000 && $4<122600000 {print $3,$7,"BAB"}' > BAB.Chr10.candidate.pi

```

##### Shape variations

```bash
# Working directory
cd /data/zhaojiantao/pepper/vg/13GWAS/BAB.BAP/Shape

# Extrac the samples
for i in $(seq -w 1 12); do
plink --file ../BAB.BAP.Chr$i --keep sample.list --make-bed --recode --out BAB.BAP.Chr$i &
done

# Calculate the effective number of SNPs
for i in $(seq -w 1 12); do
    java -jar /data/zhaojiantao/tools/GeneticTypeIError/gec.jar \
        --effect-number \
        --plink-binary BAB.BAP.Chr$i
    mv gec.sum BAB.BAP.Chr$i.gec.sum
done

# Calculate the effective number of SNPs of all chromosomes
cat *gec.sum | grep -v Effective_Number | awk '{ sum += $2 } END { print sum }' 

# Merge the chromosomal-level map/ped files
ls *map | sed 's/.map//g' > mergefiles.txt
plink --merge-list mergefiles.txt --make-bed --recode --out BAB.BAP

# Extrac the samples
for i in $(seq -w 1 12); do
plink --file ../BAB.BAP.Chr$i --keep sample.list --make-bed --recode --out BAB.BAP.Chr$i
done

# Calculate the effective number of SNPs
for i in $(seq -w 1 12); do
    java -jar /data/zhaojiantao/tools/GeneticTypeIError/gec.jar \
        --effect-number \
        --plink-binary BAB.BAP.Chr$i
    mv gec.sum BAB.BAP.Chr$i.gec.sum
done

# Merge the chromosomal-level map/ped files
ls *map | sed 's/.map//g' > mergefiles.txt
plink --merge-list mergefiles.txt --make-bed --recode --out BAB.BAP

# LD prune
plink --file BAB.BAP \
    --indep-pairwise 50 10 0.2 \
    --allow-extra-chr --double-id \
    --out BAB.BAP
# Extract the prunned SNPs
plink --file BAB.BAP \
    --extract BAB.BAP.prune.in \
    --recode --double-id \
    --out BAB.BAP.prune

# Calculate PCA
plink --file BAB.BAP.prune --pca --out BAB.BAP.prune

# Prepare the PCA file as cofactor
awk '{print $1,$2,1,$3,$4,$5,$6,$7}' BAB.BAP.prune.eigenvec > BAB.BAP.prune.PCA.txt

# Transpose the ped file
for i in $(seq -w 1 12); do
    plink --file BAB.BAP.Chr$i --recode12 --output-missing-genotype 0 --transpose --out BAB.BAP.Chr$i &
done

# Calculate kinship
plink --file BAB.BAP.prune --recode12 --output-missing-genotype 0 --transpose --out BAB.BAP.prune --allow-extra-chr

# Estimate kinship
emmax-kin -v -d 10 BAB.BAP.prune

# EMMAX
for i in $(seq -w 1 12); do
    emmax -v -d 10 -t BAB.BAP.Chr$i -p Shape.normal.txt -k BAB.BAP.prune.BN.kinf -c BAB.BAP.prune.PCA.txt -o BAB.BAP.Chr$i &
done

# Prepare the output for manhattan plot
cat BAB.BAP.Chr*.ps | cut -f1 | sed 's/_/\t/g' | sed 's/Chr0//g' | sed 's/Chr//g' | \
    paste - <(cat BAB.BAP.Chr*.ps) | \
    awk '{print $1,$2,$3,$5}' | sed '1i\CHR BP SNP P' | sed 's/ /\t/g' \
    > BAB.BAP.Shape.GWAS.txt

```

#### Seed size

##### Extract the phenotype

```bash
# Working directory
cd /data/zhaojiantao/pepper/vg/13GWAS/Seed

# Extract the samples
for i in BAB.BAP GLA.ANN CHN FRU Annuum ; do
grep -f ../../Pop.list/$i.list seed.txt | awk '{print $1,$1,$2}' | sort -k1,1 | sed 's/ /\t/g' > $i/L.txt
grep -f ../../Pop.list/$i.list seed.txt | awk '{print $1,$1,$3}' | sort -k1,1 | sed 's/ /\t/g' > $i/R.txt
grep -f ../../Pop.list/$i.list seed.txt | awk '{print $1,$1,$4}' | sort -k1,1 | sed 's/ /\t/g' > $i/S.txt
grep -f ../../Pop.list/$i.list seed.txt | awk '{print $1,$1,$5}' | sort -k1,1 | sed 's/ /\t/g' > $i/W.txt
done

```

##### GLA.ANN

```bash
# Working directory
cd /data/zhaojiantao/pepper/vg/13GWAS/Seed/GLA.ANN

# Prepare the sample list
awk '{print $1,$1}' L.txt | sed 's/ /\t/g' > GLA.ANN.list

# Convert vcf to tped
for i in $(seq -w 1 12); do
    plink --vcf ../../../04SNP/Miss_0.3/Chr$i.SNP.vcf.gz \
        --keep GLA.ANN.list \
        --allow-extra-chr \
        --recode \
        --double-id \
        --geno 0.7 \
        --maf 0.05 \
        --out GLA.ANN.Chr$i
    # Rename the chr id
    sed 's/CannZSG_//g' GLA.ANN.Chr$i.map > GLA.ANN.Chr$i.map2
    mv GLA.ANN.Chr$i.map2 GLA.ANN.Chr$i.map
    # Make bed
    plink --file GLA.ANN.Chr$i --make-bed --out GLA.ANN.Chr$i
done

# Calculate the effective number of SNPs
for i in $(seq -w 1 12); do
    java -jar /data/zhaojiantao/tools/GeneticTypeIError/gec.jar \
        --effect-number \
        --plink-binary GLA.ANN.Chr$i
    mv gec.sum GLA.ANN.Chr$i.gec.sum
done

# Calculate the effective number of SNPs of all chromosomes
cat *gec.sum | grep -v Effective_Number | awk '{ sum += $2 } END { print sum }' 

# Population structure
plink --vcf <(zcat ../../../064DTV/504.4DTV.Miss0.1.vcf.gz | grep -v Chr00) \
    --keep GLA.ANN.list \
    --allow-extra-chr \
    --recode \
    --double-id \
    --geno 0.7 \
    --maf 0.05 \
    --out GLA.ANN.ffds

# Calculate PCA
plink --file GLA.ANN.ffds --pca --out GLA.ANN.ffds --allow-extra-chr

# Remove unnecessary files
rm *log *nosex

# Prepare the PCA file as cofactor
awk '{print $1,$2,1,$3,$4,$5,$6,$7}' GLA.ANN.ffds.eigenvec > GLA.ANN.ffds.PCA.txt

# Transpose the ped file
for i in $(seq -w 1 12); do
    plink --file GLA.ANN.Chr$i --recode12 --output-missing-genotype 0 --transpose --out GLA.ANN.Chr$i
done

# Calculate kinship
plink --file GLA.ANN.ffds --recode12 --output-missing-genotype 0 --transpose --allow-extra-chr --out GLA.ANN.ffds
# Rename the chr id
cut -f1 -d ' ' GLA.ANN.ffds.tped | sed 's/CannZSG_Chr0//g' | sed 's/CannZSG_Chr//g' | \
    paste - GLA.ANN.ffds.tped | sed 's/\t/ /g' | cut -f1,3- -d ' ' > GLA.ANN.ffds.tped2
mv GLA.ANN.ffds.tped2 GLA.ANN.ffds.tped

# Estimate kinship
emmax-kin -v -d 10 GLA.ANN.ffds

# EMMAX
for i in $(seq -w 1 12); do
for j in L R S W ; do
    emmax -v -d 10 -t GLA.ANN.Chr$i -p $j.txt -k GLA.ANN.ffds.BN.kinf -c GLA.ANN.ffds.PCA.txt -o GLA.ANN.$j.Chr$i &
done
done

# Prepare the output for manhattan plot
for i in L R S W ; do
    cat GLA.ANN.$i.Chr*.ps | cut -f1 | sed 's/_/\t/g' | sed 's/Chr0//g' | sed 's/Chr//g' | \
        paste - <(cat GLA.ANN.$i.Chr*.ps) | \
        awk '{print $1,$2,$3,$5}' | sed '1i\CHR BP SNP P' | sed 's/ /\t/g' \
        > GLA.ANN.$i.GWAS.txt
done

plink --file GLA.ANN.Chr07 \
    --r2 \
    --ld-snp Chr07_13218329 \
    --ld-window-kb 1000 \
    --ld-window 1000000 \
    --ld-window-r2 0.5 \
    --out Chr07_13218329

```

###### Manhattan plot

```r
library(qqman)
library(data.table)
GLA.ANN.L.GWAS <- fread("GLA.ANN.L.GWAS.txt", header = TRUE)
jpeg("GLA.ANN.L.GWAS.jpeg", height = 8, width = 16, units="cm", res = 300)
manhattan(GLA.ANN.L.GWAS, chr = "CHR", bp = "BP",
          snp = "SNP",col = c("orange", "darkgreen"),
          p = "P", ylim=c(0,10), logp = TRUE,cex=0.3, 
          cex.axis = 0.5,cex.main = 1, 
          suggestiveline = FALSE, genomewideline = -log10(1/2.58721e+06))
dev.off()

GLA.ANN.R.GWAS <- fread("GLA.ANN.R.GWAS.txt", header = TRUE)
jpeg("GLA.ANN.R.GWAS.jpeg", height = 8, width = 16, units="cm", res = 300)
manhattan(GLA.ANN.R.GWAS, chr = "CHR", bp = "BP",
          snp = "SNP",col = c("orange", "darkgreen"),
          p = "P", ylim=c(0,30), logp = TRUE,cex=0.3, 
          cex.axis = 0.5,cex.main = 1, 
          suggestiveline = FALSE, genomewideline = -log10(1/2.58721e+06))
dev.off()

GLA.ANN.S.GWAS <- fread("GLA.ANN.S.GWAS.txt", header = TRUE)
jpeg("GLA.ANN.S.GWAS.jpeg", height = 8, width = 16, units="cm", res = 300)
manhattan(GLA.ANN.S.GWAS, chr = "CHR", bp = "BP",
          snp = "SNP",col = c("orange", "darkgreen"),
          p = "P", ylim=c(0,10), logp = TRUE,cex=0.3, 
          cex.axis = 0.5,cex.main = 1, 
          suggestiveline = FALSE, genomewideline = -log10(1/2.58721e+06))
dev.off()

GLA.ANN.W.GWAS <- fread("GLA.ANN.W.GWAS.txt", header = TRUE)
jpeg("GLA.ANN.W.GWAS.jpeg", height = 8, width = 16, units="cm", res = 300)
manhattan(GLA.ANN.W.GWAS, chr = "CHR", bp = "BP",
          snp = "SNP",col = c("orange", "darkgreen"),
          p = "P", ylim=c(0,10), logp = TRUE,cex=0.3, 
          cex.axis = 0.5,cex.main = 1, 
          suggestiveline = FALSE, genomewideline = -log10(1/2.58721e+06))
dev.off()

```

##### FRU

```bash
# Working directory
cd /data/zhaojiantao/pepper/vg/13GWAS/Seed/FRU

# Prepare the sample list
awk '{print $1,$1}' L.txt | sed 's/ /\t/g' > FRU.list

# Convert vcf to tped
for i in $(seq -w 1 12); do
    plink --vcf ../../../04SNP/Miss_0.3/Chr$i.SNP.vcf.gz \
        --keep FRU.list \
        --allow-extra-chr \
        --recode \
        --double-id \
        --geno 0.7 \
        --maf 0.05 \
        --out FRU.Chr$i
    # Rename the chr id
    sed 's/CannZSG_//g' FRU.Chr$i.map > FRU.Chr$i.map2
    mv FRU.Chr$i.map2 FRU.Chr$i.map
    # Make bed
    plink --file FRU.Chr$i --make-bed --out FRU.Chr$i
done

# Calculate the effective number of SNPs
for i in $(seq -w 1 12); do
    java -jar /data/zhaojiantao/tools/GeneticTypeIError/gec.jar \
        --effect-number \
        --plink-binary FRU.Chr$i
    mv gec.sum FRU.Chr$i.gec.sum
done

# Calculate the effective number of SNPs of all chromosomes
cat *gec.sum | grep -v Effective_Number | awk '{ sum += $2 } END { print sum }' 

# Population structure
plink --vcf <(zcat ../../../064DTV/504.4DTV.Miss0.1.vcf.gz | grep -v Chr00) \
    --keep FRU.list \
    --allow-extra-chr \
    --recode \
    --double-id \
    --geno 0.7 \
    --maf 0.05 \
    --out FRU.ffds

# Calculate PCA
plink --file FRU.ffds --pca --out FRU.ffds --allow-extra-chr

# Remove unnecessary files
rm *log *nosex

# Prepare the PCA file as cofactor
awk '{print $1,$2,1,$3,$4,$5,$6,$7}' FRU.ffds.eigenvec > FRU.ffds.PCA.txt

# Transpose the ped file
for i in $(seq -w 1 12); do
    plink --file FRU.Chr$i --recode12 --output-missing-genotype 0 --transpose --out FRU.Chr$i
done

# Calculate kinship
plink --file FRU.ffds --recode12 --output-missing-genotype 0 --transpose --allow-extra-chr --out FRU.ffds
# Rename the chr id
cut -f1 -d ' ' FRU.ffds.tped | sed 's/CannZSG_Chr0//g' | sed 's/CannZSG_Chr//g' | \
    paste - FRU.ffds.tped | sed 's/\t/ /g' | cut -f1,3- -d ' ' > FRU.ffds.tped2
mv FRU.ffds.tped2 FRU.ffds.tped

# Estimate kinship
emmax-kin -v -d 10 FRU.ffds

# EMMAX
for i in $(seq -w 1 12); do
for j in L R S W ; do
    emmax -v -d 10 -t FRU.Chr$i -p $j.txt -k FRU.ffds.BN.kinf -c FRU.ffds.PCA.txt -o FRU.$j.Chr$i &
done
done

# Prepare the output for manhattan plot
for i in L R S W ; do
    cat FRU.$i.Chr*.ps | cut -f1 | sed 's/_/\t/g' | sed 's/Chr0//g' | sed 's/Chr//g' | \
        paste - <(cat FRU.$i.Chr*.ps) | \
        awk '{print $1,$2,$3,$5}' | sed '1i\CHR BP SNP P' | sed 's/ /\t/g' \
        > FRU.$i.GWAS.txt
done

```

###### Manhattan plot

```r
library(qqman)
library(data.table)
FRU.L.GWAS <- fread("FRU.L.GWAS.txt", header = TRUE)
jpeg("FRU.L.GWAS.jpeg", height = 8, width = 16, units="cm", res = 300)
manhattan(FRU.L.GWAS, chr = "CHR", bp = "BP",
          snp = "SNP",col = c("orange", "darkgreen"),
          p = "P", ylim=c(0,10), logp = TRUE,cex=0.3, 
          cex.axis = 0.5,cex.main = 1, 
          suggestiveline = FALSE, genomewideline = -log10(1/2.58721e+06))
dev.off()

FRU.R.GWAS <- fread("FRU.R.GWAS.txt", header = TRUE)
jpeg("FRU.R.GWAS.jpeg", height = 8, width = 16, units="cm", res = 300)
manhattan(FRU.R.GWAS, chr = "CHR", bp = "BP",
          snp = "SNP",col = c("orange", "darkgreen"),
          p = "P", ylim=c(0,10), logp = TRUE,cex=0.3, 
          cex.axis = 0.5,cex.main = 1, 
          suggestiveline = FALSE, genomewideline = -log10(1/2.58721e+06))
dev.off()

FRU.S.GWAS <- fread("FRU.S.GWAS.txt", header = TRUE)
jpeg("FRU.S.GWAS.jpeg", height = 8, width = 16, units="cm", res = 300)
manhattan(FRU.S.GWAS, chr = "CHR", bp = "BP",
          snp = "SNP",col = c("orange", "darkgreen"),
          p = "P", ylim=c(0,10), logp = TRUE,cex=0.3, 
          cex.axis = 0.5,cex.main = 1, 
          suggestiveline = FALSE, genomewideline = -log10(1/2.58721e+06))
dev.off()

FRU.W.GWAS <- fread("FRU.W.GWAS.txt", header = TRUE)
jpeg("FRU.W.GWAS.jpeg", height = 8, width = 16, units="cm", res = 300)
manhattan(FRU.W.GWAS, chr = "CHR", bp = "BP",
          snp = "SNP",col = c("orange", "darkgreen"),
          p = "P", ylim=c(0,10), logp = TRUE,cex=0.3, 
          cex.axis = 0.5,cex.main = 1, 
          suggestiveline = FALSE, genomewideline = -log10(1/2.58721e+06))
dev.off()

```


##### CHN

```bash
# Working directory
cd /data/zhaojiantao/pepper/vg/13GWAS/Seed/CHN

# Prepare the sample list
awk '{print $1,$1}' L.txt | sed 's/ /\t/g' > CHN.list

# Convert vcf to tped
for i in $(seq -w 1 12); do
    plink --vcf ../../../04SNP/Miss_0.3/Chr$i.SNP.vcf.gz \
        --keep CHN.list \
        --allow-extra-chr \
        --recode \
        --double-id \
        --geno 0.7 \
        --maf 0.05 \
        --out CHN.Chr$i
    # Rename the chr id
    sed 's/CannZSG_//g' CHN.Chr$i.map > CHN.Chr$i.map2
    mv CHN.Chr$i.map2 CHN.Chr$i.map
    # Make bed
    plink --file CHN.Chr$i --make-bed --out CHN.Chr$i
done

# Calculate the effective number of SNPs
for i in $(seq -w 1 12); do
    java -jar /data/zhaojiantao/tools/GeneticTypeIError/gec.jar \
        --effect-number \
        --plink-binary CHN.Chr$i
    mv gec.sum CHN.Chr$i.gec.sum
done

# Calculate the effective number of SNPs of all chromosomes
cat *gec.sum | grep -v Effective_Number | awk '{ sum += $2 } END { print sum }' 

# Population structure
plink --vcf <(zcat ../../../064DTV/504.4DTV.Miss0.1.vcf.gz | grep -v Chr00) \
    --keep CHN.list \
    --allow-extra-chr \
    --recode \
    --double-id \
    --geno 0.7 \
    --maf 0.05 \
    --out CHN.ffds

# Calculate PCA
plink --file CHN.ffds --pca --out CHN.ffds --allow-extra-chr

# Remove unnecessary files
rm *log *nosex

# Prepare the PCA file as cofactor
awk '{print $1,$2,1,$3,$4,$5,$6,$7}' CHN.ffds.eigenvec > CHN.ffds.PCA.txt

# Transpose the ped file
for i in $(seq -w 1 12); do
    plink --file CHN.Chr$i --recode12 --output-missing-genotype 0 --transpose --out CHN.Chr$i
done

# Calculate kinship
plink --file CHN.ffds --recode12 --output-missing-genotype 0 --transpose --allow-extra-chr --out CHN.ffds
# Rename the chr id
cut -f1 -d ' ' CHN.ffds.tped | sed 's/CannZSG_Chr0//g' | sed 's/CannZSG_Chr//g' | \
    paste - CHN.ffds.tped | sed 's/\t/ /g' | cut -f1,3- -d ' ' > CHN.ffds.tped2
mv CHN.ffds.tped2 CHN.ffds.tped

# Estimate kinship
emmax-kin -v -d 10 CHN.ffds

# EMMAX
for i in $(seq -w 1 12); do
for j in L R S W ; do
    emmax -v -d 10 -t CHN.Chr$i -p $j.txt -k CHN.ffds.BN.kinf -c CHN.ffds.PCA.txt -o CHN.$j.Chr$i &
done
done

# Prepare the output for manhattan plot
for i in L R S W ; do
    cat CHN.$i.Chr*.ps | cut -f1 | sed 's/_/\t/g' | sed 's/Chr0//g' | sed 's/Chr//g' | \
        paste - <(cat CHN.$i.Chr*.ps) | \
        awk '{print $1,$2,$3,$5}' | sed '1i\CHR BP SNP P' | sed 's/ /\t/g' \
        > CHN.$i.GWAS.txt
done

# Candidate region
plink --file CHN.Chr08 \
    --r2 \
    --ld-snp Chr08_107756202 \
    --ld-window-kb 1000 \
    --ld-window 1000000 \
    --ld-window-r2 0.5 \
    --out Chr08_107756202

```

###### Manhattan plot

```r
library(qqman)
library(data.table)
CHN.L.GWAS <- fread("CHN.L.GWAS.txt", header = TRUE)
jpeg("CHN.L.GWAS.jpeg", height = 8, width = 16, units="cm", res = 300)
manhattan(CHN.L.GWAS, chr = "CHR", bp = "BP",
          snp = "SNP",col = c("orange", "darkgreen"),
          p = "P", ylim=c(0,10), logp = TRUE,cex=0.3, 
          cex.axis = 0.5,cex.main = 1, 
          suggestiveline = FALSE, genomewideline = -log10(1/2.48239e+06))
dev.off()

CHN.R.GWAS <- fread("CHN.R.GWAS.txt", header = TRUE)
jpeg("CHN.R.GWAS.jpeg", height = 8, width = 16, units="cm", res = 300)
manhattan(CHN.R.GWAS, chr = "CHR", bp = "BP",
          snp = "SNP",col = c("orange", "darkgreen"),
          p = "P", ylim=c(0,10), logp = TRUE,cex=0.3, 
          cex.axis = 0.5,cex.main = 1, 
          suggestiveline = FALSE, genomewideline = -log10(1/2.48239e+06))
dev.off()

CHN.S.GWAS <- fread("CHN.S.GWAS.txt", header = TRUE)
jpeg("CHN.S.GWAS.jpeg", height = 8, width = 16, units="cm", res = 300)
manhattan(CHN.S.GWAS, chr = "CHR", bp = "BP",
          snp = "SNP",col = c("orange", "darkgreen"),
          p = "P", ylim=c(0,10), logp = TRUE,cex=0.3, 
          cex.axis = 0.5,cex.main = 1, 
          suggestiveline = FALSE, genomewideline = -log10(1/2.48239e+06))
dev.off()

CHN.W.GWAS <- fread("CHN.W.GWAS.txt", header = TRUE)
jpeg("CHN.W.GWAS.jpeg", height = 8, width = 16, units="cm", res = 300)
manhattan(CHN.W.GWAS, chr = "CHR", bp = "BP",
          snp = "SNP",col = c("orange", "darkgreen"),
          p = "P", ylim=c(0,10), logp = TRUE,cex=0.3, 
          cex.axis = 0.5,cex.main = 1, 
          suggestiveline = FALSE, genomewideline = -log10(1/2.48239e+06))
dev.off()

```

##### Annuum

```bash
# Working directory
cd /data/zhaojiantao/pepper/vg/13GWAS/Seed/Annuum

# Prepare the sample list
awk '{print $1,$1}' L.txt | sed 's/ /\t/g' > Annuum.list

# Convert vcf to tped
for i in $(seq -w 1 12); do
    plink --vcf ../../../04SNP/Miss_0.3/Chr$i.SNP.vcf.gz \
        --keep Annuum.list \
        --allow-extra-chr \
        --recode \
        --double-id \
        --geno 0.7 \
        --maf 0.05 \
        --out Annuum.Chr$i
    # Rename the chr id
    sed 's/CannZSG_//g' Annuum.Chr$i.map > Annuum.Chr$i.map2
    mv Annuum.Chr$i.map2 Annuum.Chr$i.map
    # Make bed
    plink --file Annuum.Chr$i --make-bed --out Annuum.Chr$i
done

# Calculate the effective number of SNPs
for i in $(seq -w 1 12); do
    java -jar /data/zhaojiantao/tools/GeneticTypeIError/gec.jar \
        --effect-number \
        --plink-binary Annuum.Chr$i
    mv gec.sum Annuum.Chr$i.gec.sum
done

# Calculate the effective number of SNPs of all chromosomes
cat *gec.sum | grep -v Effective_Number | awk '{ sum += $2 } END { print sum }' 

# Population structure
plink --vcf <(zcat ../../../064DTV/504.4DTV.Miss0.1.vcf.gz | grep -v Chr00) \
    --keep Annuum.list \
    --allow-extra-chr \
    --recode \
    --double-id \
    --geno 0.7 \
    --maf 0.05 \
    --out Annuum.ffds

# Calculate PCA
plink --file Annuum.ffds --pca --out Annuum.ffds --allow-extra-chr

# Remove unnecessary files
rm *log *nosex

# Prepare the PCA file as cofactor
awk '{print $1,$2,1,$3,$4,$5,$6,$7}' Annuum.ffds.eigenvec > Annuum.ffds.PCA.txt

# Transpose the ped file
for i in $(seq -w 1 12); do
    plink --file Annuum.Chr$i --recode12 --output-missing-genotype 0 --transpose --out Annuum.Chr$i
done

# Calculate kinship
plink --file Annuum.ffds --recode12 --output-missing-genotype 0 --transpose --allow-extra-chr --out Annuum.ffds
# Rename the chr id
cut -f1 -d ' ' Annuum.ffds.tped | sed 's/CannZSG_Chr0//g' | sed 's/CannZSG_Chr//g' | \
    paste - Annuum.ffds.tped | sed 's/\t/ /g' | cut -f1,3- -d ' ' > Annuum.ffds.tped2
mv Annuum.ffds.tped2 Annuum.ffds.tped

# Estimate kinship
emmax-kin -v -d 10 Annuum.ffds

# EMMAX
for i in $(seq -w 1 12); do
for j in L R S W ; do
    emmax -v -d 10 -t Annuum.Chr$i -p $j.txt -k Annuum.ffds.BN.kinf -c Annuum.ffds.PCA.txt -o Annuum.$j.Chr$i &
done
done

# Prepare the output for manhattan plot
for i in L R S W ; do
    cat Annuum.$i.Chr*.ps | cut -f1 | sed 's/_/\t/g' | sed 's/Chr0//g' | sed 's/Chr//g' | \
        paste - <(cat Annuum.$i.Chr*.ps) | \
        awk '{print $1,$2,$3,$5}' | sed '1i\CHR BP SNP P' | sed 's/ /\t/g' \
        > Annuum.$i.GWAS.txt
done

```

###### Manhattan plot

```r
library(qqman)
library(data.table)
Annuum.L.GWAS <- fread("Annuum.L.GWAS.txt", header = TRUE)
jpeg("Annuum.L.GWAS.jpeg", height = 8, width = 16, units="cm", res = 300)
manhattan(Annuum.L.GWAS, chr = "CHR", bp = "BP",
          snp = "SNP",col = c("orange", "darkgreen"),
          p = "P", ylim=c(0,10), logp = TRUE,cex=0.3, 
          cex.axis = 0.5,cex.main = 1, 
          suggestiveline = FALSE, genomewideline = -log10(1/3.19274e+06))
dev.off()

Annuum.R.GWAS <- fread("Annuum.R.GWAS.txt", header = TRUE)
jpeg("Annuum.R.GWAS.jpeg", height = 8, width = 16, units="cm", res = 300)
manhattan(Annuum.R.GWAS, chr = "CHR", bp = "BP",
          snp = "SNP",col = c("orange", "darkgreen"),
          p = "P", ylim=c(0,10), logp = TRUE,cex=0.3, 
          cex.axis = 0.5,cex.main = 1, 
          suggestiveline = FALSE, genomewideline = -log10(1/3.19274e+06))
dev.off()

Annuum.S.GWAS <- fread("Annuum.S.GWAS.txt", header = TRUE)
jpeg("Annuum.S.GWAS.jpeg", height = 8, width = 16, units="cm", res = 300)
manhattan(Annuum.S.GWAS, chr = "CHR", bp = "BP",
          snp = "SNP",col = c("orange", "darkgreen"),
          p = "P", ylim=c(0,10), logp = TRUE,cex=0.3, 
          cex.axis = 0.5,cex.main = 1, 
          suggestiveline = FALSE, genomewideline = -log10(1/3.19274e+06))
dev.off()

Annuum.W.GWAS <- fread("Annuum.W.GWAS.txt", header = TRUE)
jpeg("Annuum.W.GWAS.jpeg", height = 8, width = 16, units="cm", res = 300)
manhattan(Annuum.W.GWAS, chr = "CHR", bp = "BP",
          snp = "SNP",col = c("orange", "darkgreen"),
          p = "P", ylim=c(0,10), logp = TRUE,cex=0.3, 
          cex.axis = 0.5,cex.main = 1, 
          suggestiveline = FALSE, genomewideline = -log10(1/3.19274e+06))
dev.off()

```

##### BAB.BAP

```bash
# Working directory
cd /data/zhaojiantao/pepper/vg/13GWAS/Seed/BAB.BAP

# Prepare the sample list
awk '{print $1,$1}' L.txt | sed 's/ /\t/g' > BAB.BAP.list

# Convert vcf to tped
for i in $(seq -w 1 12); do
    plink --vcf ../../../04SNP/Miss_0.3/Chr$i.SNP.vcf.gz \
        --keep BAB.BAP.list \
        --allow-extra-chr \
        --recode \
        --double-id \
        --geno 0.7 \
        --maf 0.05 \
        --out BAB.BAP.Chr$i
    # Rename the chr id
    sed 's/CannZSG_//g' BAB.BAP.Chr$i.map > BAB.BAP.Chr$i.map2
    mv BAB.BAP.Chr$i.map2 BAB.BAP.Chr$i.map
    # Make bed
    plink --file BAB.BAP.Chr$i --make-bed --out BAB.BAP.Chr$i
done

# Calculate the effective number of SNPs
for i in $(seq -w 1 12); do
    java -jar /data/zhaojiantao/tools/GeneticTypeIError/gec.jar \
        --effect-number \
        --plink-binary BAB.BAP.Chr$i
    mv gec.sum BAB.BAP.Chr$i.gec.sum
done

# Calculate the effective number of SNPs of all chromosomes: 3.19274e+06
cat *gec.sum | grep -v Effective_Number | awk '{ sum += $2 } END { print sum }' 

# Population structure
plink --vcf <(zcat ../../../064DTV/504.4DTV.Miss0.1.vcf.gz | grep -v Chr00) \
    --keep BAB.BAP.list \
    --allow-extra-chr \
    --recode \
    --double-id \
    --geno 0.7 \
    --maf 0.05 \
    --out BAB.BAP.ffds

# Calculate PCA
plink --file BAB.BAP.ffds --pca --out BAB.BAP.ffds --allow-extra-chr

# Remove unnecessary files
rm *log *nosex

# Prepare the PCA file as cofactor
awk '{print $1,$2,1,$3,$4,$5,$6,$7}' BAB.BAP.ffds.eigenvec > BAB.BAP.ffds.PCA.txt

# Transpose the ped file
for i in $(seq -w 1 12); do
    plink --file BAB.BAP.Chr$i --recode12 --output-missing-genotype 0 --transpose --out BAB.BAP.Chr$i
done

# Calculate kinship
plink --file BAB.BAP.ffds --recode12 --output-missing-genotype 0 --transpose --allow-extra-chr --out BAB.BAP.ffds
# Rename the chr id
cut -f1 -d ' ' BAB.BAP.ffds.tped | sed 's/CannZSG_Chr0//g' | sed 's/CannZSG_Chr//g' | \
    paste - BAB.BAP.ffds.tped | sed 's/\t/ /g' | cut -f1,3- -d ' ' > BAB.BAP.ffds.tped2
mv BAB.BAP.ffds.tped2 BAB.BAP.ffds.tped

# Estimate kinship
emmax-kin -v -d 10 BAB.BAP.ffds

# EMMAX
for i in $(seq -w 1 12); do
for j in L R S W ; do
    emmax -v -d 10 -t BAB.BAP.Chr$i -p $j.txt -k BAB.BAP.ffds.BN.kinf -c BAB.BAP.ffds.PCA.txt -o BAB.BAP.$j.Chr$i &
done
done

# Prepare the output for manhattan plot
for i in L R S W ; do
    cat BAB.BAP.$i.Chr*.ps | cut -f1 | sed 's/_/\t/g' | sed 's/Chr0//g' | sed 's/Chr//g' | \
        paste - <(cat BAB.BAP.$i.Chr*.ps) | \
        awk '{print $1,$2,$3,$5}' | sed '1i\CHR BP SNP P' | sed 's/ /\t/g' \
        > BAB.BAP.$i.GWAS.txt
done

```

###### Manhattan plot

```r
library(qqman)
library(data.table)
BAB.BAP.L.GWAS <- fread("BAB.BAP.L.GWAS.txt", header = TRUE)
jpeg("BAB.BAP.L.GWAS.jpeg", height = 8, width = 16, units="cm", res = 300)
manhattan(BAB.BAP.L.GWAS, chr = "CHR", bp = "BP",
          snp = "SNP",col = c("orange", "darkgreen"),
          p = "P", ylim=c(0,10), logp = TRUE,cex=0.3, 
          cex.axis = 0.5,cex.main = 1, 
          suggestiveline = FALSE, genomewideline = -log10(1/3.19274e+06))
dev.off()

BAB.BAP.R.GWAS <- fread("BAB.BAP.R.GWAS.txt", header = TRUE)
jpeg("BAB.BAP.R.GWAS.jpeg", height = 8, width = 16, units="cm", res = 300)
manhattan(BAB.BAP.R.GWAS, chr = "CHR", bp = "BP",
          snp = "SNP",col = c("orange", "darkgreen"),
          p = "P", ylim=c(0,20), logp = TRUE,cex=0.3, 
          cex.axis = 0.5,cex.main = 1, 
          suggestiveline = FALSE, genomewideline = -log10(1/3.19274e+06))
dev.off()

BAB.BAP.S.GWAS <- fread("BAB.BAP.S.GWAS.txt", header = TRUE)
jpeg("BAB.BAP.S.GWAS.jpeg", height = 8, width = 16, units="cm", res = 300)
manhattan(BAB.BAP.S.GWAS, chr = "CHR", bp = "BP",
          snp = "SNP",col = c("orange", "darkgreen"),
          p = "P", ylim=c(0,10), logp = TRUE,cex=0.3, 
          cex.axis = 0.5,cex.main = 1, 
          suggestiveline = FALSE, genomewideline = -log10(1/3.19274e+06))
dev.off()

BAB.BAP.W.GWAS <- fread("BAB.BAP.W.GWAS.txt", header = TRUE)
jpeg("BAB.BAP.W.GWAS.jpeg", height = 8, width = 16, units="cm", res = 300)
manhattan(BAB.BAP.W.GWAS, chr = "CHR", bp = "BP",
          snp = "SNP",col = c("orange", "darkgreen"),
          p = "P", ylim=c(0,10), logp = TRUE,cex=0.3, 
          cex.axis = 0.5,cex.main = 1, 
          suggestiveline = FALSE, genomewideline = -log10(1/3.19274e+06))
dev.off()

```

#### CHN fruit shape

```bash
# Working directory
cd /data/zhaojiantao/pepper/vg/13GWAS/CHN

# Prepare the sample list
awk '{print $1,$1}' L.txt | sed 's/ /\t/g' > CHN.list

# Convert vcf to tped
for i in $(seq -w 1 12); do
    plink --vcf ../../04SNP/Miss_0.3/Chr$i.SNP.vcf.gz \
        --keep CHN.list \
        --allow-extra-chr \
        --recode \
        --double-id \
        --geno 0.7 \
        --maf 0.05 \
        --out CHN.Chr$i
    # Rename the chr id
    sed 's/CannZSG_//g' CHN.Chr$i.map > CHN.Chr$i.map2
    mv CHN.Chr$i.map2 CHN.Chr$i.map
    # Make bed
    plink --file CHN.Chr$i --make-bed --out CHN.Chr$i
done

# Calculate the effective number of SNPs
for i in $(seq -w 1 12); do
    java -jar /data/zhaojiantao/tools/GeneticTypeIError/gec.jar \
        --effect-number \
        --plink-binary CHN.Chr$i
    mv gec.sum CHN.Chr$i.gec.sum
done

# Calculate the effective number of SNPs of all chromosomes: 3.19274e+06
cat *gec.sum | grep -v Effective_Number | awk '{ sum += $2 } END { print sum }' 

# Population structure
plink --vcf <(zcat ../../064DTV/504.4DTV.Miss0.1.vcf.gz | grep -v Chr00) \
    --keep CHN.list \
    --allow-extra-chr \
    --recode \
    --double-id \
    --geno 0.7 \
    --maf 0.05 \
    --out CHN.ffds

# Calculate PCA
plink --file CHN.ffds --pca --out CHN.ffds --allow-extra-chr

# Remove unnecessary files
rm *log *nosex

# Prepare the PCA file as cofactor
awk '{print $1,$2,1,$3,$4,$5,$6,$7}' CHN.ffds.eigenvec > CHN.ffds.PCA.txt

# Transpose the ped file
for i in $(seq -w 1 12); do
    plink --file CHN.Chr$i --recode12 --output-missing-genotype 0 --transpose --out CHN.Chr$i
done

# Calculate kinship
plink --file CHN.ffds --recode12 --output-missing-genotype 0 --transpose --allow-extra-chr --out CHN.ffds
# Rename the chr id
cut -f1 -d ' ' CHN.ffds.tped | sed 's/CannZSG_Chr0//g' | sed 's/CannZSG_Chr//g' | \
    paste - CHN.ffds.tped | sed 's/\t/ /g' | cut -f1,3- -d ' ' > CHN.ffds.tped2
mv CHN.ffds.tped2 CHN.ffds.tped

# Estimate kinship
emmax-kin -v -d 10 CHN.ffds

# EMMAX
for i in $(seq -w 1 12); do
for j in L R S W ; do
    emmax -v -d 10 -t CHN.Chr$i -p $j.txt -k CHN.ffds.BN.kinf -c CHN.ffds.PCA.txt -o CHN.$j.Chr$i &
done
done

# Prepare the output for manhattan plot
for i in L R S W ; do
    cat CHN.$i.Chr*.ps | cut -f1 | sed 's/_/\t/g' | sed 's/Chr0//g' | sed 's/Chr//g' | \
        paste - <(cat CHN.$i.Chr*.ps) | \
        awk '{print $1,$2,$3,$5}' | sed '1i\CHR BP SNP P' | sed 's/ /\t/g' \
        > CHN.$i.GWAS.txt
done

# Candidate region
plink --file CHN.Chr07 \
    --r2 \
    --ld-snp Chr07_49272382 \
    --ld-window-kb 1000 \
    --ld-window 1000000 \
    --ld-window-r2 0.5 \
    --out Chr07_49272382

```

###### Manhattan plot

```r
library(qqman)
library(data.table)
CHN.L.GWAS <- fread("CHN.L.GWAS.txt", header = TRUE)
jpeg("CHN.L.GWAS.jpeg", height = 8, width = 16, units="cm", res = 300)
manhattan(CHN.L.GWAS, chr = "CHR", bp = "BP",
          snp = "SNP",col = c("orange", "darkgreen"),
          p = "P", ylim=c(0,10), logp = TRUE,cex=0.3, 
          cex.axis = 0.5,cex.main = 1, 
          suggestiveline = FALSE, genomewideline = -log10(1/2.13784e+06))
dev.off()


```

### SweeD

#### Installation

```bash
# Installation
cd /data/zhaojiantao/tools
git clone https://github.com/alachins/sweed.git
cd sweed
# Install
make -f Makefile.gcc
make -f Makefile.PTHREADS.gcc

# Install DMTCP library
DMTCPAWARELIB=/data/zhaojiantao/tools/sweed/dmtcp/lib/dmtcp

```

#### BAP with 10k grid

```bash
# Working directory
cd /data/zhaojiantao/pepper/vg/14SweeD/BAP10K

# SweeD does not accept zipped vcf
# Extract the samples
for i in $(seq -w 1 12); do
    bcftools view -S BAP.list \
        ../../04SNP/Miss_0.3/Chr$i.SNP.vcf.gz \
        > Chr$i.SNP.vcf &
done

# Prepare the bed file containing the chr and grid numbers
grep Chr ../../PanReference/Cann_Zhangshugang.chr.fa.fai | sed 's/CannZSG_//g' | \
    awk '{print $1,int($2/10000+0.5)}' \
    > Chr.grids.10K

# Run SweeD
IFS=$'\n'
for LINE in $(cat Chr.grids.10K); do
    Chr=$(echo ${LINE} | awk '{print $1}')
    grid=$(echo ${LINE} | awk '{print $2}')
    echo " SweeD-MPFR-P -name $Chr -input $Chr.SNP.vcf -grid $grid -maf 0.05 -missing 0.3 -threads 60 \
    -s 218 -eN 0.093 0.908 -eN 0.318 1.41 -eN 0.588 3.85 -eN 1.088 3.879 "
done > sweed.sh

# Add the chromosome ids to the output
sed -n '4,$p' SweeD_Report.Chr01 | awk '{print "Chr01_"$1,"1",$1,$4,$5,$2}' | sed 's/ /\t/g' > SweeD_Report.Chr01-1
sed -n '4,$p' SweeD_Report.Chr02 | awk '{print "Chr02_"$1,"2",$1,$4,$5,$2}' | sed 's/ /\t/g' > SweeD_Report.Chr02-1
sed -n '4,$p' SweeD_Report.Chr03 | awk '{print "Chr03_"$1,"3",$1,$4,$5,$2}' | sed 's/ /\t/g' > SweeD_Report.Chr03-1
sed -n '4,$p' SweeD_Report.Chr04 | awk '{print "Chr04_"$1,"4",$1,$4,$5,$2}' | sed 's/ /\t/g' > SweeD_Report.Chr04-1
sed -n '4,$p' SweeD_Report.Chr05 | awk '{print "Chr05_"$1,"5",$1,$4,$5,$2}' | sed 's/ /\t/g' > SweeD_Report.Chr05-1
sed -n '4,$p' SweeD_Report.Chr06 | awk '{print "Chr06_"$1,"6",$1,$4,$5,$2}' | sed 's/ /\t/g' > SweeD_Report.Chr06-1
sed -n '4,$p' SweeD_Report.Chr07 | awk '{print "Chr07_"$1,"7",$1,$4,$5,$2}' | sed 's/ /\t/g' > SweeD_Report.Chr07-1
sed -n '4,$p' SweeD_Report.Chr08 | awk '{print "Chr08_"$1,"8",$1,$4,$5,$2}' | sed 's/ /\t/g' > SweeD_Report.Chr08-1
sed -n '4,$p' SweeD_Report.Chr09 | awk '{print "Chr09_"$1,"9",$1,$4,$5,$2}' | sed 's/ /\t/g' > SweeD_Report.Chr09-1
sed -n '4,$p' SweeD_Report.Chr10 | awk '{print "Chr10_"$1,"10",$1,$4,$5,$2}' | sed 's/ /\t/g' > SweeD_Report.Chr10-1
sed -n '4,$p' SweeD_Report.Chr11 | awk '{print "Chr11_"$1,"11",$1,$4,$5,$2}' | sed 's/ /\t/g' > SweeD_Report.Chr11-1
sed -n '4,$p' SweeD_Report.Chr12 | awk '{print "Chr12_"$1,"12",$1,$4,$5,$2}' | sed 's/ /\t/g' > SweeD_Report.Chr12-1

# Combine the results
cat SweeD_Report.Chr*-1 | \
    sed '1i\SNP Chr Pos Start End SweeD' | sed 's/ /\t/g' \
    > SweeD_Report.BAP.10k.txt

# Reverse sort 
grep -v SNP SweeD_Report.BAP.10k.txt | awk '$6>0 {print}' | sort -k6gr | head -3021 > SweeD_Report.BAP.10k.top1.txt
grep -v SNP SweeD_Report.BAP.10k.txt | awk '$6>0 {print}' | sort -k6gr | head -15106 > SweeD_Report.BAP.10k.top5.txt
cut -f1 SweeD_Report.BAP.10k.top5.txt | sed 's/_/\t/g' | paste - SweeD_Report.BAP.10k.top5.txt | \
    awk '{print $3,$1,int($2+0.5),int($2+1+0.5),$8}' | sed '1i\SNP Chr Start End XPCLR' | sed 's/ /\t/g' \
    > SweeD_Report.BAP.10k.top5-1.txt

# Sort the file in excel and add header

# cutoff
grep -v SNP SweeD_Report.BAP.10k.txt | awk '$6>0 {print}' | sort -k6gr | sed -n '3021p'
grep -v SNP SweeD_Report.BAP.10k.txt | awk '$6>0 {print}' | sort -k6gr | sed -n '15106p'

# Merge the sweeps with distance <100000
perl merge.pl SweeD_Report.BAP.10k.top1-1.txt | sed '1d' > SweeD_Report.BAP.10k.top1.merged.txt
perl merge.pl SweeD_Report.BAP.10k.top5-1.txt | sed '1d' > SweeD_Report.BAP.10k.top5.merged.txt
```

#### BAP with 100k size

```bash
# Working directory
cd /data/zhaojiantao/pepper/vg/14SweeD/BAP100K

# Prepare the bed file containing the chr and grid numbers
grep Chr ../../PanReference/Cann_Zhangshugang.chr.fa.fai | sed 's/CannZSG_//g' | \
    awk '{print $1,int($2/100000+0.5)}' \
    > Chr.grids.100K

# Generate command lines
IFS=$'\n'
for LINE in $(cat Chr.grids.100K); do
    Chr=$(echo ${LINE} | awk '{print $1}')
    grid=$(echo ${LINE} | awk '{print $2}')
    echo " SweeD-MPFR-P -name $Chr -input ../BAP10K/$Chr.SNP.vcf -grid $grid -threads 60 -maf 0.05 -missing 0.3 \
    -s 218 -eN 0.093 0.908 -eN 0.318 1.41 -eN 0.588 3.85 -eN 1.088 3.879 "
done > sweed.100K.sh

# Add the chromosome ids to the output
sed -n '4,$p' SweeD_Report.Chr01 | awk '{print "Chr01_"$1,"1",$1,$4,$5,$2}' | sed 's/ /\t/g' > SweeD_Report.Chr01-1
sed -n '4,$p' SweeD_Report.Chr02 | awk '{print "Chr02_"$1,"2",$1,$4,$5,$2}' | sed 's/ /\t/g' > SweeD_Report.Chr02-1
sed -n '4,$p' SweeD_Report.Chr03 | awk '{print "Chr03_"$1,"3",$1,$4,$5,$2}' | sed 's/ /\t/g' > SweeD_Report.Chr03-1
sed -n '4,$p' SweeD_Report.Chr04 | awk '{print "Chr04_"$1,"4",$1,$4,$5,$2}' | sed 's/ /\t/g' > SweeD_Report.Chr04-1
sed -n '4,$p' SweeD_Report.Chr05 | awk '{print "Chr05_"$1,"5",$1,$4,$5,$2}' | sed 's/ /\t/g' > SweeD_Report.Chr05-1
sed -n '4,$p' SweeD_Report.Chr06 | awk '{print "Chr06_"$1,"6",$1,$4,$5,$2}' | sed 's/ /\t/g' > SweeD_Report.Chr06-1
sed -n '4,$p' SweeD_Report.Chr07 | awk '{print "Chr07_"$1,"7",$1,$4,$5,$2}' | sed 's/ /\t/g' > SweeD_Report.Chr07-1
sed -n '4,$p' SweeD_Report.Chr08 | awk '{print "Chr08_"$1,"8",$1,$4,$5,$2}' | sed 's/ /\t/g' > SweeD_Report.Chr08-1
sed -n '4,$p' SweeD_Report.Chr09 | awk '{print "Chr09_"$1,"9",$1,$4,$5,$2}' | sed 's/ /\t/g' > SweeD_Report.Chr09-1
sed -n '4,$p' SweeD_Report.Chr10 | awk '{print "Chr10_"$1,"10",$1,$4,$5,$2}' | sed 's/ /\t/g' > SweeD_Report.Chr10-1
sed -n '4,$p' SweeD_Report.Chr11 | awk '{print "Chr11_"$1,"11",$1,$4,$5,$2}' | sed 's/ /\t/g' > SweeD_Report.Chr11-1
sed -n '4,$p' SweeD_Report.Chr12 | awk '{print "Chr12_"$1,"12",$1,$4,$5,$2}' | sed 's/ /\t/g' > SweeD_Report.Chr12-1

# Combine the results
cat SweeD_Report.Chr*-1 | \
    sed '1i\SNP Chr Pos Start End SweeD' | sed 's/ /\t/g' \
    > SweeD_Report.BAP.100k.txt

# Reverse sort 
grep -v SNP SweeD_Report.BAP.100k.txt | awk '$6>0 {print}' | sort -k6gr | head -302 > SweeD_Report.BAP.100k.top1.txt
# Sort the file in excel and add header

# cutoff
grep -v SNP SweeD_Report.BAP.100k.txt | awk '$6>0 {print}' | sort -k6gr | sed -n '302p'

# Merge the sweeps with distance <100000
perl merge.pl SweeD_Report.BAP.100k.top1.txt | sed '1d' > SweeD_Report.BAP.100k.top1.merged.txt


```

#### ANN with 10k grid

```bash
# Working directory
cd /data/zhaojiantao/pepper/vg/14SweeD/ANN

# Extract the samples
for i in $(seq -w 1 12); do
    bcftools view -S ANN.list \
        ../../04SNP/Miss_0.3/Chr$i.SNP.vcf.gz \
        > Chr$i.SNP.vcf
done

# Prepare the bed file containing the chr and grid numbers
grep Chr ../../PanReference/Cann_Zhangshugang.chr.fa.fai | sed 's/CannZSG_//g' | \
    awk '{print $1,int($2/10000+0.5)}' \
    > Chr.grids.10K

# Run SweeD
IFS=$'\n'
for LINE in $(cat Chr.grids.10K); do
    Chr=$(echo ${LINE} | awk '{print $1}')
    grid=$(echo ${LINE} | awk '{print $2}')
    echo " SweeD-MPFR-P -name $Chr -input $Chr.SNP.vcf -grid $grid -threads 80 -maf 0.05 -missing 0.3 \
    -s 150 -eN 0.033 0.347 -eN 0.056 0.415 -eN 0.092 0.574 -eN 0.157 1.298 -eN 0.254 3.585 -eN 0.484 5.908 "
done > sweed.sh

# Add the chromosome ids to the output
sed -n '4,$p' SweeD_Report.Chr01 | awk '{print "Chr01_"$1,"1",$1,$4,$5,$2}' | sed 's/ /\t/g' > SweeD_Report.Chr01-1
sed -n '4,$p' SweeD_Report.Chr02 | awk '{print "Chr02_"$1,"2",$1,$4,$5,$2}' | sed 's/ /\t/g' > SweeD_Report.Chr02-1
sed -n '4,$p' SweeD_Report.Chr03 | awk '{print "Chr03_"$1,"3",$1,$4,$5,$2}' | sed 's/ /\t/g' > SweeD_Report.Chr03-1
sed -n '4,$p' SweeD_Report.Chr04 | awk '{print "Chr04_"$1,"4",$1,$4,$5,$2}' | sed 's/ /\t/g' > SweeD_Report.Chr04-1
sed -n '4,$p' SweeD_Report.Chr05 | awk '{print "Chr05_"$1,"5",$1,$4,$5,$2}' | sed 's/ /\t/g' > SweeD_Report.Chr05-1
sed -n '4,$p' SweeD_Report.Chr06 | awk '{print "Chr06_"$1,"6",$1,$4,$5,$2}' | sed 's/ /\t/g' > SweeD_Report.Chr06-1
sed -n '4,$p' SweeD_Report.Chr07 | awk '{print "Chr07_"$1,"7",$1,$4,$5,$2}' | sed 's/ /\t/g' > SweeD_Report.Chr07-1
sed -n '4,$p' SweeD_Report.Chr08 | awk '{print "Chr08_"$1,"8",$1,$4,$5,$2}' | sed 's/ /\t/g' > SweeD_Report.Chr08-1
sed -n '4,$p' SweeD_Report.Chr09 | awk '{print "Chr09_"$1,"9",$1,$4,$5,$2}' | sed 's/ /\t/g' > SweeD_Report.Chr09-1
sed -n '4,$p' SweeD_Report.Chr10 | awk '{print "Chr10_"$1,"10",$1,$4,$5,$2}' | sed 's/ /\t/g' > SweeD_Report.Chr10-1
sed -n '4,$p' SweeD_Report.Chr11 | awk '{print "Chr11_"$1,"11",$1,$4,$5,$2}' | sed 's/ /\t/g' > SweeD_Report.Chr11-1
sed -n '4,$p' SweeD_Report.Chr12 | awk '{print "Chr12_"$1,"12",$1,$4,$5,$2}' | sed 's/ /\t/g' > SweeD_Report.Chr12-1

# Combine the results
cat SweeD_Report.Chr*-1 | \
    sed '1i\SNP Chr Pos Start End SweeD' | sed 's/ /\t/g' \
    > SweeD_Report.ANN.10k.txt

# Reverse sort and top1% signals
sort -k6gr SweeD_Report.ANN.10k.txt | head -3021 > SweeD_Report.ANN.10k.top1.txt
sort -k6gr SweeD_Report.ANN.10k.txt | head -15106 > SweeD_Report.ANN.10k.top5.txt
cut -f1 SweeD_Report.ANN.10k.top5.txt | sed 's/_/\t/g' | paste - SweeD_Report.ANN.10k.top5.txt | \
    awk '{print $3,$1,int($2+0.5),int($2+1+0.5),$8}' | sed '1i\SNP Chr Start End XPCLR' | sed 's/ /\t/g' \
    > SweeD_Report.ANN.10k.top5-1.txt

# cutoff
sort -k6gr SweeD_Report.ANN.10k.txt | sed -n '3021p'
sort -k6gr SweeD_Report.ANN.10k.txt | sed -n '15106p'

# Merge the sweeps with distance <100000
perl merge.pl SweeD_Report.ANN.10k.top1-1.txt | sed '1d' > SweeD_Report.ANN.10k.top1.merged.txt
perl merge.pl SweeD_Report.ANN.10k.top5-1.txt | sed '1d' > SweeD_Report.ANN.10k.top5.merged.txt

```


#### ANN with 100k grid

```bash
# Working directory
cd /data/zhaojiantao/pepper/vg/14SweeD/ANN100K

# Prepare the bed file containing the chr and grid numbers
grep Chr ../../PanReference/Cann_Zhangshugang.chr.fa.fai | sed 's/CannZSG_//g' | \
    awk '{print $1,int($2/100000+0.5)}' \
    > Chr.grids.100K

# Generate command lines
IFS=$'\n'
for LINE in $(cat Chr.grids.100K); do
    Chr=$(echo ${LINE} | awk '{print $1}')
    grid=$(echo ${LINE} | awk '{print $2}')
    echo " SweeD-MPFR-P -name $Chr -input ../ANN10K/$Chr.SNP.vcf -grid $grid -threads 80 -maf 0.05 -missing 0.3 \
    -s 150 -eN 0.033 0.347 -eN 0.056 0.415 -eN 0.092 0.574 -eN 0.157 1.298 -eN 0.254 3.585 -eN 0.484 5.908  "
done > sweed.100K.sh

# Add the chromosome ids to the output
sed -n '4,$p' SweeD_Report.Chr01 | awk '{print "Chr01_"$1,"1",$1,$4,$5,$2}' | sed 's/ /\t/g' > SweeD_Report.Chr01-1
sed -n '4,$p' SweeD_Report.Chr02 | awk '{print "Chr02_"$1,"2",$1,$4,$5,$2}' | sed 's/ /\t/g' > SweeD_Report.Chr02-1
sed -n '4,$p' SweeD_Report.Chr03 | awk '{print "Chr03_"$1,"3",$1,$4,$5,$2}' | sed 's/ /\t/g' > SweeD_Report.Chr03-1
sed -n '4,$p' SweeD_Report.Chr04 | awk '{print "Chr04_"$1,"4",$1,$4,$5,$2}' | sed 's/ /\t/g' > SweeD_Report.Chr04-1
sed -n '4,$p' SweeD_Report.Chr05 | awk '{print "Chr05_"$1,"5",$1,$4,$5,$2}' | sed 's/ /\t/g' > SweeD_Report.Chr05-1
sed -n '4,$p' SweeD_Report.Chr06 | awk '{print "Chr06_"$1,"6",$1,$4,$5,$2}' | sed 's/ /\t/g' > SweeD_Report.Chr06-1
sed -n '4,$p' SweeD_Report.Chr07 | awk '{print "Chr07_"$1,"7",$1,$4,$5,$2}' | sed 's/ /\t/g' > SweeD_Report.Chr07-1
sed -n '4,$p' SweeD_Report.Chr08 | awk '{print "Chr08_"$1,"8",$1,$4,$5,$2}' | sed 's/ /\t/g' > SweeD_Report.Chr08-1
sed -n '4,$p' SweeD_Report.Chr09 | awk '{print "Chr09_"$1,"9",$1,$4,$5,$2}' | sed 's/ /\t/g' > SweeD_Report.Chr09-1
sed -n '4,$p' SweeD_Report.Chr10 | awk '{print "Chr10_"$1,"10",$1,$4,$5,$2}' | sed 's/ /\t/g' > SweeD_Report.Chr10-1
sed -n '4,$p' SweeD_Report.Chr11 | awk '{print "Chr11_"$1,"11",$1,$4,$5,$2}' | sed 's/ /\t/g' > SweeD_Report.Chr11-1
sed -n '4,$p' SweeD_Report.Chr12 | awk '{print "Chr12_"$1,"12",$1,$4,$5,$2}' | sed 's/ /\t/g' > SweeD_Report.Chr12-1

# Combine the results
cat SweeD_Report.Chr*-1 | \
    sed '1i\SNP Chr Pos Start End SweeD' | sed 's/ /\t/g' \
    > SweeD_Report.ANN.100k.txt

# Reverse sort 
grep -v SNP SweeD_Report.ANN.100k.txt | awk '$6>0 {print}' | sort -k6gr | head -302 > SweeD_Report.ANN.100k.top1.txt
# Sort the file in excel and add header

# cutoff
grep -v SNP SweeD_Report.ANN.100k.txt | awk '$6>0 {print}' | sort -k6gr | sed -n '302p'

# Merge the sweeps with distance <100000
perl merge.pl SweeD_Report.ANN.10k.top1.txt | sed '1d' > SweeD_Report.ANN.10k.top1.merged.txt

```

#### Overlaped regions

```bash
# Working directory
cd /data/zhaojiantao/pepper/vg/14SweeD/Overlapp

# Top GLA.ANN pi ratio
grep -v SNP Pi.GLA.ANN.txt  | head -28826 | cut -f1,7 | sed 's/_/\t/g' > Pi.GLA.ANN.top10.bed
grep -v SNP Pi.GLA.ANN.txt  | head -57652 | cut -f1,7 | sed 's/_/\t/g' > Pi.GLA.ANN.top20.bed
grep -v SNP Pi.GLA.ANN.txt  | head -86477 | cut -f1,7 | sed 's/_/\t/g' > Pi.GLA.ANN.top30.bed
grep -v SNP Pi.GLA.ANN.txt  | head -115303 | cut -f1,7 | sed 's/_/\t/g' > Pi.GLA.ANN.top40.bed
grep -v SNP Pi.GLA.ANN.txt  | head -144129 | cut -f1,7 | sed 's/_/\t/g' > Pi.GLA.ANN.top50.bed

# Top BAB.BAP pi ratio
grep -v SNP Pi.BAB.BAP.txt  | head -28822 | cut -f1,7 | sed 's/_/\t/g' > Pi.BAB.BAP.top10.bed
grep -v SNP Pi.BAB.BAP.txt  | head -57644 | cut -f1,7 | sed 's/_/\t/g' > Pi.BAB.BAP.top20.bed
grep -v SNP Pi.BAB.BAP.txt  | head -86465 | cut -f1,7 | sed 's/_/\t/g' > Pi.BAB.BAP.top30.bed
grep -v SNP Pi.BAB.BAP.txt  | head -115287 | cut -f1,7 | sed 's/_/\t/g' > Pi.BAB.BAP.top40.bed
grep -v SNP Pi.BAB.BAP.txt  | head -144109 | cut -f1,7 | sed 's/_/\t/g' > Pi.BAB.BAP.top50.bed

# ANN top 5% signal bed
cut -f1 ../ANN10K/SweeD_Report.ANN.10k.top5.txt | sed 's/_/\t/g' | \
    awk '{print $1,int($2)-5000,int($2+1)+5000,$1"_"$2}' | sed 's/ /\t/g' | grep -v SNP \
    > ANN.10k.top5.bed

# BAP top 5% signal bed
cut -f1 ../BAP10K/SweeD_Report.BAP.10k.top5.txt | sed 's/_/\t/g' | \
    awk '{print $1,int($2)-5000,int($2+1)+5000,$1"_"$2}' | sed 's/ /\t/g' | grep -v SNP \
    > BAP.10k.top5.bed

paste BAP.10k.top5.bed ../BAP10K/SweeD_Report.BAP.10k.top5.txt | \
    awk '{print $4,$1,$2,$3,$10}' | sed '1i\SNP Chr Start End XPCLR' | sed 's/ /\t/g' \
    > BAP.10k.top5.bed2
perl merge.pl BAP.10k.top5.bed2 | sed '1d' > BAP.10k.top5.bed2.merged

# Extract the overlapped sweed scores for ANN
for i in $(seq -w 10 10 50); do 
bedtools intersect -a ANN.10k.top5.bed -b Pi.GLA.ANN.top$i.bed -wo  > SweeD.ANN.overlapp.top$i.txt
done

# Extract the overlapped sweed scores for BAP
for i in $(seq -w 10 10 50); do 
bedtools intersect -a BAP.10k.top5.bed -b Pi.BAB.BAP.top$i.bed -wo  > SweeD.BAP.overlapp.top$i.txt
done

for i in $(seq -w 10 10 50); do 
cut -f4 SweeD.BAP.overlapp.top$i.txt | sort | uniq | wc -l
done

# Merge the overlapped pi ratio windows
for i in $(seq -w 10 10 50); do
cut -f5-8 SweeD.ANN.overlapp.top$i.txt | sort | uniq | awk '{print $1"_"$2,$1,$2,$3,$4}' | \
    sed '1i\SNP Chr Start End XPCLR' | sed 's/ /\t/g' \
    > SweeD.ANN.overlapp.top$i.uniq.txt
done

for i in $(seq -w 10 10 50); do
perl merge.pl SweeD.ANN.overlapp.top$i.uniq.txt | sed '1d' > SweeD.ANN.overlapp.top$i.uniq2.txt
done

# Calculate the total size
for i in $(seq -w 10 10 50); do
awk '{ sum += $3-$2+1 } END { print sum }' SweeD.ANN.overlapp.top$i.uniq2.txt
done

# Merge the overlapped pi ratio windows
for i in $(seq -w 10 10 50); do
cut -f5-8 SweeD.BAP.overlapp.top$i.txt | sort | uniq | awk '{print $1"_"$2,$1,$2,$3,$4}' | \
    sed '1i\SNP Chr Start End XPCLR' | sed 's/ /\t/g' \
    > SweeD.BAP.overlapp.top$i.uniq.txt
done

for i in $(seq -w 10 10 50); do
perl merge.pl SweeD.BAP.overlapp.top$i.uniq.txt | sed '1d' > SweeD.BAP.overlapp.top$i.uniq2.txt
done

# Calculate the total size
for i in $(seq -w 10 10 50); do
awk '{ sum += $3-$2+1 } END { print sum }' SweeD.BAP.overlapp.top$i.uniq2.txt
done

# BAP Genes
for i in $(seq -w 10 10 50); do
bedtools intersect \
    -a /data/zhaojiantao/pepper/References/S8/Zhangshugang.gff \
    -b SweeD.BAP.overlapp.top$i.uniq2.txt | \
    grep gene | sed 's/;/\t/g' | awk '{print $9}' | sed 's/ID=//g' | cut -f1 -d '.' | sort | uniq \
    > SweeD.BAP.overlapp.top$i.uniq2.genes.txt
done

# ANN Genes
for i in $(seq -w 10 10 50); do
bedtools intersect \
    -a /data/zhaojiantao/pepper/References/S8/Zhangshugang.gff \
    -b SweeD.ANN.overlapp.top$i.uniq2.txt | \
    grep gene | sed 's/;/\t/g' | awk '{print $9}' | sed 's/ID=//g' | cut -f1 -d '.' | sort | uniq \
    > SweeD.ANN.overlapp.top$i.uniq2.genes.txt
done

# Overlapped genes
for i in $(seq -w 10 10 50); do
comm -12 <(sort SweeD.ANN.overlapp.top$i.uniq2.genes.txt) \
    <(sort SweeD.BAP.overlapp.top$i.uniq2.genes.txt) \
    > SweeD.ANN.BAP.overlapp.top$i.uniq2.genes.txt
done

# Overlapped size
for i in $(seq -w 10 10 50); do
bedtools intersect \
    -a SweeD.ANN.overlapp.top$i.uniq2.txt \
    -b SweeD.BAP.overlapp.top$i.uniq2.txt \
    > SweeD.ANN.BAP.overlapp.top$i.uniq2.size.txt
done

# Calculate the total size
for i in $(seq -w 10 10 50); do
awk '{ sum += $3-$2+1 } END { print sum }' SweeD.ANN.BAP.overlapp.top$i.uniq2.size.txt
done

# Overlapped region between Sweed and pi ratio
bedtools intersect \
    -a ANN.10k.top5.bed \
    -b Pi.GLA.ANN.top10.bed -wo \
    > SweeD.ANN.overlapp.top10.txt

# Overlaped bed
bedtools intersect \
    -a SweeD.ANN.overlapp.top30.uniq2.txt \
    -b SweeD.BAP.overlapp.top30.uniq2.txt \
    > Overlapped.bed

# Extract the overlapped SNPs
# ANN
bedtools intersect \
    -a <(cut -f1 ../ANN10K/SweeD_Report.ANN.10k.txt | sed 's/_/\t/g' | \
        awk '{print $1,int($2),int($2+1),$1"_"$2}' | sed 's/ /\t/g' | grep -v SNP) \
    -b Overlapped.bed | cut -f4 | deal_table.pl -transpose \
    > SweeD.ANN.overlapp.SNPs.txt

# BAP
bedtools intersect \
    -a <(cut -f1 ../BAP10K/SweeD_Report.BAP.10k.txt | sed 's/_/\t/g' | \
        awk '{print $1,int($2),int($2+1),$1"_"$2}' | sed 's/ /\t/g' | grep -v SNP) \
    -b Overlapped.bed | cut -f4 | deal_table.pl -transpose \
    > SweeD.BAP.overlapp.SNPs.txt

# GLA.ANN
bedtools intersect \
    -a <(awk '{print $2,$3,$4,$5,$1}' ../GLA.ANN.ld8/SweeD_Report.GLA.ANN.10k.ld8.top1.txt | sed 's/ /\t/g' | grep -v Pos) \
    -b Overlapp.bed | cut -f5 | deal_table.pl -transpose \
    > GLA.ANN.overlapp.SNPs.txt

bedtools intersect \
    -a ../BAB.BAP.ld8/SweeD_Report.BAB.BAP.10k.ld8.top5%0.merged.txt \
    -b ../GLA.ANN.ld8/SweeD_Report.GLA.ANN.10k.ld8.top5%0.merged.txt \
    > top5%0.merged.txt

bedtools intersect \
    -a ../BAB.BAP.ld8/SweeD_Report.BAB.BAP.10k.ld8.top1.merged.txt \
    -b ../GLA.ANN.ld8/SweeD_Report.GLA.ANN.10k.ld8.top1.merged.txt \
    > top1.merged.txt

bedtools intersect \
    -a /data/zhaojiantao/pepper/References/S8/Zhangshugang.gff \
    -b top5%0.merged.txt | \
    grep gene | sed 's/;/\t/g' | awk '{print $9}' | sed 's/ID=//g' | cut -f1 -d '.' | sort | uniq \
    > top5%0.merged.genes.txt

bedtools intersect \
    -a /data/zhaojiantao/pepper/References/S8/Zhangshugang.gff \
    -b top1.merged.txt | \
    grep gene | sed 's/;/\t/g' | awk '{print $9}' | sed 's/ID=//g' | cut -f1 -d '.' | sort | uniq \
    > top1.merged.genes.txt

```

#### Comparison with Cao et al. 2022

```bash
# Working directory
cd /data/zhaojiantao/pepper/References/Zunla-1/v2.0/Cao

# Prepare the gene bed
cut -f1,4,5,9 ../Zunla-1.genes.gff | grep ID | sed 's/ID=//g' | sed 's/;//g' > Zunla-1.genes.bed
# Reformat the CDS to one row
perl -pe '/^>/ ? print "\n" : chomp' ../Zunla-1.CDS.fa > Zunla-1.CDS.fa

# Extract the CDS
for i in A B C D ; do
grep -f Stage$i.genes -A 1 Zunla-1.CDS.fa | sed '/^--$/d' > Stage$i.cds.fa
done

# Run blastn
for i in A B C D ; do
blastn -query Stage$i.cds.fa \
    -db /data/zhaojiantao/pepper/References/S8/blastdb/S8.cds.index \
    -evalue 1e-10 -num_threads 60 -outfmt 6 \
    -out Stage$i.cds.blastn6 
done

# Extract the best matches
for i in A B C D ; do
awk '!seen[$1]++' Stage$i.cds.blastn6 > Stage$i.cds.best.blastn6
done

## Challenge: Many of the genes could not find alignment

# Use the liftoff approach
liftoff -polish \
    -copies -sc 0.97 \
    -exclude_partial -s 0.97 -a 0.95 -p 80 \
    -dir intermediate.ncbi \
    -u unmapped.ncbi \
    -o Zunla.liftoff.gff3 \
    -g /data/zhaojiantao/pepper/References/S8/Zhangshugang.gff3 \
    /data/zhaojiantao/pepper/References/Zunla-1/v2.0/Zunla-1.fasta \
    /data/zhaojiantao/pepper/References/S8/Zhangshugang.fa

# Extrtact the remaining genes
grep -v "#" Zunla.liftoff.gff3_polished | awk '$3=="gene" {print $9}' | \
    cut -f1 -d ';' | sed 's/ID=//g' | cut -f1 -d '_' | sort | uniq | \
    grep -vf - /data/zhaojiantao/pepper/References/S8/Zhangshugang.gff3 \
    > Zhangshugang.remain.gff3

# Liftoff for the remaining genes
liftoff -polish \
    -copies -sc 0.97 \
    -exclude_partial -s 0.97 -a 0.95 -p 80 \
    -dir intermediate.ncbi \
    -u unmapped.ncbi \
    -o Zunla.remain.liftoff.gff3 \
    -g Zhangshugang.remain.gff3 \
    /data/zhaojiantao/pepper/References/Zunla-1/v2.0/Zunla-1.fasta \
    /data/zhaojiantao/pepper/References/S8/Zhangshugang.fa

# Combine the results and filter the infos
cat Zunla.liftoff.gff3_polished Zunla.remain.liftoff.gff3_polished | \
    grep -v "#" | awk '$3=="gene" {print}' | sed 's/;/\t/g' | \
    awk '{ for (i=1; i<=NF; i++) if ($i ~ /coverage|sequence|valid_ORFs/) {cols[i] = $i} } {printf "%s %s %s %s %s ", $9, $1, $4, $5, $7; for (i=6; i<=NF; i++) if (cols[i]) printf "%s ", $i; printf "\n" }' | \
    sed 's/ID=gene://g' | sed 's/coverage=//g' | sed 's/sequence_ID=//g' | sed 's/ID=//g' | \
    sed 's/valid_ORFs=0/Invalid/g' | sed 's/valid_ORFs=1/Valid/g' | sed 's/ /\t/g' \
    > Zunla.liftoff.genes

# Format Zhangshugang gene info
grep mRNA ../Zunla-1.genes.gff | sed 's/;/\t/g' | \
    awk '{print $9,$1,$4,$5,$7}' | sed 's/ID=//g' | sed 's/ /\t/g' \
    > Zunla.genes

# Find the overlapped genes
bedtools intersect -wo -a <(awk '{print $2,$3,$4,$1,$5}' Zunla.genes | sed 's/ /\t/g') \
    -b <(awk '{print $2,$3,$4,$1,$5,$6,$7,$8}' Zunla.liftoff.genes | sed 's/ /\t/g') | \
        sed '1i\Chr Start End Zunla Strand Chr Start End Zunla.liftoff Strand Coverage Identity ValidORFs OverlapSize' | \
        awk '{print $4,$1,$2,$3,$5,$9,$6,$7,$8,$10,$11,$12,$13,$14}' | sed 's/ /\t/g' \
        > Zunla_Zhangshugang.genes-1

# Remove those genes located on different chromosomes
cut -f6,7 Zunla_M82.genes-1 | grep -v Zunla | sed 's/g/\t/g' | cut -f1,3 | sed 's/Solyc/Chr/g' | \
    paste - <(cut -f6 Zunla_M82.genes-1 | grep -v Zunla) | \
    awk '{if ($1!=$2) print $3}' > Zunla_M82.genes-A
# Those genes from Chr00
cut -f6,7 Zunla_M82.genes-1 | grep -v Zunla | sed 's/g/\t/g' | cut -f1,3 | sed 's/Solyc/Chr/g' | \
    paste - <(cut -f6 Zunla_M82.genes-1 | grep -v Zunla) | \
    awk '{if ($1=="Chr00" || $2=="Chr00" ) print $3}' > Zunla_M82.genes-B

# Extract the final clean results        
grep -vf Zunla_M82.genes-B Zunla_M82.genes-A | \
        grep -vf - Zunla_M82.genes-1 \
        > Zunla_M82.genes-2

# Solyc gene infos from Zunla v4.0
grep -v "#" ../../Zunla1706/ITAG4.0_gene_models.gff | sed 's/;/ /g' | \
        awk ' $3=="gene" {print $9,$1,$4,$5,$7}' | sed 's/ID=gene://g' | \
        sed 's/ch/Chr/g' | sed 's/ /\t/g' \
        > Zunla.genes

# Extract the corresponding Zunla gene infos
for i in $(grep -v M82 Zunla_M82.genes-2 | cut -f6 | cut -f1 -d "_"); do grep $i Zunla.genes ; done | \
        sed '1i\Zunla Chr Start End Strand' | sed 's/ /\t/g' \
        > Zunla.genes-1

# Combine the final output
paste Zunla_M82.genes-2 Zunla.genes-1 > Zunla_M82.genes.txt

# For those with multiple copies, only keep the original one
cut -f6 Zunla_M82.genes.txt | grep -v Zunla.liftoff | sort | uniq -c | awk '$1>1 {print $2}' > Zunla_M82.genes.multi
grep -f Zunla_M82.genes.multi Zunla_M82.genes.txt

# Extract those Zunla genes could liftoff to M82 but with no overlapps with M82 gene annotations
grep -f <(cut -f6 Zunla_M82.genes.txt | grep -v Zunla.liftoff) Zunla.liftoff.genes \
        > Zunla.liftoff.uniq.genes 

# Extract the corresponding Zunla gene infos
for i in $(cut -f1 Zunla.liftoff.uniq.genes | cut -f1 -d "_"); do grep $i Zunla.genes ; done | \
        sed 's/ /\t/g' \
        > Zunla.genes-2

# Combine the Zunla liftoff to M82 with no overlapps with M82
paste Zunla.liftoff.uniq.genes Zunla.genes-2 | awk '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13}' | \
        sed '1i\Zunla.liftoff Chr Start End Strand Coverage Identity ORFs Zunla Chr Start End Strand' | \
        sed 's/ /\t/g' \
        > Zunla.liftoff.uniq.genes.txt

# Total Zunla v4 liftoff to M82
cat Zunla.genes-1 Zunla.genes-2  | grep -v Zunla | sort | uniq



```



### Introgression

```bash
# Working directory
cd /data/zhaojiantao/pepper/vg/15ABBA-BABA

# Prepare the Pop list file
awk '{print $1,"Ca.annuum"}' ../Pop.list/ANN.list > ANN.list
awk '{print $1,"Ca.pubescens"}' ../Pop.list/PUB.list > PUB.list
cat ../Pop.list/BAB.list ../Pop.list/BAP.list | awk '{print $1,"Ca.baccatum"}' > BAC.list
cat ../Pop.list/FRU.list ../Pop.list/CHN.list | awk '{print $1,"Ca.FRU.CHN"}' > FRU.CHN.list
cat *list > BAC2FRU.CHN.list
cut -f1 -d ' ' BAC2FRU.CHN.list > BAC2FRU.CHN.list.id

# Extract the genotypes
for i in $(seq -w 1 12); do
    bcftools view -S BAC2FRU.CHN.list.id -i 'F_MISSING<0.3' --min-af 0.05 \
        ../04SNP/Miss_0.3/Chr$i.SNP.vcf.gz \
        > Chr$i.SNP.vcf &
done

# Convert vcf to ABAB-ABBA geno file format
for i in $(seq -w 1 12); do
    python /data/zhaojiantao/tools/ABBA-BABA/genomics_general-master/VCF_processing/parseVCF.py \
        -i Chr$i.SNP.vcf --minQual 30 --gtf flag=DP \
        > Chr$i.clean.geno &
done

# Combine the files
head -1 Chr01.clean.geno > header
cat Chr*.clean.geno | grep -v "#CHROM" > clean.geno
cat header clean.geno | bgzip -@ 20 > All.clean.geno.gz

# Remove unnecessary files
rm header clean.geno Chr*.clean.geno

# ABBA-BABA
ABBABABAwindows.py \
    -g All.clean.geno.gz \
    -o PUB.BAC.FRUCHN.ANN.out \
    -w 500000 -s 50000 -m 100 -f phased --minData 0.5 -T 90 \
    -P1 Ca.annuum \
    -P2 Ca.FRU.CHN \
    -P3 Ca.baccatum \
    -O Ca.pubescens \
    --popsFile BAC2FRU.CHN.list \
    --writeFailedWindows

# Select the top 5% highest fd values
cut -f1 PUB.BAC.FRUCHN.ANN.boxplot.out | sed 's/Chr0//g' | sed 's/Chr//g' | \
    paste - PUB.BAC.FRUCHN.ANN.boxplot.out | cut -f1,3- | sed '1i\Chr Start End Pos SNP fd' | sed 's/ /\t/g' \
    > PUB.BAC.FRUCHN.ANN.boxplot.vg.out.txt

```

### Introgression with LD prune

```bash
# Working directory
cd /data/zhaojiantao/pepper/vg/15ABBA-BABA/LDprune

# LD prune
for i in $(seq -w 1 12); do
    # LD estimation
    plink --gzvcf ../Chr$i.SNP.vcf \
        --indep-pairwise 50 10 0.8 \
        --allow-extra-chr --double-id \
        --out Chr$i
    # Extract the prunned SNPs
    bcftools view --include ID==@Chr$i.prune.in ../Chr$i.SNP.vcf > Chr$i.prune.vcf
done

# Convert vcf to ABAB-ABBA geno file format
for i in $(seq -w 1 12); do
    python /data/zhaojiantao/tools/ABBA-BABA/genomics_general-master/VCF_processing/parseVCF.py \
        -i Chr$i.prune.vcf --minQual 30 --gtf flag=DP \
        > Chr$i.clean.geno
done

# Combine the files
head -1 Chr01.clean.geno > header
cat Chr*.clean.geno | grep -v "#CHROM" > clean.geno
cat header clean.geno | bgzip -@ 20 > All.clean.geno.gz

# Remove unnecessary files
rm header clean.geno Chr*.clean.geno

# ABBA-BABA
ABBABABAwindows.py \
    -g All.clean.geno.gz \
    -o PUB.BAC.FRUCHN.ANN.out \
    -w 500000 -s 50000 -m 100 -f phased --minData 0.5 -T 90 \
    -P1 Ca.annuum \
    -P2 Ca.FRU.CHN \
    -P3 Ca.baccatum \
    -O Ca.pubescens \
    --popsFile BAC2FRU.CHN.list \
    --writeFailedWindows

# Select the top 5% highest fd values
cut -f1 PUB.BAC.FRUCHN.ANN.boxplot.out | sed 's/Chr0//g' | sed 's/Chr//g' | \
    paste - PUB.BAC.FRUCHN.ANN.boxplot.out | cut -f1,3- | sed '1i\Chr Start End Pos SNP fd' | sed 's/ /\t/g' \
    > PUB.BAC.FRUCHN.ANN.boxplot.vg.out.ld.txt

# Calculate the average value of fd values.
cat PUB.BAC.FRUCHN.ANN.out | sed 's/,/ /g' | awk '$9>0 && $10<1{print $0}' | grep -v nan  | awk '{s+=$10} END {print s/NR}'
# Calculate the total value of fd
cat PUB.BAC.FRUCHN.ANN.out | sed 's/,/ /g' | awk '$9>0 && $10<1{print $0}' | grep -v nan  | awk '{s+=$6*$10} END {print s}' 
# Calculate the total G
cat PUB.BAC.FRUCHN.ANN.out | sed 's/,/ /g' | awk '$9>0 && $10<1{print $0}' | grep -v nan  | awk '{s+=$6} END {print s}' 


```




## Mapping quality evaluation

### Mapping quality of the 1314 collection

```bash
# Working directory
cd /data/zhaojiantao/pepper/pepper_survey/01mapping

# Trimm the raw sequences
for i in $(cat sample.id); do
    # Remove adaptors and low-quality sequences
    trimmomatic PE ../00fastq/$i.R1.fastq.gz ../00fastq/$i.R2.fastq.gz \
        $i.R1P.fq $i.R1U.fq \
        $i.R2P.fq $i.R2U.fq \
        -threads 40 \
        ILLUMINACLIP:/data/zhaojiantao/tools/miniconda3/share/trimmomatic-0.39-2/adapters/TruSeq3-PE-2.fa:2:30:10:1:TRUE \
        SLIDINGWINDOW:4:20 LEADING:3 TRAILING:3 MINLEN:40 ;
done

# Working directory
cd /data/zhaojiantao/pepper/pepper_survey/02bam

# Map the reads to S8
for i in $(cat ../01mapping/sample.id); do    
    sentieon bwa mem -M -R "@RG\tLB:Library1\tID:FeiLab_${i}\tSM:$i\tPL:ILLUMINA" \
        -t 80 /data/zhaojiantao/pepper/References/S8/Zhangshugang.fasta \
        ../01mapping/$i.R1P.fq ../01mapping/$i.R2P.fq | \
        sentieon util sort -o $i.bam -t 80 --sam2bam -i -
done

# Raw reads pairs
grep "Input Read Pairs" trimm2.log | cut -f4 -d ' '

# Cleaned reads pairs
grep "Input Read Pairs" trimm2.log | cut -f7 -d ' '

# Count the number of paired reads from fastq
echo $(cat yourfile.fastq|wc -l)/4|bc

# Total mapped reads
for i in $(cat ../01mapping/sample2.id); do
    sentieon driver -t 20 \
        -r /data/zhaojiantao/pepper/References/S8/Zhangshugang.fasta \
        -i $i.bam \
        --algo AlignmentStat \
        $i.stat
done

# Total cleaned reads
for i in $(ls *stat); do cut -f2 $i | sed -n '5p' ; done > AA

# Total aligned reads
for i in $(ls *stat); do cut -f6 $i | sed -n '5p' ; done > BB

samtools flagstat bailajiao.bam > bailajiao.flagstat

# Total reads
for i in $(ll *flagstat | awk '{print $9}' | sed 's/.flagstat//g'); do sed -n '1p' $i.flagstat | awk '{print $1}'; done
# Total mapped Reads
for i in $(ll *flagstat | awk '{print $9}' | sed 's/.flagstat//g'); do sed -n '7p' $i.flagstat | awk '{print $1}'; done

# Unmapped singletons
for i in $(ll *flagstat | awk '{print $9}' | sed 's/.flagstat//g'); do sed -n '11p' $i.flagstat | awk '{print $1}'; done
# Total paired Reads
for i in $(ll *flagstat | awk '{print $9}' | sed 's/.flagstat//g'); do sed -n '7p' $i.flagstat | awk '{print $1}'; done
# Total reads + 0 paired in sequencing
for i in $(ll *flagstat | awk '{print $9}' | sed 's/.flagstat//g'); do sed -n '6p' $i.flagstat | awk '{print $1}'; done

# Calculate genome coverage
for i in $(ll *flagstat | awk '{print $9}' | sed 's/.flagstat//g'); do
    bedtools genomecov -ibam $i.bam -bg > $i.bed 
    awk '$4>=1 {print $3-$2+1}' $i.bed | awk '{sum+=$1} END {print 1,sum}' > $i.1.bases
    awk '$4>=2 {print $3-$2+1}' $i.bed | awk '{sum+=$1} END {print 2,sum}' > $i.2.bases
    awk '$4>=5 {print $3-$2+1}' $i.bed | awk '{sum+=$1} END {print 5,sum}' > $i.5.bases
    cat $i.1.bases $i.2.bases $i.5.bases > $i.mapping.bases
    rm $i.1.bases $i.2.bases $i.5.bases $i.bed $i.bam $i.bam.bai
done

# Clean the files
# Remove those raw files with bam file generated
for i in $(ls *bai | sed 's/.bam.bai//g'); do rm ../00fastq/$i.R1.fastq.gz ../00fastq/$i.R2.fastq.gz ; done

# Remove the cleaned paired fastq files with bam file generated
for i in $(ls *bai | sed 's/.bam.bai//g'); do rm ../01mapping/$i.R1P.fq ../01mapping/$i.R2P.fq ; done

```

### Mapping quality of the 500 core collection

```bash
# Change working directory
cd /data/zhaojiantao/pepper/pepper_mapping_ratio

for i in $(cat Sample.id); do
    # Map the reads to S8
    sentieon bwa mem -M -R "@RG\tLB:Library1\tID:FeiLab_${i}\tSM:$i\tPL:ILLUMINA" \
        -t 100 /data/zhaojiantao/pepper/References/S8/pepper_S8.fasta \
        /data/jiaochen/pepper_reseq/clean_reads/combined/$i.R1.comb.fastq.gz \
        /data/jiaochen/pepper_reseq/clean_reads/combined/$i.R2.comb.fastq.gz | \
        sentieon util sort -o $i.bam -t 100 --sam2bam -i -
done

# Calculate mapped ratios
samtools mpileup -f /data/zhaojiantao/pepper/References/S8/pepper_S8.fasta PI_406847.bam -d 100 -o PI_406847.pileup
# Base with mapped reads > 1
awk '$4>1 {print $0}' PI_406847.pileup | wc -l
# Base with mapped reads > 2
awk '$4>2 {print $0}' PI_406847.pileup | wc -l
# Total bases
cat PI_406847.pileup | wc -l

# Total mapped reads
samtools flagstat PI_439367.bam > PI_439367.flagstat

# Total reads
for i in $(ll *flagstat | awk '{print $9}' | sed 's/.flagstat//g'); do sed -n '1p' $i.flagstat | awk '{print $1}'; done
# Total mapped Reads
for i in $(ll *flagstat | awk '{print $9}' | sed 's/.flagstat//g'); do sed -n '5p' $i.flagstat | awk '{print $1}'; done
# Unmapped singletons
for i in $(ll *flagstat | awk '{print $9}' | sed 's/.flagstat//g'); do sed -n '11p' $i.flagstat | awk '{print $1}'; done
# Total paired Reads
for i in $(ll *flagstat | awk '{print $9}' | sed 's/.flagstat//g'); do sed -n '7p' $i.flagstat | awk '{print $1}'; done
# Total reads + 0 paired in sequencing
for i in $(ll *flagstat | awk '{print $9}' | sed 's/.flagstat//g'); do sed -n '6p' $i.flagstat | awk '{print $1}'; done
# Calculate genome coverage

for i in $(ll *flagstat | awk '{print $9}' | sed 's/.flagstat//g'); do
    bedtools genomecov -ibam $i.bam -bg > $i.bed 
    awk '$4>=1 {print $3-$2}' $i.bed | awk '{sum+=$1} END {print 1,sum}' > $i.1.bases
    awk '$4>=2 {print $3-$2}' $i.bed | awk '{sum+=$1} END {print 2,sum}' > $i.2.bases
    awk '$4>=5 {print $3-$2}' $i.bed | awk '{sum+=$1} END {print 5,sum}' > $i.5.bases
    cat $i.1.bases $i.2.bases $i.5.bases > $i.mapping.bases
    rm $i.1.bases $i.2.bases $i.5.bases $i.bed $i.bam $i.bam.bai
done

```

### Genotype check of C. baccatum

```bash
# Working directory
cd /data/zhaojiantao/pepper/Ca.baccatum.check

# Index the genome
bwa index /data/zhaojiantao/pepper/References/PBC81/PBC81.fasta
samtools faidx /data/zhaojiantao/pepper/References/PBC81/PBC81.fasta
java -jar /data/zhaojiantao/tools/picard/build/libs/picard.jar CreateSequenceDictionary R=PBC81.fasta O=PBC81.dict

# Create the FASTA sequence directory file of the reference to be used in gatk
gatk CreateSequencedirectory -R /data/zhaojiantao/pepper/References/PBC81/PBC81.fasta

# Call variants
for i in Grif_9354 Grif_9355 PI_257110 PI_260564 PI_260595 PI_281308 PI_281310 PI_281436 PI_337524 PI_370010; do
    # Remove adaptors and low-quality sequences
    trimmomatic PE ../00fastq/$i.R1.fastq.gz ../00fastq/$i.R2.fastq.gz \
        $i.R1P.fq $i.R1U.fq \
        $i.R2P.fq $i.R2U.fq \
        -threads 80 \
        ILLUMINACLIP:/data/zhaojiantao/tools/miniconda3/share/trimmomatic-0.39-2/adapters/TruSeq3-PE-2.fa:2:30:10:1:TRUE \
        SLIDINGWINDOW:4:20 LEADING:3 TRAILING:3 MINLEN:40 ;
    # Map the reads to S8
    sentieon bwa mem -M -R "@RG\tLB:Library1\tID:FeiLab_${i}\tSM:$i\tPL:ILLUMINA" \
        -t 80 /data/zhaojiantao/pepper/References/PBC81/PBC81.fasta \
        $i.R1P.fq $i.R2P.fq | \
        sentieon util sort -o $i.bam -t 80 --sam2bam -i -
    # Collect read information to remove or mark duplicates
    sentieon driver -t 80 -i $i.bam \
        --algo LocusCollector --fun score_info $i.score.gz
    # Remove or mark duplicates
    sentieon driver -t 80 -i $i.bam \
        --algo Dedup --score_info $i.score.gz \
        --metrics $i.dedup_metrix.txt $i.deduped.bam
    # Call gvcf
    sentieon driver -t 80 -r /data/zhaojiantao/pepper/References/PBC81/PBC81.fasta \
        -i $i.deduped.bam \
        --algo Haplotyper \
        --emit_mode gvcf \
        $i.gvcf.gz
    # Remove unnecessary files
    rm $i.R1U.fq $i.R2U.fq *.score.gz* *dedup_metrix.txt *R1P.fq *R2P.fq
done

# Joint variant call
sentieon driver -t 80 -r /data/zhaojiantao/pepper/References/PBC81/PBC81.fasta --algo GVCFtyper \
    Ca.baccatum10.vcf.gz *gvcf.gz 

# Extract the SNPs could be mapped in both references
# BAC SNPs
grep -v "NA" rand500k.map_loc.REF_BACslct | cut -f5,6 > BAC.snps
# Extract the genotypes
vcftools --gzvcf Ca.baccatum10.vcf.gz \
    --positions BAC.snps \
    --recode --recode-INFO-all \
    --out Ca.baccatum10

# ANN SNPs
grep -v "NA" rand500k.map_loc.REF_BACslct | cut -f2,3 > ANN.snps

# Extract the genotypes
for i in $(seq -w 0 12); do
    vcftools --gzvcf /data/zhaojiantao/pepper/vcf.clean/Chr$i.clean.vcf.gz \
        --keep Ca.baccatum10.list \
        --positions ANN.snps \
        --recode --recode-INFO-all \
        --out Ca.annuum10.Chr$i &
done

# Combine the vcfs for the 500 core collection
grep "#" Ca.annuum10.Chr00.recode.vcf > header
cat Ca.annuum10.Chr*.recode.vcf | grep -v "#" | cat header - > Ca.annuum10.recode.vcf

# Extract biallelic SNPs
bcftools view --types snps -m 2 -M 2 Ca.baccatum10.comm.recode.vcf.gz -O v -o Ca.baccatum10.comm.biallelic.vcf.gz

zcat Ca.baccatum10.comm.biallelic.vcf.gz | grep -v "#" | cut -f1,2 | head -50000 

# Find the common SNPs
zcat Ca.baccatum10.recode.vcf.gz | grep -v "#" | cut -f1,2 | \
    grep -Fwf - rand500k.map_loc.REF_BACslct.filter \
    > Ca.baccatum10.comm1.snps.txt

zcat Ca.baccatum10.biallelic.vcf.gz | grep -v "#" | cut -f1,2 | \
    grep -Fwf - rand500k.map_loc.REF_BACslct.filter | \
    grep -Fwf <(grep -v "#" Ca.annuum10.recode.vcf | cut -f1,2) - \
    > Ca.annuum10.comm1.snps.txt

# Use the top 50k SNPs
head -50000 Ca.baccatum10.comm1.snps.txt > Ca.baccatum10.comm50k.snps.txt
grep -Fwf <(zcat Ca.baccatum10.comm.biallelic.vcf.gz | grep -v "#" | cut -f1,2 | \
    grep -v 43775082 | grep -v 37665850 | grep -v 40250533 | \
    grep -v 151351965 | grep -v 142012760 | grep -v 62417332 | grep -v 139406871) \
    Ca.baccatum10.comm1.snps.txt | head -50000 > Ca.baccatum10.comm50k.snps.txt

# Extract the common genotypes for C. baccatum
vcftools --gzvcf Ca.baccatum10.comm.recode.vcf.gz \
    --positions <(cut -f5,6 Ca.baccatum10.comm50k.snps.txt) \
    --recode --recode-INFO-all \
    --out Ca.baccatum10.comm50k

# Extract the common genotypes for C. annuum
vcftools --gzvcf Ca.annuum10.comm.recode.vcf.gz \
    --positions <(cut -f2,3 Ca.baccatum10.comm50k.snps.txt) \
    --recode --recode-INFO-all \
    --out Ca.annuum10.comm50k

# Prepare the file for genotype consistency calculation
plink --vcf Ca.annuum10.comm50k.recode.vcf --recode --make-bed --allow-extra-chr --out Ca.annuum10.comm50k
plink --vcf Ca.baccatum10.comm50k.recode.vcf --recode --make-bed --allow-extra-chr --out Ca.baccatum10.comm50k

# Replace the C. baccatum SNP position with correspoding position in C. annuum
awk '{print "Chr0"$1,$2,$3,$4}' Ca.baccatum10.comm50k.map | \
    sed 's/Chr010/Chr10/g' | sed 's/Chr011/Chr11/g' | \
    sed 's/Chr012/Chr12/g' | sed 's/Chr0scaffold/scaffold/g' | sed 's/ /\t/g' \
    > Ca.baccatum10.comm50k2.map

grep -Fwf <(cut -f1,4 Ca.baccatum10.comm50k2.map) \
    Ca.baccatum10.comm50k.snps.txt | sort -k5,5 -k6n \
    > Ca.baccatum10.comm50k3.map

# Cross check the consistency of ids in excel
paste Ca.baccatum10.comm50k2.map Ca.baccatum10.comm50k3.map | \
    awk '{if ($4=$10) print "TRUE",$0 ; else print "FALSE",$0}'  | grep FALSE

awk '{print $2,".",0,$3}' Ca.baccatum10.comm50k3.map | sed 's/ /\t/g' > Ca.baccatum10.comm.annuum.map

cp Ca.baccatum10.comm50k.ped Ca.baccatum10.comm.annuum.ped

# Convert ped to tped
plink --file Ca.annuum10.comm50k --recode transpose --out Ca.annuum10.comm50k --allow-extra-chr
plink --file Ca.baccatum10.comm.annuum --recode transpose --out Ca.baccatum10.comm.annuum --allow-extra-chr

# 互补配对
grep -Fwf <(cut -f1 Ca.baccatum10.comm50k.snps.txt) <(zcat rand500k.map_loc.REF_BACslct.filter.gz) | \
    awk '{if ($4!=$8) print $2"_"$3, "Negative"; else print $2"_"$3, "Positive"}' | \
    grep Negative | cut -f1 -d ' ' \
    > negative.snps

grep -Fwf <(cut -f1 Ca.baccatum10.comm50k.snps.txt) <(zcat rand500k.map_loc.REF_BACslct.filter.gz) | \
    awk '{if ($4!=$8) print $2"_"$3, "Negative"; else print $2"_"$3, "Positive"}' | \
    grep Positive | cut -f1 -d ' ' \
    > positive.snps

# Negative genotype
grep -Fwf negative.snps <(awk '{print $1,"Chr0"$1"_"$4,$0}' Ca.baccatum10.comm.annuum.tped | cut -f1,2,5- -d ' ') | cut -f5- -d ' ' | \
    sed 's/A/1/g' | sed 's/T/2/g' | sed 's/C/3/g' | sed 's/G/4/g' | \
    sed 's/1/T/g' | sed 's/2/A/g' | sed 's/3/G/g' | sed 's/4/C/g' | \
    paste <(grep -Fwf negative.snps <(awk '{print $1,"Chr0"$1"_"$4,$0}' Ca.baccatum10.comm.annuum.tped | cut -f1,2,5- -d ' ') | cut -f1-4 -d ' ' ) - | sed 's/\t/ /g' \
    > Ca.baccatum10.comm.annuum.tped.negative

grep -Fvf negative.snps <(awk '{print $1,"Chr0"$1"_"$4,$0}' Ca.baccatum10.comm.annuum.tped | cut -f1,2,5- -d ' ') \
    > Ca.baccatum10.comm.annuum.tped.positive

cat Ca.baccatum10.comm.annuum.tped.negative Ca.baccatum10.comm.annuum.tped.positive | \
    sort -k1n -k4n > Ca.baccatum10.comm.annuum.tped2

# Calculate the number of consistency
cut -f5- -d ' ' Ca.baccatum10.comm.annuum.tped2 | sed 's! \([^ ]\+\)\( \|$\)!\1 !g' > Ca.baccatum10.comm.annuum.tped2-1
cut -f5- -d ' ' Ca.annuum10.comm50k.tped | sed 's! \([^ ]\+\)\( \|$\)!\1 !g' > Ca.annuum10.comm50k.tped-1
cut -f5- -d ' ' Ca.baccatum10.comm.annuum.tped2 > Ca.baccatum10.comm.annuum.tped2-2
cut -f5- -d ' ' Ca.annuum10.comm50k.tped > Ca.annuum10.comm50k.tped-2


# missing
paste Ca.annuum10.comm50k.tped-2 Ca.baccatum10.comm.annuum.tped2-2 | \
    awk '{print $1,$2,$21,$22}' | grep 0 | wc -l

# Consistency of both allele for sample 1
paste Ca.annuum10.comm50k.tped-2 Ca.baccatum10.comm.annuum.tped2-2 | \
    awk '{print $1,$2,$21,$22}' | grep -v 0 |  sort | uniq -c \
    > count.sample1

paste Ca.annuum10.comm50k.tped-2 Ca.baccatum10.comm.annuum.tped2-2 | awk '{print $1,$2,$21,$22}' | grep -v 0 |  sort | uniq -c > count.sample1

# Consistency of both alleles
# Ref=Ref; Alt=Alt
awk '{if ($2==$4 && $3==$5) print $0,"TRUE"}' count.sample1 | awk '{sum+=$1} END {print sum}' 
# Ref=Alt; Alt=Ref
awk '{if ($2!=$3 && $2==$5 && $3==$4) print $0,"TRUE"}' count.sample1 | awk '{sum+=$1} END {print sum}' 

# AA==TT
grep "A A T T" count.sample1 | awk '{sum+=$1} END {print sum}' 

# TT==AA
grep "T T A A" count.sample1 | awk '{sum+=$1} END {print sum}' 

# CC==GG
grep "C C G G" count.sample1 | awk '{sum+=$1} END {print sum}' 

# GG==CC
grep "G G C C" count.sample1 | awk '{sum+=$1} END {print sum}' 


# Consistency of one allele with Zhangshugang homo
paste Ca.annuum10.comm50k.tped-2 Ca.baccatum10.comm.annuum.tped2-2 | \
    awk '{print $1,$2,$21,$22}' | grep -v 0 | \
    awk '{if ($1==$3 && $2==$4) print $0,"TRUE1"; else if ($1==$4 && $2==$3) print $0,"TRUE2"; else if ($1==$3 || $1==$4 || $2==$3 || $2==$4) print $0,"HALF" ; else print $0,"FALSE"}' | grep HALF | \
    awk '{if ($1==$2)print $0,"Zhangshugang.home.TRUE"; else print $0,"Zhangshugang.home.FALSE"}' | \
    cut -f6 -d ' ' | sort | uniq -c

```

## SNP calling and evaluation

### SNPs calling and quality control of the 1300 pepper collection

```bash
# working directory
cd /data/zhaojiantao/pepper/pepper_survey/01mapping

# Obtain all the sample ids
ll *R1.fastq.gz | awk '{print $9}' | sed 's/.R1.fastq.gz//g' > sample.id

# Index the genome
bwa index /data/zhaojiantao/pepper/References/S8/pepper_S8.fasta
samtools faidx /data/zhaojiantao/pepper/References/S8/pepper_S8.fasta
java -jar /data/zhaojiantao/tools/picard/build/libs/picard.jar CreateSequenceDictionary R=pepper_S8.fasta O=pepper_S8.dict

# Create the FASTA sequence directory file of the reference to be used in gatk
gatk CreateSequencedirectory -R /data/zhaojiantao/pepper/References/S8/pepper_S8.fasta

# Call variants
for i in $(cat sample.id); do
    # Remove adaptors and low-quality sequences
    trimmomatic PE ../00fastq/$i.R1.fastq.gz ../00fastq/$i.R2.fastq.gz \
        $i.R1P.fq $i.R1U.fq \
        $i.R2P.fq $i.R2U.fq \
        -threads 80 \
        ILLUMINACLIP:/data/zhaojiantao/tools/miniconda3/share/trimmomatic-0.39-2/adapters/TruSeq3-PE-2.fa:2:30:10:1:TRUE \
        SLIDINGWINDOW:4:20 LEADING:3 TRAILING:3 MINLEN:40 ;
    # Map the reads to S8
    sentieon bwa mem -M -R "@RG\tLB:Library1\tID:FeiLab_${i}\tSM:$i\tPL:ILLUMINA" \
        -t 80 /data/zhaojiantao/pepper/References/S8/pepper_S8.fasta \
        $i.R1P.fq $i.R2P.fq | \
        sentieon util sort -o $i.bam -t 80 --sam2bam -i -
    # Collect read information to remove or mark duplicates
    sentieon driver -t 80 -i $i.bam \
        --algo LocusCollector --fun score_info $i.score.gz
    # Remove or mark duplicates
    sentieon driver -t 80 -i $i.bam \
        --algo Dedup --score_info $i.score.gz \
        --metrics $i.dedup_metrix.txt $i.deduped.bam
    # Call gvcf
    sentieon driver -t 80 -r /data/zhaojiantao/pepper/References/S8/pepper_S8.fasta \
        -i $i.deduped.bam \
        --algo Haplotyper \
        --emit_mode gvcf \
        $i.gvcf.gz
    # Remove unnecessary files
    rm $i.R1U.fq $i.R2U.fq *bam* *.score.gz* *dedup_metrix.txt *R1P.fq *R2P.fq
done

# Joint variant call
sentieon driver -t 80 -r /data/zhaojiantao/pepper/References/S8/pepper_S8.fasta --algo GVCFtyper \
    1314.vcf.gz *gvcf.gz 

```

### SNP calling quality check

```bash
# Working directory
cd /data/zhaojiantao/pepper/pepper_survey

# For 1314 collection
vcftools --gzvcf All.filtered.4dtv.thin.vcf.gz \
    --keep ../Pop.list/500.from.1314.list \
    --recode --recode-INFO-all \
    --out All.4dtv.1314

# For 500 collection
vcftools --gzvcf ../4DTV/ffds.vcf.gz \
    --positions <(zcat All.filtered.4dtv.thin.vcf.gz | grep -v "#" | cut -f1,2) \
    --recode --recode-INFO-all \
    --out All.4dtv.500

# Create chromosome map
bcftools view -H All.4dtv.1314.recode.vcf | \
    cut -f 1 | uniq | awk '{print $0"\t"$0}' > \
    All.4dtv.1314.chrom-map.txt

bcftools view -H All.4dtv.500.recode.vcf | \
    cut -f 1 | uniq | awk '{print $0"\t"$0}' > \
    All.4dtv.500.chrom-map.txt

# Make ped file using this chromosome map
vcftools --vcf All.4dtv.1314.recode.vcf \
    --plink --chrom-map All.4dtv.1314.chrom-map.txt \
    --out All.4dtv.1314

vcftools --vcf All.4dtv.500.recode.vcf \
    --plink --chrom-map All.4dtv.500.chrom-map.txt \
    --out All.4dtv.500

# The number of Nomissing values
for i in $(seq 4 503); do 
    cut -f1-3,$i 1314.txt | awk '$4!=00 {print $0}' | wc ; 
done | awk '{print $1}'

# The number of overlapped Nomissing values
for i in $(seq 4 503); do \
    comm -12 <(cut -f1-3,$i 1314.txt | \awk '$4!=00 {print $0}' | sort) \
    <(cut -f1-3,$i 500.txt | awk '$4!=00 {print $0}' | sort) | wc ; 
done | awk '{print $1}'

```

### SNP count for each group

```bash
# Working directory
cd /data/zhaojiantao/pepper/SNP_Number_count

# Extract the polymrphic SNPs for each group
for j in annuum bac.baccatum bac.pendulum chacoense chinense frutescens glabriusculum pubescens ; do
    for i in $(seq -w 0 12); do
        vcftools --gzvcf ../vcf.clean/Chr$i.clean.vcf.gz \
            --keep ../Pop.list/Ca.$j.list \
            --recode --recode-INFO-all \
            --mac 1 \
            --out Ca.$j.Chr$i
    done
done

# zip the file

# SNP count for each group
for i in annuum bac.baccatum bac.pendulum chacoense chinense frutescens glabriusculum pubescens ; do
    zcat Ca.$i.Chr*.recode.vcf.gz | grep -v "#" | awk '{print $1"_"$2}' > Ca.$i.snps
done

# Further filter with maf = 0.05
for j in annuum bac.baccatum bac.pendulum chacoense chinense frutescens glabriusculum pubescens ; do
    for i in $(seq -w 0 12); do
        vcftools --gzvcf Ca.$j.Chr$i.recode.vcf.gz --maf 0.05 --recode --out Ca.$j.Chr$i.maf
    done
done

# SNP count for each group
for i in annuum bac.baccatum bac.pendulum chacoense chinense frutescens glabriusculum pubescens ; do
    cat Ca.$i.Chr*.maf.recode.vcf | grep -v "#" | awk '{print $1"_"$2}' > Ca.$i.maf.snps
done

for i in annuum bac.baccatum bac.pendulum chacoense chinense frutescens glabriusculum pubescens ; do
    zcat Ca.$i.Chr*.maf.recode.vcf.gz | grep -v "#" | awk '{print $1"_"$2}' > Ca.$i.maf.snps
done

```

### SNP prunning

```bash
# Working directory
cd /data/zhaojiantao/pepper/SNP.prune

# Convert vcf to ped/map
for i in $(seq -w 1 12); do
    plink --vcf ../vcf.clean/Chr$i.clean.vcf.gz \
    --allow-extra-chr --double-id \
    --recode --make-bed \
    --out Chr$i &
done

# Prepare the map file format
for i in $(seq -w 1 12); do
    awk '{print $1,$1"_"$4,$3,$4}' Chr$i.map | sed 's/ /\t/g' > Chr$i.map-1
    mv Chr$i.map-1 Chr$i.map
done

# SNP pruning
for i in $(seq -w 1 12); do
    plink --file Chr$i \
    --indep-pairwise 500 10 0.1 \
    --allow-extra-chr --double-id \
    --out Chr$i &
    # Extract the prunned SNPs
    plink --file Chr$i \
        --extract Chr$i.prune.in \
        --recode vcf --double-id --allow-extra-chr \
        --out Chr$i.prune &
done

# Total number of SNPs:298,000
change double ids to single id
sh single.id.sh
for i in $(seq -w 1 12); do
    grep -v "#" Chr$i.prune.vcf | cat single.id.info - > Chr$i.prune.vcf-1
    mv Chr$i.prune.vcf-1 Chr$i.prune.vcf
done

# Combine the vcf file
cat <(grep "#" Chr01.prune.vcf) \
    <(cat Chr*.prune.vcf | grep -v "#") \
    > Combined.prune.vcf
 
# Further filtering missing alleles
plink --file Combined.prune \
    --geno 0 --maf 0.05 \
    --recode 12 \
    --output-missing-genotype 9 \
    --out Combined.prune.nomiss

```

### SNP effects and distribution of different subgroups.

```bash
cd /data/zhaojiantao/pepper/SNP

# Generate the Plink map/ped format for each subgroup
for i in Ca.annuum Ca.bac.baccatum Ca.bac.pendulum Ca.chacoense Ca.chinense Ca.frutescens Ca.glabriusculum Ca.pubescens; do
    for j in $(seq -w 1 12); do
    vcftools --gzvcf ../vcf.clean/Chr$j.clean.vcf.gz \
        --mac 2 --max-missing 0.3 \
        --keep ../Pop.list/$i.list \
        --chrom-map <(echo Chr$j Chr$j | sed 's/ /\t/g') \
        --plink --out $i.Chr$j
    done 
done

for i in Ca.annuum Ca.bac.baccatum Ca.bac.pendulum Ca.chacoense Ca.chinense Ca.frutescens Ca.glabriusculum Ca.pubescens; do
    cat $i.Chr*.map | awk '{print $2}' > $i.combined.map &
done

paste *combined.map | sed '1i\Ca.annuum Ca.bac.baccatum Ca.bac.pendulum Ca.chacoense Ca.chinense Ca.frutescens Ca.glabriusculum Ca.pubescens' | sed 's/ /\t/g' > SNP.venn.txt

# Evaluate the SNP effects using [SnpEff](http://pcingola.github.io/SnpEff/se_running/)
# Build the [databases](https://www.cnblogs.com/zhanmaomao/p/10964636.html)

cd /data/zhaojiantao/tools/snpEff
# Prepare and bulid the genome index
cp /data/zhaojiantao/pepper/References/S8/Zhangshugang.fasta data/Zhangshugang/sequences.fa
cp /data/zhaojiantao/pepper/References/S8/Zhangshugang.fasta data/genomes/Zhangshugang.fa

cp /data/zhaojiantao/pepper/References/S8/Zhangshugang.gff data/Zhangshugang/genes.gff

echo "Zhangshugang.genome : Zhangshugang" >> snpEff.config

# Build the database
java -jar snpEff.jar build -gff3 -v Zhangshugang

# Calculate the SNP effects
for i in $(seq -w 0 12); do
    cd /data/zhaojiantao/pepper/vcf.clean/SnpEff.Chr$i
    java -Xmx8g -jar /data/zhaojiantao/tools/snpEff/snpEff.jar \
        -c /data/zhaojiantao/tools/snpEff/snpEff.config \
        Zhangshugang \
        ../Chr$i.clean.vcf.gz \
        > Chr$i.ann.vcf &
done

# Calculate indel effects
for i in $(seq -w 0 12); do
    cd /data/zhaojiantao/pepper/Indels/SnpEff.indel.Chr$i
    java -Xmx8g -jar /data/zhaojiantao/tools/snpEff/snpEff.jar \
        -c /data/zhaojiantao/tools/snpEff/snpEff.config \
        Zhangshugang \
        ../Chr$i.fltr5.indel.vcf.gz \
        > Chr$i.indel.ann.vcf &
done

# Calculate the uniq variants by taking the first effect
cd /data/zhaojiantao/pepper/vcf.clean
for i in $(seq -w 0 12); do
    zcat SnpEff.Chr$i/Chr$i.ann.vcf.gz | grep -v "#" | cut -f8 | sed 's/|/\t/g' | cut -f2 | sort | uniq -c > SnpEff.Chr$i/Chr$i.ann.count
done

# Calculate the uniq variants by taking the first effect
cd /data/zhaojiantao/pepper/all_indels_filter5
for i in $(seq -w 0 12); do
    grep -v "#" SnpEff.indel.Chr$i/Chr$i.indel.ann.vcf | cut -f8 | sed 's/;/\t/g' | awk '{for(i=1;i<=NF;i++){if ($i ~ /ANN/){print $i}}}' | \
    sed 's/|/\t/g' | cut -f2 | sort | uniq -c > SnpEff.indel.Chr$i/Chr$i.indel.ann.count
done

# Upset distribution of cleaned SNPs
for i in $(ll *combined.map.gz | awk '{print $9}' | sed 's/.combined.map.gz//g' ); do 
    zcat $i.combined.map.gz > $i.snp ; 
done

```

#### Upset plot in R

```r
library(UpSetR)
library(data.table)

Ca.annuum <- fread("Ca.annuum.snp", header=FALSE, sep = "\n")
Ca.bac.baccatum <- fread("Ca.bac.baccatum.snp", header=FALSE, sep = "\n")
Ca.bac.pendulum <- fread("Ca.bac.pendulum.snp", header=FALSE, sep = "\n")
Ca.chacoense <- fread("Ca.chacoense.snp", header=FALSE, sep = "\n")
Ca.chinense <- fread("Ca.chinense.snp", header=FALSE, sep = "\n")
Ca.frutescens <- fread("Ca.frutescens.snp", header=FALSE, sep = "\n")
Ca.glabriusculum <- fread("Ca.glabriusculum.snp", header=FALSE, sep = "\n")
Ca.pubescens <- fread("Ca.pubescens.snp", header=FALSE, sep = "\n")

Ca.annuum <- as.vector(Ca.annuum$V1)
Ca.bac.baccatum <- as.vector(Ca.bac.baccatum$V1)
Ca.bac.pendulum <- as.vector(Ca.bac.pendulum$V1)
Ca.chacoense <- as.vector(Ca.chacoense$V1)
Ca.chinense <- as.vector(Ca.chinense$V1)
Ca.frutescens <- as.vector(Ca.frutescens$V1)
Ca.glabriusculum <- as.vector(Ca.glabriusculum$V1)
Ca.pubescens <- as.vector(Ca.pubescens$V1)

read_sets = list(C.annuum = Ca.annuum, C.baccatum.var.baccatum = Ca.bac.baccatum,
    C.baccatum.var.pendulum = Ca.bac.pendulum, C.chacoense = Ca.chacoense,
    C.chinense = Ca.chinense, C.frutescens = Ca.frutescens,
    C.annuum.var.glabriusculum = Ca.glabriusculum, C.pubescens = Ca.pubescens)

jpeg("SNP_upset.jpeg", height = 15, width = 25, units="cm", res = 300)
pdf("SNP_upset.pdf", height = 6, width = 10)
upset(fromList(read_sets), order.by = "freq",number.angles = 90,
    mainbar.y.label = "SNP Intersections", sets.x.label = "Total cleaned SNPs",
    sets = c("C.annuum", "C.baccatum.var.baccatum", 
    "C.baccatum.var.pendulum","C.chacoense",
    "C.chinense","C.frutescens",
    "C.annuum.var.glabriusculum","C.pubescens"))
dev.off()

```



## 4DTV SNPs

```bash
# Working directory
cd /data/zhaojiantao/pepper/4DTV

# Call 4dtv SNPs
4DTV_scan.pl -g ../pepper_S8_annotation/fasta \
    -f ../pepper_S8_annotation/gff \
    -o ffds.list \
    > ffds.log

# Select the ids
cat ../4DTV/ffds.list | cut -f1,2 > ../4DTV/ffds.list.clean.txt

```
## IQ-tree

### IQ-Tree for the ~1300 accessions

```bash
# Working directory
cd /data/zhaojiantao/pepper/pepper_survey

# Select the Header
grep "#" All.snp.filter.vcf > aa

# Select those PASSED SNPs 
awk '$7=="PASS"' All.snp.filter.vcf > bb

# Combine with the header and rename the chromosomal ids consistent with reference
cat aa bb | \
    sed 's/Superscaffold1_54_332676983/Chr01/g' | \
    sed 's/Superscaffold2_27_177319215/Chr02/g' | \
    sed 's/Superscaffold3_40_289790774/Chr03/g' | \
    sed 's/Superscaffold4_26_248932513/Chr04/g' | \
    sed 's/Superscaffold5_64_254874144/Chr05/g' | \
    sed 's/Superscaffold6_82_253233553/Chr06/g' | \
    sed 's/Superscaffold7_54_266382521/Chr07/g' | \
    sed 's/Superscaffold8_30_174326481/Chr08/g' | \
    sed 's/Superscaffold9_43_278428952/Chr09/g' | \
    sed 's/Superscaffold10_31_210332287/Chr10/g' | \
    sed 's/Superscaffold11_37_275185330/Chr11/g' | \
    sed 's/Superscaffold12_60_260676116/Chr12/g' | \
    sed 's/unanchor1_1_17280/Contig01/g' | sed 's/unanchor2_1_14241/Contig02/g' | \
    sed 's/unanchor3_1_14184/Contig03/g' | sed 's/unanchor4_1_14135/Contig04/g' | \
    sed 's/unanchor5_1_13935/Contig05/g' | sed 's/unanchor6_1_13927/Contig06/g' | \
    sed 's/unanchor7_1_13733/Contig07/g' | sed 's/unanchor8_1_13524/Contig08/g' | \
    sed 's/unanchor9_1_13450/Contig09/g' | sed 's/unanchor10_1_13288/Contig10/g' | \
    sed 's/unanchor11_1_13186/Contig11/g' | sed 's/unanchor12_1_12934/Contig12/g' | \
    sed 's/unanchor13_1_12859/Contig13/g' | sed 's/unanchor14_1_12026/Contig14/g' | \
    sed 's/unanchor15_1_12004/Contig15/g' | sed 's/unanchor16_1_11909/Contig16/g' | \
    sed 's/unanchor17_1_11861/Contig17/g' | sed 's/unanchor18_1_11798/Contig18/g' | \
    sed 's/unanchor19_1_11724/Contig19/g' | sed 's/unanchor20_1_11540/Contig20/g' | \
    sed 's/unanchor21_1_10681/Contig21/g' | sed 's/unanchor22_1_10492/Contig22/g' | \
    sed 's/unanchor23_1_10491/Contig23/g' | sed 's/unanchor24_1_10144/Contig24/g' | \
    sed 's/unanchor25_1_10010/Contig25/g' | sed 's/unanchor26_1_9833/Contig26/g' | \
    sed 's/unanchor27_1_9659/Contig27/g' | sed 's/unanchor28_1_9285/Contig28/g' | \
    sed 's/unanchor29_1_9105/Contig29/g' | sed 's/unanchor30_1_9048/Contig30/g' | \
    sed 's/unanchor31_1_9046/Contig31/g' | sed 's/unanchor32_1_9019/Contig32/g' | \
    sed 's/unanchor33_1_8665/Contig33/g' | sed 's/unanchor34_1_8639/Contig34/g' | \
    sed 's/unanchor35_1_8180/Contig35/g' | sed 's/unanchor36_1_8163/Contig36/g' | \
    sed 's/unanchor37_1_8109/Contig37/g' | sed 's/unanchor38_1_8032/Contig38/g' | \
    sed 's/unanchor39_1_7734/Contig39/g' | sed 's/unanchor40_1_7545/Contig40/g' | \
    sed 's/unanchor41_1_7514/Contig41/g' | sed 's/unanchor42_1_7491/Contig42/g' | \
    sed 's/unanchor43_1_7293/Contig43/g' | sed 's/unanchor44_1_7178/Contig44/g' | \
    sed 's/unanchor45_1_6978/Contig45/g' | sed 's/unanchor46_1_6489/Contig46/g' | \
    sed 's/unanchor47_1_6375/Contig47/g' | sed 's/unanchor48_1_6084/Contig48/g' | \
    sed 's/unanchor49_1_5994/Contig49/g' | sed 's/unanchor50_1_5934/Contig50/g' | \
    sed 's/unanchor51_1_5922/Contig51/g' | sed 's/unanchor52_1_5576/Contig52/g' | \
    sed 's/unanchor53_1_5065/Contig53/g' | sed 's/unanchor54_1_4952/Contig54/g' | \
    sed 's/unanchor55_1_4829/Contig55/g' | sed 's/unanchor56_1_4778/Contig56/g' | \
    sed 's/unanchor57_1_4529/Contig57/g' | sed 's/unanchor58_1_4262/Contig58/g' | \
    sed 's/unanchor59_1_3903/Contig59/g' | sed 's/unanchor60_1_3887/Contig60/g' | \
    sed 's/unanchor61_1_3635/Contig61/g' | sed 's/unanchor62_1_3597/Contig62/g' | \
    sed 's/unanchor63_1_3447/Contig63/g' | sed 's/unanchor64_1_3080/Contig64/g' | \
    sed 's/unanchor65_1_2873/Contig65/g' | sed 's/unanchor66_1_2852/Contig66/g' | \
    sed 's/unanchor67_1_2512/Contig67/g' | sed 's/unanchor68_1_2453/Contig68/g' | \
    sed 's/unanchor69_1_2175/Contig69/g' | sed 's/unanchor70_1_2025/Contig70/g' | \
    sed 's/unanchor71_1_1781/Contig71/g' | sed 's/unanchor72_1_1707/Contig72/g' | \
    sed 's/unanchor73_1_1626/Contig73/g' | sed 's/unanchor74_1_1578/Contig74/g' | \
    sed 's/unanchor75_1_1279/Contig75/g' | sed 's/unanchor76_1_1221/Contig76/g' | \
    sed 's/unanchor77_1_936/Contig77/g' | sed 's/unanchor78_1_594/Contig78/g' | \
    sed 's/unanchor79_1_50000/Contig79/g' | sed 's/unanchor80_1_36440/Contig80/g' | \
    sed 's/unanchor81_1_188396/Contig81/g' | sed 's/unanchor82_1_233858/Contig82/g' | \
    sed 's/unanchor83_1_262094/Contig83/g' | sed 's/unanchor84_1_55865/Contig84/g' | \
    sed 's/unanchor85_1_116314/Contig85/g' | sed 's/unanchor86_1_125000/Contig86/g' | \
    sed 's/unanchor87_1_547308/Contig87/g' | sed 's/unanchor88_1_126556/Contig88/g' | \
    sed 's/unanchor89_1_444622/Contig89/g' | gzip \
    > All.snp.filtered.vcf.gz

# 4DTV snps ids
zcat ../IQ-Tree504/Ca.S8.ffds.vcf.gz | \
    grep -v "#" | cut -f1,2 | grep -v Ca.S8.Chr00 | \
    sed 's/Ca.S8.Chr01/Superscaffold1_54_332676983/g' | sed 's/Ca.S8.Chr02/Superscaffold2_27_177319215/g' | \
    sed 's/Ca.S8.Chr03/Superscaffold3_40_289790774/g' | sed 's/Ca.S8.Chr04/Superscaffold4_26_248932513/g' | \
    sed 's/Ca.S8.Chr05/Superscaffold5_64_254874144/g' | sed 's/Ca.S8.Chr06/Superscaffold6_82_253233553/g' | \
    sed 's/Ca.S8.Chr07/Superscaffold7_54_266382521/g' | sed 's/Ca.S8.Chr08/Superscaffold8_30_174326481/g' | \
    sed 's/Ca.S8.Chr09/Superscaffold9_43_278428952/g' | sed 's/Ca.S8.Chr10/Superscaffold10_31_210332287/g' | \
    sed 's/Ca.S8.Chr11/Superscaffold11_37_275185330/g' | sed 's/Ca.S8.Chr12/Superscaffold12_60_260676116/g' \
    > 504.4DTV.bisnps

# Extract the SNPs
vcftools --gzvcf 1314.vcf.gz \
    --positions 504.4DTV.bisnps \
    --recode --stdout > 1314.4DTV-1.vcf

# Select the SNPs
gatk SelectVariants \
    -R /data/zhaojiantao/pepper/References/S8/pepper_S8.fasta \
    -V 1314.4DTV-1.vcf \
    -select-type SNP \
    -O 1314.4DTV.biSNPs.vcf

# Filter the SNPs
echo Chr$i Chr$i | sed 's/ /\t/g' > Chr$i.clean-map.txt ; 
# Make ped file using this chromosome map
vcftools --vcf 1314.4DTV.biSNPs.vcf \
    --max-missing 0.5 --maf 0.05 \
    --recode --stdout | bgzip -c > 1314.4DTV.biSNPs.filter.vcf.gz

# Extract the tomato accessions as outgroup
bcftools view -s ERR418106,ERR418107,ERR418103,ERR418102 \
    /data/zhaojiantao/pepper/IQ-Tree504/Ca.S8.ffds.addTomato.vcf | \
    sed 's/Ca.S8.Chr01/Superscaffold1_54_332676983/g' | sed 's/Ca.S8.Chr02/Superscaffold2_27_177319215/g' | \
    sed 's/Ca.S8.Chr03/Superscaffold3_40_289790774/g' | sed 's/Ca.S8.Chr04/Superscaffold4_26_248932513/g' | \
    sed 's/Ca.S8.Chr05/Superscaffold5_64_254874144/g' | sed 's/Ca.S8.Chr06/Superscaffold6_82_253233553/g' | \
    sed 's/Ca.S8.Chr07/Superscaffold7_54_266382521/g' | sed 's/Ca.S8.Chr08/Superscaffold8_30_174326481/g' | \
    sed 's/Ca.S8.Chr09/Superscaffold9_43_278428952/g' | sed 's/Ca.S8.Chr10/Superscaffold10_31_210332287/g' | \
    sed 's/Ca.S8.Chr11/Superscaffold11_37_275185330/g' | sed 's/Ca.S8.Chr12/Superscaffold12_60_260676116/g' \
    > Tomato.vcf

# Find the common SNPs of tomato
zcat 1314.4DTV.biSNPs.filter.vcf | grep -v "#" | cut -f1,2 > 1314.4DTV.biSNPs.filter.snps
vcftools --vcf Tomato.vcf \
    --positions 1314.4DTV.biSNPs.filter.snps \
    --recode --stdout | bgzip -c > Tomato.4DTV.vcf.gz

# Extract biSNPs
tabix  Tomato.4DTV.vcf.gz
bcftools view --types snps -m 2 -M 2 \
    Tomato.4DTV.vcf.gz -O z \
    -o Tomato.4DTV.biSNPs.vcf.gz

zcat Tomato.4DTV.biSNPs.vcf.gz | grep -v "#" | cut -f1,2 > Tomato.4DTV.biSNPs
vcftools --gzvcf 1314.4DTV.biSNPs.filter.vcf.gz \
    --positions Tomato.4DTV.biSNPs \
    --recode --stdout | bgzip -c > 1314.4DTV.biSNPs.filter2.vcf.gz

# Combine the vcf files
bcftools index 1314.4DTV.biSNPs.filter2.vcf.gz
bcftools index Tomato.4DTV.biSNPs.vcf.gz

# Merge samples
bcftools merge -Oz \
    Tomato.4DTV.biSNPs.vcf.gz 1314.4DTV.biSNPs.filter2.vcf.gz \
    > 1314.4DTV.biSNPs.filter.withtomato.vcf.gz

# Convert vcf to phylip
python /data/zhaojiantao/tools/PythonScripts/vcf2phylip.py -i 1314.4DTV.biSNPs.filter.withtomato.vcf.gz
mv 1314.4DTV.biSNPs.filter.withtomato.min4.phy 1314.withtomato.phy

# Run IQ-tree
iqtree -s 1314.withtomato.phy -nt AUTO -m MFP

```




### IQ-tree for 504 accessions

```bash
# Working directory
cd /data/zhaojiantao/pepper/IQ-Tree504

# Convert vcf to phylip
python /data/zhaojiantao/tools/PythonScripts/vcf2phylip.py -i Ca.S8.ffds.addTomato.vcf.gz
mv ffds.addTomato.min4.phy ffds.addTomato.phy

# IQ-tree
iqtree -s ffds.addTomato.phy -nt AUTO -m MFP -o ERR418106

# SNP pruning for 4DTV SNPs
## Add the SNP infos in the vcf files
grep -v "#" Ca.S8.ffds.addTomato.vcf | \
    awk '{print $1,$2,$1"_"$2,$0}' | sed 's/Ca.S8.//g' | sed 's/ /\t/g' | \
    cut -f1-3,7- > Ca.S8.ffds.addTomato.vcf-1
grep "#" Ca.S8.ffds.addTomato.vcf | sed 's/Ca.S8.//g' | \
    cat - Ca.S8.ffds.addTomato.vcf-1 \
    > Ca.S8.ffds.addTomato.vcf-2
mv Ca.S8.ffds.addTomato.vcf-2 Ca.S8.ffds.addTomato.vcf

```

### IQ-tree for 504 accessions with prunned SNPs

```bash
# Working directory
cd /data/zhaojiantao/pepper/SNP.prune

# Convert vcf to phylip
python /data/zhaojiantao/tools/PythonScripts/vcf2phylip.py -i Combined.prune.vcf
# Shorten the sample ids
cat Combined.prune.min4.phy | \
    sed 's/ERR418102_ERR418102/ERR418102/g' | sed 's/ERR418103_ERR418103/ERR418103/g' | \
    sed 's/ERR418106_ERR418106/ERR418106/g' | sed 's/ERR418107_ERR418107/ERR418107/g' | \
    sed 's/CPPSIH_1_05_CPPSIH_1_05/CPPSIH_1_05/g' | sed 's/CXJ14_CXJ14/CXJ14/g' | \
    sed 's/CXJ34_CXJ34/CXJ34/g' | sed 's/CXJ38_CXJ38/CXJ38/g' | \
    sed 's/CXJ52_CXJ52/CXJ52/g' | sed 's/CXJ69_CXJ69/CXJ69/g' | \
    sed 's/CXJ92_CXJ92/CXJ92/g' | sed 's/CXJ96_CXJ96/CXJ96/g' | \
    sed 's/Grif_12454_Grif_12454/Grif_12454/g' | sed 's/Grif_1613_Grif_1613/Grif_1613/g' | \
    sed 's/Grif_1614_Grif_1614/Grif_1614/g' | sed 's/Grif_1615_Grif_1615/Grif_1615/g' | \
    sed 's/Grif_9118_Grif_9118/Grif_9118/g' | sed 's/Grif_9122_Grif_9122/Grif_9122/g' | \
    sed 's/Grif_9138_Grif_9138/Grif_9138/g' | sed 's/Grif_9140_Grif_9140/Grif_9140/g' | \
    sed 's/Grif_9142_Grif_9142/Grif_9142/g' | sed 's/Grif_9152_Grif_9152/Grif_9152/g' | \
    sed 's/Grif_9184_Grif_9184/Grif_9184/g' | sed 's/Grif_9188_Grif_9188/Grif_9188/g' | \
    sed 's/Grif_9194_Grif_9194/Grif_9194/g' | sed 's/Grif_9197_Grif_9197/Grif_9197/g' | \
    sed 's/Grif_9199_Grif_9199/Grif_9199/g' | sed 's/Grif_9201_Grif_9201/Grif_9201/g' | \
    sed 's/Grif_9226_Grif_9226/Grif_9226/g' | sed 's/Grif_9229_Grif_9229/Grif_9229/g' | \
    sed 's/Grif_9230_Grif_9230/Grif_9230/g' | sed 's/Grif_9236_Grif_9236/Grif_9236/g' | \
    sed 's/Grif_9241_Grif_9241/Grif_9241/g' | sed 's/Grif_9249_Grif_9249/Grif_9249/g' | \
    sed 's/Grif_9256_Grif_9256/Grif_9256/g' | sed 's/Grif_9257_Grif_9257/Grif_9257/g' | \
    sed 's/Grif_9267_Grif_9267/Grif_9267/g' | sed 's/Grif_9289_Grif_9289/Grif_9289/g' | \
    sed 's/Grif_9301_Grif_9301/Grif_9301/g' | sed 's/Grif_9312_Grif_9312/Grif_9312/g' | \
    sed 's/Grif_9313_Grif_9313/Grif_9313/g' | sed 's/Grif_9315_Grif_9315/Grif_9315/g' | \
    sed 's/Grif_9352_Grif_9352/Grif_9352/g' | sed 's/Grif_9354_Grif_9354/Grif_9354/g' | \
    sed 's/Grif_9355_Grif_9355/Grif_9355/g' | sed 's/Grif_9357_Grif_9357/Grif_9357/g' | \
    sed 's/Grif_9359_Grif_9359/Grif_9359/g' | sed 's/Grif_9365_Grif_9365/Grif_9365/g' | \
    sed 's/Grif_9376_Grif_9376/Grif_9376/g' | sed 's/Grif_9381_Grif_9381/Grif_9381/g' | \
    sed 's/Grif_9388_Grif_9388/Grif_9388/g' | sed 's/Grif_9391_Grif_9391/Grif_9391/g' | \
    sed 's/Grif_9393_Grif_9393/Grif_9393/g' | sed 's/Grif_9400_Grif_9400/Grif_9400/g' | \
    sed 's/Grif_9433_Grif_9433/Grif_9433/g' | sed 's/OLJ20_OLJ20/OLJ20/g' | \
    sed 's/PI_148629_PI_148629/PI_148629/g' | sed 's/PI_159235_PI_159235/PI_159235/g' | \
    sed 's/PI_159242_PI_159242/PI_159242/g' | sed 's/PI_159246_PI_159246/PI_159246/g' | \
    sed 's/PI_159249_PI_159249/PI_159249/g' | sed 's/PI_159255_PI_159255/PI_159255/g' | \
    sed 's/PI_159260_PI_159260/PI_159260/g' | sed 's/PI_159261_PI_159261/PI_159261/g' | \
    sed 's/PI_188481_PI_188481/PI_188481/g' | sed 's/PI_193470_PI_193470/PI_193470/g' | \
    sed 's/PI_194879_PI_194879/PI_194879/g' | sed 's/PI_194880_PI_194880/PI_194880/g' | \
    sed 's/PI_194910_PI_194910/PI_194910/g' | sed 's/PI_195296_PI_195296/PI_195296/g' | \
    sed 's/PI_195770_PI_195770/PI_195770/g' | sed 's/PI_199506_PI_199506/PI_199506/g' | \
    sed 's/PI_208738_PI_208738/PI_208738/g' | sed 's/PI_209590_PI_209590/PI_209590/g' | \
    sed 's/PI_215699_PI_215699/PI_215699/g' | sed 's/PI_215731_PI_215731/PI_215731/g' | \
    sed 's/PI_215741_PI_215741/PI_215741/g' | sed 's/PI_224405_PI_224405/PI_224405/g' | \
    sed 's/PI_224416_PI_224416/PI_224416/g' | sed 's/PI_224426_PI_224426/PI_224426/g' | \
    sed 's/PI_224427_PI_224427/PI_224427/g' | sed 's/PI_224428_PI_224428/PI_224428/g' | \
    sed 's/PI_238049_PI_238049/PI_238049/g' | sed 's/PI_238051_PI_238051/PI_238051/g' | \
    sed 's/PI_238053_PI_238053/PI_238053/g' | sed 's/PI_238054_PI_238054/PI_238054/g' | \
    sed 's/PI_238055_PI_238055/PI_238055/g' | sed 's/PI_238056_PI_238056/PI_238056/g' | \
    sed 's/PI_238060_PI_238060/PI_238060/g' | sed 's/PI_241650_PI_241650/PI_241650/g' | \
    sed 's/PI_241671_PI_241671/PI_241671/g' | sed 's/PI_241673_PI_241673/PI_241673/g' | \
    sed 's/PI_241675_PI_241675/PI_241675/g' | sed 's/PI_241676_PI_241676/PI_241676/g' | \
    sed 's/PI_241679_PI_241679/PI_241679/g' | sed 's/PI_257046_PI_257046/PI_257046/g' | \
    sed 's/PI_257050_PI_257050/PI_257050/g' | sed 's/PI_257062_PI_257062/PI_257062/g' | \
    sed 's/PI_257074_PI_257074/PI_257074/g' | sed 's/PI_257080_PI_257080/PI_257080/g' | \
    sed 's/PI_257083_PI_257083/PI_257083/g' | sed 's/PI_257104_PI_257104/PI_257104/g' | \
    sed 's/PI_257109_PI_257109/PI_257109/g' | sed 's/PI_257110_PI_257110/PI_257110/g' | \
    sed 's/PI_257116_PI_257116/PI_257116/g' | sed 's/PI_257117_PI_257117/PI_257117/g' | \
    sed 's/PI_257121_PI_257121/PI_257121/g' | sed 's/PI_257124_PI_257124/PI_257124/g' | \
    sed 's/PI_257141_PI_257141/PI_257141/g' | sed 's/PI_257148_PI_257148/PI_257148/g' | \
    sed 's/PI_257157_PI_257157/PI_257157/g' | sed 's/PI_257161_PI_257161/PI_257161/g' | \
    sed 's/PI_257173_PI_257173/PI_257173/g' | sed 's/PI_257174_PI_257174/PI_257174/g' | \
    sed 's/PI_257181_PI_257181/PI_257181/g' | sed 's/PI_257182_PI_257182/PI_257182/g' | \
    sed 's/PI_260427_PI_260427/PI_260427/g' | sed 's/PI_260429_PI_260429/PI_260429/g' | \
    sed 's/PI_260430_PI_260430/PI_260430/g' | sed 's/PI_260431_PI_260431/PI_260431/g' | \
    sed 's/PI_260435_PI_260435/PI_260435/g' | sed 's/PI_260436_PI_260436/PI_260436/g' | \
    sed 's/PI_260437_PI_260437/PI_260437/g' | sed 's/PI_260478_PI_260478/PI_260478/g' | \
    sed 's/PI_260479_PI_260479/PI_260479/g' | sed 's/PI_260485_PI_260485/PI_260485/g' | \
    sed 's/PI_260486_PI_260486/PI_260486/g' | sed 's/PI_260488_PI_260488/PI_260488/g' | \
    sed 's/PI_260492_PI_260492/PI_260492/g' | sed 's/PI_260498_PI_260498/PI_260498/g' | \
    sed 's/PI_260501_PI_260501/PI_260501/g' | sed 's/PI_260505_PI_260505/PI_260505/g' | \
    sed 's/PI_260508_PI_260508/PI_260508/g' | sed 's/PI_260509_PI_260509/PI_260509/g' | \
    sed 's/PI_260515_PI_260515/PI_260515/g' | sed 's/PI_260516_PI_260516/PI_260516/g' | \
    sed 's/PI_260535_PI_260535/PI_260535/g' | sed 's/PI_260538_PI_260538/PI_260538/g' | \
    sed 's/PI_260545_PI_260545/PI_260545/g' | sed 's/PI_260549_PI_260549/PI_260549/g' | \
    sed 's/PI_260550_PI_260550/PI_260550/g' | sed 's/PI_260551_PI_260551/PI_260551/g' | \
    sed 's/PI_260564_PI_260564/PI_260564/g' | sed 's/PI_260567_PI_260567/PI_260567/g' | \
    sed 's/PI_260575_PI_260575/PI_260575/g' | sed 's/PI_260582_PI_260582/PI_260582/g' | \
    sed 's/PI_260590_PI_260590/PI_260590/g' | sed 's/PI_260595_PI_260595/PI_260595/g' | \
    sed 's/PI_260611_PI_260611/PI_260611/g' | sed 's/PI_266041_PI_266041/PI_266041/g' | \
    sed 's/PI_267729_PI_267729/PI_267729/g' | sed 's/PI_267730_PI_267730/PI_267730/g' | \
    sed 's/PI_273415_PI_273415/PI_273415/g' | sed 's/PI_273419_PI_273419/PI_273419/g' | \
    sed 's/PI_281308_PI_281308/PI_281308/g' | sed 's/PI_281309_PI_281309/PI_281309/g' | \
    sed 's/PI_281310_PI_281310/PI_281310/g' | sed 's/PI_281314_PI_281314/PI_281314/g' | \
    sed 's/PI_281341_PI_281341/PI_281341/g' | sed 's/PI_281396_PI_281396/PI_281396/g' | \
    sed 's/PI_281398_PI_281398/PI_281398/g' | sed 's/PI_281407_PI_281407/PI_281407/g' | \
    sed 's/PI_281421_PI_281421/PI_281421/g' | sed 's/PI_281422_PI_281422/PI_281422/g' | \
    sed 's/PI_281423_PI_281423/PI_281423/g' | sed 's/PI_281424_PI_281424/PI_281424/g' | \
    sed 's/PI_281428_PI_281428/PI_281428/g' | sed 's/PI_281430_PI_281430/PI_281430/g' | \
    sed 's/PI_281436_PI_281436/PI_281436/g' | sed 's/PI_281440_PI_281440/PI_281440/g' | \
    sed 's/PI_281442_PI_281442/PI_281442/g' | sed 's/PI_281444_PI_281444/PI_281444/g' | \
    sed 's/PI_290980_PI_290980/PI_290980/g' | sed 's/PI_290982_PI_290982/PI_290982/g' | \
    sed 's/PI_290983_PI_290983/PI_290983/g' | sed 's/PI_291999_PI_291999/PI_291999/g' | \
    sed 's/PI_293349_PI_293349/PI_293349/g' | sed 's/PI_310435_PI_310435/PI_310435/g' | \
    sed 's/PI_311126_PI_311126/PI_311126/g' | sed 's/PI_315007_PI_315007/PI_315007/g' | \
    sed 's/PI_315020_PI_315020/PI_315020/g' | sed 's/PI_315021_PI_315021/PI_315021/g' | \
    sed 's/PI_315023_PI_315023/PI_315023/g' | sed 's/PI_321387_PI_321387/PI_321387/g' | \
    sed 's/PI_337524_PI_337524/PI_337524/g' | sed 's/PI_338992_PI_338992/PI_338992/g' | \
    sed 's/PI_355813_PI_355813/PI_355813/g' | sed 's/PI_355815_PI_355815/PI_355815/g' | \
    sed 's/PI_355817_PI_355817/PI_355817/g' | sed 's/PI_355819_PI_355819/PI_355819/g' | \
    sed 's/PI_355820_PI_355820/PI_355820/g' | sed 's/PI_357467_PI_357467/PI_357467/g' | \
    sed 's/PI_358811_PI_358811/PI_358811/g' | sed 's/PI_358968_PI_358968/PI_358968/g' | \
    sed 's/PI_360725_PI_360725/PI_360725/g' | sed 's/PI_360732_PI_360732/PI_360732/g' | \
    sed 's/PI_368073_PI_368073/PI_368073/g' | sed 's/PI_368080_PI_368080/PI_368080/g' | \
    sed 's/PI_368081_PI_368081/PI_368081/g' | sed 's/PI_368085_PI_368085/PI_368085/g' | \
    sed 's/PI_370008_PI_370008/PI_370008/g' | sed 's/PI_370009_PI_370009/PI_370009/g' | \
    sed 's/PI_370010_PI_370010/PI_370010/g' | sed 's/PI_379193_PI_379193/PI_379193/g' | \
    sed 's/PI_387833_PI_387833/PI_387833/g' | sed 's/PI_406847_PI_406847/PI_406847/g' | \
    sed 's/PI_406948_PI_406948/PI_406948/g' | sed 's/PI_406987_PI_406987/PI_406987/g' | \
    sed 's/PI_407450_PI_407450/PI_407450/g' | sed 's/PI_410407_PI_410407/PI_410407/g' | \
    sed 's/PI_414729_PI_414729/PI_414729/g' | sed 's/PI_415094_PI_415094/PI_415094/g' | \
    sed 's/PI_420377_PI_420377/PI_420377/g' | sed 's/PI_420378_PI_420378/PI_420378/g' | \
    sed 's/PI_420379_PI_420379/PI_420379/g' | sed 's/PI_424732_PI_424732/PI_424732/g' | \
    sed 's/PI_435917_PI_435917/PI_435917/g' | sed 's/PI_438535_PI_438535/PI_438535/g' | \
    sed 's/PI_438538_PI_438538/PI_438538/g' | sed 's/PI_438567_PI_438567/PI_438567/g' | \
    sed 's/PI_438619_PI_438619/PI_438619/g' | sed 's/PI_438624_PI_438624/PI_438624/g' | \
    sed 's/PI_438630_PI_438630/PI_438630/g' | sed 's/PI_438635_PI_438635/PI_438635/g' | \
    sed 's/PI_438637_PI_438637/PI_438637/g' | sed 's/PI_438638_PI_438638/PI_438638/g' | \
    sed 's/PI_438641_PI_438641/PI_438641/g' | sed 's/PI_438642_PI_438642/PI_438642/g' | \
    sed 's/PI_438643_PI_438643/PI_438643/g' | sed 's/PI_438667_PI_438667/PI_438667/g' | \
    sed 's/PI_439245_PI_439245/PI_439245/g' | sed 's/PI_439256_PI_439256/PI_439256/g' | \
    sed 's/PI_439311_PI_439311/PI_439311/g' | sed 's/PI_439359_PI_439359/PI_439359/g' | \
    sed 's/PI_439360_PI_439360/PI_439360/g' | sed 's/PI_439363_PI_439363/PI_439363/g' | \
    sed 's/PI_439366_PI_439366/PI_439366/g' | sed 's/PI_439367_PI_439367/PI_439367/g' | \
    sed 's/PI_439369_PI_439369/PI_439369/g' | sed 's/PI_439371_PI_439371/PI_439371/g' | \
    sed 's/PI_439378_PI_439378/PI_439378/g' | sed 's/PI_439380_PI_439380/PI_439380/g' | \
    sed 's/PI_439384_PI_439384/PI_439384/g' | sed 's/PI_439386_PI_439386/PI_439386/g' | \
    sed 's/PI_439390_PI_439390/PI_439390/g' | sed 's/PI_439395_PI_439395/PI_439395/g' | \
    sed 's/PI_439398_PI_439398/PI_439398/g' | sed 's/PI_439401_PI_439401/PI_439401/g' | \
    sed 's/PI_439402_PI_439402/PI_439402/g' | sed 's/PI_439403_PI_439403/PI_439403/g' | \
    sed 's/PI_439405_PI_439405/PI_439405/g' | sed 's/PI_439407_PI_439407/PI_439407/g' | \
    sed 's/PI_439410_PI_439410/PI_439410/g' | sed 's/PI_439411_PI_439411/PI_439411/g' | \
    sed 's/PI_439414_PI_439414/PI_439414/g' | sed 's/PI_439415_PI_439415/PI_439415/g' | \
    sed 's/PI_439421_PI_439421/PI_439421/g' | sed 's/PI_439426_PI_439426/PI_439426/g' | \
    sed 's/PI_439427_PI_439427/PI_439427/g' | sed 's/PI_439431_PI_439431/PI_439431/g' | \
    sed 's/PI_439439_PI_439439/PI_439439/g' | sed 's/PI_439450_PI_439450/PI_439450/g' | \
    sed 's/PI_439463_PI_439463/PI_439463/g' | sed 's/PI_439466_PI_439466/PI_439466/g' | \
    sed 's/PI_439472_PI_439472/PI_439472/g' | sed 's/PI_439475_PI_439475/PI_439475/g' | \
    sed 's/PI_439477_PI_439477/PI_439477/g' | sed 's/PI_439481_PI_439481/PI_439481/g' | \
    sed 's/PI_439488_PI_439488/PI_439488/g' | sed 's/PI_439492_PI_439492/PI_439492/g' | \
    sed 's/PI_439498_PI_439498/PI_439498/g' | sed 's/PI_439502_PI_439502/PI_439502/g' | \
    sed 's/PI_439504_PI_439504/PI_439504/g' | sed 's/PI_439509_PI_439509/PI_439509/g' | \
    sed 's/PI_439511_PI_439511/PI_439511/g' | sed 's/PI_439512_PI_439512/PI_439512/g' | \
    sed 's/PI_439521_PI_439521/PI_439521/g' | sed 's/PI_439526_PI_439526/PI_439526/g' | \
    sed 's/PI_439528_PI_439528/PI_439528/g' | sed 's/PI_441516_PI_441516/PI_441516/g' | \
    sed 's/PI_441523_PI_441523/PI_441523/g' | sed 's/PI_441531_PI_441531/PI_441531/g' | \
    sed 's/PI_441532_PI_441532/PI_441532/g' | sed 's/PI_441536_PI_441536/PI_441536/g' | \
    sed 's/PI_441537_PI_441537/PI_441537/g' | sed 's/PI_441541_PI_441541/PI_441541/g' | \
    sed 's/PI_441547_PI_441547/PI_441547/g' | sed 's/PI_441548_PI_441548/PI_441548/g' | \
    sed 's/PI_441558_PI_441558/PI_441558/g' | sed 's/PI_441578_PI_441578/PI_441578/g' | \
    sed 's/PI_441590_PI_441590/PI_441590/g' | sed 's/PI_441592_PI_441592/PI_441592/g' | \
    sed 's/PI_441594_PI_441594/PI_441594/g' | sed 's/PI_441597_PI_441597/PI_441597/g' | \
    sed 's/PI_441602_PI_441602/PI_441602/g' | sed 's/PI_441605_PI_441605/PI_441605/g' | \
    sed 's/PI_441608_PI_441608/PI_441608/g' | sed 's/PI_441615_PI_441615/PI_441615/g' | \
    sed 's/PI_441625_PI_441625/PI_441625/g' | sed 's/PI_441631_PI_441631/PI_441631/g' | \
    sed 's/PI_441656_PI_441656/PI_441656/g' | sed 's/PI_441657_PI_441657/PI_441657/g' | \
    sed 's/PI_441674_PI_441674/PI_441674/g' | sed 's/PI_441685_PI_441685/PI_441685/g' | \
    sed 's/PI_441690_PI_441690/PI_441690/g' | sed 's/PI_441699_PI_441699/PI_441699/g' | \
    sed 's/PI_441710_PI_441710/PI_441710/g' | sed 's/PI_446900_PI_446900/PI_446900/g' | \
    sed 's/PI_446909_PI_446909/PI_446909/g' | sed 's/PI_487450_PI_487450/PI_487450/g' | \
    sed 's/PI_487451_PI_487451/PI_487451/g' | sed 's/PI_487452_PI_487452/PI_487452/g' | \
    sed 's/PI_487457_PI_487457/PI_487457/g' | sed 's/PI_497972_PI_497972/PI_497972/g' | \
    sed 's/PI_497974_PI_497974/PI_497974/g' | sed 's/PI_497984_PI_497984/PI_497984/g' | \
    sed 's/PI_497985_PI_497985/PI_497985/g' | sed 's/PI_511885_PI_511885/PI_511885/g' | \
    sed 's/PI_532990_PI_532990/PI_532990/g' | sed 's/PI_543181_PI_543181/PI_543181/g' | \
    sed 's/PI_543182_PI_543182/PI_543182/g' | sed 's/PI_543187_PI_543187/PI_543187/g' | \
    sed 's/PI_543195_PI_543195/PI_543195/g' | sed 's/PI_543197_PI_543197/PI_543197/g' | \
    sed 's/PI_543198_PI_543198/PI_543198/g' | sed 's/PI_543199_PI_543199/PI_543199/g' | \
    sed 's/PI_543204_PI_543204/PI_543204/g' | sed 's/PI_543207_PI_543207/PI_543207/g' | \
    sed 's/PI_555611_PI_555611/PI_555611/g' | sed 's/PI_555612_PI_555612/PI_555612/g' | \
    sed 's/PI_555616_PI_555616/PI_555616/g' | sed 's/PI_555620_PI_555620/PI_555620/g' | \
    sed 's/PI_555638_PI_555638/PI_555638/g' | sed 's/PI_555642_PI_555642/PI_555642/g' | \
    sed 's/PI_555645_PI_555645/PI_555645/g' | sed 's/PI_560944_PI_560944/PI_560944/g' | \
    sed 's/PI_573337_PI_573337/PI_573337/g' | sed 's/PI_585247_PI_585247/PI_585247/g' | \
    sed 's/PI_585250_PI_585250/PI_585250/g' | sed 's/PI_585254_PI_585254/PI_585254/g' | \
    sed 's/PI_585264_PI_585264/PI_585264/g' | sed 's/PI_585265_PI_585265/PI_585265/g' | \
    sed 's/PI_585266_PI_585266/PI_585266/g' | sed 's/PI_585267_PI_585267/PI_585267/g' | \
    sed 's/PI_585268_PI_585268/PI_585268/g' | sed 's/PI_585269_PI_585269/PI_585269/g' | \
    sed 's/PI_585272_PI_585272/PI_585272/g' | sed 's/PI_585273_PI_585273/PI_585273/g' | \
    sed 's/PI_585274_PI_585274/PI_585274/g' | sed 's/PI_585275_PI_585275/PI_585275/g' | \
    sed 's/PI_585276_PI_585276/PI_585276/g' | sed 's/PI_592528_PI_592528/PI_592528/g' | \
    sed 's/PI_593490_PI_593490/PI_593490/g' | sed 's/PI_593500_PI_593500/PI_593500/g' | \
    sed 's/PI_593507_PI_593507/PI_593507/g' | sed 's/PI_593517_PI_593517/PI_593517/g' | \
    sed 's/PI_593519_PI_593519/PI_593519/g' | sed 's/PI_593526_PI_593526/PI_593526/g' | \
    sed 's/PI_593527_PI_593527/PI_593527/g' | sed 's/PI_593545_PI_593545/PI_593545/g' | \
    sed 's/PI_593546_PI_593546/PI_593546/g' | sed 's/PI_593547_PI_593547/PI_593547/g' | \
    sed 's/PI_593555_PI_593555/PI_593555/g' | sed 's/PI_593557_PI_593557/PI_593557/g' | \
    sed 's/PI_593574_PI_593574/PI_593574/g' | sed 's/PI_593576_PI_593576/PI_593576/g' | \
    sed 's/PI_593577_PI_593577/PI_593577/g' | sed 's/PI_593606_PI_593606/PI_593606/g' | \
    sed 's/PI_593607_PI_593607/PI_593607/g' | sed 's/PI_593608_PI_593608/PI_593608/g' | \
    sed 's/PI_593611_PI_593611/PI_593611/g' | sed 's/PI_593613_PI_593613/PI_593613/g' | \
    sed 's/PI_593616_PI_593616/PI_593616/g' | sed 's/PI_593617_PI_593617/PI_593617/g' | \
    sed 's/PI_593618_PI_593618/PI_593618/g' | sed 's/PI_593619_PI_593619/PI_593619/g' | \
    sed 's/PI_593620_PI_593620/PI_593620/g' | sed 's/PI_593621_PI_593621/PI_593621/g' | \
    sed 's/PI_593623_PI_593623/PI_593623/g' | sed 's/PI_593624_PI_593624/PI_593624/g' | \
    sed 's/PI_593625_PI_593625/PI_593625/g' | sed 's/PI_593626_PI_593626/PI_593626/g' | \
    sed 's/PI_593627_PI_593627/PI_593627/g' | sed 's/PI_593628_PI_593628/PI_593628/g' | \
    sed 's/PI_593630_PI_593630/PI_593630/g' | sed 's/PI_593632_PI_593632/PI_593632/g' | \
    sed 's/PI_593633_PI_593633/PI_593633/g' | sed 's/PI_593634_PI_593634/PI_593634/g' | \
    sed 's/PI_593636_PI_593636/PI_593636/g' | sed 's/PI_593637_PI_593637/PI_593637/g' | \
    sed 's/PI_593638_PI_593638/PI_593638/g' | sed 's/PI_593639_PI_593639/PI_593639/g' | \
    sed 's/PI_593640_PI_593640/PI_593640/g' | sed 's/PI_593641_PI_593641/PI_593641/g' | \
    sed 's/PI_593642_PI_593642/PI_593642/g' | sed 's/PI_593921_PI_593921/PI_593921/g' | \
    sed 's/PI_593922_PI_593922/PI_593922/g' | sed 's/PI_593924_PI_593924/PI_593924/g' | \
    sed 's/PI_593926_PI_593926/PI_593926/g' | sed 's/PI_594137_PI_594137/PI_594137/g' | \
    sed 's/PI_594138_PI_594138/PI_594138/g' | sed 's/PI_594139_PI_594139/PI_594139/g' | \
    sed 's/PI_595905_PI_595905/PI_595905/g' | sed 's/PI_595907_PI_595907/PI_595907/g' | \
    sed 's/PI_596059_PI_596059/PI_596059/g' | sed 's/PI_631129_PI_631129/PI_631129/g' | \
    sed 's/PI_631135_PI_631135/PI_631135/g' | sed 's/PI_631136_PI_631136/PI_631136/g' | \
    sed 's/PI_631137_PI_631137/PI_631137/g' | sed 's/PI_631144_PI_631144/PI_631144/g' | \
    sed 's/PI_631146_PI_631146/PI_631146/g' | sed 's/PI_631151_PI_631151/PI_631151/g' | \
    sed 's/PI_631152_PI_631152/PI_631152/g' | sed 's/PI_632918_PI_632918/PI_632918/g' | \
    sed 's/PI_632923_PI_632923/PI_632923/g' | sed 's/PI_632926_PI_632926/PI_632926/g' | \
    sed 's/PI_632928_PI_632928/PI_632928/g' | sed 's/PI_632930_PI_632930/PI_632930/g' | \
    sed 's/PI_633751_PI_633751/PI_633751/g' | sed 's/PI_633755_PI_633755/PI_633755/g' | \
    sed 's/PI_633757_PI_633757/PI_633757/g' | sed 's/PI_633758_PI_633758/PI_633758/g' | \
    sed 's/PI_635818_PI_635818/PI_635818/g' | sed 's/PI_639126_PI_639126/PI_639126/g' | \
    sed 's/PI_639129_PI_639129/PI_639129/g' | sed 's/PI_639132_PI_639132/PI_639132/g' | \
    sed 's/PI_639140_PI_639140/PI_639140/g' | sed 's/PI_639639_PI_639639/PI_639639/g' | \
    sed 's/PI_639640_PI_639640/PI_639640/g' | sed 's/PI_639646_PI_639646/PI_639646/g' | \
    sed 's/PI_639648_PI_639648/PI_639648/g' | sed 's/PI_639651_PI_639651/PI_639651/g' | \
    sed 's/PI_639652_PI_639652/PI_639652/g' | sed 's/PI_639653_PI_639653/PI_639653/g' | \
    sed 's/PI_639657_PI_639657/PI_639657/g' | sed 's/PI_639659_PI_639659/PI_639659/g' | \
    sed 's/PI_639661_PI_639661/PI_639661/g' | sed 's/PI_639662_PI_639662/PI_639662/g' | \
    sed 's/PI_639682_PI_639682/PI_639682/g' | sed 's/PI_640503_PI_640503/PI_640503/g' | \
    sed 's/PI_640770_PI_640770/PI_640770/g' | sed 's/PI_640882_PI_640882/PI_640882/g' | \
    sed 's/PI_640884_PI_640884/PI_640884/g' | sed 's/PI_640889_PI_640889/PI_640889/g' | \
    sed 's/PI_640897_PI_640897/PI_640897/g' | sed 's/PI_640900_PI_640900/PI_640900/g' | \
    sed 's/PI_640901_PI_640901/PI_640901/g' | sed 's/PI_640902_PI_640902/PI_640902/g' | \
    sed 's/PI_640909_PI_640909/PI_640909/g' | sed 's/PI_642956_PI_642956/PI_642956/g' | \
    sed 's/PI_643124_PI_643124/PI_643124/g' | sed 's/PI_645493_PI_645493/PI_645493/g' | \
    sed 's/PI_645508_PI_645508/PI_645508/g' | sed 's/PI_645534_PI_645534/PI_645534/g' | \
    sed 's/PI_645555_PI_645555/PI_645555/g' | sed 's/PI_645556_PI_645556/PI_645556/g' | \
    sed 's/PI_645557_PI_645557/PI_645557/g' | sed 's/PI_645560_PI_645560/PI_645560/g' | \
    sed 's/PI_645681_PI_645681/PI_645681/g' | sed 's/PI_653659_PI_653659/PI_653659/g' | \
    sed 's/PI_653671_PI_653671/PI_653671/g' | sed 's/PI_653677_PI_653677/PI_653677/g' | \
    sed 's/PI_653679_PI_653679/PI_653679/g' | sed 's/PI_653681_PI_653681/PI_653681/g' | \
    sed 's/PI_653744_PI_653744/PI_653744/g' | sed 's/PI_653747_PI_653747/PI_653747/g' | \
    sed 's/PI_655064_PI_655064/PI_655064/g' | sed 's/PI_659105_PI_659105/PI_659105/g' | \
    sed 's/PI_659106_PI_659106/PI_659106/g' | sed 's/PI_660971_PI_660971/PI_660971/g' | \
    sed 's/PI_660972_PI_660972/PI_660972/g' | sed 's/PI_661082_PI_661082/PI_661082/g' | \
    sed 's/PI_661088_PI_661088/PI_661088/g' | sed 's/PI_665008_PI_665008/PI_665008/g' | \
    sed 's/PI_666511_PI_666511/PI_666511/g' | sed 's/PI_666550_PI_666550/PI_666550/g' | \
    sed 's/PI_666556_PI_666556/PI_666556/g' | sed 's/PI_666563_PI_666563/PI_666563/g' | \
    sed 's/PI_666566_PI_666566/PI_666566/g' | sed 's/PI_666567_PI_666567/PI_666567/g' | \
    sed 's/PI_666580_PI_666580/PI_666580/g' | sed 's/PI_666586_PI_666586/PI_666586/g' | \
    sed 's/PI_666589_PI_666589/PI_666589/g' | sed 's/PI_666591_PI_666591/PI_666591/g' | \
    sed 's/PI_666594_PI_666594/PI_666594/g' | sed 's/PI_666595_PI_666595/PI_666595/g' | \
    sed 's/PI_666596_PI_666596/PI_666596/g' | sed 's/PI_674456_PI_674456/PI_674456/g' | \
    sed 's/PI_674459_PI_674459/PI_674459/g' | sed 's/PI_675162_PI_675162/PI_675162/g' | \
    sed 's/PI_678386_PI_678386/PI_678386/g' | sed 's/PI_688675_PI_688675/PI_688675/g' | \
    sed 's/YMH38_YMH38/YMH38/g' | sed 's/Z111_Z111/Z111/g' | \
    sed 's/bailajiao_bailajiao/bailajiao/g' | sed 's/changanyingjiao_changanyingjiao/changanyingjiao/g' | \
    sed 's/chaotiaojiao_chaotiaojiao/chaotiaojiao/g' | sed 's/qixingjiao_qixingjiao/qixingjiao/g' | \
    sed 's/sanweijiao_sanweijiao/sanweijiao/g' | sed 's/sanyingjiao_sanyingjiao/sanyingjiao/g' | \
    sed 's/xiaoguojiao_xiaoguojiao/xiaoguojiao/g' | sed 's/xiaomila_xiaomila/xiaomila/g' | \
    sed 's/zhudachang_zhudachang/zhudachang/g' | sed 's/zoukejiao_zoukejiao/zoukejiao/g' \
    > Combined.prune.min4.phy2
mv Combined.prune.min4.phy2 Combined.prune.min4.phy

# IQ-tree
iqtree -s Combined.prune.min4.phy -nt 80 -m MFP

```


## PCA analysis

```bash
# Working directory
cd /data/zhaojiantao/pepper/PCA

# Create chromosome map
bcftools view -H ../4DTV/ffds.vcf.gz | cut -f 1 | uniq | awk '{print $0"\t"$0}' > ffds.chrom-map.txt

# Make ped file using this chromosome map
vcftools --gzvcf ../4DTV/ffds.vcf.gz --plink --chrom-map ffds.chrom-map.txt --out ffds

# Make bed for the ped/map plink file
plink --file ffds --make-bed --allow-extra-chr --out ffds

# Calculate PCA
plink --file ffds --allow-extra-chr --pca --out ffds

```

## Population structure analysis

```bash
# Working directory
cd /data/zhaojiantao/pepper/fastStructure

# Install fastStructure
# Be careful: fastStructure is based on python2, rathor than python3
# Copy fastStructure file
git clone https://github.com/rajanil/fastStructure
cd fastStructure
git fetch
git merge origin/master
# Build library extensions, with force extention -f, which is crucial to run
cd ~/proj/fastStructure/vars
python setup.py build_ext -f --inplace
cd ~/proj/fastStructure
python setup.py build_ext -f --inplace

# Run fastStructure
# Create chromosome map
bcftools view -H ../IQ-Tree504/ffds.thin.vcf.gz | cut -f 1 | uniq | awk '{print $0"\t"$0}' > ffds.thin.chrom-map.txt

# Make ped file using this chromosome map
vcftools --gzvcf ../IQ-Tree504/ffds.thin.vcf.gz --plink --chrom-map ffds.thin.chrom-map.txt --out ffds.thin

# Make bed for the ped/map plink file
plink --file ffds.thin --make-bed --allow-extra-chr --out ffds.thin

# Remove unnecessary files
rm ffds.thin.chrom-map.txt ffds.thin.log ffds.thin.nosex 

# fastStructure
for i in $(seq 1 20); do 
    python ~/tools/fastStructure/structure.py -K $i --input=ffds.thin --output=fastStructure/output --full --seed=100
done

# Choose the optimal K (K = 7)
python ~/tools/fastStructure/chooseK.py --input=fastStructure/output

# Combine the fastStructure results at K = 7
# Run the fastStructure for 30 times using all ffds SNPs
for i in $(seq 1 30); do 
    python ~/tools/fastStructure/structure.py -K 7 --input=ffds --output=fastStructure/K7/Run$i --full --seed=100
done

# Remove unnecessary files Run*.7.meanQ.ind 
rm id empty Run*.7.meanP Run*.7.meanQ.ind Run*.7.varP Run*.7.varQ Run*.7.log K7.miscfile K7.perm_datafile 

```

### Population structure use pruned 4DTV SNPs

```bash
# Working directory
cd /data/zhaojiantao/pepper/SNP.prune

# fastStructure use pruned 4DTV
for i in $(seq 1 10); do 
    python2 /data/zhaojiantao/tools/fastStructure/structure.py -K $i \
        --input=Combined.prune \
        --output=output \
        --full --seed=100
done

```

## Genome-wide LD

```bash
# Check the subgroup of Ca.glabriusculum
cd /data/zhaojiantao/pepper/LD/Phylogeny
# Extract the Ca.glabriusculum and one Ca.pubescens to test
vcftools --gzvcf ../../IQ-Tree504/ffds.thin.vcf.gz \
    --keep ../../Pop.list/Ca.glabriusculum.check.list \
    --recode --recode-INFO-all \
    --out Ca.glabriusculum.check

# Convert vcf to phylip
python /data/zhaojiantao/tools/PythonScripts/vcf2phylip.py -i Ca.glabriusculum.check.recode.vcf

# IQ-tree
iqtree -s Ca.glabriusculum.check.recode.min4.phy -nt AUTO -m MFP

# Quality control of genotypes: maf >= 0.05
for i in $(seq -w 0 12); do
    vcftools --gzvcf ../vcf.raw/Chr$i.vcf.gz --maf 0.05 --recode --stdout | bgzip -c > ../vcf.clean/Chr$i.clean.vcf.gz
done

# Use PopLDdecay
for i in $(seq -w 1 12); do
    for j in Ca.annuum Ca.bac.baccatum Ca.bac.pendulum Ca.chacoense Ca.chinense Ca.frutescens Ca.glabriusculum Ca.glabriusculum_II Ca.pubescens ; do
        PopLDdecay -InVCF ../vcf.clean/Chr$i.clean.vcf.gz -SubPop ../Pop.list/$j.list \
            -MaxDist 2000 -MAF 0.05 -OutStat $j.Chr$i ;
    done
done

# LD for the candidate region of chromosome 6
# Extract the region
plink --vcf ../vcf.clean/Chr06.clean.vcf.gz  \
    --allow-extra-chr --chr Chr06 --geno 0.2 --maf 0.05 \
    --from-bp 243400000 --to-bp 244400000 \
    --recode vcf --threads 40 --double-id \
    --out Chr06.region

for i in Ca.annuum Ca.bac.baccatum Ca.bac.pendulum Ca.chinense Ca.frutescens Ca.pubescens ; do
    PopLDdecay -InVCF Chr06.region.vcf -SubPop <(awk '{print $1"_"$1}' ../Pop.list/$i.list) \
        -MaxDist 1000 -MAF 0.05 -OutStat Chr06.region.$i ;
    /data/zhaojiantao/tools/PopLDdecay-master/bin/Plot_OnePop.pl -inFile Chr06.region.$i.stat.gz -output Chr06.region.$i  
done

less Chr06.region.Ca.bac.baccatum.bin.gz | awk '$1<400000 {print $1,$2,"A01"}' > Chr06.BAB.region.ld
less Chr06.region.Ca.bac.pendulum.bin.gz | awk '$1<400000 {print $1,$2,"A02"}' > Chr06.BAP.region.ld
less Chr06.region.Ca.frutescens.bin.gz | awk '$1<400000 {print $1,$2,"A03"}' > Chr06.FRU.region.ld
less Chr06.region.Ca.chinense.bin.gz | awk '$1<400000 {print $1,$2,"A04"}' > Chr06.CHN.region.ld
less Chr06.region.Ca.annuum.bin.gz | awk '$1<400000 {print $1,$2,"A05"}' > Chr06.ANN.region.ld
cat *ld 

# Plot the LD
for i in Ca.annuum Ca.bac.baccatum Ca.bac.pendulum Ca.chacoense Ca.chinense Ca.frutescens Ca.glabriusculum Ca.glabriusculum_II Ca.pubescens ; do
    for j in $(seq -w 1 12); do
        /data/zhaojiantao/tools/PopLDdecay-master/bin/Plot_OnePop.pl -inFile $i.Chr$j.stat.gz -output $i.Chr$j
    done
done

# Add chromosome information to each stat files
for i in Ca.annuum Ca.bac.baccatum Ca.bac.pendulum Ca.chacoense Ca.chinense Ca.frutescens Ca.glabriusculum Ca.glabriusculum_II Ca.pubescens ; do
    zcat $i.Chr01.bin.gz | awk '{print $1,$2,1}' > $i.Chr01.bin.gz-2 ;
    zcat $i.Chr02.bin.gz | awk '{print $1,$2,2}' > $i.Chr02.bin.gz-2 ;
    zcat $i.Chr03.bin.gz | awk '{print $1,$2,3}' > $i.Chr03.bin.gz-2 ;
    zcat $i.Chr04.bin.gz | awk '{print $1,$2,4}' > $i.Chr04.bin.gz-2 ;
    zcat $i.Chr05.bin.gz | awk '{print $1,$2,5}' > $i.Chr05.bin.gz-2 ;
    zcat $i.Chr06.bin.gz | awk '{print $1,$2,6}' > $i.Chr06.bin.gz-2 ;
    zcat $i.Chr07.bin.gz | awk '{print $1,$2,7}' > $i.Chr07.bin.gz-2 ;
    zcat $i.Chr08.bin.gz | awk '{print $1,$2,8}' > $i.Chr08.bin.gz-2 ;
    zcat $i.Chr09.bin.gz | awk '{print $1,$2,9}' > $i.Chr09.bin.gz-2 ;
    zcat $i.Chr10.bin.gz | awk '{print $1,$2,10}' > $i.Chr10.bin.gz-2 ;
    zcat $i.Chr11.bin.gz | awk '{print $1,$2,11}' > $i.Chr11.bin.gz-2 ;
    zcat $i.Chr12.bin.gz | awk '{print $1,$2,12}' > $i.Chr12.bin.gz-2 ;
done

# Combine the results for each subgroup
for i in Ca.annuum Ca.bac.baccatum Ca.bac.pendulum Ca.chacoense Ca.chinense Ca.frutescens Ca.glabriusculum Ca.pubescens ; do
cat $i.Chr*.bin.gz-2 | grep -v "#" > $i.ld
done

# Combine all subgroups together
paste Ca.annuum.ld Ca.bac.baccatum.ld Ca.bac.pendulum.ld Ca.chacoense.ld \
    Ca.chinense.ld Ca.frutescens.ld Ca.glabriusculum.ld Ca.glabriusculum_II.ld Ca.pubescens.ld | \
    awk '{print $1,$3,$2,$5,$8,$11,$14,$17,$20,$23}' | \
    sed '1i\Dist Chr Ca.annuum Ca.bac.baccatum Ca.bac.pendulum Ca.chacoense Ca.chinense Ca.frutescens Ca.glabriusculum Ca.glabriusculum_II Ca.pubescens' | sed 's/ /\t/g' > Group.ld.txt

# Calculate the average LD across 12 chromosomes
for i in Ca.annuum Ca.bac.baccatum Ca.bac.pendulum Ca.chacoense Ca.chinense Ca.frutescens \
    Ca.glabriusculum Ca.glabriusculum_II Ca.pubescens ; do
    paste $i.Chr*.bin.gz-2 | sed '1d' | \
        awk '{print $1,($2+$5+$8+$11+$14+$17+$20+$23+$26+$29+$32+$35)/12}' \
        > $i.average.ld
done

# Combine the results
paste *average.ld | awk '{print $1,$2,$4,$6,$8,$10,$12,$14,$16,$18}' | \
    sed '1i\Dist Ca.annuum Ca.bac.baccatum Ca.bac.pendulum Ca.chacoense Ca.chinense Ca.frutescens Ca.glabriusculum Ca.glabriusculum_II Ca.pubescens' | sed 's/ /\t/g' > Group.average.ld.txt


# Remove unnecessary files
rm *png *stat.gz *pdf *bin.gz *bin.gz-2 *average.ld

# Find the ld at LD =0.3
for i in Ca.annuum Ca.bac.baccatum Ca.bac.pendulum Ca.chacoense Ca.chinense Ca.frutescens \
    Ca.glabriusculum Ca.glabriusculum_II Ca.pubescens; do 
    echo $i ;
    cat $i.average.ld | awk '$2>0.2999 && $2<0.3001 {print $0}' 
done

# Find the minimum ld
for i in Ca.annuum Ca.bac.baccatum Ca.bac.pendulum Ca.chacoense Ca.chinense Ca.frutescens Ca.glabriusculum Ca.pubescens; do 
    cat $i.ld | awk '$3==1 {print $0}' | sort -nk2 | head -1 | cut -f2 -d ' ' ; 
done

```

## Fst differentiation

```bash
# Working directory
cd /data/zhaojiantao/pepper/Fst/

# Calculate Fst between subgroups
for i in annuum bac.baccatum bac.pendulum chacoense chinense frutescens glabriusculum glabriusculum_II pubescens; do
    for j in bac.baccatum bac.pendulum chacoense chinense frutescens glabriusculum glabriusculum_II pubescens; do
        for chr in $(seq -w 1 12); do
        vcftools --gzvcf ../vcf.clean/Chr$chr.clean.vcf.gz \
            --weir-fst-pop ../Pop.list/Ca.$i.list \
            --weir-fst-pop ../Pop.list/Ca.$j.list \
            --fst-window-size 100000 \
            --fst-window-step 10000 \
            --out Ca.$i.$j.Chr$chr 
        done &
    done
done

for chr in $(seq -w 1 12); do
    vcftools --gzvcf ../vcf.clean/Chr$chr.clean.vcf.gz \
        --weir-fst-pop ../Pop.list/Ca.bac.pendulum.list \
        --weir-fst-pop ../Pop.list/Ca.glabriusculum_I.list \
        --fst-window-size 100000 \
        --fst-window-step 10000 \
        --out Ca.bac.pendulum.glabriusculum_I.Chr$chr &
done


# Combine the results
for i in $(ll *weir.fst | awk '{print $9}' | sed 's/.Chr/\t/g' | cut -f1 | sort | uniq); do 
    cat $i.Chr*.windowed.weir.fst | \
        grep -v CHROM | \
        awk '{print $7=$1"_"$2"_"$3,$1,$2,$3,$8=($2+$3)/2,$4,$5,$6}' | \
        awk '$8 > 0 {print $0}' | \
        awk '{gsub(/Chr0/, "", $2)} 1' | \
        awk '{gsub(/Chr/, "", $2)} 1' | \
        sed '1i\ID CHROM BIN_START BIN_END BIN_Middle N_VARIANTS WEIGHTED_FST MEAN_FST' | \
        sed 's/ /\t/g' \
        > $i.fst.txt
done

# Remove chromosomal level files
rm *.windowed.weir.fst *log 

# Calculate the mean Fst value
for i in $(ll *fst.txt.gz | awk '{print $9}' | sed 's/.fst.txt//g'); do
    echo $i; zcat $i.fst.txt.gz | grep -v MEAN | awk '{sum+=$8} END {print sum/NR}'
done

# Remove unnecessary files
rm *fst.txt-2 common.id

# Check the total line numbers to determine the top 1% lines
zcat Ca.frutescens.annuum.fst.txt.gz | wc -l
zcat Ca.frutescens.annuum.fst.txt.gz | sort -rgk8 | head -2970 | cut -f1-4,8 | sed '1i\SNP\tChr\tStart\tEnd\tFst' \
    > Ca.frutescens.annuum.top1.fst
# Numberic sort for each chromosome
for i in $(seq -w 1 12); do
    grep Chr$i Ca.frutescens.annuum.top1.fst | sort -nk3 > Ca.frutescens.annuum.top1.fst.Chr$i
done
# Combine the sort
cat Ca.frutescens.annuum.top1.fst.Chr* > Ca.frutescens.annuum.top1.fst.sorted
# Merge the top 1% windows
perl merge.windows.pl Ca.frutescens.annuum.top1.fst.sorted | sed '1d' > Ca.frutescens.annuum.top1.fst.merged.txt
# Remove unnecessary files
rm Ca.frutescens.annuum.top1.fst Ca.frutescens.annuum.top1.fst.Chr*

# Find the number of genes 
for i in $(ll *merged.txt | awk '{print $9}' | sed 's/.txt//g') ; do
    bedtools intersect \
        -a <(sed 's/Chr0//g' ../References/S8/Ca.S8.gff | sed 's/Chr//g') \
        -b <(cut -f1-3 $i.txt) | \
        grep mRNA | sed 's/;/\t/g' | awk '{print $9}' | sed 's/ID=//g' | grep rna \
        > $i.genes.txt
done

# Uniq gene count
for i in $(ll *merged.genes.txt | awk '{print $9}' | sed 's/.txt//g') ; do
    sed 's/\./\t/g' $i.txt | cut -f1 | sort | uniq > $i.uniq
done

# Calculate the total size of differentiated regions
for i in *top1.fst.merged.txt ; do
    awk '{s+=$3-$2+1}END{print s}' $i 
done

# Prepare the file for significant test
for i in $(ls *fst.txt.gz | sed 's/.fst.txt.gz//g' | sed 's/Ca.//g'); do \
    zcat Ca.$i.fst.txt.gz | grep -v MEAN | awk '{print $8,"'$i'"}' > Ca.$i.fst
done

cat Ca*fst | sed '1i\Fst Group' | sed 's/ /\t/g' > Combined.fst.txt

```

### Fst with accessions no taxonomy reassignment

```bash
# Working directory
cd /data/zhaojiantao/pepper/Fst/No.taxonomy.correction

# Calculate Fst between subgroups
for i in annuum bac.baccatum bac.pendulum chacoense chinense frutescens glabriusculum_I glabriusculum_II; do
    for j in bac.baccatum bac.pendulum chacoense chinense frutescens glabriusculum_I glabriusculum_II pubescens; do
        for chr in $(seq -w 1 12); do
        vcftools --gzvcf ../../vcf.clean/Chr$chr.clean.vcf.gz \
            --weir-fst-pop Ca.$i.list \
            --weir-fst-pop Ca.$j.list \
            --fst-window-size 100000 \
            --fst-window-step 10000 \
            --out Ca.$i.$j.Chr$chr 
        done &
    done
done


# Combine the results
for i in $(ll *weir.fst | awk '{print $9}' | sed 's/.Chr/\t/g' | cut -f1 | sort | uniq); do 
    cat $i.Chr*.windowed.weir.fst | \
        grep -v CHROM | \
        awk '{print $7=$1"_"$2"_"$3,$1,$2,$3,$8=($2+$3)/2,$4,$5,$6}' | \
        awk '$8 > 0 {print $0}' | \
        awk '{gsub(/Chr0/, "", $2)} 1' | \
        awk '{gsub(/Chr/, "", $2)} 1' | \
        sed '1i\ID CHROM BIN_START BIN_END BIN_Middle N_VARIANTS WEIGHTED_FST MEAN_FST' | \
        sed 's/ /\t/g' \
        > $i.fst.txt
done

# Remove chromosomal level files
rm *.windowed.weir.fst *log 

# Calculate the mean Fst value
for i in $(ll *fst.txt | awk '{print $9}' | sed 's/.fst.txt//g'); do
    echo $i; cat $i.fst.txt | grep -v MEAN | awk '{sum+=$8} END {print sum/NR}'
done

# Remove unnecessary files
rm *fst.txt-2 common.id

# Prepare the file for significant test
for i in $(ls *fst.txt.gz | sed 's/.fst.txt.gz//g' | sed 's/Ca.//g'); do \
    zcat Ca.$i.fst.txt.gz | grep -v MEAN | awk '{print $8,"'$i'"}' > Ca.$i.fst
done

cat Ca*fst | sed '1i\Fst Group' | sed 's/ /\t/g' > Combined.fst.txt

```

### Fst with accessions of pure background

```bash
# Working directory
cd /data/zhaojiantao/pepper/Fst/Pure.background

# Calculate Fst between BAB and CHA
for i in $(seq -w 1 12); do
vcftools --gzvcf /data/zhaojiantao/pepper/vcf.clean/Chr$i.clean.vcf.gz \
    --weir-fst-pop /data/zhaojiantao/pepper/Pop.list/Ca.bac.baccatum.list \
    --weir-fst-pop /data/zhaojiantao/pepper/Pop.list/Ca.chacoense.list \
    --fst-window-size 100000 \
    --fst-window-step 10000 \
    --out Ca.bac.baccatum.chacoense.Chr$i &
done

# Combine the results
for i in $(ll *weir.fst | awk '{print $9}' | sed 's/.Chr/\t/g' | cut -f1 | sort | uniq); do 
    cat $i.Chr*.windowed.weir.fst | \
        grep -v CHROM | \
        awk '{print $7=$1"_"$2"_"$3,$1,$2,$3,$8=($2+$3)/2,$4,$5,$6}' | \
        awk '$8 > 0 {print $0}' | \
        awk '{gsub(/Chr0/, "", $2)} 1' | \
        awk '{gsub(/Chr/, "", $2)} 1' | \
        sed '1i\ID CHROM BIN_START BIN_END BIN_Middle N_VARIANTS WEIGHTED_FST MEAN_FST' | \
        sed 's/ /\t/g' \
        > $i.fst.txt
done

# Remove chromosomal level files
rm *.windowed.weir.fst *log 

# Calculate the mean Fst value
for i in $(ll *fst.txt | awk '{print $9}' | sed 's/.fst.txt//g'); do
    echo $i; cat $i.fst.txt | grep -v MEAN | awk '{sum+=$8} END {print sum/NR}'
done

# Remove unnecessary files
rm *fst.txt-2 common.id

# Prepare the file for significant test
for i in $(ls *fst.txt.gz | sed 's/.fst.txt.gz//g' | sed 's/Ca.//g'); do \
    zcat Ca.$i.fst.txt.gz | grep -v MEAN | awk '{print $8,"'$i'"}' > Ca.$i.fst
done

cat Ca*fst | sed '1i\Fst Group' | sed 's/ /\t/g' > Combined.fst.txt

```

### Fst permutation test

#### BAB_CHA

```bash
# Working directory
cd /data/zhaojiantao/pepper/Fst/PermutationTest/BAB_CHA

# Without permutation
for i in $(seq -w 1 12); do
vcftools --vcf /data/zhaojiantao/pepper/SNP.prune/Chr$i.prune.vcf \
    --weir-fst-pop /data/zhaojiantao/pepper/Pop.list/Ca.bac.baccatum.list \
    --weir-fst-pop /data/zhaojiantao/pepper/Pop.list/Ca.chacoense.list \
    --fst-window-size 100000 \
    --fst-window-step 10000 \
    --out Ca.bac.baccatum.chacoense.Chr$i &
done

# Calculate the mean value
cat Ca.bac.baccatum.chacoense.Chr*.windowed.weir.fst | grep -v CHROM | \
    awk '$6>=0 {print}' | awk '{sum+=$6} END {print sum/NR}' > Ca.bac.baccatum.chacoense.meanFst
# 0.411388

# Remove tempory file
rm *log *weir.fst

# Permutation test
wc -l /data/zhaojiantao/pepper/Pop.list/Ca.bac.baccatum.list /data/zhaojiantao/pepper/Pop.list/Ca.chacoense.list

# Prepare the list
cat /data/zhaojiantao/pepper/Pop.list/Ca.bac.baccatum.list /data/zhaojiantao/pepper/Pop.list/Ca.chacoense.list \
    > Ca.bac.baccatum.chacoense.list

# Randomize the samples 
for i in $(seq -w 1 100); do
    sort --random-sort Ca.bac.baccatum.chacoense.list | head -12 > Ca.bac.baccatum.Random$i.list
    grep -Fvf Ca.bac.baccatum.Random$i.list Ca.bac.baccatum.chacoense.list > Ca.chacoense.Random$i.list
done

# Run the permutation test with 100 replicates
for i in $(seq -w 1 100); do
    for j in $(seq -w 1 12); do
    vcftools --vcf /data/zhaojiantao/pepper/SNP.prune/Chr$j.prune.vcf \
        --weir-fst-pop Ca.bac.baccatum.Random$i.list \
        --weir-fst-pop Ca.chacoense.Random$i.list \
        --fst-window-size 100000 \
        --fst-window-step 10000 \
        --out Ca.bac.baccatum.chacoense.Random$i.Chr$j &
    done
done

# Calculate the mean value
for i in $(seq -w 1 100); do
    cat Ca.bac.baccatum.chacoense.Random$i.Chr*.windowed.weir.fst | \
        grep -v CHROM | awk '$6>=0 {print}' | awk '{sum+=$6} END {print sum/NR}' \
        > Ca.bac.baccatum.chacoense.Random$i.meanFst
done

# Combine the permutation result
cat *Random*.meanFst > Ca.bac.baccatum.chacoense.Random.meanFst

# Remove tempory file
rm *log *weir.fst Ca.bac.baccatum.Random*.list Ca.chacoense.Random*.list 
rm Ca.bac.baccatum.chacoense.Random0*.meanFst Ca.bac.baccatum.chacoense.Random100.meanFst

# Generate scripts for significant t test in R
echo "
data <- read.table("Ca.bac.baccatum.chacoense.Random.meanFst", header=FALSE)
AA <- t.test(x=data, mu=0.411388, conf.level=0.95)
chars <- capture.output(print(AA))
writeLines(chars, con = file("t.test.txt")) " | \
sed 's/Ca.bac.baccatum.chacoense.Random.meanFst/"Ca.bac.baccatum.chacoense.Random.meanFst"/g' |\
sed 's/t.test.txt/"t.test.txt"/g' \
> t.test.R

# Significant t test in R
Rscript t.test.R

```

## Nucleotide diversity

```bash
# Workind directory
cd /data/zhaojiantao/pepper/Pi

# Calculate per-site nucleotide diversity
for i in Ca.annuum Ca.bac.baccatum Ca.bac.pendulum Ca.chacoense Ca.chinense Ca.frutescens Ca.glabriusculum Ca.glabriusculum_II ; do
    for j in $(seq -w 1 12); do
        vcftools --gzvcf ../../vcf.clean/Chr$j.clean.vcf.gz --keep $i.list --site-pi --out $i.Chr$j &
    done 
done

# Calculate windowed pi 
for i in Ca.annuum Ca.bac.baccatum Ca.bac.pendulum Ca.chacoense Ca.chinense Ca.frutescens Ca.glabriusculum Ca.glabriusculum_II; do
    for j in $(seq -w 1 12); do
        perl PI_slide_win.pl -w 100000 -s 10000 $i.Chr$j.sites.pi
    done &
done

# Combine the results of each chromosome
for i in Ca.annuum Ca.bac.baccatum Ca.bac.pendulum Ca.chacoense Ca.chinense Ca.frutescens Ca.glabriusculum Ca.glabriusculum_II Ca.pubescens; do
    cat $i.Chr*.w100000.s10000.txt | \
        grep -v "#Chr" | \
        awk '{print $1"_"$3"_"$4,$1,$3,$4,($3+$4)/2,$6,$7}' | \
        awk '{gsub(/Chr0/, "", $2)} 1' | \
        awk '{gsub(/Chr/, "", $2)} 1' | \
        sed '1i\ID CHROM BIN_START BIN_END BIN_Middle N_VARIANTS PI' | \
        sed 's/ /\t/g' \
        > $i.windowed.pi.txt
done

# Find the common ids
awk '{print $1,$2,$3,$7,"Ca.annuum"}' Ca.annuum.windowed.pi.txt > Ca.annuum.windowed.pi.txt-2
awk '{print $1,$2,$3,$7,"Ca.bac.baccatum"}' Ca.bac.baccatum.windowed.pi.txt > Ca.bac.baccatum.windowed.pi.txt-2
awk '{print $1,$2,$3,$7,"Ca.bac.pendulum"}' Ca.bac.pendulum.windowed.pi.txt > Ca.bac.pendulum.windowed.pi.txt-2
awk '{print $1,$2,$3,$7,"Ca.chacoense"}' Ca.chacoense.windowed.pi.txt > Ca.chacoense.windowed.pi.txt-2
awk '{print $1,$2,$3,$7,"Ca.chinense"}' Ca.chinense.windowed.pi.txt > Ca.chinense.windowed.pi.txt-2
awk '{print $1,$2,$3,$7,"Ca.frutescens"}' Ca.frutescens.windowed.pi.txt > Ca.frutescens.windowed.pi.txt-2
awk '{print $1,$2,$3,$7,"Ca.glabriusculum_II"}' Ca.glabriusculum_II.windowed.pi.txt > Ca.glabriusculum_II.windowed.pi.txt-2
awk '{print $1,$2,$3,$7,"Ca.glabriusculum"}' Ca.glabriusculum.windowed.pi.txt > Ca.glabriusculum.windowed.pi.txt-2
awk '{print $1,$2,$3,$7,"Ca.pubescens"}' Ca.pubescens.windowed.pi.txt > Ca.pubescens.windowed.pi.txt-2

cat *windowed.pi.txt-2 | grep -v ID | sed '1i\SNP Chromosome Position PI Taxonomy' | sed 's/ /\t/g' > Group.pi.txt

# Calculate the mean value of Pi
for i in *windowed.pi.txt ; do
    echo $i ;
    grep -v ID $i | awk '{ total += $7 } END { print total/NR }' 
done

cat *Random*.meanFst > Ca.bac.baccatum.chacoense.Random.meanFst

# Remove unnecessary files
rm *.windowed.pi.txt-2 common.snp *windowed.pi *.w100000.s10000.txt *sites.pi *log *

```

### Nucleotide diversity with accessions no taxonomy reassignment

```bash
# Working directory
cd /data/zhaojiantao/pepper/Pi/No.taxonomy.correction

# Calculate per-site nucleotide diversity
for i in Ca.annuum Ca.bac.baccatum Ca.bac.pendulum Ca.chacoense Ca.chinense Ca.frutescens Ca.glabriusculum_I Ca.glabriusculum_II Ca.pubescens ; do
    for j in $(seq -w 1 12); do
        vcftools --gzvcf ../../vcf.clean/Chr$j.clean.vcf.gz --keep $i.list --site-pi --out $i.Chr$j &
    done &
done

# Calculate windowed pi 
for i in Ca.annuum Ca.bac.baccatum Ca.bac.pendulum Ca.chacoense Ca.chinense Ca.frutescens Ca.glabriusculum_I Ca.glabriusculum_II Ca.pubescens ; do
    for j in $(seq -w 1 12); do
        perl PI_slide_win.pl -w 100000 -s 10000 $i.Chr$j.sites.pi
    done &
done

# Combine the results of each chromosome
for i in Ca.annuum Ca.bac.baccatum Ca.bac.pendulum Ca.chacoense Ca.chinense Ca.frutescens Ca.glabriusculum_I Ca.glabriusculum_II Ca.pubescens; do
    cat $i.Chr*.w100000.s10000.txt | \
        grep -v "#Chr" | \
        awk '{print $1"_"$3"_"$4,$1,$3,$4,($3+$4)/2,$6,$7}' | \
        awk '{gsub(/Chr0/, "", $2)} 1' | \
        awk '{gsub(/Chr/, "", $2)} 1' | \
        sed '1i\ID CHROM BIN_START BIN_END BIN_Middle N_VARIANTS PI' | \
        sed 's/ /\t/g' \
        > $i.windowed.pi.txt
done

# Find the common ids
awk '{print $1,$2,$3,$7,"A"}' Ca.pubescens.windowed.pi.txt > Ca.pubescens.windowed.pi.txt-2
awk '{print $1,$2,$3,$7,"B"}' Ca.chacoense.windowed.pi.txt > Ca.chacoense.windowed.pi.txt-2
awk '{print $1,$2,$3,$7,"C"}' Ca.bac.baccatum.windowed.pi.txt > Ca.bac.baccatum.windowed.pi.txt-2
awk '{print $1,$2,$3,$7,"D"}' Ca.bac.pendulum.windowed.pi.txt > Ca.bac.pendulum.windowed.pi.txt-2
awk '{print $1,$2,$3,$7,"E"}' Ca.glabriusculum_I.windowed.pi.txt > Ca.glabriusculum_I.windowed.pi.txt-2
awk '{print $1,$2,$3,$7,"F"}' Ca.glabriusculum_II.windowed.pi.txt > Ca.glabriusculum_II.windowed.pi.txt-2
awk '{print $1,$2,$3,$7,"G"}' Ca.annuum.windowed.pi.txt > Ca.annuum.windowed.pi.txt-2
awk '{print $1,$2,$3,$7,"H"}' Ca.frutescens.windowed.pi.txt > Ca.frutescens.windowed.pi.txt-2
awk '{print $1,$2,$3,$7,"I"}' Ca.chinense.windowed.pi.txt > Ca.chinense.windowed.pi.txt-2

cat *windowed.pi.txt-2 | grep -v ID | sed '1i\SNP Chromosome Position PI Taxonomy' | sed 's/ /\t/g' > Group.pi.nocorrect.txt

# Calculate the mean value of Pi
for i in *windowed.pi.txt ; do
    echo $i ;
    grep -v ID $i | awk '{ total += $7 } END { print total/NR }' 
done

# Calculate the mean value of pi for all 
grep -v SNP Group.pi.nocorrect.txt | awk '{ total += $4 } END { print total/NR }' 

# Clean the files
rm *w100000.s10000.txt *log *-2
for i in Ca.annuum Ca.bac.baccatum Ca.bac.pendulum Ca.chacoense Ca.chinense Ca.frutescens Ca.glabriusculum_I Ca.glabriusculum_II Ca.pubescens; do
    cat $i.Chr*.sites.pi | gzip > $i.sites.pi.gz
done

```

### Nucleotide diversity with accessions of pure background

```bash
# Working directory
cd /data/zhaojiantao/pepper/Pi/Pure.background

# Calculate per-site nucleotide diversity
for i in Ca.annuum Ca.bac.baccatum Ca.bac.pendulum Ca.chacoense Ca.chinense Ca.frutescens Ca.glabriusculum Ca.glabriusculum_II Ca.pubescens ; do
    for j in $(seq -w 1 12); do
        vcftools --gzvcf ../../vcf.clean/Chr$j.clean.vcf.gz --keep $i.list --site-pi --out $i.Chr$j
    done
done

# Calculate windowed pi 
for i in Ca.annuum Ca.bac.baccatum Ca.bac.pendulum Ca.chacoense Ca.chinense Ca.frutescens Ca.glabriusculum_I Ca.glabriusculum_II Ca.pubescens ; do
    for j in $(seq -w 1 12); do
        perl PI_slide_win.pl -w 100000 -s 10000 $i.Chr$j.sites.pi
    done &
done

# Combine the results of each chromosome
for i in Ca.annuum Ca.bac.baccatum Ca.bac.pendulum Ca.chacoense Ca.chinense Ca.frutescens Ca.glabriusculum_I Ca.glabriusculum_II Ca.pubescens; do
    cat $i.Chr*.w100000.s10000.txt | \
        grep -v "#Chr" | \
        awk '{print $1"_"$3"_"$4,$1,$3,$4,($3+$4)/2,$6,$7}' | \
        awk '{gsub(/Chr0/, "", $2)} 1' | \
        awk '{gsub(/Chr/, "", $2)} 1' | \
        sed '1i\ID CHROM BIN_START BIN_END BIN_Middle N_VARIANTS PI' | \
        sed 's/ /\t/g' \
        > $i.windowed.pi.txt
done

# Find the common ids
awk '{print $1,$2,$3,$7,"A"}' Ca.pubescens.windowed.pi.txt > Ca.pubescens.windowed.pi.txt-2
awk '{print $1,$2,$3,$7,"B"}' Ca.chacoense.windowed.pi.txt > Ca.chacoense.windowed.pi.txt-2
awk '{print $1,$2,$3,$7,"C"}' Ca.bac.baccatum.windowed.pi.txt > Ca.bac.baccatum.windowed.pi.txt-2
awk '{print $1,$2,$3,$7,"D"}' Ca.bac.pendulum.windowed.pi.txt > Ca.bac.pendulum.windowed.pi.txt-2
awk '{print $1,$2,$3,$7,"E"}' Ca.glabriusculum_I.windowed.pi.txt > Ca.glabriusculum_I.windowed.pi.txt-2
awk '{print $1,$2,$3,$7,"F"}' Ca.glabriusculum_II.windowed.pi.txt > Ca.glabriusculum_II.windowed.pi.txt-2
awk '{print $1,$2,$3,$7,"G"}' Ca.annuum.windowed.pi.txt > Ca.annuum.windowed.pi.txt-2
awk '{print $1,$2,$3,$7,"H"}' Ca.frutescens.windowed.pi.txt > Ca.frutescens.windowed.pi.txt-2
awk '{print $1,$2,$3,$7,"I"}' Ca.chinense.windowed.pi.txt > Ca.chinense.windowed.pi.txt-2

cat *windowed.pi.txt-2 | grep -v ID | sed '1i\SNP Chromosome Position PI Taxonomy' | sed 's/ /\t/g' > Group.pi.pure.txt

# Calculate the mean value of Pi
for i in *windowed.pi.txt ; do
    echo $i ;
    grep -v ID $i | awk '{ total += $7 } END { print total/NR }' 
done

# Calculate the mean value of pi for all 
grep -v SNP Group.pi.pure.txt | awk '{ total += $4 } END { print total/NR }' 

```

## Treemix analysis

Reference: [treemix analysis](https://speciationgenomics.github.io/Treemix/)

```bash
# Working directory
cd /data/zhaojiantao/pepper/treemix

for i in $(seq -w 1 12); do
    # Step 16.1 Create chromosome map
    echo Chr$i Chr$i | sed 's/ /\t/g' > Chr$i.clean-map.txt ; 
    # Step 16.2 Make ped file using this chromosome map
    vcftools --gzvcf ../vcf.clean/Chr$i.clean.vcf.gz \
        --keep ../Pop.list/Pop8File.txt \
        --max-missing 1 \
        --plink --chrom-map Chr$i.clean-map.txt \
        --out Chr$i.noN
    # Step 16.3 It assumes unlinked SNPs and we are thus first going to prune the file for SNPs in high LD.
    plink --file Chr$i.noN --indep-pairwise 50 10 0.2 --allow-extra-chr --out Chr$i.noN
    plink --file Chr$i.noN --extract Chr$i.noN.prune.in --make-bed --allow-extra-chr --out Chr$i.noN.pruned
    # Step 16.4 Prepare the pop clust file
    cat ../Pop.list/Ca.annuum.list | awk '{print $1,$1,"Ca.annuum"}' > Ca.annuum.clust
    cat ../Pop.list/Ca.bac.baccatum.list | awk '{print $1,$1,"Ca.bac.baccatum"}' > Ca.bac.baccatum.clust
    cat ../Pop.list/Ca.bac.pendulum.list | awk '{print $1,$1,"Ca.bac.pendulum"}' > Ca.bac.pendulum.clust
    cat ../Pop.list/Ca.chacoense.list | awk '{print $1,$1,"Ca.chacoense"}' > Ca.chacoense.clust
    cat ../Pop.list/Ca.chinense.list | awk '{print $1,$1,"Ca.chinense"}' > Ca.chinense.clust
    cat ../Pop.list/Ca.frutescens.list | awk '{print $1,$1,"Ca.frutescens"}' > Ca.frutescens.clust
    cat ../Pop.list/Ca.glabriusculum.list | awk '{print $1,$1,"Ca.glabriusculum"}' > Ca.glabriusculum.clust
    cat ../Pop.list/Ca.pubescens.list | awk '{print $1,$1,"Ca.pubescens"}' > Ca.pubescens.clust
    cat *clust > Pop8.clust.txt
    rm *clust
    # Step 16.5 Calculate the frequency
    plink -bfile Chr$i.noN.pruned --freq --missing --allow-extra-chr --within Pop8.clust.txt --out Chr$i.noN.pruned
done

# Step 16.6 Prepare the input file for treemix
cat Chr01.noN.pruned.frq.strat | sed -n '1p' > header
cat Chr*.noN.pruned.frq.strat | grep -v CHR > noN.pruned.frq-1
cat header noN.pruned.frq-1 | gzip > noN.pruned.frq.gz

# Step 16.7 Convert frequency file to treemix file format
plink2treemix.py noN.pruned.frq.gz noN.pruned.treemix.frq.gz

# Step 16.8 Remove unnecessary files
rm *clean-map.txt *bed *bim *fam *log *map *nosex *ped *frq.strat *imiss *lmiss *prune.in *prune.out header noN.pruned.frq-1 

# Step 16.9 run treemix
cd /data/zhaojiantao/pepper/treemix/migration

for i in $(seq -w 1 10); do  
    treemix \
        -i ../Ca.S8.noN.pruned.treemix.frq.gz \
        -m $i \
        -o Ca.S8.noN.pruned.treemix.$i \
        -root Ca.pubescens \
        -bootstrap \
        -k 1400 &
done

# Find the optimal edges using [OptM R package](https://cran.r-project.org/web/packages/OptM/readme/README.html)
# The PL and BC model are parametric models that 
for m in {1..10}; do
    for i in {1..10}; do
        # Generate random seed
        s=$RANDOM
        treemix -i Ca.S8.noN.pruned.treemix.frq.gz -o OptM/S8.${i}.${m} -m ${m} -k 1400 -seed $s -bootstrap 100 -root Ca.pubescens
    done
done

test.linear = optM("OptM", method = "linear")
jpeg("treemix.migration.linear.jpeg", height = 15, width = 20, units="cm", res = 300)
plot_optM(test.linear, method = "linear")
dev.off()

for i in {1..1000}; do
    # Generate random seed
    s=$RANDOM
    treemix -i Ca.S8.noN.pruned.treemix.frq.gz -o OptM2/S8.${i}.1 -m 1 -k 1400 -seed $s -bootstrap 100 -root Ca.pubescens
done

```

## ABBA-BABA tests

```bash
# Working directory
cd /data/zhaojiantao/pepper/ABBA-BABA

# Convert vcf to ABAB-ABBA geno file format for BAC2FRU+CHN
# Extract and convert the analyzed samples
for i in $(seq -w 1 12); do
    python /data/zhaojiantao/tools/ABBA-BABA/genomics_general-master/VCF_processing/parseVCF.py \
        -i ../Chr$i.clean.vcf.gz --skipIndels --minQual 30 --gtf flag=DP |\ 
        bgzip -c > Chr$i.clean.geno.gz 
done


for i in $(seq -w 1 12); do
    python /data/zhaojiantao/tools/ABBA-BABA/genomics_general-master/VCF_processing/parseVCF.py \
        -i ../Chr$i.clean.vcf.gz --skipIndels --minQual 30 --gtf flag=DP | \
        bgzip -c > Chr$i.clean.geno.gz 
done
# Combine the files
zcat Chr01.clean.geno.gz | head -1 > header
zcat Chr*.clean.geno.gz | grep -v "#CHROM" > clean.geno
cat header clean.geno | gzip > clean.geno.gz

# Remove unnecessary files
rm header clean.geno Chr*.clean.geno.gz 
rm *clean-map.txt 

# Prepare the popfile txt
cat ../treemix/Pop8.clust.txt | cut -f 1,3 -d ' ' | sed 's/ /\t/g' > Pop8File.txt
# Introgression between BAB to FRU/CHN
cat Pop8File.txt | sed 's/Ca.bac.baccatum/Ca.baccatum/g' | sed 's/Ca.bac.pendulum/Ca.baccatum/g' | \
    sed 's/Ca.chinense/Ca.FRU.CHN/g' | sed 's/Ca.frutescens/Ca.FRU.CHN/g' > Pop.BAB2ANN.txt

# Calculate the allele frequency
python2 /data/zhaojiantao/tools/ABBA-BABA/genomics_general-master/freq.py -t 60 \
    -g Ca.S8.clean.geno.gz \
    -p Ca.annuum -p Ca.baccatum \
    -p Ca.FRU.CHN -p Ca.pubescens \
    --popsFile Pop.BAB2ANN.txt --target derived \
    -o BAC2FRUCHN.derFreq.tsv.gz

```

## Introgression

```bash
# Working directory
cd /data/zhaojiantao/pepper/ABBA-BABA

# Extract and filter the SNPs
plink --vcf ../vcf.clean/Chr06.clean.vcf.gz  \
    --keep <(awk '{print $1,$1}' Pop.PUB.BAB2ANN.txt) \
    --allow-extra-chr --chr Chr06 --geno 0.1 --maf 0.05 \
    --from-bp 241500000 --to-bp 247500000 \
    --recode -output-missing-genotype 9 --double-id --threads 40 \
    --out Chr06.region
awk '{print $1,$1"_"$4,$3,$4}' Chr06.region.map | sed 's/ /\t/g' > Chr06.region.map2
mv Chr06.region.map2 Chr06.region.map

# Calculate freq for BAC
plink --file Chr06.region \
    --keep <(cat ../Pop.list/Ca.bac*.list | awk '{print $1,$1}') --double-id \
    --freq --out Chr06.region.BAC

# Calculate frequency for PUB
plink --file Chr06.region \
    --keep <(cat ../Pop.list/Ca.pubescens.list | awk '{print $1,$1}') --double-id \
    --freq --out Chr06.region.PUB

# Find the SNPs where PUB.dominant not equal to BAC.dominant allele
paste Chr06.region.PUB.frq Chr06.region.BAC.frq | awk '{print $4,1-$5,$10,1-$11}' | \
    grep -v A2 | paste Chr06.region.map - | sed 's/ /\t/g' | \
    awk '{if ($5==$7)print $0,"TRUE "; else print $0,"FALSE "}' | grep FALSE | cut -f2 > Chr06.region.remain.snps

# Transpose to tped
plink --file Chr06.region --extract Chr06.region.remain.snps \
    --recode transpose --out Chr06.region

# Combine every two lines
cut -f5- -d ' ' Chr06.region.tped | sed 's! \([^ ]\+\)\( \|$\)!\1 !g' > Chr06.region.geno-1

# Chr06.region.BAC.dominant.allele
grep -v CHR Chr06.region.BAC.frq | awk '{print $4$4}' > Chr06.region.BAC.dominant.allele.out
paste Chr06.region.PUB.frq Chr06.region.BAC.frq | awk '{print $4,1-$5,$10,1-$11}' | \
    grep -v A2 | paste Chr06.region.map - | sed 's/ /\t/g' | \
    awk '{if ($5==$7)print $0,"TRUE "; else print $0,"FALSE "}' | grep FALSE | cut -f7 \
    > Chr06.region.BAC.dominant.allele.out1

# Check the SNP IDs
paste Chr06.region.BAC.dominant.allele.out Chr06.region.geno-1 | \
    awk '{for(i=2;i<=449;i=i+1) if ($i==$1){printf "1 "} else {printf "0 "}; printf "\n"}' \
    > Chr06.region.alleles

# Change heterozygous to 0.5
for i in $(seq -w 1 448); do
    cut -f $i -d ' ' Chr06.region.geno-1 | sed 's/./& /g' | \
    paste Chr06.region.BAC.dominant.allele.out1 - | \
    awk '{
        if ($2==$1 && $3==$1) 
        print 1 ;
        else if ($2==$1 && $3!=$1)
        print 0.5
        else 
        print 0
    }' > Accession$i.geno
done
paste Accession*.geno | sed 's/\t/ /g' > Chr06.region.alleles
rm Accession*.geno

# Transpose the genotypes
python -c "import sys; print('\n'.join(' '.join(c) for c in zip(*(l.split() for l in sys.stdin.readlines() if l.strip()))))" < Chr06.region.alleles > Chr06.region.alleles.trans

# Sort the genos
paste Chr06.sample.id Chr06.region.alleles.trans | \
    sed 's/\t/ /g' | \
    sort -k3,3 -t ' ' | \
    cut -f1,2,4- -d ' ' \
    > Chr06.region.introgression 

# Transpose the snp infos
cut -f4 -d ' ' Chr06.region.tped | sed '1i\ID\nTaxonomy' > Chr06.region.pos
python -c "import sys; print('\n'.join(' '.join(c) for c in zip(*(l.split() for l in sys.stdin.readlines() if l.strip()))))" < Chr06.region.pos > Chr06.region.pos.trans

# Combine the matrix
cat Chr06.region.pos.trans Chr06.region.introgression > Chr06.region.introgression.txt
# Transpose the data
python -c "import sys; print('\n'.join(' '.join(c) for c in zip(*(l.split() for l in sys.stdin.readlines() if l.strip()))))" < Chr06.region.introgression.txt > Chr06.region.introgression.trans.txt

grep -v ID Chr06.region.introgression.trans.txt | \
    grep -v Taxonomy | \
    cut -f2- -d ' ' \
    > Chr06.region.introgression.trans2

for i in $(seq 0 35472); do echo $i; done | \
    paste - <(grep -v ID Chr06.region.introgression.trans.txt | grep -v Taxonomy) | awk '{print 6,$1,1,$2,0}' > Chr06.region.introgression.trans3

# Cluster the result in sliding window
for i in $(seq -w 1 448); do
    perl cluster_introgression.pl \
        <(cut -f$i -d ' ' Chr06.region.introgression.trans2 | \
        paste - Chr06.region.introgression.trans3 | awk '{print $2,$3,$4,$5,$6,$1}') \
        -wind_size 1000000 -wind_step 100000 \
        > Chr06.region.w1000000.s100000.Accession$i &
done

# Combine the result
paste Chr06.region.w1000000.s100000.Accession* > Chr06.region.w1000000.s100000.Accessions.region
python -c "import sys; print('\n'.join(' '.join(c) for c in zip(*(l.split() for l in sys.stdin.readlines() if l.strip()))))" < Chr06.region.w1000000.s100000.Accessions.region | grep Avg | grep -v AdjAvg > Chr06.region.w1000000.s100000.Accessions.region.trans

for i in $(seq -w 2 448); do
    rm Chr06.region.w1000000.s100000.Accession$i
done

# Original values
python -c "import sys; print('\n'.join(' '.join(c) for c in zip(*(l.split() for l in sys.stdin.readlines() if l.strip()))))" \
    < Chr06.region.w1000000.s100000.Accessions.region.trans | \
    grep -v Avg >Chr06.region.w1000000.s100000.Accessions.region.clean-1

awk '{print $1,$2,$3,($2+$3)/2}' Chr06.region.w1000000.s100000.Accession001 | \
    grep -v chrID | paste - Chr06.region.w1000000.s100000.Accessions.region.clean-1 | sed 's/ /\t/g' \
    > Chr06.region.w1000000.s100000.Accessions.region.clean.txt

python -c "import sys; print('\n'.join(' '.join(c) for c in zip(*(l.split() for l in sys.stdin.readlines() if l.strip()))))" \
    <Chr06.region.w1000000.s100000.Accessions.region.clean.txt | sed 's/ /\t/g' \
    > Chr06.region.w1000000.s100000.Accessions.region.clean.trans.txt

```


### Introgression Chr07:4,900,001-10,000,000

```bash
# Working directory
cd /data/zhaojiantao/pepper/ABBA-BABA

# Extract and filter the SNPs
plink --vcf ../vcf.clean/Chr07.clean.vcf.gz  \
    --keep <(awk '{print $1,$1}' Pop.PUB.BAB2ANN.txt) \
    --allow-extra-chr --chr Chr07 --geno 0.1 --maf 0.05 \
    --from-bp 3900001 --to-bp 13000000 \
    --recode -output-missing-genotype 9 --double-id \
    --out Chr07.region1
awk '{print $1,$1"_"$4,$3,$4}' Chr07.region1.map | sed 's/ /\t/g' > Chr07.region1.map2
mv Chr07.region1.map2 Chr07.region1.map

# Calculate freq for BAC
plink --file Chr07.region1 \
    --keep <(cat ../Pop.list/Ca.bac*.list | awk '{print $1,$1}') --double-id \
    --freq --out Chr07.region1.BAC

# Calculate frequency for PUB
plink --file Chr07.region1 \
    --keep <(cat ../Pop.list/Ca.pubescens.list | awk '{print $1,$1}') --double-id \
    --freq --out Chr07.region1.PUB

# Find the SNPs where PUB.dominant not equal to BAC.dominant allele
paste Chr07.region1.PUB.frq Chr07.region1.BAC.frq | awk '{print $4,1-$5,$10,1-$11}' | \
    grep -v A2 | paste Chr07.region1.map - | sed 's/ /\t/g' | \
    awk '{if ($5==$7)print $0,"TRUE "; else print $0,"FALSE "}' | grep FALSE | cut -f2 > Chr07.region1.remain.snps

# Transpose to tped
plink --file Chr07.region1 --extract Chr07.region1.remain.snps \
    --recode transpose --out Chr07.region1

# Combine every two lines
cut -f5- -d ' ' Chr07.region1.tped | sed 's! \([^ ]\+\)\( \|$\)!\1 !g' > Chr07.region1.geno-1

# Chr07.region1.BAC.dominant.allele
grep -v CHR Chr07.region1.BAC.frq | awk '{print $4$4}' > Chr07.region1.BAC.dominant.allele.out
paste Chr07.region1.PUB.frq Chr07.region1.BAC.frq | awk '{print $4,1-$5,$10,1-$11}' | \
    grep -v A2 | paste Chr07.region1.map - | sed 's/ /\t/g' | \
    awk '{if ($5==$7)print $0,"TRUE "; else print $0,"FALSE "}' | grep FALSE | cut -f7 \
    > Chr07.region1.BAC.dominant.allele.out1

# Check the SNP IDs
#paste Chr07.region1.BAC.dominant.allele.out Chr07.region1.geno-1 | \
#    awk '{for(i=2;i<=449;i=i+1) if ($i==$1){printf "1 "} else {printf "0 "}; printf "\n"}' \
#    > Chr07.region1.alleles

# Change heterozygous to 0.5
for i in $(seq -w 1 448); 
    do cut -f $i -d ' ' Chr07.region1.geno-1 | sed 's/./& /g' | \
    paste Chr07.region1.BAC.dominant.allele.out1 - | \
    awk '{
        if ($2==$1 && $3==$1) 
        print 1 ;
        else if ($2==$1 && $3!=$1)
        print 0.5
        else 
        print 0
    }' > Accession$i.geno &
done
paste Accession*.geno | sed 's/\t/ /g' > Chr07.region1.alleles
rm Accession*.geno

# Transpose the genotypes
python -c "import sys; print('\n'.join(' '.join(c) for c in zip(*(l.split() for l in sys.stdin.readlines() if l.strip()))))" < Chr07.region1.alleles > Chr07.region1.alleles.trans

# Sort the genos
paste Chr07.sample.id Chr07.region1.alleles.trans | sed 's/\t/ /g' | sort -k3 | cut -f1,2,4- -d ' ' > Chr07.region1.introgression 

# Transpose the snp infos
cut -f4 -d ' ' Chr07.region1.tped | sed '1i\ID\nTaxonomy' > Chr07.region1.pos
python -c "import sys; print('\n'.join(' '.join(c) for c in zip(*(l.split() for l in sys.stdin.readlines() if l.strip()))))" < Chr07.region1.pos > Chr07.region1.pos.trans

# Combine the matrix
cat Chr07.region1.pos.trans Chr07.region1.introgression > Chr07.region1.introgression.txt
# Transpose the data
python -c "import sys; print('\n'.join(' '.join(c) for c in zip(*(l.split() for l in sys.stdin.readlines() if l.strip()))))" < Chr07.region1.introgression.txt > Chr07.region1.introgression.trans.txt

grep -v ID Chr07.region1.introgression.trans.txt | grep -v Taxonomy | cut -f2- -d ' ' > Chr07.region1.introgression.trans2

for i in $(seq 0 38089); do echo $i; done | \
    paste - <(grep -v ID Chr07.region1.introgression.trans.txt | grep -v Taxonomy) | awk '{print 6,$1,1,$2,0}' > Chr07.region1.introgression.trans3

# Cluster the result in sliding window
for i in $(seq -w 1 448); do
    perl cluster_introgression.pl \
        <(cut -f$i -d ' ' Chr07.region1.introgression.trans2 | \
        paste - Chr07.region1.introgression.trans3 | awk '{print $2,$3,$4,$5,$6,$1}') \
        -wind_size 1000000 -wind_step 100000 \
        > Chr07.region1.w1000000.s100000.Accession$i &
done

# Combine the result
paste Chr07.region1.w1000000.s100000.Accession* > Chr07.region1.w1000000.s100000.Accessions.region1
python -c "import sys; print('\n'.join(' '.join(c) for c in zip(*(l.split() for l in sys.stdin.readlines() if l.strip()))))" < Chr07.region1.w1000000.s100000.Accessions.region1 | grep Avg | grep -v AdjAvg > Chr07.region1.w1000000.s100000.Accessions.region1.trans

for i in $(seq -w 2 448); do
    rm Chr07.region1.w1000000.s100000.Accession$i
done

# Original values
python -c "import sys; print('\n'.join(' '.join(c) for c in zip(*(l.split() for l in sys.stdin.readlines() if l.strip()))))" \
    < Chr07.region1.w1000000.s100000.Accessions.region1.trans | \
    grep -v Avg >Chr07.region1.w1000000.s100000.Accessions.region1.clean-1

awk '{print $1,$2,$3,($2+$3)/2}' Chr07.region1.w1000000.s100000.Accession001 | \
    grep -v chrID | paste - Chr07.region1.w1000000.s100000.Accessions.region1.clean-1 | sed 's/ /\t/g' \
    > Chr07.region1.w1000000.s100000.Accessions.region1.clean.txt

python -c "import sys; print('\n'.join(' '.join(c) for c in zip(*(l.split() for l in sys.stdin.readlines() if l.strip()))))" \
    <Chr07.region1.w1000000.s100000.Accessions.region1.clean.txt | sed 's/ /\t/g' \
    > Chr07.region1.w1000000.s100000.Accessions.region1.clean.trans.txt

# Prepare the plot of fd
sed 's/,/\t/g' PUB.BAC.FRUCHN.ANN.out | \
    awk '$1=="Ca.S8.Chr07" && $2>=3900001 && $3<=13000000 {print $0}' | grep -v "nan" | \
    cut -f1-4,10 > PUB.BAC.FRUCHN.ANN.Chr07.region1.out.txt

```


### Introgression Chr07:258,150,001-259,300,000

```bash
# Working directory
cd /data/zhaojiantao/pepper/ABBA-BABA

# Extract and filter the SNPs
plink --vcf ../vcf.clean/Chr07.clean.vcf.gz  \
    --keep <(awk '{print $1,$1}' Pop.PUB.BAB2ANN.txt) \
    --allow-extra-chr --chr Chr07 --geno 0.1 --maf 0.05 \
    --from-bp 256150001 --to-bp 261300000 \
    --recode -output-missing-genotype 9 --double-id \
    --out Chr07.region2
awk '{print $1,$1"_"$4,$3,$4}' Chr07.region2.map | sed 's/ /\t/g' > Chr07.region2.map2
mv Chr07.region2.map2 Chr07.region2.map

# Calculate freq for BAC
plink --file Chr07.region2 \
    --keep <(cat ../Pop.list/Ca.bac*.list | awk '{print $1,$1}') --double-id \
    --freq --out Chr07.region2.BAC

# Calculate frequency for PUB
plink --file Chr07.region2 \
    --keep <(cat ../Pop.list/Ca.pubescens.list | awk '{print $1,$1}') --double-id \
    --freq --out Chr07.region2.PUB

# Find the SNPs where PUB.dominant not equal to BAC.dominant allele
paste Chr07.region2.PUB.frq Chr07.region2.BAC.frq | awk '{print $4,1-$5,$10,1-$11}' | \
    grep -v A2 | paste Chr07.region2.map - | sed 's/ /\t/g' | \
    awk '{if ($5==$7)print $0,"TRUE "; else print $0,"FALSE "}' | grep FALSE | cut -f2 > Chr07.region2.remain.snps

# Transpose to tped
plink --file Chr07.region2 --extract Chr07.region2.remain.snps \
    --recode transpose --out Chr07.region2

# Combine every two lines
cut -f5- -d ' ' Chr07.region2.tped | sed 's! \([^ ]\+\)\( \|$\)!\1 !g' > Chr07.region2.geno-1

# Chr07.region2.BAC.dominant.allele
grep -v CHR Chr07.region2.BAC.frq | awk '{print $4$4}' > Chr07.region2.BAC.dominant.allele.out
paste Chr07.region2.PUB.frq Chr07.region2.BAC.frq | awk '{print $4,1-$5,$10,1-$11}' | \
    grep -v A2 | paste Chr07.region2.map - | sed 's/ /\t/g' | \
    awk '{if ($5==$7)print $0,"TRUE "; else print $0,"FALSE "}' | grep FALSE | cut -f7 \
    > Chr07.region2.BAC.dominant.allele.out1

# Check the SNP IDs
paste Chr07.region2.BAC.dominant.allele.out Chr07.region2.geno-1 | \
    awk '{for(i=2;i<=449;i=i+1) if ($i==$1){printf "1 "} else {printf "0 "}; printf "\n"}' \
    > Chr07.region2.alleles

# Change heterozygous to 0.5
for i in $(seq -w 1 448); do
    cut -f $i -d ' ' Chr07.region2.geno-1 | sed 's/./& /g' | \
    paste Chr07.region2.BAC.dominant.allele.out1 - | \
    awk '{
        if ($2==$1 && $3==$1) 
        print 1 ;
        else if ($2==$1 && $3!=$1)
        print 0.5
        else 
        print 0
    }' > Accession$i.geno &
done

paste Accession*.geno | sed 's/\t/ /g' > Chr07.region2.alleles
rm Accession*.geno

# Transpose the genotypes
python -c "import sys; print('\n'.join(' '.join(c) for c in zip(*(l.split() for l in sys.stdin.readlines() if l.strip()))))" < Chr07.region2.alleles > Chr07.region2.alleles.trans

# Sort the genos
paste Chr07.sample.id Chr07.region2.alleles.trans | sed 's/\t/ /g' | sort -k3 | cut -f1,2,4- -d ' ' > Chr07.region2.introgression 

# Transpose the snp infos
cut -f4 -d ' ' Chr07.region2.tped | sed '1i\ID\nTaxonomy' > Chr07.region2.pos
python -c "import sys; print('\n'.join(' '.join(c) for c in zip(*(l.split() for l in sys.stdin.readlines() if l.strip()))))" < Chr07.region2.pos > Chr07.region2.pos.trans

# Combine the matrix
cat Chr07.region2.pos.trans Chr07.region2.introgression > Chr07.region2.introgression.txt
# Transpose the data
python -c "import sys; print('\n'.join(' '.join(c) for c in zip(*(l.split() for l in sys.stdin.readlines() if l.strip()))))" < Chr07.region2.introgression.txt > Chr07.region2.introgression.trans.txt

grep -v ID Chr07.region2.introgression.trans.txt | grep -v Taxonomy | cut -f2- -d ' ' > Chr07.region2.introgression.trans2

for i in $(seq 0 26308); do echo $i; done | \
    paste - <(grep -v ID Chr07.region2.introgression.trans.txt | grep -v Taxonomy) | awk '{print 6,$1,1,$2,0}' > Chr07.region2.introgression.trans3

# Cluster the result in sliding window
for i in $(seq -w 1 448); do
    perl cluster_introgression.pl \
        <(cut -f$i -d ' ' Chr07.region2.introgression.trans2 | \
        paste - Chr07.region2.introgression.trans3 | awk '{print $2,$3,$4,$5,$6,$1}') \
        -wind_size 1000000 -wind_step 100000 \
        > Chr07.region2.w1000000.s100000.Accession$i &
done

# Combine the result
paste Chr07.region2.w1000000.s100000.Accession* > Chr07.region2.w1000000.s100000.Accessions.region2
python -c "import sys; print('\n'.join(' '.join(c) for c in zip(*(l.split() for l in sys.stdin.readlines() if l.strip()))))" < Chr07.region2.w1000000.s100000.Accessions.region2 | grep Avg | grep -v AdjAvg > Chr07.region2.w1000000.s100000.Accessions.region2.trans

for i in $(seq -w 2 448); do
    rm Chr07.region2.w1000000.s100000.Accession$i
done

# Original values
python -c "import sys; print('\n'.join(' '.join(c) for c in zip(*(l.split() for l in sys.stdin.readlines() if l.strip()))))" \
    < Chr07.region2.w1000000.s100000.Accessions.region2.trans | \
    grep -v Avg >Chr07.region2.w1000000.s100000.Accessions.region2.clean-1

awk '{print $1,$2,$3,($2+$3)/2}' Chr07.region2.w1000000.s100000.Accession001 | \
    grep -v chrID | paste - Chr07.region2.w1000000.s100000.Accessions.region2.clean-1 | sed 's/ /\t/g' \
    > Chr07.region2.w1000000.s100000.Accessions.region2.clean.txt

python -c "import sys; print('\n'.join(' '.join(c) for c in zip(*(l.split() for l in sys.stdin.readlines() if l.strip()))))" \
    <Chr07.region2.w1000000.s100000.Accessions.region2.clean.txt | sed 's/ /\t/g' \
    > Chr07.region2.w1000000.s100000.Accessions.region2.clean.trans.txt

# Prepare the plot of fd
sed 's/,/\t/g' PUB.BAC.FRUCHN.ANN.out | \
    awk '$1=="Ca.S8.Chr07" && $2>=256150001 && $3<=261300000 {print $0}' | grep -v "nan" | \
    cut -f1-4,10 > PUB.BAC.FRUCHN.ANN.Chr07.region2.out.txt

```

### Introgression Chr12:250550001-255900000

```bash
# Working directory
cd /data/zhaojiantao/pepper/ABBA-BABA
# Extract and filter the SNPs
plink --vcf ../vcf.clean/Chr12.clean.vcf.gz  \
    --keep <(awk '{print $1,$1}' Pop.PUB.BAB2ANN.txt) \
    --allow-extra-chr --chr Chr12 --geno 0.1 --maf 0.05 \
    --from-bp 250550001 --to-bp 255900000 \
    --recode -output-missing-genotype 9 --double-id \
    --out Chr12.region1
awk '{print $1,$1"_"$4,$3,$4}' Chr12.region1.map | sed 's/ /\t/g' > Chr12.region1.map2
mv Chr12.region1.map2 Chr12.region1.map

# Calculate freq for BAC
plink --file Chr12.region1 \
    --keep <(cat ../Pop.list/Ca.bac*.list | awk '{print $1,$1}') --double-id \
    --freq --out Chr12.region1.BAC

# Calculate frequency for PUB
plink --file Chr12.region1 \
    --keep <(cat ../Pop.list/Ca.pubescens.list | awk '{print $1,$1}') --double-id \
    --freq --out Chr12.region1.PUB

# Find the SNPs where PUB.dominant not equal to BAC.dominant allele
paste Chr12.region1.PUB.frq Chr12.region1.BAC.frq | awk '{print $4,1-$5,$10,1-$11}' | \
    grep -v A2 | paste Chr12.region1.map - | sed 's/ /\t/g' | \
    awk '{if ($5==$7)print $0,"TRUE "; else print $0,"FALSE "}' | grep FALSE | cut -f2 > Chr12.region1.remain.snps

# Transpose to tped
plink --file Chr12.region1 --extract Chr12.region1.remain.snps \
    --recode transpose --out Chr12.region1

# Combine every two lines
cut -f5- -d ' ' Chr12.region1.tped | sed 's! \([^ ]\+\)\( \|$\)!\1 !g' > Chr12.region1.geno-1

# Chr12.region1.BAC.dominant.allele
grep -v CHR Chr12.region1.BAC.frq | awk '{print $4$4}' > Chr12.region1.BAC.dominant.allele.out
paste Chr12.region1.PUB.frq Chr12.region1.BAC.frq | awk '{print $4,1-$5,$10,1-$11}' | \
    grep -v A2 | paste Chr12.region1.map - | sed 's/ /\t/g' | \
    awk '{if ($5==$7)print $0,"TRUE "; else print $0,"FALSE "}' | grep FALSE | cut -f7 \
    > Chr12.region1.BAC.dominant.allele.out1

# Check the SNP IDs
paste Chr12.region1.BAC.dominant.allele.out Chr12.region1.geno-1 | \
    awk '{for(i=2;i<=449;i=i+1) if ($i==$1){printf "1 "} else {printf "0 "}; printf "\n"}' \
    > Chr12.region1.alleles

# Change heterozygous to 0.5
for i in $(seq -w 1 448); 
    do cut -f $i -d ' ' Chr12.region1.geno-1 | sed 's/./& /g' | \
    paste Chr12.region1.BAC.dominant.allele.out1 - | \
    awk '{
        if ($2==$1 && $3==$1) 
        print 1 ;
        else if ($2==$1 && $3!=$1)
        print 0.5
        else 
        print 0
    }' > Accession$i.geno &
done

paste Accession*.geno | sed 's/\t/ /g' > Chr12.region1.alleles
rm Accession*.geno

# Transpose the genotypes
python -c "import sys; print('\n'.join(' '.join(c) for c in zip(*(l.split() for l in sys.stdin.readlines() if l.strip()))))" < Chr12.region1.alleles > Chr12.region1.alleles.trans

# Sort the genos
paste Chr12.sample.id Chr12.region1.alleles.trans | sed 's/\t/ /g' | sort -k3 | cut -f1,2,4- -d ' ' > Chr12.region1.introgression 

# Transpose the snp infos
cut -f4 -d ' ' Chr12.region1.tped | sed '1i\ID\nTaxonomy' > Chr12.region1.pos
python -c "import sys; print('\n'.join(' '.join(c) for c in zip(*(l.split() for l in sys.stdin.readlines() if l.strip()))))" < Chr12.region1.pos > Chr12.region1.pos.trans

# Combine the matrix
cat Chr12.region1.pos.trans Chr12.region1.introgression > Chr12.region1.introgression.txt
# Transpose the data
python -c "import sys; print('\n'.join(' '.join(c) for c in zip(*(l.split() for l in sys.stdin.readlines() if l.strip()))))" < Chr12.region1.introgression.txt > Chr12.region1.introgression.trans.txt

grep -v ID Chr12.region1.introgression.trans.txt | grep -v Taxonomy | cut -f2- -d ' ' > Chr12.region1.introgression.trans2

for i in $(seq 0 26308); do echo $i; done | \
    paste - <(grep -v ID Chr12.region1.introgression.trans.txt | grep -v Taxonomy) | awk '{print 6,$1,1,$2,0}' > Chr12.region1.introgression.trans3

# Cluster the result in sliding window
for i in $(seq -w 1 448); do
    perl cluster_introgression.pl \
        <(cut -f$i -d ' ' Chr12.region1.introgression.trans2 | \
        paste - Chr12.region1.introgression.trans3 | awk '{print $2,$3,$4,$5,$6,$1}') \
        -wind_size 1000000 -wind_step 100000 \
        > Chr12.region1.w1000000.s100000.Accession$i &
done

# Combine the result
paste Chr12.region1.w1000000.s100000.Accession* > Chr12.region1.w1000000.s100000.Accessions.region1
python -c "import sys; print('\n'.join(' '.join(c) for c in zip(*(l.split() for l in sys.stdin.readlines() if l.strip()))))" < Chr12.region1.w1000000.s100000.Accessions.region1 | grep Avg | grep -v AdjAvg > Chr12.region1.w1000000.s100000.Accessions.region1.trans

for i in $(seq -w 2 448); do
    rm Chr12.region1.w1000000.s100000.Accession$i
done

# Original values
python -c "import sys; print('\n'.join(' '.join(c) for c in zip(*(l.split() for l in sys.stdin.readlines() if l.strip()))))" \
    < Chr12.region1.w1000000.s100000.Accessions.region1.trans | \
    grep -v Avg >Chr12.region1.w1000000.s100000.Accessions.region1.clean-1

awk '{print $1,$2,$3,($2+$3)/2}' Chr12.region1.w1000000.s100000.Accession001 | \
    grep -v chrID | paste - Chr12.region1.w1000000.s100000.Accessions.region1.clean-1 | sed 's/ /\t/g' \
    > Chr12.region1.w1000000.s100000.Accessions.region1.clean.txt

python -c "import sys; print('\n'.join(' '.join(c) for c in zip(*(l.split() for l in sys.stdin.readlines() if l.strip()))))" \
    <Chr12.region1.w1000000.s100000.Accessions.region1.clean.txt | sed 's/ /\t/g' \
    > Chr12.region1.w1000000.s100000.Accessions.region1.clean.trans.txt

# Prepare the plot of fd
sed 's/,/\t/g' PUB.BAC.FRUCHN.ANN.out | \
    awk '$1=="Ca.S8.Chr12" && $2>=250550001 && $3<=255900000 {print $0}' | grep -v "nan" | \
    cut -f1-4,10 > PUB.BAC.FRUCHN.ANN.Chr12.region1.out.txt

```

### PUB.BAC.FRUCHN.GLA

```bash
ABBABABAwindows.py -g Ca.S8.clean.geno.gz -o PUB.BAC.FRUCHN.ANN.out -w 500000 -s 50000 -m 100 -f phased --minData 0.5 -T 30 \
    -P1 Ca.annuum -P2 Ca.FRU.CHN -P3 Ca.baccatum -O Ca.pubescens --popsFile Pop.BAB2ANN.txt --writeFailedWindows
# PUB.BAB.CHN.ANN
ABBABABAwindows.py -g Ca.S8.clean.geno.gz -o PUB.BAB.CHN.ANN.out -w 500000 -s 50000 -m 100 -f phased --minData 0.5 -T 20 \
    -P1 Ca.annuum -P2 Ca.chinense -P3 Ca.bac.baccatum -O Ca.pubescens --popsFile Pop8File.txt --writeFailedWindows
# PUB.BAB.FRU.ANN
ABBABABAwindows.py -g Ca.S8.clean.geno.gz -o PUB.BAB.FRU.ANN.out -w 500000 -s 50000 -m 100 -f phased --minData 0.5 -T 20 \
    -P1 Ca.annuum -P2 Ca.frutescens -P3 Ca.bac.baccatum -O Ca.pubescens --popsFile Pop8File.txt --writeFailedWindows
# PUB.BAP.CHN.ANN
ABBABABAwindows.py -g Ca.S8.clean.geno.gz -o PUB.BAP.CHN.ANN.out -w 500000 -s 50000 -m 100 -f phased --minData 0.5 -T 10 \
    -P1 Ca.annuum -P2 Ca.chinense -P3 Ca.bac.pendulum -O Ca.pubescens --popsFile Pop8File.txt --writeFailedWindows 
# PUB.BAP.FRU.ANN
ABBABABAwindows.py -g Ca.S8.clean.geno.gz -o PUB.BAP.FRU.ANN.out -w 500000 -s 50000 -m 100 -f phased --minData 0.5 -T 10 \
    -P1 Ca.annuum -P2 Ca.frutescens -P3 Ca.bac.pendulum -O Ca.pubescens --popsFile Pop8File.txt --writeFailedWindows 

# Select the top 5% highest fd values
sed 's/,/\t/g' PUB.BAC.FRUCHN.ANN.out | awk '$9>0 && $10<=1 {print $1,$2,$10}' \
    > PUB.BAC.FRUCHN.ANN.boxplot.out

# Prepare the fd boxplot
sed 's/,/\t/g' PUB.BAC.FRUCHN.ANN.out | awk '$9>0 && $10<=1 {print $10,"PUB.BAC.FRUCHN.ANN"}' \
    > PUB.BAC.FRUCHN.ANN.boxplot.out
sed 's/,/\t/g' PUB.BAB.FRU.ANN.out | awk '$9>0 && $10<=1 {print $10,"PUB.BAB.FRU.ANN"}' \
    > PUB.BAB.FRU.ANN.boxplot.out
sed 's/,/\t/g' PUB.BAB.CHN.ANN.out | awk '$9>0 && $10<=1 {print $10,"PUB.BAB.CHN.ANN"}' \
    > PUB.BAB.CHN.ANN.boxplot.out
sed 's/,/\t/g' PUB.BAP.FRU.ANN.out | awk '$9>0 && $10<=1 {print $10,"PUB.BAP.FRU.ANN"}' \
    > PUB.BAP.FRU.ANN.boxplot.out
sed 's/,/\t/g' PUB.BAP.CHN.ANN.out | awk '$9>0 && $10<=1 {print $10,"PUB.BAP.CHN.ANN"}' \
    > PUB.BAP.CHN.ANN.boxplot.out
cat *boxplot.out | sed '1i\fd Group' | awk '$1>0 {print $0}' | sed 's/ /\t/g' > BAC2FRUCHN.fd.boxplot.txt

# Prepare the fd line plot on Chr06
sed 's/,/\t/g' PUB.BAC.FRUCHN.GLA.out | \
    awk '$1=="Ca.S8.Chr06" && $9>0 && $10<=1 && $4>240000000 && $4<248000000 {print $1,$4,$10,"PUB.BAC.FRUCHN.GLA"}' \
    > PUB.BAC.FRUCHN.GLA.Chr06.region.out
sed 's/,/\t/g' PUB.BAB.FRU.ANN.out | \
    awk '$1=="Ca.S8.Chr06" && $9>0 && $10<=1 && $4>240000000 && $4<248000000 {print $1,$4,$10,"PUB.BAB.FRU.ANN"}' \
    > PUB.BAB.FRU.ANN.Chr06.region.out
sed 's/,/\t/g' PUB.BAB.CHN.ANN.out | \
    awk '$1=="Ca.S8.Chr06" && $9>0 && $10<=1 && $4>240000000 && $4<248000000 {print $1,$4,$10,"PUB.BAB.CHN.ANN"}' \
    > PUB.BAB.CHN.ANN.Chr06.region.out
sed 's/,/\t/g' PUB.BAP.FRU.ANN.out | \
    awk '$1=="Ca.S8.Chr06" && $9>0 && $10<=1 && $4>240000000 && $4<248000000 {print $1,$4,$10,"PUB.BAP.FRU.ANN"}' \
    > PUB.BAP.FRU.ANN.Chr06.region.out
sed 's/,/\t/g' PUB.BAP.CHN.ANN.out | \
    awk '$1=="Ca.S8.Chr06" && $9>0 && $10<=1 && $4>240000000 && $4<248000000 {print $1,$4,$10,"PUB.BAP.CHN.ANN"}' \
    > PUB.BAP.CHN.ANN.Chr06.region.out
cat *region.out | sed '1i\Chr Pos fd Group' | sed 's/ /\t/g' > fd.BAC2FRUCHN.Chr06.region.txt

# Calculate the average value of fd values.
cat PUB.BAC.FRUCHN.GLA.out | sed 's/,/ /g' | awk '$9>0 && $10<1{print $0}' | grep -v nan  | awk '{s+=$10} END {print s/NR}'
# Calculate the total value of fd
cat PUB.BAC.FRUCHN.GLA.out | sed 's/,/ /g' | awk '$9>0 && $10<1{print $0}' | grep -v nan  | awk '{s+=$6*$10} END {print s}' 
# Calculate the total G
cat PUB.BAC.FRUCHN.GLA.out | sed 's/,/ /g' | awk '$9>0 && $10<1{print $0}' | grep -v nan  | awk '{s+=$6} END {print s}' 

# For all the other comparisons
ABBABABAwindows.py -g Ca.S8.clean.geno.gz -o PUB.CHA.BAB.ANN.out -w 500000 -s 50000 -m 100 -f phased --minData 0.5 -T 20 \
    -P1 Ca.annuum -P2 Ca.bac.baccatum -P3 Ca.chacoense -O Ca.pubescens --popsFile Pop8File.txt --writeFailedWindows 
# C.bac.pendulum as P2
ABBABABAwindows.py -g Ca.S8.clean.geno.gz -o PUB.CHA.BAP.ANN.out -w 500000 -s 50000 -m 100 -f phased --minData 0.5 -T 20 \
    -P1 Ca.annuum -P2 Ca.bac.pendulum -P3 Ca.chacoense -O Ca.pubescens --popsFile Pop8File.txt --writeFailedWindows 
# Ca.chinense as P2    
ABBABABAwindows.py -g Ca.S8.clean.geno.gz -o PUB.CHA.CHN.ANN.out -w 500000 -s 50000 -m 100 -f phased --minData 0.5 -T 20 \
    -P1 Ca.annuum -P2 Ca.chinense -P3 Ca.chacoense -O Ca.pubescens --popsFile Pop8File.txt --writeFailedWindows
# Ca.frutescens as P2
ABBABABAwindows.py -g Ca.S8.clean.geno.gz -o PUB.CHA.FRU.ANN.out -w 500000 -s 50000 -m 100 -f phased --minData 0.5 -T 20 \
    -P1 Ca.annuum -P2 Ca.frutescens -P3 Ca.chacoense -O Ca.pubescens --popsFile Pop8File.txt --writeFailedWindows
# Ca.glabriusculum as P2
ABBABABAwindows.py -g Ca.S8.clean.geno.gz -o PUB.CHA.GLA.ANN.out -w 500000 -s 50000 -m 100 -f phased --minData 0.5 -T 20 \
    -P1 Ca.annuum -P2 Ca.glabriusculum -P3 Ca.chacoense -O Ca.pubescens --popsFile Pop8File.txt --writeFailedWindows

# Step 17.5.2 Prepare the file format for fd boxplot using CHA as P3
cat PUB.CHA.BAB.ANN.out | sed 's/,/ /g' | awk '$9>0 {print $0}' | grep -v nan > PUB.CHA.BAB.ANN.clean.txt
cat PUB.CHA.BAP.ANN.out | sed 's/,/ /g' | awk '$9>0 {print $0}' | grep -v nan > PUB.CHA.BAP.ANN.clean.txt
cat PUB.CHA.CHN.ANN.out | sed 's/,/ /g' | awk '$9>0 {print $0}' | grep -v nan > PUB.CHA.CHN.ANN.clean.txt
cat PUB.CHA.FRU.ANN.out | sed 's/,/ /g' | awk '$9>0 {print $0}' | grep -v nan > PUB.CHA.FRU.ANN.clean.txt
cat PUB.CHA.GLA.ANN.out | sed 's/,/ /g' | awk '$9>0 {print $0}' | grep -v nan > PUB.CHA.GLA.ANN.clean.txt
cat PUB.CHA.BAB.ANN.clean.txt | awk '{print $1,$4,$10,"BAB"}' > PUB.CHA.BAB.ANN.clean-2
cat PUB.CHA.BAP.ANN.clean.txt | awk '{print $1,$4,$10,"BAP"}' > PUB.CHA.BAP.ANN.clean-2
cat PUB.CHA.CHN.ANN.clean.txt | awk '{print $1,$4,$10,"CHN"}' > PUB.CHA.CHN.ANN.clean-2
cat PUB.CHA.FRU.ANN.clean.txt | awk '{print $1,$4,$10,"FRU"}' > PUB.CHA.FRU.ANN.clean-2
cat PUB.CHA.GLA.ANN.clean.txt | awk '{print $1,$4,$10,"GLA"}' > PUB.CHA.GLA.ANN.clean-2
# Prepare the bxoplot
cat PUB.CHA.*.ANN.clean-2 | grep -v fd | cut -f3,4 -d ' ' | sed '1i\fd Taxonomy' | sed 's/ /\t/g' > PUB.CHA.all.ANN.fd.txt

# Step 17.5.3 Prepare the genome-wide distribution plot
# Genome-wide distribution of fd values for Annuum clade
cat PUB.CHA.BAB.ANN.out | sed 's/,/ /g' | awk '{print $1,$4,$9,$10,"BAB"}' > PUB.CHA.BAB.ANN.clean-1
cat PUB.CHA.BAP.ANN.out | sed 's/,/ /g' | awk '{print $1,$4,$9,$10,"BAP"}' > PUB.CHA.BAP.ANN.clean-1
cat PUB.CHA.CHN.ANN.out | sed 's/,/ /g' | awk '{print $1,$4,$9,$10,"CHN"}' > PUB.CHA.CHN.ANN.clean-1
cat PUB.CHA.FRU.ANN.out | sed 's/,/ /g' | awk '{print $1,$4,$9,$10,"FRU"}' > PUB.CHA.FRU.ANN.clean-1
cat PUB.CHA.GLA.ANN.out | sed 's/,/ /g' | awk '{print $1,$4,$9,$10,"GLA"}' > PUB.CHA.GLA.ANN.clean-1

# Only chromosome 11
cat PUB.CHA.*.ANN.clean-1 | grep Chr11 | \
    sed 's/nan/0/g' | sed 's///g' | \
    sed '1i\scaffold mid D fd Taxonomy' > PUB.CHA.all.ANN.fd.Chr11.txt

# All chromosomes
cat PUB.CHA.*.ANN.clean-1 | \
    sed 's/nan/0/g' | sed 's///g' | \
    sed '1i\scaffold mid D fd Taxonomy' > PUB.CHA.all.ANN.fd.txt

rm PUB.CHA.*.ANN.clean-1

# Step 17.5.4 Calculate mean fd and PGI
# Calculate the average value of fd values.
for i in PUB.CHA.*.ANN.clean.txt; do cat $i | awk '{s+=$10} END {print s/NR}' ; done
# Calculate the total value of fd
for i in PUB.CHA.*.ANN.clean.txt; do cat $i | awk '{s+=$6*$10} END {print s}' ; done
# Calculate the total G
for i in PUB.CHA.*.ANN.clean.txt; do cat $i | awk '{s+=$6} END {print s}' ; done

# Step 17.6 Run ABBA-BABA using C.bac.baccatum as P3
ABBABABAwindows.py -g Ca.S8.clean.geno.gz -o PUB.BAB.BAP.ANN.out -w 500000 -s 50000 -m 100 -f phased --minData 0.5 -T 20 \
    -P1 Ca.annuum -P2 Ca.bac.pendulum -P3 Ca.bac.baccatum -O Ca.pubescens --popsFile Pop8File.txt --writeFailedWindows
ABBABABAwindows.py -g Ca.S8.clean.geno.gz -o PUB.BAB.FRU.ANN.out -w 500000 -s 50000 -m 100 -f phased --minData 0.5 -T 20 \
    -P1 Ca.annuum -P2 Ca.frutescens -P3 Ca.bac.baccatum -O Ca.pubescens --popsFile Pop8File.txt --writeFailedWindows
ABBABABAwindows.py -g Ca.S8.clean.geno.gz -o PUB.BAB.CHN.ANN.out -w 500000 -s 50000 -m 100 -f phased --minData 0.5 -T 20 \
    -P1 Ca.annuum -P2 Ca.chinense -P3 Ca.bac.baccatum -O Ca.pubescens --popsFile Pop8File.txt --writeFailedWindows
ABBABABAwindows.py -g Ca.S8.clean.geno.gz -o PUB.BAB.GLA.ANN.out -w 500000 -s 50000 -m 100 -f phased --minData 0.5 -T 20 \
    -P1 Ca.annuum -P2 Ca.glabriusculum -P3 Ca.bac.baccatum -O Ca.pubescens --popsFile Pop8File.txt --writeFailedWindows

# Run ABBA-BABA using C.bac.pendulum as P3
ABBABABAwindows.py -g Ca.S8.clean.geno.gz -o PUB.BAP.FRU.ANN.out -w 500000 -s 50000 -m 100 -f phased --minData 0.5 -T 20 \
    -P1 Ca.annuum -P2 Ca.frutescens -P3 Ca.bac.pendulum -O Ca.pubescens --popsFile Pop8File.txt --writeFailedWindows
ABBABABAwindows.py -g Ca.S8.clean.geno.gz -o PUB.BAP.CHN.ANN.out -w 500000 -s 50000 -m 100 -f phased --minData 0.5 -T 20 \
    -P1 Ca.annuum -P2 Ca.chinense -P3 Ca.bac.pendulum -O Ca.pubescens --popsFile Pop8File.txt --writeFailedWindows
ABBABABAwindows.py -g Ca.S8.clean.geno.gz -o PUB.BAP.GLA.ANN.out -w 500000 -s 50000 -m 100 -f phased --minData 0.5 -T 20 \
    -P1 Ca.annuum -P2 Ca.glabriusculum -P3 Ca.bac.pendulum -O Ca.pubescens --popsFile Pop8File.txt --writeFailedWindows

# Step 17.8 Run ABBA-BABA using Ca.frutescens as P3
ABBABABAwindows.py -g Ca.S8.clean.geno.gz -o BAP.FRU.CHN.ANN.out -w 500000 -s 50000 -m 100 -f phased --minData 0.5 -T 10 \
    -P1 Ca.annuum -P2 Ca.chinense -P3 Ca.frutescens -O Ca.bac.pendulum --popsFile Pop8File.txt --writeFailedWindows ;
ABBABABAwindows.py -g Ca.S8.clean.geno.gz -o BAP.FRU.GLA.ANN.out -w 500000 -s 50000 -m 100 -f phased --minData 0.5 -T 60 \
    -P1 Ca.annuum -P2 Ca.glabriusculum -P3 Ca.frutescens -O Ca.bac.pendulum --popsFile Pop8File.txt --writeFailedWindows

# Run ABBA-BABA using C. chinense as P3
ABBABABAwindows.py -g Ca.S8.clean.geno.gz -o FRU.CHN.GLA.ANN.out -w 500000 -s 50000 -m 100 -f phased --minData 0.5 -T 20 \
    -P1 Ca.annuum -P2 Ca.glabriusculum -P3 Ca.chinense -O Ca.frutescens --popsFile Pop8File.txt --writeFailedWindows 

# Step 17.9 Run ABBA-BABA using ANN as P2
ABBABABAwindows.py -g Ca.S8.clean.geno.gz -o PUB.CHA.ANN.GLA.out -w 500000 -s 50000 -m 100 -f phased --minData 0.5 -T 10 \
    -P1 Ca.glabriusculum -P2 Ca.annuum  -P3 Ca.chacoense -O Ca.pubescens --popsFile Pop8File.txt --writeFailedWindows
ABBABABAwindows.py -g Ca.S8.clean.geno.gz -o PUB.BAB.ANN.GLA.out -w 500000 -s 50000 -m 100 -f phased --minData 0.5 -T 10 \
    -P1 Ca.glabriusculum -P2 Ca.annuum  -P3 Ca.bac.baccatum -O Ca.pubescens --popsFile Pop8File.txt --writeFailedWindows
ABBABABAwindows.py -g Ca.S8.clean.geno.gz -o PUB.BAP.ANN.GLA.out -w 500000 -s 50000 -m 100 -f phased --minData 0.5 -T 10 \
    -P1 Ca.glabriusculum -P2 Ca.annuum  -P3 Ca.bac.pendulum -O Ca.pubescens --popsFile Pop8File.txt --writeFailedWindows
ABBABABAwindows.py -g Ca.S8.clean.geno.gz -o BAP.FRU.ANN.GLA.out -w 500000 -s 50000 -m 100 -f phased --minData 0.5 -T 10 \
    -P1 Ca.glabriusculum -P2 Ca.annuum  -P3 Ca.frutescens -O Ca.bac.pendulum --popsFile Pop8File.txt --writeFailedWindows
ABBABABAwindows.py -g Ca.S8.clean.geno.gz -o BAP.CHN.ANN.GLA.out -w 500000 -s 50000 -m 100 -f phased --minData 0.5 -T 10 \
    -P1 Ca.glabriusculum -P2 Ca.annuum  -P3 Ca.chinense -O Ca.bac.pendulum --popsFile Pop8File.txt --writeFailedWindows

# Combine every two columns
head -1  Chr06.tped | sed 's! \([^ ]\+\)\( \|$\)!\1 !g'

```

## SweeD

### Installation
```bash
# Installation
cd /data/zhaojiantao/tools
git clone https://github.com/alachins/sweed.git
cd sweed
# Install
make -f Makefile.gcc
make -f Makefile.PTHREADS.gcc

# Install DMTCP library
DMTCPAWARELIB=/data/zhaojiantao/tools/sweed/dmtcp/lib/dmtcp

```

### BAC with 10k grid

```bash
# SweeD does not accept zipped vcf
# I will first extract the BAC accessions and filter the SNPs for a faster running
cd /data/zhaojiantao/pepper/SweeD/BAC
# Convert vcf to transposed ped
for i in $(seq -w 1 12); do
    plink --vcf ../../vcf.clean/Chr$i.clean.vcf.gz \
        --keep BAC.list \
        --double-id \
        --chr Chr$i \
        --maf 0.05 \
        --geno 0.1 \
        --recode vcf \
        --out Chr$i & 
done

# Prepare the bed file containing the chr and grid numbers
grep -v Chr00 ../../References/Zhangshugang/Zhangshugang.fa.fai | \
    awk '{print $1,int($2/10000+0.5)}' \
    > Chr.grids

# Run SweeD
IFS=$'\n'
for LINE in $(cat Chr.grids); do
    Chr=$(echo ${LINE} | awk '{print $1}')
    grid=$(echo ${LINE} | awk '{print $2}')
    SweeD-P -name $Chr -input $Chr.vcf -grid $grid -threads 10 &
done

# Add the chromosome ids to the output
sed -n '4,$p' SweeD_Report.Chr01 | awk '{print "1", $0}' > SweeD_Report.Chr01-1
sed -n '4,$p' SweeD_Report.Chr02 | awk '{print "2", $0}' > SweeD_Report.Chr02-1
sed -n '4,$p' SweeD_Report.Chr03 | awk '{print "3", $0}' > SweeD_Report.Chr03-1
sed -n '4,$p' SweeD_Report.Chr04 | awk '{print "4", $0}' > SweeD_Report.Chr04-1
sed -n '4,$p' SweeD_Report.Chr05 | awk '{print "5", $0}' > SweeD_Report.Chr05-1
sed -n '4,$p' SweeD_Report.Chr06 | awk '{print "6", $0}' > SweeD_Report.Chr06-1
sed -n '4,$p' SweeD_Report.Chr07 | awk '{print "7", $0}' > SweeD_Report.Chr07-1
sed -n '4,$p' SweeD_Report.Chr08 | awk '{print "8", $0}' > SweeD_Report.Chr08-1
sed -n '4,$p' SweeD_Report.Chr09 | awk '{print "9", $0}' > SweeD_Report.Chr09-1
sed -n '4,$p' SweeD_Report.Chr10 | awk '{print "10", $0}' > SweeD_Report.Chr10-1
sed -n '4,$p' SweeD_Report.Chr11 | awk '{print "11", $0}' > SweeD_Report.Chr11-1
sed -n '4,$p' SweeD_Report.Chr12 | awk '{print "12", $0}' > SweeD_Report.Chr12-1

# Combine the results
cat SweeD_Report.Chr*-1 | awk '{print $1,"Chr"$1"_"$2,$2,$3,$4,$5,$6}' | \
    sed '1i\Chr SNP Position Likelihoood Alpha StartPos EndPos' | sed 's/ /\t/g' \
    > SweeD_Report.BAC.10k.txt

```

### BAC with pruned SNPs

```bash
# Working directory
cd /data/zhaojiantao/pepper/SweeD/BAC.prune

# Simplify the sample id and add SNP info
for i in $(seq -w 1 12); do
    # Simplify the sample id
    grep "#" ../BAC.10k/Chr$i.vcf | \
        sed 's/Grif_9199_Grif_9199/Grif_9199/g' | sed 's/Grif_9201_Grif_9201/Grif_9201/g' | \
        sed 's/Grif_9352_Grif_9352/Grif_9352/g' | sed 's/Grif_9354_Grif_9354/Grif_9354/g' | \
        sed 's/Grif_9355_Grif_9355/Grif_9355/g' | sed 's/PI_159235_PI_159235/PI_159235/g' | \
        sed 's/PI_159242_PI_159242/PI_159242/g' | sed 's/PI_159249_PI_159249/PI_159249/g' | \
        sed 's/PI_159260_PI_159260/PI_159260/g' | sed 's/PI_188481_PI_188481/PI_188481/g' | \
        sed 's/PI_199506_PI_199506/PI_199506/g' | sed 's/PI_215699_PI_215699/PI_215699/g' | \
        sed 's/PI_215741_PI_215741/PI_215741/g' | sed 's/PI_241679_PI_241679/PI_241679/g' | \
        sed 's/PI_257110_PI_257110/PI_257110/g' | sed 's/PI_257141_PI_257141/PI_257141/g' | \
        sed 's/PI_257157_PI_257157/PI_257157/g' | sed 's/PI_257161_PI_257161/PI_257161/g' | \
        sed 's/PI_257173_PI_257173/PI_257173/g' | sed 's/PI_257174_PI_257174/PI_257174/g' | \
        sed 's/PI_260488_PI_260488/PI_260488/g' | sed 's/PI_260535_PI_260535/PI_260535/g' | \
        sed 's/PI_260538_PI_260538/PI_260538/g' | sed 's/PI_260545_PI_260545/PI_260545/g' | \
        sed 's/PI_260549_PI_260549/PI_260549/g' | sed 's/PI_260550_PI_260550/PI_260550/g' | \
        sed 's/PI_260551_PI_260551/PI_260551/g' | sed 's/PI_260564_PI_260564/PI_260564/g' | \
        sed 's/PI_260567_PI_260567/PI_260567/g' | sed 's/PI_260575_PI_260575/PI_260575/g' | \
        sed 's/PI_260582_PI_260582/PI_260582/g' | sed 's/PI_260590_PI_260590/PI_260590/g' | \
        sed 's/PI_260595_PI_260595/PI_260595/g' | sed 's/PI_266041_PI_266041/PI_266041/g' | \
        sed 's/PI_267729_PI_267729/PI_267729/g' | sed 's/PI_281308_PI_281308/PI_281308/g' | \
        sed 's/PI_281309_PI_281309/PI_281309/g' | sed 's/PI_281310_PI_281310/PI_281310/g' | \
        sed 's/PI_281398_PI_281398/PI_281398/g' | sed 's/PI_281407_PI_281407/PI_281407/g' | \
        sed 's/PI_281436_PI_281436/PI_281436/g' | sed 's/PI_290982_PI_290982/PI_290982/g' | \
        sed 's/PI_290983_PI_290983/PI_290983/g' | sed 's/PI_293349_PI_293349/PI_293349/g' | \
        sed 's/PI_315020_PI_315020/PI_315020/g' | sed 's/PI_337524_PI_337524/PI_337524/g' | \
        sed 's/PI_355813_PI_355813/PI_355813/g' | sed 's/PI_370010_PI_370010/PI_370010/g' | \
        sed 's/PI_424732_PI_424732/PI_424732/g' | sed 's/PI_439359_PI_439359/PI_439359/g' | \
        sed 's/PI_439360_PI_439360/PI_439360/g' | sed 's/PI_439363_PI_439363/PI_439363/g' | \
        sed 's/PI_439366_PI_439366/PI_439366/g' | sed 's/PI_439367_PI_439367/PI_439367/g' | \
        sed 's/PI_439369_PI_439369/PI_439369/g' | sed 's/PI_439371_PI_439371/PI_439371/g' | \
        sed 's/PI_439380_PI_439380/PI_439380/g' | sed 's/PI_439384_PI_439384/PI_439384/g' | \
        sed 's/PI_439386_PI_439386/PI_439386/g' | sed 's/PI_439390_PI_439390/PI_439390/g' | \
        sed 's/PI_439395_PI_439395/PI_439395/g' | sed 's/PI_439398_PI_439398/PI_439398/g' | \
        sed 's/PI_439401_PI_439401/PI_439401/g' | sed 's/PI_439402_PI_439402/PI_439402/g' | \
        sed 's/PI_439403_PI_439403/PI_439403/g' | sed 's/PI_439405_PI_439405/PI_439405/g' | \
        sed 's/PI_439407_PI_439407/PI_439407/g' | sed 's/PI_439410_PI_439410/PI_439410/g' | \
        sed 's/PI_439411_PI_439411/PI_439411/g' | sed 's/PI_439472_PI_439472/PI_439472/g' | \
        sed 's/PI_439528_PI_439528/PI_439528/g' | sed 's/PI_441516_PI_441516/PI_441516/g' | \
        sed 's/PI_441523_PI_441523/PI_441523/g' | sed 's/PI_441531_PI_441531/PI_441531/g' | \
        sed 's/PI_441532_PI_441532/PI_441532/g' | sed 's/PI_441536_PI_441536/PI_441536/g' | \
        sed 's/PI_441537_PI_441537/PI_441537/g' | sed 's/PI_441541_PI_441541/PI_441541/g' | \
        sed 's/PI_441547_PI_441547/PI_441547/g' | sed 's/PI_441548_PI_441548/PI_441548/g' | \
        sed 's/PI_441558_PI_441558/PI_441558/g' | sed 's/PI_441578_PI_441578/PI_441578/g' | \
        sed 's/PI_441590_PI_441590/PI_441590/g' | sed 's/PI_441592_PI_441592/PI_441592/g' | \
        sed 's/PI_441594_PI_441594/PI_441594/g' | sed 's/PI_441597_PI_441597/PI_441597/g' | \
        sed 's/PI_441656_PI_441656/PI_441656/g' | sed 's/PI_441674_PI_441674/PI_441674/g' | \
        sed 's/PI_441685_PI_441685/PI_441685/g' | sed 's/PI_441699_PI_441699/PI_441699/g' | \
        sed 's/PI_446900_PI_446900/PI_446900/g' | sed 's/PI_446909_PI_446909/PI_446909/g' | \
        sed 's/PI_497972_PI_497972/PI_497972/g' | sed 's/PI_555611_PI_555611/PI_555611/g' | \
        sed 's/PI_585247_PI_585247/PI_585247/g' | sed 's/PI_585250_PI_585250/PI_585250/g' | \
        sed 's/PI_593606_PI_593606/PI_593606/g' | sed 's/PI_594137_PI_594137/PI_594137/g' | \
        sed 's/PI_594138_PI_594138/PI_594138/g' | sed 's/PI_595905_PI_595905/PI_595905/g' | \
        sed 's/PI_596059_PI_596059/PI_596059/g' | sed 's/PI_631146_PI_631146/PI_631146/g' | \
        sed 's/PI_632923_PI_632923/PI_632923/g' | sed 's/PI_632926_PI_632926/PI_632926/g' | \
        sed 's/PI_632928_PI_632928/PI_632928/g' | sed 's/PI_633751_PI_633751/PI_633751/g' | \
        sed 's/PI_633755_PI_633755/PI_633755/g' | sed 's/PI_633758_PI_633758/PI_633758/g' | \
        sed 's/PI_639129_PI_639129/PI_639129/g' | sed 's/PI_639648_PI_639648/PI_639648/g' | \
        sed 's/PI_640882_PI_640882/PI_640882/g' | sed 's/PI_640884_PI_640884/PI_640884/g' | \
        sed 's/PI_640889_PI_640889/PI_640889/g' | sed 's/PI_643124_PI_643124/PI_643124/g' | \
        sed 's/PI_659105_PI_659105/PI_659105/g' | sed 's/PI_660972_PI_660972/PI_660972/g' | \
        sed 's/PI_688675_PI_688675/PI_688675/g' \
        > Chr$i.header
    # Add SNP info
    grep -v "#" ../BAC.10k/Chr$i.vcf | \
        awk '{print $1,$2,$1"_"$2}' | sed 's/ /\t/g' | \
        paste - <(grep -v "#" ../BAC.10k/Chr$i.vcf) | \
        cut -f1-3,7- | \
        cat Chr$i.header - \
        > Chr$i.vcf
    # Remove temporary file
    rm Chr$i.header
done

# SNP pruning
for i in $(seq -w 1 12); do
    plink --vcf Chr$i.vcf \
    --indep-pairwise 50 10 0.2 \
    --allow-extra-chr --double-id \
    --out Chr$i
    # Extract the prunned SNPs
    plink --vcf Chr$i.vcf \
        --extract Chr$i.prune.in \
        --recode vcf --double-id --allow-extra-chr \
        --out Chr$i.prune
    # Remove temporary files
    rm *log *nosex *out *in
done

# Prepare the bed file containing the chr and grid numbers
grep -v Chr00 ../../References/Zhangshugang/Zhangshugang.fa.fai | \
    awk '{print $1,int($2/10000+0.5)}' \
    > Chr.grids

# Run SweeD
IFS=$'\n'
for LINE in $(cat Chr.grids); do
    Chr=$(echo ${LINE} | awk '{print $1}')
    grid=$(echo ${LINE} | awk '{print $2}')
    SweeD-P -name $Chr -input $Chr.prune.vcf -grid $grid -threads 10 &
done

# Add the chromosome ids to the output
sed -n '4,$p' SweeD_Report.Chr01 | awk '{print "1", $0}' > SweeD_Report.Chr01-1
sed -n '4,$p' SweeD_Report.Chr02 | awk '{print "2", $0}' > SweeD_Report.Chr02-1
sed -n '4,$p' SweeD_Report.Chr03 | awk '{print "3", $0}' > SweeD_Report.Chr03-1
sed -n '4,$p' SweeD_Report.Chr04 | awk '{print "4", $0}' > SweeD_Report.Chr04-1
sed -n '4,$p' SweeD_Report.Chr05 | awk '{print "5", $0}' > SweeD_Report.Chr05-1
sed -n '4,$p' SweeD_Report.Chr06 | awk '{print "6", $0}' > SweeD_Report.Chr06-1
sed -n '4,$p' SweeD_Report.Chr07 | awk '{print "7", $0}' > SweeD_Report.Chr07-1
sed -n '4,$p' SweeD_Report.Chr08 | awk '{print "8", $0}' > SweeD_Report.Chr08-1
sed -n '4,$p' SweeD_Report.Chr09 | awk '{print "9", $0}' > SweeD_Report.Chr09-1
sed -n '4,$p' SweeD_Report.Chr10 | awk '{print "10", $0}' > SweeD_Report.Chr10-1
sed -n '4,$p' SweeD_Report.Chr11 | awk '{print "11", $0}' > SweeD_Report.Chr11-1
sed -n '4,$p' SweeD_Report.Chr12 | awk '{print "12", $0}' > SweeD_Report.Chr12-1

# Combine the results
cat SweeD_Report.Chr*-1 | awk '{print $1,"Chr"$1"_"$2,$2,$3,$4,$5,$6}' | \
    sed '1i\Chr SNP Position Likelihoood Alpha StartPos EndPos' | sed 's/ /\t/g' \
    > SweeD_Report.BAC.prune.10k.txt

```


### GLA with 10k grid

```bash
# For GLA
cd /data/zhaojiantao/pepper/SweeD/GLA.10k

# Prepare the sample list
cat /data/zhaojiantao/pepper/Pop.list.pure/Ca.glabriusculum_I.list \
    /data/zhaojiantao/pepper/Pop.list.pure/Ca.annuum.list | \
    awk '{print $1,$1}' | sed 's/ /\t/g' > GLA.list

# Convert vcf to transposed ped
for i in $(seq -w 1 12); do
    plink --vcf ../../vcf.clean/Chr$i.clean.vcf.gz \
        --keep GLA.list \
        --double-id \
        --chr Chr$i \
        --maf 0.05 \
        --geno 0.1 \
        --recode vcf \
        --out Chr$i & 
done

# Prepare the bed file containing the chr and grid numbers
grep -v Chr00 ../../References/Zhangshugang/Zhangshugang.fa.fai | \
    awk '{print $1,int($2/10000+0.5)}' \
    > Chr.grids

# Run SweeD
IFS=$'\n'
for LINE in $(cat Chr.grids); do
    Chr=$(echo ${LINE} | awk '{print $1}')
    grid=$(echo ${LINE} | awk '{print $2}')
    SweeD-P -name $Chr -input $Chr.vcf -grid $grid -threads 10 &
done

# Add the chromosome ids to the output
sed -n '4,$p' SweeD_Report.Chr01 | awk '{print "1", $0}' > SweeD_Report.Chr01-1
sed -n '4,$p' SweeD_Report.Chr02 | awk '{print "2", $0}' > SweeD_Report.Chr02-1
sed -n '4,$p' SweeD_Report.Chr03 | awk '{print "3", $0}' > SweeD_Report.Chr03-1
sed -n '4,$p' SweeD_Report.Chr04 | awk '{print "4", $0}' > SweeD_Report.Chr04-1
sed -n '4,$p' SweeD_Report.Chr05 | awk '{print "5", $0}' > SweeD_Report.Chr05-1
sed -n '4,$p' SweeD_Report.Chr06 | awk '{print "6", $0}' > SweeD_Report.Chr06-1
sed -n '4,$p' SweeD_Report.Chr07 | awk '{print "7", $0}' > SweeD_Report.Chr07-1
sed -n '4,$p' SweeD_Report.Chr08 | awk '{print "8", $0}' > SweeD_Report.Chr08-1
sed -n '4,$p' SweeD_Report.Chr09 | awk '{print "9", $0}' > SweeD_Report.Chr09-1
sed -n '4,$p' SweeD_Report.Chr10 | awk '{print "10", $0}' > SweeD_Report.Chr10-1
sed -n '4,$p' SweeD_Report.Chr11 | awk '{print "11", $0}' > SweeD_Report.Chr11-1
sed -n '4,$p' SweeD_Report.Chr12 | awk '{print "12", $0}' > SweeD_Report.Chr12-1

# Combine the results
cat SweeD_Report.Chr*-1 | awk '{print $1,"Chr"$1"_"$2,$2,$3,$4,$5,$6}' | \
    sed '1i\Chr SNP Position Likelihoood Alpha StartPos EndPos' | sed 's/ /\t/g' \
    > SweeD_Report.GLA.10k.txt

```




## XP-CLR for C. annuum

```bash
# Working directory
cd /data/zhaojiantao/pepper/XPCLR

# Blast the markers to S8 to update the physical position
blastall -p blastn -i cM.BP.fa -d ../pepper_genome/fa -m 8 -e 1e-40 -o cM.BP.S8.blast

# Prepare the SNP format for inferring the pseudo genetic position
zcat ../vcf.clean/Chr01.clean.vcf.gz | grep -v "#" | awk '{print "Chr01_"$2,$2}' > Chr01.BP.txt
zcat ../vcf.clean/Chr02.clean.vcf.gz | grep -v "#" | awk '{print "Chr02_"$2,$2}' > Chr02.BP.txt
zcat ../vcf.clean/Chr03.clean.vcf.gz | grep -v "#" | awk '{print "Chr03_"$2,$2}' > Chr03.BP.txt
zcat ../vcf.clean/Chr04.clean.vcf.gz | grep -v "#" | awk '{print "Chr04_"$2,$2}' > Chr04.BP.txt
zcat ../vcf.clean/Chr05.clean.vcf.gz | grep -v "#" | awk '{print "Chr05_"$2,$2}' > Chr05.BP.txt
zcat ../vcf.clean/Chr06.clean.vcf.gz | grep -v "#" | awk '{print "Chr06_"$2,$2}' > Chr06.BP.txt
zcat ../vcf.clean/Chr07.clean.vcf.gz | grep -v "#" | awk '{print "Chr07_"$2,$2}' > Chr07.BP.txt
zcat ../vcf.clean/Chr08.clean.vcf.gz | grep -v "#" | awk '{print "Chr08_"$2,$2}' > Chr08.BP.txt
zcat ../vcf.clean/Chr09.clean.vcf.gz | grep -v "#" | awk '{print "Chr09_"$2,$2}' > Chr09.BP.txt
zcat ../vcf.clean/Chr10.clean.vcf.gz | grep -v "#" | awk '{print "Chr10_"$2,$2}' > Chr10.BP.txt
zcat ../vcf.clean/Chr11.clean.vcf.gz | grep -v "#" | awk '{print "Chr11_"$2,$2}' > Chr11.BP.txt
zcat ../vcf.clean/Chr12.clean.vcf.gz | grep -v "#" | awk '{print "Chr12_"$2,$2}' > Chr12.BP.txt

# Prepare the real genetic maps
for i in $(seq -w 1 12); do 
    cat Linkage.clean.txt | grep Chr$i |  awk '{print $2,$5,$4}' > Chr$i.map.txt ; 
done

# Infer the genetic position 
for i in $(seq -w 1 12); do 
    perl assign_genetic.pl Chr$i.map.txt Chr$i.BP.txt
done

# Compare the relationship between the inferred and real genetic maps
for i in $(seq -w 1 12); do 
    cat Chr$i.map.txt | awk '{print $0,"Real"}' > Chr$i.map-1.txt
    cat Chr$i.snp.genetic.txt | awk '{print $1,$2,$3,"Inferred"}' | \
        cat - Chr$i.map-1.txt | sed '1i\SNP Pos cM Group' | sed 's/ /\t/g' \
        > Chr$i.compare.txt
    rm Chr$i.map-1.txt
done

# SNP position files in LinkageMap folder
for i in $(seq -w 1 12); do 
    zcat ../vcf.clean/Chr$i.clean.vcf.gz | grep -v "#" | awk '{print $1,$2,$4,$5}' > Chr$i.BP-2.txt
    paste Chr$i.snp.genetic.txt Chr$i.BP-2.txt | \
        awk '{print $1,$4,$3,$2,$6,$7}' | sed 's/ /\t/g' > ../XPCLR/Chr$i.snp
    rm Chr$i.BP-2.txt
done

# Prepare the genotypes for xpclr in XPCLR folder
# Make ped file using this chromosome map
for j in Ca.annuum Ca.bac.baccatum Ca.bac.pendulum Ca.chacoense Ca.chinense Ca.frutescens Ca.glabriusculum ; do 
    cat ../Pop.list/$j.list | awk '{print $1 "\t" $1}' > $j.list.txt
done

for i in $(seq -w 1 12); do
    for j in Ca.annuum Ca.bac.baccatum Ca.bac.pendulum Ca.chacoense Ca.chinense Ca.frutescens Ca.glabriusculum ; do 
        # Convert vcf to transposed ped
        plink --vcf ../vcf.clean/Chr$i.clean.vcf.gz \
            --keep $j.list.txt  \
            --allow-extra-chr --chr Chr$i \
            --recode 01 transpose \
            -output-missing-genotype 9 --double-id \
            --set-missing-var-ids @:# --keep-allele-order \
            --out $j.Chr$i
        # Prepare the transposed ped file to xpclr geno file format
        cut -d " " -f 5- $j.Chr$i.tped  | awk '{print $0" "}' > $j.Chr$i.geno
    done
done

# Remove unnecessary files
rm *log *nosex *tfam *tped 

# XPCLR of Ca.glabriusculum vs Ca.annuum
# Chr01
# Prepare the input files
for i in $(seq -w 500000 500000 9086854); do
    cat Ca.annuum.Chr01.geno | sed -n "1,499999p" > GLA.ANN/Ca.annuum.Chr01.0000000.geno ; 
    cat Ca.glabriusculum.Chr01.geno | sed -n "1,499999p" > GLA.ANN/Ca.glabriusculum.Chr01.0000000.geno ; 
    cat Chr01.snp | sed -n "1,499999p" > GLA.ANN/Chr01.0000000.snp
    cat Ca.annuum.Chr01.geno | sed -n "1p;${i},+499999p" > GLA.ANN/Ca.annuum.Chr01.$i.geno ; 
    cat Ca.glabriusculum.Chr01.geno | sed -n "1p;${i},+499999p" > GLA.ANN/Ca.glabriusculum.Chr01.$i.geno ; 
    cat Chr01.snp | sed -n "1p;${i},+499999p" > GLA.ANN/Chr01.$i.snp ; 
done
# Run XPCLR
for i in $(seq -w 000000 500000 9086854); do
    XPCLR -xpclr GLA.ANN/Ca.annuum.Chr01.$i.geno \
        GLA.ANN/Ca.glabriusculum.Chr01.$i.geno \
        GLA.ANN/Chr01.$i.snp GLA.ANN/GLA.ANN.Chr01.$i.out \
        -w1 0.0005 100 100 Chr01 -p0 0.7 > GLA.ANN/GLA.ANN.Chr01.$i.out.xpclr.log 
done

# Chr02
# Prepare the input files
for i in $(seq -w 500000 500000 4776664); do
    cat Ca.annuum.Chr02.geno | sed -n "1,499999p" > GLA.ANN/Ca.annuum.Chr02.0000000.geno ; 
    cat Ca.glabriusculum.Chr02.geno | sed -n "1,499999p" > GLA.ANN/Ca.glabriusculum.Chr02.0000000.geno ; 
    cat Chr02.snp | sed -n "1,499999p" > GLA.ANN/Chr02.0000000.snp
    cat Ca.annuum.Chr02.geno | sed -n "1p;${i},+499999p" > GLA.ANN/Ca.annuum.Chr02.$i.geno ; 
    cat Ca.glabriusculum.Chr02.geno | sed -n "1p;${i},+499999p" > GLA.ANN/Ca.glabriusculum.Chr02.$i.geno ; 
    cat Chr02.snp | sed -n "1p;${i},+499999p" > GLA.ANN/Chr02.$i.snp ; 
done
# Run XPCLR
    for i in $(seq -w 000000 500000 4776664); do
    XPCLR -xpclr GLA.ANN/Ca.annuum.Chr02.$i.geno \
        GLA.ANN/Ca.glabriusculum.Chr02.$i.geno \
        GLA.ANN/Chr02.$i.snp GLA.ANN/GLA.ANN.Chr02.$i.out \
        -w1 0.0005 100 100 Chr02 -p0 0.7 > GLA.ANN/GLA.ANN.Chr02.$i.out.xpclr.log 
done

# Chr03
# Prepare the input files
for i in $(seq -w 500000 500000 8369414); do
    cat Ca.annuum.Chr03.geno | sed -n "1,499999p" > GLA.ANN/Ca.annuum.Chr03.0000000.geno ; 
    cat Ca.glabriusculum.Chr03.geno | sed -n "1,499999p" > GLA.ANN/Ca.glabriusculum.Chr03.0000000.geno ; 
    cat Chr03.snp | sed -n "1,499999p" > GLA.ANN/Chr03.0000000.snp
    cat Ca.annuum.Chr03.geno | sed -n "1p;${i},+499999p" > GLA.ANN/Ca.annuum.Chr03.$i.geno ; 
    cat Ca.glabriusculum.Chr03.geno | sed -n "1p;${i},+499999p" > GLA.ANN/Ca.glabriusculum.Chr03.$i.geno ; 
    cat Chr03.snp | sed -n "1p;${i},+499999p" > GLA.ANN/Chr03.$i.snp ; 
done
# Run XPCLR
for i in $(seq -w 000000 500000 8369414); do
XPCLR -xpclr GLA.ANN/Ca.annuum.Chr03.$i.geno \
    GLA.ANN/Ca.glabriusculum.Chr03.$i.geno \
    GLA.ANN/Chr03.$i.snp GLA.ANN/GLA.ANN.Chr03.$i.out \
    -w1 0.0005 100 100 Chr03 -p0 0.7 > GLA.ANN/GLA.ANN.Chr03.$i.out.xpclr.log 
done

# Chr04
# Prepare the input files
for i in $(seq -w 500000 500000 7616000); do
    cat Ca.annuum.Chr04.geno | sed -n "1,499999p" > GLA.ANN/Ca.annuum.Chr04.0000000.geno ; 
    cat Ca.glabriusculum.Chr04.geno | sed -n "1,499999p" > GLA.ANN/Ca.glabriusculum.Chr04.0000000.geno ; 
    cat Chr04.snp | sed -n "1,499999p" > GLA.ANN/Chr04.0000000.snp
    cat Ca.annuum.Chr04.geno | sed -n "1p;${i},+499999p" > GLA.ANN/Ca.annuum.Chr04.$i.geno ; 
    cat Ca.glabriusculum.Chr04.geno | sed -n "1p;${i},+499999p" > GLA.ANN/Ca.glabriusculum.Chr04.$i.geno ; 
    cat Chr04.snp | sed -n "1p;${i},+499999p" > GLA.ANN/Chr04.$i.snp ; 
done
# Run XPCLR
for i in $(seq -w 000000 500000 7616000); do
XPCLR -xpclr GLA.ANN/Ca.annuum.Chr04.$i.geno \
    GLA.ANN/Ca.glabriusculum.Chr04.$i.geno \
    GLA.ANN/Chr04.$i.snp GLA.ANN/GLA.ANN.Chr04.$i.out \
    -w1 0.0005 100 100 Chr04 -p0 0.7 > GLA.ANN/GLA.ANN.Chr04.$i.out.xpclr.log 
done 

# Chr05
# Prepare the input files
for i in $(seq -w 500000 500000 7570910); do
    cat Ca.annuum.Chr05.geno | sed -n "1,499999p" > GLA.ANN/Ca.annuum.Chr05.0000000.geno ; 
    cat Ca.glabriusculum.Chr05.geno | sed -n "1,499999p" > GLA.ANN/Ca.glabriusculum.Chr05.0000000.geno ; 
    cat Chr05.snp | sed -n "1,499999p" > GLA.ANN/Chr05.0000000.snp
    cat Ca.annuum.Chr05.geno | sed -n "1p;${i},+499999p" > GLA.ANN/Ca.annuum.Chr05.$i.geno ; 
    cat Ca.glabriusculum.Chr05.geno | sed -n "1p;${i},+499999p" > GLA.ANN/Ca.glabriusculum.Chr05.$i.geno ; 
    cat Chr05.snp | sed -n "1p;${i},+499999p" > GLA.ANN/Chr05.$i.snp ; 
done
# Run XPCLR
for i in $(seq -w 000000 500000 7570910); do
XPCLR -xpclr GLA.ANN/Ca.annuum.Chr05.$i.geno \
    GLA.ANN/Ca.glabriusculum.Chr05.$i.geno \
    GLA.ANN/Chr05.$i.snp GLA.ANN/GLA.ANN.Chr05.$i.out \
    -w1 0.0005 100 100 Chr05 -p0 0.7 > GLA.ANN/GLA.ANN.Chr05.$i.out.xpclr.log 
done

# Chr06
# Prepare the input files
for i in $(seq -w 500000 500000 7522020); do
    cat Ca.annuum.Chr06.geno | sed -n "1,499999p" > GLA.ANN/Ca.annuum.Chr06.0000000.geno ; 
    cat Ca.glabriusculum.Chr06.geno | sed -n "1,499999p" > GLA.ANN/Ca.glabriusculum.Chr06.0000000.geno ; 
    cat Chr06.snp | sed -n "1,499999p" > GLA.ANN/Chr06.0000000.snp
    cat Ca.annuum.Chr06.geno | sed -n "1p;${i},+499999p" > GLA.ANN/Ca.annuum.Chr06.$i.geno ; 
    cat Ca.glabriusculum.Chr06.geno | sed -n "1p;${i},+499999p" > GLA.ANN/Ca.glabriusculum.Chr06.$i.geno ; 
    cat Chr06.snp | sed -n "1p;${i},+499999p" > GLA.ANN/Chr06.$i.snp ; 
done
# Run XPCLR
for i in $(seq -w 000000 500000 7522020); do
    XPCLR -xpclr GLA.ANN/Ca.annuum.Chr06.$i.geno \
        GLA.ANN/Ca.glabriusculum.Chr06.$i.geno \
        GLA.ANN/Chr06.$i.snp GLA.ANN/GLA.ANN.Chr06.$i.out \
        -w1 0.0005 100 100 Chr06 -p0 0.7 > GLA.ANN/GLA.ANN.Chr06.$i.out.xpclr.log 
done

# Chr07
# Prepare the input files
for i in $(seq -w 500000 500000 7824408); do
    cat Ca.annuum.Chr07.geno | sed -n "1,499999p" > GLA.ANN/Ca.annuum.Chr07.0000000.geno ; 
    cat Ca.glabriusculum.Chr07.geno | sed -n "1,499999p" > GLA.ANN/Ca.glabriusculum.Chr07.0000000.geno ; 
    cat Chr07.snp | sed -n "1,499999p" > GLA.ANN/Chr07.0000000.snp
    cat Ca.annuum.Chr07.geno | sed -n "1p;${i},+499999p" > GLA.ANN/Ca.annuum.Chr07.$i.geno ; 
    cat Ca.glabriusculum.Chr07.geno | sed -n "1p;${i},+499999p" > GLA.ANN/Ca.glabriusculum.Chr07.$i.geno ; 
    cat Chr07.snp | sed -n "1p;${i},+499999p" > GLA.ANN/Chr07.$i.snp ; 
done
# Run XPCLR
for i in $(seq -w 000000 500000 7824408); do
    XPCLR -xpclr GLA.ANN/Ca.annuum.Chr07.$i.geno \
        GLA.ANN/Ca.glabriusculum.Chr07.$i.geno \
        GLA.ANN/Chr07.$i.snp GLA.ANN/GLA.ANN.Chr07.$i.out \
        -w1 0.0005 100 100 Chr07 -p0 0.7 > GLA.ANN/GLA.ANN.Chr07.$i.out.xpclr.log 
done

# Chr08
# Prepare the input files
for i in $(seq -w 500000 500000 5099543); do
    cat Ca.annuum.Chr08.geno | sed -n "1,499999p" > GLA.ANN/Ca.annuum.Chr08.0000000.geno ; 
    cat Ca.glabriusculum.Chr08.geno | sed -n "1,499999p" > GLA.ANN/Ca.glabriusculum.Chr08.0000000.geno ; 
    cat Chr08.snp | sed -n "1,499999p" > GLA.ANN/Chr08.0000000.snp
    cat Ca.annuum.Chr08.geno | sed -n "1p;${i},+499999p" > GLA.ANN/Ca.annuum.Chr08.$i.geno ; 
    cat Ca.glabriusculum.Chr08.geno | sed -n "1p;${i},+499999p" > GLA.ANN/Ca.glabriusculum.Chr08.$i.geno ; 
    cat Chr08.snp | sed -n "1p;${i},+499999p" > GLA.ANN/Chr08.$i.snp ; 
done
# Run XPCLR
for i in $(seq -w 000000 500000 5099543); do
    XPCLR -xpclr GLA.ANN/Ca.annuum.Chr08.$i.geno \
        GLA.ANN/Ca.glabriusculum.Chr08.$i.geno \
        GLA.ANN/Chr08.$i.snp GLA.ANN/GLA.ANN.Chr08.$i.out \
        -w1 0.0005 100 1000 Chr08 -p0 0.7 > GLA.ANN/GLA.ANN.Chr08.$i.out.xpclr.log 
done

# Chr09
# Prepare the input files
for i in $(seq -w 500000 500000 8390222); do
    cat Ca.annuum.Chr09.geno | sed -n "1,499999p" > GLA.ANN/Ca.annuum.Chr09.0000000.geno ; 
    cat Ca.glabriusculum.Chr09.geno | sed -n "1,499999p" > GLA.ANN/Ca.glabriusculum.Chr09.0000000.geno ; 
    cat Chr09.snp | sed -n "1,499999p" > GLA.ANN/Chr09.0000000.snp
    cat Ca.annuum.Chr09.geno | sed -n "1p;${i},+499999p" > GLA.ANN/Ca.annuum.Chr09.$i.geno ; 
    cat Ca.glabriusculum.Chr09.geno | sed -n "1p;${i},+499999p" > GLA.ANN/Ca.glabriusculum.Chr09.$i.geno ; 
    cat Chr09.snp | sed -n "1p;${i},+499999p" > GLA.ANN/Chr09.$i.snp ; 
done
# Run XPCLR
for i in $(seq -w 000000 500000 8390222); do
    XPCLR -xpclr GLA.ANN/Ca.annuum.Chr09.$i.geno \
        GLA.ANN/Ca.glabriusculum.Chr09.$i.geno \
        GLA.ANN/Chr09.$i.snp GLA.ANN/GLA.ANN.Chr09.$i.out \
        -w1 0.0005 100 100 Chr09 -p0 0.7 > GLA.ANN/GLA.ANN.Chr09.$i.out.xpclr.log 
done

# Chr10
# Prepare the input files
for i in $(seq -w 500000 500000 6143906); do
    cat Ca.annuum.Chr10.geno | sed -n "1,499999p" > GLA.ANN/Ca.annuum.Chr10.0000000.geno ; 
    cat Ca.glabriusculum.Chr10.geno | sed -n "1,499999p" > GLA.ANN/Ca.glabriusculum.Chr10.0000000.geno ; 
    cat Chr10.snp | sed -n "1,499999p" > GLA.ANN/Chr10.0000000.snp
    cat Ca.annuum.Chr10.geno | sed -n "1p;${i},+499999p" > GLA.ANN/Ca.annuum.Chr10.$i.geno ; 
    cat Ca.glabriusculum.Chr10.geno | sed -n "1p;${i},+499999p" > GLA.ANN/Ca.glabriusculum.Chr10.$i.geno ; 
    cat Chr10.snp | sed -n "1p;${i},+499999p" > GLA.ANN/Chr10.$i.snp ; 
done
# Run XPCLR
for i in $(seq -w 000000 500000 6143906); do
    XPCLR -xpclr GLA.ANN/Ca.annuum.Chr10.$i.geno \
        GLA.ANN/Ca.glabriusculum.Chr10.$i.geno \
        GLA.ANN/Chr10.$i.snp GLA.ANN/GLA.ANN.Chr10.$i.out \
        -w1 0.0005 100 100 Chr10 -p0 0.7 > GLA.ANN/GLA.ANN.Chr10.$i.out.xpclr.log 
done

# Chr11
# Prepare the input files
for i in $(seq -w 500000 500000 8763242); do
    cat Ca.annuum.Chr11.geno | sed -n "1,499999p" > GLA.ANN/Ca.annuum.Chr11.0000000.geno ; 
    cat Ca.glabriusculum.Chr11.geno | sed -n "1,499999p" > GLA.ANN/Ca.glabriusculum.Chr11.0000000.geno ; 
    cat Chr11.snp | sed -n "1,499999p" > GLA.ANN/Chr11.0000000.snp
    cat Ca.annuum.Chr11.geno | sed -n "1p;${i},+499999p" > GLA.ANN/Ca.annuum.Chr11.$i.geno ; 
    cat Ca.glabriusculum.Chr11.geno | sed -n "1p;${i},+499999p" > GLA.ANN/Ca.glabriusculum.Chr11.$i.geno ; 
    cat Chr11.snp | sed -n "1p;${i},+499999p" > GLA.ANN/Chr11.$i.snp ; 
done
# Run XPCLR
for i in $(seq -w 000000 500000 8763242); do
    XPCLR -xpclr GLA.ANN/Ca.annuum.Chr11.$i.geno \
        GLA.ANN/Ca.glabriusculum.Chr11.$i.geno \
        GLA.ANN/Chr11.$i.snp GLA.ANN/GLA.ANN.Chr11.$i.out \
        -w1 0.0005 100 100 Chr11 -p0 0.7 > GLA.ANN/GLA.ANN.Chr11.$i.out.xpclr.log 
done

# Chr12
# Prepare the input files
for i in $(seq -w 500000 500000 7930008); do
    cat Ca.annuum.Chr12.geno | sed -n "1,499999p" > GLA.ANN/Ca.annuum.Chr12.0000000.geno ; 
    cat Ca.glabriusculum.Chr12.geno | sed -n "1,499999p" > GLA.ANN/Ca.glabriusculum.Chr12.0000000.geno ; 
    cat Chr12.snp | sed -n "1,499999p" > GLA.ANN/Chr12.0000000.snp
    cat Ca.annuum.Chr12.geno | sed -n "1p;${i},+499999p" > GLA.ANN/Ca.annuum.Chr12.$i.geno ; 
    cat Ca.glabriusculum.Chr12.geno | sed -n "1p;${i},+499999p" > GLA.ANN/Ca.glabriusculum.Chr12.$i.geno ; 
    cat Chr12.snp | sed -n "1p;${i},+499999p" > GLA.ANN/Chr12.$i.snp ; 
done
# Run XPCLR
for i in $(seq -w 000000 500000 7930008); do
    XPCLR -xpclr GLA.ANN/Ca.annuum.Chr12.$i.geno \
        GLA.ANN/Ca.glabriusculum.Chr12.$i.geno \
        GLA.ANN/Chr12.$i.snp GLA.ANN/GLA.ANN.Chr12.$i.out \
        -w1 0.0005 100 100 Chr12 -p0 0.7 > GLA.ANN/GLA.ANN.Chr12.$i.out.xpclr.log 
done

# Remove unnecessary files
rm *geno *snp *log

# Prepare the output format for distribution plots
# Average XPCLR score using 100 kb window at 10 kb step size sliding window
for i in $(seq -w 1 12); do
    cat GLA.ANN.Chr$i.*.out.xpclr.txt > GLA.ANN.Chr$i.combined.out.xpclr.txt
    perl ../cluster_xpclrscore.pl GLA.ANN.Chr$i.combined.out.xpclr.txt \
        -wind_size 100000 -wind_step 10000 \
        > GLA.ANN.Chr$i.averaged.out.xpclr.txt
done

awk '{print "Chr01_"$2,1,$2,$3,$7}' GLA.ANN.Chr01.averaged.out.xpclr.txt > GLA.ANN.Chr01.out.xpclr.txt-1
awk '{print "Chr02_"$2,2,$2,$3,$7}' GLA.ANN.Chr02.averaged.out.xpclr.txt > GLA.ANN.Chr02.out.xpclr.txt-1
awk '{print "Chr03_"$2,3,$2,$3,$7}' GLA.ANN.Chr03.averaged.out.xpclr.txt > GLA.ANN.Chr03.out.xpclr.txt-1
awk '{print "Chr04_"$2,4,$2,$3,$7}' GLA.ANN.Chr04.averaged.out.xpclr.txt > GLA.ANN.Chr04.out.xpclr.txt-1
awk '{print "Chr05_"$2,5,$2,$3,$7}' GLA.ANN.Chr05.averaged.out.xpclr.txt > GLA.ANN.Chr05.out.xpclr.txt-1
awk '{print "Chr06_"$2,6,$2,$3,$7}' GLA.ANN.Chr06.averaged.out.xpclr.txt > GLA.ANN.Chr06.out.xpclr.txt-1
awk '{print "Chr07_"$2,7,$2,$3,$7}' GLA.ANN.Chr07.averaged.out.xpclr.txt > GLA.ANN.Chr07.out.xpclr.txt-1
awk '{print "Chr08_"$2,8,$2,$3,$7}' GLA.ANN.Chr08.averaged.out.xpclr.txt > GLA.ANN.Chr08.out.xpclr.txt-1
awk '{print "Chr09_"$2,9,$2,$3,$7}' GLA.ANN.Chr09.averaged.out.xpclr.txt > GLA.ANN.Chr09.out.xpclr.txt-1
awk '{print "Chr10_"$2,10,$2,$3,$7}' GLA.ANN.Chr10.averaged.out.xpclr.txt > GLA.ANN.Chr10.out.xpclr.txt-1
awk '{print "Chr11_"$2,11,$2,$3,$7}' GLA.ANN.Chr11.averaged.out.xpclr.txt > GLA.ANN.Chr11.out.xpclr.txt-1
awk '{print "Chr12_"$2,12,$2,$3,$7}' GLA.ANN.Chr12.averaged.out.xpclr.txt > GLA.ANN.Chr12.out.xpclr.txt-1

# Prepare the output for manhattan plot
cat GLA.ANN.Chr*.out.xpclr.txt-1 | grep -v Wind | sed '1i\SNP Chr Pos PosEnd XPCLR' | sed 's/ /\t/g' > GLA.ANN.out.xpclr.txt

# Remove unnecessary files
rm GLA.ANN.Chr*.out.xpclr.txt-1 GLA.ANN.Chr*.combined.out.xpclr.txt GLA.ANN.Chr*.averaged.out.xpclr.txt

# Filter the selective sweep windows
bedtools merge -d 100000 \
    -i <(cat BAB.BAP.out.xpclr.top5.txt  | sed 's/_/\t/g' | awk '{print $1,$4,$5}' | grep -v XPCLR | sed 's/ /\t/g' ) \
    > BAB.BAP.out.xpclr.top5.merged.txt

bedtools merge -d 10000 \
    -i <(cat BAB.BAP.out.xpclr.top5.txt  | sed 's/_/\t/g' | awk '{print $1,$4,$5}' | grep -v XPCLR | sed 's/ /\t/g' ) \
    > BAB.BAP.out.xpclr.top5.merged.10k.txt

bedtools merge -d 10000 \
    -i <(cat GLA.ANN.out.xpclr.top5.txt  | sed 's/_/\t/g' | awk '{print $1,$4,$5}' | grep -v XPCLR | sed 's/ /\t/g' ) \
    > GLA.ANN.out.xpclr.top5.merged.10k.txt

# Alternatively, use python script
python region.py BAB.BAP.out.xpclr.top5.txt > BAB.BAP.out.xpclr.top5.merged.txt

# Find the overlapped region betwen BAB.BAP and GLA.ANN
bedtools intersect \
    -a <(cut -f1-3 BAB.BAP.out.xpclr.top5.merged.txt | grep -v Start) \
    -b <(cut -f1-3 GLA.ANN.out.xpclr.top5.merged.txt | grep -v Start) \
    > BAB.BAP.GLA.ANN.overlapped.txt

# Find the number of genes within the overlapped regions
bedtools intersect \
    -a ../References/S8/Ca.S8.gff \
    -b BAB.BAP.GLA.ANN.overlapped.txt \
    > BAB.BAP.GLA.ANN.overlapped.genes.gff3

# Find the number of genes within the specific regions
for i in $(ll *merged.txt | awk '{print $9}' | sed 's/.txt//g') ; do
    bedtools intersect \
        -a ../References/S8/Ca.S8.gff\
        -b <(cut -f1-3 $i.txt | grep -v Start) \
        > $i.genes.gff3
done

```

### XP-CLR for C. annuum with all C. annuum var. grabriusculum

```bash
# Working directory
cd /data/zhaojiantao/pepper/XPCLR/GLA.ANN.All

# Blast the markers to S8 to update the physical position
blastall -p blastn -i cM.BP.fa -d ../pepper_genome/fa -m 8 -e 1e-40 -o cM.BP.S8.blast

# Prepare the SNP format for inferring the pseudo genetic position
zcat ../vcf.clean/Chr01.clean.vcf.gz | grep -v "#" | awk '{print "Chr01_"$2,$2}' > Chr01.BP.txt
zcat ../vcf.clean/Chr02.clean.vcf.gz | grep -v "#" | awk '{print "Chr02_"$2,$2}' > Chr02.BP.txt
zcat ../vcf.clean/Chr03.clean.vcf.gz | grep -v "#" | awk '{print "Chr03_"$2,$2}' > Chr03.BP.txt
zcat ../vcf.clean/Chr04.clean.vcf.gz | grep -v "#" | awk '{print "Chr04_"$2,$2}' > Chr04.BP.txt
zcat ../vcf.clean/Chr05.clean.vcf.gz | grep -v "#" | awk '{print "Chr05_"$2,$2}' > Chr05.BP.txt
zcat ../vcf.clean/Chr06.clean.vcf.gz | grep -v "#" | awk '{print "Chr06_"$2,$2}' > Chr06.BP.txt
zcat ../vcf.clean/Chr07.clean.vcf.gz | grep -v "#" | awk '{print "Chr07_"$2,$2}' > Chr07.BP.txt
zcat ../vcf.clean/Chr08.clean.vcf.gz | grep -v "#" | awk '{print "Chr08_"$2,$2}' > Chr08.BP.txt
zcat ../vcf.clean/Chr09.clean.vcf.gz | grep -v "#" | awk '{print "Chr09_"$2,$2}' > Chr09.BP.txt
zcat ../vcf.clean/Chr10.clean.vcf.gz | grep -v "#" | awk '{print "Chr10_"$2,$2}' > Chr10.BP.txt
zcat ../vcf.clean/Chr11.clean.vcf.gz | grep -v "#" | awk '{print "Chr11_"$2,$2}' > Chr11.BP.txt
zcat ../vcf.clean/Chr12.clean.vcf.gz | grep -v "#" | awk '{print "Chr12_"$2,$2}' > Chr12.BP.txt

# Prepare the real genetic maps
for i in $(seq -w 1 12); do 
    cat Linkage.clean.txt | grep Chr$i |  awk '{print $2,$5,$4}' > Chr$i.map.txt ; 
done

# Infer the genetic position 
for i in $(seq -w 1 12); do 
    perl assign_genetic.pl Chr$i.map.txt Chr$i.BP.txt
done

# Compare the relationship between the inferred and real genetic maps
for i in $(seq -w 1 12); do 
    cat Chr$i.map.txt | awk '{print $0,"Real"}' > Chr$i.map-1.txt
    cat Chr$i.snp.genetic.txt | awk '{print $1,$2,$3,"Inferred"}' | \
        cat - Chr$i.map-1.txt | sed '1i\SNP Pos cM Group' | sed 's/ /\t/g' \
        > Chr$i.compare.txt
    rm Chr$i.map-1.txt
done

# SNP position files in LinkageMap folder
for i in $(seq -w 1 12); do 
    zcat ../vcf.clean/Chr$i.clean.vcf.gz | grep -v "#" | awk '{print $1,$2,$4,$5}' > Chr$i.BP-2.txt
    paste Chr$i.snp.genetic.txt Chr$i.BP-2.txt | \
        awk '{print $1,$4,$3,$2,$6,$7}' | sed 's/ /\t/g' > ../XPCLR/Chr$i.snp
    rm Chr$i.BP-2.txt
done

# Prepare the genotypes for xpclr in XPCLR folder
# Make ped file using this chromosome map
for j in Ca.annuum Ca.glabriusculum ; do 
    cat ../../Pop.list/$j.list | awk '{print $1 "\t" $1}' > $j.list.txt
done

for i in $(seq -w 1 12); do
    for j in Ca.glabriusculum ; do 
        # Convert vcf to transposed ped
        plink --vcf ../../vcf.clean/Chr$i.clean.vcf.gz \
            --keep $j.list.txt  \
            --allow-extra-chr --chr Chr$i \
            --recode 01 transpose \
            -output-missing-genotype 9 --double-id \
            --set-missing-var-ids @:# --keep-allele-order \
            --out $j.Chr$i
        # Prepare the transposed ped file to xpclr geno file format
        cut -d " " -f 5- $j.Chr$i.tped  | awk '{print $0" "}' > $j.Chr$i.geno
    done
done

# Remove unnecessary files
rm *log *nosex *tfam *tped 

# XPCLR of Ca.glabriusculum vs Ca.annuum
# Chr01
# Split the genotypes
# Prepare the input files
for i in $(seq -w 500000 500000 9086854); do
    cat Ca.annuum.Chr01.geno | sed -n "1,499999p" > GLA.ANN/Ca.annuum.Chr01.0000000.geno ; 
    cat Ca.glabriusculum.Chr01.geno | sed -n "1,499999p" > GLA.ANN/Ca.glabriusculum.Chr01.0000000.geno ; 
    cat Chr01.snp | sed -n "1,499999p" > GLA.ANN/Chr01.0000000.snp
    cat Ca.annuum.Chr01.geno | sed -n "1p;${i},+499999p" > GLA.ANN/Ca.annuum.Chr01.$i.geno ; 
    cat Ca.glabriusculum.Chr01.geno | sed -n "1p;${i},+499999p" > GLA.ANN/Ca.glabriusculum.Chr01.$i.geno ; 
    cat Chr01.snp | sed -n "1p;${i},+499999p" > GLA.ANN/Chr01.$i.snp ; 
done

# Run XPCLR
for i in $(seq -w 000000 500000 9086854); do
    XPCLR -xpclr GLA.ANN/Ca.annuum.Chr01.$i.geno \
        GLA.ANN/Ca.glabriusculum.Chr01.$i.geno \
        GLA.ANN/Chr01.$i.snp GLA.ANN/GLA.ANN.Chr01.$i.out \
        -w1 0.0005 100 100 Chr01 -p0 0.7 > GLA.ANN/GLA.ANN.Chr01.$i.out.xpclr.log 
done

# Chr02
# Prepare the input files
for i in $(seq -w 500000 500000 4776664); do
    cat Ca.annuum.Chr02.geno | sed -n "1,499999p" > GLA.ANN/Ca.annuum.Chr02.0000000.geno ; 
    cat Ca.glabriusculum.Chr02.geno | sed -n "1,499999p" > GLA.ANN/Ca.glabriusculum.Chr02.0000000.geno ; 
    cat Chr02.snp | sed -n "1,499999p" > GLA.ANN/Chr02.0000000.snp
    cat Ca.annuum.Chr02.geno | sed -n "1p;${i},+499999p" > GLA.ANN/Ca.annuum.Chr02.$i.geno ; 
    cat Ca.glabriusculum.Chr02.geno | sed -n "1p;${i},+499999p" > GLA.ANN/Ca.glabriusculum.Chr02.$i.geno ; 
    cat Chr02.snp | sed -n "1p;${i},+499999p" > GLA.ANN/Chr02.$i.snp ; 
done

# Run XPCLR
    for i in $(seq -w 000000 500000 4776664); do
    XPCLR -xpclr GLA.ANN/Ca.annuum.Chr02.$i.geno \
        GLA.ANN/Ca.glabriusculum.Chr02.$i.geno \
        GLA.ANN/Chr02.$i.snp GLA.ANN/GLA.ANN.Chr02.$i.out \
        -w1 0.0005 100 100 Chr02 -p0 0.7 > GLA.ANN/GLA.ANN.Chr02.$i.out.xpclr.log &
done

# Chr03
# Prepare the input files
for i in $(seq -w 500000 500000 8369414); do
    cat Ca.annuum.Chr03.geno | sed -n "1,499999p" > GLA.ANN/Ca.annuum.Chr03.0000000.geno ; 
    cat Ca.glabriusculum.Chr03.geno | sed -n "1,499999p" > GLA.ANN/Ca.glabriusculum.Chr03.0000000.geno ; 
    cat Chr03.snp | sed -n "1,499999p" > GLA.ANN/Chr03.0000000.snp
    cat Ca.annuum.Chr03.geno | sed -n "1p;${i},+499999p" > GLA.ANN/Ca.annuum.Chr03.$i.geno ; 
    cat Ca.glabriusculum.Chr03.geno | sed -n "1p;${i},+499999p" > GLA.ANN/Ca.glabriusculum.Chr03.$i.geno ; 
    cat Chr03.snp | sed -n "1p;${i},+499999p" > GLA.ANN/Chr03.$i.snp ; 
done
# Run XPCLR
for i in $(seq -w 000000 500000 8369414); do
XPCLR -xpclr GLA.ANN/Ca.annuum.Chr03.$i.geno \
    GLA.ANN/Ca.glabriusculum.Chr03.$i.geno \
    GLA.ANN/Chr03.$i.snp GLA.ANN/GLA.ANN.Chr03.$i.out \
    -w1 0.0005 100 100 Chr03 -p0 0.7 > GLA.ANN/GLA.ANN.Chr03.$i.out.xpclr.log &
done

# Chr04
# Prepare the input files
for i in $(seq -w 500000 500000 7616000); do
    cat Ca.annuum.Chr04.geno | sed -n "1,499999p" > GLA.ANN/Ca.annuum.Chr04.0000000.geno ; 
    cat Ca.glabriusculum.Chr04.geno | sed -n "1,499999p" > GLA.ANN/Ca.glabriusculum.Chr04.0000000.geno ; 
    cat Chr04.snp | sed -n "1,499999p" > GLA.ANN/Chr04.0000000.snp
    cat Ca.annuum.Chr04.geno | sed -n "1p;${i},+499999p" > GLA.ANN/Ca.annuum.Chr04.$i.geno ; 
    cat Ca.glabriusculum.Chr04.geno | sed -n "1p;${i},+499999p" > GLA.ANN/Ca.glabriusculum.Chr04.$i.geno ; 
    cat Chr04.snp | sed -n "1p;${i},+499999p" > GLA.ANN/Chr04.$i.snp ; 
done
# Run XPCLR
for i in $(seq -w 000000 500000 7616000); do
XPCLR -xpclr GLA.ANN/Ca.annuum.Chr04.$i.geno \
    GLA.ANN/Ca.glabriusculum.Chr04.$i.geno \
    GLA.ANN/Chr04.$i.snp GLA.ANN/GLA.ANN.Chr04.$i.out \
    -w1 0.0005 100 100 Chr04 -p0 0.7 > GLA.ANN/GLA.ANN.Chr04.$i.out.xpclr.log 
done 

# Chr05
# Prepare the input files
for i in $(seq -w 500000 500000 7570910); do
    cat Ca.annuum.Chr05.geno | sed -n "1,499999p" > GLA.ANN/Ca.annuum.Chr05.0000000.geno ; 
    cat Ca.glabriusculum.Chr05.geno | sed -n "1,499999p" > GLA.ANN/Ca.glabriusculum.Chr05.0000000.geno ; 
    cat Chr05.snp | sed -n "1,499999p" > GLA.ANN/Chr05.0000000.snp
    cat Ca.annuum.Chr05.geno | sed -n "1p;${i},+499999p" > GLA.ANN/Ca.annuum.Chr05.$i.geno ; 
    cat Ca.glabriusculum.Chr05.geno | sed -n "1p;${i},+499999p" > GLA.ANN/Ca.glabriusculum.Chr05.$i.geno ; 
    cat Chr05.snp | sed -n "1p;${i},+499999p" > GLA.ANN/Chr05.$i.snp ; 
done
# Run XPCLR
for i in $(seq -w 000000 500000 7570910); do
XPCLR -xpclr GLA.ANN/Ca.annuum.Chr05.$i.geno \
    GLA.ANN/Ca.glabriusculum.Chr05.$i.geno \
    GLA.ANN/Chr05.$i.snp GLA.ANN/GLA.ANN.Chr05.$i.out \
    -w1 0.0005 100 100 Chr05 -p0 0.7 > GLA.ANN/GLA.ANN.Chr05.$i.out.xpclr.log &
done

# Chr06
# Prepare the input files
for i in $(seq -w 500000 500000 7522020); do
    cat Ca.annuum.Chr06.geno | sed -n "1,499999p" > GLA.ANN/Ca.annuum.Chr06.0000000.geno ; 
    cat Ca.glabriusculum.Chr06.geno | sed -n "1,499999p" > GLA.ANN/Ca.glabriusculum.Chr06.0000000.geno ; 
    cat Chr06.snp | sed -n "1,499999p" > GLA.ANN/Chr06.0000000.snp
    cat Ca.annuum.Chr06.geno | sed -n "1p;${i},+499999p" > GLA.ANN/Ca.annuum.Chr06.$i.geno ; 
    cat Ca.glabriusculum.Chr06.geno | sed -n "1p;${i},+499999p" > GLA.ANN/Ca.glabriusculum.Chr06.$i.geno ; 
    cat Chr06.snp | sed -n "1p;${i},+499999p" > GLA.ANN/Chr06.$i.snp ; 
done
# Run XPCLR
for i in $(seq -w 000000 500000 7522020); do
    XPCLR -xpclr GLA.ANN/Ca.annuum.Chr06.$i.geno \
        GLA.ANN/Ca.glabriusculum.Chr06.$i.geno \
        GLA.ANN/Chr06.$i.snp GLA.ANN/GLA.ANN.Chr06.$i.out \
        -w1 0.0005 100 100 Chr06 -p0 0.7 > GLA.ANN/GLA.ANN.Chr06.$i.out.xpclr.log &
done

# Chr07
# Prepare the input files
for i in $(seq -w 500000 500000 7824408); do
    cat Ca.annuum.Chr07.geno | sed -n "1,499999p" > GLA.ANN/Ca.annuum.Chr07.0000000.geno ; 
    cat Ca.glabriusculum.Chr07.geno | sed -n "1,499999p" > GLA.ANN/Ca.glabriusculum.Chr07.0000000.geno ; 
    cat Chr07.snp | sed -n "1,499999p" > GLA.ANN/Chr07.0000000.snp
    cat Ca.annuum.Chr07.geno | sed -n "1p;${i},+499999p" > GLA.ANN/Ca.annuum.Chr07.$i.geno ; 
    cat Ca.glabriusculum.Chr07.geno | sed -n "1p;${i},+499999p" > GLA.ANN/Ca.glabriusculum.Chr07.$i.geno ; 
    cat Chr07.snp | sed -n "1p;${i},+499999p" > GLA.ANN/Chr07.$i.snp ; 
done
# Run XPCLR
for i in $(seq -w 000000 500000 7824408); do
    XPCLR -xpclr GLA.ANN/Ca.annuum.Chr07.$i.geno \
        GLA.ANN/Ca.glabriusculum.Chr07.$i.geno \
        GLA.ANN/Chr07.$i.snp GLA.ANN/GLA.ANN.Chr07.$i.out \
        -w1 0.0005 100 100 Chr07 -p0 0.7 > GLA.ANN/GLA.ANN.Chr07.$i.out.xpclr.log &
done

# Chr08
# Prepare the input files
for i in $(seq -w 500000 500000 5099543); do
    cat Ca.annuum.Chr08.geno | sed -n "1,499999p" > GLA.ANN/Ca.annuum.Chr08.0000000.geno ; 
    cat Ca.glabriusculum.Chr08.geno | sed -n "1,499999p" > GLA.ANN/Ca.glabriusculum.Chr08.0000000.geno ; 
    cat Chr08.snp | sed -n "1,499999p" > GLA.ANN/Chr08.0000000.snp
    cat Ca.annuum.Chr08.geno | sed -n "1p;${i},+499999p" > GLA.ANN/Ca.annuum.Chr08.$i.geno ; 
    cat Ca.glabriusculum.Chr08.geno | sed -n "1p;${i},+499999p" > GLA.ANN/Ca.glabriusculum.Chr08.$i.geno ; 
    cat Chr08.snp | sed -n "1p;${i},+499999p" > GLA.ANN/Chr08.$i.snp ; 
done
# Run XPCLR
for i in $(seq -w 000000 500000 5099543); do
    XPCLR -xpclr GLA.ANN/Ca.annuum.Chr08.$i.geno \
        GLA.ANN/Ca.glabriusculum.Chr08.$i.geno \
        GLA.ANN/Chr08.$i.snp GLA.ANN/GLA.ANN.Chr08.$i.out \
        -w1 0.0005 100 1000 Chr08 -p0 0.7 > GLA.ANN/GLA.ANN.Chr08.$i.out.xpclr.log &
done

# Chr09
# Prepare the input files
for i in $(seq -w 500000 500000 8390222); do
    cat Ca.annuum.Chr09.geno | sed -n "1,499999p" > GLA.ANN/Ca.annuum.Chr09.0000000.geno ; 
    cat Ca.glabriusculum.Chr09.geno | sed -n "1,499999p" > GLA.ANN/Ca.glabriusculum.Chr09.0000000.geno ; 
    cat Chr09.snp | sed -n "1,499999p" > GLA.ANN/Chr09.0000000.snp
    cat Ca.annuum.Chr09.geno | sed -n "1p;${i},+499999p" > GLA.ANN/Ca.annuum.Chr09.$i.geno ; 
    cat Ca.glabriusculum.Chr09.geno | sed -n "1p;${i},+499999p" > GLA.ANN/Ca.glabriusculum.Chr09.$i.geno ; 
    cat Chr09.snp | sed -n "1p;${i},+499999p" > GLA.ANN/Chr09.$i.snp ; 
done
# Run XPCLR
for i in $(seq -w 000000 500000 8390222); do
    XPCLR -xpclr GLA.ANN/Ca.annuum.Chr09.$i.geno \
        GLA.ANN/Ca.glabriusculum.Chr09.$i.geno \
        GLA.ANN/Chr09.$i.snp GLA.ANN/GLA.ANN.Chr09.$i.out \
        -w1 0.0005 100 100 Chr09 -p0 0.7 > GLA.ANN/GLA.ANN.Chr09.$i.out.xpclr.log &
done

# Chr10
# Prepare the input files
for i in $(seq -w 500000 500000 6143906); do
    cat Ca.annuum.Chr10.geno | sed -n "1,499999p" > GLA.ANN/Ca.annuum.Chr10.0000000.geno ; 
    cat Ca.glabriusculum.Chr10.geno | sed -n "1,499999p" > GLA.ANN/Ca.glabriusculum.Chr10.0000000.geno ; 
    cat Chr10.snp | sed -n "1,499999p" > GLA.ANN/Chr10.0000000.snp
    cat Ca.annuum.Chr10.geno | sed -n "1p;${i},+499999p" > GLA.ANN/Ca.annuum.Chr10.$i.geno ; 
    cat Ca.glabriusculum.Chr10.geno | sed -n "1p;${i},+499999p" > GLA.ANN/Ca.glabriusculum.Chr10.$i.geno ; 
    cat Chr10.snp | sed -n "1p;${i},+499999p" > GLA.ANN/Chr10.$i.snp ; 
done
# Run XPCLR
for i in $(seq -w 000000 500000 6143906); do
    XPCLR -xpclr GLA.ANN/Ca.annuum.Chr10.$i.geno \
        GLA.ANN/Ca.glabriusculum.Chr10.$i.geno \
        GLA.ANN/Chr10.$i.snp GLA.ANN/GLA.ANN.Chr10.$i.out \
        -w1 0.0005 100 100 Chr10 -p0 0.7 > GLA.ANN/GLA.ANN.Chr10.$i.out.xpclr.log &
done

# Chr11
# Prepare the input files
for i in $(seq -w 500000 500000 8763242); do
    cat Ca.annuum.Chr11.geno | sed -n "1,499999p" > GLA.ANN/Ca.annuum.Chr11.0000000.geno ; 
    cat Ca.glabriusculum.Chr11.geno | sed -n "1,499999p" > GLA.ANN/Ca.glabriusculum.Chr11.0000000.geno ; 
    cat Chr11.snp | sed -n "1,499999p" > GLA.ANN/Chr11.0000000.snp
    cat Ca.annuum.Chr11.geno | sed -n "1p;${i},+499999p" > GLA.ANN/Ca.annuum.Chr11.$i.geno ; 
    cat Ca.glabriusculum.Chr11.geno | sed -n "1p;${i},+499999p" > GLA.ANN/Ca.glabriusculum.Chr11.$i.geno ; 
    cat Chr11.snp | sed -n "1p;${i},+499999p" > GLA.ANN/Chr11.$i.snp ; 
done
# Run XPCLR
for i in $(seq -w 000000 500000 8763242); do
    XPCLR -xpclr GLA.ANN/Ca.annuum.Chr11.$i.geno \
        GLA.ANN/Ca.glabriusculum.Chr11.$i.geno \
        GLA.ANN/Chr11.$i.snp GLA.ANN/GLA.ANN.Chr11.$i.out \
        -w1 0.0005 100 100 Chr11 -p0 0.7 > GLA.ANN/GLA.ANN.Chr11.$i.out.xpclr.log &
done

# Chr12
# Prepare the input files
for i in $(seq -w 500000 500000 7930008); do
    cat Ca.annuum.Chr12.geno | sed -n "1,499999p" > GLA.ANN/Ca.annuum.Chr12.0000000.geno ; 
    cat Ca.glabriusculum.Chr12.geno | sed -n "1,499999p" > GLA.ANN/Ca.glabriusculum.Chr12.0000000.geno ; 
    cat Chr12.snp | sed -n "1,499999p" > GLA.ANN/Chr12.0000000.snp
    cat Ca.annuum.Chr12.geno | sed -n "1p;${i},+499999p" > GLA.ANN/Ca.annuum.Chr12.$i.geno ; 
    cat Ca.glabriusculum.Chr12.geno | sed -n "1p;${i},+499999p" > GLA.ANN/Ca.glabriusculum.Chr12.$i.geno ; 
    cat Chr12.snp | sed -n "1p;${i},+499999p" > GLA.ANN/Chr12.$i.snp ; 
done
# Run XPCLR
for i in $(seq -w 000000 500000 7930008); do
    XPCLR -xpclr GLA.ANN/Ca.annuum.Chr12.$i.geno \
        GLA.ANN/Ca.glabriusculum.Chr12.$i.geno \
        GLA.ANN/Chr12.$i.snp GLA.ANN/GLA.ANN.Chr12.$i.out \
        -w1 0.0005 100 100 Chr12 -p0 0.7 > GLA.ANN/GLA.ANN.Chr12.$i.out.xpclr.log &
done

# Remove unnecessary files
rm *geno *snp *log

# Prepare the output format for distribution plots
# Average XPCLR score using 100 kb window at 10 kb step size sliding window
for i in $(seq -w 1 12); do
    cat GLA.ANN/GLA.ANN.Chr$i.*.out.xpclr.txt > GLA.ANN.Chr$i.combined.out.xpclr.txt
    perl ../cluster_xpclrscore.pl GLA.ANN.Chr$i.combined.out.xpclr.txt \
        -wind_size 100000 -wind_step 10000 \
        > GLA.ANN.Chr$i.averaged.out.xpclr.txt
done

awk '{print "Chr01_"$2,1,$2,$3,$7}' GLA.ANN.Chr01.averaged.out.xpclr.txt > GLA.ANN.Chr01.out.xpclr.txt-1
awk '{print "Chr02_"$2,2,$2,$3,$7}' GLA.ANN.Chr02.averaged.out.xpclr.txt > GLA.ANN.Chr02.out.xpclr.txt-1
awk '{print "Chr03_"$2,3,$2,$3,$7}' GLA.ANN.Chr03.averaged.out.xpclr.txt > GLA.ANN.Chr03.out.xpclr.txt-1
awk '{print "Chr04_"$2,4,$2,$3,$7}' GLA.ANN.Chr04.averaged.out.xpclr.txt > GLA.ANN.Chr04.out.xpclr.txt-1
awk '{print "Chr05_"$2,5,$2,$3,$7}' GLA.ANN.Chr05.averaged.out.xpclr.txt > GLA.ANN.Chr05.out.xpclr.txt-1
awk '{print "Chr06_"$2,6,$2,$3,$7}' GLA.ANN.Chr06.averaged.out.xpclr.txt > GLA.ANN.Chr06.out.xpclr.txt-1
awk '{print "Chr07_"$2,7,$2,$3,$7}' GLA.ANN.Chr07.averaged.out.xpclr.txt > GLA.ANN.Chr07.out.xpclr.txt-1
awk '{print "Chr08_"$2,8,$2,$3,$7}' GLA.ANN.Chr08.averaged.out.xpclr.txt > GLA.ANN.Chr08.out.xpclr.txt-1
awk '{print "Chr09_"$2,9,$2,$3,$7}' GLA.ANN.Chr09.averaged.out.xpclr.txt > GLA.ANN.Chr09.out.xpclr.txt-1
awk '{print "Chr10_"$2,10,$2,$3,$7}' GLA.ANN.Chr10.averaged.out.xpclr.txt > GLA.ANN.Chr10.out.xpclr.txt-1
awk '{print "Chr11_"$2,11,$2,$3,$7}' GLA.ANN.Chr11.averaged.out.xpclr.txt > GLA.ANN.Chr11.out.xpclr.txt-1
awk '{print "Chr12_"$2,12,$2,$3,$7}' GLA.ANN.Chr12.averaged.out.xpclr.txt > GLA.ANN.Chr12.out.xpclr.txt-1

# Prepare the output for manhattan plot
cat GLA.ANN.Chr*.out.xpclr.txt-1 | grep -v Wind | sed '1i\SNP Chr Pos PosEnd XPCLR' | sed 's/ /\t/g' > GLA.all.ANN.out.xpclr.txt

# Remove unnecessary files
rm GLA.ANN.Chr*.out.xpclr.txt-1 GLA.ANN.Chr*.combined.out.xpclr.txt GLA.ANN.Chr*.averaged.out.xpclr.txt

# Extract the top 5% regions
grep "SNP" GLA.all.ANN.out.xpclr.txt | \
    cat - <(grep -v "SNP" GLA.all.ANN.out.xpclr.txt | sort -nrk5 | head -14870 | sort -nk2 -nk3) \
    > GLA.all.ANN.out.xpclr.top5.txt

# Merge the top 5% regions
python ../region.py GLA.all.ANN.out.xpclr.top5.txt GLA.all.ANN.out.xpclr.top5.merged.txt
# Alternatively use perl script
perl ../merge.pl GLA.all.ANN.out.xpclr.top5.txt  > GLA.all.ANN.out.xpclr.top5.merged2.txt


# Filter the selective sweep windows
bedtools merge -d 100000 \
    -i <(cat BAB.BAP.out.xpclr.top5.txt  | sed 's/_/\t/g' | awk '{print $1,$4,$5}' | grep -v XPCLR | sed 's/ /\t/g' ) \
    > BAB.BAP.out.xpclr.top5.merged.txt

# Find the overlapped region betwen BAB.BAP and GLA.ANN
bedtools intersect \
    -a <(cut -f1-3 ../BAB.BAP.out.xpclr.top5.merged.txt | grep -v Start) \
    -b <(cut -f1-3 GLA.all.ANN.out.xpclr.top5.merged.txt | grep -v Start) \
    > BAB.BAP.GLA.ANN.overlapped.txt

# Find the number of genes with C. annuum
bedtools intersect \
    -a ../../References/S8/Ca.S8.gff \
    -b <(cut -f1-3 GLA.all.ANN.out.xpclr.top5.merged.txt | grep -v Chromsome) \
    > GLA.all.ANN.genes.gff3

# Number of genes
grep ID=gene  GLA.all.ANN.genes.gff3 | wc

# Find the number of genes within the overlapped regions
bedtools intersect \
    -a ../../References/S8/Ca.S8.gff \
    -b BAB.BAP.GLA.ANN.overlapped.txt \
    > BAB.BAP.GLA.ANN.overlapped.genes.gff3

# Find the SNPs within the overlapped regions
# GLA.ANN
bedtools intersect -a <(cut -f1 GLA.all.ANN.out.xpclr.txt | \
    sed 's/_/\t/g' | paste - GLA.all.ANN.out.xpclr.txt | \
    awk '{print $1,$2,$6,$3}' | grep -v SNP | sed 's/ /\t/g') -b BAB.BAP.GLA.ANN.overlapped.txt | cut -f4 \
    > GLA.all.ANN.overlapped.snps.txt

# BAB.BAP
bedtools intersect -a <(cut -f1 ../BAB.BAP.out.xpclr.txt | \
    sed 's/_/\t/g' | paste - ../BAB.BAP.out.xpclr.txt | \
    awk '{print $1,$2,$6,$3}' | grep -v SNP | sed 's/ /\t/g') -b BAB.BAP.GLA.ANN.overlapped.txt | cut -f4 \
    > BAB.BAP.overlapped.snps.txt

# Find the number of genes within the specific regions
for i in $(ll *merged.txt | awk '{print $9}' | sed 's/.txt//g') ; do
    bedtools intersect \
        -a ../References/S8/Ca.S8.gff\
        -b <(cut -f1-3 $i.txt | grep -v Start) \
        > $i.genes.gff3
done

```

### XP-CLR for C. bac. baccatum versus C. bac. pendulum

```bash
# Working directory
cd /data/zhaojiantao/pepper/XPCLR/BAB.BAP

# Prepare the genotypes
for i in $(seq -w 1 12); do
    for j in Ca.bac.baccatum Ca.bac.pendulum; do 
        cat ../../Pop.list/$j.list | awk '{print $1 "\t" $1}' > $j.list.txt
        # Convert vcf to transposed ped
        plink --vcf ../../vcf.clean/Chr$i.clean.vcf.gz \
            --keep $j.list.txt  \
            --allow-extra-chr --chr Chr$i \
            --recode 01 transpose \
            --geno 0.5 \
            --maf 0.05 \
            -output-missing-genotype 9 --double-id \
            --set-missing-var-ids @:# --keep-allele-order \
            --out $j.Chr$i
    done
done

# Find the common SNPs between C. bac. baccatum and C. bac. pendulum
for i in $(seq -w 1 12); do
    cut -f2 -d ' ' Ca.bac.pendulum.Chr$i.tped > BAP.Chr$i.snp
    cut -f2 -d ' ' Ca.bac.baccatum.Chr$i.tped > BAB.Chr$i.snp
    comm -12 <(sort BAP.Chr$i.snp) <(sort BAB.Chr$i.snp) > BAB.BAP.Chr$i.snp
    grep -Fwf BAB.BAP.Chr$i.snp Ca.bac.pendulum.Chr$i.tped | \
        cut -d " " -f 5- | awk '{print $0" "}' > Ca.bac.pendulum.Chr$i.geno
    grep -Fwf BAB.BAP.Chr$i.snp Ca.bac.baccatum.Chr$i.tped | \
        cut -d " " -f 5- | awk '{print $0" "}' > Ca.bac.baccatum.Chr$i.geno
    grep -Fwf BAB.BAP.Chr$i.snp ../Chr$i.snp > Chr$i.snp
done

# Remove unnecessary files
rm *log *nosex *tfam *tped 

# XPCLR of Ca.bac.baccatum vs Ca.bac.pendulum
# Chr01
# Prepare the input files
for i in $(seq -w 50000 50000 740681); do
    cat Ca.bac.pendulum.Chr01.geno | sed -n "1,49999p" > Ca.bac.pendulum.Chr01.000000.geno ; 
    cat Ca.bac.baccatum.Chr01.geno | sed -n "1,49999p" > Ca.bac.baccatum.Chr01.000000.geno ; 
    cat Chr01.snp | sed -n "1,49999p" > Chr01.000000.snp
    cat Ca.bac.pendulum.Chr01.geno | sed -n "1p;${i},+49999p" > Ca.bac.pendulum.Chr01.$i.geno ; 
    cat Ca.bac.baccatum.Chr01.geno | sed -n "1p;${i},+49999p" > Ca.bac.baccatum.Chr01.$i.geno ; 
    cat Chr01.snp | sed -n "1p;${i},+49999p" > Chr01.$i.snp ; 
done
# Run XPCLR
for i in $(seq -w 00000 50000 740681); do
    XPCLR -xpclr Ca.bac.pendulum.Chr01.$i.geno \
        Ca.bac.baccatum.Chr01.$i.geno \
        Chr01.$i.snp BAB.BAP.Chr01.$i.out \
        -w1 0.0005 100 100 Chr01 -p0 0.7 > BAB.BAP.Chr01.$i.out.xpclr.log 
done

# Chr02
# Prepare the input files
for i in $(seq -w 50000 50000 291463); do
    cat Ca.bac.pendulum.Chr02.geno | sed -n "1,49999p" > Ca.bac.pendulum.Chr02.000000.geno ; 
    cat Ca.bac.baccatum.Chr02.geno | sed -n "1,49999p" > Ca.bac.baccatum.Chr02.000000.geno ; 
    cat Chr02.snp | sed -n "1,49999p" > Chr02.000000.snp
    cat Ca.bac.pendulum.Chr02.geno | sed -n "1p;${i},+49999p" > Ca.bac.pendulum.Chr02.$i.geno ; 
    cat Ca.bac.baccatum.Chr02.geno | sed -n "1p;${i},+49999p" > Ca.bac.baccatum.Chr02.$i.geno ; 
    cat Chr02.snp | sed -n "1p;${i},+49999p" > Chr02.$i.snp ; 
done
# Run XPCLR
for i in $(seq -w 00000 50000 291463); do
    XPCLR -xpclr Ca.bac.pendulum.Chr02.$i.geno \
        Ca.bac.baccatum.Chr02.$i.geno \
        Chr02.$i.snp BAB.BAP.Chr02.$i.out \
        -w1 0.0005 100 100 Chr02 -p0 0.7 > BAB.BAP.Chr02.$i.out.xpclr.log 
done

# Chr03
# Prepare the input files
for i in $(seq -w 50000 50000 521851); do
    cat Ca.bac.pendulum.Chr03.geno | sed -n "1,49999p" > Ca.bac.pendulum.Chr03.000000.geno ; 
    cat Ca.bac.baccatum.Chr03.geno | sed -n "1,49999p" > Ca.bac.baccatum.Chr03.000000.geno ; 
    cat Chr03.snp | sed -n "1,49999p" > Chr03.000000.snp
    cat Ca.bac.pendulum.Chr03.geno | sed -n "1p;${i},+49999p" > Ca.bac.pendulum.Chr03.$i.geno ; 
    cat Ca.bac.baccatum.Chr03.geno | sed -n "1p;${i},+49999p" > Ca.bac.baccatum.Chr03.$i.geno ; 
    cat Chr03.snp | sed -n "1p;${i},+49999p" > Chr03.$i.snp ; 
done
# Run XPCLR
for i in $(seq -w 00000 50000 521851); do
XPCLR -xpclr Ca.bac.pendulum.Chr03.$i.geno \
    Ca.bac.baccatum.Chr03.$i.geno \
    Chr03.$i.snp BAB.BAP.Chr03.$i.out \
    -w1 0.0005 100 100 Chr03 -p0 0.7 > BAB.BAP.Chr03.$i.out.xpclr.log 
done

# Chr04
# Prepare the input files
for i in $(seq -w 50000 50000 563213); do
cat Ca.bac.pendulum.Chr04.geno | sed -n "1,49999p" > Ca.bac.pendulum.Chr04.000000.geno ; 
cat Ca.bac.baccatum.Chr04.geno | sed -n "1,49999p" > Ca.bac.baccatum.Chr04.000000.geno ; 
cat Chr04.snp | sed -n "1,499999p" > Chr04.000000.snp
cat Ca.bac.pendulum.Chr04.geno | sed -n "1p;${i},+49999p" > Ca.bac.pendulum.Chr04.$i.geno ; 
cat Ca.bac.baccatum.Chr04.geno | sed -n "1p;${i},+49999p" > Ca.bac.baccatum.Chr04.$i.geno ; 
cat Chr04.snp | sed -n "1p;${i},+49999p" > Chr04.$i.snp ; 
done
# Run XPCLR
for i in $(seq -w 00000 50000 563213); do
XPCLR -xpclr Ca.bac.pendulum.Chr04.$i.geno \
    Ca.bac.baccatum.Chr04.$i.geno \
    Chr04.$i.snp BAB.BAP.Chr04.$i.out \
    -w1 0.0005 100 100 Chr04 -p0 0.7 > BAB.BAP.Chr04.$i.out.xpclr.log 
done 

# Chr05
# Prepare the input files
for i in $(seq -w 50000 50000 538616); do
cat Ca.bac.pendulum.Chr05.geno | sed -n "1,49999p" > Ca.bac.pendulum.Chr05.000000.geno ; 
cat Ca.bac.baccatum.Chr05.geno | sed -n "1,49999p" > Ca.bac.baccatum.Chr05.000000.geno ; 
cat Chr05.snp | sed -n "1,49999p" > Chr05.000000.snp
cat Ca.bac.pendulum.Chr05.geno | sed -n "1p;${i},+49999p" > Ca.bac.pendulum.Chr05.$i.geno ; 
cat Ca.bac.baccatum.Chr05.geno | sed -n "1p;${i},+49999p" > Ca.bac.baccatum.Chr05.$i.geno ; 
cat Chr05.snp | sed -n "1p;${i},+49999p" > Chr05.$i.snp ; 
done
# Run XPCLR
for i in $(seq -w 00000 50000 538616); do
XPCLR -xpclr Ca.bac.pendulum.Chr05.$i.geno \
    Ca.bac.baccatum.Chr05.$i.geno \
    Chr05.$i.snp BAB.BAP.Chr05.$i.out \
    -w1 0.0005 100 100 Chr05 -p0 0.7 > BAB.BAP.Chr05.$i.out.xpclr.log 
done

# Chr06
# Prepare the input files
for i in $(seq -w 50000 50000 548970); do
cat Ca.bac.pendulum.Chr06.geno | sed -n "1,49999p" > Ca.bac.pendulum.Chr06.000000.geno ; 
cat Ca.bac.baccatum.Chr06.geno | sed -n "1,49999p" > Ca.bac.baccatum.Chr06.000000.geno ; 
cat Chr06.snp | sed -n "1,49999p" > Chr06.000000.snp
cat Ca.bac.pendulum.Chr06.geno | sed -n "1p;${i},+49999p" > Ca.bac.pendulum.Chr06.$i.geno ; 
cat Ca.bac.baccatum.Chr06.geno | sed -n "1p;${i},+49999p" > Ca.bac.baccatum.Chr06.$i.geno ; 
cat Chr06.snp | sed -n "1p;${i},+49999p" > Chr06.$i.snp ; 
done
# Run XPCLR
for i in $(seq -w 00000 50000 548970); do
XPCLR -xpclr Ca.bac.pendulum.Chr06.$i.geno \
    Ca.bac.baccatum.Chr06.$i.geno \
    Chr06.$i.snp BAB.BAP.Chr06.$i.out \
    -w1 0.0005 100 100 Chr06 -p0 0.7 > BAB.BAP.Chr06.$i.out.xpclr.log 
done

# Chr07
# Prepare the input files
for i in $(seq -w 50000 50000 487786); do
cat Ca.bac.pendulum.Chr07.geno | sed -n "1,49999p" > Ca.bac.pendulum.Chr07.000000.geno ; 
cat Ca.bac.baccatum.Chr07.geno | sed -n "1,49999p" > Ca.bac.baccatum.Chr07.000000.geno ; 
cat Chr07.snp | sed -n "1,49999p" > Chr07.000000.snp
cat Ca.bac.pendulum.Chr07.geno | sed -n "1p;${i},+49999p" > Ca.bac.pendulum.Chr07.$i.geno ; 
cat Ca.bac.baccatum.Chr07.geno | sed -n "1p;${i},+49999p" > Ca.bac.baccatum.Chr07.$i.geno ; 
cat Chr07.snp | sed -n "1p;${i},+49999p" > Chr07.$i.snp ; 
done
# Run XPCLR
for i in $(seq -w 00000 50000 487786); do
XPCLR -xpclr Ca.bac.pendulum.Chr07.$i.geno \
    Ca.bac.baccatum.Chr07.$i.geno \
    Chr07.$i.snp BAB.BAP.Chr07.$i.out \
    -w1 0.0005 100 100 Chr07 -p0 0.7 > BAB.BAP.Chr07.$i.out.xpclr.log 
done

# Chr08
# Prepare the input files
for i in $(seq -w 50000 50000 238001); do
cat Ca.bac.pendulum.Chr08.geno | sed -n "1,49999p" > Ca.bac.pendulum.Chr08.000000.geno ; 
cat Ca.bac.baccatum.Chr08.geno | sed -n "1,49999p" > Ca.bac.baccatum.Chr08.000000.geno ; 
cat Chr08.snp | sed -n "1,49999p" > Chr08.000000.snp
cat Ca.bac.pendulum.Chr08.geno | sed -n "1p;${i},+49999p" > Ca.bac.pendulum.Chr08.$i.geno ; 
cat Ca.bac.baccatum.Chr08.geno | sed -n "1p;${i},+49999p" > Ca.bac.baccatum.Chr08.$i.geno ; 
cat Chr08.snp | sed -n "1p;${i},+49999p" > Chr08.$i.snp ; 
done
# Run XPCLR
for i in $(seq -w 00000 50000 238001); do
XPCLR -xpclr Ca.bac.pendulum.Chr08.$i.geno \
    Ca.bac.baccatum.Chr08.$i.geno \
    Chr08.$i.snp BAB.BAP.Chr08.$i.out \
    -w1 0.0005 100 100 Chr08 -p0 0.7 > BAB.BAP.Chr08.$i.out.xpclr.log 
done

# Chr09
# Prepare the input files
for i in $(seq -w 50000 50000 482094); do
cat Ca.bac.pendulum.Chr09.geno | sed -n "1,49999p" > Ca.bac.pendulum.Chr09.000000.geno ; 
cat Ca.bac.baccatum.Chr09.geno | sed -n "1,49999p" > Ca.bac.baccatum.Chr09.000000.geno ; 
cat Chr09.snp | sed -n "1,49999p" > Chr09.000000.snp
cat Ca.bac.pendulum.Chr09.geno | sed -n "1p;${i},+49999p" > Ca.bac.pendulum.Chr09.$i.geno ; 
cat Ca.bac.baccatum.Chr09.geno | sed -n "1p;${i},+49999p" > Ca.bac.baccatum.Chr09.$i.geno ; 
cat Chr09.snp | sed -n "1p;${i},+49999p" > Chr09.$i.snp ; 
done
# Run XPCLR
for i in $(seq -w 00000 50000 482094); do
XPCLR -xpclr Ca.bac.pendulum.Chr09.$i.geno \
    Ca.bac.baccatum.Chr09.$i.geno \
    Chr09.$i.snp BAB.BAP.Chr09.$i.out \
    -w1 0.0005 100 100 Chr09 -p0 0.7 > BAB.BAP.Chr09.$i.out.xpclr.log 
done

# Chr10
# Prepare the input files
for i in $(seq -w 50000 50000 583266); do
cat Ca.bac.pendulum.Chr10.geno | sed -n "1,49999p" > Ca.bac.pendulum.Chr10.000000.geno ; 
cat Ca.bac.baccatum.Chr10.geno | sed -n "1,49999p" > Ca.bac.baccatum.Chr10.000000.geno ; 
cat Chr10.snp | sed -n "1,49999p" > Chr10.000000.snp
cat Ca.bac.pendulum.Chr10.geno | sed -n "1p;${i},+49999p" > Ca.bac.pendulum.Chr10.$i.geno ; 
cat Ca.bac.baccatum.Chr10.geno | sed -n "1p;${i},+49999p" > Ca.bac.baccatum.Chr10.$i.geno ; 
cat Chr10.snp | sed -n "1p;${i},+49999p" > Chr10.$i.snp ; 
done
# Run XPCLR
for i in $(seq -w 00000 50000 583266); do
XPCLR -xpclr Ca.bac.pendulum.Chr10.$i.geno \
    Ca.bac.baccatum.Chr10.$i.geno \
    Chr10.$i.snp BAB.BAP.Chr10.$i.out \
    -w1 0.0005 100 100 Chr10 -p0 0.7 > BAB.BAP.Chr10.$i.out.xpclr.log 
done

# Chr11
# Prepare the input files
for i in $(seq -w 50000 50000 488040); do
cat Ca.bac.pendulum.Chr11.geno | sed -n "1,49999p" > Ca.bac.pendulum.Chr11.000000.geno ; 
cat Ca.bac.baccatum.Chr11.geno | sed -n "1,49999p" > Ca.bac.baccatum.Chr11.000000.geno ; 
cat Chr11.snp | sed -n "1,49999p" > Chr11.000000.snp
cat Ca.bac.pendulum.Chr11.geno | sed -n "1p;${i},+49999p" > Ca.bac.pendulum.Chr11.$i.geno ; 
cat Ca.bac.baccatum.Chr11.geno | sed -n "1p;${i},+49999p" > Ca.bac.baccatum.Chr11.$i.geno ; 
cat Chr11.snp | sed -n "1p;${i},+49999p" > Chr11.$i.snp ; 
done
# Run XPCLR
for i in $(seq -w 00000 50000 488040); do
XPCLR -xpclr Ca.bac.pendulum.Chr11.$i.geno \
    Ca.bac.baccatum.Chr11.$i.geno \
    Chr11.$i.snp BAB.BAP.Chr11.$i.out \
    -w1 0.0005 100 100 Chr11 -p0 0.7 > BAB.BAP.Chr11.$i.out.xpclr.log 
done

# Chr12
# Prepare the input files
for i in $(seq -w 50000 50000 507900); do
cat Ca.bac.pendulum.Chr12.geno | sed -n "1,49999p" > Ca.bac.pendulum.Chr12.000000.geno ; 
cat Ca.bac.baccatum.Chr12.geno | sed -n "1,49999p" > Ca.bac.baccatum.Chr12.000000.geno ; 
cat Chr12.snp | sed -n "1,49999p" > Chr12.000000.snp
cat Ca.bac.pendulum.Chr12.geno | sed -n "1p;${i},+49999p" > Ca.bac.pendulum.Chr12.$i.geno ; 
cat Ca.bac.baccatum.Chr12.geno | sed -n "1p;${i},+49999p" > Ca.bac.baccatum.Chr12.$i.geno ; 
cat Chr12.snp | sed -n "1p;${i},+49999p" > Chr12.$i.snp ; 
done
# Run XPCLR
for i in $(seq -w 00000 50000 507900); do
XPCLR -xpclr Ca.bac.pendulum.Chr12.$i.geno \
    Ca.bac.baccatum.Chr12.$i.geno \
    Chr12.$i.snp BAB.BAP.Chr12.$i.out \
    -w1 0.0005 100 100 Chr12 -p0 0.7 > BAB.BAP.Chr12.$i.out.xpclr.log 
done

# Remove unnecessary files
rm *geno *snp *log

# Prepare the output format for distribution plots
# Average XPCLR score using 100 kb window at 10 kb step size sliding window
for i in $(seq -w 1 12); do
    cat BAB.BAP.Chr$i.*out.xpclr.txt > BAB.BAP.Chr$i.combined.out.xpclr.txt
    perl ../cluster_xpclrscore.pl BAB.BAP.Chr$i.combined.out.xpclr.txt \
        -wind_size 100000 -wind_step 10000 \
        > BAB.BAP.Chr$i.averaged.out.xpclr.txt
done

awk '{print "Chr01_"$2,1,$2,$3,$7}' BAB.BAP.Chr01.averaged.out.xpclr.txt > BAB.BAP.Chr01.out.xpclr.txt-1
awk '{print "Chr02_"$2,2,$2,$3,$7}' BAB.BAP.Chr02.averaged.out.xpclr.txt > BAB.BAP.Chr02.out.xpclr.txt-1
awk '{print "Chr03_"$2,3,$2,$3,$7}' BAB.BAP.Chr03.averaged.out.xpclr.txt > BAB.BAP.Chr03.out.xpclr.txt-1
awk '{print "Chr04_"$2,4,$2,$3,$7}' BAB.BAP.Chr04.averaged.out.xpclr.txt > BAB.BAP.Chr04.out.xpclr.txt-1
awk '{print "Chr05_"$2,5,$2,$3,$7}' BAB.BAP.Chr05.averaged.out.xpclr.txt > BAB.BAP.Chr05.out.xpclr.txt-1
awk '{print "Chr06_"$2,6,$2,$3,$7}' BAB.BAP.Chr06.averaged.out.xpclr.txt > BAB.BAP.Chr06.out.xpclr.txt-1
awk '{print "Chr07_"$2,7,$2,$3,$7}' BAB.BAP.Chr07.averaged.out.xpclr.txt > BAB.BAP.Chr07.out.xpclr.txt-1
awk '{print "Chr08_"$2,8,$2,$3,$7}' BAB.BAP.Chr08.averaged.out.xpclr.txt > BAB.BAP.Chr08.out.xpclr.txt-1
awk '{print "Chr09_"$2,9,$2,$3,$7}' BAB.BAP.Chr09.averaged.out.xpclr.txt > BAB.BAP.Chr09.out.xpclr.txt-1
awk '{print "Chr10_"$2,10,$2,$3,$7}' BAB.BAP.Chr10.averaged.out.xpclr.txt > BAB.BAP.Chr10.out.xpclr.txt-1
awk '{print "Chr11_"$2,11,$2,$3,$7}' BAB.BAP.Chr11.averaged.out.xpclr.txt > BAB.BAP.Chr11.out.xpclr.txt-1
awk '{print "Chr12_"$2,12,$2,$3,$7}' BAB.BAP.Chr12.averaged.out.xpclr.txt > BAB.BAP.Chr12.out.xpclr.txt-1

# Prepare the output for manhattan plot
cat BAB.BAP.Chr*.out.xpclr.txt-1 | \
    grep -v Wind | \
    sed '1i\SNP Chr Pos PosEnd XPCLR' | \
    sed 's/ /\t/g' \
    > BAB.BAP.out.xpclr.txt

# Remove unnecessary files
rm BAB.BAP.Chr*.out.xpclr.txt-1 BAB.BAP.Chr*.combined.out.xpclr.txt BAB.BAP.Chr*.averaged.out.xpclr.txt

```

### Enrichment analysis

```bash
# Working directory
cd /data/zhaojiantao/pepper/References/S8

# Genome gene list
less pepper_S8.gff3.gz | grep gene | cut -f9 | sed 's/;/\t/g' | \
    cut -f1 | grep rna | sed 's/ID=//g' | sed 's/Parent=//g' \
    > S8.genome.gene.list.txt

# Prepare the gene lists
for i in BAB.BAP.0.05.divergent.genes GLA.ANN.0.05.divergent.genes \
    GLA.ANN.0.05.narrow.convergent.genes GLA.ANN.0.05.broad.convergent.genes ; do
    grep gene $i | cut -f9 | sed 's/;/.1;/g' | sed 's/;/\t/g' | \
        cut -f1 | sed 's/ID=//g' | sed 's/Parent=//g' | sed 's/gene/rna/g' \
        > $i.list.txt
done

```



## Effective population size
Software link: [SMC++](https://github.com/popgenmethods/smcpp#tips-for-using-smc)

```bash
# Working directory
cd /data/zhaojiantao/pepper/SMC++

# Install 
pip install git+https://github.com/popgenmethods/smcpp

# Working directory
cd /data/zhaojiantao/pepper/SMC++/GLA

for i in $(seq -w 1 12); do
    for j in $(cat GLA.list) ; do
        # Extract the tested samples from all vcf
        vcftools --gzvcf /data/zhaojiantao/pepper/vcf.clean/Chr$i.clean.vcf.gz \
            --keep GLA.list \
            --maf 0.05 \
            --max-missing 0.8 \
            --recode --out Chr$i.clean
        # bgzip the vcf, not gzip
        bgzip Chr$i.clean.recode.vcf
        # Index the file
        tabix Chr$i.clean.recode.vcf.gz
        # Convert vcf to smc 
        smc++ vcf2smc --cores 60 --missing-cutoff 5000 -d $j $j \
            Chr$i.clean.recode.vcf.gz Chr$i.$j.smc.gz Chr$i \
            GLA:Grif_9241,PI_555616,PI_632930,Grif_9142,PI_311126,PI_511885,PI_593557,Grif_9249,PI_593526,PI_439311
    done
done

# Run smc++
smc++ estimate -o smc.out 6.96e-9 --cores 20 *smc.gz

# Generate the plots
smc++ plot GLA.smc.jpeg -g 1 smc.out/model.final.json

```

### SMC for GLA with different geographical background

```bash
# Working directory
cd /data/zhaojiantao/pepper/SMC++/Accessions.diverse.background/GLA1

for i in $(seq -w 1 12); do
    for j in $(cat GLA.list) ; do
        # Extract the tested samples from all vcf
        vcftools --gzvcf /data/zhaojiantao/pepper/vcf.clean/Chr$i.clean.vcf.gz \
            --keep GLA.list \
            --maf 0.05 \
            --max-missing 0.8 \
            --recode --out Chr$i.clean
        # bgzip the vcf, not gzip
        bgzip Chr$i.clean.recode.vcf
        # Index the file
        tabix Chr$i.clean.recode.vcf.gz
        # Convert vcf to smc 
        smc++ vcf2smc --cores 60 --missing-cutoff 5000 -d $j $j \
            Chr$i.clean.recode.vcf.gz Chr$i.$j.smc.gz Chr$i \
            GLA:Grif_9138,Grif_9142,PI_406948,PI_438567,PI_593526,PI_593574,PI_631136,PI_555616,PI_661082,PI_674459
    done
done

# Run smc++
smc++ estimate -o smc.out 6.96e-9 --cores 20 *smc.gz

# Generate the plots
smc++ plot GLA.smc.jpeg -g 1 model.final.json

```



## GWAS analysis of fruit shape
### GWAS for _C. annuum_

```bash
# GLA.ANN working directory
cd /data/zhaojiantao/pepper/GWAS/GLA.ANN

for i in $(seq -w 1 12); do
    # Convert vcf to tped
    plink --vcf /data/zhaojiantao/pepper/vcf.clean/Chr$i.clean.vcf.gz \
        --keep GLA.ANN.list \
        --allow-extra-chr --chr Chr$i \
        --recode \
        --double-id \
        --geno 0.5 \
        --maf 0.05 \
        --out GLA.ANN.Chr$i
done

for i in $(seq -w 1 12); do
    sed 's/Chr0//g' GLA.ANN.Chr$i.map | sed 's/Chr//g' | \
    awk '{print $1,$1"_"$4,$3,$4}' | sed 's/ /\t/g' > GLA.ANN.Chr$i.map-1
done

rm *map *log *nosex

for i in $(seq -w 1 12); do
    mv GLA.ANN.Chr$i.map-1 GLA.ANN.Chr$i.map
done

# Make-bed
for i in $(seq -w 1 12); do
    plink --file GLA.ANN.Chr$i \
        --make-bed \
        --out GLA.ANN.Chr$i
done

# Calculate the effective number of SNPs
for i in $(seq -w 1 12); do
    java -jar /data/zhaojiantao/tools/GeneticTypeIError/gec.jar \
        --effect-number \
        --plink-binary GLA.ANN.Chr$i
    mv gec.sum GLA.ANN.Chr$i.gec.sum
done

# Calculate the effective number of SNPs of all chromosomes
cat *gec.sum | grep -v Effective_Number | awk '{ sum += $2 } END { print sum }' 

# Population structure
plink --vcf <(zcat ../../IQ-Tree504/Ca.S8.ffds.vcf.gz | grep -v Chr00) \
    --keep GLA.ANN.list \
    --allow-extra-chr \
    --recode \
    --double-id \
    --geno 0.2 \
    --maf 0.05 \
    --out GLA.ANN.ffds

sed 's/Chr0//g' GLA.ANN.ffds.map | sed 's/Chr//g' | awk '{print $1,$1"_"$4,$3,$4}' | sed 's/ /\t/g' > GLA.ANN.ffds.map-1
rm GLA.ANN.ffds.map
mv GLA.ANN.ffds.map-1 GLA.ANN.ffds.map

# Calculate PCA
plink --file GLA.ANN.ffds --pca --out GLA.ANN.ffds

# Remove unnecessary files
rm *log *nosex

# Prepare the PCA file as cofactor
awk '{print $1,$2,1,$3,$4,$5,$6,$7}' GLA.ANN.ffds.eigenvec > GLA.ANN.ffds.PCA.txt

# Transpose the ped file
for i in $(seq -w 1 12); do
    plink --file GLA.ANN.Chr$i --recode12 --output-missing-genotype 0 --transpose --out GLA.ANN.Chr$i
done

# Calculate kinship
plink --file GLA.ANN.ffds --recode12 --output-missing-genotype 0 --transpose --out GLA.ANN.ffds
emmax-kin -v -d 10 GLA.ANN.ffds

# EMMAX
for i in $(seq -w 1 12); do
for j in FruitLength FruitWidth FruitRatio ; do
    emmax -v -d 10 -t GLA.ANN.Chr$i -p $j.txt -k GLA.ANN.ffds.BN.kinf -c GLA.ANN.ffds.PCA.txt -o GLA.ANN.$j.Chr$i
done
done

# Prepare the output for manhattan plot
for i in FruitLength FruitWidth FruitRatio ; do
    cat GLA.ANN.$i.Chr*.ps | cut -f1 | sed 's/_/\t/g' | \
        paste - <(cat GLA.ANN.$i.Chr*.ps) | \
        awk '{print $1,$2,$3,$5}' | sed '1i\CHR BP SNP P' | sed 's/ /\t/g' \
        > GLA.ANN.$i.GWAS.txt
done

# Calculate the LD for the peak SNP
plink --file GLA.ANN.Chr02 \
    --r2 \
    --ld-snp 2_160951499 \
    --ld-window-kb 1000 \
    --ld-window 1050600 \
    --ld-window-r2 0.351 \
    --out 2_160951499

plink --file BAB.BAP.Chr10 \
    --r2 \
    --ld-snp 10_121586106 \
    --ld-window-kb 1000 \
    --ld-window 126500 \
    --ld-window-r2 0.366 \
    --out 10_121586106

plink --file BAB.BAP.Chr02 \
    --r2 \
    --ld-snp 2_172705768 \
    --ld-window-kb 1000 \
    --ld-window 126500 \
    --ld-window-r2 0.366 \
    --out 2_172705768

for i in BAP.long BAP.bell ; do
    vcftools --gzvcf ../vcf.clean/Chr10.clean.vcf.gz --keep $i.list --site-pi --out $i.Chr10
    perl PI_slide_win.pl -w 100000 -s 10000 $i.Chr10.sites.pi
done

# Step 15.2 Calculate windowed pi 
for i in Ca.annuum Ca.bac.baccatum Ca.bac.pendulum Ca.chacoense Ca.chinense Ca.frutescens Ca.glabriusculum Ca.glabriusculum_II Ca.pubescens ; do
    for j in $(seq -w 1 12); do
    perl PI_slide_win.pl -w 100000 -s 10000 $i.Chr$j.sites.pi
done
done

```

### Extract the genotypes for CaOvate

```bash
# Working directory
cd /data/zhaojiantao/pepper/Pi

# Extract the genotype and convert to plink ped
plink --vcf /data/zhaojiantao/pepper/vcf.clean/Chr10.clean.vcf.gz \
    --keep <(awk '{print $1,$1}' ../Pop.list/Ca.bac.baccatum.list | sed 's/ /\t/g' | cat - BAP.bell.list BAP.long.list) \
    --chr Chr10 --from-mb 120 --to-mb 122.6 \
    --recode --make-bed \
    --double-id \
    --maf 0.05 \
    --out Ca.baccatum.Chr10.region

# Convert to 012
plink --file Ca.baccatum.Chr10.region --recode12 transpose --out Ca.baccatum.Chr10.region

# Combine the genotypes
cut -f5- -d ' ' Ca.baccatum.Chr10.region.tped | sed 's! \([^ ]\+\)\( \|$\)!\1 !g' > Ca.baccatum.Chr10.region.tped-1

# Prepare the genotype
cut -f1 -d ' ' Ca.baccatum.Chr10.region.tfam > Ca.baccatum.Chr10.region.tfam-1

python -c "import sys; print('\n'.join(' '.join(c) for c in zip(*(l.split() for l in sys.stdin.readlines() if l.strip()))))" < Ca.baccatum.Chr10.region.tfam-1 > Ca.baccatum.Chr10.region.tfam.trans

cat Ca.baccatum.Chr10.region.tfam.trans Ca.baccatum.Chr10.region.tped-1 | sed 's/ /\t/g' > Ca.baccatum.Chr10.region.tped-2
python -c "import sys; print('\n'.join(' '.join(c) for c in zip(*(l.split() for l in sys.stdin.readlines() if l.strip()))))" < Ca.baccatum.Chr10.region.tped-2 > Ca.baccatum.Chr10.region.tped-3


# BAB geno
grep -Fwf ../Pop.list/Ca.bac.baccatum.list Ca.baccatum.Chr10.region.tped-3 > Ca.baccatum.Chr10.region.BAB.geno
grep -Fwf <(cut -f1 BAP.bell.list) Ca.baccatum.Chr10.region.tped-3 > Ca.baccatum.Chr10.region.bell.geno
grep -Fwf <(cut -f1 BAP.long.list) Ca.baccatum.Chr10.region.tped-3 > Ca.baccatum.Chr10.region.long.geno

cat Ca.baccatum.Chr10.region.*.geno | sed 's/ /\t/g' > Ca.baccatum.Chr10.region.txt

python -c "import sys; print('\n'.join(' '.join(c) for c in zip(*(l.split() for l in sys.stdin.readlines() if l.strip()))))" < Ca.baccatum.Chr10.region.txt > Ca.baccatum.Chr10.region.trans.txt


```


### Manhattan plot for GLA.ANN 

```r
library(qqman)
library(data.table)
source("manhattan.R")

GLA.ANN.GWAS.FruitLength <- fread("GLA.ANN.FruitLength.GWAS.txt", header = TRUE)
jpeg("GLA.ANN.FruitLength.jpeg", height = 8, width = 20, units="cm", res = 300)
manhattan(GLA.ANN.GWAS.FruitLength, chr="CHR", bp="BP", snp="SNP", p="P", 
          col = c("orange", "darkgreen"), cex=0.1, ylim=c(0,8),
          genomewideline = 6.29, suggestiveline = FALSE)
dev.off()

GLA.ANN.GWAS.FruitWidth <- fread("GLA.ANN.FruitWidth.GWAS.txt", header = TRUE)
jpeg("GLA.ANN.FruitWidth.jpeg", height = 8, width = 20, units="cm", res = 300)
manhattan(GLA.ANN.GWAS.FruitWidth, chr="CHR", bp="BP", snp="SNP", p="P", 
          col = c("orange", "darkgreen"), cex=0.1, ylim=c(0,8),
          genomewideline = 6.29, suggestiveline = FALSE)
dev.off()

GLA.ANN.GWAS.FruitRatio <- fread("GLA.ANN.FruitRatio.GWAS.txt", header = TRUE)
jpeg("GLA.ANN.FruitRatio.jpeg", height = 8, width = 15, units="cm", res = 300)
manhattan(GLA.ANN.GWAS.FruitRatio, chr="CHR", bp="BP", snp="SNP", p="P", 
          col = c("orange", "darkgreen"), cex=0.1, ylim=c(0,8),
          genomewideline = 6.29, suggestiveline = FALSE)
dev.off()

```

### Manhattan plot for BAB.BAP

```r
BAB.BAP.GWAS.FruitLength <- fread("BAB.BAP.FruitLength.GWAS.txt", header = TRUE)
jpeg("BAB.BAP.FruitLength.jpeg", height = 8, width = 20, units="cm", res = 300)
manhattan(BAB.BAP.GWAS.FruitLength, chr="CHR", bp="BP", snp="SNP", p="P", 
          col = c("purple", "darkorange"), cex=0.1, ylim=c(0,8),
          genomewideline = 6.39, suggestiveline = FALSE)
dev.off()

BAB.BAP.GWAS.FruitWidth <- fread("BAB.BAP.FruitWidth.GWAS.txt", header = TRUE)
jpeg("BAB.BAP.FruitWidth.jpeg", height = 8, width = 20, units="cm", res = 300)
manhattan(BAB.BAP.GWAS.FruitWidth, chr="CHR", bp="BP", snp="SNP", p="P", 
          col = c("purple", "darkorange"), cex=0.1, ylim=c(0,8),
          genomewideline = 6.39, suggestiveline = FALSE)
dev.off()

BAB.BAP.GWAS.FruitRatio <- fread("BAB.BAP.FruitRatio.GWAS.txt", header = TRUE)
jpeg("BAB.BAP.FruitRatio.jpeg", height = 8, width = 15, units="cm", res = 300)
manhattan(BAB.BAP.GWAS.FruitRatio, chr="CHR", bp="BP", snp="SNP", p="P", 
          col = c("purple", "darkorange"), cex=0.1, ylim=c(0,12),
          genomewideline = 6.39, suggestiveline = FALSE)
dev.off()

```

## Transcriptome data

```bash
# Working directory
cd /data/zhaojiantao/pepper/transcriptome

# Process the RNAseq data for fruit shape at different developing stages
# Save the above command line to Qualityclean.sh and run: nohup sh Qualityclean.sh > Qualityclean.log &
for i in $(cat fruitpeel.id); do
        # Trim adaptor and low quality
        trimmomatic PE $i/${i}_R1.fq.gz $i/${i}_R2.fq.gz \
                $i/${i}_R1P.fq $i/${i}_R1U.fq \
                $i/${i}_R2P.fq $i/${i}_R2U.fq \
                -threads 60 \
                ILLUMINACLIP:/data/zhaojiantao/tools/anaconda3/share/trimmomatic-0.39-2/adapters/TruSeq3-PE-2.fa:2:30:10:1:TRUE \
                SLIDINGWINDOW:4:20 LEADING:3 TRAILING:3 MINLEN:40
        # Remove polyA
        prinseq++ -threads 60 -VERBOSE 0 \
                -fastq $i/${i}_R1P.fq \
                -fastq2 $i/${i}_R2P.fq \
                -min_len 40 \
                -trim_tail_left 10 \
                -trim_tail_right 10 \
                -out_name $i/${i}
        # Remove rRNA
        bowtie -v 3 -k 1 -p 60 --al $i/${i}_R1P.rRNA.fq \
                /data/zhaojiantao/database/rRNA_region/rRNA_combined \
                -q $i/${i}_good_out_R1.fastq >/dev/null
        bowtie -v 3 -k 1 -p 60 --al $i/${i}_R2P.rRNA.fq \
                /data/zhaojiantao/database/rRNA_combined/rRNA_combined \
                -q $i/${i}_good_out_R2.fastq >/dev/null
        grep @ $i/${i}_R1P.rRNA.fq | cut -f1 -d ' ' | sed 's/@//g' > $i/${i}_R1P.rRNA.readsID 
        grep @ $i/${i}_R2P.rRNA.fq | cut -f1 -d ' ' | sed 's/@//g' > $i/${i}_R2P.rRNA.readsID
        cat $i/${i}_R1P.rRNA.readsID $i/${i}_R2P.rRNA.readsID | sort | uniq > $i/${i}_rRNA.sorted.readsID
        seqkit grep -vf $i/${i}_rRNA.sorted.readsID $i/${i}_good_out_R1.fastq -o $i/${i}_R1P.clean.fq.gz
        seqkit grep -vf $i/${i}_rRNA.sorted.readsID $i/${i}_good_out_R2.fastq -o $i/${i}_R2P.clean.fq.gz
        # Remove unnecessary files
        rm $i/${i}_R*U.fq $i/${i}_R*P.fq $i/${i}_R*P.rRNA.fq $i/*.readsID $i/*bad* $i/*good* $i/*single*
        # Index the reference genome. Only run this for one time.
        # hisat2-build -f /data/zhaojiantao/pepper/References/S8/Zhangshugang.fasta \
        #        /data/zhaojiantao/pepper/References/S8/Zhangshugang.fasta
        #Map the reads
        hisat2 --rna-strandness RF -p 32 \
            --summary-file $i/$i.cleaned.summary \
            -x /data/zhaojiantao/pepper/References/S8/Zhangshugang.fasta \
            -1 $i/${i}_R1P.clean.fq.gz \
            -2 $i/${i}_R2P.clean.fq.gz \
            -S $i/$i.cleaned.sam
        # Convert sam to bam
        samtools view -S -b $i/$i.cleaned.sam > $i/$i.cleaned.bam
        # Remove sam file
        rm $i/$i.cleaned.sam
done

# Prepare the gene bed file 
cat /data/zhaojiantao/pepper/References/S8/Zhangshugang.gff | \
    awk '$3=="CDS" && $9~/^ID/{sub(/;.*/,"",$9);L[substr($9,4)]+=$5-$4+1}END{for(i in L){print i,L[i]}}' | \
    sed 's/cds-//g' | sort -k1 > /data/zhaojiantao/pepper/References/S8/Zhangshugang.gff.trans.length

sed 's/;/\t/g' /data/zhaojiantao/pepper/References/S8/Zhangshugang.gff | \
    awk '$3=="mRNA" {print $1,$4,$5,$7,$9}' | sed 's/ID=//g' | sort -k5 | \
    paste - /data/zhaojiantao/pepper/References/S8/Zhangshugang.gff.trans.length | \
    sort -t ' ' -k1,1 -k5,5 |  awk '{print $1,$2,$3,$5,$7,$4}' | sed 's/ /\t/g' \
    > /data/zhaojiantao/pepper/References/S8/Zhangshugang.gff.bed

rm /data/zhaojiantao/pepper/References/S8/Zhangshugang.gff.trans.length

# Get the longest transcripts for each gene using [AGAT](https://agat.readthedocs.io/en/latest/index.html)
# agat_sp_keep_longest_isoform.pl: The script aims to filter isoforms when present. For a locus: - when all isoforms have CDS we keep the one with the longest CDS. - when some isoforms have CDS some others not, we keep the one with the longest CDS. - when none of the isoforms have CDS, we keep the one with the longest concatenated exons.
# Source: https://agat.readthedocs.io/en/latest/tools/agat_sp_keep_longest_isoform.html
agat_sp_keep_longest_isoform.pl -gff Ca.S8.gff -o Ca.S8.longestisoform.gff

# Count the reads
featureCounts -T 20 -s 1 -p \
    -a /data/zhaojiantao/pepper/References/S8/Ca.S8.longestisoform.gff \
    -t gene -g ID */*.bam -O \
    -o raw.count 

# Generate raw count
mRNAtool.pl -t count -s PS -l fr-firststrand \
    -f /data/zhaojiantao/pepper/References/S8/Zhangshugang.gff.bed */*.bam 

# Prepare the mapped reads summary file
# Total reads
for i in $(cat fruitpeel.id); do sed -n '1p' $i/$i.cleaned.summary | awk '{print $1}' ; done

# Pair-end concordantly exactly 1 time
for i in $(cat fruitpeel.id); do sed -n '4p' $i/$i.cleaned.summary | awk '{print $1}' ; done
# Pair-end concordantly >1 times
for i in $(cat fruitpeel.id); do sed -n '5p' $i/$i.cleaned.summary | awk '{print $1}' ; done
# Single concordantly exactly 1 time
for i in $(cat fruitpeel.id); do sed -n '13p' $i/$i.cleaned.summary | awk '{print $1}' ; done
# Single concordantly exactly >1 time
for i in $(cat fruitpeel.id); do sed -n '14p' $i/$i.cleaned.summary | awk '{print $1}' ; done

# Calculate the reads count
mRNAtool.pl -t count -s PS -l fr-firststrand \
    -f /data/zhaojiantao/pepper/References/S8/Zhangshugang.gff.bed */*.bam
    
# Normalize the raw count
mRNAtool.pl -t norm -f /data/zhaojiantao/pepper/References/S8/Zhangshugang.gff.bed -m 1 -u mapped_reads exp_raw_count.txt > rpkm.txt

# Step 6.3 Calculate correlation
mRNAtool.pl -t corre rpkm.txt > rpkm_corre.txt

```

## Personal scripts

### 4DTV_scan.pl

```perl
#! /usr/bin/env perl
use strict;
use warnings;
use Bio::SeqIO;
#use Getopt::Std;
use Getopt::Long;
use File::Basename;
# Author: Xin Wang   20131020
my ($genome_file,$gff_file,$fold);
GetOptions(
				"g|genome=s" => \$genome_file,
				"f|gff=s" => \$gff_file,
                "o|output=s" => \$fold,
        );
        
die(qq/Usage: Get the 4DTV site from genome \n
Options: -g|genome      The genome fasta\n
         -f|gff         The GFF file\n
         -o|output      The output file\n
/) unless ( $genome_file && $gff_file && $fold);

#=@ARGV;

#die "Usage: $0 <genome file> <gff file> <4 fold degeneration sites output>\n" if(@ARGV < 3);

my %gff;
open(I,"< $gff_file");
my $no=0;
while(<I>){
    chomp;
    next if(/^#/);
    my @a=split(/\s+/);
    next unless($a[2] eq "CDS");
    $no++;
    my ($chr,$start,$end,$strand,$phase,$name)=($a[0],$a[3],$a[4],$a[6],$a[7],$a[8]);
    # $chr=~s/chr//g;
    $name=~/Parent=([^;]+)/;
    $name=$1;
    $gff{$chr}{$name}{$no}{start}=$start;
    $gff{$chr}{$name}{$no}{end}=$end;
    $gff{$chr}{$name}{$no}{strand}=$strand;
    $gff{$chr}{$name}{$no}{phase}=$phase;
}
close I;

print STDERR "GFF reading complete!\n";

my $fa=Bio::SeqIO->new(-file=>$genome_file,-format=>'fasta');

my $control=0;
open(O,"> $fold");
while(my $seq=$fa->next_seq){
    my $chr=$seq->id;
    my $seq=$seq->seq;
    next unless(exists $gff{$chr});

    foreach my $name(keys %{$gff{$chr}}){
        my $strand="NA";
        my $line="";
        # die "debug" if($control++>100);
        foreach my $no(sort { $gff{$chr}{$name}{$a}{start} <=> $gff{$chr}{$name}{$b}{start} } keys %{$gff{$chr}{$name}}){
            if($strand eq "NA"){
	$strand=$gff{$chr}{$name}{$no}{strand};
            }
            my $start=$gff{$chr}{$name}{$no}{start};
            my $end=$gff{$chr}{$name}{$no}{end};
            my $len=$end-$start+1;
            my $subline=substr($seq,$start-1,$len);
            $line.=$subline;
        }

        my @fold_sites_location_in_genome;

        if($strand eq "+"){
            my $pep = &translate_nucl($line);

            my @fold_sites_location_in_cds = &get_location_in_cds($line);

            next if(scalar(@fold_sites_location_in_cds) == 0);

            my %corresponding_postion;
            my $accumulating_length=0;
            foreach my $no(sort { $gff{$chr}{$name}{$a}{start} <=> $gff{$chr}{$name}{$b}{start} } keys %{$gff{$chr}{$name}}){
	my $start=$gff{$chr}{$name}{$no}{start};
	my $end=$gff{$chr}{$name}{$no}{end};
	my $len=$end-$start+1;

	my $real_start=$start;
	my $real_end  =$end;
	my $cds_start =$accumulating_length+1;
	my $cds_end   =$accumulating_length+$len;

	my $real_pos=$real_start;
	for(my $i=$cds_start;$i<=$cds_end;$i++){
	    my $cds_pos=$i;
	    $corresponding_postion{$cds_pos}=$real_pos;
	    $real_pos++;
	}
	$accumulating_length+=$len;
            }
            foreach my $cds_pos(@fold_sites_location_in_cds){
	my $real_pos=$corresponding_postion{$cds_pos};
	push @fold_sites_location_in_genome,$real_pos;
            }
        }
        elsif($strand eq "-"){
            $line=reverse($line);
            $line=~tr/ATCGatcg/TAGCtagc/;
            my $pep = &translate_nucl($line);

            my @fold_sites_location_in_cds = &get_location_in_cds($line);

            next if(scalar(@fold_sites_location_in_cds) == 0);

            my %corresponding_postion;
            my $accumulating_length=0;
            foreach my $no(sort { $gff{$chr}{$name}{$b}{start} <=> $gff{$chr}{$name}{$a}{start} } keys %{$gff{$chr}{$name}}){
	my $start=$gff{$chr}{$name}{$no}{start};
	my $end=$gff{$chr}{$name}{$no}{end};
	my $len=$end-$start+1;

	my $real_start=$start;
	my $real_end  =$end;
	my $cds_start =$accumulating_length+1;
	my $cds_end   =$accumulating_length+$len;

	my $real_pos=$real_end;
	for(my $i=$cds_start;$i<=$cds_end;$i++){
	    my $cds_pos=$i;
	    $corresponding_postion{$cds_pos}=$real_pos;
	    $real_pos--;
	}
	$accumulating_length+=$len;
            }
            foreach my $cds_pos(@fold_sites_location_in_cds){
	my $real_pos=$corresponding_postion{$cds_pos};
	push @fold_sites_location_in_genome,$real_pos;
            }
        }

        foreach my $pos(@fold_sites_location_in_genome){
            print O "$chr\t$pos\t$name\n";
        }

        print STDERR "$chr\t$name\n";
    }
}
close O;

sub get_location_in_cds{
    my $line=shift;
    my @location;
    my $len=length($line);

    if($len == 0){
        return(@location);
    }

    my @bases=split("",$line);
    my @ATCG=("A","T","C","G");

    for(my $i=0;$i<$len;$i+=3){
        next if(!$bases[$i] or !$bases[$i+1] or !$bases[$i+2]);
        my @codon_bases=($bases[$i],$bases[$i+1],$bases[$i+2]);
        my $codon_bases=join "",@codon_bases;
        $codon_bases=~tr/atcg/ATCG/;
        next if($codon_bases=~/[^ATCG]/);
        my $codon_pep  = &translate_nucl($codon_bases);

        for(my $j=0;$j<3;$j++){
            my $light = 1;
            my @new_codon_bases = @codon_bases;
            foreach my $bases(@ATCG){
	$new_codon_bases[$j]=$bases;
	my $new_codon_bases = join "",@new_codon_bases;
	my $new_codon_pep   = &translate_nucl($new_codon_bases);
	if($new_codon_pep ne $codon_pep){
	    $light=0;
	    last;
	}
            }
            if($light == 1){
	my $location_in_cds=$i+$j+1;
	push @location,$location_in_cds;
            }
        }
    }
    my $number=@location;
    print STDERR " $number four-fold-degenerate-sites in this gene\n";
    return(@location);
}

sub translate_nucl{
    my $seq=shift;
    my $seq_obj=Bio::Seq->new(-seq=>$seq,-alphabet=>'dna');
    my $pro=$seq_obj->translate;
    $pro=$pro->seq;
    return($pro);
}
```

### assign_genetic.pl

```perl
#!/usr/bin/perl -w
use strict;

my ($genetic_map,$in) = @ARGV;

my %map;

open IN,"<",$genetic_map;
while(<IN>){
	chomp;
	my @line = split;
	my @infor = split/\_/,$line[0];
	my $chr = $infor[0];
	my $pos = $infor[1];
	my $genetic_distance = $line[2];
	$map{$chr}{$pos} = $genetic_distance;
}
close IN;

my %snp;
open IN,"<",$in;
while(<IN>){
	chomp;
	my @line = split;
	my @infor = split/\_/,$line[0];
	my $chr = $infor[0];
	my $pos = $infor[1];
	$snp{$chr}{$pos}{$line[0]} =0;
}

#my @test = ("SL2.50ch01");
for my $chr(keys %snp){
	my @pos = sort {$a <=> $b} keys %{$snp{$chr}};
	my $pos_length = @pos;
	my @map = sort {$a <=> $b} keys %{$map{$chr}};
	my $map_length = @map;
	open OUT,">","$chr.snp.genetic.txt";
	for my $pos(@pos){
		my $distance = "NA";
		if($pos <=$map[0]){
			$distance = $pos*$map{$chr}{$map[0]}/($map[0]);
		}
		elsif($pos >=$map[-1]){
			my $bin =(abs($map{$chr}{$map[-1]} - $map{$chr}{$map[-2]}))/($map[-1] - $map[-2]+1);
			$distance = $map{$chr}{$map[-1]} + ($pos -$map[-1])*$bin;
		}
		for my $i(0..($map_length-2)){
			if(($pos>$map[$i]) && ($pos<=$map[$i+1])){
				my $bin =(abs($map{$chr}{$map[$i+1]} - $map{$chr}{$map[$i]}))/($map[$i+1] - $map[$i]+1);
				$distance = $map{$chr}{$map[$i]} + ($pos -$map[$i])*$bin;
			}
		}
		my @id = keys %{$snp{$chr}{$pos}};
		$distance = sprintf("%.3f",$distance);
		print OUT"$id[0]\t$pos\t$distance\n";
	}
	close OUT;
}

```

### PI_slide_win.pl

```perl
#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Std;
use Set::IntervalTree;
$|=1;

our ($opt_i, $opt_p, $opt_s, $opt_d, $opt_w, $opt_o);
getopt("s:i:p:d:w:o:");

&usage if (@ARGV < 1);


print "Starting $0...\n";

# default parameters
$opt_w = $opt_w || 100000;
$opt_s = $opt_s || 10000;


for my $infile(@ARGV){
	my $time1 = time();
	open IN,"<",$infile;
	my %pi_per_site;
	my %genome_length;
	while(<IN>){
		chomp;
		next if(/^CHROM/);
		my @line = split/\t/,$_;
		my $chr = $line[0];
		my $pos = $line[1];
		my $pi = $line[2];
		if($pi eq "-nan"){
			$pi =0;
		}
		$pi_per_site{$chr}{$pos} =$pi;
		#$genome_length{$chr}{$pos} =0;
	}
	close IN;
	our %window_hash = ();	#####change the my to our to avoid the Variable "%window_hash" will not stay shared 
	foreach my $chr (sort keys %pi_per_site){
		my @pos = sort {$a <=> $b} keys %{$pi_per_site{$chr}};
		my $chr_len = $pos[-1];
	#	print "\tprocessing $opt_p$chr, length: $chr_len, ";
		# get the window number, window start, window end
		my $window_num = 0;
		if ($chr_len <= $opt_w) {               #windowsize greater than chro_length
			$window_num = 1;
			($window_hash{$chr}{$window_num}{start}, $window_hash{$chr}{$window_num}{end}) = (1, $chr_len);
		}else{                                  #windowsize less than chro_length
			$window_num = &split_genome($chr, $chr_len);
		}
	#	print "window_number: $window_num.\n";
	}
	##generate interval
	print "Step1: calculating PI per window for $infile\n";
	my %average_pi;	
	for my $chr(sort {$a cmp $b} keys %window_hash){
		my $chr_tree = $chr;
		$chr_tree = Set::IntervalTree->new;
		for my $w_num(sort {$a <=> $b} keys %{$window_hash{$chr}}){
			my $start = $window_hash{$chr}{$w_num}{'start'};
			my $end = $window_hash{$chr}{$w_num}{'end'};
			#print "$chr\t$w_num\t$window_hash{$chr}{$w_num}{'start'}\t$window_hash{$chr}{$w_num}{'end'}\t$N_num\n";
			$chr_tree->insert("$chr"."_"."$w_num"."_"."$start"."_"."$end", $start, $end);
			$average_pi{$chr}{$w_num}{$start}{$end}[0]=0;   # variation number
			$average_pi{$chr}{$w_num}{$start}{$end}[1] =0; #PI sum
		}
		for my $pos(sort {$a <=> $b} keys %{$pi_per_site{$chr}}){
			my $interval = $chr_tree->fetch($pos-0.5,$pos+0.5);
			my @interval_array = @$interval;
			for my $interval_infor(@interval_array){
				my @infor = split/\_/,$interval_infor;
				$average_pi{$infor[0]}{$infor[1]}{$infor[2]}{$infor[3]}[0]++; # variation number
				$average_pi{$infor[0]}{$infor[1]}{$infor[2]}{$infor[3]}[1]=$average_pi{$infor[0]}{$infor[1]}{$infor[2]}{$infor[3]}[1] + $pi_per_site{$chr}{$pos}; #Sum PI			
			}
		}	
	}
	
	####output average PI for one window
	print "Step2: outputing PI per window for $infile\n";
	my $whole_pi =0;
	my $whole_win_num = 0;
	#print "\#Chr\tWindow number\tStart\tEnd\tBin_Size\tNon_N_Bin_Size\tPI_Bin\tPI_Non_N_Bin\n";
	open OUT,">","$infile.w$opt_w.s$opt_s.txt";
	print OUT"\#Chr\tWindow number\tStart\tEnd\tBin_Size\tN_var\tPI\n";
	for my $chr(sort {$a cmp $b} keys %average_pi){
		#print OUT"\#Chr\tWindow number\tStart\tEnd\tBin_Size\tNon_N_Bin_Size\tPI_Bin\tPI_Non_N_Bin\n";
		for my $win_num(sort {$a <=> $b} keys %{$average_pi{$chr}}){
			print OUT"$chr\t$win_num\t";
			$whole_win_num++;
			for my $start(sort {$a <=> $b} keys %{$average_pi{$chr}{$win_num}}){
				for my $end(sort {$a <=> $b} keys %{$average_pi{$chr}{$win_num}{$start}}){
					my $bin_size = $end -$start +1;
					my $average_pi = $average_pi{$chr}{$win_num}{$start}{$end}[1]/$bin_size;
					$whole_pi = $whole_pi + $average_pi;
					#print OUT"$start\t$end\t$bin_size\t$non_n_bin_size\t$average_pi{$chr}{$win_num}{$start}{$end}[0]\t$average_pi{$chr}{$win_num}{$start}{$end}[1]\n";;
					print OUT"$start\t$end\t$bin_size\t$average_pi{$chr}{$win_num}{$start}{$end}[0]\t$average_pi\n";;		
				}
			}
		}
	}
	close OUT;
	my $average_whole_pi = $whole_pi/$whole_win_num;
	print "The whole genome PI for $infile is $average_whole_pi\n\n";	
	my $time2 = time();
	my $time  = sprintf "%.4f", ( $time2 - $time1 ) / 60;
	print "$0 finished! $time mins elapsed for $infile.\n";
	
	sub split_genome {
		my ($chro_id, $chro_length) = @_;
		my $win_num =0;
		for (my $start=1;$start < ($chro_length - $opt_w + $opt_s);$start=$start + $opt_s){
			my $end = $start + $opt_w - 1;
			$win_num++;
			if($end < $chro_length){
				$window_hash{$chro_id}{$win_num}{'start'}=$start;
				$window_hash{$chro_id}{$win_num}{'end'}=$end;
			}
			else{
				$end = $chro_length;
				$window_hash{$chro_id}{$win_num}{'start'}=$start;
				$window_hash{$chro_id}{$win_num}{'end'}=$end;
				#print "$size\n";
			}
		}
	}
}	

sub usage {
print"Usage: $0 -OPTIONS VALUES infile1 infile2 infile3 ......
This script is used to caculate the average PI using sliding window.
Options:
     -w  NO   windowsize to caculate the average PI, default = 100000
     -s  NO   step size, should be greater than 0kb and less than windowsize, default = 10000 
     \n"
}
```

### merge.pl
This script is used to merge the top 4% XP-CLR scores with distances within 100000

```perl
#!/usr/bin/perl
while(<>) {

        chomp;

        next if (/SNP/);

        @a = split "\t";

        if ($a[2] - $end > 100000 || $a[1] ne $chr || eof() ) { print "$chr\t$start\t$end\t$max\n"; $start = $a[2]; $end = $a[3]; $max = $a[4]; }

        else { $end = $a[3]; if ($a[4] > $max) { $max = $a[4]; } }

        $chr = $a[1];

}

```
