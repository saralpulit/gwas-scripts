#!/bin/bash
#$ -cwd

$1=mydata
$2=hapmap

plink=/path/to/plink1.9

### An example script for prepping PLINK-formatted data for principal component analysis using a reference dataset
### For this script, assume that your reference data is HapMap 3 and you pass the root names of the HapMap files as the second command line argument

### First, generate a high-quality of set of SNPs from your GWAS data on which you want to perform PCA
  
  ## Drop all of the SNPs in the region of the LCT locus (hg19 coordinates)
  awk ' $1==2 && $4 > 129883530 && $4 < 140283530 { print $2 } ' $mydata.bim > $mydata.pca.EXCLUDE.snps
  
  ## Drop all of the SNPs in the major histocompatibility complex (hg19 coordinates)
  awk ' $1==6 && $4 > 24092021 && $4 < 38892022 { print $2 } ' $mydata.bim >> $mydata.pca.EXCLUDE.snps
  
  ## Drop all of the SNPs in the inverted regions on chromosomes 8 and 17
  awk ' $1==8 && $4 > 6612592 && $4 < 13455629 { print $2 } ' $mydata.bim >> $mydata.pca.EXCLUDE.snps
  awk ' $1==17 && $4 > 40546474 && $4 < 44644684 { print $2 } ' $mydata.bim >> $mydata.pca.EXCLUDE.snps
  
  ## Drop all non-autosomal SNPs
  awk ' $1>22 { print $2 } ' $mydata.bim >> $mydata.pca.EXCLUDE.snps

### Select common SNPs (maf > 0.10) and low-missingness SNPs (missingness < 0.1%), drop the SNPs from above, and LD prune your data at an r2 = 0.2
$plink --bfile $mydata \
  --bfile $mydata \
  --allow-no-sex \
  --geno 0.001 \
  --maf 0.10 \
  --exclude $mydata.pca.EXCLUDE.snps \
  --indep-pairwise 50 5 0.2 \
  --out $mydata
  
## Extract the "prune.in" SNPs from your data; these are the SNPs we want for PCA
$plink --bfile $mydata \
  --allow-no-sex \
  --extract $mydata.prune.in \
  --make-bed \
  --out $mydata.pca

### Now, merge your gwas data with the reference data
$plink --bfile $mydata.pca \
  --allow-no-sex \
  --bmerge $hapmap.bed $hapmap.bim $hapmap.fam \
  --make-bed \
  --out $mydata.$hapmap.pca.merge
  
### Merging may trigger an issue with alignment of alleles
### Plink will automatically generate a list of SNPs that need to be flipped, called *.missnp
### use this block of code to flip alleles in your gwas data if that happens, and attempt merging again
$plink --bfile $mydata.pca \
  --allow-no-sex \
  --flip $mydata.$hapmap.pca.merge.missnp \
  --make-bed \
  --out $mydata.pca
  
### If flipping alleles does not resolve the issue, you can just drop the problematic SNPs (if it's a small number. a large number of problematic SNPs might signal a data formatting problem)
$plink --bfile $mydata.pca \
  --allow-no-sex \
  --exclude $mydata.$hapmap.pca.merge.missnp \
  --make-bed \
  --out $mydata.pca
  
## Lastly, once the merge is done, make sure you only keep the low missingness SNPs (i.e., present in your GWAS data and in your reference data)
$plink --bfile $mydata.$hapmap.pca.merge \
  --allow-no-sex \
  --geno 0.01 \
  --make-bed \
  --out $mydata.$hapmap.pca.merge
  
  
