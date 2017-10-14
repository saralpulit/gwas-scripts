#!/bin/sh
#$ -cwd
#$ -S /bin/bash
#$ -l h_rt=24:00:00

### To run: ./run_smartPCA.sh $PED_ROOT
### $PED_ROOT is the rootname of your plink files
### Runs PCA within your own GWAS samples (e.g., to generate PCs to be used as covariates in your GWAS)

### Pointers to programs
CONVERTF_EXEC=/home/software/EIG4.2/bin/convertf
SMARTPCA_EXEC=/home/software/EIG4.2/bin/smartpca
EIGENSTRAT_EXEC=/home/software/EIG4.2/bin/eigenstrat
PED_ROOT=$1

### MUST put smartpca bin directory in path for smartpca.perl to work
PATH="/home/software/eig4.2m/EIG4.2/bin:$PATH"

### Convert to eigenstrat format using convertf

echo "genotypename:    $PED_ROOT.ped" > $PED_ROOT.convertf.par
echo "snpname:         $PED_ROOT.map" >> $PED_ROOT.convertf.par
echo "indivname:       $PED_ROOT.ped" >> $PED_ROOT.convertf.par
echo "outputformat:    EIGENSTRAT" >> $PED_ROOT.convertf.par
echo "genotypeoutname: $PED_ROOT.geno" >> $PED_ROOT.convertf.par
echo "snpoutname:      $PED_ROOT.snp" >> $PED_ROOT.convertf.par
echo "indivoutname:    $PED_ROOT.ind" >> $PED_ROOT.convertf.par
echo "familynames:     NO" >> $PED_ROOT.convertf.par

$CONVERTF_EXEC -p $PED_ROOT.convertf.par

### Now, calculate eigenvectors
### Note that your reference samples will be labeled "ref" and your GWAS samples will be labeled "test"

echo "genotypename:    $PED_ROOT.geno" > $PED_ROOT.smartpca.par
echo "snpname:         $PED_ROOT.snp" >> $PED_ROOT.smartpca.par
echo "indivname:       $PED_ROOT.ped" >> $PED_ROOT.smartpca.par
echo "evecoutname:     $PED_ROOT.evec" >> $PED_ROOT.smartpca.par
echo "evaloutname:     $PED_ROOT.eval" >> $PED_ROOT.smartpca.par
echo "altnormstyle:    NO" >> $PED_ROOT.smartpca.par
echo "numoutevec:      100" >> $PED_ROOT.smartpca.par ## This can be reset if you want more PCs (up to N-1)
echo "numoutlieriter:  0"  >> $PED_ROOT.smartpca.par
echo "numthreads:      10" >> $PED_ROOT.smartpca.par

$SMARTPCA_EXEC -p $PED_ROOT.smartpca.par > $PED_ROOT.smartpca.log
