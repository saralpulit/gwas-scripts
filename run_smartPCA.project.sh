#!/bin/sh
#$ -cwd
#$ -S /bin/bash
#$ -l h_rt=24:00:00

### To run: ./run_smartPCA.project.sh $PED_ROOT
### $PED_ROOT is the rootname of your plink files

### Pointers to programs
CONVERTF_EXEC=/home/software/EIG4.2/bin/convertf
SMARTPCA_EXEC=/home/software/EIG4.2/bin/smartpca
EIGENSTRAT_EXEC=/home/software/EIG4.2/bin/eigenstrat
PED_ROOT=$1

### MUST put smartpca bin directory in path for smartpca.perl to work
PATH="/home/software/eig4.2m/EIG4.2/bin:$PATH"

### Convert to eigenstrat format using convertf
### note that this data should contain both your GWAS data and reference data (e.g., HapMap3, 1000 Genomes)
### the pedind file contains the first 5 columns of the PED file, and then a label for a "test" sample (your GWAS data) or a "ref" sample (your reference data)
echo "genotypename:    $PED_ROOT.ped" > $PED_ROOT.convertf.par
echo "snpname:         $PED_ROOT.map" >> $PED_ROOT.convertf.par
echo "indivname:       $PED_ROOT.pedind" >> $PED_ROOT.convertf.par
echo "outputformat:    EIGENSTRAT" >> $PED_ROOT.convertf.par
echo "genotypeoutname: $PED_ROOT.geno" >> $PED_ROOT.convertf.par
echo "snpoutname:      $PED_ROOT.snp" >> $PED_ROOT.convertf.par
echo "indivoutname:    $PED_ROOT.ind" >> $PED_ROOT.convertf.par
echo "familynames:     NO" >> $PED_ROOT.convertf.par

$CONVERTF_EXEC -p $PED_ROOT.convertf.par

### Now, calculate eigenvectors
### Note that your reference samples will be labeled "ref" and your GWAS samples will be labeled "test"

echo "ref"             > $PED_ROOT.POPlist
echo "genotypename:    $PED_ROOT.geno" > $PED_ROOT.smartpca.par
echo "snpname:         $PED_ROOT.snp" >> $PED_ROOT.smartpca.par
echo "indivname:       $PED_ROOT.pedind" >> $PED_ROOT.smartpca.par
echo "evecoutname:     $PED_ROOT.evec" >> $PED_ROOT.smartpca.par
echo "evaloutname:     $PED_ROOT.eval" >> $PED_ROOT.smartpca.par
echo "altnormstyle:    NO" >> $PED_ROOT.smartpca.par
echo "numoutevec:      100" >> $PED_ROOT.smartpca.par ## This can be reset if you want more PCs (up to N-1)
echo "numoutlieriter:  0"  >> $PED_ROOT.smartpca.par
echo "numthreads:      10" >> $PED_ROOT.smartpca.par
echo "poplistname:     $PED_ROOT.POPlist" >> $PED_ROOT.smartpca.par

$SMARTPCA_EXEC -p $PED_ROOT.smartpca.par > $PED_ROOT.smartpca.log
