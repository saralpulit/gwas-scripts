## gwas-scripts
A repository containing various scripts useful for performing quality control on data from genome-wide association studies and visualizing results

### Scripts contained in this repository:

#### (1) qqplot.R
A script for generating a standard quantile-quantile (QQ) plot and calculating lambda for your data. The script takes in an input file containing the column "P" (with your GWAS p-values). 

The script is run as follows:

```Rscript qqplot.R inputfile stat_type output_file color xmamx```  

   *stat_type: P (for p-values), Z (for z-scores), or CHISQ (for chi-square values)*  
   *color: color you want your points plotted in*  
   *xmax: maximum value for the x-axis (and y-axis, as the plot is square)*  
   
#### (2) qq_by_maf.R
A script for generating a QQ plot that plots data stratified by allele frequency. The input data must contain the following columns: P (p-values from your GWAS), FRQ (the frequency of the tested allele from your GWAS). 

The script runs as follows:

```Rscript qq_by_maf.R inputfile outputfile```

#### (3) qq_by_info.R
A script for generating a QQ plot that plots data stratified by imputation quality. Does exactly the same as the QQ by maf script, but stratifying by imputation quality (INFO) score (typically a number between 0 and 1). The input data must contain the following columns: P (p-values from your GWAS), INFO (imputation quality score).

The script runs as follows:

```Rscript qq_by_info.R inputfile outputfile```

#### (4) manhattan.R
A script for plotting a Manhattan plot. Note that this depends on having the "qqman" library (written by Stephen Turner; please see his excellent GitHub for more information https://github.com/stephenturner/qqman). The input data must contain the following columns and be in this order: SNP (rsID), CHR (chromosome), BP (base position), P (p-value)

#### (5) qq.three_plots.R
Combines scripts 1-3 into a single script to produce three QQ plots at once. The columns of your data must include: P, FRQ, and INFO.

The script runs as follows:

```Rscript qqplot.R inputfile stat_type output_file color xmamx```  

   *stat_type: P (for p-values), Z (for z-scores), or CHISQ (for chi-square values)*  
   *color: color you want your points plotted in*  
   *xmax: maximum value for the x-axis (and y-axis, as the plot is square)*  
   
#### (6) lambda.R
A quick script to calculate lambda in your data. Your data must have column "P" containing p-values.

The script runs as follows:

```Rscript lambda.R inputfile stat_type```

   *stat_type: P (for p-values), Z (for z-scores), or CHISQ (for chi-square values)*
   
#### (7) prep.PCA.sh
A script for preparing (Plink-formatted) data to run principal component analysis.

#### (8) run_smartPCA.project.sh
A script that makes a call to EIGENSTRAT to plot  
   (a) your GWAS samples, projected onto  
   (b) a reference dataset (for example, HapMap 3)  
   
Input data should be:  
   (a) A Plink-formatted PED file containing the genotypes from your GWAS samples and the samples from your reference data.
   (b) A Plink-formatted MAP file, corresponding with the PED file, containing marker information.
   (c) A PEDIND file. Contains the first 5 columns of the PED file. The 6th column indicates "test" (i.e., a sample from your GWAS data) or "ref" (i.e., a sample from your reference data.

The script runs as follows:

```./run_smartPCA.project.sh file-root-name```

where file-root-name is the root name of your PED, MAP and PEDIND files.

#### (9) run_smartPCA.sh
Exactly as above to run principal component analysis in EIGENSTRAT, but only performing the analysis within your own GWAS data samples. You therefore only need a PED and MAP file.

The script runs as follows:

```./run_smartPCA.sh file-root-name```

where file-root-name is the root name of your PED and MAP files.


   
