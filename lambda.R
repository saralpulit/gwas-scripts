## calculate genomic inflation (lambda) in your data set
## Usage: Rscript inputfile stat_type
## stat_type is P, Z, or CHISQ

input=commandArgs(trailingOnly=T)[1]
stat_type=commandArgs(trailingOnly=T)[2] ## P, Z, CHISQ

## Reads data
## Header should have column "P"

S <- read.table(input,header=T)
if (stat_type == "Z")
z=S$P

if (stat_type == "CHISQ")
z=sqrt(S$P)

if (stat_type == "P")
z=qnorm(S$P/2)

## calculates lambda
lambda = round(median(z^2,na.rm=T)/qchisq(0.5,df=1),3)
print(lambda)
