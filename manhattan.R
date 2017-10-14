## To run: Rscript manhattan.R inputfile outputfile color1 color2
## color1 and color2 are the alternating color for the chromosomes

rm(list=ls())
library(qqman)

### Columns need to be SNP, CHR, BP, P (in that order)

input=commandArgs(trailingOnly=T)[1]
output=commandArgs(trailingOnly=T)[2]
color1=commandArgs(trailingOnly=T)[3]
color2=commandArgs(trailingOnly=T)[4]

S <- read.table(input,header=T)

### Manhattan plot
### Change the header of the file

### Update the column names to work with manhattan script
colnames(S) <- c("SNP","CHR","BP","P")

### Generate a QQ plot
png(file=paste(output,".QQ.png",sep=""),height=480,width=480)
qq(S$P)
dev.off()

### Drop everything with p > 0.5
S <- subset(S,S$P<0.5)

### Plot the manhattan plot
chr <- as.numeric(max(S$CHR))

png(file=paste(output,".manhattan.png",sep=""),height=480,width=960)
manhattan(S, pt.bg=c(color1,color2),pt.col=c(color1,color2), pch=21, genomewideline=T, suggestiveline=F, pt.cex=1.1,limitchromosomes=1:chr)
lines(x=c(0,10000000000),y=c(7.30103,7.30103),col="grey50",lty="dashed",lwd=0.5)

dev.off()
