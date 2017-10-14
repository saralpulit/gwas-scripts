### Generates a QQ plot from GWAS data
### Input file must contain column "P"
### To run: Rscript qqplot.R inputfile stat_type outputfile color xmax

## stat_type: P, Z, CHISQ
## xmax: dimensions of the plot

rm(list=ls())

input=commandArgs(trailingOnly=T)[1]
stat_type=commandArgs(trailingOnly=T)[2]
output=commandArgs(trailingOnly=T)[3]
color=commandArgs(trailingOnly=T)[4]
xmax=commandArgs(trailingOnly=T)[5]

xmax <- as.numeric(xmax)
print(xmax)

## Plot function ##
   plotQQ <- function(z,color){
   p <- 2*pnorm(-abs(z))
   p <- sort(p)
   expected <- c(1:length(p))
   lobs <- -(log10(p))
   lexp <- -(log10(expected / (length(expected)+1)))

   # plots all points with p < 0.05
   p_sig = subset(p,p<0.05)
   points(lexp[1:length(p_sig)], lobs[1:length(p_sig)], pch=23, cex=0.8, col=color, bg=color)

   # samples 5,000 points from p > 0.05 (to keep file size down)
   n=5001
   i<- c(length(p)- c(0,round(log(2:(n-1))/log(n)*length(p))),1)
   lobs_bottom=subset(lobs[i],lobs[i] <= 3)
   lexp_bottom=lexp[i[1:length(lobs_bottom)]]
   points(lexp_bottom, lobs_bottom, pch=23, cex=0.8, col=color, bg=color)

}


## Reads data
## Header should have column "P"

S <- read.table(input,header=T)

if (stat_type == "Z")
z=S$P

if (stat_type == "CHISQ")
z=sqrt(S$P)

if (stat_type == "PVAL")
z=qnorm(S$P/2)

## calculates lambda
lambda = round(median(z^2,na.rm=T)/qchisq(0.5,df=1),3)
print(lambda)

## Plots axes and null distribution
png(file=paste(output,".QQ.png",sep="."), width=8, height=8, unit="cm", pointsize=4, res=300)

plot(c(0,xmax+0.2), c(0,xmax+0.2), col="black", lwd=3, type="l", xlab="Expected Distribution (-log10 of P value)", ylab="Observed Distribution (-log10 of P value)", las=1, xaxs="i", yaxs="i", axes=F, bty="l",main=c(substitute(paste("",lambda," = ", lam),list(lam = lambda)),expression()))

axis(1,at=seq(0,xmax,1),labels=T)
axis(2,at=seq(0,xmax,1),labels=T,las=2)

## plots data
plotQQ(z,color);

## provides legend
legend(x="topleft",legend=c("Expected","Observed"),pch=23,cex=1,pt.bg=c("black",color),bty="n")

rm(z)

dev.off()
