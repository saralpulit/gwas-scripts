## Generates (a) a standard QQ plot, (b) a QQ plot stratified by frequency, and (c) a QQ plot stratified by info score
## To run: Rscript inputfile stat_type output file color1 xmax

## stat_type: P, Z, CHISQ
## File must contain the columns: P, INFO, FRQ

rm(list=ls())

### The input file needs to contain the following headers (for each SNP): FRQ INFO P

input=commandArgs(trailingOnly=T)[1]
stat_type=commandArgs(trailingOnly=T)[2]
output=commandArgs(trailingOnly=T)[3]
color1=commandArgs(trailingOnly=T)[4]
xmax=commandArgs(trailingOnly=T)[5]

xmax <- as.numeric(xmax)
print(xmax)

## Plot function ##
   plotQQ <- function(z,color,cex){
   p <- 2*pnorm(-abs(z))
   p <- sort(p)
   expected <- c(1:length(p))
   lobs <- -(log10(p))
   lexp <- -(log10(expected / (length(expected)+1)))

   ## plots all points
   ## can be adjusted to plot only a subset of points and then randomly sample for non-significant SNPs (see lines below)
   p_sig = subset(p,p<1)
   points(lexp[1:length(p_sig)], lobs[1:length(p_sig)], pch=23, cex=cex, col=color, bg=color)

   ## samples 5,000 points from p > 0.01 (to minimize number of data points)
   ## n=5001
   ## o <- c(length(p)- c(0,round(log(2:(n-1))/log(n)*length(p))),1)
   ## lobs_bottom=subset(lobs[i],lobs[i] <= 2)
   ## lexp_bottom=lexp[i[1:length(lobs_bottom)]]
   ## points(lexp_bottom, lobs_bottom, pch=23, cex=0.8, col=color, bg=color)

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
png(file=paste(output, ".QQ.png", sep=""), width=8, height=8, unit="cm", res=600, pointsize=7)

plot(c(0,xmax+0.2), c(0,xmax+0.2), col="black", lwd=3, type="l", xlab="Expected Distribution (-log10 of P value)", ylab="Observed Distribution (-log10 of P value)", las=1, xaxs="i", yaxs="i", axes=F, bty="l",main=c(substitute(paste("",lambda," = ", lam),list(lam = lambda)),expression()))

axis(1,at=seq(0,xmax,1),labels=T)
axis(2,at=seq(0,xmax,1),labels=T,las=2)

## plots data
plotQQ(z,color1,0.8);

## provides legend
legend(x="topleft",legend=c("Expected","Observed"),pch=23,cex=0.5,pt.bg=c("black",color1),bty="n")

rm(z)
dev.off()

#########################
#### SPLIT UP BY FRQ ####
#########################

## Plot 1 - by coded allele frequency

pvals_lo1=subset(S,(S$FRQ > 0.20 & S$FRQ < 0.8))
pvals_lo2=subset(S,((S$FRQ < 0.20 & S$FRQ > 0.05) | (S$FRQ > 0.8 & S$FRQ < 0.95)))
pvals_lo3=subset(S,((S$FRQ < 0.05 & S$FRQ > 0.01) | (S$FRQ > 0.95 & S$FRQ < 0.99)))
pvals_lo4=subset(S,((S$FRQ < 0.01 & S$FRQ > 0.001) | (S$FRQ > 0.99 & S$FRQ < 0.999)))
pvals_lo5=subset(S,(S$FRQ < 0.001 | S$FRQ > 0.999))

z=qnorm(S$P/2)
z_lo1=qnorm(pvals_lo1$P/2)
z_lo2=qnorm(pvals_lo2$P/2)
z_lo3=qnorm(pvals_lo3$P/2)
z_lo4=qnorm(pvals_lo4$P/2)
z_lo5=qnorm(pvals_lo5$P/2)


## calculates lambda overall and for subsets
lambda = round(median(z^2,na.rm=T)/qchisq(0.5,df=1),3)
l1 = round(median(z_lo1^2,na.rm=T)/qchisq(0.5,df=1),3)
l2 = round(median(z_lo2^2,na.rm=T)/qchisq(0.5,df=1),3)
l3 = round(median(z_lo3^2,na.rm=T)/qchisq(0.5,df=1),3)
l4 = round(median(z_lo4^2,na.rm=T)/qchisq(0.5,df=1),3)
l5 = round(median(z_lo5^2,na.rm=T)/qchisq(0.5,df=1),3)


## Plots axes and null distribution

png(file=paste(output,"freqQQ.png",sep="."), width=8, height=8, unit="cm", res=600, pointsize=7)
plot(c(0,xmax+0.2), c(0,xmax+0.2), col="gray25", lwd=3, type="l", xlab="Expected Distribution (-log10 of P value)", ylab="Observed Distribution (-log10 of P value)", xlim=c(0,8), ylim=c(0,8), las=1, xaxs="i", yaxs="i", bty="l",main=c(substitute(paste("QQ plot: ",lambda," = ", lam),list(lam = lambda)),expression()))

## plots data

plotQQ(z,"navy",0.4);
plotQQ(z_lo5,"royalblue2",0.3);
plotQQ(z_lo4,"purple",0.3);
plotQQ(z_lo3,"forestgreen",0.3);
plotQQ(z_lo2,"aquamarine",0.3);
plotQQ(z_lo1,"olivedrab2",0.3);

## provides legend
legend(x="toplef",legend=c("Expected (null)","Observed",
substitute(paste("MAF > 0.20 [", lambda," = ", lam, "]"),list(lam = l1)),expression(),
substitute(paste("0.05 < MAF < 0.20 [", lambda," = ", lam, "]"),list(lam = l2)),expression(),
substitute(paste("0.01 < MAF < 0.05 [", lambda," = ", lam, "]"),list(lam = l3)),expression(),
substitute(paste("0.001 < MAF < 0.01 [", lambda," = ", lam, "]"),list(lam = l4)),expression(),
substitute(paste("MAF < 0.001 [", lambda," = ", lam, "]"),list(lam = l5)),expression()),
pch=c((vector("numeric",6)+1)*23), cex=0.5, pt.cex=1, pt.bg=c("grey25","navy","olivedrab2","aquamarine","forestgreen","purple","royalblue2"))

rm(z)
dev.off()

##########################
#### SPLIT UP BY INFO ####
##########################

pvals_lo1=subset(S,(S$INFO > 0.9 ))
pvals_lo2=subset(S,((S$INFO <= 0.9 & S$INFO > 0.8)))
pvals_lo3=subset(S,((S$INFO <= 0.8 & S$INFO > 0.7)))
pvals_lo4=subset(S,((S$INFO <= 0.7 & S$INFO > 0.6)))
pvals_lo5=subset(S,(S$INFO <= 0.6))

z=qnorm(S$P/2)
z_lo1=qnorm(pvals_lo1$P/2)
z_lo2=qnorm(pvals_lo2$P/2)
z_lo3=qnorm(pvals_lo3$P/2)
z_lo4=qnorm(pvals_lo4$P/2)
z_lo5=qnorm(pvals_lo5$P/2)


## calculates lambda
lambda = round(median(z^2,na.rm=T)/qchisq(0.5,df=1),3)
l1 = round(median(z_lo1^2,na.rm=T)/qchisq(0.5,df=1),3)
l2 = round(median(z_lo2^2,na.rm=T)/qchisq(0.5,df=1),3)
l3 = round(median(z_lo3^2,na.rm=T)/qchisq(0.5,df=1),3)
l4 = round(median(z_lo4^2,na.rm=T)/qchisq(0.5,df=1),3)
l5 = round(median(z_lo5^2,na.rm=T)/qchisq(0.5,df=1),3)


## Plots axes and null distribution
png(file=paste(output,"impqualQQ.png",sep="."), width=8, height=8,unit="cm", res=600, pointsize=7)
plot(c(0,xmax+0.2), c(0,xmax+0.2), col="black", lwd=3, type="l", xlab="Expected Distribution (-log10 of P value)", ylab="Observed Distribution (-log10 of P value)", xlim=c(0,8), ylim=c(0,8), las=1, xaxs="i", yaxs="i", bty="l",main=c(substitute(paste("QQ plot: ",lambda," = ", lam),list(lam = lambda)),expression()))

## plots data

plotQQ(z,"navy",0.4);
plotQQ(z_lo5,"royalblue2",0.3);
plotQQ(z_lo4,"purple",0.3);
plotQQ(z_lo3,"forestgreen",0.3);
plotQQ(z_lo2,"aquamarine",0.3);
plotQQ(z_lo1,"olivedrab2",0.3);

## provides legend
legend(x="topleft",legend=c("Expected (null)","Observed",
substitute(paste("Imp Qual > 0.9 [", lambda," = ", lam, "]"),list(lam = l1)),expression(),
substitute(paste("0.8 < Imp Qual < 0.9 [", lambda," = ", lam, "]"),list(lam = l2)),expression(),
substitute(paste("0.7 < Imp qual < 0.8 [", lambda," = ", lam, "]"),list(lam = l3)),expression(),
substitute(paste("0.6 < Imp qual < 0.7 [", lambda," = ", lam, "]"),list(lam = l4)),expression(),
substitute(paste("Imp qual < 0.6 [", lambda," = ", lam, "]"),list(lam = l5)),expression()),
pch=c((vector("numeric",6)+1)*23), cex=0.5, pt.cex=1, pt.bg=c("black","navy","olivedrab2","aquamarine","forestgreen","purple","royalblue2"))

rm(z)
dev.off()
