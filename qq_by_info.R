## To run: Rscript qq_by_info.R inputfile outputfile

rm(list=ls())

input=commandArgs(trailingOnly=T)[1]
output=commandArgs(trailingOnly=T)[2]

## QQ Plot function ##
plotQQ <- function(z,color,cex){
    p <- 2*pnorm(-abs(z))
    p <- sort(p)
    expected <- c(1:length(p))
    lobs <- -(log10(p))
    lexp <- -(log10(expected / (length(expected)+1)))

    ## plots all points with p < 1
    p_sig = subset(p,p<1)
    points(lexp[1:length(p_sig)], lobs[1:length(p_sig)], pch=23, cex=.3, col=color, bg=color)

    ## samples 5000 points from p > 0.01
    n=5001
    i<- c(length(p)- c(0,round(log(2:(n-1))/log(n)*length(p))),1)
    lobs_bottom=subset(lobs[i],lobs[i] <= 2)
    lexp_bottom=lexp[i[1:length(lobs_bottom)]]
    points(lexp_bottom, lobs_bottom, pch=23, cex=cex, col=color, bg=color)
}


## Reads data
## Header must contain P INFO

S <- read.table(input,header=T)

### Plo - by imputation quality (INFO)

pvals_lo0=subset(S,(S$INFO > 1.1))
pvals_lo1=subset(S,(S$INFO > 0.9 & S$INFO < 1.1))
pvals_lo2=subset(S,((S$INFO <= 0.9 & S$INFO > 0.8)))
pvals_lo3=subset(S,((S$INFO <= 0.8 & S$INFO > 0.7)))
pvals_lo4=subset(S,((S$INFO <= 0.7 & S$INFO > 0.6)))
pvals_lo5=subset(S,(S$INFO <= 0.6))

z=qnorm(S$P/2)
z_lo0=qnorm(pvals_lo0$p/2)
z_lo1=qnorm(pvals_lo1$P/2)
z_lo2=qnorm(pvals_lo2$P/2)
z_lo3=qnorm(pvals_lo3$P/2)
z_lo4=qnorm(pvals_lo4$P/2)
z_lo5=qnorm(pvals_lo5$P/2)


## calculates lambda
lambda = round(median(z^2,na.rm=T)/qchisq(0.5,df=1),3)
l0 = round(median(z_lo0^2,na.rm=T)/qchisq(0.5,df=1),3)
l1 = round(median(z_lo1^2,na.rm=T)/qchisq(0.5,df=1),3)
l2 = round(median(z_lo2^2,na.rm=T)/qchisq(0.5,df=1),3)
l3 = round(median(z_lo3^2,na.rm=T)/qchisq(0.5,df=1),3)
l4 = round(median(z_lo4^2,na.rm=T)/qchisq(0.5,df=1),3)
l5 = round(median(z_lo5^2,na.rm=T)/qchisq(0.5,df=1),3)


## Plots axes and null distribution
png(file=paste(output,"impqualQQ.png",sep="."), width=8, height=8,unit="cm", pointsize=4, res=300)
plot(c(0,8), c(0,8), col="black", lwd=3, type="l", xlab="Expected Distribution (-log10 of P value)", ylab="Observed Distribution (-log10 of P value)", xlim=c(0,8), ylim=c(0,8), las=1, xaxs="i", yaxs="i", bty="l",main=c(substitute(paste("QQ plot: ",lambda," = ", lam),list(lam = lambda)),expression()))

## plots data
plotQQ(z,rgb(255,79,0,maxColorValue=255),0.4);
plotQQ(z_lo5,"lightpink2",0.3);
plotQQ(z_lo4,"purple",0.3);
plotQQ(z_lo3,"deepskyblue1",0.3);
plotQQ(z_lo2,"slateblue3",0.3);
plotQQ(z_lo1,"olivedrab2",0.3);
plotQQ(z_lo0,"palevioletred3",0.3);

## provides legend
legend(.25,8,legend=c("Expected (null)","Observed",
substitute(paste("Imp Qual > 1.1 [", lambda," = ", lam, "]"),list(lam = l0)),expression(),
substitute(paste("0.9 < Imp Qual < 1.1 [", lambda," = ", lam, "]"),list(lam = l1)),expression(),
substitute(paste("0.8 < Imp Qual < 0.9 [", lambda," = ", lam, "]"),list(lam = l2)),expression(),
substitute(paste("0.7 < Imp qual < 0.8 [", lambda," = ", lam, "]"),list(lam = l3)),expression(),
substitute(paste("0.6 < Imp qual < 0.7 [", lambda," = ", lam, "]"),list(lam = l4)),expression(),
substitute(paste("Imp qual < 0.6 [", lambda," = ", lam, "]"),list(lam = l5)),expression()),
pch=c((vector("numeric",6)+1)*23), cex=1.1, pt.cex=1.5, pt.bg=c("black",rgb(255,79,0,maxColorValue=255),"palevioletred3","olivedrab2","slateblue3","deepskyblue1","purple","lightpink2"))

rm(z)
dev.off()
