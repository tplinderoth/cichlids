#!/usr/bin/env Rscript

chr <- c(1:20, 22, 23)
qcdir <- '/home/tl483/rds/rds-rd109-durbin-group/projects/cichlid/cichlid_g.vcf/main_set_2019-07-10/fAstCal1.2/GenotypeCorrected/depth_filter/malawi_variants/qc/'
outplot <- '/home/tl483/rds/rds-rd109-durbin-group/projects/cichlid/cichlid_g.vcf/main_set_2019-07-10/fAstCal1.2/GenotypeCorrected/depth_filter/malawi_variants/qc/malawi_callset_DP_plot.pdf'
pdf(outplot)

x <- data.frame(CHROM=NULL, POS=NULL, DP=NULL)

par(mfrow=c(4,4))

for (i in chr) {
	file <- paste0(qcdir,"malawi_cichlid_variants_chr",i,"_DP.txt")
	y <- read.table(file, head=TRUE, na.string=".")
	x <- rbind(x,y)
	y[,3][which(is.infinite(y[,3]))] <- NA # mask inf and -inf values
	cat(paste0("plotting chr",i,"\n"))
	hist(y$DP, main=paste0("chr",i), ylab="Number sites", xlab="DP", xlim=c(0,quantile(y$DP,0.99,na.rm=TRUE)+sqrt(var(y$DP,na.rm=TRUE))),breaks=100000)
}

par(mfrow=c(1,1))

# Genome-wide stats and plot

cat("Calculating genome-wide DP stats and plotting\n")

dpmed <- median(x$DP, na.rm=TRUE)
dpa <- round(dpmed - (0.25*dpmed))
dpb <- round(dpmed + (0.25*dpmed))

cat("Genome median DP: ",dpmed,"\n")
cat("Genome DP lower cutoff: ",dpa,"\n")
cat("Genome DP upper cutoff: ",dpb,"\n")

hist(x$DP, main='Genome', ylab="Number sites", xlab="DP", xlim=c(0,quantile(x$DP,0.99,na.rm=TRUE)+sqrt(var(x$DP,na.rm=TRUE))),breaks=100000)
abline(v=dpa, col="red")
abline(v=dpb, col="red")
legend('topright', c(paste0("lower: ",dpa), paste0("upper: ",dpb)), bty='n',cex=0.75)

invisible(dev.off())
