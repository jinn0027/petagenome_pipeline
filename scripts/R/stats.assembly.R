
args <- commandArgs(trailingOnly=TRUE)
in_ <- args[1]; # file containing length information
col <- as.numeric(args[2]) # column that contains length information
df.in <- read.delim(in_, header=FALSE, sep="\t", quote="")
v.len <- sort(as.numeric(df.in[,col]), decreasing=TRUE)
total <- sum(v.len)
idx90 <- sum(cumsum(v.len) <= total * 0.9) + 1
idx80 <- sum(cumsum(v.len) <= total * 0.8) + 1
idx70 <- sum(cumsum(v.len) <= total * 0.7) + 1
idx60 <- sum(cumsum(v.len) <= total * 0.6) + 1
idx50 <- sum(cumsum(v.len) <= total * 0.5) + 1
idx40 <- sum(cumsum(v.len) <= total * 0.4) + 1
idx30 <- sum(cumsum(v.len) <= total * 0.3) + 1
idx20 <- sum(cumsum(v.len) <= total * 0.2) + 1
idx10 <- sum(cumsum(v.len) <= total * 0.1) + 1
n90 <- v.len[idx90]
n80 <- v.len[idx80]
n70 <- v.len[idx70]
n60 <- v.len[idx60]
n50 <- v.len[idx50]
n40 <- v.len[idx40]
n30 <- v.len[idx30]
n20 <- v.len[idx20]
n10 <- v.len[idx10]
max_l <- max(v.len)

#print(paste("Total length:",total))
#print(paste("# contigs:",length(v.len)))
#print(paste("N50:",n50))
v.header <- c("total_length", "N_contigs", "N90", "N80", "N70", "N60", "N50", "N40", "N30", "N20", "N10", "max")
v.out <- c(total, length(v.len), n90, n80, n70, n60, n50, n40, n30, n20, n10, max_l)
cat(t(v.header), sep="\t"); cat("\n");
cat(t(v.out), sep="\t"); cat("\n");

# N50
# N50 statistic defines assembly quality. Given a set of contigs, each with its
# own length, the N50 length is defined as the shortest sequence length at 50%
# of the genome. It can be thought of as the point of half of the mass of the
# distribution; the number of bases from all contigs shorter than the N50 will
# be close to the number of bases from all contigs longer than the N50.
# For example, 9 contigs with the lengths 2,3,4,5,6,7,8,9,10, their sum is 54,
# the size of the genome also happens to be 54. 50% of this assembly would be
# 2+3+4+5+6+7=27 (>25) Thus the N50=7 which is that contig along with the
# smaller contigs that contains half of sequence of a particular genome.
# Note: When comparing N50 values from different assemblies, the assembly sizes
# must be the same size in order for N50 to be meaningful.
