
v.args <- commandArgs(trailingOnly=TRUE)
dir_ <- v.args[1] # directory of input and output
in_ <- v.args[2] # input file
label <- v.args[3] # label for output files

# input
df.in <- read.delim(in_, row.names=1, header=T, check.names=F)

# total nucleotide length in each sample
v.total_length <- c()
for(j in 1:ncol(df.in)){
	v.total_length <- c(v.total_length, sum(as.numeric(df.in[,j])))
}
# avoid dividing by 0
v.total_length[ v.total_length == 0 ] <- 1
names(v.total_length) <- colnames(df.in)
# print
cat(paste("\ntotal length of reads mapped to contigs in each sample", label, "\n"))
print(v.total_length)

# normalization factor decided by the minimum total read length among samples
nf <-  1000000000 # 1e9
# scaling values
v.scale <- v.total_length / nf
cat(paste("\nscaling factors in each sample", label, "\n"))
print(v.scale)

v.scale[1] <- 1 # for length
# normalized value (mapped read nucleotides)
df.norm <- scale(df.in, center=F, scale=v.scale)
# cover rate (mapped read nucleotides / contig length)
df.norm.rate <- cbind(df.norm[,1], t( scale(t(df.norm[,-1]), center=F, scale=as.vector(df.norm[,1])) ))
colnames(df.norm.rate) <- colnames(df.norm)

#------------------------------
# output
# total length
out_tot_len_ <- paste(dir_, "/total_length.",label,".txt",sep="")
write.table(v.total_length, file=out_tot_len_, sep="\t", row.names=T, col.names=F, quote=F)

# scaling factors used for normalization
out_sf_ <- paste(dir_, "/scaling_factor.",label,".txt",sep="")
write.table(v.scale, file=out_sf_, sep="\t", row.names=T, col.names=F, quote=F)

# normalized read lengths covering contigs
out_norm_ <- paste(dir_, "/normalized.length.",label,".txt",sep="")
write.table(cbind(rownames(df.norm), df.norm), file=out_norm_, sep="\t", row.names=F, col.names=T, quote=F)

# normalized read coverages (rate) of contigs
out_rate_ <- paste(dir_, "/normalized.rate.",label,".txt",sep="")
write.table(cbind(rownames(df.norm.rate), df.norm.rate), file=out_rate_, sep="\t", row.names=F, col.names=T, quote=F)


