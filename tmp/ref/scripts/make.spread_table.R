library(tidyr)

v.args <- commandArgs(trailingOnly=TRUE)
in_ <- v.args[1] # input file
out_ <- v.args[2] # output file
col_name_1 <- v.args[3] # name of the 1st column

df.in <- read.delim(in_, header=F, stringsAsFactors=F, check.names=F)
df.spread <- df.in %>% spread(V2, V3)
df.spread[is.na(df.spread)] <- 0
colnames(df.spread)[1] <- col_name_1
write.table(df.spread, file=out_, row.names=F, col.names=T, quote=F, sep="\t")

