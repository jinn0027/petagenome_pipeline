library(tidyr)
library(plyr)

v.args <- commandArgs(trailingOnly=TRUE)
in_ <- v.args[1] # input file
dir_out_ <- v.args[2] # output directory
samples <- v.args[3] # samples
date <- v.args[4] # date of KEGG data
db_label <- v.args[5] # KEGG db used
v.samples <- unlist(strsplit(samples, ","))

# input
pathway_ko_ <- paste("/home/user/kegg.20180916/genes/ko/pathway_ko.txt", sep="")
ko_definition_ <-paste("/home/user/kegg.20180916/genes/ko/ko.NAME.DEFINITION.txt", sep="")

df.in <- read.delim(in_, header=T, sep="\t", quote="", stringsAsFactors=FALSE, check.names=F)
df.t <- df.in[,v.samples]

# make id (pathway:KO)
v.id <- paste(df.in$pathway, df.in$KO, sep=":")

#v.thre <- c(0.75, 0.8, 0.85, 0.9, 0.95, 1)
v.thre <- c(0.85)
for(thre in v.thre){
	print(paste("start",thre))

	dir_out_each_ <- paste(dir_out_,"/",db_label,".ko.orf_",thre,sep="")
	if(! file.exists(dir_out_each_)){
		dir.create(dir_out_each_)
	}

	df.0_1 <- df.t
	df.0_1[df.0_1 >= thre] <- 1
	df.0_1[df.0_1 < thre] <- 0

	df.x <- cbind(v.id, df.0_1)
	colnames(df.x)[1] <- "id"

	# Number of KO covered length of KO
	df.y <- ddply(df.x, .(id), colwise(sum))

	# split id (pathway:KO)
	v.id_split <- unlist(strsplit(as.character(df.y$id), ":"))
	v.pathway <- as.numeric(v.id_split[seq(1,length(v.id_split)-1,by=2)])
	v.ko <- v.id_split[seq(2,length(v.id_split),by=2)]
	df.z <- cbind(v.pathway, v.ko, df.y[,v.samples])
	colnames(df.z)[c(1,2)] <- c("pathway", "ko")

	# output 0 / 1 files
	df.pathway_ko <- read.delim(pathway_ko_, header=F, sep="\t", quote="", stringsAsFactors=FALSE, check.names=F)
	df.ko_definition <- read.delim(ko_definition_, header=F, sep="\t", quote="", stringsAsFactors=FALSE, check.names=F)
	v.name <- df.ko_definition[,2]
	v.definition <- df.ko_definition[,3]
	names(v.name) <- df.ko_definition[,1]
	names(v.definition) <- df.ko_definition[,1]
	v.pathway.uniq <- unique(v.pathway)
	for(pathway in v.pathway.uniq){
		# all KOs in a pathway
		v.ko <- unlist(strsplit(df.pathway_ko[df.pathway_ko[,1]==pathway, 3], ","))
		# matrix of KOs detected
		df.out.exist <- df.z[df.z$pathway==pathway,]
		v.ko.exist <- df.out.exist$ko
		# undetected KOs
		v.ko.non <- setdiff(v.ko, v.ko.exist)
		n_non <- length(v.ko.non)
		# matrix of undetected KOs
		mx.non <- cbind(rep(pathway, n_non), v.ko.non, matrix(0, nrow=n_non, ncol=length(v.samples)))
		colnames(mx.non) <- colnames(df.out.exist)
		# final matrix
		df.out <- rbind(df.out.exist, mx.non)
		df.out.def <- transform(df.out, name=v.name[as.character(df.out$ko)], definition=v.definition[as.character(df.out$ko)])
		colnames(df.out.def) <- c("pathway","ko",v.samples,"name","definition")
		df.out.def <- df.out.def[, c("pathway","ko","name","definition",v.samples)]
		write.table(df.out, file=paste(dir_out_each_,"/sample.kegg_ko.",pathway,".count.txt",sep=""), row.names=F, col.names=T, quote=F, sep="\t")
		write.table(df.out.def, file=paste(dir_out_each_,"/sample.kegg_ko.",pathway,".count.definition.txt",sep=""), row.names=F, col.names=T, quote=F, sep="\t")
	}

	# calculate ratio
	df.ko_count <- data.frame(matrix(rep(NA, 1 + length(v.samples)), nrow=1))[numeric(0), ]
	for(pathway in v.pathway.uniq){
		n_ko <- df.pathway_ko[df.pathway_ko[,1]==pathway, 2]
		v.count_ko <- apply(df.z[df.z$pathway==pathway, v.samples], 2, sum)
		df.ko_count <- rbind(df.ko_count, c(pathway, v.count_ko))
	}
	colnames(df.ko_count) <- c("pathway", v.samples)
	write.table(df.ko_count, file=paste(dir_out_,"/",db_label,".pathway.ko_count.orf_",thre,".txt",sep=""), row.names=F, col.names=T, quote=F, sep="\t")
}

