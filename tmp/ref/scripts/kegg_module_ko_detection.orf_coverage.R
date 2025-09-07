library(tidyr)
library(plyr)

v.args <- commandArgs(trailingOnly=TRUE)
in_ <- v.args[1] # input file
dir_out_ <- v.args[2] # output directory
samples <- v.args[3] # samples
date <- v.args[4] # date of KEGG data
db_label <- v.args[5] # KEGG db used
thresholds <- v.args[6] # thresholds for ORF coverage
v.samples <- unlist(strsplit(samples, ","))
v.thre <- as.numeric(unlist(strsplit(thresholds, ",")))

# input
module_ko_ <- paste("/home/user/kegg.20180916/genes/ko/module_ko.txt", sep="")
ko_definition_ <-paste("/home/user/kegg.20180916/genes/ko/ko.definition.txt", sep="")


df.in <- read.delim(in_, header=T, sep="\t", quote="", stringsAsFactors=FALSE, check.names=F)
df.t <- df.in[,v.samples]

# make id (module:KO)
v.id <- paste(df.in$module, df.in$KO, sep=":")
df.x <- cbind(v.id, df.t)
colnames(df.x)[1] <- "id"

# max covered length of KO
df.y <- ddply(df.x, .(id), colwise(max))

#v.thre <- c(0.75, 0.8, 0.85, 0.9, 0.95, 1)
#v.thre <- c(0.85)
for(thre in v.thre){
	print(paste("start",thre))

	dir_out_each_ <- paste(dir_out_,"/",db_label,".ko.orf_",thre,sep="")
	if(! file.exists(dir_out_each_)){
		dir.create(dir_out_each_)
	}
	df.0_1 <- df.y[,v.samples]
	df.0_1[df.0_1 >= thre] <- 1
	df.0_1[df.0_1 < thre] <- 0

	# split id (module:KO)
	v.id_split <- unlist(strsplit(as.character(df.y$id), ":"))
	v.module <- v.id_split[seq(1,length(v.id_split)-1,by=2)]
	v.ko <- v.id_split[seq(2,length(v.id_split),by=2)]
	df.z <- cbind(v.module, v.ko, df.0_1)
	colnames(df.z)[c(1,2)] <- c("module", "ko")

	# output 0 / 1 files
	df.module_ko <- read.delim(module_ko_, header=F, sep="\t", quote="", stringsAsFactors=FALSE, check.names=F)
	df.ko_definition <- read.delim(ko_definition_, header=F, sep="\t", quote="", stringsAsFactors=FALSE, check.names=F)
	v.definition <- df.ko_definition[,2]
	names(v.definition) <- df.ko_definition[,1]
	v.module.uniq <- unique(v.module)
	for(module in v.module.uniq){
		# all KOs in a module
		v.ko <- unlist(strsplit(df.module_ko[df.module_ko[,1]==module, 3], ","))
		# matrix of KOs detected
		df.out.exist <- df.z[df.z$module==module,]
		v.ko.exist <- df.out.exist$ko
		# undetected KOs
		v.ko.non <- setdiff(v.ko, v.ko.exist)
		n_non <- length(v.ko.non)
		# matrix of undetected KOs
		mx.non <- cbind(rep(module, n_non), v.ko.non, matrix(0, nrow=n_non, ncol=length(v.samples)))
		colnames(mx.non) <- colnames(df.out.exist)
		# final matrix
		df.out <- rbind(df.out.exist, mx.non)
		df.out.def <- transform(df.out, definition=v.definition[as.character(df.out$ko)])
		colnames(df.out.def) <- c("module","ko",v.samples,"definition")
		df.out.def <- df.out.def[, c("module","ko","definition",v.samples)]
		write.table(df.out, file=paste(dir_out_each_,"/sample.kegg_ko.",module,".0_1.txt",sep=""), row.names=F, col.names=T, quote=F, sep="\t")
		write.table(df.out.def, file=paste(dir_out_each_,"/sample.kegg_ko.",module,".0_1.definition.txt",sep=""), row.names=F, col.names=T, quote=F, sep="\t")
	}
	
	# calculate ratio
	mx.ko_ratio <- c()
	for(module in v.module.uniq){
		n_ko <- df.module_ko[df.module_ko[,1]==module, 2]
		v.count_ko <- apply(df.z[df.z$module==module, v.samples], 2, sum)
		v.ratio_ko <- as.numeric(v.count_ko / n_ko)
		mx.ko_ratio <- rbind(mx.ko_ratio, v.ratio_ko)
	}
	df.ko_ratio <- cbind(v.module.uniq, mx.ko_ratio)
	colnames(df.ko_ratio) <- c("module", v.samples)
	write.table(df.ko_ratio, file=paste(dir_out_,"/",db_label,".module.ko_ratio.orf_",thre,".txt",sep=""), row.names=F, col.names=T, quote=F, sep="\t")
}

