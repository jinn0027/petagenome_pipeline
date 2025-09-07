
args <- commandArgs(trailingOnly = TRUE)

# input file and parameters
in_ <- args[1]
out <- args[2] # output base
col_id <- as.numeric(args[3]) # columun of id
lineage <- args[4] # lineage to consider (eg: "kingdom,phylum,class,order,family,subfamily,genus,species")
col_lineage <- args[5] # colmun of lineage (eg: "5,6,7,8,9,10,11,12")
thre <- as.numeric(args[6]) # threshold of fraction to assign taxonomy
bln_taxid <- args[7] # blank: no taxid, any: showing taxid

# output
out_ <- paste(out, ".txt", sep="")
out_uc_ <- paste(out, ".uc_.txt", sep="")

# parse arguments of lineage specification
v.lineage <- unlist(strsplit(lineage, ","))
v.col_lineage <- as.numeric(unlist(strsplit(col_lineage, ",")))
# check arguments
if(length(v.lineage) == length(v.col_lineage)){
	cat("specified lineage: ")
	cat(paste(v.col_lineage, v.lineage, sep=":"))
	cat("\n")
}else{
	stop(message="different number of lineages were specified.")
}

# read data
df.in.all <- read.delim(in_, header=F, sep="\t", quote="", colClasses="character", stringsAsFactors=FALSE)
n_col <- ncol(df.in.all)
df.in <- df.in.all[, c(col_id, v.col_lineage)]
v.id <- unique(df.in.all[, col_id])

# header of output (usually: "Contig\tkingdom\tphylum\tclass\torder\tfamily\tsubfamily\tgenus\tspecies")
#v.header <- c("Contig", v.lineage)
v.header <- c("contig", "lca", "rank", v.lineage)
df.out <- data.frame(matrix(rep(NA, length(v.header)), nrow=1))[numeric(0), ]
df.out.uc <- df.out
colnames(df.out) <- v.header
for(id in v.id){
	df.x <- df.in[df.in[,1]==id, ]
	n <- nrow(df.x)
	v.majority <- c()
	lca <- "0:undefined"
	lca_col <- 0
	for(i_col in length(v.lineage):1){
		v.tax <- df.x[, i_col + 1]
		v.table.tax <- sort(table(v.tax), decreasing=T)
		if(v.table.tax[1] > thre * n){
			# taxonomy is affiliated
			tax <- names(v.table.tax)[1]
			if(lca == "0:undefined"){
				lca <- tax
				lca_rank <- v.lineage[i_col]
				lca_col <- i_col
			}
		}else{
			# taxonomy cannot be determined
			tax <- "0:undefined"
		}
		v.majority <- c(tax, v.majority)
	}

	if(is.na(bln_taxid)){
		# remove taxid
		v.majority <- sub("\\d+:", "", v.majority)
		lca <- sub("\\d+:", "", lca)
		uc_lca <- paste("uc_",lca,sep="")
	}else{
		# keep taxid
		uc_lca <- sub(":",":uc_",lca)
	}

	# store LCA in dataframe
	v.out <- c(id, lca, lca_rank, v.majority)
	names(v.out) <- v.header
	df.out <- rbind(df.out, t(as.data.frame(v.out)))

	# uc_
	# store LCA with uc_ information
	v.majority.uc <- sub("undefined","unknown",v.majority)
	if(lca_col > 0 & lca_col < length(v.lineage)){
		# set "uc_" + lca name under the rank of lca
		for(i_col in (lca_col + 1):length(v.lineage)){
			v.majority.uc[i_col] <- uc_lca
		}
	}
	v.out.uc <- c(id, lca, lca_rank, v.majority.uc)
	df.out.uc <- rbind(df.out.uc, t(as.data.frame(v.out.uc)))
}

#----- write output data -----
write.table(df.out, file=out_, sep="\t", row.names=F, col.names=F, quote=F)
write.table(df.out.uc, file=out_uc_, sep="\t", row.names=F, col.names=F, quote=F)

