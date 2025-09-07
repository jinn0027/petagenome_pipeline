
library(ggplot2)
library(grid)
library(RColorBrewer)
library(dplyr)
library(tidyr)
library(scales)
library(vegan)

v.args <- commandArgs(trailingOnly=TRUE)
in_ <- v.args[1] # input file
dir_out_ <- v.args[2] # output directory
out_pdf_ <- v.args[3] # output pdf file
host <- v.args[4] # host organism
samples <- v.args[5] # samples
max_taxa <- as.numeric(v.args[6]) # max taxa to show
base_size <- as.numeric(v.args[7]) # base_size in ggplot2

v.samples <- unlist(strsplit(samples, ","))
v.out_bac_level <- c("class", "order", "family")

# load data
df.all.in <- read.delim(in_, header=TRUE, sep="\t", quote="", stringsAsFactors=F, check.names=F) # check.names=F is to avoid automatic coversion of column names that are not syntactically valid names.
df.phylum <- df.all.in[df.all.in[,1]=="Bacteria" & df.all.in[,2]=="phylum", -c(1,2)]
df.class <- df.all.in[df.all.in[,1]=="Bacteria" & df.all.in[,2]=="class", -c(1,2)]
df.order <- df.all.in[df.all.in[,1]=="Bacteria" & df.all.in[,2]=="order", -c(1,2)]
df.family <- df.all.in[df.all.in[,1]=="Bacteria" & df.all.in[,2]=="family", -c(1,2)]
df.genus <- df.all.in[df.all.in[,1]=="Bacteria" & df.all.in[,2]=="genus", -c(1,2)]
df.species <- df.all.in[df.all.in[,1]=="Bacteria" & df.all.in[,2]=="species", -c(1,2)]

#------------------------------
# diversity and richness
#------------------------------
v.div.shan <- c()
v.div.simp <- c()
v.div.invsimp <- c()
v.richness <- c()
for(sample in (v.samples)){
	v.x <- df.species[,sample]
	names(v.x) <- df.species[,"lineage"]
	# diversity and richness
	shan <- diversity(v.x, "shannon")
	simp <- diversity(v.x, "simpson")
	invsimp <- diversity(v.x, "inv")
	richness <- specnumber(v.x)
	# diversity and richness
	v.div.shan <- c(v.div.shan, shan)
	v.div.simp <- c(v.div.simp, simp)
	v.div.invsimp <- c(v.div.invsimp, invsimp)
	v.richness <- c(v.richness, richness)
}
names(v.div.shan) <- v.samples
names(v.div.simp) <- v.samples
names(v.div.invsimp) <- v.samples
names(v.richness) <- v.samples

# write table of diversity and richness
out_vegan_ <- paste(dir_out_,"/diversity_indeces.richness.txt",sep="")
df.vegan <- data.frame(sample=v.samples, shannon=v.div.shan, simpson=v.div.simp, invsimp=v.div.invsimp, richness=v.richness)
write.table(df.vegan, file=out_vegan_, sep="\t", row.names=F, col.names=T, quote=F)

# output
pdf(file=out_pdf_, colormodel="cmyk", pointsize=8, paper="a4", width=0, height=0)
tryCatch({
	#--------------------------------------------------
	# phylum level composition
	#--------------------------------------------------
	print("01 pylum level composition")
	if(host == "human"){
		v.fix_p <- c("p__Actinobacteria", "p__Bacteroidetes", "p__Firmicutes", "p__Proteobacteria", "p__Verrucomicrobia", "p__Fusobacteria")
	}else if(host == "mouse"){
		v.fix_p <- c("p__Actinobacteria", "p__Bacteroidetes", "p__Firmicutes", "p__Proteobacteria", "p__Verrucomicrobia", "p__Deferribacteres")
	}

	df.x <- df.phylum[, c("name", v.samples)]
	df.x$name <- paste("p__", df.x$name, sep="")
	print(paste("#phyla: ",nrow(df.phylum),sep=""))

	# normalize
	v.total <- apply(df.x[, -1], 2, sum)
	df.x[,-1] <- scale(df.x[,-1], center=FALSE, scale=v.total)
	# write table
	out_table_ <- paste(dir_out_,"/relative_abundance.phylum.txt",sep="")
	write.table(df.x, file=out_table_, sep="\t", row.names=F, col.names=T, quote=F)

	# set pylum other
	v.name <- df.x[,"name"]
	for(name_other in v.name[!v.name %in% v.fix_p]){
		df.x[df.x$name == name_other, "name"] <- "Other"
	}

	# set phylum colors
	v.out_color <- c(brewer.pal(length(v.fix_p), "Set1"), "#AAAAAA")
	names(v.out_color) <- c(v.fix_p, "Other")

	# check existing phylum names
	v.phylum <- c()
	for(fix_name in c(v.fix_p,"Other")){
		if(sum(df.x$name == fix_name) > 0){
			v.phylum <- c(v.phylum, fix_name)
		}
	}
	# update color code
	v.out_color <- v.out_color[v.phylum]

	# page layout
	grid.newpage() # make new blank figure
	pushViewport(viewport(layout=grid.layout(20, 5)))
	vp.1 <- viewport(layout.pos.row=1:6, layout.pos.col=1:5)
	vp.2 <- viewport(layout.pos.row=8:12, layout.pos.col=1:5)
	vp.3 <- viewport(layout.pos.row=14:17, layout.pos.col=1:5)
	vp.4 <- viewport(layout.pos.row=19:20, layout.pos.col=1:5)

	# barplot
	df.bar <- df.x %>% gather(sample, abundance, -name)
	p <- ggplot(df.bar, aes(x=factor(sample, levels=v.samples), y=abundance, fill=factor(name, levels=rev(v.phylum))))
	p <- p + theme_bw(base_size=base_size)
	p <- p + geom_bar(stat="identity", position = "fill")
	p <- p + ggtitle(paste("Bacterial pylum relative abundance")) + xlab("Samples") + ylab("Relative abundance")
	p <- p + theme(axis.text.x = element_text(angle=90, hjust=1, vjust =0.5, size=base_size*0.9))
	p <- p + scale_fill_manual(values=rev(v.out_color))
	p <- p + guides(fill = guide_legend(ncol=1, reverse=FALSE, title="Bacterial phylum"))
	print(p, vp=vp.1)
	print(p, vp=vp.2)
	print(p, vp=vp.3)
	print(p, vp=vp.4)


	#--------------------------------------------------
	# class/order/family level composition
	#--------------------------------------------------
	for(out_bac_level in v.out_bac_level){
		if(out_bac_level == "class"){
			print("02 class level composition")
			df.x <- df.class[, c("lineage", v.samples)]
			print(paste("#classes: ",nrow(df.class),sep=""))
		}else if(out_bac_level == "order"){
			print("03 order level composition")
			df.x <- df.order[, c("lineage", v.samples)]
			print(paste("#orders: ",nrow(df.order),sep=""))
		}else if(out_bac_level == "family"){
			print("04 family level composition")
			df.x <- df.family[, c("lineage", v.samples)]
			print(paste("#families: ",nrow(df.family),sep=""))
		}
	
		#------------------------------
		# select classes/families to show
		#------------------------------
		print("  .1 select taxa to show")
		# normalize
		v.total <- apply(df.x[, -1], 2, sum)
		df.x[,-1] <- scale(df.x[,-1], center=FALSE, scale=v.total)
		# sort taxonomic labels by their abundance
		v.taxonomy.sum <- apply(df.x[, -1], 1, sum)
		names(v.taxonomy.sum) <- df.x[,1]
		df.x <- df.x[rev(order(v.taxonomy.sum)),]
		v.taxonomy.sort <- df.x[,1]
		# write table
		out_table_ <- paste(dir_out_,"/relative_abundance.",out_bac_level,".txt",sep="")
		write.table(df.x, file=out_table_, sep="\t", row.names=F, col.names=T, quote=F)

		# select taxa to show
		v.taxonomy.select <- c()
		mx.tax <- c() # matrix of selected taxonomy labels (phylum, class and class-level labels)
		for(i in seq(v.taxonomy.sort)){
			tax <- v.taxonomy.sort[i]
			v.split.tax <- strsplit(tax, "|", fixed=TRUE)[[1]]
			if(v.split.tax[2] %in% v.fix_p){
				if(length(v.taxonomy.select) < max_taxa){
					v.taxonomy.select <- c(v.taxonomy.select, tax)
					if(out_bac_level == "class"){
						mx.tax <- rbind(mx.tax, v.split.tax[c(2,3)]) # 2:phulum, 3:class
					}else if(out_bac_level == "order"){
						mx.tax <- rbind(mx.tax, v.split.tax[c(2,3,4)]) # 2:phulum, 3:class, 4:order
					}else if(out_bac_level == "family"){
						mx.tax <- rbind(mx.tax, v.split.tax[c(2,3,5)]) # 2:phulum, 3:class, 5:family
					}
				}
			}
		}
	
		# put label "Other" for unselected classes
		for(taxonomy_other in v.taxonomy.sort[!v.taxonomy.sort %in% v.taxonomy.select]){
			df.x[df.x[,1] == taxonomy_other, 1] <- "Other"
		}
		df.x.select <- df.x[df.x[,1]!="Other", ]
		v.x.other <- c(0, apply(df.x[df.x[,1]=="Other", -1], 2, sum))
	
		#------------------------------
		# sort data
		#------------------------------
		print("  .2 sort data")
		# sum for selected taxonomy
		v.sum.select <- apply(df.x.select[,-1], 1, sum)
		names(v.sum.select) <- df.x.select[, 1]
	
		# order of phyla shown is fixed
		n_phylum <- length(v.fix_p)
		v.order.p <- 1:n_phylum
		names(v.order.p) <- v.fix_p
	
		# order of taxa other than "Other"
		v.sum.2nd <- tapply(v.sum.select, mx.tax[,2], sum)
		v.order <- order(v.order.p[mx.tax[,1]], -v.sum.2nd[mx.tax[,2]], -v.sum.select)
		# put "Other" at the last
		df.x.sort <- rbind(df.x.select[v.order, ], v.x.other)
		mx.tax.sort <- mx.tax[v.order,]
		# taxonomy label to show
		if(out_bac_level == "class"){
			v.out_taxonomy <- c(mx.tax.sort[,2], "Other") # class
		}else{
			v.out_taxonomy <- c(mx.tax.sort[,3], "Other") # order / family
		}
		df.x.sort[, 1] <- v.out_taxonomy
		colnames(df.x.sort)[1] <- "taxonomy"
	
		#------------------------------
		# set color code
		#------------------------------
		print("  .3 set color code")
		v.p_color <- brewer.pal(n_phylum, "Set1"); v.p_color[3] <- "#4DCF4A" # color for phylum
		color_base <- "white" # gradient base color
		v.out_color <- c() # final color code vector
		for(p in 1:n_phylum){
			v.class <- unique(mx.tax.sort[mx.tax.sort[,1]==v.fix_p[p], 2])
			# color for ramping among the same phylum
			if(p == 1){
				v.c_color <- v.p_color[c(-4)]
			}else if(p == 2){
				v.c_color <- v.p_color[c(-3,-4)][c(p:7,1:(p-1))]
			}else if(p == 3){
				v.c_color <- v.p_color[c(-4)][c(p:7,1:(p-1))]
			}else if(p == 4){
				v.c_color <- v.p_color[c(-7,-8)][c(p:6,1:(p-1))]
			}else if(p > 4){
				v.c_color <- v.p_color[c(p:8,1:(p-1))]
			}
			for(c in seq(v.class)){
				# number of family in a class
				n_family <- sum(mx.tax.sort[,1]==v.fix_p[p] & mx.tax.sort[,2]==v.class[c])
				# set gradient colors for classes in the class
				colfunc_p <- colorRampPalette(c(v.p_color[p], v.c_color[c]))
				p_color_base <- colfunc_p(10)[4]
				colfunc <- colorRampPalette(c(p_color_base, color_base))
				v.color <- colfunc(n_family + 2)[1:n_family]
				# set to the final color code vector
				v.out_color <- c(v.out_color, v.color)
			}
		}
		v.out_color <- c(v.out_color, "#AAAAAA") # grey for "Other"
	
		#------------------------------
		# draw barplot
		#------------------------------
		print("  .4 draw barplot")
		# page layout
		grid.newpage() # make new blank figure
		pushViewport(viewport(layout=grid.layout(20, 5)))
		vp.1 <- viewport(layout.pos.row=1:6, layout.pos.col=1:5)
		vp.2 <- viewport(layout.pos.row=8:12, layout.pos.col=1:5)
		vp.3 <- viewport(layout.pos.row=14:17, layout.pos.col=1:5)
		vp.4 <- viewport(layout.pos.row=19:20, layout.pos.col=1:5)

		df.bar <- df.x.sort %>% gather(sample, abundance, -taxonomy)
		# barplot (fill)
		p <- ggplot(df.bar, aes(x=factor(sample, levels=v.samples), y=abundance, fill=factor(taxonomy, levels=rev(v.out_taxonomy))))
		p <- p + theme_bw(base_size=base_size)
		p <- p + geom_bar(stat="identity", position = "fill")
		p <- p + scale_fill_manual(values=rev(v.out_color))
		p <- p + ggtitle(paste("Bacterial ",out_bac_level," relative abundance",sep="")) + xlab("Samples") + ylab("Relative abundance")
		p <- p + guides(fill = guide_legend(ncol=1, reverse=FALSE, title=paste("Bacterial ",out_bac_level,sep="")))
		p <- p + theme(axis.text.x = element_text(angle=90, hjust=1, vjust =0.5, size=base_size*0.9))
		p <- p + theme(legend.key.size = unit(0.75,"line"), legend.text = element_text(size=base_size*0.9))
		print(p, vp=vp.1)
		print(p, vp=vp.2)
		print(p, vp=vp.3)
		print(p, vp=vp.4)
	}

}, finally = { dev.off() }
) 
