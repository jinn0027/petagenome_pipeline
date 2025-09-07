#Install these packages if they are not already installed
#install.packages(c("coin","exactRankTests","BiocManager","ggplot2","grid","RColorBrewer","plyr","ggpubr","ggrepel"))

#install package qvalue if it's not already installed
#library(BiocManager)
#BiocManager::install("qvalue") 

library(coin)
library(exactRankTests)
library(qvalue)

library(ggplot2)
library(grid)
library(RColorBrewer)
library(plyr)
library(ggpubr)
library(ggrepel)

# Define control_number by number of samples of the control (In data_table, control samples must be placed sequentially first, followed by target samples. 
#Otherwise, v.samples1 and v.samples2 must be defined manually, such as v.samples1<-c("sample001", "sample002", ..). etc. )

control_number <- 19
target_name <- "E.faecalis"
control_name<- "E.faecium"
significant_target <-paste("Significant in",target_name)
significant_control <-paste("Significant in",control_name)

#Define base line and threshold mean
base_line<- 0.0001
thre_mean <- 0.0001

output_rank <- "ko"
dir_project<-"./"
merge_label<-"GVHD"

data_table <-paste0("./Figure1d_dataset_GVHD_Nature.txt")
out_all_ <- paste0("./",output_rank,".KEGG.coverage.all.p_value.two_sided.",merge_label,".add_base_line_",base_line,".txt",sep="")

# Extracting column names
df.in.all <-read.table(data_table, header=T, sep="\t", quote="", stringsAsFactors=FALSE, check.name=FALSE)
col_names <- colnames(df.in.all)

# Creating v.samples1 vector (control)
v.samples1 <- col_names[1:control_number+1]

# Creating v.samples2 vector (target)
v.samples2 <- col_names[(control_number + 2):length(col_names)]

# Printing the vectors
print(v.samples1)
print(v.samples2)

df.ra<- df.in.all
colnames(df.ra) <- c("KO", v.samples1, v.samples2)
print(colnames(df.ra))

v.p.wilcox <- c()
v.p.t <- c()
v.fold_change <- c()
v.mean_1 <- c()
v.mean_2 <- c()
for(i in 1:nrow(df.ra)){
        # add base_line abundance
        v.x_1 <- as.numeric(df.ra[i, v.samples1]) + base_line
        v.x_2 <- as.numeric(df.ra[i, v.samples2]) + base_line
        p.wilcox <- wilcox.exact(v.x_1, v.x_2, paired = FALSE, alternative = "two.sided")$p.value
        fold_change <- mean(v.x_2) / mean(v.x_1)
        v.p.wilcox <- c(v.p.wilcox, p.wilcox)
        v.fold_change <- c(v.fold_change, fold_change)
        v.mean_1 <- c(v.mean_1, mean(v.x_1))
        v.mean_2 <- c(v.mean_2, mean(v.x_2))
}

# Q-value
v.q.wilcox <- qvalue(v.p.wilcox)$qvalues
df.out <- data.frame(id=rownames(df.ra), name=df.ra$"KO", p_wilcox=v.p.wilcox, q_wilcox=v.q.wilcox, fold_change=v.fold_change, mean_1=v.mean_1, mean_2=v.mean_2)
write.table(df.out, file=out_all_, row.names=F, col.names=T, quote=F, sep="\t")


#adjust these threshold values
base_size <- 36
thre_log_q <- -log10(0.05)
thre_log_fc <- log2(1.5)
ggrepel.max.overlaps = Inf

df.in<-read.table(out_all_, header=T, sep="\t", quote="", stringsAsFactors=FALSE, check.name=FALSE)

# output
out_pdf_ <- paste("./volcano.",output_rank,".coverage.",merge_label,".add_base_line_",base_line,".two.pdf",sep="")
out_eps_ <- paste("./volcano.",output_rank,".coverage.",merge_label,".add_base_line_",base_line,".two.eps",sep="")

v.out_type <- c("pdf","eps")
for(out_type in v.out_type){
        if(out_type == "pdf"){
# PDF
                pdf(file=out_pdf_, colormodel="cmyk", pointsize=30,  width=20, height=20)
}else if(out_type == "eps"){
                # EPS
                postscript(file=out_eps_, colormodel="cmyk", pointsize=10,  width=50, height=50, horizontal=F)
                setEPS()
        }
        tryCatch({

                grid.newpage() # make new blank figure
                pushViewport(viewport(layout=grid.layout(4, 8)))
                vp.i  <- viewport(layout.pos.row=1:4, layout.pos.col=1:8)

                #------------------------------
                # volcano plot
                #------------------------------
               df.x <- df.in[df.in[,"mean_1"] > thre_mean | df.in[,"mean_2"] > thre_mean, ]
                print(nrow(df.x))
                v.x <- log2(df.x[,"fold_change"])
                v.y <- -log10(df.x[,"q_wilcox"])
#                v.name <- df.x[,"name"]
                v.name <-"" 
               v.select <- rep("Not significant", nrow(df.x))
                v.select[v.x > thre_log_fc & v.y > thre_log_q] <- significant_target
                v.select[v.x < -thre_log_fc & v.y > thre_log_q] <- significant_control
#               v.name<-ifelse(v.name == "K01791","K01791: wecB",
#                       ifelse(v.name == "K08605","K08605: gelE",
#                       ifelse(v.name == "K05946","K05946: tagA, tarA",
#                       ifelse(v.name == "K07173","K07173: luxS",""))))  
                x_lab_ <- paste("log2(",target_name,"/",control_name,")")
                max_value_x <- ceiling(max(log2(df.x[,"fold_change"])))
                min_value_x <-floor(min(log2(df.x[,"fold_change"])))
                max_value_y <-ceiling(max(-log10(df.x[,"q_wilcox"])))
                min_value_y <-floor(min(-log10(df.x[,"q_wilcox"])))
                df.plot <- data.frame(x=v.x, y=v.y, name=v.name, significance=factor(v.select, levels=c("Not significant", significant_target,significant_control)))
                p <- ggplot(df.plot, aes(x=x, y=y))
                p <- p + theme_bw(base_size=base_size)
                p <- p + ggtitle("Volcano plot of gene (KO) expression")+theme(plot.title = element_text(hjust = 0.5)) 
                p <- p + geom_point(alpha=1.0, aes(colour=significance), shape=16, size = 4)+theme(axis.text.x = element_text(color = "black"))+theme(axis.text.y = element_text(color = "black"))
                p <- p + scale_colour_manual(values=c("gray40","#ff0000","#0080ff"))
                p <- p + geom_text_repel(size=12, nudge_x=0.5,nudge_y=-0.5,color="black", aes(label=name),max.overlaps = Inf)
                p + geom_segment(aes(x = x, y = y, xend = x, yend = y - 1), color = "black") + geom_segment(aes(x = x, y = y - 1, xend = x, yend = y - 2), color = "black")
                p <- p + xlab(x_lab_) + ylab("-log10(q-value)")
                p <- p + xlim(c(min_value_x,max_value_x))
                p <- p + ylim(c(0,max_value_y))
                p <- p + geom_hline(yintercept = thre_log_q, color="gray", size=1, linetype="longdash")
                p <- p + geom_vline(xintercept = thre_log_fc, color="gray", size=1, linetype="longdash")
                p <- p + geom_vline(xintercept = -thre_log_fc, color="gray", size=1, linetype="longdash")
                print(p, vp=vp.i)

        }, finally = { dev.off() }
        )

}



