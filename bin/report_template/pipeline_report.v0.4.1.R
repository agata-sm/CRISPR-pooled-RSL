#!/usr/bin/env Rscript

# code for CRISPR pipeline report
# use with crispr_pipeline_report_v0.4.Rmd

#### IMPORTANT!

#### VARIABLES TO DEFINE IN CLI
#proj.dir
#data.type
#proj.name.pref
#organism



#Rscript report_launcher.R  proj.dir proj.name.pref data.type organism &> reads.report.stderrout
#OBS! positional arguments


## ---- dirs


# types of data for chunk evaluation
is.RSL=isTRUE(data.type=="RSL")


#metadata
comparisons.file=file.path(proj.dir,"metadata","comparisons.txt")
samples.file=file.path(proj.dir,"metadata","metadata.txt")

#indata dirs
mageck_datadir=file.path(proj.dir,"results",data.type,"rra_mageck")

if(is.RSL){
  ctable_datadir=file.path(proj.dir,"results",data.type,"count_table_filtered")
}else{
  ctable_datadir=file.path(proj.dir,"results",data.type,"count_table_mageck")
}



#file seq stats >  chunks only executed for data type RSL
file.seq.stats="counter.stdout.parsed.txt"
file.seqstats=file.path(proj.dir,"results/RSL/Counter",file.seq.stats)
seqstatsSTAUS=file.exists(file.seqstats)

#file summary stats by mageck > chunks only executed for data type reads


#results dir
resdir=file.path("results",data.type,"rra_annotation")
dir.create(resdir, recursive = TRUE)

plotdir=file.path(resdir,"plots")
dir.create(plotdir)

plotdir_qc=file.path(plotdir,"QC_plots")
dir.create(plotdir_qc)

## ---- prep_environment
options(connectionObserver = NULL) # https://github.com/rstudio/rstudio/issues/9219


library(tidyverse)
library(ggplot2)
library(DT)
library(plotly)
library(viridis)
library(ggrepel)
library(dplyr)
library(reshape2)
library(htmltools)
library(magrittr)
library(MASS)
library(cowplot)

library(MAGeCKFlute)
#library(idr)
library(edgeR)

library(fgsea)
library(clusterProfiler)
library(enrichplot)
library(DOSE)
library(biomaRt)
library(ReactomePA)
library(reactome.db)

require(org.Hs.eg.db)
require(org.Mm.eg.db)

# for plots
theme_pca=theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),
    panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
    strip.background=element_blank(),axis.text.x=element_text(colour="black"),
    axis.text.y=element_text(colour="black"),axis.ticks=element_line(colour="black"),plot.margin=unit(c(1,1,1,1),"line"),
    legend.text=element_text(size=10) )

fontsize =theme(axis.text=element_text(size=16), axis.title=element_text(size=20))

fig_n=0
tab_n=0
tab_n_sv=0

#set here - significance cutoff for RRA (mainly used for plot labeling)
FDR.CO=0.05
mycolour = c('ns'="gray80", 'dn'= "#377eb8", 'up'="#e41a1c")


source("./crispr_pipeline_report_functions.R")


#######################################################
########## descriptive stats for read summarisation

## ---- data_countingstats_rsl

samples.tab=read.table(samples.file, sep="\t", header=TRUE, blank.lines.skip=TRUE)
samples.tab$library=samples.tab$sample


if(is.RSL){
    if(seqstatsSTAUS){

  seqstats=read.delim(file.seqstats, header = TRUE, sep = "\t", quote = "\"", dec = ".", fill = TRUE, row.names=NULL)
  seqstats$missing_sgRNA=seqstats$sgRNA_total-seqstats$sgRNA_file

  seqstats$file_fastq=factor(seqstats$file_fastq, levels=samples.tab$file)
  seqstats=seqstats[match(samples.tab$file, seqstats$file_fastq),]
  seqstats$sample=samples.tab$sample

  seq_plot=ggplot(data=seqstats, aes(x=sample, y=fraction_reads_assigned, fill=sample)) + geom_bar(stat="identity") +
      scale_fill_viridis(discrete=TRUE, alpha=0.35,option="turbo") +
      theme_bw() +
      geom_text(aes(label=fraction_reads_assigned), vjust=1.6, color="black", size=3.5)


  seq_plot2=ggplot(data=seqstats, aes(x=sample, y=reads_assigned_file, fill=sample)) + geom_bar(stat="identity") +
      scale_fill_viridis(discrete=TRUE, alpha=0.35,option="turbo") +
      theme_bw() +
      geom_text(aes(label=fraction_reads_assigned), vjust=1, hjust=1, angle=90, color="black", size=4) +
      theme(axis.text.x = element_text(angle = 90, vjust = 1, size = 10, hjust = 1), 
        axis.text.y = element_text(size=10)) +
      theme(axis.title.x = element_blank()) +
      labs(y="Reads assigned to sgRNA")

  tot_sgRNAs=seqstats[1,3]

  seqstats$sgRNA_file=format(as.numeric(seqstats$sgRNA_file), nsmall=0, big.mark=",")
  seqstats$reads_assigned_file=format(as.numeric(seqstats$reads_assigned_file), nsmall=0, big.mark=",")
  seqstats$reads_total=format(as.numeric(seqstats$reads_total), nsmall=0, big.mark=",")

  ##change
  seqstats=seqstats[,c(8,2,7,4,5,6)]
  colnames(seqstats)=c("sample","sgRNA detected","sgRNAs not detected","reads assigned","reads total","fraction reads assigned")


  }else{
    tot_sgRNAs="na"
  }
}

## ---- data_countingstats_reads
if(!is.RSL){
  file.summary.stats=file.path(ctable_datadir,paste0(proj.name.pref,".countsummary.txt"))
  summary.stats=read.delim(file.summary.stats, header = TRUE, sep = "\t", quote = "\"", dec = ".", fill = TRUE, row.names=NULL)
  summary.stats$sample=factor(summary.stats$Label, levels=samples.tab$library)
  summary.stats$detected=summary.stats$TotalsgRNAs-summary.stats$Zerocounts

  seq_plot2=ggplot(data=summary.stats.table, aes(x=sample, y=Mapped, fill=sample)) + geom_bar(stat="identity") +
      scale_fill_viridis(discrete=TRUE, alpha=0.35,option="turbo") +
      theme_bw() +
      geom_text(aes(label=Percentage), vjust=1, hjust=1, angle=90, color="black", size=4) +
      theme(axis.text.x = element_text(angle = 90, vjust = 1, size = 10, hjust = 1), 
        axis.text.y = element_text(size=10)) +
      theme(axis.title.x = element_blank()) +
      labs(y="Reads assigned to sgRNA")

  tot_sgRNAs=summary.stats[1,6]

  summary.stats$detected=format(as.numeric(summary.stats$detected), nsmall=0, big.mark=",")
  summary.stats$Reads=format(as.numeric(summary.stats$Reads), nsmall=0, big.mark=",")
  summary.stats$Mapped=format(as.numeric(summary.stats$Mapped), nsmall=0, big.mark=",")

  summary.stats.table=summary.stats[,c(14,15,7,4,3,5,8)]
  summary.stats.table=summary.stats.table[match(samples.tab$library, summary.stats.table$sample),]

  colnames(summary.stats.table)=c("sample","sgRNA detected","sgRNAs not detected","reads assigned","reads total","fraction reads assigned","Gini index")


}

## ---- data-countingstats-table
if(is.RSL){

if(seqstatsSTAUS){

knitr::kable(seqstats, row.names = FALSE, caption = "Summary statistics of sgRNA quantification.")
}
}

if(!is.RSL){
knitr::kable(summary.stats.table, row.names = FALSE, caption = "Summary statistics of sgRNA quantification.")
}



## ---- count-stats-plot
seq_plot2 + 
	theme(aspect.ratio = 1) + theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))  + 
	theme(legend.text=element_text(size=10),legend.title=element_text(size=0)) 


## ---- count-stats-plot-save
fig_n=fig_n+1
fname=paste("Figure",fig_n,"sgRNAassignStatistics.pdf",sep=".")
pdf(file.path(plotdir,fname))
seq_plot2 + 
  theme(aspect.ratio = 1) + theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))  + 
  theme(legend.text=element_text(size=10),legend.title=element_text(size=0)) 
dev.off()

## ---- data_countingstats-save
if(is.RSL){
if(seqstatsSTAUS){

fname="sgRNAassignStatistics.tsv"
save.table(df=seqstats, file=fname, dir=resdir)
}
}

## ---- reads-per-guide

if (data.type == "RSL"){
  
  readsperGuide_cnt.fname=file.path(ctable_datadir,paste(proj.name.pref,"RSL.perguide.tsv",sep="."))
  frequencies_reads.fname=file.path(ctable_datadir,paste(proj.name.pref,"frequencies.filt_reads.tsv",sep="."))
  ylab.txt="RSL counts"

  print(frequencies_reads.fname)
  print(readsperGuide_cnt.fname)

}else{
  readsperGuide_cnt.fname=file.path(ctable_datadir,paste(proj.name.pref,"count_normalized.txt",sep="."))
  readsperGuide_cnt_nonnorm.fname=file.path(ctable_datadir,paste(proj.name.pref,"count.txt",sep="."))
  ylab.txt="Read counts"

}



reads.pG=read.delim(readsperGuide_cnt.fname, header = TRUE, sep = "\t", quote = "\"",dec = ".", fill = TRUE, row.names=NULL)
dat_reads_h = gather(reads.pG[,-2], variable, value,-sgRNA)

dat_reads_h_ord=dat_reads_h
dat_reads_h_ord$variable=factor(dat_reads_h_ord$variable, levels=samples.tab$library[order(samples.tab$library)])


box_rawreads_perguide=ggplot(dat_reads_h_ord, aes(x=variable, y=value, color=variable))  + geom_boxplot()+
	theme_bw()+theme(legend.position="none")+
	scale_color_viridis(discrete=TRUE, alpha=0.35,option="turbo") +
	theme(axis.text.x=element_text(angle=90,hjust=1)) +
  	labs(y=paste(ylab.txt,"per guide",sep=" "),x="") +
  	coord_flip()

# box_rawreads_perguide2=ggplot(dat_reads_h, aes(x=variable, y=value))  + geom_boxplot()+
# 	theme_bw()+ scale_y_sqrt() + 
# 	theme(axis.text.x=element_text(angle=90,hjust=1,size = 10), axis.text.y = element_text(size=10)) +
#   	labs(y=paste("Square root of",ylab.txt,"per guide",sep=" "),x="") +
#   	coord_flip()


# box_rawreads_perguide3=ggplot(dat_reads_h, aes(x=variable, y=value))  + geom_boxplot()+
# 	theme_bw()+ scale_y_continuous(trans='log2')+ theme(legend.position="none")+
# 	theme(axis.text.x=element_text(angle=90,hjust=1,size = 10), axis.text.y = element_text(size=10)) +
#   	labs(y=paste("Log2 of",ylab.txt,"per guide",sep=" "),x="") +
#   	coord_flip()

# #outliers are removed
# box_rawreads_perguide4=ggplot(dat_reads_h, aes(x=variable, y=value))  + geom_boxplot(outlier.shape=NA)+
# 	theme_bw()+ theme(legend.position="none")+
# 	theme(axis.text.x=element_text(angle=90,hjust=1)) + scale_y_continuous(limits = quantile(dat_reads_h$value, c(0.1, 0.9))) +
#   	labs(y=paste(ylab.txt,"per guide (outliers are removed)",sep=" "),x="") +
#   	coord_flip()


# violin_rawreads_perguide=ggplot(dat_reads_h, aes(x=variable, y=value))  + geom_violin()+
# 	theme_bw()+ theme(legend.position="none")+
# 	theme(axis.text.x=element_text(angle=90,hjust=1)) +
#     labs(y=paste(ylab.txt,"per guide",sep=" "),x="") +
#   	coord_flip()

# #with box
# violin_rawreads_perguide2=ggplot(dat_reads_h, aes(x=variable, y=value,color=variable))  +  theme(legend.position="none")+
# 	geom_violin(width=0.8)+ geom_boxplot(width=0.1, color="black", alpha=0.2) +
# 	theme_bw()+ scale_color_viridis(discrete=TRUE, alpha=0.75,option="turbo") +
# 	theme(axis.text.x=element_text(angle=90,hjust=1)) +
#     labs(y=paste(ylab.txt,"per guide",sep=" "),x="") +
#   	coord_flip()


# log scale
box_rawreads_perguide_log=ggplot(dat_reads_h_ord, aes(x=variable, y=value, color=variable))  + geom_boxplot()+
  theme_bw()+theme(legend.position="none")+
  scale_y_continuous(trans='log10')+
  scale_color_viridis(discrete=TRUE, alpha=0.35,option="turbo") +
  theme(axis.text.x=element_text(angle=90,hjust=1)) +
    labs(y=paste("log10 (",ylab.txt,") per guide",sep=""),x="") +
    coord_flip()



#outliers removed
violin_rawreads_perguide3=ggplot(dat_reads_h_ord, aes(x=variable, y=value,color=variable))  + theme(legend.position="none") +
	geom_violin(width=0.8)+ geom_boxplot(width=0.1, color="black", alpha=0.2) +
	theme_bw()+ scale_color_viridis(discrete=TRUE, alpha=0.75,option="turbo") +
	scale_y_continuous(limits = quantile(dat_reads_h$value, c(0.1, 0.9))) +
	theme(axis.text.x=element_text(angle=90,hjust=1)) +
    labs(y=paste(ylab.txt,"per guide (outliers are removed)",sep=" "),x="") +
  	coord_flip()

# #log scale
# violin_rawreads_perguide4=ggplot(dat_reads_h, aes(x=variable, y=value,color=variable))  + 
# 	geom_violin(width=0.8)+ geom_boxplot(width=0.1, color="black", alpha=0.2) +
# 	theme_bw()+ scale_color_viridis(discrete = TRUE) +
# 	scale_y_continuous(trans='log2')+
# 	theme(axis.text.x=element_text(angle=90,hjust=1)) +
#   	labs(y="Log2 raw reads per guide",x="") +
#   	coord_flip()



## ---- reads-per-guide-boxplot
box_rawreads_perguide + theme(legend.position="none")
#box_rawreads_perguide3
#box_rawreads_perguide4

## ---- reads-per-guide-boxplot-log
box_rawreads_perguide_log + theme(legend.position="none")
#box_rawreads_perguide3
#box_rawreads_perguide4


## ---- reads-per-guide-violinplot
#violin_rawreads_perguide2 + theme(legend.position="none")
violin_rawreads_perguide3 + theme(legend.position="none") 



## ---- boxplots-save
fig_n=fig_n+1
fname=paste("Figure",fig_n,"ReadssgRNA.pdf",sep=".")
pdf(file.path(plotdir,fname))

box_rawreads_perguide + theme(legend.position="none")+ 
  theme(aspect.ratio = 1) + theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))

box_rawreads_perguide_log + theme(legend.position="none")+ 
  theme(aspect.ratio = 1) + theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))

violin_rawreads_perguide3 + theme(legend.position="none") + 
  theme(aspect.ratio = 1) + theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))

dev.off()


## ---- histogram-rslguide
freq_data=read.delim(frequencies_reads.fname, header = TRUE, sep = "\t", quote = "\"",dec = ".", fill = TRUE, row.names=NULL)

raw_freq_r = gather(freq_data, variable, value,-readcount)
colnames(raw_freq_r)=c("readcount","sample","frequency")


title_plot1=paste("Reads per RSL-guide after input-based filtering for RSL analysis.",sep=(" "))


plot1=ggplot(data=raw_freq_r, aes(x=readcount, y=frequency, colour=sample, fill=sample)) + geom_bar(stat="identity") +
  facet_wrap(~sample)+
    scale_fill_viridis(discrete=TRUE,option="turbo")+
    scale_colour_viridis(discrete=TRUE,option="turbo")+
  theme_bw() + theme(legend.position="bottom") +
  labs(y="count",x="reads per RSL-guide")


# x_lim=40
# y_lim=500000

# better lim specification
x_lim=40 
freqs.lim.i=raw_freq_r%>%filter(readcount>=1, readcount<=x_lim)
y_lim=round(max(freqs.lim.i$frequency)+0.1*(max(freqs.lim.i$frequency)))


plot2=ggplot(data=raw_freq_r, aes(x=readcount, y=frequency, colour=sample, fill=sample)) + geom_bar(stat="identity") +
  facet_wrap(~sample)+
    scale_fill_viridis(discrete=TRUE,option="turbo")+
    scale_colour_viridis(discrete=TRUE,option="turbo")+
  theme_bw() + theme(legend.position="bottom") +
  coord_cartesian(xlim=c(0,x_lim),ylim=c(0,y_lim) )+
  labs(y="count",x="reads per RSL-guide")


# fig_n=fig_n+1
# fname=paste("Figure",fig_n,"ReadCount_per_RSLguide.pdf",sep=".")
# pdf(file.path(plotdir,fname))

# plot1 + theme(legend.position="none")+ 
#   theme(aspect.ratio = 1) + theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")) +
#   labs(title=title_plot1, y="count",x="reads per RSL-guide")
# plot2 + theme(legend.position="none") + 
#   theme(aspect.ratio = 1) + theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")) +
#   labs(title=title_plot1, y="count",x="reads per RSL-guide")

# dev.off()


# detected guides table
total_sgRNA=nrow(reads.pG)

reads.pG.mat=as.matrix(reads.pG[,-c(1,2)])
detect_guides=colSums(reads.pG.mat != 0)

detguides_table=as.data.frame(cbind(names(detect_guides),detect_guides))
detguides_table=rbind(detguides_table, c("total present",total_sgRNA))
colnames(detguides_table)=c("sample","guides")
detguides_table=detguides_table[match(samples.tab$library, detguides_table$sample),]


## ---- histogram-rslguide-plot1
plot1
plot1=plot1 + theme(legend.position="none")+ 
  theme(aspect.ratio = 1) + theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")) +
  labs(title=title_plot1, y="count",x="reads per RSL-guide")
fig_n=fig_n+1
ggsave(filename=paste("Figure",fig_n,"ReadCount_per_RSLguide_full_scale.pdf",sep="."),path=plotdir,device="pdf")



## ---- histogram-rslguide-plot2
plot2
plot2=plot2 + theme(legend.position="none")+ 
  theme(aspect.ratio = 1) + theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")) +
  labs(title=title_plot1, y="count",x="reads per RSL-guide")
fig_n=fig_n+1
ggsave(filename=paste("Figure",fig_n,"ReadCount_per_RSLguide_redused_scale.pdf",sep="."),path=plotdir,device="pdf")


## ---- detected-guides-table
knitr::kable(detguides_table, row.names = FALSE, caption = "Number of detected guides after filtering for RSL analysis.")



#######################################################
####################### RRA


## ---- read-contrasts
contrasts.tab=read.table(comparisons.file, sep="\t", header=TRUE, blank.lines.skip=TRUE)
my.contrasts=contrasts.tab$name
my.contrasts=unique(my.contrasts)
n.cont=nrow(contrasts.tab)



## ---- contrasts-table
knitr::kable(contrasts.tab, row.names = FALSE, caption = "Comparisons analysed in this report.")


## ---- collect-rra-res
all.res.rra=vector(mode = "list", length = nrow(contrasts.tab))
all.res=vector(mode = "list", length = nrow(contrasts.tab))


## ---- contrasts-loop ## not in the report; instead this code is given directly in the chunk
res <- vector(mode = "list", length = n.cont)
options(knitr.duplicate.label = "allow")

for (i in my.contrasts) {
  res[[i]] <- knitr::knit_child("crispr_report_rra.Rmd", quiet = TRUE, envir = environment())
}

cat(unlist(res), sep = '\n')




## ---- pca_prep

samples.tab=read.table(samples.file, sep="\t", header=TRUE, blank.lines.skip=TRUE)
samples.tab=samples.tab[,-1]
colnames(samples.tab)[1]="library"


## PCA on scaled data  after log2 transformation + pseudocount

## reads: scaling by mageck
if(!is.RSL){
  reads.pG.unscaled=read.delim(readsperGuide_cnt_nonnorm.fname, header = TRUE, sep = "\t", quote = "\"",dec = ".", fill = TRUE, row.names=NULL)

  #sgRNA
  dat_pca=reads.pG.unscaled
  rownames(dat_pca)=dat_pca$sgRNA
  dat_pca=dat_pca[,-c(1,2)]


  pca.sG=plot_pca_TMM(dat_pca=dat_pca, annot=samples.tab)

  dge = DGEList(counts=dat_pca) # TMM normalisation only
  dge = calcNormFactors(dge)
  dat_pca_lognorm_sgRNA <- cpm(dge, normalized.lib.sizes = TRUE,log = TRUE, prior.count = 1)

  #genes
  df=reads.pG.unscaled

  data_umi_bygene= df %>%
  group_by(Gene) %>%
  summarise(across(where(is.numeric), sum))

  dat_pca=as.data.frame(data_umi_bygene)
  dat_pca$Gene[is.na(dat_pca$Gene)]="na"
  rownames(dat_pca)=dat_pca$Gene
  dat_pca=dat_pca[,-c(1)]

  pca.G=plot_pca_TMM(dat_pca=dat_pca, annot=samples.tab)

  dge = DGEList(counts=dat_pca) # TMM normalisation only
  dge = calcNormFactors(dge)
  dat_pca_lognorm_gene <- cpm(dge, normalized.lib.sizes = TRUE,log = TRUE, prior.count = 1)

}

## RSL: TMM scaling

if(is.RSL){
  #sgRNA
  dat_pca=reads.pG
  rownames(dat_pca)=dat_pca$sgRNA
  dat_pca=dat_pca[,-c(1,2)]

  pca.sG=plot_pca_TMM(dat_pca=dat_pca, annot=samples.tab)

  dge = DGEList(counts=dat_pca) # TMM normalisation only
  dge = calcNormFactors(dge)
  dat_pca_lognorm_sgRNA <- cpm(dge, normalized.lib.sizes = TRUE,log = TRUE, prior.count = 1)


  #genes
  df=reads.pG

  data_umi_bygene= df %>%
  group_by(Gene) %>%
  summarise(across(where(is.numeric), sum))

  dat_pca=as.data.frame(data_umi_bygene)
  dat_pca$Gene[is.na(dat_pca$Gene)]="na"
  rownames(dat_pca)=dat_pca$Gene
  dat_pca=dat_pca[,-c(1)]


  pca.G=plot_pca_TMM(dat_pca=dat_pca, annot=samples.tab)

  dge = DGEList(counts=dat_pca) # TMM normalisation only
  dge = calcNormFactors(dge)
  dat_pca_lognorm_gene <- cpm(dge, normalized.lib.sizes = TRUE,log = TRUE, prior.count = 1)
}


#get legend
mock_pca=pca.sG
legend <- cowplot::get_legend(mock_pca)


#pca.plot=cowplot::plot_grid(pca1,pca2,mock_pca, ncol=3, labels=LETTERS[1:2])

pca1=pca.sG+
    theme(aspect.ratio = 1) + theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))  + 
    theme(legend.text=element_text(size=10),legend.title=element_text(size=0)) +
    theme(legend.position="none")

pca2=pca.G+
    theme(aspect.ratio = 1) + theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))  + 
    theme(legend.text=element_text(size=10),legend.title=element_text(size=0)) +
    theme(legend.position="none")

plot_pca_combined=plot_grid(plot_grid(pca1, pca2, ncol=2, labels=LETTERS[1:2]),
                 plot_grid(legend, ncol=1), ncol=1,
                  rel_widths=c(1, 0.2))


fig_n=fig_n+1
ggsave(filename=paste("Figure",fig_n,"PCA.pdf",sep="."),path=plotdir,device="pdf", plot=plot_pca_combined)


## ---- pca-plot
plot_pca_combined



## ---- smpl_correlations
smpl_corr_sgRNA=cor(dat_pca_lognorm_sgRNA,use="pairwise.complete.obs", method="spearman")
smpl_corr_gene=cor(dat_pca_lognorm_gene,use="pairwise.complete.obs", method="spearman")


cor_hm_sgRNA=plot_corr_hm(corr_matirx=smpl_corr_sgRNA)
cor_hm_gene=plot_corr_hm(corr_matirx=smpl_corr_gene)


fig_n=fig_n+1
fname=paste("Figure",fig_n,"SpearmanCorrHeatmaps.pdf",sep=".")
pdf(file.path(plotdir,fname))

cor_hm_sgRNA +theme(aspect.ratio = 1) + theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))  + 
  theme(legend.text=element_text(size=8),legend.title=element_text(size=0))

cor_hm_sgRNA +theme(aspect.ratio = 1) + theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))  + 
  theme(legend.text=element_text(size=8),legend.title=element_text(size=0))

dev.off()


## ---- smpl-correlations-gene
cor_hm_gene +theme(aspect.ratio = 1) + theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))  + 
  theme(legend.text=element_text(size=8),legend.title=element_text(size=0))


## ---- smpl-correlations-sgRNA
cor_hm_sgRNA +theme(aspect.ratio = 1) + theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))  + 
  theme(legend.text=element_text(size=8),legend.title=element_text(size=0))




## ---- contrast_scatters


# here check if it is a legit pair for scatter plotting
groups_to_plot=unique(contrasts.tab$group)

pairs_to_plot=list()


for (g in groups_to_plot){
  comparisons_tab=contrasts.tab[contrasts.tab$group==g,]
  comparisons=comparisons_tab$name

  if( length(comparisons) >=2 ){

    comparison_pairs_g=combn(comparisons,2)

    contrasts.pairs.number=ncol(comparison_pairs_g)

    for (i in c(1:contrasts.pairs.number)){

        contr.pair.i=comparison_pairs_g[,i]

        name.gi=paste(g,i,sep=".")

        pairs_to_plot[[name.gi]]=contr.pair.i
      }
  }
}




#this has been defined before, has to be unique now
#my.contrasts=contrasts.tab$name


#contrasts combinations
contrasts.pairs.mtx=combn(my.contrasts,2)

contrasts.pairs.number=ncol(contrasts.pairs.mtx)

#produce the plots
#scatters=vector(mode = "list", length = contrasts.pairs.number)
scatters=vector(mode = "list")


for (i in c(1:contrasts.pairs.number)){

#for (i in c(1:pairs_to_plot)){

  contr.pair.i=contrasts.pairs.mtx[,i]

  if ( any(pairs_to_plot %in% list(contr.pair.i)) ){

    res.df=data.frame()
    for (j in c( contr.pair.i[1], contr.pair.i[2])){
   
      if(!is.RSL){
        res.df.i=as.data.frame(cbind(all.res[[j]]$neg.lfc,all.res[[j]]$id ))
        colnames(res.df.i)=c("lfc","id")
        res.df.i$comparison=j
        res.df=rbind(res.df,res.df.i)
      }
      if(is.RSL){
        res.df.i=as.data.frame(cbind(all.res[[j]]$median.logFC,all.res[[j]]$id ))
        colnames(res.df.i)=c("lfc","id")
        res.df.i$comparison=j
        res.df=rbind(res.df,res.df.i)
      }
    }

    res.df$lfc=as.numeric(res.df$lfc)
    df_scatter=spread(res.df,comparison,lfc)

    #add this to avoid error in density distribution
    #https://stackoverflow.com/questions/53075331/error-using-geom-density-2d-in-r-computation-failed-in-stat-density2d-b
    #Error in MASS::kde2d(x, y, ...) : 
    #  missing or infinite values in the data are not allowed
    # OBS! this still does not work for RSL data, so will only be run for reads data
    pseudocount=0.01
    df_scatter[,4]=df_scatter[,2]+pseudocount
    df_scatter[,5]=df_scatter[,3]+pseudocount

    df_scatter$density <- get_density(df_scatter[,2], df_scatter[,3], n = 100)

    pl2=ggplot(df_scatter, aes(x=df_scatter[,4],y=df_scatter[,5], text=paste(id,"; lfc",colnames(df_scatter)[2],round(df_scatter[,2],digits=3),"; lfc",colnames(df_scatter)[3],round(df_scatter[,3],digits=3)) ))+
      geom_point() +
      theme_bw() +
       xlab(paste0("log2 Fold Change in ",colnames(df_scatter)[2])) + ylab(paste0("log2 Fold Change in ",colnames(df_scatter)[3]))


    pl2.d=pl2+geom_point(aes(df_scatter[,4], df_scatter[,5], color = density)) + scale_color_viridis()

    fig_n=fig_n+1
    ggsave(filename=paste("Figure",fig_n,"log2FCscatterplot",colnames(df_scatter)[2],colnames(df_scatter)[3],"comparison_group",g,i,"pdf",sep="."),path=plotdir,device="pdf")


    pl2_int=ggplotly(pl2.d, tooltip=c("text"))

    scatters[[i]]=pl2_int

  }
}


## ---- contrast-scatters-plot

for (i in c(1:contrasts.pairs.number)){
  print(scatters[[i]])
}


## ---- nn

