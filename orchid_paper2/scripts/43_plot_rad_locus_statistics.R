#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(data.table); library(ggplot2)})
options(bitmapType="cairo")
args <- commandArgs(trailingOnly=TRUE)
base <- if (length(args)) args[[1]] else getwd()
indir <- file.path(base,"results","rad_locus_statistics")
outdir <- file.path(indir,"plots"); dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
stats <- fread(file.path(indir,"rad_locus_statistics.tsv"),na.strings="NA")
snps <- fread(cmd=paste("gzip -cd",shQuote(file.path(indir,"rad_snp_plot_table.tsv.gz"))))

# Nuclear landscape only. Species colors are fixed across every pair and plot.
stats <- stats[chrom!="NC_017609.1"]
snps <- snps[chrom!="NC_017609.1"]
pair_levels <- c("PAM-PAP","PAM-PSA","PAP-PSA")
metric_levels <- c("FST","dXY","pi","Tajima's D")
species_cols <- c(PAM="#0072B2",PAP="#D55E00",PSA="#009E73")

# Fragmented-reference Manhattan coordinate: occupied scaffolds in first-
# appearance order with fixed visual gaps. It is not a chromosome coordinate.
ct <- unique(snps[,.(chrom)]); ct[,contig_index:=.I]
span <- snps[,.(lo=min(pos),hi=max(pos)),by=chrom][ct,on="chrom"]
span[,width:=pmax(hi-lo+1L,1L)]; span[,offset:=shift(cumsum(width+1000),fill=0)]
snps <- span[snps,on="chrom"]; snps[,genome_pos:=offset+(pos-lo)]

# Make one pairwise series for FST/dXY and two species-specific series for pi/D.
make_locus_long <- function(d) {
 out <- vector("list",nrow(d)*6L); k <- 0L
 for(i in seq_len(nrow(d))) {
  z <- d[i]; taxa <- strsplit(z$pair,"-",fixed=TRUE)[[1]]
  vals <- list(c("FST","pair statistic",z$hudson_fst),
               c("dXY","pair statistic",z$dxy),
               c("pi",taxa[1],z$pi_group1),c("pi",taxa[2],z$pi_group2),
               c("Tajima's D",taxa[1],z$tajima_d_group1),c("Tajima's D",taxa[2],z$tajima_d_group2))
  for(a in vals) {k<-k+1L;out[[k]]<-data.table(pair=z$pair,locus_id=z$locus_id,metric=a[1],series=a[2],value=as.numeric(a[3]))}
 }
 rbindlist(out[seq_len(k)])
}
locus_long <- make_locus_long(stats)
locus_long[,metric:=factor(metric,levels=metric_levels)]
locus_long[,pair:=factor(pair,levels=pair_levels)]

# Each species-specific pi/D distribution gets its own lower-5% threshold.
# Pairwise FST/dXY retain upper-5% thresholds.
thresholds <- locus_long[,.(q05=quantile(value,.05,na.rm=TRUE),q95=quantile(value,.95,na.rm=TRUE)),by=.(pair,metric,series)]
thresholds[,tail:=ifelse(metric %in% c("FST","dXY"),"upper","lower")]
thresholds[,threshold:=ifelse(tail=="upper",q95,q05)]
locus_long <- thresholds[,.(pair,metric,series,tail,threshold)][locus_long,on=.(pair,metric,series)]
locus_long[,outlier:=fifelse(tail=="upper",value>=threshold,value<=threshold,na=FALSE)]

# Repeat each locus value at all SNPs belonging to that RAD locus.
coords <- unique(snps[,.(pair,locus_id,chrom,pos,genome_pos)])
man <- coords[locus_long,on=.(pair,locus_id),allow.cartesian=TRUE,nomatch=0]
man[,panel:=factor(paste(pair,metric,sep=" | "),levels=as.vector(t(outer(pair_levels,metric_levels,paste,sep=" | "))))]
man[,display:=fifelse(series=="pair statistic",fifelse(outlier,"Directional 5%","Pair statistic"),series)]
man[,point_alpha:=ifelse(outlier,1,.48)]
plot_cols <- c("Pair statistic"="grey25","Directional 5%"="#CC0000",species_cols)
thresholds[,panel:=factor(paste(pair,metric,sep=" | "),levels=levels(man$panel))]
thresholds[,display:=ifelse(series=="pair statistic","Directional 5%",series)]

p <- ggplot(man,aes(genome_pos,value))+
 geom_point(aes(color=display,alpha=point_alpha),shape=16,size=.60,stroke=0,na.rm=TRUE)+
 geom_hline(data=thresholds,aes(yintercept=threshold,color=display),linetype="dashed",linewidth=.4,show.legend=FALSE)+
 facet_wrap(~panel,ncol=4,scales="free_y")+
 scale_color_manual(values=plot_cols,name=NULL)+scale_alpha_identity()+
 labs(x="Occupied scaffolds in reference order",y="Locus statistic",title="RAD-locus population-genetic statistics",subtitle="Independent panel scales; upper 5% FST/dXY and species-specific lower 5% pi/Tajima's D")+
 theme_bw(base_size=9)+theme(panel.grid.minor=element_blank(),panel.grid.major.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),legend.position="top",strip.background=element_rect(fill="grey92"))
ggsave(file.path(outdir,"rad_locus_statistics_top5_manhattan.pdf"),p,width=14,height=8.5,useDingbats=FALSE,bg="white")
ggsave(file.path(outdir,"rad_locus_statistics_top5_manhattan.png"),p,width=14,height=8.5,dpi=240,bg="white")

# Violin distributions use the same species-specific pi and Tajima's D values.
locus_long[,xgroup:=interaction(pair,series,drop=TRUE,lex.order=TRUE)]
summ <- locus_long[!is.na(value),.(n=.N,mean=mean(value),sd=sd(value),median=median(value),threshold=unique(threshold)),by=.(pair,metric,series)]
pv <- ggplot(locus_long[!is.na(value)],aes(pair,value,group=xgroup,fill=series,color=series))+
 geom_violin(position=position_dodge(width=.82),scale="width",trim=TRUE,linewidth=.3,alpha=.58,na.rm=TRUE)+
 geom_boxplot(position=position_dodge(width=.82),width=.16,outlier.shape=NA,linewidth=.25,alpha=.75,na.rm=TRUE)+
 facet_wrap(~metric,ncol=2,scales="free_y")+
 scale_fill_manual(values=c("pair statistic"="grey55",species_cols),name=NULL)+
 scale_color_manual(values=c("pair statistic"="grey25",species_cols),name=NULL)+
 labs(x=NULL,y="Locus statistic",title="RAD-locus genetic-index distributions",subtitle="Independent panel scales; pi and Tajima's D are shown separately for each species")+
 theme_bw(base_size=10)+theme(panel.grid.minor=element_blank(),axis.text.x=element_text(angle=25,hjust=1),legend.position="top",strip.background=element_rect(fill="grey92"))
ggsave(file.path(outdir,"rad_locus_statistics_violins.pdf"),pv,width=10,height=8,useDingbats=FALSE,bg="white")
ggsave(file.path(outdir,"rad_locus_statistics_violins.png"),pv,width=10,height=8,dpi=240,bg="white")

fwrite(summ,file.path(indir,"rad_locus_statistics_summary.tsv"),sep="\t")
fwrite(thresholds[,.(pair,metric,series,q05,q95,tail,threshold)],file.path(indir,"top5_thresholds_species_specific.tsv"),sep="\t")
