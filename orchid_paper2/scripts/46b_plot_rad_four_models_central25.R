#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(data.table);library(ggplot2)})
args=commandArgs(trailingOnly=TRUE);base=if(length(args))args[[1]]else"orchid_paper2"
d=file.path(base,"results","rad_locus_statistics");pdir=file.path(d,"plots");dir.create(pdir,recursive=TRUE,showWarnings=FALSE)
x=fread(file.path(d,"rad_locus_statistics.four_models_central25_dxy.tsv"));s=fread(file.path(d,"four_model_central25_dxy_summary.tsv"))
lev=c("divergence_with_gene_flow","allopatric_selection","recurrent_selection","balancing_selection","unclassified")
lab=c("Divergence with gene flow","Allopatric selection","Recurrent selection","Balancing selection","Unclassified")
pal=c("#D73027","#E69F00","#0072B2","#7B3294","grey78")
x[,model:=factor(four_model_class_central25_dxy,levels=lev,labels=lab)]
x[,pair:=factor(pair,levels=c("PAM-PAP","PAM-PSA","PAP-PSA"))]
p1=ggplot(x,aes(dxy,hudson_fst,color=model))+geom_point(shape=16,size=1.05,alpha=.72)+facet_wrap(~pair,scales="free")+
 scale_color_manual(values=setNames(pal,lab),name=NULL)+labs(x="dXY",y="Hudson FST",title="Four-model classification with central-25% allopatric dXY",subtitle="Allopatric selection requires dXY between Q37.5 and Q62.5; FST and species-specific pi retain quartile criteria")+
 theme_bw(base_size=10)+theme(panel.grid.minor=element_blank(),legend.position="top")
ggsave(file.path(pdir,"rad_four_models_central25_dxy_fst_dxy_scatter.pdf"),p1,width=10.5,height=4.5,useDingbats=FALSE,bg="white")
ggsave(file.path(pdir,"rad_four_models_central25_dxy_fst_dxy_scatter.png"),p1,width=10.5,height=4.5,dpi=240,bg="white")
s[,model:=factor(model,levels=lev,labels=lab)];s[,pair:=factor(pair,levels=c("PAM-PAP","PAM-PSA","PAP-PSA"))]
p2=ggplot(s[model!="Unclassified"],aes(pair,fraction_all_loci,fill=model))+geom_col(width=.7)+scale_fill_manual(values=setNames(pal[1:4],lab[1:4]),name=NULL)+scale_y_continuous(labels=function(z)sprintf("%d%%",round(100*z)))+labs(x=NULL,y="Fraction of all retained loci",title="Four-model assignments with central-25% allopatric dXY",subtitle="Unclassified loci remain in the denominator")+theme_bw(base_size=10)+theme(panel.grid.minor=element_blank(),legend.position="top",axis.text.x=element_text(angle=20,hjust=1))
ggsave(file.path(pdir,"rad_four_models_central25_dxy_proportions.pdf"),p2,width=9,height=5,useDingbats=FALSE,bg="white")
ggsave(file.path(pdir,"rad_four_models_central25_dxy_proportions.png"),p2,width=9,height=5,dpi=240,bg="white")
p3=ggplot(s[model!="Unclassified"],aes(pair,fraction_classified_loci,fill=model))+geom_col(width=.7)+scale_fill_manual(values=setNames(pal[1:4],lab[1:4]),name=NULL)+scale_y_continuous(labels=function(z)sprintf("%d%%",round(100*z)))+labs(x=NULL,y="Fraction of classified loci",title="Composition of central-25% dXY-classified RAD loci",subtitle="Allopatric dXY occupies Q37.5-Q62.5")+theme_bw(base_size=10)+theme(panel.grid.minor=element_blank(),legend.position="top",axis.text.x=element_text(angle=20,hjust=1))
ggsave(file.path(pdir,"rad_four_models_central25_dxy_composition_classified.pdf"),p3,width=9,height=5,useDingbats=FALSE,bg="white")
ggsave(file.path(pdir,"rad_four_models_central25_dxy_composition_classified.png"),p3,width=9,height=5,dpi=240,bg="white")
