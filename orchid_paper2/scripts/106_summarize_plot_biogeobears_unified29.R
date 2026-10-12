#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(BioGeoBEARS); library(ape)})
root <- normalizePath(Sys.getenv("ORCHID_PAPER2_ROOT"), mustWork=TRUE)
base_env <- Sys.getenv("BGB_BASE",unset="")
base <- if(nzchar(base_env)) normalizePath(base_env,mustWork=TRUE) else
 file.path(root,"orchid_paper2/results/biogeobears_unified29")
max_range_size <- as.integer(readLines(file.path(base,"inputs","max_range_size.txt"),n=1))
dataset_label <- if(file.exists(file.path(base,"inputs","dataset_label.txt")))
 trimws(readLines(file.path(base,"inputs","dataset_label.txt"),n=1)) else "Unified 29-tip BioGeoBEARS"
dataset_tag <- if(file.exists(file.path(base,"inputs","dataset_tag.txt")))
 trimws(readLines(file.path(base,"inputs","dataset_tag.txt"),n=1)) else "unified29"
figdir <- file.path(base,"figures"); dir.create(figdir,recursive=TRUE,showWarnings=FALSE)
files <- list.files(file.path(base,"runs"), pattern="[.]summary[.]tsv$", recursive=TRUE, full.names=TRUE)
if (length(files) != 606L) stop("Expected 606 completed summaries, found ",length(files))
x <- do.call(rbind,lapply(files,read.delim,check.names=FALSE))
x$AIC <- -2*x$log_likelihood + 2*x$n_parameters
x$AICc <- x$AIC + 2*x$n_parameters*(x$n_parameters+1)/(x$n_tips-x$n_parameters-1)
x$at_lower_bound <- x$d <= 1.01e-12 | x$e <= 1.01e-12
x$tree_order <- ifelse(x$tree_id=="mcc",0,as.integer(x$tree_id))
x <- x[order(x$tree_order,x$AICc),]
x$delta_AICc <- ave(x$AICc,x$tree_id,FUN=function(z) z-min(z))
x$AICc_weight <- ave(x$delta_AICc,x$tree_id,FUN=function(z){w<-exp(-z/2);w/sum(w)})
write.table(x,file.path(base,"all_model_fits.tsv"),sep="\t",quote=FALSE,row.names=FALSE)

winners <- do.call(rbind,lapply(split(x,x$tree_id),function(z) z[which.min(z$AICc),]))
write.table(winners,file.path(base,"tree_model_winners.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
postwin <- subset(winners,tree_id!="mcc")
freq <- as.data.frame(table(postwin$model),stringsAsFactors=FALSE)
names(freq) <- c("model","winning_tree_count"); freq$winning_tree_fraction <- freq$winning_tree_count/100
write.table(freq,file.path(base,"posterior_model_winner_frequencies.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
mcc <- subset(x,tree_id=="mcc"); best <- mcc$model[which.min(mcc$AICc)]
writeLines(best,file.path(base,"selected_mcc_model.txt"))

load(file.path(base,"runs","mcc",paste0(best,".RData")))
tr <- read.tree(file.path(base,"inputs","mcc_myr.newick"))
# Pruning PEQ leaves the former PEQ-to-ingroup stem as root.edge. It is not
# part of the 29-tip crown chronogram and must not be drawn as extra time.
tr$root.edge <- NULL
tipranges <- getranges_from_LagrangePHYLIP(lgdata_fn=file.path(base,"inputs","geography.data"))
scriptdir <- normalizePath(system.file("extdata/a_scripts",package="BioGeoBEARS"))
probs <- as.matrix(res$ML_marginal_prob_each_state_at_branch_top_AT_node)
popmeta <- read.delim(file.path(root,"orchid_paper2/17_herodotools/inputs/bioregions/tip_to_thesis_bioregion.tsv"),
 stringsAsFactors=FALSE)
species_by_tip <- setNames(ifelse(popmeta$lineage=="PAM","PAM","PAP"),popmeta$tip)
species_by_tip["MDA"] <- "PSA"
species_cols <- c(PAM="#D55E00",PAP="#009E73",PSA="#7B2CBF")
area_order <- read.delim(file.path(base,"inputs","area_order.tsv"),stringsAsFactors=FALSE)
region_ids <- area_order$region_id
global_region_cols <- c(taiwan="#E41A1C",north_luzon="#FF7F00",south_luzon="#F0E442",
 greater_negros_panay="#A65628",mindoro="#4DAF4A",palawan="#00BFC4",
 mindanao="#56B4E9",borneo="#0057D9",java="#6A00FF",
 sulawesi_moluccas_lesser_sunda="#F000C8",new_guinea_queensland="#FF3B9D")
region_cols <- global_region_cols[region_ids]
if(anyNA(region_cols)) stop("Missing region color")
region_by_tip <- setNames(popmeta$region_id,popmeta$tip)
# BioGeoBEARS' default state order is null followed by all combinations of
# one through four of the eleven areas. Mix constituent area colors for
# multi-area ranges, retaining distinct colors for single-area states.
if(!is.null(res$inputs$states_list)) {
 states <- lapply(res$inputs$states_list,function(s) {
   if(!length(s) || all(is.na(s))) return(integer(0))
   as.integer(s)+1L
 })
} else {
 states <- list(integer(0))
 for(k in seq_len(max_range_size)) states <- c(states,combn(seq_along(region_ids),k,simplify=FALSE))
}
stopifnot(length(states)==ncol(probs))
state_cols <- vapply(states,function(s) {
 if(length(s)==0) return("#222222")
 z <- col2rgb(region_cols[region_ids[s]])/255
 rgb(rowMeans(z)[1],rowMeans(z)[2],rowMeans(z)[3])
},character(1))
sea <- read.delim(file.path(root,"orchid_paper2/results/biogeobears/sea_level/deboer2014-noaa.txt"),comment.char="#")
sea <- sea[sea$age_calkaBP <= max(node.depth.edgelength(tr))*1000,]
sea$Ma <- sea$age_calkaBP/1000
# One-kyr plotting resolution preserves Quaternary cycles without raster overdraw.
sea <- sea[seq(1,nrow(sea),by=10),]
write.table(sea[,c("Ma","sealev")],file.path(base,"sea_level_deboer2014_plot_data.tsv"),sep="\t",quote=FALSE,row.names=FALSE)

draw_timescale <- function(root_age, fine=FALSE) {
  if (!fine) {
    b <- data.frame(start=c(5.333,2.58,0.0117),end=c(2.58,0.0117,0),
      label=c("Pliocene","Pleistocene","Holocene"),col=c("#F4D35E","#B7E4C7","#FFF3B0"))
  } else {
    b <- data.frame(start=c(3.60,2.58,1.80,0.774,0.129,0.0117),end=c(2.58,1.80,0.774,0.129,0.0117,0),
      label=c("Piacenzian","Gelasian","Calabrian","Chibanian","Late Pleist.","Holocene"),
      col=c("#F6BD60","#A8DADC","#80B1D3","#B3CDE3","#CCEBC5","#FFF3B0"))
  }
  plot(NA,xlim=c(root_age,0),ylim=c(0,1),axes=FALSE,xlab="",ylab="",xaxs="i",yaxs="i")
  for(i in seq_len(nrow(b))) {
    lo <- max(0,b$end[i]); hi <- min(root_age,b$start[i]); if(hi<=lo) next
    rect(hi,0,lo,1,col=b$col[i],border="white")
    if(abs(hi-lo)>0.08*root_age) text((hi+lo)/2,.5,b$label[i],cex=1.08)
  }
  if(fine) axis(1,at=pretty(c(0,root_age)),labels=pretty(c(0,root_age)),cex.axis=1.08,mgp=c(1.4,.35,0),tcl=-.25)
}

draw <- function() {
  root_age <- max(node.depth.edgelength(tr))
  nt <- length(tr$tip.label)
  # A dedicated label column keeps the time-zero boundary identical in the
  # tree, sea-level and geological panels while leaving room for tip boxes.
  layout(matrix(c(1,2,3,0,4,0,5,0),nrow=4,byrow=TRUE),
    widths=c(5.0,1.18),heights=c(9.3,2.0,.72,.82))
  par(mar=c(.25,5.2,2.8,0),xpd=FALSE)
  plot(tr,show.tip.label=FALSE,direction="rightwards",use.edge.length=TRUE,
    x.lim=c(0,root_age),xaxs="i",y.lim=c(.5,length(tr$tip.label)+.5),
    edge.width=.85,edge.color="#333333",no.margin=FALSE)
  pp <- get("last_plot.phylo",envir=.PlotPhyloEnv)
  tree_ylim <- par("usr")[3:4]
  inn <- (nt+1):(nt+tr$Nnode)
  # Extend every terminal branch to the common present boundary so its
  # dashed leader, region box and population name form one aligned row.
  segments(pp$xx[seq_len(nt)],pp$yy[seq_len(nt)],root_age,pp$yy[seq_len(nt)],
    lty=2,lwd=.65,col="grey45")
  par(xpd=NA)
  nodelabels(node=inn,pie=probs[inn,,drop=FALSE],piecol=state_cols,
    cex=.48,frame="none")
  par(xpd=FALSE)
  title(main=paste0(dataset_label,": ",best),cex.main=1.28,line=1.0)

  par(mar=c(.25,0,2.8,.4),xpd=FALSE)
  # Reuse the exact plotted y range from plot.phylo. This avoids the subtle
  # vertical stretching caused by constructing an independent 1:nt axis.
  plot(NA,xlim=c(0,1),ylim=tree_ylim,xaxs="i",yaxs="i",axes=FALSE,xlab="",ylab="")
  tip_y <- pp$yy[seq_len(nt)]
  tip_prob <- probs[seq_len(nt),,drop=FALSE]
  tip_state <- max.col(tip_prob,ties.method="first")
  for(i in seq_len(nt)) {
    segments(0,tip_y[i],.12,tip_y[i],lty=2,lwd=.85,col="grey40")
    rect(.12,tip_y[i]-.19,.19,tip_y[i]+.19,
      col=region_cols[region_by_tip[tr$tip.label[i]]],border="black",lwd=.7)
    sp <- species_by_tip[tr$tip.label[i]]
    text(.215,tip_y[i],tr$tip.label[i],adj=c(0,.5),cex=1.02,font=2,col=species_cols[sp])
  }

  par(mar=c(.2,4.8,.7,.25),xpd=FALSE)
  plot(sea$Ma,sea$sealev,type="l",xlim=c(root_age,0),xaxs="i",xaxt="n",
    col="#1769AA",lwd=.8,xlab="",ylab="Sea level (m) relative to present",
    cex.lab=1.12,cex.axis=1.05)
  abline(h=0,lty=2,col="grey45",lwd=.7)
  mtext("de Boer et al. (2014) global reconstruction",side=3,adj=0,cex=1.02,line=.08)
  par(mar=c(0,4.8,0,.25),xpd=FALSE)
  draw_timescale(root_age,FALSE)
  par(mar=c(4.5,5.2,0,.6),xpd=FALSE)
  draw_timescale(root_age,TRUE)
  mtext("Millions of years before present",side=1,line=2.45,cex=1.18)
}

draw_color_legend <- function() {
  max_prob <- apply(probs,2,max,na.rm=TRUE)
  single_idx <- which(lengths(states)==1L)
  # Include every multi-area state that contributes at least 1% to any pie.
  multi_idx <- which(lengths(states)>1L & max_prob>=0.01)
  shown <- c(single_idx,multi_idx)
  state_code <- vapply(states,function(s) {
    if(!length(s)) return("null")
    paste(area_order$code[s],collapse=" + ")
  },character(1))
  state_regions <- vapply(states,function(s) {
    if(!length(s)) return("null range")
    paste(area_order$region_name[s],collapse=" + ")
  },character(1))
  legend_table <- data.frame(state_index=shown,state_code=state_code[shown],
    regions=state_regions[shown],color=state_cols[shown],
    maximum_node_probability=max_prob[shown],
    state_type=ifelse(lengths(states)[shown]==1L,"single region","multiple-region combination"))
  write.table(legend_table,file.path(figdir,paste0(dataset_tag,"_biogeobears_state_color_legend.tsv")),
    sep="\t",quote=FALSE,row.names=FALSE)

  plot.new(); plot.window(xlim=c(0,1),ylim=c(0,1),xaxs="i",yaxs="i")
  title(main="BioGeoBEARS state and species color key",cex.main=1.55,line=.35)
  text(.02,.93,"Population-name colors",adj=c(0,.5),font=2,cex=1.18)
  sx <- c(.22,.38,.54)
  points(sx,rep(.93,3),pch=15,col=species_cols,cex=1.7)
  text(sx+.018,rep(.93,3),names(species_cols),adj=c(0,.5),cex=1.15,
    col=species_cols,font=2)
  text(.02,.865,"Single-region colors (letters are used below for combinations)",
    adj=c(0,.5),font=2,cex=1.12)
  for(i in seq_along(single_idx)) {
    rr <- (i-1L) %/% 3L; cc <- (i-1L) %% 3L
    xx <- .02+cc*.325; yy <- .81-rr*.065
    rect(xx,yy-.016,xx+.022,yy+.016,col=state_cols[single_idx[i]],border="black")
    text(xx+.029,yy,paste0(state_code[single_idx[i]]," = ",state_regions[single_idx[i]]),
      adj=c(0,.5),cex=.92)
  }
  yhead <- .53
  text(.02,yhead,"Multiple-region colors visible in ancestral pies (maximum probability >= 0.01)",
    adj=c(0,.5),font=2,cex=1.12)
  if(length(multi_idx)) {
    nr <- ceiling(length(multi_idx)/3)
    dy <- min(.042,.43/max(1,nr))
    for(i in seq_along(multi_idx)) {
      rr <- (i-1L) %% nr; cc <- (i-1L) %/% nr
      xx <- .02+cc*.325; yy <- yhead-.055-rr*dy
      rect(xx,yy-.012,xx+.020,yy+.012,col=state_cols[multi_idx[i]],border="black")
      text(xx+.027,yy,state_code[multi_idx[i]],adj=c(0,.5),cex=.82)
    }
  } else text(.02,yhead-.06,"None",adj=c(0,.5),cex=1)
}
pdf(file.path(figdir,paste0(dataset_tag,"_biogeobears_pies_geology_sealevel.pdf")),width=15,height=14,useDingbats=FALSE)
draw(); dev.off()
png(file.path(figdir,paste0(dataset_tag,"_biogeobears_pies_geology_sealevel.png")),width=4500,height=4200,res=300)
draw(); dev.off()
pdf(file.path(figdir,paste0(dataset_tag,"_biogeobears_state_color_legend.pdf")),width=15,height=10,useDingbats=FALSE)
draw_color_legend(); dev.off()
png(file.path(figdir,paste0(dataset_tag,"_biogeobears_state_color_legend.png")),width=4500,height=3000,res=300)
draw_color_legend(); dev.off()

writeLines(c(
  paste("Selected MCC model:",best),
  paste("MCC AICc:",format(min(mcc$AICc),digits=8)),
  paste("Posterior-tree winners:",paste(freq$model,freq$winning_tree_count,sep="=",collapse=", ")),
  paste("Maximum ancestral range size:",max_range_size),
  paste0("Tips: ",Ntip(tr),"; regions: ",length(region_ids),"; posterior trees: 100")
),file.path(base,"analysis_summary.txt"))
