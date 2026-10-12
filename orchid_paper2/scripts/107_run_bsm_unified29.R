#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(BioGeoBEARS); library(ape)})
root <- normalizePath(Sys.getenv("ORCHID_PAPER2_ROOT"), mustWork=TRUE)
base_env <- Sys.getenv("BGB_BASE",unset="")
base <- if(nzchar(base_env)) normalizePath(base_env,mustWork=TRUE) else
 file.path(root,"orchid_paper2/results/biogeobears_unified29")
out <- file.path(base,"bsm_100"); dir.create(out,recursive=TRUE,showWarnings=FALSE)
model <- trimws(readLines(file.path(base,"selected_mcc_model.txt"),n=1))
load(file.path(base,"runs","mcc",paste0(model,".RData")))
tipranges <- getranges_from_LagrangePHYLIP(lgdata_fn=file.path(base,"inputs","geography.data"))
areanames <- names(tipranges@df)

if (file.exists(file.path(out,"bsm_100_complete.RData"))) {
  load(file.path(out,"bsm_100_complete.RData"))
} else {
  inputs <- get_inputs_for_stochastic_mapping(res=res)
  save(inputs,file=file.path(out,"stochastic_mapping_inputs.RData"))
  set.seed(20260807)
  bsm <- runBSM(res,stochastic_mapping_inputs_list=inputs,maxnum_maps_to_try=200,
   nummaps_goal=100,maxtries_per_branch=40000,save_after_every_try=TRUE,
   savedir=out,seedval=20260807,wait_before_save=0.01,master_nodenum_toPrint=0)
  save(bsm,file=file.path(out,"bsm_100_complete.RData"))
}

src <- simulate_source_areas_ana_clado(res,bsm$RES_clado_events_tables,
 bsm$RES_ana_events_tables,areanames)
clado <- src$clado_events_tables; ana <- src$ana_events_tables
counts <- count_ana_clado_events(clado,ana,areanames,areanames)
save(counts,src,file=file.path(out,"bsm_event_counts.RData"))

# Write every canonical two-dimensional count/direction table returned by BioGeoBEARS.
for(nm in names(counts)) {
 obj <- counts[[nm]]
 if(length(dim(obj))==2L && length(obj)>0L)
   write.table(obj,file.path(out,paste0(nm,".tsv")),sep="\t",quote=FALSE,col.names=NA)
}

# Independent per-map audit using the actual event histories.
clado_levels <- c("founder event (j)","sympatry (y)","subset sympatry (s)","vicariance (v)")
ana_levels <- c("a","d","e")
permap <- lapply(seq_along(clado),function(i) {
 ctab <- clado[[i]]
 atab <- ana[[i]]
 cvals <- if(is.data.frame(ctab) && "clado_event_type" %in% names(ctab)) ctab$clado_event_type else character(0)
 avals <- if(is.data.frame(atab) && "event_type" %in% names(atab)) atab$event_type else character(0)
 ct <- table(factor(cvals,levels=clado_levels))
 at <- table(factor(avals,levels=ana_levels))
 data.frame(map=i,j=unname(ct[1]),y=unname(ct[2]),s=unname(ct[3]),v=unname(ct[4]),
   a=unname(at[1]),d=unname(at[2]),e=unname(at[3]))
})
permap <- do.call(rbind,permap)
write.table(permap,file.path(out,"event_counts_per_map.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
labels <- c(j="Founder-event/jump dispersal",a="Range-switching dispersal",
 d="Range-expansion dispersal",e="Local extinction/range contraction",
 v="Vicariance",y="Narrow sympatry",s="Subset sympatry")
cols <- names(labels)
# Use BioGeoBEARS' canonical summary as the authoritative event totals. This
# avoids dependence on version-specific text labels in raw cladogenetic tables.
canon <- counts$summary_counts_BSMs
canon_col <- c(j="founder",a="a",d="d",e="e",v="vicariance",y="sympatry",s="subset")
summary <- data.frame(event_code=cols,event_type=unname(labels),
 mean_count=as.numeric(canon["means",canon_col[cols]]),
 sd_count=as.numeric(canon["stdevs",canon_col[cols]]),stringsAsFactors=FALSE)
summary$percent_all_events <- 100*summary$mean_count/sum(summary$mean_count)
summary$event_class <- ifelse(summary$event_code %in% c("a","d","e"),"anagenetic","cladogenetic")
summary$percent_within_class <- ave(summary$mean_count,summary$event_class,
 FUN=function(z)100*z/sum(z))
write.table(summary,file.path(out,"bsm_event_mean_sd_and_ratios.tsv"),sep="\t",quote=FALSE,row.names=FALSE)

# Compact publication plot analogous to the snook BSM table.
cols_plot <- c("j"="#D55E00","a"="#56B4E9","d"="#0072B2","e"="#999999",
 "v"="#CC79A7","y"="#009E73","s"="#F0E442")
plot_ylim <- c(0,max(summary$percent_all_events,na.rm=TRUE)*1.18)
pdf(file.path(out,"bsm_event_ratios.pdf"),width=11,height=9,useDingbats=FALSE)
par(mar=c(17,5.8,3.2,2.2),xpd=FALSE); bp<-barplot(summary$percent_all_events,names.arg=summary$event_type,
 las=2,cex.names=.92,cex.axis=1.0,cex.lab=1.08,
 col=cols_plot[summary$event_code],border=NA,ylab="Percentage of all events",
 main=paste0("100 biogeographic stochastic maps: ",model),cex.main=1.15,ylim=plot_ylim)
text(bp,summary$percent_all_events,sprintf("%.1f%%",summary$percent_all_events),pos=3,cex=.9)
dev.off()
png(file.path(out,"bsm_event_ratios.png"),width=3300,height=2700,res=300)
par(mar=c(17,5.8,3.2,2.2),xpd=FALSE); bp<-barplot(summary$percent_all_events,names.arg=summary$event_type,
 las=2,cex.names=.92,cex.axis=1.0,cex.lab=1.08,
 col=cols_plot[summary$event_code],border=NA,ylab="Percentage of all events",
 main=paste0("100 biogeographic stochastic maps: ",model),cex.main=1.15,ylim=plot_ylim)
text(bp,summary$percent_all_events,sprintf("%.1f%%",summary$percent_all_events),pos=3,cex=.9)
dev.off()

writeLines(c(paste("Selected MCC model:",model),"BSM replicates: 100",
 "Reported values: mean, SD, percent of all events, and percent within anagenetic/cladogenetic class",
 "Directionality matrices: canonical count_ana_clado_events() TSV outputs"),
 file.path(out,"README.txt"))
