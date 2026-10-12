.libPaths(c(normalizePath("orchid_paper2/22_diversification/Rlib"),.libPaths()))
suppressPackageStartupMessages({library(ape);library(BAMMtools);library(data.table)})
base <- normalizePath("orchid_paper2/22_diversification")
out <- file.path(base,"summary")
jobs <- fread(file.path(base,"inputs","bamm_run_registry.tsv"))[shift_prior==1]
objects <- list(); audit <- list()
for(cal in c("zhang","givnish")) {
  tr <- read.tree(file.path(base,"inputs",paste0(cal,"_mcc.nwk")))
  pieces <- list()
  for(runid in 1:4) {
    j <- jobs[calibration==cal & run==runid]
    stopifnot(nrow(j)==1)
    ev <- fread(file.path(j$directory,"event_data.txt"))
    trace <- fread(file.path(j$directory,"mcmc_out.txt"))
    retained <- trace$generation[(floor(nrow(trace)*.25)+1):nrow(trace)]
    selected <- retained[round(seq(1,length(retained),length.out=500))]
    stopifnot(length(unique(selected))==500,all(selected %in% ev$generation))
    audit[[length(audit)+1]] <- data.table(calibration=cal,run=runid,
      original_generation=selected,pooled_sample=(runid-1)*500+seq_len(500))
    z <- ev[generation %in% selected]
    z[,generation:=(runid-1)*500+match(generation,selected)]
    pieces[[runid]] <- z
  }
  pooled <- rbindlist(pieces)
  ed <- getEventData(tr,as.data.frame(pooled),burnin=0,nsamples=NULL,verbose=FALSE)
  stopifnot(length(ed$eventData)==2000)
  objects[[cal]] <- dtRates(ed,tau=.005)
}
fwrite(rbindlist(audit),file.path(out,"bamm_phylogeny_selected_samples.tsv"),sep="\t")
saveRDS(objects,file.path(out,"bamm_primary_pooled_phylorate_objects.rds"))
ratevector <- function(ed,k) {
  if(k=="s") ed$dtrates$rates[[1]] else if(k=="e") ed$dtrates$rates[[2]] else ed$dtrates$rates[[1]]-ed$dtrates$rates[[2]]
}
kinds <- c("s","e","netdiv")
titles <- c(s="Branching rate",e="Loss parameter",netdiv="Net branching rate")
breaks <- lapply(kinds,function(k) {
  r <- range(unlist(lapply(objects,ratevector,k=k)),finite=TRUE)
  stopifnot(diff(r)>0)
  seq(r[1],r[2],length.out=65)
});names(breaks)<-kinds
fwrite(rbindlist(lapply(kinds,function(k)data.table(rate=k,minimum=min(breaks[[k]]),maximum=max(breaks[[k]])))),file.path(out,"bamm_phylogeny_color_scales.tsv"),sep="\t")
drawpanel <- function(cal,k) {
  ed <- objects[[cal]]; height <- max(ed$end)
  p <- plot(ed,spex=k,method="phylogram",labels=TRUE,legend=FALSE,
    lwd=2.4,cex=.85,pal=c("#2166AC","#67A9CF","#D0D0D0","#EF8A62","#B2182B"),
    colorbreaks=breaks[[k]],logcolor=FALSE,breaksmethod="linear",
    xlim=c(-.025*height,1.17*height),ylim=c(-9,29),mar=c(.5,.4,2,.4),par.reset=FALSE)
  title(main=paste(tools::toTitleCase(cal),titles[[k]],sep=" | "),cex.main=1.05)
  age <- pretty(c(0,height),n=4);age<-age[age>=0 & age<=height]
  axis(1,at=height-age,labels=format(age,trim=TRUE),pos=-.9,tcl=-.2,cex.axis=.85)
  text(height/2,-3.25,"Age (Ma)",cex=.85)
  addBAMMlegend(p,location=c(.08*height,.88*height,-5.9,-5.3),
    direction="horizontal",side=1,nTicks=2,cex.axis=.85,labelDist=.5)
  text(height*.48,-8.1,"Rate per lineage per Myr",cex=.8)
}
render <- function(name,cals,kinds,w,h) {
  draw <- function() {
    par(mfrow=c(length(cals),length(kinds)),oma=c(.2,.2,2.4,.2),family="sans")
    for(cal in cals) for(k in kinds) drawpanel(cal,k)
    mtext("BAMM: posterior mean rates on population-lineage trees",outer=TRUE,line=1,cex=1.15,font=2)
    mtext("Expected-shift prior = 1; four runs, 500 post-burn-in samples per run",outer=TRUE,line=-.1,cex=.85)
  }
  pdf(file.path(out,paste0(name,".pdf")),width=w,height=h,useDingbats=FALSE);draw();dev.off()
  png(file.path(out,paste0(name,".png")),width=w,height=h,units="in",res=300,type="cairo");draw();dev.off()
}
render("bamm_primary_phylogeny_rates",c("zhang","givnish"),kinds,15,12)
for(cal in names(objects)) render(paste0("bamm_",cal,"_phylogeny_rates"),cal,kinds,15,6.5)
render("bamm_primary_branching_phylogenies",c("zhang","givnish"),"s",6,12)
cat("Completed four PDF/PNG figure pairs; 2,000 samples per calibration.\n")
