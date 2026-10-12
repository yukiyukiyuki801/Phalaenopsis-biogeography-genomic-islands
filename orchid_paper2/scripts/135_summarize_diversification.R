.libPaths(c(normalizePath("orchid_paper2/22_diversification/Rlib"),.libPaths()))
suppressPackageStartupMessages({library(ape);library(coda);library(BAMMtools);library(ggplot2);library(data.table)})
base <- normalizePath("orchid_paper2/22_diversification")
out <- file.path(base,"summary")
dir.create(out,showWarnings=FALSE)
saveplot <- function(p,name,w=10,h=7) {
  ggsave(file.path(out,paste0(name,".pdf")),p,width=w,height=h,bg="white")
  ggsave(file.path(out,paste0(name,".png")),p,width=w,height=h,dpi=300,bg="white")
}
mode <- commandArgs(trailingOnly=TRUE)[1]
if(mode=="rpanda") {
  files <- list.files(file.path(base,"rpanda"),pattern="_summary.tsv$",recursive=TRUE,full.names=TRUE)
  stopifnot(length(files)>0)
  d <- rbindlist(lapply(files,fread))
  d[,delta_AICc:=AICc-min(AICc),by=.(calibration,tree_id)]
  d[,Akaike_weight:=exp(-delta_AICc/2)/sum(exp(-delta_AICc/2)),by=.(calibration,tree_id)]
  d[,n_fitted_models:=.N,by=.(calibration,tree_id)]
  d[,complete_candidate_set:=n_fitted_models==9]
  d[,rank:=frank(AICc,ties.method="min"),by=.(calibration,tree_id)]
  setorder(d,calibration,tree_id,rank)
  fwrite(d,file.path(out,"rpanda_model_comparison.tsv"),sep="\t")
  # Missing model fits are conspicuous and not treated as zero support.
  if(any(!d$complete_candidate_set)) warning("Incomplete candidate sets; weights are provisional.")
  p <- ggplot(d[tree_id=="mcc"],aes(model,Akaike_weight,fill=family))+
    geom_col()+facet_wrap(~calibration,ncol=1)+
    scale_fill_manual(values=c(N="grey55",E="#0072B2",T="#E69F00"),labels=c(E="Paleotemperature",N="Constant null",T="Time-only"))+
    labs(x="Candidate model",y="Akaike weight",fill=NULL,title="RPANDA: branching among sampled population lineages",subtitle="Nine-model comparison; N1 = constant birth-death null")+theme_bw(base_size=12)
  saveplot(p,"rpanda_mcc_model_weights")
  winners <- d[tree_id=="mcc",.SD[which.min(AICc)],by=calibration]
  curves <- rbindlist(lapply(seq_len(nrow(winners)),function(i){
    w<-winners[i]; z<-fread(file.path(base,"rpanda",w$calibration,w$tree_id,paste0(w$model,"_rate_curve.tsv")))
    z[,`:=`(calibration=w$calibration,model=w$model)];z
  }))
  fwrite(curves,file.path(out,"rpanda_best_mcc_rate_curves.tsv"),sep="\t")
  long<-melt(curves,id.vars=c("age_Ma","calibration","model"),variable.name="rate",value.name="value")
  p<-ggplot(long,aes(age_Ma,value,color=rate))+geom_line()+facet_wrap(~calibration,scales="free_y",ncol=1)+scale_x_reverse()+labs(x="Age (Ma)",y="Model rate per lineage per Myr",title="Best MCC models: conditional sampled-lineage rate estimates")+theme_bw(base_size=12)
  saveplot(p,"rpanda_best_mcc_rate_curves")
} else if(mode=="bamm") {
  jobs<-fread(file.path(base,"inputs","bamm_run_registry.tsv"))
  ess<-list();shifts<-list();rates<-list();chains<-list()
  for(i in seq_len(nrow(jobs))) {
    j<-jobs[i]; path<-file.path(j$directory,"mcmc_out.txt")
    if(!file.exists(path)) next
    z<-fread(path)
    z<-z[(floor(nrow(z)*.25)+1):nrow(z)]
    numeric_cols<-names(z)[vapply(z,is.numeric,logical(1))]
    numeric_cols<-setdiff(numeric_cols,c("generation"))
    for(k in numeric_cols) {
      v<-z[[k]]
      e<-if(sd(v)==0) NA_real_ else as.numeric(effectiveSize(v))
      ess[[length(ess)+1]]<-data.table(calibration=j$calibration,shift_prior=j$shift_prior,run=j$run,parameter=k,ESS=e,constant_trace=sd(v)==0,n_samples=length(v))
    }
    if("N_shifts" %in% names(z)) {
      s<-z[,.(n=.N),by=N_shifts];s[,`:=`(probability=n/sum(n),calibration=j$calibration,shift_prior=j$shift_prior,run=j$run)]
      shifts[[length(shifts)+1]]<-s
    }
    chains[[paste(j$calibration,j$shift_prior,j$run,sep="_")]]<-z
    bf<-try(computeBayesFactors(path,expectedNumberOfShifts=j$shift_prior,burnin=.25),silent=TRUE)
    saveRDS(bf,file.path(j$directory,"shift_count_Bayes_factors.rds"))
    tr<-read.tree(file.path(base,"inputs",paste0(j$calibration,"_mcc.nwk")))
    ed<-getEventData(tr,eventdata=file.path(j$directory,"event_data.txt"),burnin=.25,nsamples=500,verbose=FALSE)
    rm<-getRateThroughTimeMatrix(ed,nslices=201)
    for(k in c("lambda","mu","net")) {
      mat<-if(k=="net") rm$lambda-rm$mu else rm[[k]]
      qs<-apply(mat,2,quantile,probs=c(.025,.5,.975),na.rm=TRUE)
      rates[[length(rates)+1]]<-data.table(calibration=j$calibration,shift_prior=j$shift_prior,run=j$run,rate=k,age_Ma=max(branching.times(tr))-rm$times,lower=qs[1,],median=qs[2,],upper=qs[3,])
    }
    saveRDS(getBestShiftConfiguration(ed,expectedNumberOfShifts=j$shift_prior),file.path(j$directory,"best_shift_configuration.rds"))
  }
  stopifnot(length(ess)>0)
  fwrite(rbindlist(ess),file.path(out,"bamm_all_parameter_ESS.tsv"),sep="\t")
  saveRDS(chains,file.path(out,"bamm_postburn_scalar_traces.rds"))
  psrf<-list()
  for(cal in c("zhang","givnish")) for(pr in c(.5,1,2)) {
    keys<-paste(cal,pr,1:4,sep="_")
    zz<-chains[intersect(keys,names(chains))]
    if(length(zz)!=4) next
    common<-Reduce(intersect,lapply(zz,names))
    common<-setdiff(common,c("generation"))
    common<-common[vapply(common,function(k) all(vapply(zz,function(z) is.numeric(z[[k]])&&is.finite(sd(z[[k]]))&&sd(z[[k]])>0,logical(1))),logical(1))]
    n<-min(vapply(zz,nrow,integer(1)))
    if(!length(common)) next
    mm<-mcmc.list(lapply(zz,function(z) mcmc(as.matrix(z[seq_len(n),..common]))))
    gd<-try(gelman.diag(mm,autoburnin=FALSE,multivariate=FALSE)$psrf,silent=TRUE)
    if(!inherits(gd,"try-error")) psrf[[length(psrf)+1]]<-data.table(calibration=cal,shift_prior=pr,parameter=rownames(gd),PSRF=gd[,1],PSRF_upper=gd[,2])
  }
  if(length(psrf)) fwrite(rbindlist(psrf),file.path(out,"bamm_between_run_PSRF.tsv"),sep="\t")
  if(length(shifts)) {
    s<-rbindlist(shifts);fwrite(s,file.path(out,"bamm_shift_count_probabilities_by_run.tsv"),sep="\t")
    p<-ggplot(s,aes(N_shifts,probability,color=factor(run)))+geom_point()+geom_line()+facet_grid(calibration~shift_prior)+labs(x="Number of rate shifts",y="Posterior probability",color="Independent run")+theme_bw(base_size=12)
    saveplot(p,"bamm_shift_count_support")
  }
  r<-rbindlist(rates);fwrite(r,file.path(out,"bamm_rate_curves_by_run.tsv"),sep="\t")
  p<-ggplot(r[shift_prior==1],aes(age_Ma,median,color=factor(run)))+geom_line()+facet_grid(calibration~rate,scales="free_y")+scale_x_reverse()+labs(x="Age (Ma)",y="Model rate per lineage per Myr",color="Independent run",title="BAMM: exploratory rates among sampled population lineages")+theme_bw(base_size=12)
  saveplot(p,"bamm_primary_rate_curves_by_run",12,7)
}
