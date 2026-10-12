.libPaths(c(normalizePath("orchid_paper2/22_diversification/Rlib"),.libPaths()))
suppressPackageStartupMessages({library(ggplot2);library(data.table);library(pspline)})
base<-normalizePath("orchid_paper2/22_diversification");out<-file.path(base,"summary")
cmp<-fread(file.path(out,"rpanda_model_comparison.tsv"))[tree_id=="mcc"]
env<-fread(file.path(base,"inputs","InfTemp_0_7Ma_standardized.tsv"))
meta<-fread(file.path(base,"inputs","environment_metadata.tsv"))
spl<-sm.spline(env$Age,env$E,df=30)
ratecols<-c(lambda="#D55E00",mu="#0072B2",net="#7B3294")
calcols<-c(zhang="#D55E00",givnish="#0072B2")
savefig<-function(p,name,w=11,h=8) {
  ggsave(file.path(out,paste0(name,".pdf")),p,width=w,height=h,bg="white")
  ggsave(file.path(out,paste0(name,".png")),p,width=w,height=h,dpi=300,bg="white")
}
theme_set(theme_bw(base_size=12)+theme(panel.grid.minor=element_blank(),
  plot.caption=element_text(hjust=0,size=9),strip.text=element_text(size=10)))
response<-list();temporal<-list();climate<-list();diagnostics<-list()
for(cal in c("zhang","givnish")) {
  age<-unique(cmp[calibration==cal]$crown_Ma)
  times<-seq(0,age,length.out=1001)
  E<-as.numeric(predict(spl,times))
  temp<-E*meta$temperature_sd+meta$temperature_mean
  climate[[cal]]<-data.table(calibration=cal,age_Ma=times,temperature_C=temp)
  for(model in c("N1",paste0("E",1:4),paste0("T",1:4))) {
    f<-readRDS(file.path(base,"rpanda",cal,"mcc",paste0(model,"_fit.rds")))
    fam<-substr(model,1,1)
    vl<-model %in% c("E1","E2","E4","T1","T2","T4")
    vm<-model %in% c("E3","E4","T3","T4")
    zm<-model %in% c("E1","T1")
    rates<-function(x) {
      l<-exp(f$lamb_par[1]+if(vl) f$lamb_par[2]*x else 0)+0*x
      m<-if(zm) 0*x else exp(f$mu_par[1]+if(vm) f$mu_par[2]*x else 0)+0*x
      data.table(lambda=l,mu=m,net=l-m)
    }
    z<-cbind(data.table(calibration=cal,model=model,age_Ma=times,temperature_C=temp),rates(if(fam=="E") E else times))
    temporal[[length(temporal)+1]]<-z
    if(fam=="E") {
      gridE<-seq(min(E),max(E),length.out=501)
      r<-cbind(data.table(calibration=cal,model=model,temperature_C=gridE*meta$temperature_sd+meta$temperature_mean),rates(gridE))
      response[[length(response)+1]]<-r
    }
    model_id<-model
    s<-cmp[calibration==cal & model==model_id]
    stopifnot(nrow(s)==1)
    diagnostics[[length(diagnostics)+1]]<-data.table(calibration=cal,model=model,
      beta_lambda=s$beta_lambda,beta_mu=s$beta_mu,delta_AICc=s$delta_AICc,
      successful_starts=s$successful_starts,near_best_starts=s$starts_within_0.01_logLik,
      lambda_min=min(z$lambda),lambda_max=max(z$lambda),mu_min=min(z$mu),mu_max=max(z$mu),
      net_min=min(z$net),net_max=max(z$net),
      provisional=s$starts_within_0.01_logLik<3 || max(z$mu)>100)
  }
}
temporal<-rbindlist(temporal);response<-rbindlist(response);climate<-rbindlist(climate);diagnostics<-rbindlist(diagnostics)
fwrite(temporal,file.path(out,"rpanda_all_models_temporal_predictions.tsv"),sep="\t")
fwrite(response,file.path(out,"rpanda_environment_response_predictions.tsv"),sep="\t")
fwrite(climate,file.path(out,"rpanda_smoothed_temperature_for_plots.tsv"),sep="\t")
fwrite(diagnostics,file.path(out,"rpanda_response_diagnostics.tsv"),sep="\t")
labels<-c(lambda="Branching (lambda)",mu="Loss (mu)",net="Net (lambda - mu)")
makepanels<-function(d,name,xlab,title,reverse=FALSE) {
  long<-melt(d,id.vars=c("calibration","model","age_Ma","temperature_C")[c("calibration","model","age_Ma","temperature_C") %in% names(d)],
    measure.vars=c("lambda","mu","net"),variable.name="rate",value.name="value")
  flags<-diagnostics[calibration==unique(d$calibration)]
  long[,model_label:=paste0(model,ifelse(model %in% flags[provisional==TRUE]$model," *",""))]
  long[,rate:=factor(rate,levels=c("lambda","mu","net"),labels=labels)]
  x<-if(reverse) "age_Ma" else "temperature_C"
  p<-ggplot(long,aes(.data[[x]],value,color=rate))+geom_line(linewidth=.75)+
    facet_wrap(~model_label+rate,ncol=3,scales="free_y")+
    scale_color_manual(values=setNames(ratecols,labels),guide="none")+
    labs(x=xlab,y="Fitted rate per lineage per Myr",title=title,
      subtitle="Model predictions; each panel has its own y scale",
      caption="* Provisional optimization/extreme-rate fit: do not interpret literally.\nNo uncertainty intervals available; curve shape is imposed by the fitted model.")
  if(reverse) p<-p+scale_x_reverse()
  savefig(p,name,11,10)
}
for(cal in c("zhang","givnish")) {
  makepanels(response[calibration==cal],paste0("rpanda_",cal,"_environment_rate_responses"),
    "Isotope-derived paleotemperature (degrees C)",paste(tools::toTitleCase(cal),"| E1-E4: fitted temperature responses"))
  makepanels(temporal[calibration==cal & substr(model,1,1)=="E"],paste0("rpanda_",cal,"_environment_rates_through_time"),
    "Age (Ma)",paste(tools::toTitleCase(cal),"| E1-E4: fitted rates through time"),TRUE)
  makepanels(temporal[calibration==cal & substr(model,1,1)=="T"],paste0("rpanda_",cal,"_time_rate_responses"),
    "Age (Ma)",paste(tools::toTitleCase(cal),"| T1-T4: fitted temporal responses"),TRUE)
}
# Figure-5 analogue: show E1 for both calibrations, not a selected extreme mu fit.
main<-rbindlist(list(
  climate[,.(calibration,panel="A. Paleotemperature through time",x=-age_Ma,y=temperature_C)],
  response[model=="E1",.(calibration,panel="B. E1: net rate versus temperature",x=temperature_C,y=net)],
  temporal[model=="E1",.(calibration,panel="C. E1: net rate through time",x=-age_Ma,y=net)]))
main[,calibration:=factor(calibration,levels=c("zhang","givnish"),labels=c("Zhang","Givnish"))]
p<-ggplot(main,aes(x,y,color=calibration))+geom_line(linewidth=.9)+
  facet_wrap(~calibration+panel,ncol=3,scales="free")+
  scale_color_manual(values=setNames(calcols,c("Zhang","Givnish")),guide="none")+
  scale_x_continuous(labels=function(v) format(abs(v),trim=TRUE))+
  labs(x="Panels A/C: age (Ma), past to present; panel B: paleotemperature (degrees C)",
    y="Panel A: temperature (degrees C); panels B/C: net rate per lineage per Myr",
    title="RPANDA: paleotemperature and fitted population-lineage branching",
    subtitle="E1 allows temperature-dependent branching and fixes loss to zero: net rate = branching rate",
    caption="E1 is competitive, not uniquely supported. Its temperature slope changes sign between calibrations.\nThese are fitted responses, not an independent correlation test; no confidence bands are implied.")
savefig(p,"rpanda_temperature_net_rate_main",13,8)
cat("Saved seven PDF/PNG figure pairs and prediction/diagnostic tables.\n")
