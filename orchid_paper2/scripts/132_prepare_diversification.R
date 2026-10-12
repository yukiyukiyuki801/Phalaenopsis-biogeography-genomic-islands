suppressPackageStartupMessages(library(ape))
root <- normalizePath("orchid_paper2")
base <- file.path(root,"22_diversification")
.libPaths(c(file.path(base,"Rlib"), .libPaths()))
dir.create(file.path(base,"inputs"), recursive=TRUE, showWarnings=FALSE)
load(file.path(base,"inputs","InfTemp.rda"))
raw <- aggregate(Temperature~Age, InfTemp, mean)
raw <- raw[order(raw$Age),]
write.table(raw, file.path(base,"inputs","InfTemp_unique_ages.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
# Common fixed 0-7 Ma interval covers all current posterior crowns.
stopifnot(min(raw$Age)==0, max(raw$Age)>7)
grid <- seq(0,7,by=.001)
temp <- approx(raw$Age,raw$Temperature,xout=grid)$y
center <- mean(temp); spread <- sd(temp)
env <- data.frame(Age=grid, E=(temp-center)/spread, Temperature=temp)
write.table(env,file.path(base,"inputs","InfTemp_0_7Ma_standardized.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
write.table(data.frame(source="RPANDA InfTemp", raw_rows=nrow(InfTemp),unique_ages=nrow(raw),min_Ma=min(raw$Age),max_Ma=max(raw$Age),reference_max_Ma=7,temperature_mean=center,temperature_sd=spread,spline_df=30),file.path(base,"inputs","environment_metadata.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
registry <- list(); bamm_jobs <- list()
for (cal in c("zhang","givnish")) {
  indir <- file.path(root,"16_rad_beast","final_unconstrained_60m_summary")
  mcc <- read.tree(file.path(indir,paste0(cal,"_unconstrained_final_mcc_29tip.newick")))
  post <- read.tree(file.path(indir,paste0(cal,"_unconstrained_100_posterior_29tip.newick")))
  trees <- c(list(mcc),as.list(post))
  for (i in seq_along(trees)) {
    tr <- trees[[i]]; depths <- node.depth.edgelength(tr)[seq_len(Ntip(tr))]
    delta <- max(depths)-depths
    stopifnot(Ntip(tr)==29,is.binary(tr),all(tr$edge.length>=0),max(delta)<5e-5)
    # Correct only terminal Newick rounding, at most 30 years on Myr scale.
    term <- which(tr$edge[,2]<=Ntip(tr))
    tr$edge.length[term] <- tr$edge.length[term]+delta[tr$edge[term,2]]
    tr$root.edge <- NULL
    stopifnot(is.ultrametric(tr,tol=1e-8))
    id <- if(i==1) "mcc" else sprintf("post%03d",i-1)
    path <- file.path(base,"inputs",paste0(cal,"_",id,".nwk"))
    write.tree(tr,file=path,digits=16)
    registry[[length(registry)+1]] <- data.frame(calibration=cal,tree_id=id,path=path,tips=Ntip(tr),crown_Ma=max(depths),max_terminal_rounding_adjustment_Myr=max(delta))
    if(i==1) {
      prior <- BAMMtools::setBAMMpriors(tr,outfile=NULL,suppressWarning=TRUE)
      pv <- if(is.data.frame(prior)) setNames(as.numeric(prior$value),prior$param) else prior
      template <- readLines(file.path(base,"software","bamm","examples","diversification","whales","divcontrol.txt"))
      for(shift in c(.5,1,2)) for(run in 1:4) {
        out <- file.path(base,"bamm",cal,paste0("shiftprior_",shift),paste0("run",run))
        dir.create(out,recursive=TRUE,showWarnings=FALSE)
        settings <- list(treefile=path,numberOfGenerations=10000000,mcmcWriteFreq=2000,eventDataWriteFreq=2000,printFreq=100000,expectedNumberOfShifts=shift,lambdaInitPrior=pv[["lambdaInitPrior"]],lambdaShiftPrior=pv[["lambdaShiftPrior"]],muInitPrior=pv[["muInitPrior"]],numberOfChains=4,globalSamplingFraction=1,overwrite=0,lambdaInit0=log(29/2)/max(depths),muInit0=.1*log(29/2)/max(depths))
        text <- template
        for(k in names(settings)) {
          value <- if(is.numeric(settings[[k]])) format(settings[[k]],scientific=FALSE,trim=TRUE,digits=16) else settings[[k]]
          text <- sub(paste0("^",k,"[[:space:]]*=.*"),paste(k,"=",value),text)
        }
        seed <- 132000+1000*match(cal,c("zhang","givnish"))+100*shift+run
        text <- c(text,paste("seed =",as.integer(seed)))
        writeLines(text,file.path(out,"control.txt"))
        bamm_jobs[[length(bamm_jobs)+1]] <- data.frame(calibration=cal,shift_prior=shift,run=run,directory=out,seed=seed)
      }
    }
  }
}
write.table(do.call(rbind,registry),file.path(base,"inputs","tree_registry.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
write.table(do.call(rbind,bamm_jobs),file.path(base,"inputs","bamm_run_registry.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
capture.output(sessionInfo(),file=file.path(base,"inputs","preparation_sessionInfo.txt"))
