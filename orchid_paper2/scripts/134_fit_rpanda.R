.libPaths(c(normalizePath("orchid_paper2/22_diversification/Rlib"),.libPaths()))
suppressPackageStartupMessages({library(ape);library(RPANDA)})
args <- commandArgs(trailingOnly=TRUE)
task <- as.integer(args[1]); mode <- if(length(args)>1) args[2] else "mcc"
base <- normalizePath("orchid_paper2/22_diversification")
registry <- read.delim(file.path(base,"inputs","tree_registry.tsv"))
registry <- if(mode=="mcc") registry[registry$tree_id=="mcc",] else registry[registry$tree_id!="mcc",]
models <- c("N1","E1","E2","E3","E4","T1","T2","T3","T4")
row <- registry[(task-1)%/%length(models)+1,]; model <- models[(task-1)%%length(models)+1]
stopifnot(!is.na(row$path))
tr <- read.tree(row$path); age <- max(branching.times(tr))
env0 <- read.delim(file.path(base,"inputs","InfTemp_0_7Ma_standardized.tsv"))
env <- env0[,c("Age","E")]
stopifnot(age<=max(env$Age),is.ultrametric(tr))
out <- file.path(base,"rpanda",row$calibration,row$tree_id)
dir.create(out,recursive=TRUE,showWarnings=FALSE)
family <- substr(model,1,1)
vary_l <- model %in% c("E1","E2","E4","T1","T2","T4")
vary_m <- model %in% c("E3","E4","T3","T4")
# N1 now denotes the constant birth-death null (formerly N2).
zero_m <- model %in% c("E1","T1")
# Log-rate intercepts guarantee positive rates without reflected negatives.
lc <- function(t,p) exp(p[1])+0*t
lv <- function(t,p) exp(p[1]+p[2]*t)
mc <- function(t,p) exp(p[1])+0*t
mv <- function(t,p) exp(p[1]+p[2]*t)
mz <- function(t,p) 0*t
elc <- function(t,x,p) exp(p[1])+0*t
elv <- function(t,x,p) exp(p[1]+p[2]*x)
emc <- function(t,x,p) exp(p[1])+0*t
emv <- function(t,x,p) exp(p[1]+p[2]*x)
emz <- function(t,x,p) 0*t
lambda0 <- log(Ntip(tr)/2)/age
starts <- expand.grid(mult=c(.5,1.5),beta=c(-.5,0,.5))
fits <- list(); records <- list()
for(i in seq_len(nrow(starts))) {
  lp <- c(log(lambda0*starts$mult[i]),if(vary_l) starts$beta[i])
  mp <- if(zero_m) numeric(0) else c(log(lambda0*.2),if(vary_m) starts$beta[i])
  messages <- character(0)
  f <- tryCatch(withCallingHandlers({
    if(family=="E") fit_env(tr,env,age,if(vary_l) elv else elc,if(zero_m) emz else if(vary_m) emv else emc,
       lp,mp,df=30,f=1,cst.lamb=!vary_l,cst.mu=!vary_m,fix.mu=zero_m,dt=.001,cond="crown")
    else fit_bd(tr,age,if(vary_l) lv else lc,if(zero_m) mz else if(vary_m) mv else mc,
       lp,mp,f=1,cst.lamb=!vary_l,cst.mu=!vary_m,expo.lamb=vary_l,expo.mu=vary_m,fix.mu=zero_m,dt=.001,cond="crown")
  },warning=function(w){messages<<-c(messages,conditionMessage(w));invokeRestart("muffleWarning")}),error=function(e)e)
  ok <- !inherits(f,"error") && is.finite(f$LH) && is.finite(f$aicc)
  if(ok) fits[[as.character(i)]] <- f
  records[[i]] <- data.frame(start=i,success=ok,logLik=if(ok) f$LH else NA_real_,AICc=if(ok) f$aicc else NA_real_,message=if(inherits(f,"error")) conditionMessage(f) else paste(unique(messages),collapse=";"))
}
attempts <- do.call(rbind,records)
write.table(attempts,file.path(out,paste0(model,"_optimization_starts.tsv")),sep="\t",quote=TRUE,row.names=FALSE)
if(!length(fits)) stop("All optimization starts failed: ",model)
bestid <- names(which.max(vapply(fits,function(z) z$LH,numeric(1))))
best <- fits[[bestid]]
saveRDS(best,file.path(out,paste0(model,"_fit.rds")))
k <- length(best$lamb_par)+if(zero_m) 0 else length(best$mu_par)
summary <- data.frame(calibration=row$calibration,tree_id=row$tree_id,model=model,family=family,tips=Ntip(tr),crown_Ma=age,k=k,logLik=best$LH,AICc=best$aicc,beta_lambda=if(vary_l) best$lamb_par[2] else NA_real_,beta_mu=if(vary_m) best$mu_par[2] else NA_real_,successful_starts=sum(attempts$success),starts_within_0.01_logLik=sum(abs(attempts$logLik-best$LH)<.01,na.rm=TRUE),sampling_fraction=1,spline_df=if(family=="E") 30 else NA_integer_,dt_Myr=.001)
write.table(summary,file.path(out,paste0(model,"_summary.tsv")),sep="\t",quote=FALSE,row.names=FALSE)
times <- seq(0,age,length.out=501)
if(family=="E") {
  spl <- pspline::sm.spline(env$Age,env$E,df=30)
  values <- predict(spl,times)
  lamb <- (if(vary_l) elv else elc)(times,values,best$lamb_par)
  mu <- if(zero_m) rep(0,length(times)) else (if(vary_m) emv else emc)(times,values,best$mu_par)
} else {
  lamb <- (if(vary_l) lv else lc)(times,best$lamb_par)
  mu <- if(zero_m) rep(0,length(times)) else (if(vary_m) mv else mc)(times,best$mu_par)
}
write.table(data.frame(age_Ma=times,lambda=lamb,mu=mu,net=lamb-mu),file.path(out,paste0(model,"_rate_curve.tsv")),sep="\t",quote=FALSE,row.names=FALSE)
print(summary)
