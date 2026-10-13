#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(BioGeoBEARS); library(ape)})
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 2L || length(args) > 3L) stop("usage: tree_id model [max_range_size]")
tree_id <- args[1]; model <- args[2]
max_range_size <- if(length(args)==3L) as.integer(args[3]) else 3L
if(!max_range_size %in% 2:4) stop("max_range_size must be 2, 3, or 4")
models <- c("DEC", "DECJ", "DIVALIKE", "DIVALIKEJ", "BAYAREALIKE", "BAYAREALIKEJ")
stopifnot(model %in% models)
root <- normalizePath(Sys.getenv("ORCHID_PAPER2_ROOT"), mustWork=TRUE)
base_env <- Sys.getenv("BGB_BASE",unset="")
if(!nzchar(base_env)) stop("Set BGB_BASE to a calibration-specific unconstrained_29tip/<calibration>/max_range_3 directory")
base <- normalizePath(base_env,mustWork=TRUE)
inp <- file.path(base, "inputs")
tree_file <- if (tree_id == "mcc") file.path(inp, "mcc_myr.newick") else
  file.path(inp, sprintf("posterior_%03d_myr.newick", as.integer(tree_id)))
out <- file.path(base, "runs", tree_id)
dir.create(out, recursive=TRUE, showWarnings=FALSE)

run <- define_BioGeoBEARS_run()
run$trfn <- tree_file
run$geogfn <- file.path(inp, "geography.data")
run$max_range_size <- max_range_size
run$min_branchlength <- 1e-6
run$include_null_range <- TRUE
run$on_NaN_error <- -1e50
run$speedup <- TRUE
run$use_optimx <- TRUE
run$num_cores_to_use <- 1
run$force_sparse <- FALSE
run <- readfiles_BioGeoBEARS_run(run)
run$return_condlikes_table <- TRUE
run$calc_TTL_loglike_from_condlikes_table <- TRUE
run$calc_ancprobs <- TRUE
# An explicit state list permits any range whose geographic graph is
# connected, including chains in which the terminal areas are not directly
# adjacent. Codes are converted to BioGeoBEARS' zero-based state indexes.
allowed_fn <- file.path(inp,"allowed_states.tsv")
if(file.exists(allowed_fn)) {
  allowed <- read.delim(allowed_fn,stringsAsFactors=FALSE)$state_code
  area <- read.delim(file.path(inp,"area_order.tsv"),stringsAsFactors=FALSE)
  parse_state <- function(code) {
    if(code=="_") return(NA_integer_)
    idx <- match(strsplit(code,"+",fixed=TRUE)[[1]],area$code)
    if(anyNA(idx)) stop("Unknown area code in allowed state: ",code)
    as.integer(idx-1L)
  }
  run$states_list <- lapply(allowed,parse_state)
}
# Native BioGeoBEARS adjacency requires every pair of areas in a range to be
# directly adjacent. Retain it only for legacy inputs without an explicit
# graph-connected state list.
adj_fn <- file.path(inp,"areas_adjacency.tsv")
if(file.exists(adj_fn) && !file.exists(allowed_fn)) {
  adj <- as.matrix(read.delim(adj_fn,row.names=1,check.names=FALSE))
  storage.mode(adj) <- "numeric"
  run$list_of_areas_adjacency_mats <- list(adj)
}
p <- run$BioGeoBEARS_model_object@params_table
if (grepl("DIVALIKE", model)) {
  p["s", c("type","init","est")] <- list("fixed",0,0)
  p["ysv","type"] <- "2-j"; p["ys","type"] <- "ysv*1/2"
  p["y","type"] <- "ysv*1/2"; p["v","type"] <- "ysv*1/2"
  p["mx01v", c("type","init","est")] <- list("fixed",0.5,0.5)
}
if (grepl("BAYAREALIKE", model)) {
  p["s", c("type","init","est")] <- list("fixed",0,0)
  p["v", c("type","init","est")] <- list("fixed",0,0)
  p["ysv","type"] <- "1-j"; p["ys","type"] <- "ysv*1/1"; p["y","type"] <- "1-j"
  p["mx01y", c("type","init","est")] <- list("fixed",0.9999,0.9999)
}
if (grepl("J$", model)) {
  p["j", c("type","init","est")] <- list("free",0.0001,0.0001)
  p["j","min"] <- 0.00001
  p["j","max"] <- if (grepl("DIVALIKE",model)) 1.99999 else 0.99999
}
run$BioGeoBEARS_model_object@params_table <- p
check_BioGeoBEARS_run(run)
res <- bears_optim_run(run)
save(res, run, file=file.path(out, paste0(model,".RData")))
pt <- res$outputs@params_table
tab <- data.frame(tree_id=tree_id, model=model, log_likelihood=res$total_loglikelihood,
 n_parameters=sum(pt$type=="free"), d=as.numeric(pt["d","est"]),
 e=as.numeric(pt["e","est"]), j=as.numeric(pt["j","est"]),
 n_tips=Ntip(read.tree(tree_file)))
write.table(tab, file.path(out,paste0(model,".summary.tsv")), sep="\t", quote=FALSE, row.names=FALSE)
