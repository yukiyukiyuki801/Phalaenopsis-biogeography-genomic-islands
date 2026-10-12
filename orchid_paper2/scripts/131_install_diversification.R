root <- normalizePath("orchid_paper2")
lib <- file.path(root, "22_diversification", "Rlib")
dir.create(lib, recursive=TRUE, showWarnings=FALSE)
.libPaths(c(lib, .libPaths()))
options(repos=c(CRAN="https://cloud.r-project.org"), Ncpus=8)
Sys.setenv(MAKEFLAGS="-j8")
needed <- c("ape", "coda", "BAMMtools", "RPANDA", "ggplot2", "data.table")
missing <- setdiff(needed[!vapply(needed, requireNamespace, logical(1), quietly=TRUE)],"RPANDA")
if (length(missing)) install.packages(missing, lib=lib)
if (!requireNamespace("RPANDA",quietly=TRUE)) {
  deps <- c("ape","picante","BB","bipartite","cluster","coda","deSolve","fields","fpc","geiger","glassoFast","igraph","Matrix","mvtnorm","mvMORPH","GUniFrac","ParallelLogger","phangorn","phytools","pracma","pspline","pvclust","raster","RColorBrewer","R.utils","Rmpfr","TESS","vegan")
  missing_deps <- deps[!vapply(deps,requireNamespace,logical(1),quietly=TRUE)]
  missing_cran <- setdiff(missing_deps,"TESS")
  if(length(missing_cran)) install.packages(missing_cran,lib=lib)
  if(!requireNamespace("TESS",quietly=TRUE)) {
    tess_src <- file.path(root,"22_diversification","software","TESS_cran_mirror.tar.gz")
    download.file("https://codeload.github.com/cran/TESS/tar.gz/refs/heads/master",tess_src,mode="wb")
    install.packages(tess_src,repos=NULL,type="source",lib=lib)
  }
  stopifnot(all(vapply(deps,requireNamespace,logical(1),quietly=TRUE)))
  src <- file.path(root,"22_diversification","software","RPANDA_cran_mirror.tar.gz")
  download.file("https://codeload.github.com/cran/RPANDA/tar.gz/refs/heads/master",src,mode="wb")
  install.packages(src, repos=NULL, type="source",lib=lib)
}
stopifnot(all(vapply(needed, requireNamespace, logical(1), quietly=TRUE)))
versions <- data.frame(package=needed, version=vapply(needed, function(p) as.character(packageVersion(p)), character(1)))
write.table(versions, file.path(root,"22_diversification","software_versions.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
capture.output(sessionInfo(), file=file.path(root,"22_diversification","R_sessionInfo.txt"))
