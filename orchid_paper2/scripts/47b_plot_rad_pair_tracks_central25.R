#!/usr/bin/env Rscript
# Reuse the established three-track plotting code with the stricter central-25%
# dXY classification and a separate output directory.
x=readLines("orchid_paper2/scripts/47_plot_rad_pair_tracks_models.R",warn=FALSE)
x=gsub('rad_locus_statistics.four_models_q25.tsv',
       'rad_locus_statistics.four_models_central25_dxy.tsv',x,fixed=TRUE)
x=gsub('plots","pair_tracks', 'plots","pair_tracks_central25_dxy',x,fixed=TRUE)
x=gsub('four_model_class_q25', 'four_model_class_central25_dxy',x,fixed=TRUE)
x=gsub('quartile-classified models', 'central-25% dXY model classes',x,fixed=TRUE)
x=gsub('model_bands.lines', 'central25_dxy.model_bands.lines',x,fixed=TRUE)
x=gsub('model_bands.pdf', 'central25_dxy.model_bands.pdf',x,fixed=TRUE)
x=gsub('model_bands.png', 'central25_dxy.model_bands.png',x,fixed=TRUE)
eval(parse(text=x),envir=globalenv())
