          seed = SEED_VALUE
       seqfile = ../../../matrix/strict_single_copy_concat_aa.phy
      treefile = ../../../config/TREE_FILE
       outfile = out.txt
      mcmcfile = mcmc.txt

         ndata = 1
       seqtype = 2
       usedata = 2 in.BV
         clock = 2
       RootAge = '<2.0'

         model = 3
    aaRatefile = /vast/palmer/apps/avx2/software/PAML/4.10.7-GCCcore-12.2.0/dat/wag.dat
         alpha = 0.5
         ncatG = 5
     cleandata = 0

       BDparas = 1 1 0.1
   rgene_gamma = 2 20 1
  sigma2_gamma = 1 10 1
      finetune = 1: .1 .1 .1 .1 .1 .1
         print = 1
        burnin = 100000
      sampfreq = 20
       nsample = 100000
