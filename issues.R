## issues

# Friday 4th October 2019 - when doing xcmsSet of small sample dataset

# Attaching package: 'BiocGenerics'
# 
# The following objects are masked from 'package:parallel':
#   
#   clusterApply, clusterApplyLB, clusterCall, clusterEvalQ, clusterExport, clusterMap, parApply, parCapply,
# parLapply, parLapplyLB, parRapply, parSapply, parSapplyLB
# 
# The following objects are masked from 'package:stats':
#   
#   IQR, mad, sd, var, xtabs
# 
# The following objects are masked from 'package:base':
#   
#   anyDuplicated, append, as.data.frame, basename, cbind, colMeans, colnames, colSums, dirname, do.call,
# duplicated, eval, evalq, Filter, Find, get, grep, grepl, intersect, is.unsorted, lapply, lengths, Map,
# mapply, match, mget, order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank, rbind, Reduce,
# rowMeans, rownames, rowSums, sapply, setdiff, sort, table, tapply, union, unique, unsplit, which,
# which.max, which.min
# 
# Welcome to Bioconductor
# 
# Vignettes contain introductory material; view with 'browseVignettes()'. To cite Bioconductor, see
# 'citation("Biobase")', and for packages 'citation("pkgname")'.
# 
# Loading required package: BiocParallel
# Loading required package: MSnbase
# Loading required package: mzR
# Loading required package: Rcpp
# Loading required package: S4Vectors
# Loading required package: stats4
# 
# Attaching package: 'S4Vectors'
# 
# The following object is masked from 'package:base':
#   
#   expand.grid
# 
# Loading required package: ProtGenerics
# 
# This is MSnbase version 2.8.3 
# Visit https://lgatto.github.io/MSnbase/ to get started.
# 
# 
# Attaching package: 'MSnbase'
# 
# The following object is masked from 'package:stats':
#   
#   smooth
# 
# The following object is masked from 'package:base':
#   
#   trimws
# 
# 
# This is xcms version 3.4.4 
# 
# 
# Attaching package: 'xcms'
# 
# The following object is masked from 'package:stats':
#   
#   sigma
# 
# Detecting mass traces at 10 ppm ... OK
# Detecting chromatographic peaks in 71341 regions of interest ... OK: 7294 found.
# Warning messages:
#   1: In class(object) <- "environment" :
#   Setting class(x) to "environment" sets attribute to NULL; result will no longer be an S4 object
# 2: In serialize(data, node$con) :
#   'package:stats' may not be available when loading
# 3: In serialize(data, node$con) :
#   'package:stats' may not be available when loading
# 4: In serialize(data, node$con) :
#   'package:stats' may not be available when loading