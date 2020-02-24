## Resort CDF files
## Written on: 1 March 2011
## Written By: Paul Benton
## Updated on: 27 Oct 2011
## Notes: This program will check the m/z vector for a single file for sort violation
##	If a sort violation is found the m/z vector is reorganised along with the Intensity 
##	vector. Finally a new CDF file is made which is fixed.
## !!!NB!!! : Parallel version does not report progress

require(xcms)
library(caTools)

checkAllcdfs<-function(nSlaves=1){
  AllCDFs<-list.files(recursive=TRUE, pattern="cdf", ignore.case=TRUE, full.names=TRUE)
  if(nSlaves >1){
    rmpi = "Rmpi"
    if(require(rmpi, character.only=TRUE,quietly=TRUE)){
      cl <- makeCluster(nSlaves, type = "MPI") 
    } else if(require(snow)){
      cl <- makeCluster(nSlaves, type = "SOCK")
    }
    clusterEvalQ(cl, library(xcms))
    unlist(clusterApply(cl, AllCDFs, checkCDFfile))
    stopCluster(cl)
  } else{
    sapply(AllCDFs, checkCDFfile)
    cat("\n")
  }
}

checkCDFfile<-function(file, type=".cdf"){
  cat("\n")
  cat(paste("Loading File:", file, sep=""))
  xr<-xcmsRaw(file, profstep=0)
  for(i in 1:length(xr@scanindex)){
    scan<-getScan(xr, scan=i)
    if(is.unsorted(scan[,"mz"]) == TRUE){
      cat(" x ")
      newfile<-sub(type, ".mzdata", file, ignore.case=TRUE)
      write.mzdata(xr, newfile)
      file.copy(file, sub(type, ".OLD", file, ignore.case=TRUE))
      unlink(file)
      rm(list=ls())
      gc()
      return(1)
    }
    if(i == length(xr@scanindex)){
      cat(" O ")
      rm(list=ls())
      gc()
      return(0)
    }
  }
}


checkAllcdfs()