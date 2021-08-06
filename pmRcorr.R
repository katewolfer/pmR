myCorrelation   <- function(inputMat){
  ## return Pearson correlation and p-value
  # $corr for correlation
  # $pval for pvalue
  
  corr <- matrix(NA, nrow=ncol(inputMat), ncol=ncol(inputMat))
  pval <- matrix(NA, nrow=ncol(inputMat), ncol=ncol(inputMat))
  colnames(corr) <- colnames(inputMat)
  colnames(pval) <- colnames(inputMat)
  rownames(corr) <- colnames(inputMat)
  rownames(pval) <- colnames(inputMat)
  
  # peason correlation and correlation p-value
  for(icol in 1:ncol(inputMat)){
    for(jrow in 1:ncol(inputMat)){
      tmp             <- cor.test( inputMat[,jrow], inputMat[,icol], method="pearson")
      corr[jrow,icol] <- tmp$estimate
      pval[jrow,icol] <- tmp$p.value
    }
  }
  
  return(list(corr=corr, pval=pval))
}         