whittleMetadata <- function(metadata, metadataSampleList, subsetSamples) {
  
  findOverlap <- is.element(metadataSampleList, subsetSamples)
  lookAtOverlap <- metadataSampleList[which(findOverlap == TRUE)]

  # initialise subsetted dataframe
  findSample <- grep(subsetSamples[1], metadataSampleList)
  siphonMetadata <- metadata[findSample,]
  
  # collect other rows
  for (i in 1:length(subsetSamples)){
    checkXL <- grep("XL", subsetSamples[i])
    if (length(checkXL) >= 1) {
      snipGLA <- substr(subsetSamples[i], 1, 10)
      findSample <- grep(snipGLA, metadataSampleList)
      siphonMetadata[i,] <- metadata[findSample,]
    } else {findSample <- grep(subsetSamples[i], metadataSampleList)
    siphonMetadata[i,] <- metadata[findSample,]
    }
  }
  
  return(siphonMetadata)
  
}
  
