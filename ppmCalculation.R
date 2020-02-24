##################################################
## ppm function for calculating mass accuracy   ##
## Kate Wolfer, Institute of Hepatology, London ##
## 06 December 2018                             ##
##################################################

## INPUT
## targetValues: a vector of target mass (m/z) values representing a list of MS features of interest
## observedValues: a vector of the observed/measured m/z values

## OUTPUT
## appendedTable: dataframe of target m/z values, lower mass bound, upper mass bound
##                measured mass value, calculated ppm value

ppmCalculation <- function(ppm, targetValues, observedValues) {
  lowerMass <- NULL
  upperMass <- NULL
  
  for (i in 1:length(targetValues)){
    massRange <- targetValues[i]*(ppm/1000000)
    lowerMass[i] <- targetValues[i] - massRange
    upperMass[i] <- targetValues[i] + massRange
  }
  
  
  checkPPM <- NULL
  filterTable <- NULL
  
  for (i in 1:nrow(importMasses)){
    for (j in 1:nrow(matchList)){
      checkPPM <- importMasses$mz[i] > matchList$lowerMass[j] && importMasses$mz[i] < matchList$upperMass[j]
      if (checkPPM == TRUE) {
        slotRow <- cbind(importMasses[i,],matchList[j,])
        filterTable <- rbind(filterTable, slotRow) 
      }
    }
  }
  
  
  appendedTable <- data.frame(targetValues,lowerMass,upperMass,observedValues,calculated_ppm)
  
  return(appendedTable)
  
}