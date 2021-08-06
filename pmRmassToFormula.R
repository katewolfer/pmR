## pmRmassToFormula

## the calculating formula from a given mass is a recursive subset sum problem
## https://link.springer.com/content/pdf/10.1007/s00453-007-0162-8.pdf

pmRmassToFormula <- function(mass, charge, ppm){

  library(reshape2)
  source("pmRrecursiveFormula.R")

  mass <- 613.1592 ## oxidised glutathione
  charge <- "pos"

  ## mass of electron
  e <- 0.00054857990907

  ## mass of neutron
  n <- 1.008664915885

  halides <- readline(prompt="Does the molecule contain halides (e.g. Cl) ")

  ## select general small molecule combinations
  elements <- c("C","H","O","N","S","P","Na","Cl")
  elementMasses <- c(12.000000,
                     1.007825,
                     15.994915,
                     14.003074,
                     31.972072,
                     30.973763,
                     22.989770,
                     34.968853)

  maxAtoms <- c(50,60,15,10,6,4,1,1)

  ## recalculate to accommodate the specified charge
  if(charge == "pos"){
    neutralMass <- mass + e
  } else if (charge == "neg"){
    neutralMass <- mass - e
  } else {
    stop("incorrect charge input: use pos or neg")}

  ## get min and max mass for ppm
  #ppm = ((m1-m2)/m2)/100000
  upper_ppm <- neutralMass + 0.005
  lower_ppm <- neutralMass - 0.005

  findSmallMaxAtoms <- which(maxAtoms < 20)
  getMaxAltSubs <- sum(elementMasses[findSmallMaxAtoms]*maxAtoms[findSmallMaxAtoms])

  ## sort of dynamic programming recursive approach to get values
  carbonCount <- c(0:maxAtoms[(1)])
  carbonCombos <- as.data.frame(cbind(paste0(elements[1], carbonCount), elementMasses[1]*carbonCount))
  carbonCombos[1,1] <- elements[1]
  carbonCombos[,2] <- as.numeric(carbonCombos[,2])
  subsetCarbon <- carbonCombos
  #subsetCarbon <- carbonCombos[-which(carbonCombos[,2] > upper_ppm),]
  CHTable <- pmRrecursiveFormula(upper_ppm, elements[2], elementMasses[2], maxAtoms[2], subsetCarbon)
  C1Table <- pmRrecursiveFormula(upper_ppm, elements[3], elementMasses[3], maxAtoms[3], CHTable)

  finalTable <- C1Table[which((C1Table$value) > lower_ppm),]

  # deal with each in turn
  C2Table <- C1Table# data.frame(matrix(NA, nrow = 0, ncol = 2))
  for (s in 4:6){
    atomTable <- NULL
    for (t in 1:maxAtoms[s]){
      NCheck <- C2Table$value + (elementMasses[s]*t)
      NRemovePos <- which(NCheck > (upper_ppm))
      tempTable <- C2Table[-NRemovePos,]
      tempTable$value <- tempTable$value + (elementMasses[s]*t)
      tempTable$atom <- paste0(tempTable$atom, elements[s], t)
      tempTable$atom[1] <- paste0(tempTable$atom, elements[s])
      atomTable <- rbind(atomTable, tempTable)
    }
    C2Table <- rbind(C2Table, atomTable)
  }

  ## adducts
  maxAdductMass <- sum(elementMasses[c(7:8)]*maxAtoms[c(7:8)])
  adducts <- c("Na","Cl")
  adductMasses <- c(22.989770,
                     34.968853)
  maxAdducts <- c(2,2)

  finalTable <- rbind(finalTable, C2Table)

  refineTable <- finalTable[which((finalTable$value) < upper_ppm),]
  refineTable <- refineTable[which((refineTable$value) > lower_ppm),]

  keepSpaceAdductHits <- finalTable[which((finalTable$value + maxAdductMass) < upper_ppm), ]

  C3Table <- C2Table# data.frame(matrix(NA, nrow = 0, ncol = 2))
  for (s in 1:2){
    atomTable <- NULL
    for (t in 1:maxAdducts[s]){
      NCheck <- C3Table$value + (adductMasses[s]*t)
      NRemovePos <- which(NCheck > (upper_ppm))
      tempTable <- C3Table[-NRemovePos,]
      tempTable$value <- tempTable$value + (adductMasses[s]*t)
      tempTable$atom <- paste0(tempTable$atom, adducts[s], t)
      tempTable$atom[1] <- paste0(tempTable$atom, elements[s])
      atomTable <- rbind(atomTable, tempTable)
    }
    C3Table <- rbind(C3Table, atomTable)
  }

  C3Table  <- C3Table[which((C3Table$value) < upper_ppm),]
  C3Table  <- C3Table[which((C3Table$value) > lower_ppm),]

  totalTable <- rbind(refineTable, C3Table)
  totalTable$ppm <- round(((totalTable$value-neutralMass)/neutralMass)*1000000,3)
  totalTable$ppmAbs <- abs(totalTable$ppm)

totalTable <- totalTable[order(totalTable$ppmAbs, decreasing = FALSE),]




}
