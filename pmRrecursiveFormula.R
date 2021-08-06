pmRrecursiveFormula <- function(mass, element, elementMass, maxAtoms, formulaTable){

  ## add hydrogen
  atomCount <- c(0:maxAtoms)
  atomCombos <- as.data.frame(cbind(paste0(element, atomCount), elementMass*atomCount))
  atomCombos[1,1] <- element
  atomCombos[,2] <- as.numeric(atomCombos[,2])
  subsetAtom <- atomCombos

  addAtom <- data.frame(matrix(NA, nrow = nrow(formulaTable), ncol = nrow(subsetAtom)))

  Sys.time()
  for (j in 1:nrow(addAtom)){
    addAtom[j,] <- as.numeric(subsetAtom[,2]) + as.numeric(formulaTable[j,2])
  }
  Sys.time()

  colnames(addAtom) <- atomCombos[,1]
  addAtom <- cbind(formulaTable[,1], addAtom)
  colnames(addAtom)[1] <- "atom"

  meltCH <- melt(addAtom, id.vars = "atom")
  meltCH[,1] <- paste0(meltCH[,1], meltCH[,2])
  meltCH <- meltCH[,-2]
  meltCH <- meltCH[-which(as.numeric(meltCH$value) > mass),]

  return(meltCH)

}
