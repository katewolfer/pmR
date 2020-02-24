formulaToMass <- function(formulas) {

  ## Function to calculate any formula containing C H N O S only

  options(warn = -1)

  ## select mass to calculate
  s <- formulas

  ## monoisotopic masses of elements
  H <- 1.007825
  C <- 12.000000
  O <- 15.994915
  N <- 14.003074
  S <- 31.972072
  P <- 30.973763

  ## mass of electron
  e <- 0.00054857990907

  ## mass of neutron
  n <- 1.008664915885

  ## set up mass calculations for each element
  tabulateFormula <- data.frame(matrix(NA, nrow = 1, ncol = 6))
  tabulateFormula[1,1] <- C
  tabulateFormula[1,2] <- H
  tabulateFormula[1,3] <- N
  tabulateFormula[1,4] <- O
  tabulateFormula[1,5] <- S
  tabulateFormula[1,6] <- P

  ## find the elements in the formula
  tabulateFormula[2,1] <- as.numeric(unlist(gregexpr(pattern ='C', s)))
  tabulateFormula[2,2] <- as.numeric(unlist(gregexpr(pattern ='H', s)))
  tabulateFormula[2,3] <- as.numeric(unlist(gregexpr(pattern ='N', s)))
  tabulateFormula[2,4] <- as.numeric(unlist(gregexpr(pattern ='O', s)))
  tabulateFormula[2,5] <- as.numeric(unlist(gregexpr(pattern ='S', s)))
  tabulateFormula[2,6] <- as.numeric(unlist(gregexpr(pattern ='P', s)))

  ## if any elements not represented, ensure this is flagged
  for (n in  1:6){
    if(tabulateFormula[2,n] < 0){
      tabulateFormula[2,n] <- 0
    }
  }

  ## find each eement in the formula and extract the relevant numbers
  listElements <- c("C","H","N","O","S","P")

  for (g in 1:6){

    locateE <- as.numeric(unlist(gregexpr(pattern =listElements[g], s)))

    if(locateE < 0){
      tabulateFormula[2,g] <- 0
    } else {

      findE <- as.numeric(substr(s,locateE+1,locateE+3))

      if (is.na(findE) == TRUE){
        findE <- as.numeric(substr(s,locateE+1,locateE+2))
      }

      if (is.na(findE) == TRUE){
        findE <- as.numeric(substr(s,locateE+1,locateE+1))
      }

      if (is.na(findE) == TRUE){
        findE <- 1
      }

      tabulateFormula[2,g] <- as.numeric(findE)
    }
  }

  ## add the numbers into the table
  for (n in  1:6){
    tabulateFormula[3,n] <- tabulateFormula[1,n]*tabulateFormula[2,n]
  }

  ## calculate the total mass of each element
  getMass <- sum(tabulateFormula[3,])
  return(getMass)

}
