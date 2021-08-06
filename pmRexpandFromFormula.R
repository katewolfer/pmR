pmRexpandFromFormula <- function(df, ppm, polarityMode) {

  #####################################
  ## pmR package                     ##
  ## Kate Wolfer, Universitaet Basel ##
  ## v1.0, 03 October 2019           ##
  #####################################

  ## monoisotopic masses of elements
  H <- 1.007825
  C <- 12.000000
  O <- 15.994915
  N <- 14.003074
  S <- 31.972072
  P <- 30.973763
  Na <- 22.98977
  Li <- 7.016005
  K <- 38.963708
  Cl <- 34.968853
  Fl <- 18.998403

  ## mass of electron
  e <- 0.00054857990907

  ## mass of neutron
  n <- 1.008664915885

  ## list of adducts
  if (polarityMode == "pos") {

    polarity <- c("[M+H]+",
                  "[M]+",
                  "[M+2H]+",
                  "[M+NH4]+",
                  "[M+Na]+",
                  "[M+CH3OH+H]+",
                  "[M+K]+",
                  "[M+ACN+H]+",
                  #"[M+IsoProp+H]+",
                  "[M+H+H2O]+",
                  #"[M+Li]+",
                  "[2M+H]+",
                  "[M+H-H2O]+")

  } else if (polarityMode == "neg") {

    polarity <- c("[M-H]-",
                  "[M]-",
                  "[M-H2O-H]-",
                  "[M+Cl]-",
                  "[M+FA-H]-",
                  "[M+CH3COO]-",
                  #"[M+F]-",
                  "[2M-H]-")

  }

  ## number of basic adducts to calculate
  rowsToAdd <- length(polarity)
  dimsFreshDB <- nrow(df)*rowsToAdd

  ## create the new database
  newColNames <- c("monoisotopic","adduct","adduct_mass","lower_ppm","upper_ppm")
  expDF <- data.frame(matrix(NA, nrow = dimsFreshDB, ncol = length(newColNames)+ncol(df)))
  colnames(expDF)[c(1:ncol(df))] <- colnames(df)
  colnames(expDF)[c((ncol(df)+1):(ncol(df)+(length(newColNames))))] <- newColNames

  ## populate database
  rowCounter <- 1
  source("formulaToMass.R")

  for (i in 1:nrow(df)) {

    df$monoisotopic[i] <- formulaToMass(df$formula[i])

    ## select the calculated mass
    M <- df$monoisotopic[i]

    ## row counter for the final table population
    toAdd <- rowsToAdd-1

    ## populate with feature details
    expDF$adduct[c(rowCounter:(rowCounter+toAdd))] <- polarity
    expDF[c(rowCounter:(rowCounter+toAdd)), c(1:ncol(df))] <- df[i,]

    if (polarityMode == "pos") {

      ## calculate adduct masses
      adductMass <- c(M+H-e,
                      M-e,
                      ((M+(2*H))/2)-e,
                      M+N+(H*4)-e,
                      M+Na-e,
                      M+C+(H*4)+O-e,
                      M+K-e,
                      M+(C*2)+(H*3)+N-e,
                      #M+(C*3)+(H*8)+O-e,
                      M+(H*3)+O-e,
                      #M+Li-e,
                      (2*M+H)-e,
                      M-H-O-e)

    } else if (polarityMode == "neg") {

      ## calculate adduct masses
      adductMass <- c(M-H+e,
                      M+e,
                      M-((2*H) + (2*O))+e,
                      M+Cl+e,
                      M+(C+H+(2*O))+e,
                      M+((2*C)+(3*H)+(2*O))+e,
                      #M+Fl+e,
                      (2*M-H)+e)

    }

    ## populate adduct masses
    expDF$adduct[c(rowCounter:(rowCounter+toAdd))] <- polarity
    expDF$adduct_mass[c(rowCounter:(rowCounter+toAdd))] <- adductMass


    ## add ppm range values
    massRange <- adductMass*(ppm/1000000)
    expDF$lower_ppm[c(rowCounter:(rowCounter+toAdd))] <- adductMass - massRange
    expDF$upper_ppm[c(rowCounter:(rowCounter+toAdd))] <- adductMass + massRange
    rowCounter <- rowCounter+rowsToAdd

  }

  cat(df$check == df$monoisotopic)
  return(expDF)

}

#fetchValues <- unlist(regmatches(s,gregexpr("[[:digit:]]+\\.*[[:digit:]]*",s)))
