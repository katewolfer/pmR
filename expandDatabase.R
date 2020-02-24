expandDatabase <- function(df, ppm, polarityMode) {

  #df <- read.csv("small_metabolites_database_v1_07Jan2019.csv", check.names = FALSE, stringsAsFactors = FALSE)

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
                  "[M+2H]+",
                  "[M+NH4]+",
                  "[M+Na]+",
                  "[M+CH3OH+H]+",
                  "[M+K]+",
                  "[M+ACN+H]+",
                  "[M+IsoProp+H]+",
                  "[M+H+H2O]+",
                  "[M+Li]+",
                  "[2M+H]+",
                  "[M+H-H2O]+")

  } else if (polarityMode == "neg") {

    polarity <- c("[M-H]-",
                  "[M-H2O-H]-",
                  "[M+Cl]-",
                  "[M+FA-H]-",
                  "[M+CH3COO]-",
                  "[M+F]-",
                  "[2M-H]-")

  }

  ## masses of neutral conjugates
  sulf <- S+(O*3)
  gluc <- (C*6)+(H*8)+(O*6)

  ## number of basic adducts to calculate
  rowsToAdd <- length(polarity)
  dimsFreshDB <- nrow(df)*rowsToAdd

  ## add number of rows for features which need sulf/gluc added
  findConjugates <- length(which(df$type == "bacterial"))
  getConjugates <- rowsToAdd*(2*findConjugates)
  dimsFreshDB <- dimsFreshDB + getConjugates

  print(dimsFreshDB)

  ## create the new database
  newColNames <- c("metabolite","fullName","name","monoisotopic","formula","type",
                   "Metlin","notes","adduct","conjugate", "adduct_mass","lower_ppm","upper_ppm")

  expDF <- data.frame(matrix(NA, nrow = dimsFreshDB, ncol = length(newColNames)))
  colnames(expDF) <- newColNames

  ## populate database

  rowCounter <- 1

  for (i in 1:nrow(df)) {
    ## select mass to calculate
    M <- df$monoisotopic[i]

    if (df$type[i] != "bacterial"){

      toAdd <- rowsToAdd-1

      ## populate with feature details
      expDF[c(rowCounter:(rowCounter+toAdd)),c(1:8)] <- df[i,]
      expDF$adduct[c(rowCounter:(rowCounter+toAdd))] <- polarity

      if (polarityMode == "pos") {

        ## calculate adduct masses
        adductMass <- c(M+H-e,
                        M+(2*H)-e,
                        M+N+(H*4)-e,
                        M+Na-e,
                        M+C+(H*4)+O-e,
                        M+K-e,
                        M+(C*2)+(H*3)+N-e,
                        M+(C*3)+(H*8)+O-e,
                        M+(H*3)+O-e,
                        M+Li-e,
                        (2*M)+H-e,
                        M-H-O-e)

      } else if (polarityMode == "neg") {

        ## calculate adduct masses
        adductMass <- c(M-H+e,
                        M-((2*H) + (2*O))+e,
                        M+Cl+e,
                        M+(C+H+(2*O))+e,
                        M+((2*C)+(3*H)+(2*O))+e,
                        M+Fl+e,
                        (2*M)-H+e)

      }

      ## populate adduct masses
      expDF$adduct_mass[c(rowCounter:(rowCounter+toAdd))] <- adductMass

      ## add ppm range values
      massRange <- adductMass*(ppm/1000000)
      expDF$lower_ppm[c(rowCounter:(rowCounter+toAdd))] <- adductMass - massRange
      expDF$upper_ppm[c(rowCounter:(rowCounter+toAdd))] <- adductMass + massRange
      expDF$conjugate[c(rowCounter:(rowCounter+toAdd))] <- "none"
      # lower_ppm <- adductMass - massRange
      # upper_ppm <- adductMass  + massRange
      rowCounter <- rowCounter+toAdd

    } else {toAdd <- 3*rowsToAdd-1

    ## populate with feature details
    expDF[c(rowCounter:(rowCounter+toAdd)),c(1:8)] <- df[i,]
    expDF$adduct[c(rowCounter:(rowCounter+toAdd))] <- c(polarity,polarity,polarity)

    ## calculate adduct masses
    if (polarityMode == "pos") {

      ## calculate adduct masses
      adductMass <- c(M+H-e,
                      M+(2*H)-e,
                      M+N+(H*4)-e,
                      M+Na-e,
                      M+C+(H*4)+O-e,
                      M+K-e,
                      M+(C*2)+(H*3)+N-e,
                      M+(C*3)+(H*8)+O-e,
                      M+(H*3)+O-e,
                      M+Li-e,
                      (2*M)+H-e,
                      M-H-O-e)
      ## calculate sulfate conjuate masses
      sulfateMass <- c(M+sulf+H-e,
                       M+sulf+(2*H)-e,
                       M+sulf+N+(H*4)-e,
                       M+sulf+Na-e,
                       M+sulf+C+(H*4)+O-e,
                       M+sulf+K-e,
                       M+sulf+(C*2)+(H*3)+N-e,
                       M+sulf+(C*3)+(H*8)+O-e,
                       M+sulf+(H*3)+O-e,
                       M+sulf+Li-e,
                       (2*M)+sulf+H-e,
                       M+sulf-H-O-e)

      ## calculate glucuronide conjuate masses
      glucuronideMass <- c(M+gluc+H-e,
                           M+gluc+(2*H)-e,
                           M+gluc+N+(H*4)-e,
                           M+gluc+Na-e,
                           M+gluc+C+(H*4)+O-e,
                           M+gluc+K-e,
                           M+gluc+(C*2)+(H*3)+N-e,
                           M+gluc+(C*3)+(H*8)+O-e,
                           M+gluc+(H*3)+O-e,
                           M+gluc+Li-e,
                           (2*M)+gluc+H-e,
                           M+gluc-H-O-e)

    } else if (polarityMode == "neg") {

      ## calculate adduct masses
      adductMass <- c(M-H+e,
                      M-((2*H) + (2*O))+e,
                      M+Cl+e,
                      M+(C+H+(2*O))+e,
                      M+((2*C)+(3*H)+(2*O))+e,
                      M+Fl+e,
                      (2*M)-H+e)

      ## calculate sulfate conjuate masses
      sulfateMass <- c(M+sulf-H+e,
                       M+sulf-((2*H) + (2*O))+e,
                       M+sulf+Cl+e,
                       M+sulf+(C+H+(2*O))+e,
                       M+sulf+((2*C)+(3*H)+(2*O))+e,
                       M+sulf+Fl+e,
                       (2*(M+sulf))-H+e)

      ## calculate glucuronide conjuate masses
      glucuronideMass <- c(M+gluc-H+e,
                           M+gluc-((2*H) + (2*O))+e,
                           M+gluc+Cl+e,
                           M+gluc+(C+H+(2*O))+e,
                           M+gluc+((2*C)+(3*H)+(2*O))+e,
                           M+gluc+Fl+e,
                           (2*(M+gluc))-H+e)

    }

    ## string together and populate
    collectMass <- c(adductMass,sulfateMass,glucuronideMass)
    expDF$adduct_mass[c(rowCounter:(rowCounter+toAdd))] <- collectMass

    ## add ppm range values
    massRange <- adductMass*(ppm/1000000)
    expDF$lower_ppm[c(rowCounter:(rowCounter+toAdd))] <- collectMass - massRange
    expDF$upper_ppm[c(rowCounter:(rowCounter+toAdd))] <- collectMass + massRange
    conjugateList <- c(replicate(rowsToAdd, "none"),replicate(rowsToAdd, "sulfate"),replicate(rowsToAdd, "glucuronide"))
    expDF$conjugate[c(rowCounter:(rowCounter+toAdd))] <- conjugateList
    rowCounter <- rowCounter+toAdd
    }
  }

  expDF <- expDF[-which(is.na(expDF$metabolite) == TRUE),]

  return(expDF)
}
