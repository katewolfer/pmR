pmRXCMSparameters <- function() {

  #####################################
  ## pmR package                     ##
  ## Kate Wolfer, Universitaet Basel ##
  ## v1.0, 03 October 2019           ##
  #####################################

  pm <- list()
  pm[1] <- 1 # CPWmin = 1
  pm[2] <- 25 # CPWmax = 25
  pm[3] <- 3 # ppm = 3
  pm[4] <- 3 # xsnthresh = 3
  pm[5] <- F # LM = F
  pm[6] <- 2 # integrate = 2
  pm[7] <- 3 # nSlaves = 3
  pm[8] <- 10 # RTerror = 10
  pm[9] <- 0.01 # MZerror = 0.01
  pm[10] <- c(3,1000) # prefilter = c(3,1000)
  pm[11] <- 6 # bw = 6
  pm[12] <- 4 # nCores = 4
  pm[13] <- SnowParam(as.numeric(unlist(pm[12]))) # BPPARAM <- SnowParam(nCores)
  pm[14] <- 0.5 # minFrac.val = 0.50
  pm[15] <- 0.2 # cv.threshold = 0.2

  return(pm)

}
