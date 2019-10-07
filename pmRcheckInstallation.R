pmRcheckInstallation <- function() {
  
  #https://stackoverflow.com/questions/9341635/check-for-installed-packages-before-running-install-packages
  
  pkgTest <- function(x)
  {
    if (!require(x,character.only = TRUE))
    {
      install.packages(x,dep=TRUE)
      if(!require(x,character.only = TRUE)) stop("Package not found")
    }
  }
  
}