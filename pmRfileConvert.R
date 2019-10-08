pmRfileConvert <- function(proteowizardPath,targetFolder, saveFolder) {
  
  #####################################
  ## pmR package                     ##
  ## Kate Wolfer, Universitaet Basel ##
  ## v1.0, 03 October 2019           ##
  #####################################
  
  ## function to automatically convert raw mass spec files to mzML format
  ## using Proteowizard
  
  #dir.create(paste0(targetFolder, '\\converted_data'))
  #saveFolder       <- targetFolder
  #saveFolder        <- paste0(targetFolder, '\\converted_data')
  
  # Files to convert
  files       <- dir(targetFolder, pattern='*.raw')
  targetFiles <- c()
  for (f in files){
    targetFiles <- c(targetFiles, normalizePath(file.path(targetFolder, f)))
  }
  
  ## Set command line
  inputCmd     <- c()
  for (filePath in targetFiles){
    inputCmd     <- c(inputCmd,     paste('msconvert.exe', 
                                          filePath, 
                                          '-o', 
                                          saveFolder, 
                                          #'--zlib --filter "scanEvent 1"', 
                                          '--zlib --filter "msLevel2-"',
                                          sep=' '))
  }
  
  #print(proteowizardPath)
  #print(targetFiles)
  #print(saveFolder)
  #print(inputCmd)
  #print(inputCmd_100)
  
  ## Run
  setwd(proteowizardPath)
  i   <- 0
  nb  <- length(inputCmd)
  
  for (cmd in inputCmd){
    system2('cmd.exe', wait=TRUE, input=cmd)
    i <- i+1
    print('\n')
    print(paste('------------', i, '/', nb, 'done ------------'))
    print('\n')
  }
  
  
}