## Expanding lists of features - calculates monoisotopic masses of given formulas,
## and seleted adduct masses with ppm brackets, then ouputs the list

source("expandList.R")

df <- read.csv("criegees_19Feb2020.csv", check.names = FALSE,
               stringsAsFactors = FALSE)

polarityMode <- "pos"
ppm <- 5
getList <- expandList(df, ppm, polarityMode)
