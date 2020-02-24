xr <- xcmsRaw("07Feb2020_SC_PBN_pos_full_18.mzML",
              profstep = 1,
              profmethod = "bin",
              #profparam = list(),
              includeMSn=TRUE,
              mslevel=NULL,
              #scanrange=NULL
              )

bpis <- chromatogram("07Feb2020_SC_PBN_pos_full_18.mzML", aggregationFun = "max")
peaks <- findPeaks(xr, method="MS1")



## xcmsSet
xfrag <- xcmsFragments(xs)
