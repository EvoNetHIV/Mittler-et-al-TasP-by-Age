this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)
tiff("TASP_Exp48_47_absdiff_doubled.tif", res=600, compression = "lzw", height=8, width=7.322917, units="in")