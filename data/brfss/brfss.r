# codebook is here: https://www.cdc.gov/brfss/annual_data/2017/pdf/codebook17_llcp-v2-508.pdf
# data is here (weirdly the filename has spaces in it, so I renamed it)
# https://www.cdc.gov/brfss/annual_data/2017/files/LLCP2017XPT.zip

library(foreign)
df = read.xport("LLCP2017.XPT")
hist(df$SLEPTIM1[df$SLEPTIM1<=24])

