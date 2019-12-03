## This file produces the snapshot matrix from the SST data 
## provided by STRbook.  it numbers the months by cardinal index,
## starting with 0 for the first month included in the data.

library("dplyr")
library("STRbook")

#data("SST_df", package = "STRbook")
#SST_snapshots <- filter(SST_df, Year <= 1997)  # filter the NA dates
#write.csv(SST_snapshots,"SST_df.csv", row.names = FALSE)


data("NOAA_df_1990", package = "STRbook")
NOAA_snapshots <- filter(NOAA_df_1990, year == 1993, month == 7, proc=="Tmax")  # filter the NA dates
write.csv(NOAA_snapshots,"NOAA_df.csv", row.names = FALSE)