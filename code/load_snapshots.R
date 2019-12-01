## This file produces the snapshot matrix from the SST data 
## provided by STRbook.  it numbers the months by cardinal index,
## starting with 0 for the first month included in the data.

library("dplyr")
library("STRbook")

data("SST_df", package = "STRbook")
SST_snapshots <- filter(SST_df, Year <= 1997)  # filter the NA dates
write.csv(SST_snapshots,"SST_df.csv", row.names = FALSE)
