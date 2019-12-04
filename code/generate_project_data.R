# This file filters data from the STRbook package for use 
# in the project.  I've separated this so one does not have
# to install the STRbook package to do repeat the analysis.

library("dplyr")
library("STRbook")

# Load the data and filter to extract July 1993 maximum temperatures
data("NOAA_df_1990", package = "STRbook")
NOAA_snapshots <- filter(NOAA_df_1990, year==1993, month==7, proc=="Tmax")

# Now, we need to clean up the data to include only those stations that have
# measurements for each day of the month.  This greatly simplifies the 
# DMD analysis, which is based on a dense [space, time] snapshot matrix.

# There are 31 days in July, so search for all those station id's that
# have 31 measurements
keep <- c()
for (station in unique(NOAA_snapshots$id)){
  if (length(filter(NOAA_snapshots,id==station)$julian)==31){
    keep <- append(keep, station)
  }
}

# Filter the snapshots again, keeping just those stations with measurements
# on each day.  
NOAA_snapshots2 <- filter(NOAA_snapshots,id%in%keep)

# Save to CSV for processing in the tutorial.
write.csv(NOAA_snapshots2,"NOAA_df.csv", row.names = FALSE)