library("dplyr")
library("ggplot2")

source("dmd.R")


####################
# Do some analysis #
####################

# Load the data.
df <- read.csv(file="NOAA_df.csv", header=TRUE)

# Construct the snapshot matrix Z.
# Rows are locations and columns are times.
num_locs <- length(unique(df$id))
num_time <- length(df$id)/num_locs
Z <- t(matrix(as.numeric(df$z), num_time, num_locs, byrow = FALSE))

# First several days to use for fitting and basis generation for 
# testing forecasting
t.part <- 28

# Do the DMD with all of the snapshots and assume a delta t of 1 day.
# Note, time is assumed to start at 0 in this simple implementation.
dmd_full <- compute_dmd(Z, rank=5, delta.t=1.0)

# Do the DMD part of the snapshots and assume a delta t of 1 day.
# Note, time is assumed to start at 0 in this simple implementation.
dmd_part <- compute_dmd(Z[, 1:t.part], rank=5, delta.t=1.0)

# Add translated time to data frame
df <- mutate(df, t=df$julian-min(df$julian))

# Helper function to produce a dataframe for the DMD model
make_dmd_df <- function(df, dmd){
  times <- sort(unique(df$t))
  # Generate the design matrix and add to dataframe
  X <- compute_dmd_basis(dmd, times)
  #X <- generate_design_matrix(dmd$Phi, dmd$omega, times)
  df_dmd <- mutate(df, a=X[,1], b=X[,2], c=X[,3])
  return(df_dmd)
}

# Make DMD model based on all data set and fit against all data
model_dmd_full <- lm(z~a+b+c-1, data=make_dmd_df(df, dmd_full))
# Make DMD model based on first df21 days and fit against those same days
model_dmd_part <- lm(z~a+b+c-1, data=make_dmd_df(filter(df, t<=t.part), dmd_part))
# Make second-order polynomial model based on all data
model_poly_full <- lm(z~(t+lat+lon)^2, data=filter(df, t<=t.part))
# Make second-order polynomial model based on half data and fit
model_poly_part <- lm(z~(t+lat+lon)^2, data=df)

# Predictions for all of July at spatial point
p_dmd_full <- cbind(df, predict(model_dmd_full, newdata=make_dmd_df(df, dmd_full), interval="prediction"))
p_dmd_part <- cbind(df, predict(model_dmd_part, newdata=make_dmd_df(df, dmd_part), interval="prediction"))
p_poly_full <- cbind(df, predict(model_poly_full, newdata=df, interval="prediction"))
p_poly_part <- cbind(df, predict(model_poly_part, newdata=df, interval="prediction"))

library("ggplot2")

make_prediction_plot <- function(pred, location){

  myid <- unique(pred$id)[location]
  plot_df <- filter(pred, id==myid)
  p <- ggplot(plot_df, aes(t, z)) +
       geom_point() +
       geom_line(aes(y = fit), color="blue") +
       geom_line(aes(y = upr), color = "red", linetype = "dashed") +
       geom_line(aes(y = lwr), color = "red", linetype = "dashed") 

  return(p)
}

location <- 99
plot_dmd_full <-make_prediction_plot(p_dmd_full, location)
plot_dmd_part <-make_prediction_plot(p_dmd_part, location)
plot_poly_full <-make_prediction_plot(p_poly_full, location)
plot_poly_part <-make_prediction_plot(p_poly_part, location)

ggsave('plot_dmd_full.pdf', plot_dmd_full)
ggsave('plot_dmd_part.pdf', plot_dmd_part)
ggsave('plot_poly_full.pdf', plot_poly_full)
ggsave('plot_poly_part.pdf', plot_poly_part)




