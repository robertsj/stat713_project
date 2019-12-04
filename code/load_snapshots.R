library("dplyr")
library("ggplot2")

# Compute the dynamic mode decomposition of a snapshot matrix Z and
# keep, at most, rank modes.
compute_dmd <- function(Z, rank, delta.t=1.0){
  
  # Simplify things by separating snapshots into + and -
  num_time <- dim(Z)[2]
  Z.plus <- Z[,2:num_time]
  Z.minus <- Z[,1:(num_time-1)]
  
  # Perform the SVD and extract the desired columns/values
  USV <- svd(Z.minus)
  U.r <- USV$u[,1:rank]
  V.r <- USV$v[,1:rank]
  S.r <- USV$d[1:rank]
  
  # Compute the low-rank operator and its spectrum
  A.tilde <- Conj(t(U.r)) %*% Z.plus %*% V.r %*% diag(1.0/S.r)
  WL <- eigen(A.tilde)
  W <- WL$vectors
  Lambda <- WL$values
  
  # Compute the DMD modes and frequencies assuming a fixed time step 
  # in the appropriate units.
  Phi <- Z.plus %*% V.r %*% diag(1.0/S.r) %*% W
  delta.t <- 1.0
  omega <- log(Lambda)/delta.t
  
  # Keep only one of a conjugate pairs by looking for nonnegative
  # imaginary ccomponents of eigenvalues
  keep <- which(Im(omega)>=0)
  omega <- omega[keep]
  Phi <- Phi[,keep]
  dmd <- list(Phi=Phi, omega=omega, rank=dim(omega))
  return(dmd)
}

# Generate the basis vectors for arbitrary times given the DMD
# modes and frequencies.
generate_design_matrix <- function(Phi, omega, times){
  X <- matrix(0.0,  dim(Phi)[1]*length(times), dim(Phi)[2])
  for (r in 1:dim(Phi)[2]){
    temporal <- exp(omega[r]*times)
    spatial <- Phi[,r]
    X[, r] <- Re(kronecker(spatial, temporal))
    X[, r] <- X[, r] / sqrt(sum(X[, r]^2))
  }
  return(X)
}

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
t.part <- 21

# Do the DMD with all of the snapshots and assume a delta t of 1 day.
# Note, time is assumed to start at 0 in this simple implementation.
dmd_full <- compute_dmd(Z, rank=5, delta.t=1.0)

# Do the DMD part of the snapshots and assume a delta t of 1 day.
# Note, time is assumed to start at 0 in this simple implementation.
dmd_part <- compute_dmd(Z[, 1:15], rank=5, delta.t=1.0)

# Add translated time to data frame
df <- mutate(df, t=df$julian-min(df$julian))

# Helper function to produce a dataframe for the DMD model
make_dmd_df <- function(df, dmd){
  times <- sort(unique(df$t))
  # Generate the design matrix and add to dataframe
  X <- generate_design_matrix(dmd$Phi, dmd$omega, times)
  df_dmd <- mutate(df, a=X[,1], b=X[,2], c=X[,3])
  return(df_dmd)
}

# Make DMD model based on all data set and fit against all data
model_dmd_full <- lm(z~a+b+c-1, data=make_dmd_df(df, dmd_full))
# Make DMD model based on first 21 days and fit against those same days
model_dmd_part <- lm(z~a+b+c-1, data=make_dmd_df(filter(df, t<=t.part), dmd_part))
# Make second-order polynomial model based on all data
model_poly_full <- lm(z~(t+lat+lon)^2, data=filter(df, t<=t.part))
# Make second-order polynomial model based on half data and fit
model_poly_part <- lm(z~(t+lat+lon)^2, data=df)

# Predictions for all of July at spatial point
location <- 30
id <- unique(df$id)[location]
p_dmd_full <- predict(model_dmd_full, newdata=filter(make_dmd_df(df, dmd_full), id==id), interval="prediction")
p_dmd_part <- predict(model_dmd_part, newdata=filter(make_dmd_df(df, dmd_part), id==id), interval="prediction")
p_poly_full <- predict(model_poly_full, newdata=filter(df, id==id), interval="prediction")
p_poly_part <- predict(model_poly_part, newdata=filter(df, id==id), interval="prediction")

library("ggplot2")

make_prediction_plot <- function(df, p, location){
  id <- unique(df$id)[location]
  p <- ggplot(cbind(df, p), aes(t, z)) +
    geom_point() +
    stat_smooth(method = lm)
  p + geom_line(aes(y = lwr), color = "red", linetype = "dashed") +
    geom_line(aes(y = upr), color = "red", linetype = "dashed")
  return(p)
}

location <- 1
plot_dmd_full <-make_prediction_plot(filter(df, id==id), p_dmd_full, location)
print(plot_dmd_full)





