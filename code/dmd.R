# This script contains functions to perform Dynamic Mode Decomposition
# and to produce basis vectors from that decomposition.

# Compute the singular value decomposition (SVD) of a data matrix
# according to a user-specified rank.
compute_svd <- function(Z, rank){
  
  # Perform the SVD and extract the desired columns/values
  USV <- svd(Z)
  
  # If the rank given is not positive, define it.
  if (rank == 0){
    # Using the (approximate) optimal truncation rule of omega(beta) * median(sigma)
    # from Davish and Gonoho, where beta = num_rows/num_columns.
    # Basicaly, this keeps the "right" number of singular values and vectors
    # the signal-to-noise ratio.
    beta <- min(dim(Z)) / max(dim(Z))
    omega <- 0.56*beta^3 - 0.95*beta^2 + 1.82*beta + 1.43
    rank = sum(USV$d >= omega*median(USV$d))
  }
  else if (rank < 0 || rank > min(dim(Z))) {
    rank = min(dim(Z))
  }
  stopifnot(rank > 0)
  # the matrix/vector calls handled the case when rank==1 since we still
  # need explicit matrix objects
  USV$u <- as.matrix(USV$u[,1:rank])
  USV$v <- as.matrix(USV$v[,1:rank])
  USV$d <- USV$d[1:rank]
  if (rank==1){
    USV$d <- as.matrix(USV$d)
  }
  
  return(USV)
}

# Compute the dynamic mode decomposition (DMD) of a snapshot matrix Z and
# return the DMD modes Phi and frequencies omega.
compute_dmd <- function(Z, rank=0, t.0=0.0, delta.t=1.0, keep_pairs=0){
  
  # Simplify things by separating snapshots into + and -
  num_time <- dim(Z)[2]
  Z.plus <- Z[,2:num_time]
  Z.minus <- Z[,1:(num_time-1)]
  
  # Compute the SVD 
  USV <- compute_svd(Z.minus, rank)
  U.r <- USV$u
  V.r <- USV$v
  S.r <- USV$d

  # Compute the low-rank operator and its spectrum
  A.tilde <- Conj(t(U.r)) %*% Z.plus %*% V.r %*% diag(1.0/S.r)
  WL <- eigen(A.tilde)
  W <- WL$vectors
  Lambda <- WL$values
  
  # Compute the DMD modes and frequencies assuming a fixed time step 
  # in the appropriate units.
  Phi <- Z.plus %*% V.r %*% diag(1.0/S.r) %*% W
  omega <- log(Lambda)/delta.t
  
  # By default, keep only one of a conjugate pair by looking for nonnegative
  # imaginary components of eigenvalues.
  if (keep_pairs==0){
    keep <- which(Im(omega)>=0)
    omega <- omega[keep]
    Phi <- Phi[,keep]
  }
  Phi <- as.matrix(Phi)
  
  dmd <- list(Phi=Phi, omega=omega, rank=length(omega), t.0=t.0, delta.t=delta.t, keep_pairs=keep_pairs)
  return(dmd)
}

# Produce the basis Psi 
compute_dmd_basis <- function(dmd, times){
  Phi <- dmd$Phi
  omega <- dmd$omega
  # Construct the basis but separate the real and imaginary components.  Here,
  # create a list of 1's and 0's to indicate which should be real and which
  # should be imaginary.  The reals will be any mode for which Im(omega)>=0
  if (dmd$keep_pairs==2) {
    real <- Im(omega) >= 0.0
  }
  else {
    real <- rep(TRUE, length(omega))
  }
  times <- times - dmd$t.0
  Psi <- matrix(0.0,  dim(Phi)[1]*length(times), dim(Phi)[2])
  for (r in 1:dim(Phi)[2]){
    temporal <- exp(omega[r]*times)
    spatial <- Phi[,r]
    if (real[r]){
      Psi[, r] <- Re(kronecker(spatial, temporal))
    }
    else {
      Psi[, r] <- Im(kronecker(spatial, temporal))
    }
    Psi[, r] <- Psi[, r] / sqrt(sum(Psi[, r])^2) # normalize for good measure
  }
  return(Psi)
}

# Extend an existing dataframe with the DMD basis.  This assumes that the
# observables are on a fixed space-time grid in which the spatial coordinates are
# consistent with those originally used for the DMD basis.  Moreover, the data should be 
# organized in the same, logically time-wide format, i.e, all observations for a given
# location form a contiguous block of rows in the same temporal order.  
# Note, the times included may differ from
# those originally used to construct the DMD modes and frequencies.  However, all times
# must be defined relative to the original times (and here, everything is assumed to
# start at zero).
#
# This returns the dataframe with additional columns named Psi.1, Psi.2, ..., Psi.R, where
# R is the number of basis vectors.
df_add_dmd <- function(df, dmd){
  times <- unique(df$t)
  Psi <- compute_dmd_basis(dmd, times)
  df <- cbind(df, Psi=Psi)
  return(df) 
}

# Create linear model based on DMD of rank r
lm_with_dmd <- function(df, rank=0, keep_pairs=0){
  num_locs <- length(unique(df$id))
  num_time <- length(df$id)/num_locs
  times <- sort(unique(df$t))
  stopifnot(length(times) ==  num_time)
  delta.t <- times[2]-times[1]
  t.0 <- times[1]
  Z <- t(matrix(as.numeric(df$z), num_time, num_locs, byrow = FALSE))
  dmd <- compute_dmd(Z, rank, t.0, delta.t, keep_pairs)
  df <- df_add_dmd(subset(df, select=c(z, t)), dmd)
  return(list(model=lm(z~.-t-1, data=df), dmd=dmd)) # remove t and intercept
}

# Create a prediction based on a linear model constructed with DMD basis vectors.
predict_with_dmd <- function(df, model, dmd){
  df <- df_add_dmd(subset(df, select=c(z, t)), dmd)
  pred <- predict(model, newdata=df, interval="prediction")
  return(pred)
}
