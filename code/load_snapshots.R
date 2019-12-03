# Load the data
df <- read.csv(file="NOAA_df.csv", header=TRUE)

# Construct the snapshot matrix Z.
# Rows are locations and columns are times.
num_locs <- length(unique(df$id))
num_time <- length(df$id)/num_locs
Z <- t(matrix(as.numeric(df$z), num_time, num_locs, byrow = FALSE))


# Compute the dynamic mode decomposition of a snapshot matrix Z and
# keep, at most, rank modes.
compute_dmd <- function(Z, rank, delta.t=){
  
  # Simplify things by separating snapshots into + and -
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
  
  # Compute the DMD modes and frequencies
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

#
times <- 0:30

