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
    beta <- dim(Z)[1] / dim(Z)[2]
    omega <- 0.56*beta^3 - 0.95*beta^2 + 1.82*beta + 1.43
    rank = sum(USV$d >= omega*median(USV$d))
  }
  else if (rank < 0) {
    rank = max(dim(Z))
  }
  USV$u <- USV$u[,1:rank]
  USV$v <- USV$v[,1:rank]
  USV$d <- USV$d[1:rank]
  return(USV)
}

# Compute the dynamic mode decomposition (DMD) of a snapshot matrix Z and
# return the DMD modes Phi and frequencies omega.
compute_dmd <- function(Z, rank=0, delta.t=1.0, discard_pairs=TRUE){
  
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
  
  # By default, keep only one of a conjugate pairs by looking for nonnegative
  # imaginary components of eigenvalues.
  if (discard_pairs){
    keep <- which(Im(omega)>=0)
    omega <- omega[keep]
    Phi <- Phi[,keep]
  }
  
  dmd <- list(Phi=Phi, omega=omega, rank=dim(omega), delta.t = delta.t)
  return(dmd)
}

# Produce the basis Psi 
compute_dmd_basis <- function(dmd, times){
  Phi <- dmd$Phi
  omega <- dmd$omega
  Psi <- matrix(0.0,  dim(Phi)[1]*length(times), dim(Phi)[2])
  for (r in 1:dim(Phi)[2]){
    temporal <- exp(omega[r]*times)
    spatial <- Phi[,r]
    Psi[, r] <- Re(kronecker(spatial, temporal))
    Psi[, r] <- Psi[, r] / sqrt(sum(Psi[, r]^2)) # normalize for good measure
  }
  return(Psi)
}


