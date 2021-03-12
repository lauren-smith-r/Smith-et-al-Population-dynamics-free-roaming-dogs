periodic_decay_kernal <- function(dmat, etasq, decaytime, period, sigma){
  out <- matrix(NA, ncol=ncol(dmat), nrow=nrow(dmat))
  for(i in 1:(ncol(dmat)-1)){
    for(j in (i+1):(ncol(dmat))){
      out[i,j] = etasq * exp(-0.5 * dmat[i,j]/decaytime - (2*sin(pi*(dmat[i,j]))^2)/period)
      out[j,i] = out[i,j]
    }
  }
  for(k in 1:ncol(dmat)){
    out[k,k] <- etasq + sigma
  }
  out
}

sqexp_cov_kernel <- function(dmat, etasq, rhosq, sigma){
  out <- matrix(NA, ncol=ncol(dmat), nrow = nrow(dmat))
  for(i in 1:(ncol(dmat)-1)){
    for(j in (i+1):(ncol(dmat))){
      out[i,j] = etasq * exp(-rhosq * (dmat[i,j]^2));
      out[j,i] = out[i,j];
    }
  }
  for(k in 1:ncol(dmat)){
    out[k,k] <- etasq + sigma
  }
  return(out)
}

k_cov_corr <- function(K){
  round(cov2cor(K),3)
}

