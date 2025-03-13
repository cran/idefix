# Derr_ucpp using InfoDes_cpp and det_cpp (for Modfed) 
Derr_ucpp <- function(par, des, n.alts) {
  info.des <- InfoDes_cpp(par, des, n.alts)
  detinfo <- det_cpp(info.des)
  if(is.nan(detinfo)){
    return(NaN)
  } else {
    ifelse((detinfo <= 0), return(NA), return(detinfo^(-1 / length(par))))
  }
}

# DerrS_P using InfoDes_cpp and det_cpp
DerrS.P_ucpp <- function(par, des, n.alts, i.cov) {
  info <- InfoDes_cpp(par = par, des = des, n_alts = n.alts)
  d.error <- det_cpp(info + i.cov)^(-1 / length(par))
  return(d.error)
}

# DerrC using cpp functions
DerrC_ucpp <- function(par, des, n.alts, i.cov) {
  info.des <- InfoDes_cpp(par, des, n.alts)
  detinfo <- det_cpp(info.des + i.cov)
  ifelse((detinfo <= 0), return(NA), return(detinfo^(-1 / length(par))))
}

# DBerrS.P using DerrS.P_cpp
DBerrS.P_ucpp <- function(des, par.draws, n.alts, i.cov, weights) {
  # Add alternative specific constants if necessary
  # For each draw calculate D-error.
  d.errors <- apply(par.draws, 1, DerrS.P_ucpp, des, n.alts, i.cov)
  w.d.errors <- d.errors * weights
  # DB-error. 
  db.error <- mean(w.d.errors, na.rm = TRUE)
  return(db.error)
}

# Aerror (for Modfed and CEA) 
Aerr_ucpp <- function(par, des, n.alts) {
  info.des <- round(InfoDes_cpp(par, des, n.alts), digits = 5) #rounding is to avoid floating point (infinitesimal number but above zero)
  # using the determinant as a check for the info. matrix invertibility
  if (det(info.des) < .Machine$double.eps^0.5) {
  # if (qr(info.des)$rank == ncol(info.des)) {
    aerror <- NaN
  } else {
    temp.error <- sum(diag(solve(info.des)))/length(par)
    if (temp.error < .Machine$double.eps) { #to avoid potential floating point problem (could result in an incorrect huge negative aErrors)
      aerror <- NA
    } else {
      aerror <- temp.error
    }
  }
  return(aerror)
}

