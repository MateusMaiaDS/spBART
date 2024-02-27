# Script to calculate a likely range of tau_beta
# and then approximate this as a gamma distribution
# Model is y ~ N (B beta, tau^-1)
# with beta ~ N(0, tau_beta^-1)
# so y ~ N(0, tau_beta^-1 * BB^T + I / tau)
# but y is also scaled to be in the range (-0.5, 0.5)
# and x is in the range (0, 1)

# Want to find the largest tau_beta so that var(y) < eps
# and the smallest tau_beta that maximises var(y)

# First create the basis functions
bspline <- function(x, xl, xr, ndx, bdeg){
  dx <- (xr - xl) / ndx
  knots <- seq(xl-bdeg*dx, xr + bdeg*dx, by = dx)
  B <- splines::spline.des(knots,x,bdeg + 1, 0*x,
                           outer.ok = TRUE)$design
}

# # Getting the maximum and the minimum of x
# set.seed(123)
# x <- sort(runif(1000))
# xr <- max(x)
# xl <- min(x)
# bdeg <- 3 # For a cubic spline
# nIknots_sim <- 30
# ndx <- nIknots_sim + 1 # Number of intervals within the range, i.e: number of internal knots + 1.
# B <- bspline(x = x,xl = xl,xr = xr, ndx = ndx, bdeg = bdeg)

# Now create a function which takes B
# and returns the parameters of a gamma distribution
return_gamma <- function(B, d = 0, calculate_range = FALSE) {

  diff2 <- function(x, differences = 0) {
    if (differences == 0) {
      return(x)
    } else {
      return(diff(x, differences = differences))
    }
  }
  # First find the minimum value of tau_beta
  # I think this is equivalent to the problem:
  # Find max(var(beta))
  # subject to range(B beta) in (-0.5,0.5)
  # Need to create a minimisation function for this
  min_opt <- function(beta) {
    lambda <- 10
    neg_var <- -var(beta)
    yhat <- B %*% beta
    range_yhat <- range(yhat)
    # Create a quadratic penalty away from the range
    if (range_yhat[1] < -0.5) {
      neg_var <- neg_var + lambda * (range_yhat[1] + 0.5)^2
    }
    if (range_yhat[2] > 0.5) {
      neg_var <- neg_var + lambda * (range_yhat[2] - 0.5)^2
    }
    return(neg_var)
  }
  min_opt(runif(ncol(B)))

  ans <- optim(par = rnorm(ncol(B)),
               # Not sure what these two should be set at
               lower = -1, upper = 1,
               fn = min_opt,
               method = "L-BFGS-B",
               control = list(maxit = 20000))
  beta <- ans$par
  # plot(x, B %*% beta, type = 'l')

  # So the minimum value of tau_beta
  min_tau_beta <- 1 / var(diff2(beta, differences = d))

  # Second find the max value of tau_beta
  # finding
  min_opt2 <- function(beta) {
    eps <- 0.01
    lambda <- 1 / eps^2
    var_beta <- var(beta)
    yhat <- B %*% beta
    var_yhat <- var(yhat)
    # Create a quadratic penalty away from the range
    if (var_yhat < eps) {
      var_beta <- var_beta + lambda * (var_yhat - eps)^2
    }
    return(var_beta)
  }
  # min_opt2(runif(ncol(B)))

  ans <- optim(par = rnorm(ncol(B)),
               # Not sure what these two should be set at
               lower = -1, upper = 1,
               fn = min_opt2,
               method = "L-BFGS-B",
               control = list(maxit = 20000))
  beta <- ans$par
  # plot(x, B %*% beta, type = 'l', ylim = range(B%*%beta))
  max_tau_beta <- 1 / var(diff2(beta, differences = d))

  # Finally use use min_tau_beta and max_tau_beta as 95% confidence
  # limits to obtain the parameters of a gamma distribution
  min_opt3 <- function(par) {
    shape <- par[1]
    rate <- par[2]
    alpha <- 0.95
    ucl <- qgamma(alpha + (1 - alpha)/2, shape = shape, rate = rate)
    lcl <- qgamma((1 - alpha)/2, shape = shape, rate = rate)
    return((ucl - max_tau_beta)^2 + (lcl - min_tau_beta)^2)
  }
  # min_opt3(runif(2))

  # Couldn't find a way of running B-FGS-B on this one
  # without it failing so gives some warnings
  sh_ra <- optim(par = c(1, 1),
                 fn = min_opt3)
  # This model gonna retunr actually the minimum value for \tau_beta
  min_tau_beta
  return(sh_ra$par)
}



# This function takes B and returns the parameters such as that I have a gamma distribution
#around the min. value for tau_beta
# and returns the parameters of a gamma distribution
return_min_tau_gamma <- function(B, d = 0, calculate_range = FALSE) {

  diff2 <- function(x, differences = 0) {
    if (differences == 0) {
      return(x)
    } else {
      return(diff(x, differences = differences))
    }
  }
  # First find the minimum value of tau_beta
  # I think this is equivalent to the problem:
  # Find max(var(beta))
  # subject to range(B beta) in (-0.5,0.5)
  # Need to create a minimisation function for this
  min_opt <- function(beta) {
    lambda <- 10
    neg_var <- -var(beta)
    yhat <- B %*% beta
    range_yhat <- range(yhat)
    # Create a quadratic penalty away from the range
    if (range_yhat[1] < -0.5) {
      neg_var <- neg_var + lambda * (range_yhat[1] + 0.5)^2
    }
    if (range_yhat[2] > 0.5) {
      neg_var <- neg_var + lambda * (range_yhat[2] - 0.5)^2
    }
    return(neg_var)
  }
  min_opt(runif(ncol(B)))

  ans <- optim(par = rnorm(ncol(B)),
               # Not sure what these two should be set at
               lower = -1, upper = 1,
               fn = min_opt,
               method = "L-BFGS-B",
               control = list(maxit = 20000))
  beta <- ans$par
  # plot(x, B %*% beta, type = 'l')

  # So the minimum value of tau_beta
  min_tau_beta <- 1 / var(diff2(beta, differences = d))

  # Finally use use min_tau_beta and max_tau_beta as 95% confidence
  # limits to obtain the parameters of a gamma distribution
  if(calculate_range){
    min_opt3 <- function(par) {
      shape <- par[1]
      rate <- par[2]
      alpha <- 0.95
      ucl <- qgamma(alpha + (1 - alpha)/2, shape = shape, rate = rate)
      lcl <- qgamma((1 - alpha)/2, shape = shape, rate = rate)
      return((ucl - min_tau_beta)^2 + (lcl - min_tau_beta)^2)
    }
    # min_opt3(runif(2))

    # Couldn't find a way of running B-FGS-B on this one
    # without it failing so gives some warnings
    sh_ra <- optim(par = c(1, 1),
                   fn = min_opt3)


    # This model gonna retunr actually the minimum value for \tau_beta
    if(identical(sh_ra$par,c(1,1))){
      return(NULL)
    }
  }

  return(list(gamma_parameters = c(1,1),
              min_tau_beta = min_tau_beta))
}

