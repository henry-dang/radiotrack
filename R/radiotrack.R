library("plyr")

triangulate <- function(df, x, y, bearings, group, method = mle,
                        iterations = 999, threshold = 0.0001){
  df <- plyr::ddply(.data=df, .variables=group, .fun=method,
              x, y, bearings, iterations = iterations, threshold = threshold)
}

# MLE method ------------------------------------------

mle <- function(df, x, y, bearings, group, iterations, threshold){

  y <- df[, y, drop = FALSE]
  x <- df[, x, drop = FALSE]
  bearings <- df[, bearings, drop = FALSE]

  if(max(bearings) - min(bearings) == 0){
   iterations <- 0
    counter <- iterations + 1
    x <- "Bearings cannot all be the same"
    y <- "Bearings cannot all be the same"
  }else{

  ### Using Eq. 2.6, obtain an initial value for x and y (coordinates)
  ### by ignoring the asteriks on s* and c* and assuming all d.i are equal.
  theta <- (pi / 180 * (90 - bearings))
  s.i <- sin(theta)
  c.i <- cos(theta)
  x.i <- x
  y.i <- y
  z.i <- s.i*x.i - c.i*y.i
  a <- array(c(sum(s.i^2), -sum(s.i*c.i), -sum(c.i*s.i),
               sum(c.i^2)), dim = c(2,2))
  b <- c(sum(s.i*z.i), -sum(c.i*z.i))
  xy <- solve(a, b)
  x.0 <- xy[1]
  y.0 <- xy[2]
  x <- as.numeric(0)
  y <- as.numeric(0)

  counter <- as.numeric(0)


  ### Use the initial x and y values to calculate interim values. Then iterate
  ### until x and y change by a negligible amount
  while(abs(abs(x)-abs(x.0)) > threshold & abs(abs(y)-abs(y.0)) > threshold &
        counter <= iterations){
    if(x != 0 & y != 0){
      x.0 <- x
      y.0 <- y
    }else{
      x <- x.0
      y <- y.0
    }
    ### Use the initial x and y values to calculate interim values for d_i, s*, and
    ### c*. Then iterate until x and y change by a negligible amount
    d.i <- ((x-x.i)^2 + (y-y.i)^2)^(1/2)
    s.star <- (y-y.i) / d.i^3
    c.star <- (x-x.i) / d.i^3
    a <- array(c(sum(s.i*s.star), -sum(s.i*c.star),
                 -sum(c.i*s.star), sum(c.i*c.star)), dim = c(2,2))
    b <- c(sum(s.star*z.i), -sum(c.star*z.i))
    xy <- solve(a, b)
    x <- xy[1]
    y <- xy[2]

    counter <- counter + 1

    if(is.nan(x) | is.nan(y) | counter == iterations){
      x <- "Failed"
      y <- "Failed"
      break()
    }
  }
  }
  results <- data.frame("X_Coordinate" = x,
                        "Y_Coordinate" = y,
                        "Iterations" = counter)
}


# Huber method ------------------------------------------

huber <- function(df, x, y, bearings, group, iterations, threshold){

  y <- df[, y, drop = FALSE]
  x <- df[, x, drop = FALSE]
  bearings <- df[, bearings, drop = FALSE]

  ### Using Eq. 4.7, obtain an initial value for x and y (coordinates)
  ### by ignoring the asteriks on s* and c* and assuming w.i = 1 and all d.i
  ### are equal.
  theta.i <- (pi / 180 * (90 - bearings))
  s.i <- sin(theta.i)
  c.i <- cos(theta.i)
  x.i <- x
  y.i <- y
  z.i <- s.i*x.i - c.i*y.i
  w.i <- as.numeric(1)
  a <- array(c(sum(w.i*s.i^2), -sum(w.i*s.i*c.i), -sum(w.i*c.i*s.i),
               sum(w.i*c.i^2)), dim = c(2,2))
  b <- c(sum(w.i*s.i*z.i), -sum(w.i*c.i*z.i))
  xy <- solve(a, b)
  x.0 <- xy[1]
  y.0 <- xy[2]
  x <- as.numeric(0)
  y <- as.numeric(0)

  counter <- as.numeric(0)

  ### Use the initial x and y values to calculate interim values. Then iterate
  ### until x and y change by a negligible amount
  while(abs(abs(x)-abs(x.0)) > threshold & abs(abs(y)-abs(y.0)) > threshold &
        counter <= iterations){
    if(x != 0 & y != 0){
      x.0 <- x
      y.0 <- y
    }else{
      x <- x.0
      y <- y.0
    }
    mu.i <- atan2((y - y.i), (x - x.i))
    C.w <- sum((w.i)*cos(theta.i - mu.i))/sum(w.i)

    ### Assuming that C.w is not too small, use Eq. 2.10 to approximate kappa.
    ### Then use kappa to calculate t
    kappa.inv <- 2 * (1 - C.w) + (1 - C.w)^2 *
      (0.48794 - 0.82905*C.w - 1.3915*C.w^2) / C.w
    kappa <- 1 / kappa.inv
    phi <- theta.i - mu.i
    ### Although not mentioned in the paper, there needs to be an absolute value
    ### sign in the square root.
    t <- abs((2 * kappa * (1 - cos(phi))))^(1/2)
    r <- phi %% (2*pi)   ### might not needed (see below)
    t <- ifelse(0 <= r & r <= pi, t, -t) ### might not be needed (see below)

    ### For the Huber method, calculate psi.H according to Eq. 4.5 and by assigning
    ### the tuning constant (c) a value of 1.5, as suggested by Lenth. Use psi.H
    ### to calculate interim values of weight (w.i).
    ### The sign of t is probably not important since 1) t is a square root
    ### function and should always be positive and 2) psi.H and t always have the
    ### same signs anyway, so w.i is always positive.
    psi.H <- sign(t) * pmin(abs(t), 1.5)
    w.i <- psi.H / t

    d.i <- ((x - x.i)^2 + (y - y.i)^2)^(1/2)
    s.star <- (y - y.i) / d.i^3
    c.star <- (x - x.i) / d.i^3
    a <- array(c(sum(w.i*s.i*s.star), -sum(w.i*s.i*c.star),
                 -sum(w.i*c.i*s.star), sum(w.i*c.i*c.star)), dim = c(2,2))
    b <- c(sum(w.i*s.star*z.i), -sum(w.i*c.star*z.i))
    xy <- solve(a, b)
    x <- xy[1]
    y <- xy[2]

    counter <- counter + 1

    if(is.nan(x)| is.nan(y) | counter == iterations){
      x <- "Failed"
      y <- "Failed"
      break()
    }
  }
  results <- data.frame("X_Coordinate"= x,
                        "Y_Coordinate" = y,
                        "Iterations" = counter)
}


# Andrews method ------------------------------------------

andrews <- function(df, x, y, bearings, group, iterations, threshold){

  y <- df[, y, drop = FALSE]
  x <- df[, x, drop = FALSE]
  bearings <- df[, bearings, drop = FALSE]

  ### Using Eq. 4.7, obtain an initial value for x and y (coordinates)
  ### by ignoring the asteriks on s* and c* and assuming w.i = 1 and all d.i
  ### are equal.
  theta.i <- (pi / 180 * (90 - bearings))
  s.i <- sin(theta.i)
  c.i <- cos(theta.i)
  x.i <- x
  y.i <- y
  z.i <- s.i*x.i - c.i*y.i
  w.i <- as.numeric(1)
  a <- array(c(sum(w.i*s.i^2), -sum(w.i*s.i*c.i), -sum(w.i*c.i*s.i),
               sum(w.i*c.i^2)), dim = c(2,2))
  b <- c(sum(w.i*s.i*z.i), -sum(w.i*c.i*z.i))
  xy <- solve(a, b)
  x.0 <- xy[1]
  y.0 <- xy[2]
  x <- as.numeric(0)
  y <- as.numeric(0)

  counter <- as.numeric(0)

  ### Use the initial x and y values to calculate interim values for d_i, s*, and
  ### c*. Then iterate until x and y change by a negligible amount
  while(abs(abs(x)-abs(x.0)) > threshold & abs(abs(y)-abs(y.0)) > threshold &
        counter <= iterations){
    if(x != 0 & y != 0){
      x.0 <- x
      y.0 <- y
    }else{
      x <- x.0
      y <- y.0
    }
    mu.i <- atan2((y - y.i), (x - x.i))
    C.w <- sum((w.i)*cos(theta.i - mu.i))/sum(w.i)

    ### Assuming that C.w is not too small, use Eq. 2.10 to approximate kappa.
    ### Then use kappa to calculate t
    kappa.inv <- 2*(1 - C.w) + (1 - C.w)^2 *
      (0.48794 - 0.82905*C.w - 1.3915*C.w^2) / C.w
    kappa <- 1 / kappa.inv
    phi <- theta.i - mu.i
    t <- abs((2 * kappa * (1 - cos(phi))))^(1/2)
    r <- phi %% (2*pi) ### might not be needed
    t <- ifelse(0 <= r & r <= pi, t, -t)  ### might not be needed

    ### For the Andrews method, calculate psi.A according to Eq. 4.6 (Robust
    ### Measures of Location for Directional Data, Lenth 1981) and by
    ### assigning  the tuning constant (c) a value of 1.5, as suggested by Lenth.
    ### Then use psi.A to calculate an interim values of weight (w.i).
    psi.A <- ifelse(abs(t) <= 1.5*pi, 1.5*sin(t/1.5), 0)
    w.i <- psi.A / t

    ### If all values of psi.A are 0, an error will result
    if(sum(psi.A == 0)){
      x <- "Failed"
      y <- "Failed"
      break()
    }

    d.i <- ((x - x.i)^2 + (y - y.i)^2)^(1/2)
    s.star <- (y - y.i) / d.i^3
    c.star <- (x - x.i) / d.i^3
    a <- array(c(sum(w.i*s.i*s.star), -sum(w.i*s.i*c.star),
                 -sum(w.i*c.i*s.star), sum(w.i*c.i*c.star)), dim = c(2,2))
    b <- c(sum(w.i*s.star*z.i), -sum(w.i*c.star*z.i))
    xy <- solve(a, b)
    x <- xy[1]
    y <- xy[2]

    counter <- counter + 1

    if(is.nan(x)| is.nan(y) | counter == iterations){
      x <- "Failed"
      y <- "Failed"
      break()
    }
  }
  results <- data.frame("X_Coordinate"= x,
                        "Y_Coordinate" = y,
                        "Iterations" = counter)
}


