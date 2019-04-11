library("plyr")


###### make triangulate function

triangulate <- function(df, x, y, bearings, group, method = mle,
                        iterations = 999, threshold = 0.0001){
  df <- plyr::ddply(.data=df, .variables=group, .fun=method,
                    x, y, bearings, iterations = iterations, threshold = threshold)
}

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





############## test triangulate function

tracking <- read.csv("tracking.csv")

# error message prints
tracking1 <- tracking[(tracking$Triangulation.Series==1178),]
errormes <- triangulate(tracking1, "X_Coordinate", "Y_Coordinate", "Direction",
                        .(Triangulation.Series, Bird.ID))

# error message doesn't print
no.errormes <- triangulate(tracking, "X_Coordinate", "Y_Coordinate", "Direction",
                       .(Triangulation.Series, Bird.ID))


