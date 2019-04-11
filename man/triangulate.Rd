\name{triangulate}
\alias{triangulate}
\title{Estimate location of a signal from tracking data}
\usage{
triangulate(df, x, y, bearings, group, method,
    iterations, threshold)
}
\arguments{
  \item{df}{Dataframe}

  \item{x}{Name of the X-coordinate variable, as a string.}

  \item{y}{Name of the Y-coordinate variable, as a string.}

  \item{bearings}{Name of the bearings variable, as a string.}

  \item{group}{Name(s) of the unique identifer or combinations of identifers to split dataframe by, as strings.}

  \item{method}{Method (mle, huber, andrews) used to perform calculations. Default is mle.}

  \item{iterations}{Maximum number of iterations to perform for each set of calculations before aborting. Default is 999.}

  \item{threshold}{Numeric specifying the precision of the calculations. Default is 0.0001.}
}
\description{
This function uses tracking data to compute the locaton of a radiotransmitter using either the maximum likelihood estimator, Huber, or Andrews method, as specified by the user.  The user will need to also specify the variable(s) that the function will use to group the tracking data and perform calculations on.

A dataframe containing the X- and Y-Coordinates for the location of the radiotransmitter will be returned, along with the group identifier(s) and number of iterations. If the function fails to compute a result for a group, an error message will be displayed in the X- and Y-Coordinate columns of the returned dataframe.
}
\examples{
\dontrun{
tracking <- data.frame(
  "Date" = c(rep(1/1/2018, 3), rep(1/2/2018, 3),
  rep(1/1/2018, 3), rep(1/2/2018, 3)),
  "Triangulation.Series" = c(rep(1, 6), rep(2, 6)),
  "X.Coordinate" = c(-103.099, -103.098, -103.098, -103.080,
                     -103.080, -103.080, -99.874, -99.874,
                     -99.867, -99.873, -99.867, -99.874),
  "Y.Coordinate" = c(rep(33.774, 2), 33.775, 33.795, 33.796,
                     33.795, 36.506, 36.500, rep(36.5077, 3),
                     36.5042),
  "Bearings" = c(322, 296, 290, 69, 95, 51, 137, 51, 215,
                 97, 220, 50))

results <- triangulate(
  df = tracking,
  x = "X.Coordinate",
  y = "Y.Coordinate",
  bearings = "Bearings",
  groupp = .(Date, Triangulation.Series))
)
}
}

\references{
Lenth, R.V. (1981). On Finding the Source of a Signal. Technometrics, 23(2), 149-154.
}
