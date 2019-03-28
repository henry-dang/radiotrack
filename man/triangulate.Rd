\name{triangulate}
\alias{triangulate}
\title{Estimate location of a signal from tracking data}
\usage{
triangulate(df, x, y, bearings, group, method,
    iterations, threshold)
}
\arguments{
  \item{df}{Dataframe}

  \item{x}{Name of the column containing the X-coordinate, in quotations.}

  \item{y}{Name of the column containing the Y-coordinate, in quotations.}

  \item{bearings}{Name of the column containing the bearings, in quotations.}

  \item{group}{Name of column containing a unique identifer for each group to perform
  calculations on, in quotations.}

  \item{method}{Method (mle, huber, andrews) used to perform calculations. Default is mle.}

  \item{iterations}{Maximum number of iterations to perform for each calculation before aborting. Default is 999.}

  \item{threshold}{Numeric specifying the precision of the calculations. Default is 0.0001.}
}
\description{
This function uses tracking data to compute the locaton of a radiotransmitter using either the maximum likelihood estimator, Huber, or Andrews method, as specified by the user.  The user will need to also specify a column containing the unique identifer that the function will use to group the tracking data and perform calculations on. A dataframe containing the X- and Y-Coordinate for the location of the radiotransmitter will be returned, along with the group identifier and method used. If the function fails to compute a result for a group, an error message will be displayed in the X- and Y-Coordinate columns of the returned datasframe.
}
\examples{
# results <- triangulate(df = "data", x = "X_Coordinate",
#     y = "Y_Coordinate", bearings = "Directions",
#     group = "Triangulation_Series")

# print(results)
}

\references{
Lenth, R.V. (1981). On Finding the Source of a Signal. Technometrics, 23(2), 149-154.
}
