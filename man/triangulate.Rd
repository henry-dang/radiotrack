\name{triangulate}
\alias{triangulate}
\title{Estimate location of a point using tracking data}
\usage{
hello()
}
\arguments{
  \item{df}{Dataframe}

  \item{x}{Name of the column containing the x-coordinate, in quotations.}

  \item{y}{Name of the column containing the y-coordinate, in quotations.}

  \item{bearings}{Name of the column containing the bearings, in quotations.}

  \item{group}{Name of column containing a unique identifer for each group to perform
  calculations on, in quotations.}

  \item{method}{Method (mle, huber, andrews) used to perform calculations. Default is mle.}

  \item{iterations}{Maximum number of iterations to perform for each calculation before aborting. Default is 999.}

  \item{threshold}{Numeric specifying the precision of the calculations. Default is 0.0001.}
}
\description{
Prints 'Hello, world!'.
}
\examples{
triangulate(df = "data", x = "X_Coordinate", y = "Y_Coordinate", bearings = "Directions",
group = "Triangulation_Series")
}
