% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/compareCurves.r
\name{compareCurves}
\alias{compareCurves}
\title{Compare all growth curves to each other}
\usage{
compareCurves(
  mydata,
  filename = "test.csv",
  runs = 1000,
  raw = F,
  order = c(1, 3, 6, 9, 2, 5, 8, 4, 7, 10)
)
}
\arguments{
\item{mydata}{The output file from baTable or amaTable.}

\item{runs}{The number of permutations performed by compareGrowthCurves default is 1000.}

\item{raw}{If true, returns the raw letters and groups to be checked against the final table.}

\item{order}{The final reordering of your rows so that it goes from highest to lowest conc. By default set-up for basic antimicrbial growth assay.}
}
\value{
Returns a table that shows the letter codes of each of the corresponding groups. Writes the table to an excel file titled "letters_test.csv" by default.
}
\description{
Wrapper for statmod's compareGrowthCurves
}
