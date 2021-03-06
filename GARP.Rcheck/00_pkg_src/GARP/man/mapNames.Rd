% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/aheatmap.R
\name{mapNames}
\alias{mapNames}
\title{mapping a vec to another name space}
\usage{
mapNames(vec, map, fac2string = TRUE)
}
\arguments{
\item{vec}{the input vector}

\item{map}{2-column data frame where first column is input names; 2nd column for output names}

\item{fac2string}{whether to treat factor as string}
}
\value{
a vector same length as vec
}
\description{
given a vector and a 2-column map, map names for vector based on the map.
This is useful if we want to update labels according to a table.
}
\details{
This function returns a new label vector with instructions in the map if a hit is found in the map
}
