% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/aheatmap.R
\name{prepare4aheatmap}
\alias{prepare4aheatmap}
\title{given a mat and phenoDat and barColList, calculates sampleAnno and annoCols that can
 be directly passed to aheatmap}
\usage{
prepare4aheatmap(mat, phenoDat, barColList = NULL)
}
\arguments{
\item{mat}{the main mat for heatmap}

\item{phenoDat}{pheno data frame}

\item{barColList}{column bar color list. This should be a list of color vectors (the color pallete). By default, this is set as NULL,
which means getDefaultBarColList() will be called to set the default values.}
}
\value{
a list for aheatmap visualization
}
\description{
given a mat and phenoDat and barColList, calculates sampleAnno and annoCols that can
 be directly passed to aheatmap
}
