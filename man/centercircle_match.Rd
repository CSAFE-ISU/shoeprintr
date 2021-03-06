% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{centercircle_match}
\alias{centercircle_match}
\title{Matches a center circle input and a reference shoeprint}
\usage{
centercircle_match(input, reference, output)
}
\arguments{
\item{input}{The input print}

\item{reference}{The reference print}

\item{output}{The output obtained from the match_print function}
}
\description{
Function to perform matching of center footprints by matching all the candidate circles with input circles. This function performs initial matching between 3 circles on the input print and 27 candidate circles on the reference shoeprint. After the initial cliques, it select 2 best circles for each input circle and performs reinforcement matching on them. In reinforcement matching, full points on reference circle are considered instead of user defined. The function returns the final circle parameters and the statistics of the triangle formed.
}
