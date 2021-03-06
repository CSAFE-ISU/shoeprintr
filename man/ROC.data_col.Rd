% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{ROC.data_col}
\alias{ROC.data_col}
\title{Calculate sensitivity and specificity at different threshold for drawing ROC curves}
\usage{
ROC.data_col(pred_prob, nsame, ndiff, thres, Method)
}
\arguments{
\item{nsame}{number of mates}

\item{ndiff}{number of non-mates}

\item{thres}{Threshold for ROC curve}

\item{Method}{method name}
}
\description{
Using prediction by some methods, calculate sensitivity and specificity at different threshold for drawing ROC curves.
}
