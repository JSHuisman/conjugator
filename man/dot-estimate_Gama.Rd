% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/estimate_conj_rate.R
\name{.estimate_Gama}
\alias{.estimate_Gama}
\title{Gama estimate from data}
\usage{
.estimate_Gama(data)
}
\arguments{
\item{data}{Dataframe. Required columns vary depending on
the method used to estimate the conjugation rates. The following
columns are typically required: D.t, R.t, T.t are final population sizes;
D.0, R.0 are initial pop sizes;
psi.D, psi.R, psi.T are growth rates;
t is the time of measurement (assumes start at 0).
Optional columns: T.0, psi.max (for the Simonsen method).}
}
\value{
estimate vector
}
\description{
\code{.estimate_Gama} returns the Gama
estimate: log10(T/sqrt(DR)).
This is the log10 of the Dionisio estimate.
}
\details{
João Alves Gama, Rita Zilhão, and Francisco Dionisio.
Multiple plasmid interference - Pledging allegiance to my enemy’s enemy.
Plasmid, 93(August):17–23, 2017.
}
