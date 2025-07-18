% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tspost_O2.R
\name{tspost_O2}
\alias{tspost_O2}
\title{\code{tspost_O2} Calculate posterior probablility of two-station metabolism model.}
\usage{
tspost_O2(
  MET,
  tempup,
  tempdown,
  oxyup,
  oxydown,
  light,
  tt,
  z,
  osatup,
  osatdown,
  K600mean,
  K600sd,
  gas,
  n,
  lag
)
}
\arguments{
\item{MET}{Dataframe name of cleaned raw two station data (ex. "TS_S1S2")}

\item{tempup}{Temperature data from upstream station, deg C}

\item{tempdown}{Temperature data from downstream station, deg C}

\item{oxyup}{Oxygen data from upstream station, mg-O2/L}

\item{oxydown}{Oxygen data from downstream station, mg-O2/L}

\item{light}{any light unit}

\item{tt}{travel time, days}

\item{z}{Description}

\item{osatup}{Description}

\item{osatdown}{Description}

\item{Kmean}{Description}

\item{Ksd}{Description}
}
\value{
Populate here
}
\description{
This function calculates the posterior probability of the the two station metabolism model given
parameters. Internal function, called within twostationpostsum
}
\references{
Populate here
}
