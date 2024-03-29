% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/make_opf.R
\name{make_opf}
\alias{make_opf}
\title{Make OPF}
\usage{
make_opf(parameter, opf, AMAPStudio, overwrite = T, verbose = F, java = NULL)
}
\arguments{
\item{parameter}{The VPalm parameter file path and name (see \code{\link[=write_progeny]{write_progeny()}}
or \code{\link[=write_tree]{write_tree()}})}

\item{opf}{The target OPF file path and name}

\item{AMAPStudio}{The root path to AMAPStudio}

\item{overwrite}{Boolean. Should pre-existing OPF files overwriten ?}

\item{verbose}{Should the VPalm writting informations printed to the console ?}

\item{java}{Java path (optionnal, see details).}
}
\value{
Writes an OPF file, and return \code{TRUE} if the file was successfully written.
}
\description{
Use VPalm from AMAPStudio to make an Open Plant Format file using
a parameter file generally created by \code{\link[=write_progeny]{write_progeny()}} or \code{\link[=write_tree]{write_tree()}}.
}
\details{
The \code{java} argument can be a path to the java
executable if the user needs a particular version (for example if the default Java used is the Open JDK,
because ARCHIMED is only compatible with the Oracle version).
}
