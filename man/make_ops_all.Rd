% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/make_ops.R
\name{make_ops_all}
\alias{make_ops_all}
\title{Format and write all OPS}
\usage{
make_ops_all(Progeny, design, map, path, overwrite, average = FALSE, ...)
}
\arguments{
\item{Progeny}{The progenies}

\item{design}{The design of experiment, generally computed using \code{\link[=design_plot]{design_plot()}}}

\item{map}{The age of the plantation in month after planting}

\item{path}{The path of the target folder}

\item{overwrite}{Boolean. Should pre-existing OPS files overwriten ?}

\item{average}{Boolean. Use average tree instead of the sampled ones.}

\item{...}{Further arguments to pass to \code{\link[=format_ops]{format_ops()}}}
}
\value{
An OPS for each progeny, each design and each map.
}
\description{
Format and write OPS from a set of design experiment(s),
map and progenies, by applying \code{\link[=format_ops]{format_ops()}} and \code{\link[=write_ops]{write_ops()}} sequentially
}
\details{
The function uses \code{\link[base:mapply]{base::mapply()}} to apply both \code{\link[=format_ops]{format_ops()}} and
\code{\link[=write_ops]{write_ops()}} to any number of progeny, design, or map. So if these arguments
are provided with equal length, they will be applied in parallel (\emph{i.e.} in a
multivariate mode) as \code{\link[base:mapply]{base::mapply()}} do. Each can also have length one. in this
case it will be re-cycled for each combination.
}
