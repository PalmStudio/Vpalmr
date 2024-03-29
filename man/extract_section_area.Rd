% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/extract_raw.R
\name{extract_section_area}
\alias{extract_section_area}
\title{Extract section}
\usage{
extract_section_area(df_archi, Row, Col)
}
\arguments{
\item{df_archi}{The architectural data.frame}

\item{Row}{The row index of the section}

\item{Col}{The column index of the section}
}
\value{
A data.frame with the section form data
}
\description{
Extract a section from the form
}
\details{
The row and column of the section is the one from the upper left corner
}
\examples{
\dontrun{
extract_section_area(5,1)
}
}
