\name{get_kegg_map}
\alias{get_kegg_map}
\title{get_kegg_map}
\description{
This function allows users to color individiul objects in KEGG pathway maps
and generate the output as an PNG image.
}
\usage{
get_kegg_map(keggids, pathwayid, bgcolor, fgcolor)
}
\arguments{
  \item{keggids}{
KEGG IDs to be colored. e.g.: c("C05345","C05378")
}
  \item{pathwayid}{
KEGG pathway ID: e.g.: "map01230"
}
  \item{bgcolor}{
Vector of the size of the keggids argument to indicate
background color for each object. e.g.: c("red", "green")
}
  \item{fgcolor}{
Vector of the size of the keggids argument to indicate
border and text color for each object. e.g.: c("black", "blue")
}
}
\value{
PNG image for the KEGG pathway with colored objects
}
\author{
Karan Uppal
}
\keyword{ ~pathway }

