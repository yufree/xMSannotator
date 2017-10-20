\name{get_brite_res}
\alias{get_brite_res}
\title{get_brite_res}
\description{
Get KEGG BRITE and pathway information for specific
KEGG compound IDs.
}
\usage{
get_brite_res(keggid)
}
\arguments{
  \item{keggid}{
Input KEGG compound ID. e.g.: "C00073"
}
}
\value{
Return an R object with KEGG ID mapped to KEGG
Brite category, pathway IDs, and links to external 
databases.
}
\author{
Karan Uppal
}
\keyword{ ~KEGG BRITE }