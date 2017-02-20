\name{get_mz_by_KEGGpathwayIDs}
\alias{get_mz_by_KEGGpathwayIDs}
\title{get_mz_by_KEGGpathwayIDs}
\description{
Generate list of expected m/z based on adducts for compounds in 
a given KEGG pathway.
}
\usage{
get_mz_by_KEGGpathwayIDs(kegg_pathway_ids, 
queryadductlist = c("M+H"), syssleep = 0.01)
}
\arguments{
  \item{kegg_pathway_ids}{
Vector of KEGG pathway IDs.  
e.g: c("map00270","map00966")
}
  \item{queryadductlist}{
List of adducts to be used for searching. 
eg: c("M+H","M+Na","M+K"), 

c("all") for all adducts
}
  \item{syssleep}{
Wait time between queries to prevent overloading the 
KEGG REST interface. e.g.: 0.1
}
}
\value{
Returns an R object with expected m/z
for compounds in the input list of KEGG pathways.
}
\author{
Karan Uppal <kuppal2@emory.edu>
}
\keyword{ ~KEGG pathways }
\keyword{ ~KEGG compounds }
