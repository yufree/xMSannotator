\name{get_mz_by_KEGGcompoundIDs}
\alias{get_mz_by_KEGGcompoundIDs}
\title{get_mz_by_KEGGcompoundIDs}
\description{This function generates list of 
expected m/z based on adducts for KEGG compounds.}
\usage{
get_mz_by_KEGGcompoundIDs(keggIDs, queryadductlist = c("M+H"), 
syssleep = 0.01)
}
\arguments{
  \item{keggIDs}{Vector of KEGG compound IDs. 
  e.g.: c("C00073", "C00019")}
  \item{queryadductlist}{
  List of adducts to be used for searching. 
eg: c("M+H","M+Na","M+K"), 
c("positive") for all positive adducts,
c("negative") for all negative adducts,
c("all") for all adducts
}
  \item{syssleep}{
Wait time between queries to prevent overloading the 
KEGG REST interface. e.g.: 0.1
}

}
\value{
Returns an R object with a list of expected m/z
for input list of KEGG compounds.
}
\author{
Karan Uppal <kuppal2@emory.edu>
}
\keyword{ ~KEGG compounds}
\keyword{ ~get_mz_by_KEGGcompoundIDs}
