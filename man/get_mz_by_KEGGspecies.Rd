\name{get_mz_by_KEGGspecies}
\alias{get_mz_by_KEGGspecies}
\title{get_mz_by_KEGGspecies}
\description{
This function generates list of 
expected m/z based on adducts for all compounds associated
with a species.
}
\usage{
get_mz_by_KEGGspecies(kegg_species_code = "hsa",
queryadductlist = c("M+H"),syssleep=0.01)
}
\arguments{
  \item{kegg_species_code}{
KEGG organism code. e.g.: "hsa"
for homo sapiens.
}
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
Generates a text file with a list of expected m/z
for all compounds associated with a given species
in KEGG.
}
\author{
Karan Uppal <kuppal2@emory.edu>
}
\keyword{ ~KEGG species }
