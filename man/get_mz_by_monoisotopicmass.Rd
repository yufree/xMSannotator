\name{get_mz_by_monoisotopicmass}
\alias{get_mz_by_monoisotopicmass}
\title{get_mz_by_monoisotopicmass}
\description{
Generate list of expected m/z for a specific monoisotopic mass
}
\usage{
get_mz_by_monoisotopicmass(monoisotopicmass, dbid = NA,
name = NA, formula = NA, queryadductlist = c("M+H"), 
syssleep = 0.01)
}
\arguments{
  \item{monoisotopicmass}{
Monoisotopic mass. e.g.: 149.051
}
  \item{dbid}{
Database or user-defined ID. e.g.: "M001"
}
  \item{name}{
Metabolite name. e.g.: "Methionine"
}
  \item{formula}{
Chemical formula. e.g.: "C5H11NO2S"
}
  \item{queryadductlist}{
List of adducts to be used for searching. 
eg: c("M+H","M+Na","M+K"), 
c("positive") for positive adducts,
c("negative") for negative adducts
c("all") for all adducts
}
  \item{syssleep}{
Wait time between queries to prevent overloading the 
KEGG REST interface. e.g.: 0.1
}
}
\value{
Returns an R object with a list of expected m/z
for the input monoisotopic mass.
}
\author{
Karan Uppal
}