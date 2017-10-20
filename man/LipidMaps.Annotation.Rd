\name{LipidMaps.Annotation}
\alias{LipidMaps.Annotation}
\title{LipidMaps.Annotation}
\description{
Annotation of features using LipidMaps.
}
\usage{
LipidMaps.Annotation(dataA, max.mz.diff = 10, num_nodes = 2,
queryadductlist = c("M+2H", "M+H+NH4", "M+ACN+2H",
"M+2ACN+2H", "M+H", "M+NH4", "M+Na", "M+ACN+H", 
"M+ACN+Na", "M+2ACN+H", "2M+H", "2M+Na", "2M+ACN+H"),
mode = "pos", outloc)
}
\arguments{
  \item{dataA}{
Feature table from apLCMS or XCMS. The first column should be m/z.
}
  \item{max.mz.diff}{
Mass tolerance in ppm for database matching. eg: 5
}
  \item{num_nodes}{
Number of cores to be used for parallel processing.
e.g.: 2
}
  \item{queryadductlist}{
List of adducts to be used for searching. 
eg: c("M+H","M+Na","M+K"), 
c("positive") for positive adducts,
c("negative") for negative adducts
c("all") for all adducts
}
  \item{mode}{
Ionization mode. e.g.: "pos" or "neg"
}
  \item{outloc}{
Output folder location
}
}
\details{}
\value{An object with annotation results}
\references{}
\author{
Karan Uppal <kuppal2@emory.edu>
}
\keyword{ ~LipidMaps }