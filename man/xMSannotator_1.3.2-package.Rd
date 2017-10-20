\name{xMSannotator_1.3.2-package}
\alias{xMSannotator_1.3.2-package}
\alias{xMSannotator}
\docType{package}
\title{xMSannotator}
\description{This package includes functions to perform metabolite annotation using HMDB, KEGG, LipidMaps, T3DB,
and selected data sources in ChemSpider. 

The multilevelannotation() function uses a multi-criteria approach based on data-driven modules, retention time, 
adducts, isotopes, mass defect, and various chemical rules and ratio checks for assigning annotations into different
confidence levels: high (3), medium (2), low (1), and none (0). The function is designed to work with KEGG, HMDB, 
T3DB, LipidMaps or a custom database provided by the user. You can use status, origin,
and biofluid location arguments to select subset of metabolites in HMDB for annotation. 
The "customIDs" argument can be used to only get confidence levels for a subset of metabolites.

KEGG.annotation, T3DB. annotation, HMDB.annotation, LipidMaps.annotation, and ChemSpider.annotation functions perform simple mz based
annotaion based on the user defined adducts and mass search threshold.

Additional utilites are provided for generating list of expected mz for compounds based on KEGG pathway IDs, species IDs, and compound IDs.
Users can generate KEGG pathway map images and color code specific metabolites/compounds.
}
\details{
\tabular{ll}{
Package: \tab xMSannotator\cr
Type: \tab Package\cr
Version: \tab 1.3.2\cr
Date: \tab 2017-06-29\cr
License: \tab gpl2.0\cr
LazyLoad: \tab yes\cr
}
}
\author{
Karan Uppal
Maintainer: <kuppal2@emory.edu>
}
\references{
http://pubs.acs.org/doi/abs/10.1021/acs.analchem.6b01214
}
\keyword{ xMSannotator }