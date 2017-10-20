\name{ChemSpider.annotation}
\alias{ChemSpider.annotation}
\title{ChemSpider.annotation}
\description{
The ChemSpider annotation utility maps the input mz to known compounds in
biological, chemical, and environmental databases. Users have the option to
specify which data sources in ChemSpider should be used for annotaion. 
NOTE: Please obtain a security token from www.chemspider.com/Register.aspx
before using this function. You will have to specify your role as a 
Service Subscriber.
}
\usage{
ChemSpider.annotation(dataA, max.mz.diff = 10, queryadductlist = c("M+H"), 
xMSannotator.outloc, numnodes = 1, 
datasources = c("KEGG"), tokenstr = NA,
maxhits = 30,syssleep=5)
}
\arguments{
 \item{dataA}{Feature table from apLCMS or XCMS. 
 The first column should be m/z.}
  \item{max.mz.diff}{Metlin matching m/z threshold in ppm. eg: 5}
  \item{queryadductlist}{List of adducts to be used for searching. 
  eg: c("M+H","M+Na","M+K"),
  c("positive") for all positive ion adducts, or c("negative") for all negative ion adducts, or 
  c("all") for all adducts as defined in Metlin.}
  \item{xMSannotator.outloc}{Output location where the 
  HTML and text reports will be written. 
  eg: "C:/experiment1/xMSanalyzeroutput/"}
  \item{numnodes}{Number of computing nodes to use. Default: 1}
  \item{datasources}{Provide a list of data sources to search against. 
The following data sources are supported: KEGG, Human Metabolome Database, 
LipidMAPS, MassBank, PubChem, BioCyc, SMPDB Small Molecule Pathway Database,
EPA DSSTox, EPA Toxcast, Pesticide Common Names, ChEMBL, ChEBI, and NIST. 
The names of the data sources should be provided exactly as listed above.
eg: c("KEGG", "PubChem", "EPA Toxcast", "Human Metabolome Database", "BioCyc")
}
  \item{tokenstr}{
  Please register on ChemSpider website to obtain a security token. 
  www.chemspider.com/Register.aspx
  Provide that security token here. Eg: "mysecuritycode"
}
  \item{maxhits}{
  Maximum number of result hits to retrieve. Default: 30
}
 \item{syssleep}{
  Wait time between queries. Default: 5
}
}
\details{
This utility uses the ChemSpider web services to search for compounds based on the 
input mz, mass search threshold (+/- ppm), and adduct. 
The output is generated as an HTML report and a text file that includes
pathway and compound annotations with links to several biological, chemical, and 
environmental databases. The function takes as input a data frame with a list of 
input m/z, a user-defined m/z threshold (ppm)
to define the minimum and maximum mass range, list of adducts; 
eg: c(“M+H”, “M+H-H2O”), the output folder location, list of data sources to be
searched, and a security token string. 
ChemSpider_annotation_results.html
}
\value{A list is returned.  \item{html.res }{Annotation report in HTML format}
\item{text.res }{Text delimited annotation report}
}
\author{Karan Uppal <kuppal2@emory.edu>}
\keyword{ ~annotation }
\keyword{ ~pesticides }
\keyword{ ~pathway }
\keyword{ ~ChemSpider }