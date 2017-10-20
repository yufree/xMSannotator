\name{multilevelannotation}
\alias{multilevelannotation}
\title{multilevelannotation}
\description{
The function uses a multi-level scoring algorithm to annotate features using HMDB, KEGG,
T3DB, and LipidMaps.
}
\usage{
multilevelannotation(dataA, max.mz.diff = 10, max.rt.diff = 10, 
cormethod = "pearson", 
num_nodes = 2, queryadductlist = c("all"), mode = "pos",
outloc, db_name = "HMDB", adduct_weights = NA, num_sets = 3000,
allsteps = TRUE, corthresh = 0.7, NOPS_check = TRUE, customIDs = NA,
missing.value = NA, deepsplit = 2, 
networktype = "unsigned",
minclustsize = 10, module.merge.dissimilarity = 0.2, filter.by = c("M+H"),
redundancy_check = TRUE,
min_ions_perchem = 1, biofluid.location = NA, origin = NA,
status = NA, boostIDs = NA, max_isp = 5,
MplusH.abundance.ratio.check = FALSE, 
customDB = NA, HMDBselect = "union", mass_defect_window = 0.01,
mass_defect_mode = "pos",  pathwaycheckmode = "pm")
}
\arguments{
  \item{dataA}{
Peak intensity matrix. The first two columns should be "mz" and "time"
followed by sample intensities.
}
  \item{max.mz.diff}{
Mass tolerance in ppm for database matching. e.g.: 10
}
  \item{max.rt.diff}{
Retention time (s) tolerance for finding peaks derived
from the same parent metabolite. e.g.: 10
}
  \item{cormethod}{
Method for correlation. e.g.: "pearson" or "spearman".
The "pearson" implementation is computationally faster.
}
  \item{num_nodes}{
Number of computing cores to be used for parallel processing.
e.g.:2
}
  \item{queryadductlist}{
Adduct list to be used for database matching. 
e.g. c("all") for all possible positive or negative adducts
or c("M+H","M+Na","M+ACN+H") for specifying subset of adducts. 
Run data(adduct_table) for list of all adducts.
}
  \item{mode}{
Ionization mode. e.g.:"pos" or "neg"
}
  \item{num_sets}{
How many subsets should the total number of metabolites
in a database should be divided into to faciliate
parallel processing or prevent memory overload?
e.g.: 1000
}
  \item{outloc}{
Output folder location. 
e.g.: "C:\\Documents\\ProjectX\\"
}
  \item{db_name}{
Database to be used for annotation.
e.g.: "HMDB", "KEGG", "T3DB", "LipidMaps"
}
  \item{adduct_weights}{
Adduct weight matrix. Run
data(adduct_weights) to see an example
adduct weight matrix.
}
  \item{allsteps}{
If FALSE, only step 1 that involves module and
retention time based clustering is performed. e.g.: TRUE
}
  \item{corthresh}{
Minimum correlation threshold between peaks to qualify
as adducts/isotopes of the same metabolite.
}
  \item{NOPS_check}{
Should elemental ratio checks be performed as outlined
in Fiehn 2007? e.g. TRUE
}
  \item{customIDs}{
Custom list of select database IDs (HMDB, KEGG, etc.) to be used
for annotation. This should be a data frame. 
Run data(customIDs) to see an example.
}
  \item{missing.value}{
How are missing values represented in the input
peak intensity matrix? e.g.: NA
}
  \item{deepsplit}{
How finely should the clusters be split? e.g.: 2
Please see WGCNA for reference. 
}
  \item{networktype}{
Please see WGCNA for reference: e.g: 
"unsigned" or "signed". 
}
  \item{minclustsize}{
Minimum cluster size. e.g: 10
}
  \item{module.merge.dissimilarity}{
Maximum dissimilarity measure (i.e., 1-correlation)
to be used for merging modules in WGCNA. 
e.g.:0.2 
}
  \item{filter.by}{
Require the presence of certain adducts for a 
high confidence match. e.g.: c("M+H")
}
  \item{redundancy_check}{
Should stage 5 that involves redundancy based
filtering be performed? e.g.: TRUE or FALSE
}
  \item{min_ions_perchem}{
Minimum number of adducts/isotopes to be 
present for a match to be considered high
confidence. e.g.:2
}
  \item{biofluid.location}{
Used only for HMDB. e.g.: NA, "Blood" ,"Urine", "Saliva"
Set to NA to turn off this option.
Please see http://www.hmdb.ca/metabolites or 
run data(hmdbAllinf); head(hmdbAllinf)
for more details.
}
  \item{origin}{
Used only for HMDB. e.g.: NA, "Endogenous", 
"Exogenous", etc. Set to NA to turn off this
option.
Please see http://www.hmdb.ca/metabolites 
or run data(hmdbAllinf); head(hmdbAllinf)
for more details.
}
  \item{status}{
Used for HMDB. e.g.: NA, "Detected", 
"Detected and Quantified", "Detected and Not Quantified",
"Expected and Not Quantified". Set to NA to turn off
this option. Please see
http://www.hmdb.ca/metabolites or run 
data(hmdbAllinf) for more details.
}
  \item{boostIDs}{
Databased IDs of previously validated metabolites.
e.g.: c("HMDB00696"). Set to NA to turn off
this option.
}
  \item{max_isp}{
Maximum number of expected isotopes. e.g.: 5
}
  \item{MplusH.abundance.ratio.check}{
Should MplusH be the most abundant adduct? e.g. TRUE
or FALSE
}
  \item{customDB}{
Custom database. Run: data(custom_db);
head(custom_db) to see more details on
formatting. Set to NA to turn off
this option
}
  \item{HMDBselect}{
How to select metabolites based on HMDB origin, 
biolfuid, and status filters? e.g.: "all" to take union,
"intersect" to take intersection
}
  \item{mass_defect_window}{
Mass defect window in daltons for finding isotopes.
e.g.: 0.01
}
  \item{mass_defect_mode}{
"pos" for finding positive isotopes; "neg" for
finding unexpected losses/fragments; "both" for
finding isotopes and unexpected losses/fragments
}
  \item{pathwaycheckmode}{
How to perform pathway based evaluation? 
"pm": boosts the scores if there are other metabolites
from the same pathway that are also assigned to the
same module.
"p": boosts the scores if there are other metabolites
from the same pathway without accounting for module
membership.
}
}
\details{
Multistage clustering algorithm based on intensity profiles, 
retention time characteristics, mass defect, isotope/adduct
patterns and correlation with signals for metabolic
precursors and products. The algorithm uses high-resolution
mass spectrometry data for a series of samples with
common properties and publicly available chemical,
metabolic and environmental databases to assign
confidence levels to annotation results.  
}
\value{
The function generates output at each stage:
Stage 1 includes modules and retention time based clustering
of features without any annotation
Stage 2 includes  modules and retention time based
clustering of features along with simple m/z based
database matching
Stage 3 includes scores for annotations assigned in stage 2
Stages 4 and 5 include the confidence levels before and 
after redundancy (multiple matches) filtering, respectively
}
\references{
Langfelder P, Horvath S. WGCNA: an R package for weighted
correlation network analysis. BMC Bioinformatics. 2008 Dec 29;9:559.
Wishart DS et al. HMDB 3.0--The Human Metabolome Database in 2013. 
Nucleic Acids Res. 2013 Jan;41(Database issue):D801-7. 
Kanehisa M. The KEGG database. Novartis Found Symp. 2002;
247:91-101;discussion 101-3, 119-28, 244-52. Review.
Lim E, et al.  T3DB: a comprehensively annotated database of 
common toxins and their targets. Nucleic Acids Res.
2010 Jan;38(Database issue):D781-6.
Sleno L. The use of mass defect in modern mass spectrometry.
J Mass Spectrom.  2012 Feb;47(2):226-36.
Zhang H, Zhang D, Ray K, Zhu M. Mass defect filter technique 
and its applications to drug metabolite identification by high-resolution
mass spectrometry. J Mass Spectrom. 2009 Jul;44(7):999-1016. 
}
\author{
Karan Uppal <kuppal2@emory.edu>
}
\keyword{ ~xMSannotator }
\keyword{ ~multilevelannotation }
\keyword{ ~annotation }
