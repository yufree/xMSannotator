#' KEGG.Annotation
#' 
#' Annotation of features using KEGG.
#' 
#' 
#' @param dataA Feature table from apLCMS or XCMS. The first column should be
#' m/z.
#' @param max.mz.diff Mass tolerance in ppm for database matching. eg: 5
#' @param num_nodes Number of cores to be used for parallel processing. e.g.: 2
#' @param queryadductlist List of adducts to be used for searching.  eg:
#' c("M+H","M+Na","M+K"), c("positive") for positive adducts, c("negative") for
#' negative adducts c("all") for all adducts
#' @param mode Ionization mode. e.g.: "pos" or "neg"
#' @param outloc Output folder location
#' @return An object with annotation results
#' @author Karan Uppal <kuppal2@@emory.edu>
#' @keywords ~KEGG
KEGG.Annotation <- function(dataA, max.mz.diff = 10, num_nodes = 2, 
    queryadductlist = c("M+2H", "M+H+NH4", "M+ACN+2H", "M+2ACN+2H", 
        "M+H", "M+NH4", "M+Na", "M+ACN+H", "M+ACN+Na", "M+2ACN+H", 
        "2M+H", "2M+Na", "2M+ACN+H"), gradienttype = "Acetonitrile", 
    mode = "pos", outloc) {
    
    res <- simpleAnnotation(dataA, max.mz.diff = max.mz.diff, 
        num_nodes = num_nodes, queryadductlist = queryadductlist, 
        gradienttype = gradienttype, mode = mode, outloc, 
        db_name = "KEGG")
    
    return(res)
}
