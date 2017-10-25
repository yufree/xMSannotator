#' get_mz_by_KEGGspecies
#' 
#' This function generates list of expected m/z based on adducts for all
#' compounds associated with a species.
#' 
#' 
#' @param kegg_species_code KEGG organism code. e.g.: "hsa" for homo sapiens.
#' @param queryadductlist List of adducts to be used for searching.  eg:
#' c("M+H","M+Na","M+K"), c("positive") for all positive adducts, c("negative")
#' for all negative adducts, c("all") for all adducts
#' @param syssleep Wait time between queries to prevent overloading the KEGG
#' REST interface. e.g.: 0.1
#' @return Generates a text file with a list of expected m/z for all compounds
#' associated with a given species in KEGG.
#' @author Karan Uppal <kuppal2@@emory.edu>
#' @keywords ~KEGG species
get_mz_by_KEGGspecies <- function(kegg_species_code = "hsa", 
    outloc = NA, queryadductlist = c("M+H"), syssleep = 0.01) {
    
    
    data(adduct_table)
    adduct_table <- as.data.frame(adduct_table)
    
    
    path_list <- keggList(kegg_species_code, database = "pathway")
    
    path_ids <- names(path_list)
    
    path_ids <- gsub(path_ids, pattern = "path:", replacement = "")
    
    map_res <- get_mz_by_KEGGpathwayIDs(kegg_pathway_ids = path_ids, 
        outloc = outloc, queryadductlist = queryadductlist, 
        syssleep = syssleep)
    
    if (length(map_res) > 1) {
        
        
        colnames(map_res) <- c("mz", "chemical_ID", "Name", 
            "Formula", "MonoisotopicMass", "Adduct", "AdductMass")
    }
    
    return(map_res)
}
