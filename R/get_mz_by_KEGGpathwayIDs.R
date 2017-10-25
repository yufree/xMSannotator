#' get_mz_by_KEGGpathwayIDs
#' 
#' Generate list of expected m/z based on adducts for compounds in a given KEGG
#' pathway.
#' 
#' 
#' @param kegg_pathway_ids Vector of KEGG pathway IDs.  e.g:
#' c("map00270","map00966")
#' @param queryadductlist List of adducts to be used for searching.  eg:
#' c("M+H","M+Na","M+K"),
#' 
#' c("all") for all adducts
#' @param syssleep Wait time between queries to prevent overloading the KEGG
#' REST interface. e.g.: 0.1
#' @return Returns an R object with expected m/z for compounds in the input
#' list of KEGG pathways.
#' @author Karan Uppal <kuppal2@@emory.edu>
#' @keywords ~KEGG pathways ~KEGG compounds
get_mz_by_KEGGpathwayIDs <- function(kegg_pathway_ids, outloc = NA, 
    queryadductlist = c("M+H"), syssleep = 0.01) {
    
    # if(is.na(adduct_table)==TRUE)
    {
        data(adduct_table)
        adduct_table <- as.data.frame(adduct_table)
    }
    # print(adduct_table[1:4,])
    kegg_comp_list <- {
    }
    map_res <- {
    }
    kegg_module_list <- {
    }
    for (kegg_pathway_id in kegg_pathway_ids) {
        
        Sys.sleep(syssleep)
        k1 <- keggGet(dbentries = kegg_pathway_id)
        
        kegg_comp_list <- c(kegg_comp_list, k1[[1]]$COMPOUND)
        
        if (length(kegg_comp_list) < 1) {
            
            
            kegg_module_list <- c(kegg_module_list, k1[[1]]$MODULE)
            
            
        }
    }
    
    kegg_module_ids <- names(kegg_module_list)
    
    for (kegg_pathway_id in kegg_module_ids) {
        
        Sys.sleep(syssleep)
        k1 <- keggGet(dbentries = kegg_pathway_id)
        
        kegg_comp_list <- c(kegg_comp_list, k1[[1]]$COMPOUND)
        
    }
    
    keggIDs <- names(kegg_comp_list)
    if (length(keggIDs) > 0) {
        map_res <- get_mz_by_KEGGcompoundIDs(keggIDs, queryadductlist, 
            syssleep)
        # get_mz_by_KEGGcompoundIDs(keggIDs,queryadductlist,max.mz.diff=10,adduct_table)
    }
    map_res <- unique(map_res)
    map_res <- as.data.frame(map_res)
    if (length(map_res) > 1) {
        
        
        colnames(map_res) <- c("mz", "chemical_ID", "Name", 
            "Formula", "MonoisotopicMass", "Adduct", "AdductMass")
    }
    
    return(map_res)
}
