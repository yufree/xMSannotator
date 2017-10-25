#' get_allinf
#' 
#' This function allows users to map metabolite/chemical IDs to additional
#' fields such as KEGG BRITE IDs, KEGG Pathway IDs, SMPDB IDs, HMDB Status,
#' HMDB source, etc. depending upon the database.
#' 
#' 
#' @param dataA data matrix with the first column corresponding to KEGG or HMDB
#' IDs
#' @param dbname Database name. e.g.: "KEGG", "HMDB", "T3DB"
#' @return Returns an object with chemical/metabolite IDs merged with external
#' IDs, BRITE categories, pathway IDs, etc.
#' @author Karan Uppal
get_allinf <- function(dataA, dbname) {
    
    dataA <- as.data.frame(dataA)
    cnames <- colnames(dataA)
    
    cnames[1] <- "chemical_ID"
    colnames(dataA) <- as.character(cnames)
    
    
    # if(gregexpr(curated_res$chemical_ID[1],pattern='HMDB[0-9]*')[1]==1){
    
    if (dbname == "HMDB") {
        
        data(hmdbAllinf)
        curated_res <- merge(dataA, hmdbAllinf, by.x = "chemical_ID", 
            by.y = "HMDBID")
        
        curated_res <- curated_res[, -c(40:41)]
        rm(hmdbAllinf)
    } else {
        
        if (dbname == "KEGG") {
            
            data(keggotherinf)
            curated_res <- merge(dataA, keggotherinf, by.x = "chemical_ID", 
                by.y = "KEGGID")
            rm(keggotherinf)
        } else {
            if (dbname == "T3DB") {
                
                data(t3dbotherinf)
                curated_res <- merge(dataA, t3dbotherinf, 
                  by.x = "chemical_ID", by.y = "KEGGID")
                rm(t3dbotherinf)
            }
            
        }
    }
    
    return(curated_res)
}
