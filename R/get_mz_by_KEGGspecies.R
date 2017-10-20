get_mz_by_KEGGspecies <-
function(kegg_species_code="hsa",outloc=NA,queryadductlist=c("M+H"),syssleep=0.01){
    
    
    data(adduct_table)
    adduct_table<-as.data.frame(adduct_table)
    
    
    path_list<-keggList(kegg_species_code,database="pathway")
    
    path_ids<-names(path_list)
    
    path_ids<-gsub(path_ids,pattern="path:",replacement="")
    
    map_res<-get_mz_by_KEGGpathwayIDs(kegg_pathway_ids=path_ids,outloc=outloc,queryadductlist=queryadductlist,syssleep=syssleep)
    
    if(length(map_res)>1){
        
    
        colnames(map_res)<-c("mz","chemical_ID","Name","Formula","MonoisotopicMass","Adduct","AdductMass")
    }
    
    return(map_res)
}
