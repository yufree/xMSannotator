get_mz_by_KEGGpathwayIDs <-
function(kegg_pathway_ids,outloc=NA,queryadductlist=c("M+H"),syssleep=0.01){
    
    #if(is.na(adduct_table)==TRUE)
    {
        data(adduct_table)
        adduct_table<-as.data.frame(adduct_table)
    }
    #print(adduct_table[1:4,])
    kegg_comp_list<-{}
    map_res<-{}
    kegg_module_list<-{}
    for(kegg_pathway_id in kegg_pathway_ids){
        
        Sys.sleep(syssleep)
        k1<-keggGet(dbentries=kegg_pathway_id)
        
        kegg_comp_list<-c(kegg_comp_list,k1[[1]]$COMPOUND)
        
        if(length(kegg_comp_list)<1){
            
            
            kegg_module_list<-c(kegg_module_list,k1[[1]]$MODULE)
            
            
        }
    }
    
    kegg_module_ids<-names(kegg_module_list)
    
    for(kegg_pathway_id in kegg_module_ids){
        
        Sys.sleep(syssleep)
        k1<-keggGet(dbentries=kegg_pathway_id)
        
        kegg_comp_list<-c(kegg_comp_list,k1[[1]]$COMPOUND)
        
    }
    
    keggIDs<-names(kegg_comp_list)
    if(length(keggIDs)>0){
        map_res<-get_mz_by_KEGGcompoundIDs(keggIDs,queryadductlist,syssleep)
        #get_mz_by_KEGGcompoundIDs(keggIDs,queryadductlist,max.mz.diff=10,adduct_table)
    }
    map_res<-unique(map_res)
    map_res<-as.data.frame(map_res)
    if(length(map_res)>1){
        
        
        colnames(map_res)<-c("mz","chemical_ID","Name","Formula","MonoisotopicMass","Adduct","AdductMass")
    }
    
    return(map_res)
}
