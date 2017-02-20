get_allinf <-
function(dataA,dbname){
    
    dataA<-as.data.frame(dataA)
    cnames<-colnames(dataA)
    
    cnames[1]<-"chemical_ID"
    colnames(dataA)<-as.character(cnames)
    
    
    #if(gregexpr(curated_res$chemical_ID[1],pattern="HMDB[0-9]*")[1]==1){
    
    if(dbname=="HMDB"){
        
        data(hmdbAllinf)
        curated_res<-merge(dataA,hmdbAllinf,by.x="chemical_ID",by.y="HMDBID")
        
        curated_res<-curated_res[,-c(40:41)]
        rm(hmdbAllinf)
    }else{
        
        if(dbname=="KEGG"){
            
            data(keggotherinf)
            curated_res<-merge(dataA,keggotherinf,by.x="chemical_ID",by.y="KEGGID")
            rm(keggotherinf)
        }else{
            if(dbname=="T3DB"){
                
                data(t3dbotherinf)
                curated_res<-merge(dataA,t3dbotherinf,by.x="chemical_ID",by.y="KEGGID")
                rm(t3dbotherinf)
            }
            
        }
    }

    return(curated_res)
}
