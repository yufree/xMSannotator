get_mz_by_monoisotopicmass <-
function(monoisotopicmass,dbid=NA,name=NA,formula=NA,queryadductlist=c("M+H"),syssleep=0.01,adduct_table=NA){
    cnames<-c("mz", "ID","Name","Formula","MonoisotopicMass","Adduct","AdductMass")
    
    
    if(is.na(adduct_table)==TRUE){
        
        try(rm(adduct_table),silent=TRUE)
    }
   # if(is.na(adduct_table)==TRUE)
    {
        data(adduct_table)
        adduct_table<-as.data.frame(adduct_table)
    }
    
    #adduct_table<-read.table("/Users/karanuppal/Documents/Emory/JonesLab/Projects/xMSannotator/adduct_table.txt",sep="\t",header=TRUE)
    #adduct_table<-adduct_table[c(which(adduct_table[,6]=="S"),which(adduct_table[,6]=="Acetonitrile")),]
    
    adduct_names<-as.character(adduct_table[,1])
    adductlist<-adduct_table[,4]
    mult_charge<-adduct_table[,3]
    num_mol<-adduct_table[,2]
    names(adductlist)<-as.character(adduct_names)
    names(mult_charge)<-as.character(adduct_names)
    names(num_mol)<-as.character(adduct_names)
    alladducts<-adduct_names
    
    
    if(queryadductlist[1]=="positive")
    {
        queryadductlist<-adduct_names[which(adduct_table[,5]=="positive")]
        
    }else{
        if(queryadductlist[1]=="negative")
        {
            
            queryadductlist<-adduct_names[which(adduct_table[,5]=="negative")]
            
        }else{
            if(queryadductlist[1]=="all"){
                
                
                queryadductlist<-alladducts
                
                
            }else{
                if(length(which(queryadductlist%in%alladducts==FALSE))>0){
                    
                    errormsg<-paste("Adduct should be one of:",sep="")
                    for(i in alladducts){errormsg<-paste(errormsg, i,sep=" ; ")}
                    stop(errormsg, "\n\nUsage: feat.batch.annotation.KEGG(dataA,max.mz.diff=10, queryadductlist=c(\"M+H\", \"M+Na\"), xMSannotator.outloc, numnodes=1)",
                    "\n\n OR  feat.batch.annotation.KEGG(dataA,max.mz.diff=10, queryadductlist=c(\"positive\"), xMSannotator.outloc, numnodes=1)",
                    "\n\n OR  feat.batch.annotation.KEGG(dataA,max.mz.diff=10, queryadductlist=c(\"negative\"), xMSannotator.outloc, numnodes=1)",
                    "\n\n OR  feat.batch.annotation.KEGG(dataA,max.mz.diff=10, queryadductlist=c(\"all\"), xMSannotator.outloc, numnodes=1)"
                    )
                }
                
            }
        }
    }
				
                map_res<-{}
                
                if(is.na(dbid)==TRUE){
                    dbid="-"
                }
                
                if(is.na(name)==TRUE){
                    name="-"
                }
                if(is.na(formula)==TRUE){
                    formula="-"
                }
                
                exact_mass<-monoisotopicmass
                
                
                for(adnum in 1:length(queryadductlist))
                {
                    adductname=queryadductlist[adnum]
                    adductmass=adductlist[as.character(adductname)]
                    adductcharge=mult_charge[as.character(adductname)]
                    adductnmol=num_mol[as.character(adductname)]
                    
                    # mz=((adductnmol*exact_mass)+(adductmass))/(adductcharge))
                    
                    mz=((exact_mass*adductnmol)+(adductmass))/adductcharge
                    
                   # delta_ppm=(max.mz.diff)*(mz/1000000)
                  #  min_mz=round((mz-delta_ppm),5)
                   # max_mz=round((mz+delta_ppm),5)
                    res={} 
                    mzorig=round(exact_mass,5)
                   # delta_ppm=round(delta_ppm,5)
                    
                    syssleep1<-(syssleep/5)
                    Sys.sleep(syssleep1)
                    
                    
                    cur_map_res<-cbind(mz,dbid,name,formula,adductname,adductmass,mzorig)
                    cur_map_res<-as.data.frame(cbind(mz,as.character(dbid),as.character(name),as.character(formula),mzorig,adductname,adductmass))
                    #print(cur_map_res)
                    cur_map_res<-as.data.frame(cur_map_res)
                    map_res<-rbind(map_res,cur_map_res)
                    
                }
                
                
                colnames(map_res)<-cnames
                
                map_res<-unique(map_res)
                map_res<-as.data.frame(map_res)
                
                return(map_res)
}
