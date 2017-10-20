simpleAnnotation <-
function(dataA,max.mz.diff=10,num_nodes=2,queryadductlist=c("M+2H","M+H+NH4","M+ACN+2H","M+2ACN+2H","M+H","M+NH4","M+Na","M+ACN+H","M+ACN+Na","M+2ACN+H","2M+H","2M+Na","2M+ACN+H"),
gradienttype="Acetonitrile",mode="pos",outloc,db_name="KEGG")
{
    
    #db_name="KEGG"
    if(db_name=="KEGG"){
        
        
        data(keggCompMZ)
        chemCompMZ<-keggCompMZ
        suppressWarnings(rm(keggCompMZ))
        
    }else{
        if(db_name=="HMDB"){
            data(hmdbCompMZ)
            chemCompMZ<-hmdbCompMZ
            suppressWarnings(rm(hmdbCompMZ))
            
        }else{
            if(db_name=="T3DB"){
                data(t3dbCompMZ)
                chemCompMZ<-t3dbCompMZ
                suppressWarnings(rm(t3dbCompMZ))
                
            }else{
                
                if(db_name=="LipidMaps"){
                    data(lipidmapsCompMZ)
                    chemCompMZ<-lipidmapsCompMZ
                    suppressWarnings((lipidmapsCompMZ))
                }
            }
            
        }
    }
    NOPS_check=TRUE
    
    data(adduct_table)
    
    
    adduct_table<-as.data.frame(adduct_table)
    
    allowWGCNAThreads(nThreads=num_nodes)
    #rm(adduct_table)
    #data(adduct_table)
    
    if(FALSE){
    adduct_table$Adduct<-gsub(adduct_table$Adduct,pattern="2M\\+2H",replacement="M+H")
    adduct_table$Adduct<-gsub(adduct_table$Adduct,pattern="3M\\+3H",replacement="M+H")
    
    adduct_table$Adduct<-gsub(adduct_table$Adduct,pattern="2M\\+2Na",replacement="M+Na")
    adduct_table$Adduct<-gsub(adduct_table$Adduct,pattern="3M\\+3Na",replacement="M+Na")
    
    adduct_table$Adduct<-gsub(adduct_table$Adduct,pattern="2M\\+2K",replacement="M+K")
    adduct_table$Adduct<-gsub(adduct_table$Adduct,pattern="3M\\+3K",replacement="M+K")
    }
    
    
    adduct_table<-unique(adduct_table)
    
    suppressWarnings(dir.create(outloc))
    
    suppressWarnings(
    if(queryadductlist=="all" & mode=="pos"){
        
        adduct_names<-adduct_table$Adduct[(adduct_table$Type=="S" & adduct_table$Mode=="positive") | (adduct_table$Type==gradienttype & adduct_table$Mode=="positive")]
        
        adduct_table<-adduct_table[which(adduct_table$Adduct%in%adduct_names),]
        
    }else{
        if(queryadductlist=="all" & mode=="neg"){
            
            adduct_names<-adduct_table$Adduct[(adduct_table$Type=="S" & adduct_table$Mode=="negative") | (adduct_table$Type==gradienttype & adduct_table$Mode=="negative")]
            adduct_table<-adduct_table[which(adduct_table$Adduct%in%adduct_names),]
        }else{
            
            adduct_names<-adduct_table$Adduct[which(adduct_table$Adduct%in%queryadductlist)]
            
            adduct_table<-adduct_table[which(adduct_table$Adduct%in%adduct_names),]
            
            if(length(adduct_names)<1){
                
                stop("Invalid adducts selected.")
            }
        }
    }
    )
    adduct_names<-unique(adduct_names)
    
    
    
    chemCompMZ<-chemCompMZ[which(chemCompMZ$Adduct%in%adduct_names),]
    
    print("Dimension of database matrix:")
    print(dim(chemCompMZ))
    
    
    cl<-makeSOCKcluster(num_nodes)
    
    clusterEvalQ(cl, library(XML))
    clusterEvalQ(cl, library(R2HTML))
    clusterEvalQ(cl, library(RCurl))
    clusterEvalQ(cl, library(SSOAP))
    clusterEvalQ(cl, library(limma))
    
    clusterEvalQ(cl, library(plyr))
    
    clusterEvalQ(cl, "processWSDL")
    clusterEvalQ(cl, library(png))
    clusterExport(cl, "Annotationbychemical_IDschild")
    
    clusterExport(cl, "find.Overlapping.mzs")
    
    clusterExport(cl, "find.Overlapping.mzsvparallel")
     clusterExport(cl, "overlapmzchild")
	   clusterExport(cl, "getVenn")
       clusterExport(cl, "adduct_table")
       
       s1<-seq(1,length(adduct_names))
       print("Mapping m/z to metabolites:")
#save(list=ls(),file="matching_data.Rda")

dataA<-dataA[,c(1:2)]

adduct_names<-as.character(adduct_names)
#print(adduct_names)



	if(length(adduct_names)>1){
       l2<-parLapply(cl,s1,Annotationbychemical_IDschild,dataA=dataA, queryadductlist=c(adduct_names),adduct_type=c("S",gradienttype),
       max.mz.diff=max.mz.diff,outloc=outloc,keggCompMZ=chemCompMZ,otherdbs=FALSE,otherinfo=FALSE,adduct_table=adduct_table,num_nodes=num_nodes)
       
       stopCluster(cl)
       
	
       
       levelB_res<-{};
       for(j in 1:length(l2)){
           if(length(l2[[j]])>1){
               levelB_res<-rbind(levelB_res,l2[[j]])
           }
       }
       
       rm(l2)
       
	}else{

		levelB_res<-Annotationbychemical_IDschild(adduct_index=1,dataA=dataA, queryadductlist=c(adduct_names),adduct_type=c("S",gradienttype),
       max.mz.diff=max.mz.diff,outloc=outloc,keggCompMZ=chemCompMZ,otherdbs=FALSE,otherinfo=FALSE,adduct_table=adduct_table,num_nodes=num_nodes)
	}
    
       
       if(length(levelB_res)<2){
           
           print("No matches found.")
           try(stopCluster(cl),silent=TRUE)
           try(close(cl),silent=TRUE)
           return(-1)
       }
	
       levelB_res$mz<-as.numeric(as.character(levelB_res$mz))
       
       levelB_res$time<-as.numeric(as.character(levelB_res$time))
       
       levelB_res<-as.data.frame(levelB_res)
       
       
       uniq_formula<-as.character(unique(levelB_res$Formula))
       
       bad_formula<-which(is.na(uniq_formula)==TRUE)
       if(length(bad_formula)>0){
           uniq_formula<-uniq_formula[-c(bad_formula)]
       }
       
       cl<-makeSOCKcluster(num_nodes)
       
       
       clusterExport(cl, "check_golden_rules")
       clusterExport(cl, "check_element")
       #clusterExport(cl, "uniq_formula")
       #clusterExport(cl, "NOPS_check")
       
       levelB_res_check<-parLapply(cl,1:length(uniq_formula),function(j,uniq_formula,NOPS_check){
           
           curformula<-as.character(uniq_formula[j])
           return(check_golden_rules(curformula,NOPS_check=NOPS_check))
           
       },uniq_formula=uniq_formula,NOPS_check=NOPS_check)
       suppressWarnings(stopCluster(cl))
       
       #save(levelB_res_check,file="xMSannotator_levelB_check.Rda")
       levelB_res_check2<-ldply(levelB_res_check,rbind)
       
       levelB_res_check3<-levelB_res_check2[which(levelB_res_check2[,2]==1),]
       
       
       levelB_res<-levelB_res[which(levelB_res$Formula%in%as.character(levelB_res_check3[,1])),]
       
       water_adducts<-c("M+H-H2O","M+H-2H2O","M-H2O-H")
       
       water_adduct_ind<-which(levelB_res$Adduct%in%water_adducts)
       
       cl<-makeSOCKcluster(num_nodes)
       
       
       clusterExport(cl, "check_element")
       
       
       
       if(length(water_adduct_ind)>0){
           levelB_res2<-levelB_res[c(water_adduct_ind),]
           
           levelB_res<-levelB_res[-c(water_adduct_ind),]
           
           sind1<-seq(1:dim(levelB_res2)[1])
           
           levelB_res_check3<-parLapply(cl,sind1,function(j){
               
               adduct<-as.character(levelB_res2$Adduct[j])
               curformula<-as.character(levelB_res2$Formula[j])
               
               numoxygens<-check_element(curformula,"O")
               
               if(numoxygens>0){
                   bool_check<-1
               }else{
                   bool_check<-0
               }
               
               res<-cbind(curformula,bool_check)
               res<-as.data.frame(res)
               return(res)
               
               
           })
           
           levelB_res_check4<-ldply(levelB_res_check3,rbind)
           
           valid_form<-{}
           
           if(length(which(levelB_res_check4[,2]==1))>0){
               levelB_res_check4<-levelB_res_check4[which(levelB_res_check4[,2]==1),]
               
               
               valid_form<-which(levelB_res2$Formula%in%as.character(levelB_res_check4[,1]))
           }
           if(length(valid_form)>0){
               levelB_res2<-levelB_res2[valid_form,]
               levelB_res<-rbind(levelB_res,levelB_res2)
           }
           
       }
       multiresmat<-levelB_res
       rm(levelB_res)
       dupmz<-multiresmat$mz[which(duplicated(multiresmat$mz)==TRUE)]
       
       
       MatchCategory<-rep("Multiple",dim(multiresmat)[1])
       
       MatchCategory[-which(multiresmat$mz%in%dupmz)]<-"Unique"
       
       levelB_res<-cbind(MatchCategory,multiresmat)
       
       rownames(levelB_res)<-NULL
       #save(levelB_res,file="xMSannotator_levelB.Rda")
       return(levelB_res)
}
