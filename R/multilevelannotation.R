multilevelannotation <-
function(dataA,max.mz.diff=10,max.rt.diff=10,cormethod="pearson",num_nodes=2,queryadductlist=c("all"),
gradienttype="Acetonitrile",mode="pos",outloc,db_name="HMDB", adduct_weights=NA,num_sets=3000,allsteps=TRUE,
corthresh=0.7,NOPS_check=TRUE,customIDs=NA,missing.value=NA,deepsplit=2,networktype="unsigned",minclustsize=10,
module.merge.dissimilarity=0.2,filter.by=c("M+H"),redundancy_check=TRUE,min_ions_perchem=1,biofluid.location=NA,origin=NA,
status=NA,boostIDs=NA,max_isp=5,MplusH.abundance.ratio.check=FALSE,customDB=NA,HMDBselect="union",
mass_defect_window=0.01,mass_defect_mode="pos",dbAllinf=NA,pathwaycheckmode="pm")
{

options(warn=-1)

allowWGCNAThreads(nThreads=num_nodes)


dataA<-as.data.frame(dataA)
dataA$mz<-round(dataA$mz,5)
dataA$time<-round(dataA$time,1)

dataA[,-c(1:2)]<-round(dataA[,-c(1:2)],1)

if(is.na(module.merge.dissimilarity)==TRUE){
    
    module.merge.dissimilarity=1-corthresh
}



if(is.na(customIDs)==FALSE){
    
    customIDs=as.data.frame(customIDs)
}



data(adduct_table)


print(paste("Annotating using ",db_name," database:",sep=""))

max_diff_rt<-max.rt.diff
cutheight=1-corthresh #module.merge.dissimilarity
time_step=1
step1log2scale=FALSE

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
setwd(outloc)

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
        
    }
}

adduct_names<-unique(adduct_names)

res_list<-new("list")

db_name_list=db_name

outloc_allres<-outloc

#for(db_index in 1:length(db_name_list))
{
    
    #   db_name=db_name_list[db_index]
    
    #outloc<-paste(outloc_allres,"/",db_name,"results",sep="")
 

l1<-list.files(outloc_allres)

check_step1<-which(l1=="step1_results.Rda")


if(length(check_step1)<1){

print("Status 1: Computing modules using WGCNA")

check_levelA<-which(l1=="xMSannotator_levelA_modules.Rda")

if(is.na(missing.value)==FALSE){
dataA<-replace(as.matrix(dataA),which(dataA==missing.value),NA)
dataA<-as.data.frame(dataA)
}

dataA<-unique(dataA)


mzid<-paste(dataA$mz,dataA$time,sep="_")

if(length(which(duplicated(mzid)==TRUE))>0){
	dataA<-dataA[-which(duplicated(mzid)==TRUE),]
}

#dataA<-dataA[order(dataA$mz,dataA$time),]

mzid<-paste(dataA$mz,dataA$time,sep="_")

system.time(global_cor<-WGCNA::cor(t(dataA[,-c(1:2)]),nThreads=num_nodes,method=cormethod,use = 'p'))

global_cor<-round(global_cor,2)

dataA<-unique(dataA)
setwd(outloc_allres)

mzid<-paste(dataA$mz,dataA$time,sep="_")
colnames(global_cor)<-mzid
rownames(global_cor)<-mzid

save(global_cor,file="global_cor.Rda")

mycl_metabs=NA
hr = flashClust(as.dist(1-global_cor),method="complete");
  dissTOMCormat<-(1-global_cor)
 
a1<-apply(global_cor,1,function(x){length(which(x>corthresh))})

  fname_c<-paste("NumberOfCorrelationsPerFeature_cor",corthresh,".csv",sep="")

  write.csv(a1,file=fname_c)
 
  rm(global_cor)

  #using dynamic
  mycl_metabs <-cutreeDynamic(hr,distM= dissTOMCormat,deepSplit=1, minClusterSize=minclustsize, pamRespectsDendro = FALSE, pamStage=TRUE,verbose=0)
   

    
   #mycl_metabs <-cutree(hr, h=max(hr$height)/2)
   #save(mycl_metabs,file="mycl_metabs.Rda")

clustmethod="WGCNA"
if(length(check_levelA)<1){

	if(clustmethod=="WGCNA"){
		

	levelA_res<-get_peak_blocks_modulesvhclust(dataA=dataA,simmat=NA,adjacencyfromsimilarity=FALSE,time_step=time_step,max.rt.diff=max_diff_rt,outloc,column.rm.index=NA,cor.thresh=NA,deepsplit=deepsplit,
minclustsize=minclustsize,cutheight=cutheight,cormethod=cormethod,networktype = networktype,num_nodes=num_nodes,step1log2scale=step1log2scale,mycl_metabs=mycl_metabs)


	setwd(outloc_allres)
	levelA_res<-levelA_res[,-c(1:4)]
	save(levelA_res,file="xMSannotator_levelA_modules.Rda")
	


	}else{
		
		if(clustmethod=="graph")
		{
			
			levelA_res<-get_peak_blocks_graph(dataA,simmat=global_cor,adjacencyfromsimilarity=TRUE,time_step=3,max.rt.diff=max_diff_rt,outloc,column.rm.index=NA,cor.thresh=NA,cormethod=cormethod,networktype = networktype,num_nodes=num_nodes)
						
		}else{


			if(clustmethod=="hclust"){
			levelA_res<-get_peak_blocks_hclust(dataA,time_step=time_step,max.rt.diff=max_diff_rt,outloc,column.rm.index=NA,cor.thresh=NA,deepsplit=deepsplit,minclustsize=minclustsize,cutheight=cutheight,cormethod=cormethod)
			}

		}
	
	}

	setwd(outloc_allres)
	
    
	write.csv(levelA_res,file="Stage1.csv",row.names=FALSE)
	#write.table(levelA_res,file="Stage1.txt",sep="\t",row.names=FALSE)

}else{
	setwd(outloc_allres)
	load("xMSannotator_levelA_modules.Rda")
	
	#alldegrees<-levelA_res[,c(1:4)]
	#levelA_res<-levelA_res[,-c(1:4)]

	
	
}

setwd(outloc)

levelA_res<-levelA_res[order(levelA_res$mz,levelA_res$time),]

levelA_res<-levelA_res[,c(1:3)]	

dataA<-dataA[order(dataA$mz,dataA$time),]

mean_int_vec<-apply(dataA[,-c(1:2)],1,function(x){mean(x,na.rm=TRUE)})
	


dataA<-dataA[,c(1:2)]


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
		
	}
}

adduct_names<-unique(adduct_names)



if(db_name=="HMDB"){
    


data(hmdbAllinf)

hmdbAllinf$Name<-gsub(hmdbAllinf$Name,pattern="[\\\"\']",replacement="")


hmdbAllinfv3.6=hmdbAllinf


rm(hmdbAllinf)
rm(hmdbAllinf, envir = .GlobalEnv )




suppressWarnings(if(is.na(customIDs)==TRUE){
    
    customIDs<-hmdbAllinfv3.6[,c(1,20)]
    

    
    
if(is.na(biofluid.location)==FALSE & is.na(origin)==TRUE){


	gres<-gregexpr(hmdbAllinfv3.6$BioFluidLocation,pattern=biofluid.location,ignore.case=TRUE)
	
	if(is.na(status)==FALSE){
			
			gres3<-gregexpr(hmdbAllinfv3.6$HMDBStatus,pattern=status,ignore.case=TRUE)

			gres3_good<-which(gres3==1)
			gres_good<-which(gres==1)
			if(HMDBselect=="intersect"){
			gres<-intersect(gres_good,gres3_good)
			}else{

			gres<-c(gres_good,gres3_good)
			gres<-unique(gres)

			}
			

		}else{

			gres<-which(gres==1)
			
		}
	

	customIDs<-hmdbAllinfv3.6[gres,c(1,20)]
	
	rm(hmdbAllinfv3.6)


}else{


		if(is.na(biofluid.location)==FALSE & is.na(origin)==FALSE){


        gres<-gregexpr(hmdbAllinfv3.6$BioFluidLocation,pattern=biofluid.location,ignore.case=TRUE)

	gres2<-gregexpr(hmdbAllinfv3.6$Origin,pattern=origin,ignore.case=TRUE)

	gres_good<-which(gres==1)
	gres2_good<-which(gres2==1)
	

	
	
		if(is.na(status)==FALSE){
			
			gres3<-gregexpr(hmdbAllinfv3.6$HMDBStatus,pattern=status,ignore.case=TRUE)

			gres3_good<-which(gres3==1)


				if(HMDBselect=="intersect"){
			gres<-intersect(gres_good,gres2_good,gres3_good)
			}else{

			gres<-c(gres_good,gres2_good,gres3_good)
			gres<-unique(gres)

			}
		}else{


			if(HMDBselect=="intersect"){
			gres<-intersect(gres_good,gres2_good)
			}else{

			gres<-c(gres_good,gres2_good)
			gres<-unique(gres)

			}
		}



       

                customIDs<-hmdbAllinfv3.6[gres,c(1,20)]
		

	
	}else{


				if(is.na(biofluid.location)==TRUE & is.na(origin)==FALSE){


        gres<-gregexpr(hmdbAllinfv3.6$Origin,pattern=origin,ignore.case=TRUE)


		if(is.na(status)==FALSE){
			
			gres3<-gregexpr(hmdbAllinfv3.6$HMDBStatus,pattern=status,ignore.case=TRUE)

			gres3_good<-which(gres3==1)
			gres_good<-which(gres==1)
			
			if(HMDBselect=="intersect"){
			gres<-intersect(gres_good,gres3_good)
			}else{

			gres<-c(gres_good,gres3_good)
			gres<-unique(gres)

			}
		}else{

			gres<-which(gres==1)
			
		}
	
	
	
	  customIDs<-hmdbAllinfv3.6[gres,c(1,20)]

		}else{


				if(is.na(status)==FALSE){
			
            
			gres3<-gregexpr(hmdbAllinfv3.6$HMDBStatus,pattern=status,ignore.case=TRUE)

			gres3_good<-which(gres3==1)
			
			gres<-gres3_good
			
		
	
		
	customIDs<-hmdbAllinfv3.6[gres,c(1,20)]
	  
                }else{
                    
                    
                    customIDs<-hmdbAllinfv3.6[,c(1,20)]
                }
			

		}


	}
}

})

rm(hmdbAllinfv3.6,envir = .GlobalEnv)
data(hmdbCompMZ)

hmdbCompMZ$mz<-round(as.numeric(as.character(hmdbCompMZ$mz)),5)

hmdbCompMZ$Name<-gsub(hmdbCompMZ$Name,pattern="[\\\"\']",replacement="")

#hmdbbad<-c("HMDB29244","HMDB29245","HMDB29246")

#if(length(which(hmdbCompMZ$HMDBID%in%hmdbbad))>0){
#	hmdbCompMZ<-hmdbCompMZ[-which(hmdbCompMZ$HMDBID%in%hmdbbad),]

#}

suppressWarnings(if(is.na(customIDs)==FALSE){
    
    customIDs<-unique(customIDs)
    
    
    hmdbCompMZ<-hmdbCompMZ[which(hmdbCompMZ$HMDBID%in%customIDs[,1]),]
    
    hmdbCompMZ<-hmdbCompMZ[which(hmdbCompMZ$Adduct%in%adduct_names),]
})

hmdbCompMZ<-hmdbCompMZ[which(hmdbCompMZ$Adduct%in%adduct_names),]

chemCompMZ<-hmdbCompMZ


print("Dimension of precomputed HMDB m/z database")
print(dim(chemCompMZ))

hmdbCompMZ<-1
hmdbAllinf<-1
try(rm(hmdbCompMZ),silent=TRUE)

try(rm(hmdbCompMZ, envir = .GlobalEnv ),silent=TRUE)
try(rm(hmdbAllinf, envir = .GlobalEnv ),silent=TRUE)

}else{

	if(db_name=="KEGG"){
	
    
    data(keggCompMZ)

    keggCompMZ$mz<-round(as.numeric(as.character(keggCompMZ$mz)),5)

    keggCompMZ$Name<-gsub(keggCompMZ$Name,pattern="[\\\"\']",replacement="")   


	
	suppressWarnings(if(is.na(customIDs)==FALSE){

		customIDs<-unique(customIDs)	

	
		keggCompMZ<-keggCompMZ[which(keggCompMZ$KEGGID%in%customIDs[,1]),]
	
		keggCompMZ<-keggCompMZ[which(keggCompMZ$Adduct%in%adduct_names),]
	})

	keggCompMZ<-keggCompMZ[which(keggCompMZ$Adduct%in%adduct_names),]
	chemCompMZ<-keggCompMZ
	
	print("Dimension of precomputed KEGG m/z database")
	print(dim(chemCompMZ))

	rm(keggCompMZ)
	
	rm(keggCompMZ, envir = .GlobalEnv )

	}else{
	
		if(db_name=="LipidMaps"){
		

                data(lipidmapsCompMZ)
				lipidmapsCompMZ<-lipidmapsCompMZ[which(lipidmapsCompMZ$Adduct%in%adduct_names),]
				
			lipidmapsCompMZ$mz<-round(as.numeric(as.character(lipidmapsCompMZ$mz)),5)

			lipidmapsCompMZ$Name<-gsub(lipidmapsCompMZ$Name,pattern="[\\\"\']",replacement="")

				chemCompMZ<-lipidmapsCompMZ

				print("Dimension of precomputed LipidMaps m/z database")
				print(dim(chemCompMZ))

				rm(lipidmapsCompMZ)
				rm(lipidmapsCompMZ, envir = .GlobalEnv )

	
		}else{
					if(db_name=="T3DB"){
		
                   

					data(t3dbCompMZ)
					t3dbCompMZ<-t3dbCompMZ[which(t3dbCompMZ$Adduct%in%adduct_names),]

					t3dbCompMZ$mz<-round(as.numeric(as.character(t3dbCompMZ$mz)),5)

t3dbCompMZ$Name<-gsub(t3dbCompMZ$Name,pattern="[\\\"\']",replacement="")
					chemCompMZ<-t3dbCompMZ

					print("Dimension of precomputed T3DB m/z database")
					print(dim(chemCompMZ))

					rm(t3dbCompMZ)
					rm(t3dbCompMZ, envir = .GlobalEnv )

	
				}else{
				
						if(db_name=="Custom"){
						
                            data(adduct_table)
                            
                       
                            
							inputmassmat<-customDB

							if(length(which(duplicated(inputmassmat[,1])==TRUE))>0){
							inputmassmat<-inputmassmat[-which(duplicated(inputmassmat[,1])==TRUE),]
                            
							}
                            
							mz_search_list<-{}
							mz_search_list<-lapply(1:dim(inputmassmat)[1],function(m){
                                
					#print(head(adduct_table))
					adduct_names<-as.character(adduct_names)
					#print(adduct_names)

					mz_search_list<-get_mz_by_monoisotopicmass(monoisotopicmass=as.numeric(as.character(inputmassmat[m,4])),dbid=inputmassmat[m,1],name=as.character(inputmassmat[m,2]),formula=as.character(inputmassmat[m,3]),queryadductlist = adduct_names,adduct_table=adduct_table)
	
					return(mz_search_list)
										})

							mz_search_list<-ldply(mz_search_list,rbind)
							save(mz_search_list,file="mz_search_list.Rda")
							chemCompMZ<-mz_search_list
							rm(mz_search_list)
							rm(customDB)
							rm(customDB, envir = .GlobalEnv )

							

                        }else{
                            
                            stop("db_name should be: KEGG, HMDB, T3DB, LipidMaps, or Custom")
                            
                        }
						
				
				}
			
		
		}
	}
}
data(adduct_table)

if(is.na(adduct_weights)==TRUE){
data(adduct_weights)
adduct_weights1<-matrix(nrow=2,ncol=2,0)
                adduct_weights1[1,]<-c("M+H",1)
                adduct_weights1[2,]<-c("M-H",1)
                adduct_weights<-as.data.frame(adduct_weights1)
                colnames(adduct_weights)<-c("Adduct","Weight")
adduct_weights<-as.data.frame(adduct_weights)

}

    chemCompMZ$Name<-gsub("([-:();])|[[:punct:]]","\\1",chemCompMZ$Name)
	chemCompMZ$Formula<-gsub("([-:();])|[[:punct:]]","\\1",chemCompMZ$Formula)
	
	
cnames<-colnames(chemCompMZ)

cnames[2]<-"chemical_ID"
colnames(chemCompMZ)<-cnames


	

randsetsize<-1000
if(dim(dataA)[1]<1000){

	randsetsize<-dim(dataA)[1]
}



chemCompMZ<-as.data.frame(chemCompMZ)

setwd(outloc)
l1<-list.files(outloc)
check_levelB<-which(l1=="xMSannotator_levelB.Rda")

save(chemCompMZ,file="chemCompMZ.Rda")
#save(list=ls(),file="matching_data.Rda")

if(length(check_levelB)<1){

cl<-makeSOCKcluster(num_nodes)
		
		clusterEvalQ(cl, library(XML))
		clusterEvalQ(cl, library(R2HTML))
		clusterEvalQ(cl, library(RCurl))
		clusterEvalQ(cl, library(SSOAP))
 clusterEvalQ(cl, library(limma))	
 
 clusterEvalQ(cl, library(plyr))
 	
			clusterEvalQ(cl, "processWSDL")
		clusterEvalQ(cl, library(png))
		clusterExport(cl, "Annotationbychemical_IDschild_multilevel")
        
        clusterExport(cl, "Annotationbychemical_IDschild")

		clusterExport(cl, "find.Overlapping.mzs")
        clusterExport(cl, "find.Overlapping.mzsvparallel")

     clusterExport(cl, "overlapmzchild")
	   clusterExport(cl, "getVenn")

if(length(which(duplicated(chemCompMZ$Formula)==TRUE))>0){

if(db_name=="Custom"){
#chemCompMZ_unique_formulas<-chemCompMZ

#chemCompMZ_unique_formulas<-as.data.frame(chemCompMZ_unique_formulas)

#chemCompMZ_unique_formulas$mz<-as.numeric(as.character(chemCompMZ_unique_formulas$mz))


chemCompMZ$mz<-as.numeric(as.character(chemCompMZ$mz))


}
chemCompMZ_unique_formulas<-chemCompMZ[-which(duplicated(chemCompMZ$Formula)==TRUE),]
chemCompMZ_unique_formulas<-rbind(chemCompMZ_unique_formulas,chemCompMZ[which(chemCompMZ$chemical_ID%in%chemCompMZ_unique_formulas$chemical_ID),])
chemCompMZ_unique_formulas<-unique(chemCompMZ_unique_formulas)
chemCompMZ_unique_formulas<- chemCompMZ_unique_formulas[order(chemCompMZ_unique_formulas$chemical_ID),]


}else{

chemCompMZ_unique_formulas<-chemCompMZ
}


#save(chemCompMZ_unique_formulas,file="chemCompMZ_unique_formulas.Rda")
save(chemCompMZ,file="chemCompMZ.Rda")
rm(chemCompMZ)

#formula_ID<-paste("Formula",seq(1:dim(chemCompMZ_unique_formulas)[1]),sep="_")
#chemCompMZ_unique_formulas$chemical_ID<-formula_ID



formula_table<-table(chemCompMZ_unique_formulas$Formula)
uniq_formulas<-names(formula_table)
formula_ID<-paste("Formula",seq(1:length(uniq_formulas)),sep="_")

formula_id_mat<-cbind(formula_ID,uniq_formulas)
formula_id_mat<-as.data.frame(formula_id_mat)
colnames(formula_id_mat)<-c("Formula_ID","Formula")

chemCompMZ_unique_formulas<-merge(chemCompMZ_unique_formulas,formula_id_mat,by="Formula")

chemCompMZ_unique_formulas$chemical_ID<-chemCompMZ_unique_formulas$Formula_ID
chemCompMZ_unique_formulas<-chemCompMZ_unique_formulas[,c("mz","chemical_ID","Name","Formula","MonoisotopicMass","Adduct","AdductMass")]

chemCompMZ_unique_formulas<-as.data.frame(chemCompMZ_unique_formulas)

chemIDs<-unique(chemCompMZ_unique_formulas$chemical_ID)  #unique(chemCompMZ$Formula) #


#save(chemCompMZ_unique_formulas,file="chemCompMZ_unique_formulasB.Rda")




#s1<-seq(1,length(chemIDs))
s1<-seq(1,length(adduct_names))
print("Status 2: Mapping m/z to metabolites:")

#write.csv(chemCompMZ,file="chemCompMZ.csv")

adduct_names<-as.character(adduct_names)


#annotation_multilevel_
#save(list=ls(),file="test.Rda")
#system.time(l2<-parLapply(cl,s1,Annotationbychemical_IDschild_multilevel,dataA=dataA, queryadductlist=c(adduct_names),adduct_type=c("S",gradienttype),
#adduct_table=adduct_table,max.mz.diff=max.mz.diff,outloc=outloc,keggCompMZ=chemCompMZ_unique_formulas,otherdbs=FALSE,otherinfo=FALSE,chemIDs=chemIDs,num_nodes=num_nodes))

#a1<-Annotationbychemical_IDschild_multilevel(1,dataA=dataA, queryadductlist=c(adduct_names),adduct_type=c("S",gradienttype),
#adduct_table=adduct_table,max.mz.diff=max.mz.diff,outloc=outloc,keggCompMZ=chemCompMZ_unique_formulas,otherdbs=FALSE,otherinfo=FALSE,chemIDs=chemIDs,num_nodes=num_nodes)

#a1<-Annotationbychemical_IDschild(1,dataA=dataA, queryadductlist=c(adduct_names),adduct_type=c("S",gradienttype),
#adduct_table=adduct_table,max.mz.diff=max.mz.diff,outloc=outloc,keggCompMZ=chemCompMZ_unique_formulas,otherdbs=FALSE,otherinfo=FALSE,num_nodes=num_nodes)


l2<-parLapply(cl,s1,Annotationbychemical_IDschild,dataA=dataA, queryadductlist=c(adduct_names),adduct_type=c("S",gradienttype),
adduct_table=adduct_table,max.mz.diff=max.mz.diff,outloc=outloc,keggCompMZ=chemCompMZ_unique_formulas,otherdbs=FALSE,otherinfo=FALSE,num_nodes=num_nodes)

stopCluster(cl)

#print(format(object.size(l2),unit="auto"))
#save(l2,file="l2.Rda")

rm(chemCompMZ)
	levelB_res<-{};
	for(j in 1:length(l2)){
		if(length(l2[[j]])>1){	
		levelB_res<-rbind(levelB_res,l2[[j]])
		}
	}


	rm(l2)

	if(nrow(levelB_res)<1){
		stop("No matches found.")
	}

	MatchCategory=rep("Multiple",dim(levelB_res)[1])
	
	
	
	levelB_res$mz<-as.numeric(as.character(levelB_res$mz))

	levelB_res$time<-as.numeric(as.character(levelB_res$time))

	levelB_res<-as.data.frame(levelB_res)
	
	levelB_res<-cbind(levelB_res,MatchCategory)

    

	print("DB matches")
	levelB_res<-unique(levelB_res)
	print(dim(levelB_res))

	

	uniq_formula<-as.character(unique(levelB_res$Formula))

	bad_formula<-which(is.na(uniq_formula)==TRUE)
	if(length(bad_formula)>0){
	uniq_formula<-uniq_formula[-c(bad_formula)]
	}
	
	cl<-makeSOCKcluster(num_nodes)
	

        clusterExport(cl, "check_golden_rules")
	clusterExport(cl, "check_element")

	

	levelB_res_check<-parLapply(cl,1:length(uniq_formula),function(j,uniq_formula,NOPS_check){
		
		curformula<-as.character(uniq_formula[j])
		return(check_golden_rules(curformula,NOPS_check=NOPS_check))
		
		},uniq_formula=uniq_formula,NOPS_check=NOPS_check)
	stopCluster(cl)
	

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
	#levelB_res<-levelB_res[,-c(1)]
	
	
	colnames(levelB_res)<-c("theoretical.mz","chemical_ID","Name","Formula","MonoisotopicMass","Adduct","AdductMass","mz","time","MatchCategory")
	
	save(levelB_res,file="xMSannotator_levelB.Rda")


	
}else{

	print("Status 2: Using existing m/z mapping results:")
	load("xMSannotator_levelB.Rda")
}

try(rm(hmdbAllinf, envir = .GlobalEnv ),silent=TRUE)

try(rm(hmdbAllinfv3.6),silent=TRUE)
try(rm(hmdbCompMZ),silent=TRUE)



levelA_res$mz<-round(levelA_res$mz,5)
levelB_res$mz<-round(levelB_res$mz,5)

#print("here")
#no_match<-dataA[-which(dataA$mz%in%levelB_res$mz),1:2]

#no_match_annot<-matrix(nrow=dim(no_match)[1],ncol=8,"-")

#no_match_1<-cbind(no_match[,1],no_match_annot,no_match[,-c(1)])

levelB_res <-levelB_res #[,1:10]

#colnames(no_match_1)<-colnames(levelB_res )

#levelB_res <-rbind(levelB_res ,no_match_1)

last_col_ind<-dim(levelB_res )[2]

levelB_res <-levelB_res #[,c(1,10,2:8)]


levels_A<-levels(levelA_res$Module_RTclust)

levels_A<-table(levelA_res$Module_RTclust)
levels_A<-names(levels_A)
#dataA_orig<-dataA

levelA_res<-levelA_res[order(levelA_res$mz,levelA_res$time),]

if(length(which(duplicated(mzid)==TRUE))>0){
levelA_res<-levelA_res[-which(duplicated(mzid)==TRUE),]

dataA<-dataA[-which(duplicated(mzid)==TRUE),]
}

levelA_res<-levelA_res[order(levelA_res$Module_RTclust),]
levels_A<-levels(unique(levelA_res$Module_RT))



module_num<-gsub(levelA_res$Module_RTclust,pattern="_[0-9]{1,}",replacement="")


levelA_res_all<-levelA_res[,c(1:2)]

orig_module_labels<-levelA_res$Module_RTclust

levelA_res_all$Module_RTclust<-module_num




mzdefect<-1*((levelA_res$mz-floor(levelA_res$mz)))



d1<-density(mzdefect,bw="nrd",from=min(mzdefect),to=(0.01+max(mzdefect)))


 
 t1<-d1$x
 
 	   
	   gid<-{}
	   
	   gid<-paste("ISgroup",dim(dataA)[1],sep="")
  
 	   levelA_res2<-cbind(mzdefect,levelA_res)
  

	massdefect_cor_groups<-sapply(list(myData1=levelA_res2),function(x)  split(x,cut(levelA_res2$mzdefect,breaks=seq(0,1,0.01))))
	
	if(length(which(levelA_res2$mzdefect==0)>0)){
		massdefect_cor_groups[[length(massdefect_cor_groups)+1]]<-levelA_res2[which(levelA_res2$mzdefect==0),]
	}
	diffmatB<-{}
	#gnum in 
	diffmatB<-lapply(1:length(massdefect_cor_groups),function(gnum){
		
		cur_group<-{}
		
	 	if(dim(massdefect_cor_groups[[gnum]])[1]>0){
	 		 	
	 	ISgroup<-paste("ISgroup",massdefect_cor_groups[[gnum]]$Module_RTclust,gnum,sep="_")
	 	
	 	#print(massdefect_cor_groups[[gnum]])
	 	cur_group<-as.data.frame(massdefect_cor_groups[[gnum]])
	 	cur_group<-cbind(ISgroup,cur_group)
	 	
	 
	 	if(length(cur_group)>0){
	 		
	 		cur_group<-cur_group[order(cur_group$mz,cur_group$time),]
	 		
	 		
	 		if(length(diffmatB)<1){
	 			#diffmatB<-cur_group
	 			}else{
	 	
	 
	 		if(dim(cur_group)[1]>0){
			#diffmatB<-rbind(diffmatB,cur_group)
	 		}
		}
	 
	 	}
	 	}else{
			print(gnum)
		}
		return(cur_group)
	 })
	 
	 
	 diffmatB<-ldply(diffmatB,rbind)
	 diffmatB<-as.data.frame(diffmatB)
 

	 
	 if(dim(diffmatB)[1]>0){
	 cnames<-colnames(diffmatB)
	 cnames[1]<-"ISGroup"
	
	 colnames(diffmatB)<-cnames
	}

  
if(dim(diffmatB)[1]<dim(dataA)[1]){
diffmatC<-levelA_res2[-which(levelA_res$mz%in%diffmatB$mz),]

if(nrow(diffmatC)>0){
t1<-table(diffmatB[,1])
isop_last<-paste("ISgroup_",diffmatC$Module_RTclust,"_",1,sep="")

diffmatC<-cbind(isop_last,diffmatC)
colnames(diffmatC)<-colnames(diffmatB)


diffmatD<-rbind(diffmatB[,c(1:10)],diffmatC[,c(1:10)])

rm(diffmatC)
rm(diffmatB)
}else{
	diffmatD<-diffmatB#[,c(1:10)]
	rm(diffmatB)
}
}else{
	
    diffmatD<-diffmatB#[,c(1:10)]
	rm(diffmatB)
}
rm(levelA_res2)
diffmatD<-as.data.frame(diffmatD)
diffmatD$mz<-as.numeric(diffmatD$mz)

diffmatD<-diffmatD[order(diffmatD$mz,diffmatD$time),]

levelA_res<-levelA_res[order(levelA_res$mz,levelA_res$time),]

levelA_res1<-cbind(diffmatD[,1],levelA_res)

###this didn't work
#print(head(levelA_res1))

 #mean_int_vec<-apply(levelA_res1[,-c(1:4)],1,function(x){mean(x,na.rm=TRUE)})

isop_res_md<-cbind(diffmatD[,c(4,5,1,3)],mean_int_vec,diffmatD[,c(2)])
	 
colnames(isop_res_md)<-c("mz","time","ISgroup","Module_RTclust","AvgIntensity","MD")

#print(isop_res_md[1:3,])





MD<-diffmatD[,c(2)]
levelA_res1<-cbind(levelA_res1[,c(1:4)],mean_int_vec,MD)
rm(MD)

cnames<-colnames(levelA_res1)
cnames[1]<-"ISgroup"

colnames(levelA_res1)<-cnames



multiresmat<-merge(levelB_res,levelA_res1,by="mz")

multiresmat<-multiresmat[,c("mz","time.x","MatchCategory","theoretical.mz","chemical_ID","Name","Formula","MonoisotopicMass","Adduct","ISgroup","Module_RTclust","time.y","mean_int_vec","MD")]

colnames(multiresmat)<-c("mz","time","MatchCategory","theoretical.mz","chemical_ID","Name","Formula","MonoisotopicMass","Adduct","ISgroup","Module_RTclust","time.y","mean_int_vec","MD")

rm(levelB_res)

multiresmat<-multiresmat[order(multiresmat$Module_RTclust),]

rm(m1)



t2<-table(multiresmat$mz)

same1<-which(t2==1)

uniquemz<-names(same1)

multiresmat$MatchCategory=rep("Multiple",dim(multiresmat)[1])

multiresmat$MatchCategory[which(multiresmat$mz%in%uniquemz)]<-"Unique"



setwd(outloc)




multiresmat<-multiresmat[order(multiresmat$Module_RTclust,multiresmat$chemical_ID),]



dupmz<-multiresmat$mz[which(duplicated(multiresmat$mz)==TRUE)]

tablemz<-table(multiresmat$Module_RTclust,multiresmat$mz)

multiresmat$MatchCategory<-rep("Multiple",dim(multiresmat)[1])

multiresmat$MatchCategory[-which(multiresmat$mz%in%dupmz)]<-"Unique"

if(length(which(multiresmat$chemical_ID=="-"))>0){
multiresmat<-multiresmat[-which(multiresmat$chemical_ID=="-"),]



}

tall<-table(multiresmat$chemical_ID,multiresmat$mz)

tall_mznames<-colnames(tall)



tall_checkmultichems<-apply(tall,2,sum)

tall_multichems<-tall[,which(tall_checkmultichems>1)]

names_tall_multichems<-colnames(tall_multichems)

tall_checkmultimz<-apply(tall,1,sum)

tall_multimzperchem<-tall[which(tall_checkmultimz>1),]
tall_unimzperchem<-tall[which(tall_checkmultimz==1),]

names_tall_multimzperchem<-rownames(tall_multimzperchem)

chemids<-names_tall_multimzperchem #chemids[which(t1chemids>0)]

single_mz_chemids<-rownames(tall_unimzperchem)


chemids<-chemids
chem_score<-new("list")
i<-1


del_mz<-{}


multiresmat_filt<-{}


cnames<-colnames(dataA)
 cnames[2]<-"time"
 colnames(dataA)<-as.character(cnames)
dataA<-unique(dataA)

cnames<-colnames(multiresmat)
cnames[10]<-"ISgroup"
colnames(multiresmat)<-cnames




level_module_isop_annot<-levelA_res1 


chem_ids<-table(multiresmat$chemical_ID)
chem_ids_1<-chem_ids[which(chem_ids>=1)]
chem_ids_1<-names(chem_ids_1)

chem_ids_2<-chem_ids_1



chemids<-chem_ids_2


cnames<-colnames(multiresmat)
cnames[2]<-"time"
 colnames(multiresmat)<-cnames


rm(levelA_res)
rm(levelB_res)
rm(m2)

mchemdata<-multiresmat   

rm(multiresmat)

mchemdata<-as.data.frame(mchemdata)


mchemdata$mz<-as.numeric(as.character(mchemdata$mz))
mchemdata$time<-as.numeric(as.character(mchemdata$time))


bad_rows<-which(abs(mchemdata$time-mchemdata$time.y)>0)

if(length(bad_rows)>0){

	mchemdata<-mchemdata[-bad_rows,]

}


chemids<-unique(mchemdata$chemical_ID)


 
chemids<-unique(chemids)


chemscoremat<-{}

if(length(chemids)>10){
num_sets<-length(chemids)/2
	
}else{

num_sets<-1
}

list_winsize<-num_sets

list_size<-round(length(chemids)/list_winsize)



if(length(chemids)>list_winsize){
	
	g<-seq(1,length(chemids),list_size)
	g<-factor(g)
	chemids_split<-split(1:length(chemids),f=g)
	split_size<-1:list_winsize
	
}else{
	chemids_split<-split(1:length(chemids),f=length(chemids))
	split_size<-c(1)
	}


num_sets=length(chemids_split)

#if(FALSE)
{
mchemdata$mz<-round(as.numeric(as.character(mchemdata$mz)),5)
mchemdata$time<-round(as.numeric(as.character(mchemdata$time)),1)
mchemdata$MD<-round(as.numeric(as.character(mchemdata$MD)),3)
mchemdata$mean_int_vec<-round(as.numeric(as.character(mchemdata$mean_int_vec)),1)
mchemdata$time.y<-round(as.numeric(as.character(mchemdata$time.y)),1)
}



if(max_diff_rt>=9999){
	module_num<-gsub(mchemdata$Module_RTclust,pattern="_[0-9]{1,}",replacement="")
	module_num<-paste(module_num,"_0",sep="")
	mchemdata$Module_RTclust<-module_num
}
#write.table(mchemdata,file="Stage2.txt",sep="\t",row.names=FALSE)
write.csv(mchemdata,file="Stage2.csv",row.names=FALSE)

#print("Memory used after step1")
#print(mem_used())

save(list=c("outloc","adduct_weights","boostIDs","pathwaycheckmode","adduct_table","max_diff_rt","corthresh","filter.by",
"max.rt.diff","max_isp","min_ions_perchem","max.mz.diff","db_name","allsteps","redundancy_check","num_sets"),file="tempobjects.Rda")

save(list=c("mchemdata","chemids","adduct_table","mzid","max_diff_rt","isop_res_md","corthresh","level_module_isop_annot",
"chemids_split","corthresh","max.mz.diff","outloc","num_sets","db_name","num_nodes","num_sets","adduct_weights","filter.by",
"max.rt.diff","max_isp","MplusH.abundance.ratio.check","mass_defect_window","mass_defect_mode","allsteps"),file="step1_results.Rda")

rm(mchemdata)
rm(chemids)
rm(mzid)
try(rm(global_cor),silent=TRUE)
rm(isop_res_md)
rm(level_module_isop_annot)
rm(dataA)

rm(tall_unimzperchem)
rm(tablemz)
rm(levelB_res2)
rm(levelA_res1)

rm(list=ls())
load("tempobjects.Rda")
}else{


print("Status 1: Skipping step 1.")
print("Status 2: Using existing step1_results.Rda file.")

    allsteps_temp<-allsteps
	load("tempobjects.Rda")
    allsteps<-allsteps_temp
}





if(allsteps==TRUE){

print("Status 3: Calculating scores for individual chemicals/metabolites")



if(num_sets>num_nodes)
{


cl<-makeSOCKcluster(num_nodes)

clusterEvalQ(cl, "multilevelannotationstep2")

clusterExport(cl, "multilevelannotationstep2")

clusterEvalQ(cl, "library(Rdisop)")
clusterEvalQ(cl, "library(plyr)")
#clusterEvalQ(cl, "library(pryr)")
#clusterEvalQ(cl, "library(profmem)")
#clusterEvalQ(cl, "library(gdata)")
clusterExport(cl, "get_chemscorev1.6.71")
clusterExport(cl, "getMolecule")
clusterExport(cl, "ldply")
#clusterExport(cl, "mem_used")
clusterExport(cl, "get_confidence_stage2")

#clusterExport(cl, "ll")
clusterExport(cl, "check_element")
clusterExport(cl, "group_by_rt_histv2")
clusterExport(cl, "adduct_table")
clusterExport(cl, "adduct_weights")

    parLapply(cl,1:num_sets,function(arg1){
        
        cur_fname<-paste(outloc,"/stage2/chem_score",arg1,".Rda",sep="")
        check_if_exists<-suppressWarnings(try(load(cur_fname)))
        
        if(is(check_if_exists,"try-error")){
        
      
	#suppressWarnings(multilevelannotationstep2(outloc1=outloc,list_number=arg1))
        multilevelannotationstep2(outloc1=outloc,list_number=arg1)
        }else{
            
            print(paste("List ",arg1," already exists.",sep=""))
        }
        
    })
  
        #max.time.diff=max.rt.diff,filter.by=filter.by,max_isp=max_isp,numnodes=1,MplusH.abundance.ratio.check=MplusH.abundance.ratio.check,mass_defect_window=mass_defect_window,mass_defect_mode=mass_defect_mode
    
    stopCluster(cl)
}else{

	
	for(arg1 in 1:num_sets){
		multilevelannotationstep2(outloc1=outloc,list_number=arg1)
	}

}


setwd(outloc)

#print("Memory used after step2 before removing all objects")
#print(mem_used())
rm(list=ls())

try(rm(hmdbCompMZ),silent=TRUE)

load("tempobjects.Rda")

#print("Memory used after step2 and removing all objects")
#print(mem_used())
#print(ll()[order(ll()$KB),])

#if(allsteps==TRUE)
{
print("Status 4: Pathway evaluation")
#suppressWarnings(annotres<-multilevelannotationstep3(outloc=outloc,adduct_weights=adduct_weights,boostIDs=boostIDs,pathwaycheckmode=pathwaycheckmode))
annotres<-multilevelannotationstep3(outloc=outloc,adduct_weights=adduct_weights,boostIDs=boostIDs,pathwaycheckmode=pathwaycheckmode)

#print("Memory used after step3")
#print(mem_used())

setwd(outloc)
rm(list=ls())

#print("Memory used before step4")
#print(mem_used())

try(rm(hmdbCompMZ),silent=TRUE)
try(rm(hmdbCompMZ,env=.GlobalEnv),silent=TRUE)

load("tempobjects.Rda")




#	size_objects<-sort(sapply(ls(), function(x) format(object.size(get(x)), unit = 'auto')))
#	print(size_objects)

print("Status 5: Assigning confidence levels")
#suppressWarnings(annotresstage4<-multilevelannotationstep4(outloc=outloc,max.rt.diff=max_diff_rt,filter.by=filter.by,adduct_weights=adduct_weights,min_ions_perchem=min_ions_perchem,boostIDs=boostIDs,max_isp=max_isp,max.mz.diff=max.mz.diff))
annotresstage4<-multilevelannotationstep4(outloc=outloc,max.rt.diff=max_diff_rt,filter.by=filter.by,adduct_weights=adduct_weights,
min_ions_perchem=min_ions_perchem,boostIDs=boostIDs,max_isp=max_isp,
max.mz.diff=max.mz.diff)

#print("Memory used after step4")
#print(mem_used())

rm(list=ls())

load("tempobjects.Rda")




try(rm(hmdbAllinf,env=.GlobalEnv),silent=TRUE)



if(redundancy_check==TRUE){
    
    print("Status 6:Redundancy filtering")
    
	rm(list=ls())

load("tempobjects.Rda")

    suppressWarnings(annotresstage5<-multilevelannotationstep5(outloc=outloc,max.rt.diff=max_diff_rt,filter.by=filter.by,adduct_weights=adduct_weights,
    min_ions_perchem=min_ions_perchem,boostIDs=boostIDs,max_isp=max_isp,db_name=db_name,max.mz.diff=max.mz.diff))
    
    #print("Stage 5 confidence level distribution for unique chemical/metabolite IDs")
    #print(table(annotresstage5$Confidence[-which(duplicated(annotresstage5$chemical_ID)==TRUE)]))

    }else{
    
        annotresstage5<-{}
    }
}

}else{

annotres<-{} #mchemdata
}



}


print("################")
print("Final: Processing complete")


print("Output files description:")
print("Stage 1 includes clustering of features based on intensity and retention time without annotations")
print("Stage 2 includes clustering results along with simple m/z based database matching")

if(allsteps==TRUE){
    print("Stage 3 includes scores for annotations assigned in stage 2 based on multiple criteria")
    print("Stages 4 and 5 include the confidence levels before and after redundancy (multiple matches) filtering, respectively")
}

suppressWarnings(sink(file=NULL))

setwd(outloc)
sink(file="Readme.txt")
print("Output files description:")
print("Stage1.csv includes clustering of features based on intensity and retention time without annotations")
print("DBresults/Stage2.csv includes clustering results along with simple m/z based database matching")

if(allsteps==TRUE){
    print("DBresults/Stage3.csv includes scores for annotations assigned in stage 2 based on multiple criteria")
    print("DBresults/Stage4.csv and DBresults/Stage5.csv include the confidence levels before and after redundancy (multiple matches) filtering, respectively")
}
suppressWarnings(sink(file=NULL))

if(allsteps==TRUE){
    #return(list("stage4"=annotresstage4,"stage5"=annotresstage5))
}

}
