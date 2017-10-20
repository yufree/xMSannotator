multilevelannotationstep3 <-
function(outloc1,adduct_weights=NA,num_sets=NA,boostIDs=NA,pathwaycheckmode="p",dbAllinf=NA,scorethresh=0.1){
	
	setwd(outloc1)


    if(is.na(num_sets)==FALSE){
        
        num_sets_1=num_sets
    }else{
        
        num_sets_1=NA
    }
	load("step1_results.Rda")
	load("chemCompMZ.Rda")
	rm(global_cor)
	
	rm(global_cor,env=.GlobalEnv)
	#rm(dataA)
	setwd(outloc1)

	num_sets<-length(chemids_split)

    if(is.na(num_sets_1)==FALSE){
        num_sets=num_sets_1
    }

	
	rm(dataA)
	#rm(levelA_res1)

	#print(ls())
	
		if(num_sets>=length(chemids_split)){
			num_sets<-length(chemids_split)
			
				
		}

		

	outloc<-outloc1
  

if(is.na(adduct_weights)==TRUE){
 data(adduct_weights)
 adduct_weights<-as.data.frame(adduct_weights)
#print(dim(adduct_weights))

adduct_weights1<-matrix(nrow=2,ncol=2,0)
		adduct_weights1[1,]<-c("M+H",1)
		adduct_weights1[2,]<-c("M-H",1)
		adduct_weights<-as.data.frame(adduct_weights1)
		colnames(adduct_weights)<-c("Adduct","Weight")

}	
outloc1<-paste(outloc,"/stage2/",sep="")                 
suppressWarnings(dir.create(outloc1))
setwd(outloc1)

chemscoremat<-{}

#print("num_sets")
#print(num_sets)
#print(length(chemids_split))

#for(sind in seq(1,num_sets))
chemscoremat<-lapply(1:num_sets,function(sind)
{
	cur_fname<-paste("chem_score",sind,".Rda",sep="")

	try(load(cur_fname),silent=TRUE)

	curchemscoremat<-as.data.frame(curchemscoremat)
	
	curchemscoremat$Formula<-gsub(curchemscoremat$Formula,pattern="_.*",replacement="")
	
	#colnames(curchemscoremat)<-c("cur_chem_score","Module_RTclust","mz","time", "MatchCategory","theoretical.mz","chemical_ID","Name","Formula","MonoisotopicMass","Adduct","ISgroup","time.y","mean_int_vec","MD")
	#curchemscoremat<-curchemscoremat[,1:14]
	#print(dim(curchemscoremat))
	#print(colnames(curchemscoremat))
	#chemscoremat<-rbind(chemscoremat,curchemscoremat)
	return(curchemscoremat)
})

chemscoremat<-ldply(chemscoremat,rbind)
chemscoremat<-as.data.frame(chemscoremat)

tempadduct<-chemscoremat$Adduct

chemscoremat$Adduct<-gsub(chemscoremat$Adduct,pattern="_.*",replacement="")

chemscoremat<-cbind(chemscoremat,tempadduct)
	
chemscoremat<-merge(chemscoremat,chemCompMZ[,c(2:4,6)],by=c("Formula","Adduct"))


#y because we want chemCompMZ ID and Name
chemscoremat<-chemscoremat[,c("cur_chem_score","Module_RTclust","mz","time","MatchCategory","theoretical.mz","chemical_ID.y","Name.y","Formula","MonoisotopicMass","tempadduct","ISgroup","mean_int_vec","MD")]



colnames(chemscoremat)<-c("cur_chem_score","Module_RTclust","mz","time","MatchCategory","theoretical.mz","chemical_ID","Name","Formula","MonoisotopicMass","Adduct","ISgroup","mean_int_vec","MD")


hmdbbad<-c("HMDB29244","HMDB29245","HMDB29246")

if(length(which(chemscoremat$chemical_ID%in%hmdbbad))>0){
	chemscoremat<-chemscoremat[-which(chemscoremat$chemical_ID%in%hmdbbad),]
}

otherchems<-mchemdata[which(chemscoremat$mz%in%chemscoremat)]

cnames<-colnames(chemscoremat)
cnames[1]<-"score"
colnames(chemscoremat)<-cnames




 #corthresh*(1/(max_diff_rt))*2

pthresh=0.05

if(db_name=="KEGG")
{
		data(keggotherinf)

		
        	#keggotherinf<-dbAllinf

		m1<-apply(keggotherinf,1,function(x){
			chemid<-x[1];
			g1<-gregexpr(x[4],pattern="map")
			regexp_check<-attr(g1[[1]],"match.length")
			if(regexp_check[1]<0){pathid="-";
				return(cbind(chemid,pathid))
				}else{
			pathid<-strsplit(x=x[4],split=";");
			pathid<-unlist(pathid)
			return(cbind(chemid,pathid))
			}
		}
	)
		
	m2<-ldply(m1,rbind)
		
	chemscoremat<-merge(chemscoremat,keggotherinf,by.x="chemical_ID",by.y="KEGGID")
			

	#write.table(chemscoremat,file=paste("xMSannotator_",db_name,"scorematotherinf_stage1.txt",sep=""),sep="\t",row.names=FALSE)
#chemscoremat<-read.table("xMSannotator_KEGGscorematotherinf_stage1.txt",sep="\t",header=TRUE)



	
	chemids<-as.character(chemscoremat$chemical_ID)
	
	pathway_ids<-as.character(m2[,2])
	
	pathway_ids<-unique(pathway_ids)

	module_num<-gsub(chemscoremat$Module_RTclust,pattern="_[0-9]*",replacement="")
	
	chemscoremat<-cbind(chemscoremat,module_num)
	
	chemscoremat_orig<-chemscoremat
	
	chemscoremat<-chemscoremat_orig

total_chem_count<-length(unique(m2$chemid))

if(is.na(pathwaycheckmode)==FALSE){	
	

for(path_id in pathway_ids)
	{
			
			if(path_id!="-" & path_id!="map01100"){
			pathway_chemicals<-m2[which(m2[,2]%in%path_id),1] #m2[which(m2[,2]%in%pathwayscurchemical),1]
			
            curmchemicaldata1<-chemscoremat[which(chemscoremat$chemical_ID%in%pathway_chemicals & chemscoremat$score>=scorethresh & chemscoremat$Adduct%in%as.character(adduct_weights[,1])),]
            
            #molecules of interest in pathway (a)
            num_chems_inpath<-length(unique(curmchemicaldata1$chemical_ID)) #^5
            
            #total number of chemicals in pathway
            all_cur_path_numchem<-length(unique(pathway_chemicals))
            
            #non-focus molecules associated with pathway (c)
            num_chem_inpath_notinterest<-all_cur_path_numchem-num_chems_inpath
            
            
            curmchemicaldata2<-chemscoremat[which(chemscoremat$score>=scorethresh & chemscoremat$Adduct%in%as.character(adduct_weights[,1])),]
            
            curmchemicaldata2<-curmchemicaldata2[-which(curmchemicaldata2$chemical_ID%in%pathway_chemicals),]
            
            #focus molecules not associated with this pathway (b)
            num_chems_notinpath<-length(unique(curmchemicaldata2$chemical_ID))
            
            
            all_notcurpath_numchem<-length(m2[-which(m2[,1]%in%curmchemicaldata2$chemical_ID),1]) #length(m2[-which(m2[,2]%in%path_id),1])
            
            #counts=matrix(data=c(num_chems_inpath,num_chems_notinpath,all_cur_path_numchem,all_notcurpath_numchem),nrow=2)
            
            counts=matrix(data=c(num_chems_inpath,num_chem_inpath_notinterest,num_chems_notinpath,all_notcurpath_numchem),nrow=2)
            
	    rm(curmchemicaldata2)
            p1<-fisher.test(counts)
            p1<-p1$p.value
            
            if(p1>pthresh){
                
                next;
            }else{


			t1<-table(curmchemicaldata1$module_num)
			
			module_counts<-t1[which(t1>0)]
			
			module_names<-names(module_counts)
		
			pathway_chemicals_1<-curmchemicaldata1$chemical_ID
	
			for(c in pathway_chemicals_1){
				curmchemicaldata<-curmchemicaldata1[which(as.character(curmchemicaldata1$chemical_ID)==c),]
                		chem_path_data<-m2[which(m2$chemid==c),]
                		total_num_chem<-0
				t2<-table(curmchemicaldata$module_num)	
				cur_module<-names(t2[which(t2==max(t2)[1])])
               			total_num_chem<-length(pathway_chemicals)
                
             			mzid_cur<-paste(curmchemicaldata$mz,curmchemicaldata$time,sep="_")

#dweights<-alldegrees[which(mzid%in%mzid_cur),1]

				if(nrow(curmchemicaldata)>0){
					pathwayscurchemical<-m2[which(m2[,1]==c),2]
					#for(cur_module in cur_modules)
					{					
			
					if(pathwaycheckmode=="pm"){		
					num_chems<-t1[as.character(cur_module)]
					}
                    
                    #	num_chems<-100*(num_chems/total_num_chem)
                    
					#if(num_chems>1 && length(which(dweights>min(alldegrees[,1])))>=1){
					#	num_chems=10
					#}

					num_chems<-round(num_chems,0)
                    
                    num_chems_inmodule<-length(unique(curmchemicaldata1$chemical_ID[which(curmchemicaldata$module_num==cur_module)])) #^5
                    
                    cur_module_data<-chemscoremat[which(chemscoremat$module_num==cur_module & chemscoremat$chemical_ID%in%pathway_chemicals & chemscoremat$score>=scorethresh & chemscoremat$Adduct%in%as.character(adduct_weights[,1])),]
		    a<-length(unique(cur_module_data$chemical_ID))
		     rm(cur_module_data)
		    
                    cur_module_data2<-chemscoremat[which(chemscoremat$module_num==cur_module & chemscoremat$score>=scorethresh & chemscoremat$Adduct%in%as.character(adduct_weights[,1])),]
                     b<-length(unique(cur_module_data2$chemical_ID))-a
		     rm(cur_module_data2)
		     
                    other_module_data<-chemscoremat[which(chemscoremat$module_num!=cur_module & chemscoremat$chemical_ID%in%pathway_chemicals & chemscoremat$score>=scorethresh & chemscoremat$Adduct%in%as.character(adduct_weights[,1])),]
                     c<-length(unique(other_module_data$chemical_ID))
		     rm(other_module_data)
		     
                    other_module_data2<-chemscoremat[which(chemscoremat$module_num!=cur_module & chemscoremat$score>=scorethresh & chemscoremat$Adduct%in%as.character(adduct_weights[,1])),]
                    
                  
                   
                   
                    d<-length(unique(other_module_data2$chemical_ID))-c
                    rm(other_module_data2)
		   
		    
                    counts=matrix(data=c(a,c,b,d),nrow=2)
                    if(a>1){
                        p1<-fisher.test(counts)
                        p1<-p1$p.value
                    }else{
                        
                        p1=1
                    }

                    if(p1>0.2){
                        
                        next;

                    }else{
						
					if(num_chems<3){

						num_chems<-0
					}else{	
					if(is.na(curmchemicaldata$score[1])==TRUE){

						diff_rt<-max(curmchemicaldata$time)-min(curmchemicaldata$time)
						
						if(diff_rt>max_diff_rt){
						if(length(which(t2>1))==1){
							curmchemicaldata$score<-rep(0.1,length(curmchemicaldata$score))
						}else{
							curmchemicaldata$score<-rep(0,length(curmchemicaldata$score))
						}
						}else{
							curmchemicaldata$score<-rep(0,length(curmchemicaldata$score))
						}	
					}
					
					if(curmchemicaldata$score[1]<scorethresh){
								chemscoremat$score[which(as.character(chemscoremat$chemical_ID)==c & chemscoremat$Adduct%in%as.character(adduct_weights[,1]))]=as.numeric(chemscoremat$score[which(as.character(chemscoremat$chemical_ID)==c & chemscoremat$Adduct%in%as.character(adduct_weights[,1]))][1])+num_chems
						
					}else{
						

						chemscoremat$score[which(as.character(chemscoremat$chemical_ID)==c)]=chemscoremat$score[which(as.character(chemscoremat$chemical_ID)==c)][1]+num_chems


					}				
					}
                    }
					}
					#chemscoremat[which(as.character(chemscoremat$chemical_ID)==c),]
				}
            }
			}
			
		   }
		   rm(curmchemicaldata1)
	}
}			

		cnames<-colnames(chemscoremat)

		cnames<-gsub(cnames,pattern=".x",replacement="")

		colnames(chemscoremat)<-cnames
	 
		#write.table(chemscoremat,file=paste("xMSannotator_",db_name,"scorematotherinf_stage2.txt",sep=""),sep="\t",row.names=FALSE)
}
else{
    
    if(db_name=="HMDB")
    {
	data(hmdbAllinf)

    hmdbAllinfv3.5<-hmdbAllinf[,-c(26:27)] #dbAllinf[,-c(26:27)]
	rm(hmdbAllinf,envir=.GlobalEnv)

	#hmdbAllinfv3.5<-gsub(x=hmdbAllinfv3.5,pattern="\n[\\s]+",replacement="") 
	
	#replace(as.matrix(hmdbAllinfv3.5),which(hmdbAllinfv3.5=="\n"),"")

	m1<-apply(hmdbAllinfv3.5,1,function(x){
			chemid<-x[1];
			g1<-gregexpr(x[14],pattern="SM")
			regexp_check<-attr(g1[[1]],"match.length")
			if(regexp_check[1]<0){pathid="-";
				return(cbind(chemid,pathid))
				}else{
			pathid<-strsplit(x=x[14],split=";");
			pathid<-unlist(pathid)
			return(cbind(chemid,pathid))
			}
		}
	)
		
	m2<-ldply(m1,rbind)
		
	chemscoremat<-merge(chemscoremat,hmdbAllinfv3.5,by.x="chemical_ID",by.y="HMDBID")
	
	#write.table(chemscoremat,file=paste("xMSannotator_",db_name,"scoremat_stage1.txt",sep=""),sep="\t",row.names=FALSE)

	rm(hmdbAllinfv3.5)
		
	chemids<-as.character(chemscoremat$chemical_ID)
	
	pathway_ids<-as.character(m2[,2])
	
	pathway_ids<-unique(pathway_ids)

	module_num<-gsub(chemscoremat$Module_RTclust,pattern="_[0-9]*",replacement="")
	
	chemscoremat<-cbind(chemscoremat,module_num)
	
	chemscoremat_orig<-chemscoremat
	
	chemscoremat<-chemscoremat_orig

	total_chem_count<-length(unique(m2$chemid))


if(is.na(pathwaycheckmode)==FALSE){	

	for(path_id in pathway_ids)
	{
    
    if(path_id!="-"){
        pathway_chemicals<-m2[which(m2[,2]%in%path_id),1] #m2[which(m2[,2]%in%pathwayscurchemical),1]
        
        
        curmchemicaldata1<-chemscoremat[which(chemscoremat$chemical_ID%in%pathway_chemicals & chemscoremat$score>=scorethresh & chemscoremat$Adduct%in%as.character(adduct_weights[,1])),]
        
        #molecules of interest in pathway (a)
        num_chems_inpath<-length(unique(curmchemicaldata1$chemical_ID)) #^5
        
        #total number of chemicals in pathway
        all_cur_path_numchem<-length(unique(pathway_chemicals))
        
        #non-focus molecules associated with pathway (c)
        num_chem_inpath_notinterest<-all_cur_path_numchem-num_chems_inpath
        
        
        curmchemicaldata2<-chemscoremat[which(chemscoremat$score>=scorethresh & chemscoremat$Adduct%in%as.character(adduct_weights[,1])),]
        
        curmchemicaldata2<-curmchemicaldata2[-which(curmchemicaldata2$chemical_ID%in%pathway_chemicals),]
        
        #focus molecules not associated with this pathway (b)
        num_chems_notinpath<-length(unique(curmchemicaldata2$chemical_ID))
        
        
        all_notcurpath_numchem<-length(m2[-which(m2[,1]%in%curmchemicaldata2$chemical_ID),1]) #length(m2[-which(m2[,2]%in%path_id),1])
        
        #counts=matrix(data=c(num_chems_inpath,num_chems_notinpath,all_cur_path_numchem,all_notcurpath_numchem),nrow=2)
        
        counts=matrix(data=c(num_chems_inpath,num_chem_inpath_notinterest,num_chems_notinpath,all_notcurpath_numchem),nrow=2)
        
	rm(curmchemicaldata2)
        p1<-fisher.test(counts)
        p1<-p1$p.value
        
        
        
        
        if(p1>pthresh){
            
            next;
        }else{
            
            
            #num_chems<-length(unique(curmchemicaldata$chemical_ID))^5
            
            t1<-table(curmchemicaldata1$module_num)
            
            module_counts<-t1[which(t1>0)]
            
            module_names<-names(module_counts)
            
            pathway_chemicals_1<-curmchemicaldata1$chemical_ID
            
            for(chemname in pathway_chemicals)
            {
                
                
                #  curmchemicaldata<-curmchemicaldata #chemscoremat[which(as.character(chemscoremat$chemical_ID)==c),] #& chemscoremat$Adduct%in%as.character(adduct_weights[,1])),]
                
                
                curmchemicaldata<-curmchemicaldata1[which(as.character(curmchemicaldata1$chemical_ID)==chemname),]
                
                if(nrow(curmchemicaldata)>0){
                    chem_path_data<-m2[which(m2$chemid==chemname),]
                    
                    
                    
                    #cur_module<-max(curmchemicaldata$module_num)[1]
                    t2<-table(curmchemicaldata$module_num)
                    
                    cur_module<-names(t2[which(t2==max(t2)[1])])
                    
                    mzid_cur<-paste(curmchemicaldata$mz,curmchemicaldata$time,sep="_")
                    
                    #dweights<-alldegrees[which(mzid%in%mzid_cur),1]
                    
                    
                    pathwayscurchemical<-m2[which(m2[,1]==chemname),2]
                    #for(cur_module in cur_modules)
                    {
                        if(pathwaycheckmode=="pm"){
                            num_chems<-t1[as.character(cur_module)]
                        }
                        # if(num_chems>1 && length(which(dweights>min(alldegrees[,1])))>=1){
                        #        num_chems=10
                        #}
                        total_num_chem<-length(pathway_chemicals)
                        
                        
                        # num_chems<-100*(num_chems/total_num_chem)
                        
                        num_chems<-round(num_chems,0)
                        
                        num_chems_inmodule<-length(unique(curmchemicaldata1$chemical_ID[which(curmchemicaldata$module_num==cur_module)])) #^5
                        
                        cur_module_data<-chemscoremat[which(chemscoremat$module_num==cur_module & chemscoremat$chemical_ID%in%pathway_chemicals & chemscoremat$score>=scorethresh & chemscoremat$Adduct%in%as.character(adduct_weights[,1])),]
                        a<-length(unique(cur_module_data$chemical_ID))
			rm(cur_module_data)
			
                        cur_module_data2<-chemscoremat[which(chemscoremat$module_num==cur_module & chemscoremat$score>=scorethresh & chemscoremat$Adduct%in%as.character(adduct_weights[,1])),]
                         b<-length(unique(cur_module_data2$chemical_ID))-a
			 rm(cur_module_data2)
			 
                        other_module_data<-chemscoremat[which(chemscoremat$module_num!=cur_module & chemscoremat$chemical_ID%in%pathway_chemicals & chemscoremat$score>=scorethresh & chemscoremat$Adduct%in%as.character(adduct_weights[,1])),]
                        c<-length(unique(other_module_data$chemical_ID))
			rm(other_module_data)
			
                        other_module_data2<-chemscoremat[which(chemscoremat$module_num!=cur_module & chemscoremat$score>=scorethresh & chemscoremat$Adduct%in%as.character(adduct_weights[,1])),]
                        
                        d<-length(unique(other_module_data2$chemical_ID))-c
			rm(other_module_data2)
                        
                        counts=matrix(data=c(a,c,b,d),nrow=2)
                        
                        if(a>1){
                            p1<-fisher.test(counts)
                            p1<-p1$p.value
                        }else{
                            
                            p1=1
                        }
                        
                        
                        # p1=0
                        
                        if(p1>0.2){
                            
                            next;
                            
                        }else{
                            
                            if(num_chems<3){
                                
                                num_chems<-0
                            }else{
                                
                                if(is.na(curmchemicaldata$score[1])==TRUE){
                                    
                                    diff_rt<-max(curmchemicaldata$time)-min(curmchemicaldata$time)
                                    
                                    if(diff_rt>max_diff_rt){
                                        if(length(which(t2>1))==1){
                                            curmchemicaldata$score<-rep(0.1,length(curmchemicaldata$score))
                                        }else{
                                            curmchemicaldata$score<-rep(0,length(curmchemicaldata$score))
                                        }
                                    }else{
                                        curmchemicaldata$score<-rep(0,length(curmchemicaldata$score))
                                    }
                                }
                                
                                #  print("here")
                                #print(curmchemicaldata$score)
                                #print(num_chems)
                                
                                
                                if(curmchemicaldata$score[1]<scorethresh){
                                    
                                    #chemscoremat$score[which(as.character(chemscoremat$chemical_ID)==c & chemscoremat$module_num==cur_module & chemscoremat$Adduct%in%as.character(adduct_weights[,1]))]=as.numeric(chemscoremat$score[which(as.character(chemscoremat$chemical_ID)==c & chemscoremat$Adduct%in%as.character(adduct_weights[,1]))][1])+num_chems
                                    
                                    chemscoremat$score[which(as.character(chemscoremat$chemical_ID)==chemname & chemscoremat$Adduct%in%as.character(adduct_weights[,1]))]=as.numeric(chemscoremat$score[which(as.character(chemscoremat$chemical_ID)==chemname & chemscoremat$Adduct%in%as.character(adduct_weights[,1]))][1])+num_chems
                                    
                                }else{
                                    #chemscoremat$score[which(as.character(chemscoremat$chemical_ID)==c & chemscoremat$module_num==cur_module)]=chemscoremat$score[which(as.character(chemscoremat$chemical_ID)==c)][1]+num_chems
                                    
                                    chemscoremat$score[which(as.character(chemscoremat$chemical_ID)==chemname)]=chemscoremat$score[which(as.character(chemscoremat$chemical_ID)==chemname)][1]+num_chems
                                    
                                    
                                }
                            }
                        }
                    }
                }
                
            }
        }
        
    }
    rm(curmchemicaldata1)
}



	}
	

		cnames<-colnames(chemscoremat)

		cnames<-gsub(cnames,pattern=".x",replacement="")

		colnames(chemscoremat)<-cnames
	 
		#write.table(chemscoremat,file=paste("xMSannotator_",db_name,"scoremat_stage2.txt",sep=""),sep="\t",row.names=FALSE)
    }
    else{
        
        # if(db_name=="T3DB"){
            
            #load("/Users/karanuppal/Documents/Emory/JonesLab/Projects/xMSannotator/t3dbotherinf.Rda")
            #chemscoremat_otherinf<-merge(chemscoremat,t3dbotherinf,by.x="chemical_ID",by.y="T3DB.ID")
            
            #write.table(chemscoremat_otherinf,file=paste("xMSannotator_",db_name,"scoremat_stage1.txt",sep=""),sep="\t",row.names=FALSE)
            #rm(chemscoremat_otherinf)

            
            #}
        
    }
}

 rm(curmchemicaldata1)

cnames<-colnames(chemscoremat)

		cnames<-gsub(cnames,pattern=".x",replacement="")

		colnames(chemscoremat)<-cnames
multiresmat<-chemscoremat #mchemdata_orig

cnames<-colnames(chemscoremat)

#chemscoremat$score[which(chemscoremat$chemical_ID%in%boostIDs)]<-max(chemscoremat$score,na.rm=TRUE)*100



good_ind<-which(chemscoremat$score>=scorethresh)

chemscoremat_highconf<-{}
if(length(good_ind)>0){

chemscoremat_highconf<-chemscoremat #[which(chemscoremat$score>=scorethresh),]
rm(chemscoremat)
}


chemscoremat_highconf<-chemscoremat_highconf[,c("chemical_ID","score","Module_RTclust","mz", "time" ,"MatchCategory","theoretical.mz","Name","Formula","MonoisotopicMass","Adduct","ISgroup",
"mean_int_vec","MD")]

#chemscoremat_highconf1<-chemscoremat_highconf[,c(1:12,14:15)]


write.csv(chemscoremat_highconf,file="../Stage3.csv",row.names=FALSE)

rm(chemscoremat_orig)

rm("mchemdata","chemids","adduct_table","global_cor","mzid","max_diff_rt","isop_res_md","corthresh","level_module_isop_annot",
"chemids_split","dataA","corthresh","max.mz.diff","outloc","num_sets","db_name","num_nodes","num_sets","adduct_weights","filter.by")
rm("chemCompMZ","mchemdata","hmdbAllinf","hmdbAllinfv3.6","dbAllinf")

	
return(chemscoremat_highconf)
}
