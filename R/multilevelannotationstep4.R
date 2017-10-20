multilevelannotationstep4 <-
function(outloc,max.mz.diff=5,max.rt.diff=30,adduct_weights=NA,filter.by=NA,min_ions_perchem=1,boostIDs=NA,max_isp=5,dbAllinf=NA,num_nodes=2){
	
	setwd(outloc)


	max_diff_rt=max.rt.diff
	
#chemscoremat_highconf<-read.table("Stage3.txt",sep="\t",header=TRUE) #read.table("Stage3.txt",sep="\t",header=TRUE)

chemscoremat_highconf<-read.csv("Stage3.csv")
	
chemscoremat_highconf<-as.data.frame(chemscoremat_highconf)
chemscoremat_highconf$mz<-as.numeric(chemscoremat_highconf$mz)

if(FALSE){
chemscoremat_highconf$Name<-gsub(chemscoremat_highconf$Name,pattern="[\\\"\']",replacement="")

chemscoremat_highconf$Formula<-gsub(chemscoremat_highconf$Formula,pattern="[\\\"\']",replacement="")

chemscoremat_highconf$Adduct<-gsub(chemscoremat_highconf$Adduct,pattern="[\\\"\']",replacement="")
}
cnames<-colnames(chemscoremat_highconf)

cnames<-gsub(cnames,pattern=".x",replacement="")

colnames(chemscoremat_highconf)<-cnames

chemids<-chemscoremat_highconf$chemical_ID

chemids<-unique(chemids)



chemscoremat_conf_levels<-rep("High",length(chemids))


data(adduct_table)
if(is.na(adduct_weights)==TRUE){
data(adduct_weights)
adduct_weights1<-matrix(nrow=2,ncol=2,0)
		adduct_weights1[1,]<-c("M+H",1)
		adduct_weights1[2,]<-c("M-H",1)
		adduct_weights<-as.data.frame(adduct_weights1)
		colnames(adduct_weights)<-c("Adduct","Weight")
	}	
adduct_table<-adduct_table[order(adduct_table$Adduct),]


chemscoremat_conf_levels<-{}

#cl<-makeSOCKcluster(num_nodes)

#if(length(chemids)>1000)
{

winsize=500

num_itrs<-round(length(chemids)/winsize)

#s1<-seq(1,length(chemids),num_itrs)

#for(i in s1)
{

#if(FALSE)
{
cl<-makeSOCKcluster(num_nodes)


clusterEvalQ(cl, "get_confidence_stage4")
clusterEvalQ(cl, "check_element")
clusterEvalQ(cl, "group_by_rt")

clusterExport(cl, "multilevelannotationstep2")

clusterEvalQ(cl, "library(Rdisop)")
clusterEvalQ(cl, "library(plyr)")
#clusterEvalQ(cl, "library(pryr)")
#clusterEvalQ(cl, "library(profmem)")
#clusterEvalQ(cl, "library(gdata)")

clusterExport(cl, "getMolecule",envir=environment())
clusterExport(cl, "ldply",envir=environment())
clusterExport(cl, "max_isp",envir=environment())
clusterExport(cl, "get_confidence_stage4",envir=environment())

clusterExport(cl, "min_ions_perchem",envir=environment())
clusterExport(cl, "check_element",envir=environment())
clusterExport(cl, "group_by_rt",envir=environment())
clusterExport(cl, "adduct_table",envir=environment())
clusterExport(cl, "adduct_weights",envir=environment())
clusterExport(cl,"filter.by",envir=environment())
clusterExport(cl,"max_diff_rt",envir=environment())
clusterExport(cl,"chemscoremat_highconf",envir=environment())
clusterExport(cl,"chemids",envir=environment())
}
   # parLapply(cl,1:num_sets,function(arg1){
#chemscoremat_conf_levels2<-lapply(i:(i+winsize),function(c)
#chemscoremat_conf_levels2<-sapply(1:length(chemids),function(c)

 #chemscoremat_conf_levels2<-parLapply(cl,1:length(chemids),function(c)
chemscoremat_conf_levels<-foreach(c=1:length(chemids), .combine=rbind) %dopar%
  {
        
        cur_chemid<-chemids[c]
        
        curdata<-chemscoremat_highconf[which(chemscoremat_highconf$chemical_ID==cur_chemid),]
       
	bool_check=1;


        curdata<-curdata[order(curdata$Adduct),]
        
        if((is.na(filter.by)==FALSE) && (bool_check==1)){
            
            check_adduct<-which(curdata$Adduct%in%filter.by)
            if(length(check_adduct)>0){
                bool_check=1;
            }else{
		bool_check=0;
	    }	
      
        }
        
        if(bool_check==1){
        
	    #print(cur_chemid)
	    #print(curdata)
	    final_res<-get_confidence_stage4(curdata,max_diff_rt,adduct_weights=adduct_weights,filter.by=filter.by,max_isp=max_isp,min_ions_perchem=min_ions_perchem)
          
		Confidence<-0 
		#print(final_res) 
        	if(final_res!="None"){    
            if(is.na(final_res[1,1])==FALSE){
                
                
                Confidence<-as.numeric(as.character(final_res[,1]))
                
                curdata<-final_res #[,-c(1)]
                
                rm(final_res)
                if(Confidence<2){
                    
                    if(length(which(curdata$Adduct%in%adduct_weights[which(as.numeric(adduct_weights[,2])>0),1]))>0){
                        
                        if(curdata$score>10){
        		               
					mnum<-max(as.numeric(as.character(adduct_weights[which(adduct_weights[,1]%in%curdata$Adduct),2])))[1] 
				   curdata<-curdata[which(curdata$Adduct%in%adduct_weights[which(as.numeric(as.character(adduct_weights[,2]))>=mnum),1]),]
				    Confidence<-2
                        }
                    }
                }
                
            }
		}
        }else{
            Confidence<-0
			if(length(which(curdata$Adduct%in%adduct_weights[,1]))>0){

                                if(curdata$score>=10){

                                        #curdata<-curdata[which(curdata$Adduct%in%adduct_weights[,1]),]
                                        mnum<-max(as.numeric(as.character(adduct_weights[which(adduct_weights[,1]%in%curdata$Adduct),2])))[1]
                                   #curdata<-curdata[which(curdata$Adduct%in%adduct_weights[which(as.numeric(as.character(adduct_weights[,2]))>=mnum),1]),]
                                    	if(length(which(curdata$Adduct%in%filter.by))>0){
						curdata<-curdata[which(curdata$Adduct%in%filter.by),]
						Confidence<-2
					}

                                }
                        }

        }

	if(nrow(curdata)>1){
	 if(curdata$score<10){
        if(length(unique(curdata$Adduct))<2){
        
                Confidence<-0
        }else{
		if(FALSE){
		 if(Confidence<2){

                    if(length(which(curdata$Adduct%in%adduct_weights[which(adduct_weights[,2]>1),1]))>0){

                        if(curdata$score>10){
                            #Confidence<-2
                        
				mnum<-max(as.numeric(as.character(adduct_weights[which(adduct_weights[,1]%in%curdata$Adduct),2])))[1]
                                   curdata<-curdata[which(curdata$Adduct%in%adduct_weights[which(as.numeric(as.character(adduct_weights[,2]))>=mnum),1]),]
                                    Confidence<-2

			}
                    }
                }
		}
	}
        }
	}
        
	
        curdata<-cbind(Confidence,curdata)
        curdata<-as.data.frame(curdata)
        
	curdata<-curdata[,c("Confidence","chemical_ID")]
	curdata<-unique(curdata)
	  #curdata<-as.matrix(curdata)
	#print(mem_used())
        #rm(curdata)
        return(curdata)
	#return(Confidence)
    } #)
#chemids<-c("HMDB00277","HMDB00222","HMDB00043")

#stopCluster(cl)

#save(list=c("chemscoremat_conf_levels1","chemids"),file="stage4conf_levels.Rda")

chemscoremat_conf_levels<-as.data.frame(chemscoremat_conf_levels)
#save(chemscoremat_conf_levels2,file="stage4conf_levelsB.Rda")
#print(mem_used())
 
 # stop("done saving")
   chemscoremat_conf_levels<-chemscoremat_conf_levels #ldply(chemscoremat_conf_levels2,rbind) #unlist(chemscoremat_conf_levels)
    #rm(chemscoremat_conf_levels2)
    #chemscoremat_conf_levels<-rbind(chemscoremat_conf_levels,chemscoremat_conf_levels_temp)
}


}
    #chemscoremat_highconf<-as.data.frame(chemscoremat_highconf)
    
    chemscoremat_highconf<-unique(chemscoremat_highconf)
    chemscoremat_conf_levels<-as.data.frame(chemscoremat_conf_levels) #[,c(1,3:13,15:16)])
    #print(dim(chemscoremat_conf_levels))
    #save(chemscoremat_conf_levels,file="chemscoremat_conf_levels.Rda")
    #save(chemscoremat_highconf,file="chemscoremat_highconf.Rda")
    
    #chemconf_levels<-cbind(chemscoremat_conf_levels,chemids)
    #chemconf_levels<-as.data.frame(chemconf_levels)
    #colnames(chemconf_levels)<-c("Confidence","chemids") 
    #  write.table(chemconf_levels,file="confidence_levels_chemicals.txt",sep="\t",row.names=FALSE)
    
    #curated_res<-cbind(chemscoremat_conf_levels[,1],chemscoremat_highconf)  
    curated_res<-merge(chemscoremat_conf_levels,chemscoremat_highconf,by="chemical_ID")
    
    cnames<-colnames(curated_res)
    #cnames[3]<-"score" #"Confidence"
    colnames(curated_res)<-as.character(cnames)
    
    rm(chemscoremat_highconf)
  
	curated_res_isp_check<-gregexpr(text=curated_res$Adduct,pattern="(_\\[(\\+|\\-)[0-9]*\\])")

isp_ind<-which(curated_res_isp_check>0)

if(FALSE){
if(length(isp_ind)>0){

        formula_vec<-curated_res$Formula  #aregexpr(text=curated_res$Formula,pattern="(_\\[(\\+|\\-)[0-9]*\\])")
        formula_vec<-gsub(x=formula_vec,pattern="[a-z|+|-|0-9]*_",replacement="",ignore.case=T)
        #print(formula_vec)
        #print(formula_vec[isp_ind])
        #temp_adducts<-as.character(paste("M","_",formula_vec[isp_ind],sep=""))
        curated_res$Adduct<-as.character(curated_res$Adduct)
        #print(curated_res$Adduct)
        #print(temp_adducts)
        #print(curated_res$Adduct[isp_ind])
        curated_res$Adduct[isp_ind]<-as.character(temp_adducts)
        #print(curated_res$Adduct[isp_ind])
}
}

#write.table(curated_res,file="confidence_levels_chemicals.txt",sep="\t",row.names=FALSE)

cnames<-colnames(curated_res)
    
    
    outloc3<-outloc
  suppressWarnings(dir.create(outloc3))
    setwd(outloc3)



if(FALSE)
{
    d2<-read.table("Stage2.txt",sep="\t",header=TRUE)
    
  d3pres<-d2[-which(d2$mz%in%curated_res$mz),]    
    d3pres<-d3pres[-which(d3pres$chemical_ID%in%curated_res$chemical_ID),]
    d3pres<-d3pres[which(d3pres$Adduct%in%adduct_weights[,1]),]

    conf_vec<-rep(1,length(d3pres[,1]))
    score_vec<-rep(1,length(d3pres[,1]))
    d6pres<-cbind(conf_vec,d3pres[,5],score_vec,d3pres[,11],d3pres[,1:4],d3pres[,6:10],d3pres[,13:14])

    print(curated_res[1:2,1:15])
    print(d6pres[1:2,1:15])
	print(dim(d6pres))
	print(dim(curated_res))

   curated_res<-curated_res[,c(1:15)]

   colnames(d6pres)<-colnames(curated_res)    
   curated_res<-rbind(curated_res,d6pres)
}
   curated_res<-as.data.frame(curated_res)

  
    
    cnames1<-colnames(curated_res)
    



      curated_res<-as.data.frame(curated_res)

     curated_res$mz<-as.numeric(as.character(curated_res$mz))

    
    curated_res$theoretical.mz<-as.numeric(as.character(curated_res$theoretical.mz))
    
    curated_res_temp<-curated_res[,c("mz","theoretical.mz")]


    curated_res_temp<-apply(curated_res_temp,1,as.numeric)
    curated_res_temp<-t(curated_res_temp)
	curated_res_temp<-as.data.frame(curated_res_temp)

    delta_ppm<-apply(curated_res_temp,1,function(x){
    
        ppmerror=10^6*abs(x[2]-x[1])/(x[2]);
        return(ppmerror);
    })
    delta_ppm<-round(delta_ppm,2)
    
    
    
    curated_res<-cbind(curated_res[,1:8],delta_ppm,curated_res[,9:dim(curated_res)[2]])

    curated_res<-curated_res[order(curated_res$Confidence,decreasing=TRUE),]
   
    if(is.na(boostIDs)==FALSE){
    
    cnames_boost<-colnames(boostIDs)
    
    if(length(cnames_boost)>1){
    curated_res_mzrt<-curated_res[,c("mz","time")]
    validated_mzrt<-boostIDs[,c("mz","time")]

ghilicpos<-getVenn(curated_res_mzrt, name_a="exp", validated_mzrt, name_b="boost", mz.thresh = max.mz.diff, time.thresh=max_diff_rt,
alignment.tool=NA, xMSanalyzer.outloc=getwd(),use.unique.mz=FALSE,plotvenn=FALSE)

save(ghilicpos,file="ghilicpos.Rda")

g1<-ghilicpos$common
rm(ghilicpos)

if(is.na(max_diff_rt)==FALSE){
#g1<-g1[order(g1$index.B,g1$time.difference),]

t1<-table(g1$index.B)

ind_names<-names(t1)

parent_bad_ind<-{}

if(FALSE){
g2<-{}

for(i1 in ind_names){

	temp1<-g1[which(g1$index.B==i1),]
	bad_ind<-which(temp1$time.difference>min(temp1$time.difference)[1])  #which(g1$time.difference[temp1]>min(g1$time.difference[temp1])[1])
	if(length(bad_ind)>0){
	g2<-rbind(g2,temp1[-bad_ind,])
	}else{
		g2<-rbind(g2,temp1)
	}
	#parent_bad_ind<-c(parent_bad_ind,bad_ind)
}
}

#dup_index_g1<-parent_bad_ind #which(duplicated(g1$index.B)==TRUE)
#if(length(dup_index_g1)>0){
#	g1<-g1[-dup_index_g1,]
#}

}

t1<-table(curated_res$Confidence,curated_res$chemical_ID)
cnames<-colnames(t1)
cnames<-cnames[which(cnames%in%boostIDs$ID)]

good_ind_1<-{}

for(ind2 in 1:dim(g1)[1]){

	temp_ind1<-g1$index.A[ind2]
	temp_ind2<-g1$index.B[ind2]
	
	
	if(curated_res$chemical_ID[temp_ind1]%in%boostIDs$ID[temp_ind2]){
	
		good_ind_1<-c(good_ind_1,g1$index.A[ind2])
	}
}

overlap_mz_time_id<-good_ind_1 #which(bool_vec3==2) #mz_time_index #which(curated_res$chemical_ID%in%good_ids)

curated_res$Confidence[overlap_mz_time_id]<-4
curated_res$score[overlap_mz_time_id]<-curated_res$score[overlap_mz_time_id]*100
t1<-table(curated_res$Confidence[overlap_mz_time_id],curated_res$chemical_ID[overlap_mz_time_id])

cnames1<-colnames(t1)
cnames2<-cnames1[which(t1>0)]
good_ind<-{} #which(curated_res$chemical_ID%in%cnames2)
if(length(good_ind)>0){
curated_res$Confidence[good_ind]<-4
curated_res$score[good_ind]<-curated_res$score[good_ind]*100
}

}else{
	
		good_ind<-which(curated_res$chemical_ID%in%boostIDs)
		if(length(good_ind)>0){
			curated_res$Confidence[good_ind]<-4
			curated_res$score[good_ind]<-curated_res$score[good_ind]*100
		}
}

}
  t2<-table(curated_res$mz)
    
    same1<-which(t2==1)
    
    uniquemz<-names(same1)
    
    curated_res$MatchCategory=rep("Multiple",dim(curated_res)[1])
    
    curated_res$MatchCategory[which(curated_res$mz%in%uniquemz)]<-"Unique"
    
    #write.table(curated_res,file="Stage4.txt",sep="\t",row.names=FALSE)

    write.csv(curated_res,file="Stage4.csv",row.names=FALSE)

    
    
    if(FALSE){
    fname=paste("Stage4_annotation_results",sep="")
    unlink(fname)
    
    HTMLInitFile(filename=fname,Title="Stage 4 annotation results", outdir=outloc)
    fname=paste(outloc,"/Stage4_annotation_results.html",sep="")
    HTML(curated_res,file=fname,Border=1,innerBorder=1,useCSS=TRUE)
    HTMLEndFile(file=fname)
    }
    
#print(head(curated_res))
 
curated_res<-as.data.frame(curated_res)
curated_res<-curated_res[order(curated_res$Confidence,decreasing=TRUE),]
   
print("Stage 4 confidence level distribution for unique chemical/metabolite IDs")
print(table(curated_res$Confidence[-which(duplicated(curated_res$chemical_ID)==TRUE)]))

print("Stage 4 confidence level distribution for unique chemical/metabolite formulas")
print(table(curated_res$Confidence[-which(duplicated(curated_res$Formula)==TRUE)]))



#print(table(curated_res$Confidence))   
    return(curated_res)
    
	
}
