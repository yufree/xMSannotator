group_bymassdefect_old <-
function(dataA,pct.mz.diff=10){
	
	
	mzdefect<-dataA$mz-floor(dataA$mz)
	

	#print(mzdefect[1:10])
	#dataA<-cbind(dataA)
	
	 #Step 1 Group features by m/zdim(data_a)[1]
        mz_groups<-lapply(1:dim(dataA)[1],function(j){

                                commat={}
                                commzA=new("list")
                                commzB=new("list")
                                
                                getbind_same<-c(j)
            
                                	
                                ppmb=(pct.mz.diff)*(mzdefect[j]/100)

                                #getbind_same<-c(getbind_same,which(abs(mzdefect-mzdefect[j])<=ppmb))
                                
                                 getbind_same<-which(abs(mzdefect-mzdefect[j])<=ppmb)
                                
 	
        					# gid<-paste("parent",getbind_same,sep="")
        					
        					return(getbind_same)
                           
		     })
		     
		     
	mz_groups<-unique(mz_groups)

  del_list<-{}
	   
	   gid<-{}
	   
	   gid<-paste("MDgroup",dim(dataA)[1],sep="")
	   
	   #length(mz_groups)
	   	  for(m in 1:length(mz_groups)){
	  
	   if((m%in%del_list)==FALSE){
		for(n in (m+1):length(mz_groups)){
		
			if(n>length(mz_groups)){
				break;
			}
			com1<-intersect(mz_groups[[m]],mz_groups[[n]])
			if(length(com1)>0){
			
				mz_groups[[m]]<-c(mz_groups[[m]],mz_groups[[n]])
				del_list<-c(del_list,n)
				mz_groups[[n]]<-c(0)
			}
		}
		
		#mz_groups[[m]][[1]]<-unique(mz_groups[[m]][[1]])
		}
		}
		if(length(del_list)>0){
		 mz_groups<-mz_groups[-del_list]
		}
	  
	  gidmat<-{}
	  for(k in 1:length(mz_groups)){
	  	
	  	for(m in 1:length(mz_groups[[k]])){
	  		
	  	
	  		gid[mz_groups[[k]][m]]<-paste("MDgroup",mz_groups[[k]][1],sep="")
	  		
	  		curres<-cbind(mz_groups[[k]][m],gid[mz_groups[[k]][m]])
	  		
	  		gidmat<-rbind(gidmat,curres)
	  		
	  		}
	  	
	  	}
	  	
colnames(gidmat)<-c("mzindex","MDgroup")
	 gidmat<-as.data.frame(gidmat)
	  
return(gidmat)

#return(mz_groups)
}
