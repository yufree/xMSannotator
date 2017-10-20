group_bymzisotopes <-
function(dataA,max.mz.diff=30,numisotopes=15,mdgnum=""){

dataA<-as.data.frame(dataA)
	 #Step 1 Group features by m/zdim(data_a)[1]
        mz_groups<-lapply(1:dim(dataA)[1],function(j){

                                commat={}
                                commzA=new("list")
                                commzB=new("list")
                                
                                getbind_same<-c(j)
                                
                                
                                
                                for(isotopepattern in c(1:numisotopes)){
                                	
                                	isotopemass=dataA$mz[j]+as.numeric(isotopepattern)
                                	
                                ppmb=(max.mz.diff)*(isotopemass/1000000)

                                getbind_same<-c(getbind_same,which(abs(dataA$mz-isotopemass)<=ppmb))
                              
        					}
        					
        					
        					# gid<-paste("parent",getbind_same,sep="")
        					
        					return(getbind_same)
                           
		     })
		     
		     
	

  del_list<-{}
	   
	   gid<-{}
	   
	   gid<-rep(paste("M",sep=""),dim(dataA)[1])
	   
	   
	   #if(FALSE)
	   {
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
		
		#del_list<-{}
		if(length(del_list)>0){
		 mz_groups<-mz_groups[-del_list]
		}
	  
	  }
	  
	  #ignore_list<-{}
	  
	  if(length(mz_groups)>1){
	  for(k in 1:length(mz_groups)){
	  	
	  	for(m in 1:length(mz_groups[[k]]))
	  	{
	  		#mz_groups[[k]]<-unique(mz_groups[[k]])
	  	
	  		#gid[mz_groups[[k]][m]]<-paste("ISgroup_",mdgnum,"_",mz_groups[[k]][1],sep="")
	  		
	  		gid[mz_groups[[k]][m]]<-paste("ISgroup_",mdgnum,"_",mz_groups[[k]][1],sep="")
	  		
	  		}
	  	
	  	}
	  	
	}else{
		
		for(m in 1:length(mz_groups[[1]]))
	  	{
	  		#mz_groups[[1]]<-unique(mz_groups[[1]])
	  	
	  		#gid[mz_groups[[k]][m]]<-paste("ISgroup_",mdgnum,"_",mz_groups[[k]][1],sep="")
	  		
	  		gid[mz_groups[[1]][m]]<-paste("ISgroup_",mdgnum,"_",mz_groups[[1]][1],sep="")
	  		
	  		}
	  		

		}
	  
	  gid<-cbind(gid,dataA[,c(1:2)])
	  
return(gid)

}
