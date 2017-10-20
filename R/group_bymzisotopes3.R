group_bymzisotopes3 <-
function(dataA,max.mz.diff=30,numisotopes=4){
	
	
	
	 #Step 1 Group features by m/zdim(data_a)[1]
        mz_groups<-lapply(1:dim(dataA)[1],function(j){

                                commat={}
                                commzA=new("list")
                                commzB=new("list")
                                
                                getbind_same<-c(j)
                                
                                
                                
                                for(isotopepattern in c(1:numisotopes)){
                                	
                                	isotopemass=dataA$mz[j]+isotopepattern
                                	
                                ppmb=(max.mz.diff)*(isotopemass/1000000)

                                getbind_same<-c(getbind_same,which(abs(dataA$mz-isotopemass)<=ppmb))
                                
                                
                                
        					}
        					
        					
        					# gid<-paste("parent",getbind_same,sep="")
        					
        					return(getbind_same)
                           
		     })
		     
		     
	   del_list<-{}
	   
	   gid<-{}
	   
	   #length(mz_groups)
	   for(k in 1:length(mz_groups)){
	  
	  
	   if((k%in%del_list)==FALSE){
		for(n in (k+1):length(mz_groups)){
		
			if(n>length(mz_groups)){
				break;
			}
			com1<-intersect(mz_groups[[k]],mz_groups[[n]])
			if(length(com1)>0){
			
				mz_groups[[k]]<-c(mz_groups[[k]],mz_groups[[n]])
				del_list<-c(del_list,n)
				mz_groups[[n]]<-c(0)
			}
		}
		
		#mz_groups[[m]][[1]]<-unique(mz_groups[[m]][[1]])
		}
		
		if(length(del_list)>0){
		 mz_groups<-mz_groups[-del_list]
		}
	  }
	  
	  for(k in 1:length(mz_groups)){
	  	
	  	for(m in 1:length(mz_groups[[k]])){
	  		
	  	
	  		gid[mz_groups[[k]][m]]<-paste("M",mz_groups[[k]][1],sep="")
	  		
	  		}
	  	
	  	}
	  
	  
#return(gid)

return(mz_groups)
}
