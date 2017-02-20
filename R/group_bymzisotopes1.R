group_bymzisotopes1 <-
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
        					 gid<-paste("parent",getbind_same,sep="")
        					
        					return(getbind_same)
                           
		     })

return(mz_groups)
}
