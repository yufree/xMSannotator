find.Overlapping.mzs <-
function(dataA, dataB, mz.thresh=10, time.thresh=NA, alignment.tool=NA)
{

        data_a<-as.data.frame(dataA)
        data_b<-as.data.frame(dataB)
	#data_a<-unique(data_a)
	rm(dataA)
	rm(dataB)
        
     #   data_b<-unique(data_b)
        
        com_mz_num=1
        unique_mz={}
        ppm_v={}
        rt_v={}

        commat={}

	col.names.dataA=colnames(data_a)
	col.names.dataB=colnames(data_b)

	if(is.na(alignment.tool)==FALSE){
	 if(alignment.tool=="apLCMS")
        {
              sample.col.start=5
        }
        else
        {
              if(alignment.tool=="XCMS")
              {
                    sample.col.start=9
                    col.names.dataA[1]="mz"
                    col.names.dataA[2]="time"
		    col.names.dataB[1]="mz"
                    col.names.dataB[2]="time"
                    colnames(data_a)=col.names.dataA
                    colnames(data_b)=col.names.dataB
              }
	      
	}}else{
                    #stop(paste("Invalid value for alignment.tool. Please use either \"apLCMS\" or \"XCMS\"", sep=""))
		   
		    col.names.dataA[1]="mz"
		   
		    col.names.dataB[1]="mz"
		    if(is.na(time.thresh)==FALSE){
		     col.names.dataA[2]="time"
		     col.names.dataB[2]="time"
		      #print("Using the 1st columns as \"mz\" and 2nd columsn as \"retention time\"")
		     }else{
		      #print("Using the 1st columns as \"mz\"")
		     }
		    colnames(data_a)=col.names.dataA
                    colnames(data_b)=col.names.dataB
	}

       #data_a<-data_a[order(data_a$mz),]
       #data_b<-data_b[order(data_b$mz),]
       data_a<-as.data.frame(data_a)
	data_b<-as.data.frame(data_b)
	colnames(data_a)=col.names.dataA
	colnames(data_b)=col.names.dataB
        #create header for the matrix with common features
	if(is.na(time.thresh)==FALSE){
	mznames=c("index.A","mz.data.A", "time.data.A", "index.B","mz.data.B","time.data.B", "time.difference") 
        }else{
	mznames=c("index.A","mz.data.A", "index.B","mz.data.B") 
	}

        #Step 1 Group features by m/zdim(data_a)[1]
        mz_groups<-lapply(1:dim(data_a)[1],function(j){

                                commat={}
                                commzA=new("list")
                                commzB=new("list")
                                ppmb=(mz.thresh)*(data_a$mz[j]/1000000)

                                getbind_same<-which(abs(data_b$mz-data_a$mz[j])<=ppmb)

                                if(is.na(time.thresh)==FALSE){
                                  if(length(getbind_same)>0)
                                  {
                                          nearest_time_diff=10000
                                          bestmatch={}
					  rnames={}
					  temp={}
					  commat={}
                                          for (comindex in 1:length(getbind_same))
                                          {
											  tempA=cbind(j,data_a[j,c(1,2)])
											  tempB=cbind(getbind_same[comindex],data_b[getbind_same[comindex],c(1,2)])
					                                                  temp=cbind(tempA,tempB)
											
					                                                  timediff=abs(data_a[j,2]-data_b[getbind_same[comindex],2])
					                                                  
											  temp<-cbind(temp,timediff)
											 
					                                                  if(timediff<time.thresh && timediff<=nearest_time_diff)
					                                                  {
					                                                          bestmatch=as.data.frame(temp)
					                                                          nearest_time_diff=timediff
												  
												   
													
													 
												   
												   temp<-as.data.frame(temp)
												   commat<-rbind(commat,temp)
												    rnamestemp<-paste("mz",j,"_",comindex,sep="")
												    rnames<-c(rnames,rnamestemp)
												}
										}
					       
                                          
						#commat=as.data.frame(bestmatch)
					
						if(length(commat)>=4){
							rownames(commat)=rnames
					        }


                                  }
                                }
                                else
                                {
                                    if(length(getbind_same)>0)
                                    {
                                    	temp1<-{}
                                    for (comindex in 1:length(getbind_same))
                                          {
						  tempA=cbind(j,data_a[j,c(1)])
						  tempB=cbind(getbind_same[comindex],data_b[getbind_same[comindex],c(1)])
                                                  temp=cbind(tempA,tempB)
                                                  #temp=cbind(data_a[j,c(1)],data_b[getbind_same[comindex],c(1)])
                                                 temp1<-rbind(temp1,temp)
                                          }
					  commat=as.data.frame(temp1)
					  rnames<-paste("mz",j,"_",seq(1,length(getbind_same)),sep="")
					  rownames(commat)=rnames
					  
                                    }
                                }
                                return(as.data.frame(commat))


                })
		
	#Step 2 Sub-group features from Step 1 by Retention time
        #find the features with RT values within the defined range as compared to the query feature

        uniqueinA={}
        uniqueinB={}
        commat=data.frame()

	
	if(length(mz_groups)>0){
        for(j in 1:length(mz_groups))
        {
                temp_diff={}
		
		if(is.list(mz_groups)==TRUE)
		{
			tempdata=mz_groups[[j]]
		}
		else
		{
			tempdata=mz_groups[j]
		}
		
                if(length(tempdata)>1)
                {
			
			colnames(tempdata)=mznames
			tempdata=as.data.frame(t(tempdata))
			
                        temp=tempdata
                        
                        temp=as.data.frame(temp)
			
                        if(is.null(commat)==TRUE)
                        {
                                commat=t(temp)

                        }
                        else
                        {

                                commat=rbind(commat,t(temp))

                        }
			

                }
	

        }

        if(is.null(dim(commat))==FALSE)
        {
                commat=as.data.frame(commat)
        }
	}
	
	
        return(commat)
}
