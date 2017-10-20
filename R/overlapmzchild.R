overlapmzchild <-
function(j,mz.thresh=10,time.thresh=NA,data_a,data_b){
    
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
            #for (comindex in 1:length(getbind_same))
            #chemscoremat_conf_levels<-foreach(c=1:length(chemids), .combine=rbind) %dopar%
            commat<-lapply(1:length(getbind_same),function(comindex)
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
                    #  commat<-rbind(commat,temp)
                    #rnamestemp<-paste("mz",j,"_",comindex,sep="")
                    #rnames<-c(rnames,rnamestemp)
                }
             return(temp)
            })
					     
                         commat<-ldply(commat,rbind)
                           
                           #commat=as.data.frame(bestmatch)
                           
                           #  if(length(commat)>=4){
                               #    rownames(commat)=rnames
                               #}
                           
                           
        }
    }
    else
    {
        if(length(getbind_same)>0)
        {
            temp1<-{}
            #for (comindex in 1:length(getbind_same))
            commat<-lapply(1:length(getbind_same),function(comindex)
            {
                tempA=cbind(j,data_a[j,c(1)])
                tempB=cbind(getbind_same[comindex],data_b[getbind_same[comindex],c(1)])
                temp=cbind(tempA,tempB)
               temp<-as.data.frame(temp)
		#temp=cbind(j,data_a[j,c(1)],getbind_same[comindex],data_b[getbind_same[comindex],c(1)])
                #temp1<-rbind(temp1,temp)
                return(temp)
            })
	    
            commat<-ldply(commat,rbind)
            commat=as.data.frame(commat)
            #rnames<-paste("mz",j,"_",seq(1,length(getbind_same)),sep="")
            #rownames(commat)=rnames
            
        }
    }
    return(as.data.frame(commat))
    


}
