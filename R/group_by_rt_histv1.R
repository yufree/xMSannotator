group_by_rt_histv1 <-
function(mchemicaldata,time_step=1,max_diff_rt=10,groupnum){
    
    
    range1<-range(mchemicaldata$time)
    diff1<-abs(range1[1]-range1[2])
    
    if(nrow(mchemicaldata)>1 && diff1>1){
        d1<-try(density(mchemicaldata$time,bw="nrd",from=min(mchemicaldata$time),to=(max_diff_rt+max(mchemicaldata$time,na.rm=TRUE))),silent=TRUE)
        if (is(d1, "try-error")){
            bw_i<-max_diff_rt
        }else{
            
            
            bw_i<-d1$bw
        }
        
    }else{
        bw_i<-max_diff_rt
        
    }
    
    
    
    #	if(diff1>max_diff_rt)
    {
        if(bw_i>6*max_diff_rt){
            
            hist_inc=1*max_diff_rt #min(d1$bw,max_diff_rt)
        }else{
            
            hist_inc<-min(bw_i,max_diff_rt)
        }
        
        
        # 	h1<-hist(mchemicaldata$time,breaks=seq(min(mchemicaldata$time),to=(max_diff_rt),by=hist_inc))
        
        
        h1<-hist(mchemicaldata$time,breaks=seq(min(mchemicaldata$time)-max_diff_rt,to=(max_diff_rt+max(mchemicaldata$time)),by=hist_inc))
        
        
        
        
        t1<-cbind(h1$mids,h1$density)
        t1<-as.data.frame(t1)
        colnames(t1)<-c("mids","density")
        
        
        
        
        # time_cor_groups<-sapply(list(myData1=h1$density),function(x)  split(x,cut(h1$mids,breaks=seq(min(h1$mids)-max_diff_rt,max(h1$mids)+max_diff_rt,10))))
        clusternum<-1;clusterlabels={};
        t1prev=min(mchemicaldata$time)
        for(i in 1:dim(t1)[1]){
            
            if(t1[i,2]==0 && (t1[i,1]-t1prev)>max_diff_rt){
                clusternum<-clusternum+1
                
            }
            
            if(t1[i,2]>0){
                clusterlabels<-c(clusterlabels,clusternum)
            }else{
                clusterlabels<-c(clusterlabels,0)
            }
            
        }
        t2<-cbind(t1,clusterlabels)
        t2<-as.data.frame(t2)
        if(length(which(t2[,3]==0))>0){
            t2<-t2[-which(t2[,3]==0),]
        }
        t2[1,1]<-min(mchemicaldata$time)
        clusterlabels2<-{}
        clusterlabels3<-lapply(1:length(mchemicaldata$time),function(k){
            cl1<-t2$clusterlabels[which(abs(t2$mids-mchemicaldata$time[k])==min(abs(t2$mids-mchemicaldata$time[k])))[1]]
            return(cl1)
        })
        clusterlabels3<-unlist(clusterlabels3)
        
        clusterlabels3<-replace(clusterlabels3,which(is.na(clusterlabels3)==TRUE),0)
        
    }
    #else{
    
    #	clusterlabels3<-rep(1,dim(mchemicaldata)[1])
    #}
    levelBnum<-paste(groupnum,clusterlabels3,sep="_")
    
    diffmatB<-cbind(levelBnum,mchemicaldata)
    
    
    
    diffmatB<-as.data.frame(diffmatB)
    
    
    cnames<-colnames(diffmatB)
    cnames[1]<-"Module_RTclust"
    colnames(diffmatB)<-cnames
    
    return(diffmatB)
    
    
}
