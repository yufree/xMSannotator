group_by_rt_density <-
function(subdata,time_step=1,max.rt.diff=5,groupnum){
	
	diffmatB<-{}
	d1<-density(subdata$time,bw="nrd",from=min(subdata$time),to=(0.01+max(subdata$time)))
	
	time_step<-(d1$x[2]-d1$x[1])*time_step
		
	time_cor_groups<-sapply(list(myData1=subdata),function(x)  split(x,cut(subdata$time,breaks=d1$x)))
	
	#for(gnum in 1:length(time_cor_groups)){
	diffmatB<-lapply(1:length(time_cor_groups),function(gnum){
		curmat<-time_cor_groups[[gnum]]
		
		if(dim(curmat)[1]>0){
		
		cur_mod_clust<-paste(curmat$Module_RTclust,gnum,sep="_")
		
		curmat$Module_RTclust<-cur_mod_clust
		#diffmatB<-rbind(diffmatB,curmat)
		}
		return(curmat)
	})
	
	
	diffmatB<-ldply(diffmatB,rbind)
	
	#diffmatB[which(diffmatB$mz>221 & diffmatB$mz<255),1:3]
	return(diffmatB)

}
