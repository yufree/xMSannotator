group_by_isotopes_vsplit <-
function(dataA,step.mz.diff=1,mdgnum=-1){
diffmatB<-{}	
 
 

time_mult_fact<-1
diff_mz_num<-1

t1<-dataA$time

mzdefect<-1000*((dataA$mz-floor(dataA$mz)))


d1<-density(mzdefect,bw="nrd",from=min(mzdefect),to=(0.01+max(mzdefect)))
 
 time_step<-step.mz.diff*1
 time_step<-abs(d1$x[1]-d1$x[time_step+1])
 
 t1<-d1$x
 
 	   
	   gid<-{}
	   
	   gid<-paste("ISgroup",dim(dataA)[1],sep="")
  
 
  
#time_cor_groups<-sapply(list(myData1=dataA),function(x)  #print(x);split(x,cut(mzdefect,breaks=d1$x)))
 	
#time_cor_groups<-sapply(list(myData1=dataA),function(x)  split(x,cut(mzdefect,breaks=seq(0,1000,1))))
 
 time_cor_groups<-sapply(list(myData1=dataA),function(x)  split(x,cut(mzdefect,breaks=seq(0,1000,time_step))))	
 	
 	#length(time_cor_groups)
	diffmatB<-{}
	for(gnum in 1:length(time_cor_groups)){
		
		cur_group<-{}
		
	 	if(length(time_cor_groups[[gnum]])>0){
	 		 	
	 	ISgroup<-paste("ISgroup",mdgnum,gnum,sep="_")
	 	
	 	cur_group<-as.data.frame(time_cor_groups[[gnum]])
	 	cur_group<-cbind(ISgroup,cur_group)
	 	
	 	
	 	if(length(cur_group)>0){
	 		cur_group<-as.data.frame(cur_group)
	 		cur_group<-cur_group[order(cur_group$mz),]
	 	diffmatB<-rbind(diffmatB,cur_group)
	 	#groupnum<-groupnum+1
	 	}
	 	}
	 }
	 
	 diffmatB<-as.data.frame(diffmatB)
	 
	 
#table(diffmatB$MDgroup)
	 
	 
	 cnames<-colnames(diffmatB)
	 cnames[1]<-"ISGroup"
	 colnames(diffmatB)<-cnames

	return(diffmatB) #$MDgroup)
}
