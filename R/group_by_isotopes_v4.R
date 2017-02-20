group_by_isotopes_v4 <-
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
 
 time_cor_groups<-sapply(list(myData1=dataA),function(x)  split(x,cut(mzdefect,breaks=seq(0,1000,d1$x))))	
 	
 	#length(time_cor_groups)
	diffmatB<-{}
	for(gnum in 1:length(time_cor_groups)){
		
		cur_group<-{}
		
		#print(dim(time_cor_groups[[gnum]]))
	 	if(dim(time_cor_groups[[gnum]])[1]>0){
	 		 	
	 	ISgroup<-paste("ISgroup",mdgnum,gnum,sep="_")
	 	
	 	#print(time_cor_groups[[gnum]])
	 	cur_group<-as.data.frame(time_cor_groups[[gnum]])
	 	cur_group<-cbind(ISgroup,cur_group)
	 	
	 	#print("cur group")
	 	
	 	#print(cur_group)
	 	
	 	if(length(cur_group)>0){
	 		#cur_group<-as.data.frame(cur_group)
	 		cur_group<-cur_group[order(cur_group$mz),]
	 		
	 		#print(length(diffmatB))
	 		if(length(diffmatB)<1){
	 			diffmatB<-cur_group
	 			}else{
	 	
	 	#print("here")
	 	#print(dim(diffmatB))
	 	#print(dim(cur_group))
	 	diffmatB<-rbind(diffmatB,cur_group)
	 	}
	 	#groupnum<-groupnum+1
	 	}
	 	}
	 }
	 
	 diffmatB<-as.data.frame(diffmatB)
	 
	 
#table(diffmatB$MDgroup)
	 
	 if(dim(diffmatB)[1]>0){
	 cnames<-colnames(diffmatB)
	 cnames[1]<-"ISGroup"
	 #print(cnames) 
	 colnames(diffmatB)<-cnames
	}

	return(diffmatB) #$MDgroup)
}
