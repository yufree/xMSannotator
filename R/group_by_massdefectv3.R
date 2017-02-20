group_by_massdefectv3 <-
function(dataA,pct.mz.diff=10,mdgnum=-1){
diffmatB<-{}	
 
 

time_mult_fact<-1
diff_mz_num<-1

t1<-dataA$time

mzdefect<-dataA$mz-floor(dataA$mz)


d1<-density(mzdefect,bw="nrd",from=min(mzdefect),to=(0.01+max(mzdefect)))
 
 time_step<-2
 time_step<-abs(d1$x[1]-d1$x[time_step+1])
 
 t1<-d1$x
 
 #time_step<-step.mz.diff*time_step

	 mz_groups<-sapply(1:dim(dataA)[1],function(j)
 #levelB_groups<-sapply(1:length(t1),function(j)
 {
                                commat={}
                                diffmz=new("list")
                                
                                #cur_group<-diffmat[which(diffmat$levelAnum==group_labels_level1[j]),]
                                
                                #b1<-hist(cur_group$time,breaks=seq(min(cur_group$time)-max.rt.diff,max(cur_group$time)+max.rt.diff,2*max.rt.diff))
                                
                                #for(i in 1:dim(cur_group)[1])
                                {
                                #getbind_same<-which(abs(subdata$time-subdata$time[j])<=2*time_step)
                                
                                
                                getbind_same<-which((100*abs(mzdefect-mzdefect[j])/(mzdefect[j]))<=(pct.mz.diff))
                                
                                diffmz[[diff_mz_num]]=getbind_same #dataA[getbind_same,]
                                diff_mz_num=diff_mz_num+1
                                }
                               
                                return(diffmz)
  })
  
  
		     
		     
	mz_groups<-unique(mz_groups)

  del_list<-{}
	   
	   gid<-{}
	   
	   gid<-paste("MDgroup",dim(dataA)[1],sep="")
  
  levelB_groups<-mz_groups
	  
	# if(FALSE)
	 {
	   del_list<-{}
	   levelB_groups_clean<-mz_groups
	   #length(levelB_groups)
	  
	  	
	  	for(m in 1:length(levelB_groups)){
	  
	  	if(length(levelB_groups[[m]])>0){
	   if((m%in%del_list)==FALSE){
		for(n in (m+1):length(levelB_groups)){
		
			if(n>length(levelB_groups)){
				break;
			}
			
			com1<-intersect(levelB_groups[[m]],levelB_groups[[n]])
			if(length(com1)>0){
			
				levelB_groups[[m]]<-c(levelB_groups[[m]],levelB_groups[[n]])
				del_list<-c(del_list,n)
				#levelB_groups[[n]]<-{}
			}
		}
		levelB_groups[[m]]<-unique(levelB_groups[[m]])
		#levelB_groups[[m]][[1]]<-unique(levelB_groups[[m]][[1]])
		}else{
				del_list<-c(del_list,m)
			}
		}
		
		
		
		if(length(del_list)>0){
		 levelB_groups_clean<-levelB_groups[-del_list]
		}
		}
		time_cor_groups<-unique(levelB_groups_clean)
	}
	  
	groupnum<-1
	  
	#time_cor_groups<-unique(levelB_groups)
	
	diffmatB<-{}
	for(gnum in 1:length(time_cor_groups)){
		
		cur_group<-{}
		
	 	if(length(time_cor_groups[[gnum]])>0){
	 	cur_group_ind<-unique(time_cor_groups[[gnum]])
	 	
	 	cur_group<-dataA[cur_group_ind,]
	 	
	 	MDgroup<-paste("MDgroup",mdgnum,gnum,sep="_")
	 	
	 	cur_group<-cbind(MDgroup,cur_group)
	 	
	 	
	 	#print(cur_group)
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
	 cnames[1]<-"MDGroup"
	 colnames(diffmatB)<-cnames

	return(diffmatB) #$MDgroup)
}
