group_metabs_prev <-
function(dataA, max.rt.diff=10, alignment.tool="apLCMS", clust.method="correlation",cor.method="pearson",corthresh=0.7,corpvaluethresh=0.05,mult.test.cor=FALSE)
{	
        diff_mz_num=1
        unique_mz={}
        ppm_v={}
        rt_v={}
        cnames=colnames(dataA)
        dataA<-as.data.frame(dataA)
	
		#dataA<-as.data.frame(dataA,2,as.numeric)
         if(alignment.tool=="apLCMS")
        {
              sample.col.start=6
        }
        else
        {
              if(alignment.tool=="XCMS")
              {
                    sample.col.start=9
                    cnames[1]="mz"

                    cnames[4]="time"
                    colnames(dataA)=cnames
                    
              }
              else
              {
              		 sample.col.start=3
                    cnames[1]="mz"

                    cnames[2]="time"
                    colnames(dataA)=cnames
                    #stop(paste("Invalid value for alignment.tool. Please use either \"apLCMS\" or \"XCMS\"", sep=""))
              }
        }

	
	if(clust.method=="correlation"){
				 cl<-makeCluster(10)

        #clusterExport(cl, "getCorchild")

        #goodfeats[mnum,-c(1:2)]),data_mt

        #system.time(pearson_res<-parRapply(cl,dataA[,-c(1:sample.col.start)],getCorchild,t(dataA[,-c(1:sample.col.start)]),cor.method))


        #stopCluster(cl)
		 pearson_resmat<-{}
		 
		 
		 #return(list(cormat=cormat,complete_pearsonpvalue_mat=complete_pearsonpvalue_mat,complete_pearsonqvalue_mat=complete_pearsonqvalue_mat))
	
			pearson_Res_all<-{}

			cormat<-WGCNA::cor(t(dataA[,-c(1:sample.col.start)]), use = 'p')
		    cormat<-as.data.frame(cormat)
        cormat<-as.matrix(cormat)
        
           
      		cor_list<-{}
		
		#p1<-cor2pcor(cormat)

		cor_groups<-new("list")
		
		#for(cm in 1:dim(cormat)[1]){
			
			cor_groups<-lapply(1:dim(cormat)[1],function(cm){
			corind_list=new("list")
			cor_ind<-which(cormat[,cm]>=corthresh)
			corind_list[[diff_mz_num]]=cor_ind #dataA[getbind_same,]
            diff_mz_num=diff_mz_num+1
			return(corind_list)
			})
		
		cor_group_list<-{}
		#
	  for(m in 1:length(cor_groups)){
	  
	   if((m%in%cor_group_list)==FALSE){
		for(n in (m+1):length(cor_groups)){
		
			if(n>length(cor_groups)){
				break;
			}
			com1<-intersect(cor_groups[[m]][[1]],cor_groups[[n]][[1]])
			if(length(com1)>0){
				
				cor_groups[[m]][[1]]<-c(cor_groups[[m]][[1]],cor_groups[[n]][[1]])
				cor_groups[[m]][[1]]<-unique(cor_groups[[m]][[1]])
				cor_group_list<-c(cor_group_list,n)
				cor_groups[[n]][[1]]<-c(0)
			}
		}
			cor_groups[[m]][[1]]<-unique(cor_groups[[m]][[1]])
		
		}
	  }
	  
	  if(length(cor_group_list)>0){
		cor_groups<-cor_groups[-cor_group_list]
		}
	 
	  diffmat<-{}
	cor_groups<-unique(cor_groups)
	levelAnum<-1
	 for(m in 1:length(cor_groups)){
	 	if((m%in%cor_group_list)==FALSE){
	 	
	 	
	 	cur_group_ind<-unique(cor_groups[[m]][[1]])
	 	
	 	cur_group<-dataA[cur_group_ind,]
	 	cur_group<-cbind(levelAnum,cur_group)
	 	
	 	diffmat<-rbind(diffmat,cur_group)
	 	levelAnum<-levelAnum+1
	 	}
	 	}
	

	#cor.test(as.numeric(cur_group[1,-c(1:10)]),as.numeric(cur_group[10,-c(1:10)]))

	}else{
		
		if(clust.method=="wgcna"){
			
				levelAnum<-mlabels$colors
				diffmat<-cbind(levelAnum,dataA)
				diffmat<-diffmat[order(diffmat$levelAnum),]
				
				for(m1 in mlabels){
					subdata<-diffmat[which(mlabels==m1),]
					
				}
				
		}
		
	}
	#b1<-hist(diffmat$time,breaks=seq(0,max(diffmat$time),max.rt.diff))
	
	
	diffmat<-as.data.frame(diffmat)
	
	
	
	diffmat<-diffmat[order(diffmat$time),]
	
	group_labels_level1<-levels(as.factor(diffmat$levelAnum))
	
	mz_groups<-new("list")
	
	#group_labels_level1<-names()
  #Step 1 Group features by m/z
        mz_groups<-lapply(1:length(group_labels_level1),function(j){
                                commat={}
                                diffmz=new("list")
                                
                                cur_group<-diffmat[which(diffmat$levelAnum==group_labels_level1[j]),]
                                
                                #b1<-hist(cur_group$time,breaks=seq(min(cur_group$time)-max.rt.diff,max(cur_group$time)+max.rt.diff,2*max.rt.diff))
                                
                                for(i in 1:dim(cur_group)[1]){
                                getbind_same<-which(abs(diffmat$time-cur_group$time[i])<=max.rt.diff & diffmat$levelAnum==group_labels_level1[j])
                                
                                diffmz[[diff_mz_num]]=getbind_same #dataA[getbind_same,]
                                diff_mz_num=diff_mz_num+1
                                }
                                return(diffmz)
          })
	  
	   del_list<-{}
	   #length(mz_groups)
	   for(k in 1:length(mz_groups)){
	  for(m in 1:length(mz_groups[[k]])){
	  
	   if((m%in%del_list)==FALSE){
		for(n in (m+1):length(mz_groups[[k]])){
		
			if(n>length(mz_groups[[k]])){
				break;
			}
			com1<-intersect(mz_groups[[k]][[m]],mz_groups[[k]][[n]])
			if(length(com1)>0){
			
				mz_groups[[k]][[m]]<-c(mz_groups[[k]][[m]],mz_groups[[k]][[n]])
				del_list<-c(del_list,n)
				mz_groups[[k]][[n]]<-c(0)
			}
		}
		
		#mz_groups[[m]][[1]]<-unique(mz_groups[[m]][[1]])
		}
		}
		if(length(del_list)>0){
		 mz_groups[[k]]<-mz_groups[[k]][-del_list]
		}
	  }
	  
	 
	  diffmatB<-{}
	  groupnum<-1
	time_cor_groups<-unique(mz_groups)
	
	for(m in 1:length(time_cor_groups)){
		
		for(n in 1:length(time_cor_groups[[m]])){
	 	levelBnum<-groupnum
	 	
	 	cur_group_ind<-unique(time_cor_groups[[m]][[n]])
	 	
	 	cur_group<-diffmat[cur_group_ind,]
	 	cur_group<-cbind(levelBnum,cur_group)
	 	if(length(cur_group)>0){
	 	diffmatB<-rbind(diffmatB,cur_group)
	 	groupnum<-groupnum+1
	 	}
	 	}
	 }
	
	  	
	diff_mz_num=1
	mz_groups<-lapply(1:length(mz_groups),function(j){
                                commat={}
                                diffmz=new("list")
                                getbind_same=mz_groups[[j]][[1]]
                                #curdata<-cbind(j,dataA[getbind_same,])
                                diffmz[[diff_mz_num]]=dataA[getbind_same,]
                                diff_mz_num=diff_mz_num+1
                                return(diffmz)
          })

	if(mult.test.cor==TRUE){
		mz_group_size<-lapply(1:length(mz_groups),function(j){
		num_rows<-dim(mz_groups[[j]][[1]])
		n=num_rows[[1]][[1]]
		num_comp=(n*(n-1))/2
		})

	
		num_comparisons<-sum(unlist(mz_group_size))
	}else{
		num_comparisons=1
	}

	
	
       
                               												
			        return(diffmat)
}
