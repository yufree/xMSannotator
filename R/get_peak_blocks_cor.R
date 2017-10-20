get_peak_blocks_cor <-
function(dataA, max.rt.diff=10, alignment.tool="apLCMS", clust.method="correlation",cor.method="pearson",corthresh=0.7,corpvaluethresh=0.05,time_step=3,mult.test.cor=FALSE)
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
              sample.col.start=2
        }
        else
        {
              if(alignment.tool=="XCMS")
              {
                    sample.col.start=2
                    cnames[1]="mz"

                    cnames[2]="time"
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
		
		cormat<-cor2pcor(cormat)

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
	
	
	mod_list<-diffmat$levelAnum
	
	t1<-table(mod_list)
mod_names<-names(t1)
mod_names<-as.numeric(mod_names)


 
time_mult_fact<-1


 
 diffmatC<-{}

#t1<-seq(from=1,to=max(dataA$time),by=time_step/2)
 #t1<-hist(subdata$time,breaks=seq(from=1,to=max(subdata$time),by=10/2))
 
# sigma<-min(sd(subdata$time),(quantile(subdata$time)[4]-quantile(subdata$time)[2])/1.34)
 #s1<-min(sqrt(var(x)), h)
 
 #bw1<-1.06 * sigma * length(subdata$time)^(-1/5)
 
 d1<-density(dataA$time,bw="nrd",from=min(dataA$time),to=(10+max(dataA$time)))
 
 
 time_step<-abs(d1$x[1]-d1$x[time_step+1])
 
 t1<-d1$x
 
 time_step<-1*time_step

diffmatB<-{}
for(i in 1:length(mod_names)){

#groupB_res<-sapply(1:length(mod_names),function(i){
	
	groupA_num<-mod_names[i]
	
	subdata<-dataA[which(m2$colors==groupA_num),]
	
 subdata<-subdata[order(subdata$time),]

 groupB<-group_by_rt(subdata,time_step,max.rt.diff=max.rt.diff,groupnum=groupA_num)
	 
	  diffmatB<-rbind(diffmatB,groupB)
	  
	  }
	  
	 return(diffmatB)	
	
       
                               												
			        #return(diffmat)
}
