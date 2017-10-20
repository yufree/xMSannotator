get_peak_blocks_graph <-
function(dataA,simmat,adjacencyfromsimilarity,time_step=3,max.rt.diff=5,outloc,column.rm.index=NA,cor.thresh=NA,cormethod="pearson",networktype = "unsigned",num_nodes=2){

#fname<-"/Users/karanuppal/Documents/Emory/JonesLab/Projects/NIST_QSTD_60K/apLCMS_with_xMSanalyzer_merged_data/apLCMS_feature_list_at_p1_U_p2cor0.7_CV100.txt"

#setwd("/Users/karanuppal/Documents/Emory/JonesLab/Projects/NIST_QSTD_60K/apLCMS_with_xMSanalyzer_merged_data/")

#fname<-"/Users/karanuppal/Documents/Emory/JonesLab/Projects/NIST_QSTD_60K/XCMS/xcmsCW_snthresh3step0.1mzdiff-0.001max50bw10ppm10.txt"

#setwd("/Users/karanuppal/Documents/Emory/JonesLab/Projects/NIST_QSTD_60K/XCMS/")

#setwd(outloc)

#time_step<-3

#dataA<-read.table(fname,sep="\t",header=TRUE)

#dataexpA_comp<-dataA

#1:3,6:11
if(is.na(column.rm.index)==FALSE){
dataA<-dataA[,-c(column.rm.index)]
}


cnames<-colnames(dataA)

cnames[1]<-"mz"
cnames[2]<-"time"

colnames(dataA)<-as.character(cnames)

data_m<-t(dataA)



feat_inf<-paste(dataA[,1],dataA[,2],sep="_")

sample.col.start<-2
if(FALSE)
{
if(cormethod=="pearson"){
ADJdataOne<-adjacency(data_m, 
                         type = "signed hybrid", 
                         power = 1,corOptions = "use = 'p'") 
                         }else{
                         	ADJdataOne<-adjacency(data_m, 
                         type = "signed hybrid", 
                         power = 1,corOptions = "use = 'p',method='spearman'")
                         	}
}
                         powers = c(c(1:10), seq(from = 12, to=20, by=2))
                         
                         
                         
                         if(adjacencyfromsimilarity==FALSE)
                         {
                             sft = pickSoftThreshold(data=data_m, dataIsExpr=TRUE,powerVector = powers, verbose = 0)
                             power_val=sft$powerEstimate
                             
                             if(is.na(power_val)==TRUE){
                                 power_val=6
                             }
                             
                             if(cormethod=="pearson"){
                                 ADJdataOne<-adjacency(datExpr=data_m,
                                 type = networktype,
                                 power = power_val,corOptions = "use = 'p'")
                             }else{
                                 ADJdataOne<-adjacency(datExpr=data_m,
                                 type = networktype, 
                                 power = power_val,corOptions = "use = 'p',method='spearman'")
                             }      
                         }else{
                             
                             sft = sft = pickSoftThreshold.fromSimilarity(similarity=simmat, powerVector = powers, verbose = 0)
                             power_val=sft$powerEstimate
                             
                             if(is.na(power_val)==TRUE){
                                 power_val=6
                             }
                             
                             ADJdataOne<-adjacency.fromSimilarity(similarity=simmat,power=power_val,type=networktype)
                         }
                        
                         
                         
#cormat<-cor2pcor(cormat)

if(is.na(cor.thresh)==FALSE){
ADJdataOne[ ADJdataOne<cor.thresh ] <- 0
ADJdataOne[ ADJdataOne>=cor.thresh ] <- 1
}

#rownames(ADJdataOne) <- sample(letters, nrow(ADJdataOne))
#colnames(ADJdataOne) <- seq(ncol(ADJdataOne))

g1<-graph.adjacency(ADJdataOne,weighted=TRUE,mode="undirected")
#summary(g1)

s1<-E(g1)

e1<-get.edgelist(g1)

l1<-levels(as.factor(e1[,1]))

set1<-new("list")

for(l2 in 1:length(l1)){
	
	#set1[l2]<-e1[which(e1[,1]==l1[l2]),2]
	
	set1[[l2]]<-(e1[which(e1[,1]==l1[l2]),2])
	
	set1[[l2]]<-as.numeric(gsub(set1[[l2]],pattern="V",replacement=""))
}


cor_groups<-set1
cor_group_list<-{}

cor_groups_clust<-new("list")
		#
	  for(m in 1:length(cor_groups)){
	  
	  #cor_groups_clust[[m]]<-unique(cor_groups[[m]])
	  if((m%in%cor_group_list)==FALSE)
	 {
		for(n in (m+1):length(cor_groups)){
		
			if(n>length(cor_groups)){
				break;
			}
			com1<-intersect(cor_groups[[m]],cor_groups[[n]])
			if(length(com1)>0){
				
				cor_groups[[m]]<-c(cor_groups[[m]],cor_groups[[n]])
				
				cor_group_list<-c(cor_group_list,n)
				cor_groups[[n]]<-c(0)
				
				
			}
			
		}
		
			cor_groups[[m]]<-unique(cor_groups[[m]])
		
		}
	  }
	  
	  if(length(cor_group_list)>0){
		cor_groups<-cor_groups[-cor_group_list]
		}
	 
	  diffmat<-{}
	#cor_groups<-cor_groups_clust
	levelAnum<-1
	diffmat<-lapply(1:length(cor_groups),function(m){
	 if((m%in%cor_group_list)==FALSE)
	 	if(cor_groups[[m]]!=0)
	 	{
	 	levelAnum<-m
	 	
	 	cur_group_ind<-unique(cor_groups[[m]])
	 	
	 	cur_group_ind<-as.numeric(gsub(cur_group_ind,pattern="V",replacement=""))
	 	cur_group<-dataA[cur_group_ind,]
	 	cur_group<-cbind(levelAnum,cur_group)
	 	
	 	#diffmat<-rbind(diffmat,cur_group)
	 	return(cur_group)
	 	#levelAnum<-levelAnum+1
	 	}
	 	})
	 	
diffmat<-do.call(rbind,diffmat)


diffmat<-as.data.frame(diffmat)
mod_list<-diffmat$levelAnum

t1<-table(mod_list)
mod_names<-names(t1)
#mod_names<-as.numeric(mod_names)


 
time_mult_fact<-1


 

time_mult_fact<-1


 
 diffmatC<-{}

dataA<-cbind(data_mzrt,dataA)
dataA<-as.data.frame(dataA)

 d1<-density(dataA$time,bw="nrd",from=min(dataA$time),to=(10+max(dataA$time,na.rm=TRUE)))
 
# time_step<-3
 time_step<-abs(d1$x[1]-d1$x[time_step+1])
 
 t1<-d1$x
 
 time_step<-1*time_step

diffmatB<-{}
for(i in 1:length(mod_names)){

#groupB_res<-sapply(1:length(mod_names),function(i){
	
	groupA_num<-mod_names[i]
	
	subdata<-dataA[which(moduleLabels==groupA_num),]
	#subdata<-dataA[which(m1==groupA_num),]

 subdata<-subdata[order(subdata$time),]

 groupB<-group_by_rt(subdata,time_step,max.rt.diff=max.rt.diff,groupnum=groupA_num)
	 
	  diffmatB<-rbind(diffmatB,groupB)
	  
	  }
	  
	 return(diffmatB)
}
