get_peak_blocks_modules <-
function(dataA,simmat,adjacencyfromsimilarity=FALSE,time_step=3,max.rt.diff=10,outloc,column.rm.index=NA,cor.thresh=NA,
deepsplit=2,minclustsize=30,cutheight=0.2,networktype="unsigned",num_nodes=2){

#fname<-"/Users/karanuppal/Documents/Emory/JonesLab/Projects/NIST_QSTD_60K/apLCMS_with_xMSanalyzer_merged_data/apLCMS_feature_list_at_p1_U_p2cor0.7_CV100.txt"

#setwd("/Users/karanuppal/Documents/Emory/JonesLab/Projects/NIST_QSTD_60K/apLCMS_with_xMSanalyzer_merged_data/")

#fname<-"/Users/karanuppal/Documents/Emory/JonesLab/Projects/NIST_QSTD_60K/XCMS/xcmsCW_snthresh3step0.1mzdiff-0.001max50bw10ppm10.txt"

#setwd("/Users/karanuppal/Documents/Emory/JonesLab/Projects/NIST_QSTD_60K/XCMS/")

#setwd(outloc)

#time_step<-3

#dataA<-read.table(fname,sep="\t",header=TRUE)

#dataexpA_comp<-dataA


cnames<-colnames(dataA)

cnames[1]<-"mz"
cnames[2]<-"time"

colnames(dataA)<-as.character(cnames)

data_mzrt<-dataA[,c(1:2)]


#1:3,6:11
if(is.na(column.rm.index)==FALSE){
dataA<-dataA[,-c(column.rm.index)]
}


feat_inf<-paste(dataA[,1],dataA[,2],sep="_")
dataA<-dataA[,-c(1:2)]

data_m<-t(dataA)

allowWGCNAThreads()
multiExpr = vector(mode = "list", length = 1)
multiExpr[[1]] = list(data = as.data.frame(data_m));

powers = c(c(1:10), seq(from = 12, to=20, by=2))

sft = pickSoftThreshold(data=data_m, dataIsExpr=TRUE,powerVector = powers, verbose = 0)
power_val=sft$powerEstimate

net = blockwiseModules(
datExpr=data_m, power = power_val, minModuleSize = minclustsize, deepSplit = deepsplit,
pamRespectsDendro = FALSE,
numericLabels = TRUE,saveTOMs = TRUE, verbose = 0,
corType = "pearson",networkType = networktype,nThreads=num_nodes)




consMEs = net$multiMEs;
moduleLabels = net$colors;
# Convert the numeric labels to color labels
moduleColors = labels2colors(moduleLabels)
consTree = net$dendrograms[[1]];

save(consMEs, moduleLabels, moduleColors, consTree, file = "Consensus-NetworkConstruction-auto_min30.RData")

mod_list<-moduleLabels
t1<-table(mod_list)
mod_names<-names(t1)
mod_names<-as.numeric(mod_names)


 
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
