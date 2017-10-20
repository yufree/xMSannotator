get_peak_blocks_hclust <-
function(dataA,time_step=1,max.rt.diff=5,outloc,column.rm.index=NA,cor.thresh=NA,deepsplit=4,minclustsize=2,cutheight=0.2,cormethod="spearman"){

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



if(cormethod=="pearson"){
ADJdataOne<-adjacency(data_m, 
                         type = "signed hybrid", 
                         power = 1,corOptions = "use = 'p'") 
                         }else{
                         	ADJdataOne<-adjacency(data_m, 
                         type = "signed hybrid", 
                         power = 1,corOptions = "use = 'p',method='spearman'")
                         	}      
#dissTOMCormat=TOMdist(cormat,TOMType="signed")

dissTOMCormat=TOMdist(ADJdataOne)


#if(FALSE)
{
hierTOMCormat = hclust(as.dist(dissTOMCormat),method="complete");
par(mfrow=c(1,1))
pdf("plot.pdf")
plot(hierTOMCormat,labels=F,main="Dendrogram")

#m1 <- cutree(hierTOMCormat, h=max(hierTOMCormat$height)/2)


m1=cutreeDynamic(hierTOMCormat,distM=dissTOMCormat,deepSplit=4, minClusterSize=1, pamRespectsDendro = TRUE)

#print(m1)

mod_list<-m1

#save(list=ls(),file="hier_analysis.Rda")

#g1<-graph.adjacency(ADJdataOne,weighted=TRUE)
#summary(g1)

if(FALSE)
{
colorhdataOne=cutreeDynamic(hierTOMCormat,distM= dissTOMCormat,deepSplit=deepsplit, cutHeight=cutheight,pamRespectsDendro = TRUE)

l1colors<-levels(as.factor(colorhdataOne))

if(length(l1colors)>100){
	#print(l1colors)
m1=mergeCloseModules(data_m,colors=colorhdataOne)
}
else{
	
	m1=colorhdataOne
	}

length(table(m1$colors))
}

if(FALSE){
#colorhdataOne2=cutreeDynamic(hierTOMCormat,distM= dissTOMCormat,deepSplit=1, minClusterSize=2, pamRespectsDendro = FALSE)

colorhdataOne2=cutreeHybrid(hierTOMCormat,distM= dissTOMCormat,deepSplit=deepsplit, minClusterSize=minclustsize, pamRespectsDendro = TRUE,pamStage=TRUE)

l2colors<-levels(as.factor(colorhdataOne2$labels))

#print(length(l2colors))
if(length(l2colors)>1){
#m2=mergeCloseModules(data_m,colors=colorhdataOne2,cutHeight=0.05)

m2=mergeCloseModules(data_m,colors=colorhdataOne2$labels,cutHeight=cutheight)

#print(length(table(m2$colors)))

}else{
	m2=colorhdataOne2

}
#max.rt.diff<-5


 mod_list<-as.numeric(m2$colors)
}

}

if(FALSE){
#hr <- hclust(as.dist(1-cor(t(good_metabs[,-c(1:2)]),use="pairwise.complete.obs"))) #metabolites
		hc <- hclust(as.dist(dissTOMCormat)) #,use="pairwise.complete.obs") #samples
#h73<-heatmap.2(as.matrix(good_metabs[,-c(1:2)]), Rowv=NULL, Colv=as.dendrogram(hc),  col=topo.colors(256), scale="none",key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=0.5, cexCol=0.85,xlab="",ylab="", main="",ColSideColors=patientcolors)

#mycl_samples <- cutree(hc, h=max(hc$height)/2)
#mycl_metabs <- cutree(hr, h=max(hr$height)/2)


mycl_metabs=cutreeDynamic(hc,distM=dissTOMCormat,deepSplit=4, minClusterSize=5, pamRespectsDendro = TRUE)
}

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
	
	#subdata<-dataA[which(m2$colors==groupA_num),]
	subdata<-dataA[which(m1==groupA_num),]

 subdata<-subdata[order(subdata$time),]

 groupB<-group_by_rt(subdata,time_step,max.rt.diff=max.rt.diff,groupnum=groupA_num)
	 
	  diffmatB<-rbind(diffmatB,groupB)
	  
	  }
	  
	 return(diffmatB)
}
