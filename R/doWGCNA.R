doWGCNA <-
function(data_m,deepsplit=4,minclustsize=1,cutheight=0.05){
	
	#library(WGCNA)
#fname<-"/Users/karanuppal/Documents/Emory/JonesLab/Projects/NIST_QSTD_60K/apLCMS_with_xMSanalyzer_merged_data/apLCMS_feature_list_at_p1_U_p2cor0.7_CV100.txt"
#setwd("/Users/karanuppal/Documents/Emory/JonesLab/Projects/NIST_QSTD_60K/apLCMS_with_xMSanalyzer_merged_data/")
#dataexprA<-read.table(fname,sep="\t",header=TRUE)

#data_m<-t(dataexprA[,-c(1:9)])
#feat_inf<-paste(dataexprA[,1],dataexprA[,2],sep="_")

ADJdataOne=adjacency(datExpr= data_m, power=1) #, power = beta1)
dissTOMdataOne=TOMdist(ADJdataOne)
hierTOMdataOne = hclust(as.dist(dissTOMdataOne),method="complete");
par(mfrow=c(1,1))
pdf("plot.pdf")
plot(hierTOMdataOne,labels=F,main="Dendrogram in MB20,MB50")



#save(list=ls(),file="hier_analysis.Rda")

#colorhdataOne=cutreeDynamic(hierTOMdataOne,distM= dissTOMdataOne,deepSplit=2, pamRespectsDendro = FALSE)

#m1=mergeCloseModules(data_m,colors=colorhdataOne)

#length(table(m1$colors))

colorhdataOne2=cutreeDynamic(hierTOMdataOne,distM= dissTOMdataOne,deepSplit=deepsplit, minClusterSize=minclustsize, pamRespectsDendro = FALSE)

m2=mergeCloseModules(data_m,colors=colorhdataOne2)
length(table(m2$colors))

return(m2)
	
}
