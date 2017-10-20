get_peak_blocks_modulesvhclust <-
function(dataA=NA,simmat=NA,adjacencyfromsimilarity=FALSE,time_step=3,max.rt.diff=10,outloc,column.rm.index=NA,cor.thresh=NA,deepsplit=2,
minclustsize=20,cutheight=0.2,cormethod="spearman",networktype = "unsigned",num_nodes=2,step1log2scale=TRUE,mycl_metabs=NA)
{


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

colnames(data_m)<-feat_inf


if(step1log2scale==TRUE){

		data_m<-2^(data_m)
}


sample.col.start<-2

 #which(dataA$mz>145.9 & dataA$mz<146)
 #cor.test(t(as.numeric(dataA[303,-c(1:2)])),t(as.numeric(dataA[774,-c(1:2)])))

powers = c(c(1:10), seq(from = 12, to=20, by=2))


   
do_stepwise=0 
    #print(dim(data_m))
    
    powers = c(c(1:10), seq(from = 12, to=20, by=2))
    sft = try(pickSoftThreshold(data=data_m, dataIsExpr=TRUE,powerVector = powers, verbose = 0),silent=TRUE)

    
   if(is(sft,"try-error")){

	power_val=6

   }else{
	power_val=sft$powerEstimate
    
    if(is.na(power_val)==TRUE){
        power_val=6
    }
    
   }
   
   #if(FALSE)
   {
   
   
   
   netclassA =try(blockwiseModules(datExpr=(data_m), checkMissingData = FALSE, blocks = mycl_metabs, maxBlockSize = 5000,
   blockSizePenaltyPower = 100, randomSeed = 12345, loadTOM = FALSE,
   corType = "pearson", maxPOutliers = 1, quickCor = 0, pearsonFallback = "individual",
   cosineCorrelation = FALSE, power = power_val, networkType = "unsigned",
   TOMType = "signed", TOMDenom = "min", getTOMs = NULL, saveTOMs = FALSE,
   saveTOMFileBase = "blockwiseTOM", deepSplit = deepsplit, detectCutHeight = NULL,
   minModuleSize = minclustsize, maxCoreScatter = NULL,
   minGap = NULL, maxAbsCoreScatter = NULL, minAbsGap = NULL,
   minSplitHeight = NULL, minAbsSplitHeight = NULL, useBranchEigennodeDissim = FALSE,
   minBranchEigennodeDissim = mergeCutHeight, stabilityLabels = NULL,
   minStabilityDissim = NULL, pamStage = TRUE, pamRespectsDendro = FALSE,
   reassignThreshold = 1e-06, minCoreKME = 0.5, minCoreKMESize = minclustsize/3,
   minKMEtoStay = 0.3, mergeCutHeight = cutheight, impute = TRUE,
   trapErrors = FALSE, numericLabels = FALSE, nThreads = num_nodes,
   verbose = 0, indent = 0),silent=TRUE)
   
   
   #save(netclassA,file="netclassA.Rda")
   #do_stepwise=1
    if(is(netclassA,"try-error")){ 
   
		do_stepwise=1
    }else{ 
    n1<-unlist(netclassA$colors)
    n2<-as.data.frame(n1)
    mod_list<-as.numeric(n2[,1])
 
    Alldegrees1<-softConnectivity(datExpr=data_m,power=power_val,minNSamples=2)
    #Alldegrees1<-
    Alldegrees1<-cbind(Alldegrees1,Alldegrees1,Alldegrees1,Alldegrees1)
    }
}



if(do_stepwise==1)
{ 

if(adjacencyfromsimilarity==FALSE)
{

	powers = c(c(1:10), seq(from = 12, to=20, by=2))
    sft = try(pickSoftThreshold(data=data_m, dataIsExpr=TRUE,powerVector = powers, verbose = 0),silent=TRUE)



	if(is(sft,"try-error")){
	power_val=6
	}else{
		power_val=sft$powerEstimate
	}

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

	sft = try(pickSoftThreshold.fromSimilarity(similarity=simmat, powerVector = powers, verbose = 0),silent=TRUE)
	#power_val=sft$powerEstimate

	if(is(sft,"try-error")){
        	power_val=6
        }else{
                power_val=sft$powerEstimate
        }

	if(is.na(power_val)==TRUE){
		power_val=6
	}


ADJdataOne<-adjacency.fromSimilarity(similarity=simmat,power=power_val,type=networktype)
}



rnames_simmat<-rownames(ADJdataOne)

dup_names<-which(duplicated(rnames_simmat)==TRUE)

if(length(dup_names)>0){

ADJdataOne<-ADJdataOne[-c(dup_names),-c(dup_names)]
}

#dissTOMCormat=TOMdist(cormat,TOMType="signed")

dissTOMCormat=TOMdist(ADJdataOne)

#save topological overlap dissimilarity
#save(dissTOMCormat,file="TOMdist.Rda")

#save topological overlap similarity
simTOMCormat=1-dissTOMCormat

#save(simTOMCormat,file="TOMsim.Rda")

#rm(simTOMCormat)

#if(FALSE)
{
hierTOMCormat = flashClust(as.dist(dissTOMCormat),method="complete");
par(mfrow=c(1,1))
pdf("plot.pdf")
plot(hierTOMCormat,labels=F,main="Dendrogram")
dev.off()

#save(list=ls(),file="hier_analysis.Rda")

#if(FALSE)
{


colorhdataOne2=cutreeDynamic(hierTOMCormat,distM= dissTOMCormat,deepSplit=deepsplit, minClusterSize=minclustsize, pamRespectsDendro = FALSE, pamStage=TRUE)

#colorhdataOne2=cutreeHybrid(hierTOMCormat,distM= dissTOMCormat,deepSplit=deepsplit, minClusterSize=minclustsize, pamRespectsDendro = FALSE,pamStage=TRUE)

l2colors<-levels(as.factor(colorhdataOne2)) #$labels))

#print(length(l2colors))
if(length(l2colors)>1){

m2=try(mergeCloseModules(data_m,colors=colorhdataOne2,cutHeight=cutheight),silent=TRUE)

if(is(m2,"try-error")){

	#m2<-doWGCNA(data_m,deepsplit=deepsplit,minclustsize=minclustsize,cutheight=cutheight)
	mod_list<-colorhdataOne2 #as.numeric(m2$colors)
}else{

	mod_list<-as.numeric(m2$colors)
}
#m2=mergeCloseModules(data_m,colors=colorhdataOne2$labels,cutHeight=cutheight)

#print(length(table(m2$colors)))

}else{
	m2=colorhdataOne2
	mod_list<-m2 #as.numeric(m2)
	#mod_list<-rep(0,dim(dataA)[1])
	#m2<-doWGCNA(data_m,deepsplit=deepsplit,minclustsize=minclustsize,cutheight=cutheight)
	#mod_list<-as.numeric(m2$colors)
}
#max.rt.diff<-5


 #mod_list<-as.numeric(m2$colors)
}

}

#print(length(mod_list))
#print(dim(ADJdataOne))

Alldegrees1=intramodularConnectivity(ADJdataOne,mod_list)

#rownames(Alldegrees1)<-as.character(feat_inf)
#print("mod_list")
#print(mod_list[1:4])

#print("here")

save(mod_list,file="mod_list.Rda")
}



t1<-table(mod_list)
mod_names<-names(t1)
mod_names<-as.numeric(mod_names)

time_mult_fact<-1

 
 diffmatC<-{}
rm(data_m)


dataA<-cbind(data_mzrt,dataA)
dataA<-as.data.frame(dataA)

diffmatB<-{}
diffmatB<-lapply(1:length(mod_names),function(i){

	groupA_num<-mod_names[i]
	
	subdata<-dataA[which(mod_list==groupA_num),]
	subdata<-subdata[order(subdata$time),]

 	groupB<-group_by_rt_histv1(subdata,time_step,max_diff_rt=max.rt.diff,groupnum=groupA_num)

	rownames(groupB)<-NULL	 
	
	  
	return(groupB)
	  })



    diffmatB<-ldply(diffmatB,rbind)
   # save(diffmatB,file="diffmatB.Rda")
    
	rm(dataA)
 #   print(dim(diffmatB))
#rownames(diffmatB)<-as.character(feat_inf)
	diffmatB<-cbind(Alldegrees1[,c(1:4)],diffmatB)
	
	rm(dataA)
	rm(dissTOMCormat)
	rm(ADJdataOne)

	 return(diffmatB)
}
