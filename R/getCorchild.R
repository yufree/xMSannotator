getCorchild <-
function(cur_mzdata,data_mt,cor.method){
	
	
	
	pearson_res<-lapply(1:dim(data_mt)[2],function(j){
		 return(WGCNA::cor(as.numeric(cur_mzdata),data_mt[,j],method=cor.method,use="pairwise.complete.obs"))
		
	
	})
	
	pvalues_list<-{}
	pearson_resmat<-{}
	pearson_list<-{}
	
	for(i in 1:length(pearson_res))
	{
		pearson_list<-c(pearson_list,pearson_res[[i]][1])
		nval<-length(which(is.na(data_mt[,i])==FALSE))
		
		pvalues_list<-c(pvalues_list,docortest(nval,pearson_res[[i]][[1]]))
	}
	
		
	return(list(cormat=pearson_list,complete_pearsonpvalue_mat=pvalues_list)) #,complete_pearsonqvalue_mat=qvalues_list))
}
