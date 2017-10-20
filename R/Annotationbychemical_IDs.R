Annotationbychemical_IDs <-
function(dataA,queryadductlist=c("M+H"),adduct_type=c("S","Acetonitrile"),adduct_table,max.mz.diff=10,outloc,numnodes=10, otherdbs=FALSE,otherinfo=FALSE,keggCompMZ,syssleep=0.01){

	adduct_names<-as.character(adduct_table[,1])
	adductlist<-adduct_table[,4]
	mult_charge<-adduct_table[,3]
	num_mol<-adduct_table[,2]
	names(adductlist)<-as.character(adduct_names)
	names(mult_charge)<-as.character(adduct_names)
	names(num_mol)<-as.character(adduct_names)
	alladducts<-adduct_names


	#load("~/Documents/Emory/JonesLab/Projects/xMSannotator/keggCompMZ.Rda")
	
	
	#adduct_table<-adduct_table[which(adduct_table$Type%in%adduct_type),] #=="S" | adduct_table$Type=="Acetonitrile",]

#adduct_table<-adduct_table[which(adduct_table$Mode%in%adduct_mode),] #=="positive",]

#adduct_table<-adduct_table[which(adduct_table$Adduct%in%queryadductlist),] #adduct_table$Adduct=="M+H" | adduct_table$Adduct=="M+Na",]


suppressWarnings(dir.create(outloc))

setwd(outloc)
#keggres<-KEGG.annotation(dataA=mz_search_list,queryadductlist = c("positive"),xMSannotator.outloc)


#cur_fname<-"/Users/karanuppal/Documents/Emory/JonesLab/Projects/MaHPIC/Exp2/c18/apLCMS_with_xMSanalyzer_merged_data/apLCMS_feature_list_at_p1_U_p2cor0.7_CV100.txt"
#dataA<-read.table(cur_fname,sep="\t",header=TRUE)

mz_search_list_1<-as.data.frame(keggCompMZ[which(keggCompMZ$Adduct%in%adduct_table$Adduct),c(1,7)])

mz_search_list_1<-apply(mz_search_list_1,2,as.numeric)

#print(dim(adduct_table))
	cl<-makeSOCKcluster(numnodes)
		
		clusterEvalQ(cl, library(XML))
		clusterEvalQ(cl, library(R2HTML))
		clusterEvalQ(cl, library(RCurl))
		clusterEvalQ(cl, library(SSOAP))
		

			clusterEvalQ(cl, "processWSDL")
		clusterEvalQ(cl, library(png))
		clusterEvalQ(cl, "Annotationbychemical_IDschild")
		clusterEvalQ(cl, "keggCompMZ")


clusterExport(cl, "find.Overlapping.mzs")

 clusterExport(cl, "find.Overlapping.mzsvparallel")
     clusterExport(cl, "overlapmzchild")
     
#clusterEvalQ(cl, "chemCompMZ")		
	   clusterExport(cl, "getVenn")
	   		clusterEvalQ(cl, library(limma))
		#print(paste("Query adduct: ",adductname,sep=""))
		
		mz.annot.res<-new("list")
		
		mzlist<-dataA[,1]
		min_mz<-min(mzlist)
		max_mz<-max(mzlist)

		#mz_group<-ceiling(max_mz/min_mz)
		mz_group<-10
		#mz_group<-round(max_mz/10)
		#length(mzlist)
		num_mz<-1
		
		
		s1<-seq(1,length(queryadductlist))
		
		
		mz.annot.res<-parLapply(cl,s1,Annotationbychemical_IDschild,dataA,queryadductlist,adduct_type,adduct_table,max.mz.diff,outloc, otherdbs,otherinfo,keggCompMZ)
				if(is(mz.annot.res,"try-error")){
					Sys.sleep(10)
					mz.annot.res<-parLapply(cl,s1,Annotationbychemical_IDschild,dataA, queryadductlist,adduct_type,adduct_table,max.mz.diff,outloc, otherdbs,otherinfo,keggCompMZ)					
				}
		
		if(FALSE){
		for(mind in seq(1,length(mzlist),mz_group)){
		
			stopind<-mind+mz_group-1
			if(stopind>length(mzlist)){
			
			stopind<-length(mzlist)
			}
			#s1<-mzlist[mind:stopind]
			s1<-dataA[mind:stopind,]
			s1<-unique(s1)
			num_mz<-dim(s1)[1]
			if(num_mz%%50>0){
			Sys.sleep(syssleep)	
			}else{
			Sys.sleep(syssleep)	
			}
			
			if(num_mz>100){
			
				repeat{
					
				cur.annot.res<-parLapply(cl,s1,Annotationbychemical_IDschild,queryadductlist,adduct_type,adduct_table,max.mz.diff,outloc, otherdbs,otherinfo,keggCompMZ)
				if(is(cur.annot.res,"try-error")){
					Sys.sleep(10)
					cur.annot.res<-parLapply(cl,s1,Annotationbychemical_IDschild,queryadductlist,adduct_type,adduct_table,max.mz.diff,outloc, otherdbs,otherinfo,keggCompMZ)					
				}else{
				break
				}
	
				}			
			
				mz.annot.res<-c(mz.annot.res,cur.annot.res)
			
			}else{
				#for(i in 1:length(s1))
				{
					#print("mz is")
					#print(s1[i])
					rescur<-Annotationbychemical_IDschild(s1,queryadductlist,adduct_type,adduct_table,max.mz.diff,outloc, otherdbs,otherinfo,keggCompMZ)
					#print(length(rescur))
					if(length(rescur)>0){
						rescur<-as.matrix(rescur)
						#print(dim(rescur))
						if(dim(rescur)[2]==1){
							rescur<-t(rescur)
							rescur<-as.data.frame(rescur)
						} 
						rescur<-as.data.frame(rescur)
						#print(dim(rescur))
					mz.annot.res[[length(mz.annot.res)+1L]]<-rescur
					}
				}
			}
			if(mind%%10>0){
			Sys.sleep(syssleep/2)
			}else{
				Sys.sleep(1)
				stopCluster(cl)
				cl<-makeSOCKcluster(numnodes)
				
			
				clusterEvalQ(cl, library(XML))
				clusterEvalQ(cl, library(RCurl))
				clusterEvalQ(cl, library(SSOAP))
				clusterEvalQ(cl, library(png))
			clusterEvalQ(cl, "processWSDL")
				clusterEvalQ(cl, "Annotationbychemical_IDschild")
			}
			
				fname<-paste("cur_res.Rda",sep="")
				#save(mz.annot.res,file=fname)
		}
		}
		stopCluster(cl)
		
		res={}
		
		#print(length(mz.annot.res))
		if(length(mz.annot.res)>0){
		#print(adductname)
		for(mzl in 1:length(mz.annot.res))
		{
			#print(dim(mz.annot.res[[mzl]]))
			res=rbind(res,mz.annot.res[[mzl]])
			
		}
		}
		#print(mz.annot.res)
		res<-unique(res)
		
		#print("length of res")
		#print(length(res))

	return(res)

}
