KEGG.annotationvold <-
function(dataA,max.mz.diff=10, queryadductlist=c("M+H"), xMSannotator.outloc, numnodes=1,syssleep=1,adduct_table)
{
	data_a<-as.data.frame(dataA)
	 
	#print("Using the 1st column as \"mz\" for annotation.")
		    
	mzlist<-data_a[,1]
	
	dir.create(xMSannotator.outloc,showWarnings=FALSE)
	setwd(xMSannotator.outloc)
	
	#adduct_table<-read.table("/Users/karanuppal/Documents/Emory/JonesLab/Projects/xMSannotator/adduct_table.txt",sep="\t",header=TRUE)
	#adduct_table<-adduct_table[c(which(adduct_table[,6]=="S"),which(adduct_table[,6]=="Acetonitrile")),]
	
	adduct_names<-as.character(adduct_table[,1])
	adductlist<-adduct_table[,4]
	mult_charge<-adduct_table[,3]
	num_mol<-adduct_table[,2]
	names(adductlist)<-as.character(adduct_names)
	names(mult_charge)<-as.character(adduct_names)
	names(num_mol)<-as.character(adduct_names)
	alladducts<-adduct_names
	
	if(queryadductlist[1]=="positive")
	{
			queryadductlist<-adduct_names[which(adduct_table[,5]=="positive")]
		
	}else{
		if(queryadductlist[1]=="negative")
		{
			
			queryadductlist<-adduct_names[which(adduct_table[,5]=="negative")]
			
		}else{
		if(queryadductlist[1]=="all"){
		
			
		queryadductlist<-alladducts
		
		
		}else{
			if(length(which(queryadductlist%in%alladducts==FALSE))>0){
			
				errormsg<-paste("Adduct should be one of:",sep="")
				for(i in alladducts){errormsg<-paste(errormsg, i,sep=" ; ")}
				stop(errormsg, "\n\nUsage: feat.batch.annotation.KEGG(dataA,max.mz.diff=10, queryadductlist=c(\"M+H\", \"M+Na\"), xMSannotator.outloc, numnodes=1)", 
				"\n\n OR  feat.batch.annotation.KEGG(dataA,max.mz.diff=10, queryadductlist=c(\"positive\"), xMSannotator.outloc, numnodes=1)",
				"\n\n OR  feat.batch.annotation.KEGG(dataA,max.mz.diff=10, queryadductlist=c(\"negative\"), xMSannotator.outloc, numnodes=1)",
				"\n\n OR  feat.batch.annotation.KEGG(dataA,max.mz.diff=10, queryadductlist=c(\"all\"), xMSannotator.outloc, numnodes=1)"
				)
			}
		
		}
	}	
		}
		
	parentres={}
	#cnames<-c("","Adduct","Query.m/z", "Search mass \n tolerance range (+/-)", "Metlin", "Compound.Name", "CASID", "KEGG.Compound.ID", "KEGG.Pathway.name", "KEGG.Pathway.ID", "HMDB.ID", "PubChem.Substance.ID", "PubChem.Compound.ID", "ChEBI.ID", "LIPID.MAPS.ID")
  cnames<-c("","Adduct","Query.m/z", "Search mass \n tolerance range (+/-)", "Metlin", "Compound.Name",  "Chemical.Formula", "Exact Mass", "CASID", "KEGG.Compound.ID",
  "KEGG.Pathway.ID", "KEGG.Pathway.name","HMDB.ID", "PubChem.Substance.ID", "ChEBI.ID", "LIPID.MAPS.ID","KEGG.Brite.Category",
  "KEGG.Disease.ID","KEGG.Drug.ID","KEGG.Environ.ID")
 
	
	for(adnum in 1:length(queryadductlist))
	{
		adductname=queryadductlist[adnum]
		adductmass=adductlist[as.character(adductname)]
		
		cl<-makeSOCKcluster(numnodes)
		
		clusterEvalQ(cl, library(XML))
		clusterEvalQ(cl, library(RCurl))
		clusterEvalQ(cl, "feat.batch.annotation.child")
		
		#print(paste("Query adduct: ",adductname,sep=""))
		
		mz.annot.res<-new("list")
		min_mz<-min(mzlist)
		max_mz<-max(mzlist)

		mz_group<-ceiling(max_mz/min_mz)
		
		#mz_group<-round(max_mz/10)
		#length(mzlist)
		num_mz<-1
		
		for(mind in seq(1,length(mzlist),mz_group)){
		
			stopind<-mind+mz_group
			if(stopind>length(mzlist)){
			
			stopind<-length(mzlist)
			}
			s1<-mzlist[mind:stopind]
			s1<-unique(s1)
			num_mz<-num_mz+length(s1)
			
			if(num_mz%%50>0){
			Sys.sleep((syssleep/2))	
			}else{
			Sys.sleep(syssleep)	
			}
			if(length(s1)>1){
				
				mz.annot.res<-c(mz.annot.res,parLapply(cl,s1,feat.batch.annotation.child,max.mz.diff=max.mz.diff,adductname=adductname,adduct_table=adduct_table,syssleep=syssleep))
			
			}else{
				for(i in 1:length(s1)){
					rescur<-feat.batch.annotation.child(mz.val=s1[i],max.mz.diff=max.mz.diff,adductname=adductname,adduct_table=adduct_table,syssleep=syssleep)
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
					mz.annot.res<-c(mz.annot.res,rescur)
					}
				}
			}
			if(mind%%10>0){
				Sys.sleep((syssleep/10))			
			}else{
				Sys.sleep(syssleep)
				stopCluster(cl)
				cl<-makeSOCKcluster(numnodes)
				
			
				clusterEvalQ(cl, library(XML))
				clusterEvalQ(cl, library(RCurl))
				clusterEvalQ(cl, "feat.batch.annotation.child")
			}
			
		}
		
		stopCluster(cl)
		res={}
		#print(adductname)
			if(length(mz.annot.res)>0){
		for(mzl in 1:length(mz.annot.res))
		{
			res=rbind(res,mz.annot.res[[mzl]])
			
		}
		}
		res<-unique(res)
		
		text_res<-{}
		
		if(length(res)>0){
			
		adductname=c(rep(adductname,dim(res)[1]))
		
		temp_res<-cbind(adductname,res)
		
		temp_res<-as.matrix(temp_res)
		
		
		
		if(dim(temp_res)[2]==1){
		
			temp_res<-t(temp_res)
			temp_res<-as.data.frame(temp_res)
		
		}
		
		bad_rows<-which(temp_res[,2]=="1")
		
		if(length(bad_rows)>0){
		temp_res<-temp_res[-bad_rows,]
		temp_res<-as.matrix(temp_res)
		
		#temp_res<-t(temp_res)
			if(dim(temp_res)[2]==1){
		
			temp_res<-t(temp_res)
			
		
		}

		
		}
		#temp_res<-as.data.frame(temp_res)
		colnames(temp_res)=NULL
		
			
		#text_resindex<-c(1,2,5,6,7,8,11,10,13,15,17,19,21)
		#text_resindex<-c(1,2,5,6,7,8,9,12,11,14,16,18,20,22)
		
		#text_resindex<-c(1,2,5,6:8,9,11,13,14,16,18,20,21)
		#text_resindex<-c(1,2,5,6:7,4,8,10,12,13,17,19,21)
		
		text_resindex<-c(1,2,5,6:7,4,8,9,11,13,14,16,18,20,22,26:28)
		text_resindex<-text_resindex+1
		#print(dim(temp_res))
		text_res<-temp_res[,c(1,text_resindex)]
		
		text_res<-as.matrix(text_res)
		
		if(dim(text_res)[2]==1){
		
			text_res<-t(text_res)
			
		
		}
		text_res<-as.data.frame(text_res)
		bad_rows<-which(text_res[,2]=="1")
		
		
		if(length(bad_rows)>0){
		text_res<-text_res[-bad_rows,]
		text_res<-as.matrix(text_res)
		text_res<-t(text_res)
		}
		text_res<-as.data.frame(text_res)
		
		sernum=seq(1,dim(text_res)[1])
		text_res<-cbind(sernum,text_res)
		colnames(text_res)=cnames
		text_res<-text_res[,-c(5)]
		
		parentres=rbind(parentres,temp_res)
		rm(temp_res)
		colnames(parentres)=NULL
		
		}
		#num_cols<-dim(text_res)[2]
		#text_res<-cbind(text_res[,c(1:10)],text_res[,c(num_cols)],text_res[,c(11:(num_cols-1))])
		fname=paste(xMSannotator.outloc,"/KEGG_annotation_results_",queryadductlist[adnum],".txt",sep="")
		write.table(text_res,file=fname,sep="\t",row.names=FALSE)
		
		
		
		Sys.sleep(syssleep)
	
	}
	
	html_res<-{}
	text_res<-{}
	
		
	if(length(parentres)>0){
		
		res<-parentres[order(parentres[,2]),]

	#html_resindex<-c(1,2,5,4,6:7,9,11:12,14,16,18,20,22)
	res<-as.matrix(res)
	
	if(dim(res)[2]==1){res<-t(res)}
	
	
	
	#html_resindex<-c(1,2,5,6:7,9,11:12,14,16,18,20,22)
	
	#html_resindex<-c(1,2,5,6:7,4,8,10,12,13,15,17,19,21)
	html_resindex<-c(1,2,5,6:7,4,8,10,12,13,15,17,19,21,22,23:25)
	
	html_resindex<-html_resindex+1
	html_res<-res[,c(1,html_resindex)]
		html_res<-as.matrix(html_res)
		
	
		if(dim(html_res)[2]==1){html_res<-t(html_res)}
		
		sernum=seq(1,dim(html_res)[1])
		
		#html_res<-cbind(sernum,html_res)
		html_res<-as.data.frame(html_res)
		
		cnames<-c("Adduct","Query.m/z", "Search mass \n tolerance range (+/-)","Metlin", "Compound.Name",  "Chemical.Formula",  "Exact Mass", "CASID", "KEGG.Compound.ID",
		"KEGG.Pathway.ID", "KEGG.Pathway.Name", "HMDB.ID", "PubChem.Substance.ID", "ChEBI.ID", "LIPID.MAPS.ID","KEGG.Brite.Category",
		"KEGG.Disease.ID","KEGG.Drug.ID","KEGG.Environ.ID")
 
		#"KEGG.Gene.ID",
		
		colnames(html_res)<-cnames
		
		html_res<-html_res[,-c(4)]
		#fname=paste("Annotation_results",sep="")
		#num_cols<-dim(html_res)[2]
		#text_res<-cbind(html_res[,c(1:10)],html_res[,c(num_cols)],html_res[,c(11:(num_cols-1))])
		
		fname=paste("KEGG_annotation_results",sep="")
		unlink(fname)
		#fname=paste(xMSannotator.outloc,"/KEGG_annotation_results.html",sep="")
		HTMLInitFile(filename=fname,Title="KEGG annotation results", outdir=xMSannotator.outloc)
		fname=paste(xMSannotator.outloc,"/KEGG_annotation_results.html",sep="")
		HTML(html_res,file=fname,Border=1,innerBorder=1,useCSS=TRUE)
		HTMLEndFile(file=fname)
		
		
		
		cnames<-c("","Adduct","Query.m/z", "Search mass \n tolerance range (+/-)","Metlin", "Compound.Name",  "Chemical.Formula",  "Exact Mass", "CASID", "KEGG.Compound.ID",
		"KEGG.Pathway.ID", "KEGG.Pathway.Name", "HMDB.ID", "PubChem.Substance.ID", "ChEBI.ID", "LIPID.MAPS.ID","KEGG.Brite.Category",
		"KEGG.Disease.ID","KEGG.Drug.ID","KEGG.Environ.ID")
		
		#text_resindex<-c(1,2,5,6:7,4,8,9,11,13,14,16,18,20,22) #c(1,2,5,6,7,8,9,12,11,14,16,18,20,22)
		
		#text_resindex<-c(1,2,5,6:7,4,8,10,12,13,17,19,21)
		
		#text_resindex<-c(1,2,5,6:7,4,8,9,11,13,14,16,18,20,22,29:33)
		text_resindex<-c(1,2,5,6:7,4,8,9,11,13,14,16,18,20,22,26:28)
		text_resindex<-text_resindex+1
		text_res<-res[,c(1,text_resindex)]
		
		text_res<-as.matrix(text_res)
		if(dim(text_res)[2]==1){text_res<-t(text_res)}
		
		if(length(text_res)>0){
		sernum=seq(1,dim(text_res)[1])
		}else{
			sernum={}
			}
		text_res<-cbind(sernum,text_res)
		
		text_res<-as.data.frame(text_res)
		colnames(text_res)=cnames
		text_res<-text_res[,-c(5)]
		#num_cols<-dim(text_res)[2]
		#text_res<-cbind(text_res[,c(1:10)],text_res[,c(num_cols)],text_res[,c(11:(num_cols-1))])
		fname=paste(xMSannotator.outloc,"/KEGG_annotation_results_alladducts.txt",sep="")
		write.table(text_res,file=fname,sep="\t",row.names=FALSE)
		}
	return(list("text.res"=text_res,"html.res"=html_res))
}
