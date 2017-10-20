Metlin.annotation <-
function(dataA,max.mz.diff=10, queryadductlist=c("M+H"), xMSannotator.outloc, tokenstr=NA)
{
	data_a<-as.data.frame(dataA)
	 
	#print("Using the 1st column as \"mz\" for annotation.")
		    
	if(is.na(tokenstr)==TRUE){
	
	stop("Please specify a valid METLIN security token. \nCreate a METLIN account to obtain a token: \n http://metlin.scripps.edu/soap/register.php")
	
 	}
	mzlist<-data_a[,1]
	
	#=c("PubChem","MassBank", "EPA DSSTox","EPA Toxcast","NIST Chemistry WebBook","KEGG",
	#"Human Metabolome Database", "ChEMBL", "ChEBI", "NIAID","Pesticide Common Names","SMPDB Small Molecule Pathway Database",
	#"MeSH","LipidMAPS","ChemBank","BioCyc")
	
	dir.create(xMSannotator.outloc,showWarnings=FALSE)
	setwd(xMSannotator.outloc)
	adductlist=c(1.00727,22.989171,38.963171,-35.012729,-17.002729,0.0227,7.01597,18.033871,
33.033471,42.033871,44.971171,64.015771,1.00727,11.998247,22.989171,1.00727,
8.334590,15.6618,-19.01839,-1.00727,18.998371,20.974671,34.969371,36.948571,
44.998171,59.013871,78.918885,-1.00727,-1.00727)
alladducts<-c("M+H","M+Na","M+K","M+H-2H2O","M+H-H2O", "M-H2O+NH4", "M+Li","M+NH4",
"M+CH3OH+H","M+ACN+H","M+2Na-H","M+ACN+Na","M+2H", "M+H+Na","M+2Na","M+3H",
"M+2H+Na","M+2Na+H","M-H2O-H", "M-H", "M+F","M+Na-2H","M+Cl","M+K-2H",
"M+FA-H","M+CH3COO","M+Br","M-2H","M-3H")
names(adductlist)<-c("M+H","M+Na","M+K","M+H-2H2O","M+H-H2O", "M-H2O+NH4", "M+Li","M+NH4",
"M+CH3OH+H","M+ACN+H","M+2Na-H","M+ACN+Na","M+2H", "M+H+Na","M+2Na","M+3H",
"M+2H+Na","M+2Na+H","M-H2O-H", "M-H", "M+F","M+Na-2H","M+Cl","M+K-2H",
"M+FA-H","M+CH3COO","M+Br","M-2H","M-3H")

mult_charge<-c(1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,3,3,3,1,1,1,1,1,1,1,1,1,2,3)

names(mult_charge)<-c("M+H","M+Na","M+K","M+H-2H2O","M+H-H2O", "M-H2O+NH4", "M+Li","M+NH4",
"M+CH3OH+H","M+ACN+H","M+2Na-H","M+ACN+Na","M+2H", "M+H+Na","M+2Na","M+3H",
"M+2H+Na","M+2Na+H","M-H2O-H", "M-H", "M+F","M+Na-2H","M+Cl","M+K-2H",
"M+FA-H","M+CH3COO","M+Br","M-2H","M-3H")
	
	
	
	if(queryadductlist[1]=="positive")
	{
		queryadductlist<-c("M+H","M+Na","M+K","M+H-2H2O","M+H-H2O", "M+Li","M+NH4",
		"M+CH3OH+H","M+ACN+H","M+2Na-H","M+ACN+Na","M+2H", "M+H+Na","M+2Na","M+3H",
		"M+2H+Na","M+2Na+H")
	}else{
		if(queryadductlist[1]=="negative")
	{
		queryadductlist<-c("M-H2O-H", "M-H", "M+F","M+Na-2H","M+Cl","M+K-2H",
		"M+FA-H","M+CH3COO","M-2H","M-3H")
	}else{
		if(queryadductlist[1]=="all"){
		
			
		queryadductlist<-alladducts
		
		
		}else{
			if(length(which(queryadductlist%in%alladducts==FALSE))>0){
			
				errormsg<-paste("Adduct should be one of:",sep="")
				for(i in alladducts){errormsg<-paste(errormsg, i,sep=" ; ")}
				stop(errormsg, "\n\nUsage: feat.batch.annotation.Metlin(dataA,max.mz.diff=10, queryadductlist=c(\"M+H\", \"M+Na\"), xMSannotator.outloc, numnodes=1)", 
				"\n\n OR feat.batch.annotation.Metlin(dataA,max.mz.diff=10, queryadductlist=c(\"positive\"), xMSannotator.outloc, numnodes=1)",
				"\n\n OR  feat.batch.annotation.Metlin(dataA,max.mz.diff=10, queryadductlist=c(\"negative\"), xMSannotator.outloc, numnodes=1)",
				"\n\n OR  feat.batch.annotation.Metlin(dataA,max.mz.diff=10, queryadductlist=c(\"all\"), xMSannotator.outloc, numnodes=1)"
				)
			}
		
		}
	}	
		}
	parentres={}

	res<-{}
	cnames<-c("","Adduct","Query.m/z", "MetlinID", "Mass", "Name", "Molecular.Formula")
		
		
	for(adnum in 1:length(queryadductlist))
	{
		adductname=queryadductlist[adnum]
		#adductmass=adductlist[as.character(adductname)]
		adductname<-gsub(x=adductname,pattern="\\+",replacement="%2B",perl=TRUE)
		
		#tokenstr<-"A30GCsqF5VkBoHlH"
		
		
		
		mz_seq<-seq(1,length(mzlist),500)
		temp_res<-{}
			
		for(mzind in 1:length(mzlist)){
		
		url<-paste("http://metlin.scripps.edu/REST/search/index.php?token=",tokenstr,sep="")
		suburl<-{}
			suburl<-paste(suburl,"&mass[]=",mzlist[mzind],sep="")
		
		
		
		url<-paste(url,suburl,"&adduct[]=",adductname,"&tolunits=ppm&tolerance=",max.mz.diff,sep="")
		
		#print(url)
		#195.0877&mass[]=181.0702&adduct[]=M%2BH&tolunits=ppm&tolerance=30",sep="")

		r1<-getURL(url)

		

	

		parser<-newJSONParser()
		parser$addData(r1)
		#metlinres<-parser$getObject()
		metlinres<-try(parser$getObject(),silent=TRUE)
		#print(length(metlinres))
		if (is(metlinres, "try-error")){
			cur_res<-{}
		}else{
		for(m1 in 1:length(metlinres)){
			
			metlin_res_m1<-metlinres[[m1]]
			#print(metlin_res_m1)
			for(m2 in 1:length(metlin_res_m1)){
				cur_res<-metlin_res_m1[[m2]]
				#print(length(cur_res))
				if(length(cur_res)>1){
				html_link=paste("<a href=http://metlin.scripps.edu/metabo_info.php?molid=",cur_res$molid,">",cur_res$molid,"</a>",sep="")
				temp_res<-rbind(temp_res,c(mzlist[mzind],html_link,cur_res$molid,cur_res$mass,cur_res$name,cur_res$formula))
				}
			}
	
			}
		}
		}
		
		res<-unique(temp_res)
		
		text_res<-{}
		
		if(length(res)>0){

			adductname<-gsub(x=adductname,pattern="%2B",replacement="\\+",perl=TRUE)
			adductname=rep(adductname,dim(res)[1])
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
				if(dim(temp_res)[2]==1){
					temp_res<-t(temp_res)
				}
			}
			#temp_res<-as.data.frame(temp_res)
			colnames(temp_res)=NULL
			text_resindex<-c(1,2,4:7)
			text_res<-temp_res[,c(text_resindex)]
			text_res<-as.matrix(text_res)
			text_res<-as.data.frame(text_res)
			sernum=seq(1,dim(text_res)[1])
			text_res<-cbind(sernum,text_res)
			colnames(text_res)=cnames
			parentres=rbind(parentres,temp_res)
			#rm(temp_res)
			colnames(parentres)=NULL
			
		
		fname=paste(xMSannotator.outloc,"/Metlin_annotation_results_",queryadductlist[adnum],".txt",sep="")
		write.table(text_res,file=fname,sep="\t",row.names=FALSE)
		}
		Sys.sleep(2)
	
	}
	
	html_res<-{}
	text_res<-{}
	
		
	if(length(parentres)>0){
	
		res<-parentres[order(parentres[,2]),]
		res<-as.matrix(res)
	
		if(dim(res)[2]==1){res<-t(res)}
		html_resindex<-c(1,2,3,5:7)
	
		html_res<-res[,c(html_resindex)]
	
		html_res<-as.matrix(html_res)
		if(dim(html_res)[2]==1){html_res<-t(html_res)}
		
		#sernum=seq(1,dim(html_res)[1])
		
		#html_res<-cbind(sernum,html_res)
		html_res<-as.data.frame(html_res)
		
			cnames<-c("Adduct","Query.m/z", "MetlinID", "Mass", "Name", "Molecular.Formula")
		colnames(html_res)<-cnames

		fname=paste("Metlin_annotation_results",sep="")
		unlink(fname)
	
		HTMLInitFile(filename=fname,Title="Metlin annotation results", outdir=xMSannotator.outloc)
		fname=paste(xMSannotator.outloc,"/Metlin_annotation_results.html",sep="")
		HTML(html_res,file=fname,Border=1,innerBorder=1,useCSS=TRUE)
		HTMLEndFile(file=fname)
	
			text_resindex<-c(1,2,4:7)
		
			text_res<-res[,c(text_resindex)]
			text_res<-as.matrix(text_res)

			text_res<-as.data.frame(text_res)
			
			sernum=seq(1,dim(text_res)[1])
			text_res<-cbind(sernum,text_res)
				cnames<-c("","Adduct","Query.m/z", "MetlinID", "Mass", "Name", "Molecular.Formula")
			colnames(text_res)=cnames
			
			fname=paste(xMSannotator.outloc,"/Metlin_annotation_results_alladducts.txt",sep="")
			write.table(text_res,file=fname,sep="\t",row.names=FALSE)
		}
	return(list("text.res"=text_res,"html.res"=html_res))
}
