get_mz_by_KEGGdrugIDs <-
function(keggIDs,queryadductlist=c("M+H"),max.mz.diff=10,syssleep=0.01,adduct_table=NA){
	cnames<-c("mz","chemical_ID", "Name", "Formula", "MonoisotopicMass", "Adduct","AdductMass")
	
	#rm(adduct_table)
	data(adduct_table)
	#print(adduct_table)
	
	#adduct_table<-read.table("/Users/karanuppal/Documents/Emory/JonesLab/Projects/xMSannotator/adduct_table.txt",sep="\t",header=TRUE)
	#adduct_table<-adduct_table[c(which(adduct_table[,6]=="S"),which(adduct_table[,6]=="Acetonitrile")),]
	if(is.na(adduct_table)==TRUE){
		data(adduct_table)
	}else{
	
		#print(adduct_table)
	
	}
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
		
	map_res<-{}
	for(i in 1:length(keggIDs)){
		
		Sys.sleep(syssleep)
		
		keggID<-keggIDs[i]
		#print(keggID)
		keggID<-paste(keggID,collapse="")
		kegglink<-paste("<a href=http://www.genome.jp/dbget-bin/www_bget?dr:",keggID,">",keggID,"</a>",sep="")
				
				
				
		search_link=paste("http://rest.kegg.jp/get/dr:",as.character(keggID),sep="")
		
		d2<-try(read.delim(search_link,header=FALSE),silent=TRUE)

	if(is(d2,"try-error")){
		next
	}else{
		
			dm<-keggGet(paste("dr:",keggID,sep=""))
			
			formula<-dm[[1]]$FORMULA
			chem_name<-paste(dm[[1]]$NAME,collapse="")
			
			exact_mass<-as.numeric(dm[[1]]$EXACT_MASS)
			
			if(length(exact_mass)>0){
			for(adnum in 1:length(queryadductlist))
			{
				adductname=queryadductlist[adnum]
				adductmass=adductlist[as.character(adductname)]
                                adductcharge=mult_charge[as.character(adductname)]
								adductnmol=num_mol[as.character(adductname)]
								
                                 				
				  #mono_mass=((mz.val*adduct_charge)-(adductmass))/(adductnmol)
				 
				 mz=((exact_mass*adductnmol)+(adductmass))/adductcharge
				 
				 
				delta_ppm=(max.mz.diff)*(mz/1000000)
				min_mz=round((mz-delta_ppm),5)
				max_mz=round((mz+delta_ppm),5)
				
				
				res={} 
				mzorig=round(exact_mass,5)
				delta_ppm=round(delta_ppm,5)
				
				syssleep1<-(syssleep/5)
				Sys.sleep(syssleep1)
				
				
				cur_map_res<-c(mz,keggID,chem_name,formula,mzorig,adductname,adductmass)
				map_res<-rbind(map_res,cur_map_res)
				rownames(map_res)<-NULL
		
			}
			}else{
			cur_map_res<-c("-",keggID,chem_name,"-","-","-","-")
			map_res<-rbind(map_res,cur_map_res)
			}
	
		}
			
	
			
	
	#colnames(map_res)<-cnames
	
	}
	map_res<-unique(map_res)
	map_res<-as.data.frame(map_res)
	cnames<-c("mz","KEGGID", "Name", "Formula", "MonoisotopicMass", "Adduct","AdductMass")


	#colnames(map_res)<-cnames

	#rownames(map_res)<-NULL

	return(map_res)
	

}
