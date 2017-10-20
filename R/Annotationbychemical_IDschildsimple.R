Annotationbychemical_IDschildsimple <-
function(adduct_index=NA,dataA,queryadductlist=c("M+H"),adduct_type=c("S","Acetonitrile"),adduct_table,max.mz.diff=10,outloc, otherdbs=FALSE,otherinfo=FALSE,keggCompMZ,chemIDs){

	#load("~/Documents/Emory/JonesLab/Projects/xMSannotator/keggCompMZ.Rda")
	dataA<-as.data.frame(dataA)
	adduct_names<-as.character(adduct_table[,1])
	adductlist<-adduct_table[,4]
	mult_charge<-adduct_table[,3]
	num_mol<-adduct_table[,2]
	names(adductlist)<-as.character(adduct_names)
	names(mult_charge)<-as.character(adduct_names)
	names(num_mol)<-as.character(adduct_names)
	alladducts<-adduct_names
	
    #print("here")
    queryadductlist<-as.character(queryadductlist)
        chemid_index<-adduct_index
	
	chemid_query<-chemIDs[chemid_index]
	
	chemid_index_in_compmz<-which(keggCompMZ$chemical_ID==chemid_query)
	keggCompMZ<-keggCompMZ[chemid_index,]
	
	adduct_index=NA
	if(is.na(adduct_index)==FALSE){
		
		queryadductlist=queryadductlist[adduct_index]
	}

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
	adduct_table<-as.data.frame(adduct_table)
    
  	
	adduct_table<-adduct_table[which(adduct_table$Type%in%adduct_type),] #=="S" | adduct_table$Type=="Acetonitrile",]
    
   
    adduct_table$Adduct<-as.character(adduct_table$Adduct)


suppressWarnings(dir.create(outloc))



setwd(outloc)

mz_search_list_1<-as.data.frame(keggCompMZ) #[which(keggCompMZ$Adduct%in%adduct_table$Adduct),c(1,7)])


#save(adduct_table,file="adduct_table.Rda")
#save(dataA,file="dataA.Rda")
#save(mz_search_list_1,file="mz_search_list_1.Rda")
#mz_search_list_1<-apply(mz_search_list_1,2,as.numeric)

#print(head(mz_search_list_1))

gcur<-getVenn(dataA=dataA,name_a="Experimental",name_b="DB",dataB=mz_search_list_1,mz.thresh=max.mz.diff,time.thresh=NA,
xMSanalyzer.outloc=outloc,alignment.tool=NA,plotvenn=FALSE)

#save(gcur,file="gcur.Rda")

#print("here 2")

if(length(gcur$common)>0){

mcur<-merge(keggCompMZ[which(keggCompMZ$Adduct%in%adduct_table$Adduct),],gcur$common,by.x="mz",by.y="mz.data.B")


mcur_2<-merge(mcur,dataA,by.x="mz.data.A",by.y="mz")

mcur_3<-mcur_2[order(mcur_2$Name),]

mcur_3<-mcur_3[,-c(9,10)]

cnames<-colnames(mcur_3)
cnames[1]<-"mz"
cnames[2]<-"theoretical.mz"
colnames(mcur_3)<-as.character(cnames)

mcur_4<-as.data.frame(mcur_3)

rm(keggCompMZ)
rm(dataA)
rm(mcur_2)
rm(mcur_3)
rm(mz_search_list_1)


if(dim(mcur)[1]>1){
	
	h1<-table(mcur_4$mz) #Adduct

if(length(h1)>0){
#u1<-c(u1,which(h1<=1))
u1<-which(h1<=1)
}
match_status<-rep("Multiple",dim(h1)[1])


uniq_kegg_matches<-names(u1)



match_status[u1]<-"Unique"

	
	
match_status_mat<-cbind(rownames(h1),match_status)



colnames(match_status_mat)<-c("mz","MatchCategory")
match_status_mat<-as.data.frame(match_status_mat)

mcur_5<-merge(match_status_mat,mcur_4,by="mz")

print(dim(mcur_5))

rm(mcur_4)

#h1<-table(mcur_5$chemical_ID,mcur_5$mz)



#s2<-apply(h1,2,sum)


mcur_5<-as.data.frame(mcur_5)

#mcur_5$mz<-as.numeric(mcur_5$mz)

mcur_5<-mcur_5[order(mcur_5$mz),]

}else{
	MatchCategory<-"Unique"
	cnames1<-colnames(mcur_4)
	cnames1<-c(cnames1[1],"MatchCategory",cnames[-c(1)])
	mcur_5<-cbind(mcur_4[1,1],MatchCategory,mcur_4[,-c(1)])
	mcur_5<-as.data.frame(mcur_5)
	colnames(mcur_5)<-as.character(cnames1)	


#mcur_5$mz<-as.numeric(mcur_5$mz)

mcur_5<-mcur_5[order(mcur_5$mz),]

	}
if(otherinfo==TRUE){
info_mat<-sapply(1:dim(mcur_5)[1],function(j){
	
	
	b1<-keggGet(paste("cpd:",mcur_5[j,1],sep=""))
	brite_inf<-paste(b1[[1]]$BRITE,collapse=";")
	path_inf<-paste(b1[[1]]$PATHWAYS,collapse=";")
	otherdb_inf<-paste(b1[[1]]$DBLINKS,collapse=";")
	r1<-c(as.character(mcur_5[j,1]),as.character(brite_inf),as.character(path_inf),as.character(otherdb_inf))

	return(r1)
})


info_mat_1<-as.data.frame(t(info_mat))
colnames(info_mat_1)<-c("chemical_ID","BriteCategory","Pathways","ExternalLinks")


mcur_6<-merge(info_mat_1,mcur_5,by="chemical_ID")

mcur_7<-unique(mcur_6)

rm(mcur_6)

if(otherdbs==TRUE){
info_mat_2<-sapply(1:dim(mcur_7)[1],function(j){
	
	b1<-keggLink(paste("cpd:",mcur_7[j,1],"+-e",sep=""))
	hmdbID<-"-"
	lipidmapsID<-"-"

		link_text<-b1[,2]
		
			t2<-gregexpr(pattern="hmdb:",perl=FALSE,text=link_text)
						
						if(length(t2)>1){
							g_ind<-which(t2==1)
							
							if(length(g_ind)>0){
								if(length(g_ind)>1){
								for(g in g_ind){
						t3=t2[[g]]
						
						hmdbID<-c(hmdbID,gsub(b1[g,2],pattern="hmdb:",replacement=""))
						}
						if(length(g_ind)>1){hmdbID<-paste(hmdbID,collapse=";")}
						}else{
							
							hmdbID<-gsub(b1[g_ind,2],pattern="hmdb:",replacement="")
						}
						}
						}
						
						
					t2<-gregexpr(pattern="lipidmaps:",perl=FALSE,text=link_text)
						
						if(length(t2)>1){
							g_ind<-which(t2==1)
							
							if(length(g_ind)>0){
								
								if(length(g_ind)>1){
								for(g in g_ind){
						t3=t2[[g]]
						
						lipidmapsID<-c(lipidmapsID,gsub(b1[g,2],pattern="lipidmaps:",replacement=""))
						
						
						}
						lipidmapsID<-paste(lipidmapsID,collapse=";")
						
						}else{lipidmapsID<-gsub(b1[g_ind,2],pattern="lipidmaps:",replacement="")}
						
						}
						
						}
							
					
					return(list(keggid=as.character(mcur_7[j,1]),hmdb=hmdbID,lipidmaps=lipidmapsID))
					c1<-c(as.character(mcur_7[j,1]),hmdbID,lipidmapsID)
					c1<-as.data.frame(c1)
					return(c1)		
	})
	
	info_mat_3<-{}
	#for(i in 1:dim(info_mat_2)[1]){
		
		cdata<-rbind(info_mat_2[1,],info_mat_2[2,])
		cdata<-rbind(cdata,info_mat_2[3,])
		cdata<-as.data.frame(cdata)
		info_mat_3<-rbind(info_mat_3,cdata)
		
		
	#}
	
#info_mat_3<-as.data.frame(t(info_mat_2))
info_mat_3<-t(info_mat_3)
colnames(info_mat_3)<-c("chemical_ID","HMDBID","LIPIDMAPS")

mcur_7<-as.data.frame(mcur_7)

mcur_8<-cbind(info_mat_3,mcur_7) #,by="chemical_ID")
mcur_8<-unique(mcur_8)
rownames(mcur_8)<-NULL
return(mcur_8)
}else{
	mcur_7<-as.data.frame(mcur_7)
	
	rownames(mcur_7)<-NULL
	return(mcur_7)
}
}else{
	mcur_5<-unique(mcur_5)
	return(mcur_5)
	
	}

	}else{return("no match found.")}
	#}else{return("no match found.")}
}
