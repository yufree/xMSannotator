feat.batch.annotation.child <-
function(mz.val,max.mz.diff, adductname, syssleep,adduct_table)
{
	
	#adduct_table<-read.table("/Users/karanuppal/Documents/Emory/JonesLab/Projects/xMSannotator/CAMERA_Fiehnlab_adducts_calculated_masses_vjuly1014.txt",sep="\t",header=TRUE)
	#adduct_table<-adduct_table[c(which(adduct_table[,6]=="S"),which(adduct_table[,6]=="Acetonitrile")),]
	
	adduct_names<-as.character(adduct_table[,1])
	adductlist<-adduct_table[,4]
	mult_charge<-adduct_table[,3]
	num_mol<-adduct_table[,2]
	names(adductlist)<-as.character(adduct_names)
	names(mult_charge)<-as.character(adduct_names)
	names(num_mol)<-as.character(adduct_names)
	alladducts<-adduct_names
	
	
	#print(mz.val)
	#convert to neutral mass
	#mz=mz.val-adductmass

                                adductmass=adductlist[as.character(adductname)]
                                adductcharge=mult_charge[as.character(adductname)]
								adductnmol=num_mol[as.character(adductname)]
								#mz=(mz.val-adductmass)*adductcharge
                                #mz=(exact_mass/adductcharge)+adductmass	
                                
                                #mz=
                                
                                #mz=#((nmol*M)/charge+adductMass))
                                #M=((mz-adductMass)*charge)/nmol
                                
                                #mz=#((nmol*M)+adductMass))/charge
                                #M=((mz*charge)-adductMass)/nmol
                                				

                                #reverse
                                #mono_mass=((mz.val-adductmass)*adduct_charge)/(adductnmol)
                                mono_mass=((mz.val*adduct_charge)-(adductmass))/(adductnmol)
                                
                          mz=mono_mass
        delta_ppm=(max.mz.diff)*(mz/1000000)
        min_mz=round((mz-delta_ppm),5)
        max_mz=round((mz+delta_ppm),5)
	
	
	res={} #c("-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-")
	mzorig=round(mz.val,5)
	delta_ppm=round(delta_ppm,5)
	
	syssleep1<-(syssleep/5)
	Sys.sleep(syssleep1)
	
	write.table(mz.val,file="mzval.txt",sep="\t",row.names=FALSE)

	html_link<-"-"
	search_link=paste("http://rest.kegg.jp/find/compound/",min_mz,"-",max_mz,"/exact_mass",sep="")
	
	#print(search_link)
	d1<-try(readLines(search_link),silent=TRUE)

	if(is(d1,"try-error")){

	res<-c(mz.val,rep("NC",27))
	
	 write.table(mz.val,file="kegg_bad_mzs.txt",sep="\t",row.names=FALSE,append=TRUE)

	}else{
	
	
	
	#print(dim(d1))
	cnames<-c("ENTRY","NAME","FORMULA","EXACT_MASS","REACTION","PATHWAY","ENZYME","PubChem","ChEBI","PDB")
	
	#pattern_list<-c("C[0-9]{3,5}","[:blank:]{2,}[0-9|A-Z|:punct:|(|:print:][[:punct:]|[:alnum:]]*{3,}", "FORMULA","EXACT_MASS","CAS:","PubChem:","KNApSAcK:","PDB-CCD")
	
	#pattern_list<-c("C[0-9]{3,5}","NAME", "FORMULA","EXACT_MASS","ko[0-9]{5}","CAS:","ChEBI:","LIPIDMAPS:","PubChem:", "KNApSAcK:","PDB-CCD:")
	
	#pattern_list<-c("C[0-9]{3,5}","NAME", "FORMULA","EXACT_MASS","ko[0-9]{5}","CAS:","ChEBI:","LIPIDMAPS:","PubChem:", "KNApSAcK:","PDB-CCD:", "map")
	
	pattern_list<-c("EXACT_MASS","NAME", "FORMULA","CAS:","PubChem:","ChEBI:","LIPIDMAPS:", "BRITE", "map")
	
	pattern_keggid<-"C[0-9]{3,5}"
	
	#if(dim(d1)[1]>0)
	id_list<-"-"
		CName<-"-"
		mass<-"-"
		casID<-"-"
		keggID<-"-"
		kegglink<-"-"
		keggpathid<-"-"
		keggpathname<-"-"
		keggpathlink<-"-"
		hmdbID<-"-"
		hmdblink<-"-"
		pubchemsid<-"-"
		pubchemslink<-"-"
		pubchemcid<-"-"
		pubchemclink<-"-"
		chebiid<-"-"
		chebilink<-"-"
		lipidmapsid<-"-"
		lipidmapslink<-"-"
		chemformula<-"-"
		
	if(length(d1)>0){ 
	for(i in 1:length(d1))
	{
		if(i%%5>0){
		syssleep1<-(syssleep/5)
		Sys.sleep(syssleep1)
		}else{
		syssleep1<-(syssleep/3)
		Sys.sleep(syssleep1)
		}
		id_list<-"-"
		CName<-"-"
		mass<-"-"
		casID<-"-"
		keggID<-"-"
		kegglink<-"-"
		keggpathid<-"-"
		keggpathname<-"-"
		keggpathlink<-"-"
		hmdbID<-"-"
		hmdblink<-"-"
		pubchemsid<-"-"
		pubchemslink<-"-"
		pubchemcid<-"-"
		pubchemclink<-"-"
		chebiid<-"-"
		chebilink<-"-"
		lipidmapsid<-"-"
		lipidmapslink<-"-"
		chemformula<-"-"
		keggpathinf<-{}
		
		#l1<-grep(d1[i],pattern=pattern_list[5])
		str_text=d1[i]
		t2<-gregexpr(pattern=pattern_keggid,perl=FALSE,text=str_text)
			if(t2[[1]][1]>0)
			{
				t3=t2[[1]]
				strlength=attr(t3,"match.length")-1
				t4=strsplit(as.character(str_text),"")
				
				keggID<-t4[[1]][t3[1]:(t3[1]+strlength)]
			
				
				keggID<-paste(keggID,collapse="")
				kegglink<-paste("<a href=http://www.genome.jp/dbget-bin/www_bget?cpd:",keggID,">",keggID,"</a>",sep="")
				
				
				#html_res=readHTMLTable(kegglink)			
				search_link1=paste("http://rest.genome.jp/link/cpd:",keggID,"+-e",sep="")
				
				#dlink<-readLines(search_link1)
				dlink<-getURL(search_link1)
	
				if(dlink!=""){
					dlink<-read.delim(search_link1,header=FALSE)
					dlink2<-as.data.frame(dlink)
					if(dim(dlink2)[2]>0){
					for(l in 1:dim(dlink2)[1]){
						link_text=dlink2[l,2]
						t2<-gregexpr(pattern="HMDB[0-9]{2,}",perl=FALSE,text=link_text)
						t3=t2[[1]]
						strlength=attr(t3,"match.length")-1
						t4=strsplit(as.character(link_text),"")
						if(strlength>0)
						{
							hmdbID<-t4[[1]][t3[1]:(t3[1]+strlength)]
							hmdbID<-paste(hmdbID,collapse="")
							hmdblink<-paste("<a href=http://www.hmdb.ca/metabolites/",hmdbID,">",hmdbID,"</a>",sep="")
						}
					 }
					}
				}
				
				#keggID<-"C00392"
				
				#keggID<-"C00157"
		#keggID<-"C00082"
		search_link=paste("http://rest.kegg.jp/get/cpd:",as.character(keggID),sep="")
		d2<-read.delim(search_link,header=FALSE)
		d3<-as.data.frame(d2)
		
		
		
		
		if(length(d3)>0){
		pat.res<-{}
		url_vec<-{}
			url_strs<-c("-","-","-","-","<a href=http://pubchem.ncbi.nlm.nih.gov/summary/summary.cgi?sid=","<a href=http://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:",
			"<a href=http://www.lipidmaps.org/data/get_lm_lipids_dbgif.php?LM_ID=","-")

		url_vec<-{}
		
		if(length(d3)>0){
		for(j in 1:(length(pattern_list)-1))
			{
				l1<-grep(d3[,1],pattern=pattern_list[j])
				if(length(l1)>0){
				for(ind1 in 1:length(l1)){
				
				str1<-gsub(as.character(d3[l1[ind1],1]),pattern=" ",replacement="_")
				s1<-strsplit(str1," ")
				#print(s1)
				if(j==8){
				p1<-paste("(DBLINKS)|[_]{2,}|;*",pattern_list[j],sep="")
				s2<-gsub(as.character(s1[[1]]),pattern=p1,replacement="")
				s2<-gsub(s2,pattern="_",replacement=" ")
				}else{
				p1<-paste("(DBLINKS)|[_]*|:*|;*",pattern_list[j],sep="")
				s2<-gsub(as.character(s1[[1]]),pattern=p1,replacement="")
				
				}
				#p1<-paste("([DBLINKS])|[_]|:|;",pattern_list[j],sep="")
				
				
				
				
				#paste("<a href=http://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:",chebiid,">",chebiid,"</a>",sep="")
				
				url_str_cur<-paste(url_strs[j],s2,">",s2,"</a>",sep="")
				
				#if(ind1>1){s2<-paste(s2,";",sep="")}
				
				pat.res<-c(pat.res,s2)
				pat.res<-c(pat.res,url_str_cur)		
				}
				}else{
				pat.res<-c(pat.res,rep("-",2))
				}
			}
				l1<-grep(d3[,1],pattern=pattern_list[length(pattern_list)])
				if(length(l1)>0){
				
				keggpathid<-""
				keggpathname<-""
				keggpathlink<-""
				
				
				for(ind1 in 1:length(l1)){
				temp.pat.res<-{}
				str1<-gsub(as.character(d3[l1[ind1],1]),pattern=" ",replacement="_")
				s1<-strsplit(str1," ")
				
				
				p1<-paste("(DBLINKS)|[_]{3,}|:*|;*|PATHWAY",sep="")
				#p1<-paste("([DBLINKS])|[_]|:|;",pattern_list[j],sep="")
				s2<-gsub(as.character(s1[[1]]),pattern=p1,replacement="")
				
				s3<-strsplit(s2,"__")
				#s2<-gsub(s2,"__",";",sep="")
				
				#temp.pat.res<-c(temp.pat.res,s3[[1]][1])
				keggpathurl<-paste("<a href=http://www.genome.jp/kegg-bin/show_pathway?",s3[[1]][1],"+",keggID,">",s3[[1]][1],"</a>",sep="")
				#temp.pat.res<-c(temp.pat.res,keggpathlink)				
				s4<-gsub(as.character(s3[[1]][2]),pattern="_",replacement=" ")
				
				#temp.pat.res<-c(temp.pat.res,s4)
				keggpathid<-paste(keggpathid,paste(s3[[1]][1],";",sep=""),sep="")
				keggpathlink<-paste(keggpathlink,paste(keggpathurl,";",sep=""),sep="<br>")
				keggpathname<-paste(keggpathname,paste(s4,";",sep=""),sep="<br>")
				}
				
				pat.res<-c(pat.res,keggpathid,keggpathlink,keggpathname)
				
				}else{
				pat.res<-c(pat.res,"-","-","-")
				}
				
				search_link2=paste("http://rest.kegg.jp/link/disease/cpd:",as.character(keggID),sep="")
				d4<-try(read.delim(search_link2,header=FALSE),silent=TRUE)
				if(is(d4,"try-error")){
					diseaselink<-"-"
					diseaseitemlist<-"-"
				}else{
				d5<-as.data.frame(d4)
						diseaseitemlist<-d4[,2]
				maxnum<-10
				if(length(diseaseitemlist)>maxnum)
				{
						diseaseitemlist<-diseaseitemlist[1:maxnum]
					
				}
				if(length(d5)>0){
				diseaselink<-paste("<a href=http://www.genome.jp/dbget-bin/www_bget?",diseaseitemlist,">",diseaseitemlist,"</a>",sep="")
				diseaselink<-paste(as.vector(diseaselink),";<br>",collapse="")
				diseaseitemlist<-paste(diseaseitemlist,";",collapse="")
				}else{
				diseaselink<-"-"
				diseaseitemlist<-"-"
				}
				
				}
		
				
				search_link2=paste("http://rest.kegg.jp/link/environ/cpd:",as.character(keggID),sep="")
				d4<-try(read.delim(search_link2,header=FALSE),silent=TRUE)
				if(is(d4,"try-error")){
				environlink<-"-"
				environitemlist<-"-"
				}else{
				d5<-as.data.frame(d4)
					if(length(d5)>0){
								environitemlist<-d4[,2]
				maxnum<-10
				if(length(environitemlist)>maxnum)
				{
						environitemlist<-environitemlist[1:maxnum]
					
				}
				environlink<-paste("<a href=http://www.genome.jp/dbget-bin/www_bget?",environitemlist,">",environitemlist,"</a>",sep="")
				environlink<-paste(as.vector(environlink),";<br>",collapse="")
				environitemlist<-paste(environitemlist,";",collapse="")
				}else{
				environlink<-"-"
				environitemlist<-"-"
				}
				}
				
								
				search_link2=paste("http://rest.kegg.jp/link/drug/cpd:",as.character(keggID),sep="")
				d4<-try(read.delim(search_link2,header=FALSE),silent=TRUE)
				if(is(d4,"try-error")){
				druglink<-"-"
				drugitemlist<-"-"
				}else{
				d5<-as.data.frame(d4)
				if(length(d5)>0){
					drugitemlist<-d4[,2]
				maxnum<-10
				if(length(drugitemlist)>maxnum)
				{
						drugitemlist<-drugitemlist[1:maxnum]
					
				}
				druglink<-paste("<a href=http://www.genome.jp/dbget-bin/www_bget?",drugitemlist,">",drugitemlist,"</a>",sep="")
				druglink<-paste(as.vector(druglink),";<br>",collapse="")
				drugitemlist<-paste(drugitemlist,";",collapse="")
				}else{
				druglink<-"-"
				drugitemlist<-"-"
				}
				}
				
		}
		

	#pattern_list<-c("EXACT_MASS","NAME", "FORMULA","CAS:","C[0-9]{3,5}","PubChem:","ChEBI:","LIPIDMAPS:", "map")
	
		#res<-rbind(res,c(mzorig,delta_ppm,as.character(id_list), mass, html_link, CName,chemformula,casID,keggID,kegglink,keggpathid,keggpathname,keggpathlink,hmdbID,hmdblink,pubchemsid,pubchemslink, pubchemcid,pubchemclink,chebiid,chebilink, lipidmapsid, lipidmapslink))
				
	res<-rbind(res,c(mzorig,delta_ppm,as.character(id_list), pat.res[1], html_link, pat.res[3],pat.res[5],pat.res[7],keggID,kegglink,pat.res[17],pat.res[18],pat.res[19],hmdbID,hmdblink,pat.res[9],pat.res[10],pat.res[11],pat.res[12],pat.res[13],pat.res[14],pat.res[15],
	diseaselink,druglink,environlink,diseaseitemlist,drugitemlist,environitemlist))
			
	
	
     
	
			
		
		}
		}
		
		}
	}
	metres<-html_link
	#write.table(res,file="kegg_cur_res.txt",sep="\t",append=TRUE,row.names=FALSE)
	}
	syssleep1<-(syssleep/5)
	Sys.sleep(syssleep1)

	return(res)
}
