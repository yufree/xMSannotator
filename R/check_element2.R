check_element2 <-
function(curformula,elementname){
	
	curformula<-as.character(curformula)
	strsplit_var<-strsplit(curformula,split="")
	

	g3<-gregexpr(curformula,pattern=paste(elementname,"[0-9]*",sep=""))
	
	regexp_len3<-attr(g3[[1]],"match.length")
		
	if(regexp_len3>0){
			
		if(regexp_len3==1){
			numelement<-1
			
		}else{
			numelement<-paste(strsplit_var[[1]][(g3[[1]][1]+1):(g3[[1]][1]+regexp_len3-1)],collapse="")
	
			numelement<-as.numeric(numelement)

		}
	}else{
		numelement<-0
		
	}

	return(numelement)
}
