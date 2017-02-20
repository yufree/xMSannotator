check_element <-
function(curformula,elementname){

        curformula<-as.character(curformula)
        strsplit_var<-strsplit(curformula,split="")

        strsplit_elem<-strsplit(elementname,split="")
        elem_len<-length(strsplit_elem[[1]])

        g3<-gregexpr(curformula,pattern=paste(elementname,"[0-9]*",sep=""))

	numelement<-0
	
	if(length(g3)>0){
	if(elem_len>1){
		
		 regexp_len3<-attr(g3[[1]],"match.length")[1]
		if(regexp_len3>2){
		regexp_len3<-attr(g3[[1]],"match.length")[1]-elem_len+1
		}else{
		regexp_len3<-attr(g3[[1]],"match.length")[1]-elem_len	
		}
	}else{
		 regexp_len3<-attr(g3[[1]],"match.length")[1]-elem_len
	}
	
        if(regexp_len3>=0){

                if(regexp_len3==0){
                        numelement<-1

                }else{
                         numelement<-paste(strsplit_var[[1]][(g3[[1]][1]+elem_len):(g3[[1]][1]+regexp_len3)],collapse="")

                       numelement<-as.numeric(numelement)

                }
        }else{
                numelement<-0

        }

	}
        return(numelement)
}
