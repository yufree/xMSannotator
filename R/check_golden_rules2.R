check_golden_rules2 <-
function(curformula,NOPS_check=FALSE){
	
	
	#curformula<-d2$Formula[1]
	curformula<-as.character(curformula)
	strsplit_var<-strsplit(curformula,split="")
	
	g1<-gregexpr(curformula,pattern="N[0-9]*")
	
	regexp_len<-attr(g1[[1]],"match.length")
	
	if(regexp_len>0){
		
	if(regexp_len==1){
		numnitrogens<-1
		
	}else{
	numnitrogens<-paste(strsplit_var[[1]][(g1[[1]][1]+1):(g1[[1]][1]+regexp_len-1)],collapse="")
	numnitrogens<-as.numeric(numnitrogens)
	}
	}else{
		numnitrogens<-0
		
	}
	
	
	g1<-gregexpr(curformula,pattern="C[0-9]*")
	
	regexp_len<-attr(g1[[1]],"match.length")
	numcarbons<-0
	if(regexp_len>0){
		if(regexp_len==1){
			numcarbons<-1
			
		}else{
			numcarbons<-paste(strsplit_var[[1]][(g1[[1]][1]+1):(g1[[1]][1]+regexp_len-1)],collapse="")
			numcarbons<-as.numeric(numcarbons)
		}
	}else{
		numcarbons<-0
		
	}
	
		
	
	g1<-gregexpr(curformula,pattern="O[0-9]*")
	regexp_len<-attr(g1[[1]],"match.length")
	
	if(regexp_len>0){
		if(regexp_len==1){
			numoxygens<-1
			
		}else{
			numoxygens<-paste(strsplit_var[[1]][(g1[[1]][1]+1):(g1[[1]][1]+regexp_len-1)],collapse="")
			numoxygens<-as.numeric(numoxygens)
		}
	}else{
		numoxygens<-0
		
	}

	g3<-gregexpr(curformula,pattern="H[0-9]*")
	regexp_len3<-attr(g3[[1]],"match.length")
numhydrogens<-0
	if(regexp_len3>0){
		if(regexp_len3==1){
			numhydrogens<-1
			
		}else{
			numhydrogens<-paste(strsplit_var[[1]][(g3[[1]][1]+1):(g3[[1]][1]+regexp_len3-1)],collapse="")
			numhydrogens<-as.numeric(numhydrogens)
		}
	}else{
		numhydrogens<-0
		
	}

	
	
	g3<-gregexpr(curformula,pattern="P[0-9]*")
	
	regexp_len3<-attr(g3[[1]],"match.length")
	
		if(regexp_len3>0){
			
		if(regexp_len3==1){
			numphos<-1
			
		}else{
			numphos<-paste(strsplit_var[[1]][(g3[[1]][1]+1):(g3[[1]][1]+regexp_len3-1)],collapse="")
	
			numphos<-as.numeric(numphos)

		}
	}else{
		numphos<-0
		
	}
	
	
	
	g3<-gregexpr(curformula,pattern="S[0-9]*")
	
	regexp_len3<-attr(g3[[1]],"match.length")
		
	if(regexp_len3>0){
			
		if(regexp_len3==1){
			numsulphur<-1
			
		}else{
			numsulphur<-paste(strsplit_var[[1]][(g3[[1]][1]+1):(g3[[1]][1]+regexp_len3-1)],collapse="")
	
			numsulphur<-as.numeric(numsulphur)

		}
	}else{
		numsulphur<-0
		
	}


	
	if(numcarbons<1){
		
		bool_check<-0
	}else{
	
		max_hydrogens<-(2*numcarbons)+numnitrogens+2
		
		nitrogen_to_carbon_ratio<-numnitrogens/numcarbons
		
		oxygen_to_carbon_ratio<-numoxygens/numcarbons
		
		phosphorus_to_carbon_ratio<-numphos/numcarbons
		
		sulphur_to_carbon_ratio<-numsulphur/numcarbons
	
		hydrogens_to_carbon_ratio<-numhydrogens/numcarbons
	
		bool_check<-1
		
		
		
		if(hydrogens_to_carbon_ratio<0.1 | hydrogens_to_carbon_ratio>6){
			
			bool_check=0
		}
		
		if(nitrogen_to_carbon_ratio>4){
			
			bool_check=0
		}
		
		if(oxygen_to_carbon_ratio>3){
			
			bool_check=0
		}
	
		if(phosphorus_to_carbon_ratio>2){
			
			bool_check=0
		}
		
		if(sulphur_to_carbon_ratio>3){
			
			bool_check=0
		}
		
		if(NOPS_check==TRUE){
		#NOPS>1
		if(numnitrogens>1 & numoxygens>1 & numphos>1 & numsulphur>1){
			if(numnitrogens>10 | numoxygens>20 | numphos>4 | numsulphur>3){
				
					bool_check<-0
				}
			
		}
		

		#NOP>3
		if(numnitrogens>3 & numoxygens>3 & numphos>3){
			if(numnitrogens>11 | numoxygens>22 | numphos>6){
				
					bool_check<-0
				}
			
		}
		
		#OPS>1
		if(numoxygens>1 & numphos>1 & numsulphur>1){
			if(numoxygens>14 | numphos>3 | numsulphur>3){
				
					bool_check<-0
				}
			
		}
		
		#PSN>1
		if(numnitrogens>1 & numphos>1 & numsulphur>1){
			if(numnitrogens>4 | numphos>3 | numsulphur>3){
				
					bool_check<-0
				}
			
		}

		#NOS>6
		if(numnitrogens>6 & numoxygens>6 & numsulphur>6){
			if(numnitrogens>19 | numoxygens>14 | numsulphur>8){
				
					bool_check<-0
				}
			
		}
		
		}
		
	
	}
	res<-cbind(curformula,bool_check)
	res<-as.data.frame(res)
	return(res)
	
}
