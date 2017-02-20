check_golden_rules <-
function(curformula,NOPS_check=FALSE){
	
	
	numnitrogens<-check_element(curformula,"N")
	numcarbons<-check_element(curformula,"C")
		
	numoxygens<-check_element(curformula,"O")
	
	numhydrogens<-check_element(curformula,"H")
	
	numphos<-check_element(curformula,"P")
	numsulphur<-check_element(curformula,"S")
			
		
	
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
