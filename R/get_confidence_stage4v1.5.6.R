get_confidence_stage4v1.5.6 <-
function(curdata,max_diff_rt,adduct_weights=NA,filter.by=NA){


		curdata<-curdata[order(curdata$Adduct),]
	
	cur_adducts_with_isotopes<-curdata$Adduct
	
	cur_adducts<-gsub(cur_adducts_with_isotopes,pattern="(_\\[\\+[0-9]\\])",replacement="")
	
if(is.na(adduct_weights)==TRUE){
data(adduct_table)
data(adduct_weights)
adduct_weights1<-matrix(nrow=2,ncol=2,0)
		adduct_weights1[1,]<-c("M+H",1)
		adduct_weights1[2,]<-c("M-H",1)
		adduct_weights<-as.data.frame(adduct_weights1)
		colnames(adduct_weights)<-c("Adduct","Weight")
	}
adduct_table<-adduct_table[order(adduct_table$Adduct),]
	
	curdata<-cbind(curdata,cur_adducts)
	
	curdata<-merge(curdata,adduct_table,by.x="cur_adducts",by.y="Adduct")
	
	curdata<-curdata[,-c(1)]
	
	formula_vec<-curdata$Formula

	temp1<-curdata
	temp1$Module_RTclust<-gsub(temp1$Module_RTclust,pattern="(_[0-9]*)",replacement="")
	
	table_modules<-table(temp1$Module_RTclust)

	rm(temp1)
	module_names<-names(table_modules[which(table_modules>0)])

	
	chemscoremat_conf_levels<-"High"
	
	if(length(which(adduct_table$Adduct%in%cur_adducts))<1)
    {
	
	chemscoremat_conf_levels<-"None"
	return(chemscoremat_conf_levels)	
	}

	if(length(unique(curdata$Adduct))<2 && curdata$score<10){

        chemscoremat_conf_levels<-"None"

	return(chemscoremat_conf_levels)
        }
    
    if(length(unique(curdata$Adduct))<2 && nrow(curdata)>1){
        
        chemscoremat_conf_levels<-"None"
        
        return(chemscoremat_conf_levels)
    }
	
	curdata$time<-as.numeric(as.character(curdata$time))
	delta_rt<-max(curdata$time)-min(curdata$time)
	
	delta_rt<-round(delta_rt)
	
	if(delta_rt>max_diff_rt)
	{
		#chemscoremat_conf_levels<-"Medium"
		
		#if(delta_rt>10){
		chemscoremat_conf_levels<-"None"
			
		#}
		
		#print("here")
		if(length(which(adduct_weights[,1]%in%cur_adducts))>0 && curdata$score>10){
			#curdata<-curdata[which(cur_adducts%in%adduct_weights[,1]),]
		
			#chemscoremat_conf_levels<-"Medium"		
	
			if(is.na(filter.by)==TRUE){
			mnum<-max(as.numeric(as.character(adduct_weights[which(adduct_weights[,1]%in%curdata$Adduct),2])))[1]
                                   curdata<-curdata[which(curdata$Adduct%in%adduct_weights[which(as.numeric(as.character(adduct_weights[,2]))>=mnum),1]),]
                        
				
			     	chemscoremat_conf_levels<-"Medium" 
			}else{

				if(length(which(cur_adducts%in%filter.by))>0){
					
					curdata<-curdata[which(curdata$Adduct%in%filter.by),]
					chemscoremat_conf_levels<-"Medium"
				}
			}
		}else{
			chemscoremat_conf_levels<-"None"

		}
	}
	else
	{
		#change in 1.5.3
		if(length(module_names)>1){
		chemscoremat_conf_levels<-"Medium"
		}
		else
		{
			if(nrow(curdata)<2){
				#chemscoremat_conf_levels<-"Low"
				if(length(which(cur_adducts%in%adduct_weights[,1]))<1){
                                        chemscoremat_conf_levels<-"Low"
                                }else{
							#matches an M+H and score is greater than 10
							chemscoremat_conf_levels<-"Medium"
				}
			}else{	
			min_nummol<-min(curdata$num_molecules)
			min_nummol_ind<-which(curdata$num_molecules==min_nummol)
			#num molecules check
			if(min_nummol>1){
					 chemscoremat_conf_levels<-"Low"
			}else{
		#Multiple molecules abundance check
		check1<-gregexpr(text=cur_adducts,pattern="([2-3]+M)")
		if(length(check1)>0)
		{
			#adduct_charges<-adduct_table[which(adduct_table$Adduct%in%cur_adducts),2]
			min_mol<-min(curdata$num_molecules)
			min_mol_ind<-which(curdata$num_molecules==min_mol)
	
			max_int_min_mol<-max(curdata$mean_int_vec[min_mol_ind])
	
			#max_int_min_charge<-max(curdata$mean_int_vec[min_mol_ind])

			bad_ind_status<-rep(0,length(check1))	
			
			min_mol<-min(adduct_table[which(adduct_table$Adduct%in%cur_adducts),2])
					
			#print("Max int")
			#print(min_mol)
			#print(max_int_min_mol)
			
			if(min_mol<2)
			{
			#check pattern of each adduct
			for(a1 in 1:length(check1))
			{
			
					strlength<-attr(check1[[a1]],"match.length")
					
					if(strlength[1]>(-1))
					{	
						
						#if(min_charge<2)
						{
							abund_ratio<-curdata$mean_int_vec[a1]/max_int_min_mol
			
							#print("abund ratio is")
							#print(abund_ratio)
							
							if(is.na(abund_ratio)==FALSE){
							if(abund_ratio>1){
	
									bad_ind_status[a1]<-1
									
									#chemscoremat_conf_levels<-"Low"
							}
							
							}else{
								#chemscoremat_conf_levels<-"Low"
								bad_ind_status[a1]<-1
							}
							
						}
					}
					
			}
			
			
			
				if(length(which(bad_ind_status==1))>0){
				
				#print(curdata)
				#print(bad_ind_status)
				bad_adducts<-cur_adducts[which(bad_ind_status==1)]
				
				#curdata<-curdata[-which(bad_ind_status==1),]
				#curdata<-curdata[-which(cur_adducts%in%bad_adducts),]
				
				#print(curdata)
				}
	
				if(length(nrow(curdata)>0)){	
				cur_adducts_with_isotopes<-curdata$Adduct

			        cur_adducts<-gsub(cur_adducts_with_isotopes,pattern="(_\\[\\+[0-9]\\])",replacement="")
				}else{
					
						chemscoremat_conf_levels<-"None"
					
				
				}
			}else{
			chemscoremat_conf_levels<-"Low"
			}

		
	

			#Multiply charged ions
			#if(FALSE)
			{
			
			adduct_charges<-curdata$charge # adduct_table[which(adduct_table$Adduct%in%cur_adducts),3]
			
			min_charge<-min(curdata$charge)
			
			min_charge_ind<-which(curdata$charge==min_charge)
	
			max_int_min_charge<-max(curdata$mean_int_vec) #[min_charge_ind])
	
			#max_int_min_charge<-max(curdata$mean_int_vec[min_charge_ind])

			bad_ind_status<-rep(1,length(check1))	
			
					
			
			if(min_charge>1)
			{
			
			
			chemscoremat_conf_levels<-"Low"
			}
			
			}
			
			#isotope based check for +1 and +2 to make sure that the abundant form is present
			check2<-gregexpr(text=cur_adducts_with_isotopes,pattern="(_\\[\\+[0-9]\\])")
				if(length(check2)>0){
																																			
										for(a1 in 1:length(check2)){
											strlength<-attr(check2[[a1]],"match.length")
											
											if(strlength[1]>(-1)){
											
											count_abundant_form<-length(which(cur_adducts%in%cur_adducts[a1]))
					
																																	
											if(count_abundant_form<2)
											{																	
														
														chemscoremat_conf_levels<-"None"
														
														curdata<-curdata[-c(a1),]
												
											}
										
											
											}
										}
																	
																	
						}
			
			
		

																	
		}
		
		
		
		
	}
		
		

	}
		
	
	}
	
}
	curformula<-as.character(formula_vec[1])
			curformula<-gsub(curformula,pattern="Ca",replacement="")
			curformula<-gsub(curformula,pattern="Cl",replacement="")
			curformula<-gsub(curformula,pattern="Cd",replacement="")
			curformula<-gsub(curformula,pattern="Cr",replacement="")
				
			numcarbons<-check_element(curformula,"C")
			if(numcarbons<1){
				
				chemscoremat_conf_levels<-"None"
			}

	if(FALSE){
	if(length(unique(curdata$Adduct))<2){

        chemscoremat_conf_levels<-"None"

        return(chemscoremat_conf_levels)
        }

	if(nrow(curdata)<2){
                                #chemscoremat_conf_levels<-"Low"
                                if(length(which(cur_adducts%in%adduct_weights[,1]))<1){
                                        chemscoremat_conf_levels<-"None"
                                }else{
				chemscoremat_conf_levels<-"Low"
				}
                        }
	}

	if(nrow(curdata)<1){	
	score_level<-0
	curdata<-c(-1)
	}else{
	
	#3: High
#2: medium
#1: Low
#0: None

score_level<-as.character(chemscoremat_conf_levels)

score_level<-gsub(score_level,pattern="High",replacement="3")
score_level<-gsub(score_level,pattern="Medium",replacement="2")
score_level<-gsub(score_level,pattern="Low",replacement="1")
score_level<-gsub(score_level,pattern="None",replacement="0")


	
	}
	
	final_res<-cbind(score_level,curdata)
	return(final_res)

}
