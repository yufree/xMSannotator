group_by_MDF <-
function(chemicalid,mchemicaldata,corthresh,global_cor,mzid,max_diff_rt=10,level_module_isop_annot,
adduct_table,adduct_weights)
{
	
	
	mchemicaldata$mz<-as.numeric(as.character(mchemicaldata$mz))
	
	mchemicaldata$time<-as.numeric(as.character(mchemicaldata$time))


	chemical_score<-(0)
	#print(chemicalid)
	#print(mchemicaldata)
	
	fname<-paste(chemicalid,"data.txt",sep="_")

	#write.table(mchemicaldata,file=fname,sep="\t",row.names=FALSE)
	
	#print("heress 1")
	data(adduct_table)
	
	if(length(mchemicaldata$mz)>0){
		
		if(length(which(duplicated(mchemicaldata$mz)==TRUE))>0){
			mchemicaldata<-mchemicaldata[-which(duplicated(mchemicaldata$mz)==TRUE),]
		}
		
	#mchemicaldata$Module_RTclust<-gsub(mchemicaldata$Module_RTclust,pattern="_[0-9]*",replacement="")

	if(is.na(adduct_weights)==TRUE){
	adduct_weights1<-matrix(nrow=2,ncol=2,0)
		adduct_weights1[1,]<-c("M+H",1)
		adduct_weights1[2,]<-c("M-H",1)
		adduct_weights<-as.data.frame(adduct_weights1)
		colnames(adduct_weights)<-c("Adduct","Weight")
	}
	mchemicaldata<-unique(mchemicaldata)

	bool_check<-0
	final_isp_annot_res<-mchemicaldata
	
	
	for(m in  1:length(mchemicaldata$mz))
	{
							
							isp_group<-as.character(mchemicaldata$ISgroup[m])
							
							module_rt_group<-as.character(mchemicaldata$Module_RTclust[m])
							
							query_md<-mchemicaldata$mz[m]-round(mchemicaldata$mz[m])
							
							#put_isp_masses_curmz_data<-level_module_isop_annot[which(abs((level_module_isop_annot$MD)-(query_md))<0.005),] #
							#put_isp_masses_curmz_data<-level_module_isop_annot[which(level_module_isop_annot$ISgroup==isp_group),]
							
							
							put_isp_masses_curmz_data<-level_module_isop_annot[which(abs((level_module_isop_annot$MD)-(query_md))<0.01 & level_module_isop_annot$Module_RTclust==module_rt_group),]
							
							put_isp_masses_curmz_data<-as.data.frame(put_isp_masses_curmz_data)
							
							put_isp_masses_curmz_data$mz<-as.numeric(as.character(put_isp_masses_curmz_data$mz))
							put_isp_masses_curmz_data$time<-as.numeric(as.character(put_isp_masses_curmz_data$time))
							mchemicaldata<-as.data.frame(mchemicaldata)
							
							#put_isp_masses_curmz_data$mz<-as.numeric(put_isp_masses_curmz_data$mz)
							
				#put_isp_masses_curmz_data<-put_isp_masses_curmz_data[which(put_isp_masses_curmz_data$mz>=mchemicaldata$mz[m] & put_isp_masses_curmz_data$mz<(mchemicaldata$mz[m]+3)),]
							
							#put_isp_masses_curmz_data<-put_isp_masses_curmz_data[,c(1:2,5)]
							put_isp_masses_curmz_data<-unique(put_isp_masses_curmz_data)
							
							put_isp_masses_curmz_data<-put_isp_masses_curmz_data[order(put_isp_masses_curmz_data$mz),]
							
								if(length(put_isp_masses_curmz_data)>0){
												if(isnum<0){
																			isp_sign<-"-"
																		}else{
																			isp_sign<-"+"
																			}
																		
																		#form_name<-as.character(paste(mchemicaldata[m,7],"_[+",(isnum),"]",sep=""))
																		
																		form_name<-as.character(paste(mchemicaldata[m,7],"_[",isp_sign,(abs(isnum)),"]",sep=""))
																		
																		
																		other_inf<-cbind(rep("-",7))
																		temp_var<-cbind(put_isp_masses_curmz_data[isp_v,c(1:2)],t(other_inf),put_isp_masses_curmz_data[isp_v,c(3:4)],put_isp_masses_curmz_data[isp_v,c(2,5)])
																		
																		temp_var<-as.data.frame(temp_var)
																		
																		colnames(temp_var)<-colnames(mchemicaldata)
																		
																		temp_var$Formula<-form_name
																		temp_var$Name<-as.character(mchemicaldata[m,6])
																		temp_var$chemical_ID<-as.character(mchemicaldata[m,5])
																		
																		#temp_var$Adduct<-paste(mchemicaldata[m,9],"_[+",isnum,"]",sep="")
																		
																		temp_var$Adduct<-paste(mchemicaldata[m,9],"_[",isp_sign,(abs(isnum)),"]",sep="")
																		
																		temp_var<-as.data.frame(temp_var)
																		
																		final_isp_annot_res<-as.data.frame(final_isp_annot_res)
																		
																		#write.table(temp_var,file="temp_var.txt",sep="\t")
																		#write.table(final_isp_annot_res,file="final_isp_annot_res.txt",sep="\t")
																		if(nrow(temp_var)>0){
																			
																				check_mz<-which(temp_var$mz%in%final_isp_annot_res)
																				
																				if(length(check_mz)>0){
																					temp_var<-temp_var[-c(check_mz),]
																				}
																				if(nrow(temp_var)>0){
																				final_isp_annot_res<-rbind(final_isp_annot_res,temp_var)
																				}
																			}
										
								}
	}
						
												
						
						#write.table(temp_var,file="finaltempvar.txt",sep="\t",row.names=FALSE)
						mchemicaldata<-final_isp_annot_res #[,-c(12)]
						
						#write.table(final_isp_annot_res,file="final_isp_annot_res3.txt",sep="\t")
																	
						
						dupmz<-which(duplicated(mchemicaldata$mz)==TRUE)
						
						if(length(dupmz)>0){
								mchemicaldata<-mchemicaldata[-c(dupmz),]
							
						}

	mchemicaldata_orig<-mchemicaldata
	
	table_mod<-table(mchemicaldata$Module_RTclust)
	
	table_iso<-table(mchemicaldata$ISgroup)
	
	table_mod<-table_mod[table_mod>0]
	
	table_mod<-table_mod[order(table_mod,decreasing=TRUE)]
	
	mchemicaldata_orig<-mchemicaldata
	
	#corthresh<-0.5
	top_mod<-names(table_mod)
	
	bool_check<-0

	#print(table_mod)
	
	topquant_cor<-0
	
	best_conf_level<-(-100)
	
	k_power<-1
	if(length(which(table_mod>1))>0)
	{
					best_chemical_score<-(-100)
			
					for(i in 1:length(which(table_mod>1)))
					{
					
						#print("here")
						
						mchemicaldata<-mchemicaldata_orig[which(mchemicaldata_orig$Module_RTclust==top_mod[i]),]
						
						mchemicaldata<-mchemicaldata[order(mchemicaldata$mz),]
						
						#mchemicaldata<-mchemicaldata[-which(is.na(mchemicaldata$mz)==TRUE),]
						
						#final_isp_annot_res<-mchemicaldata
			
						#bool_check<-0
											
						mzid_cur<-paste(mchemicaldata$mz,mchemicaldata$time,sep="_")
					
						cor_mz<-global_cor[which(mzid%in%mzid_cur),which(mzid%in%mzid_cur)]
						
						cor_mz<-abs(cor_mz)
						
						#print("length")
						#print(length(cor_mz))
						#print(dim(cor_mz))	
						if(length(cor_mz)>1){
									check_cor<-sapply(1:dim(cor_mz)[1],function(k){
										
										count_mz<-length(which(cor_mz[k,-k]>=corthresh))
										return(count_mz)
									})
						topquant_cor<-max(cor_mz[upper.tri(cor_mz)])
						
						#print(check_cor)
						
						check_cor2<-check_cor/length(check_cor)

						#print("check_cor2")
						#print(check_cor2)				
		
						#at least score of 2
								if(length(which(check_cor>0)==TRUE)>0)
								{
							
																hub_mz<-which(check_cor==max(check_cor)[1])
																
																hub_mz<-hub_mz[1]
														
																sub_cor_mz<-cor_mz[hub_mz,which(cor_mz[hub_mz,]>=corthresh)]
																
																#topquant_cor<-quantile(sub_cor_mz,0.9,na.rm=TRUE)[1]
																
																 if(is.na(topquant_cor)==TRUE){
                                                                                                                                	topquant_cor=0
                                                                                                                                }
                                                                                                                                
                                                               #  topquant_cor<-(topquant_cor)^6
                                                                 
																
																
														
																if(length(which(check_cor2>=0.25))>0){
																
																mchemicaldata<-mchemicaldata[which(cor_mz[hub_mz,]>=corthresh & check_cor2>=0.25),]
																#print(mchemicaldata)
																
																}else{
																
																mchemicaldata<-mchemicaldata[which(cor_mz[hub_mz,]>=corthresh),]
																}
																
																
																if(nrow(mchemicaldata)<2){
																next
																}
																
																mchemicaldata<-as.data.frame(mchemicaldata)
												
															
						
												
	#print("get confidence stage2")
		

		if(length(which(is.na(mchemicaldata$time))==TRUE)>0){
		mchemicaldata<-mchemicaldata[-which(is.na(mchemicaldata$time)==TRUE),]
		}	 
		
        #conf_level<-get_confidence_stage2(curdata=mchemicaldata,adduct_weights=adduct_weights)
	
    #print("conf_level")
    #	#print(conf_level)
		
        #	#print("here")
			
																							if(dim(mchemicaldata)[1]>1)
																							{
																									#chemical_score<-3				
				
																									k_power<-1					
																									
																								
#print("setting score 1")	
						chemical_score<-100*length(which(mchemicaldata$Adduct%in%adduct_weights[,1]))+1*(topquant_cor)*length(unique(mchemicaldata$Adduct))						
																								
																							}
																							names(chemical_score)<-chemicalid[1]
																							
								
									
									#return(list("chemical_score"=chemical_score,"filtdata"=mchemicaldata))
			

								}else{
								
							
									chemical_score<-0
												
									}
														
							}
							
							
							mchemicaldata<-mchemicaldata[order(mchemicaldata$Adduct),]
	
							

							
							#chemical_score<-chemical_score*(conf_level^conf_level)
							
							
							if(chemical_score>best_chemical_score)  #| conf_level>=best_conf_level)
							{
							
							
								
								best_chemical_score<-chemical_score
								
								best_conf_level<-conf_level
								
								
								best_mod_ind<-i
								best_data<-mchemicaldata
							
							
							}
							
							#if(FALSE)
							{
							#print("i is")
							#print(i)
							#print(mchemicaldata)
							#print(top_mod[i])
							#print("score is")
							#print(chemical_score)
							#print("topquant_cor")
							#print(topquant_cor)
							#print(diff_rt)
							#print(k_power)
							#print(conf_level)
							
							#print(length(mchemicaldata$Adduct))
							}
					}	
					
					#print("best mod")
				#	#print(best_mod_ind)
						if(best_chemical_score>0){
							chemical_score<-best_chemical_score
							best_mod<-i
							
							#print(head(best_data))
							#print(best_mod)
							
							mchemicaldata<-best_data
							names(chemical_score)<-chemicalid[1]

                                                        #return(list("chemical_score"=chemical_score,"filtdata"=mchemicaldata))
							}else{
							chemical_score<-0
							
							}
						
			}else{

						mchemicaldata<-mchemicaldata_orig
						
						mchemicaldata<-mchemicaldata[order(mchemicaldata$mz),]
						
						#mchemicaldata<-mchemicaldata[-which(is.na(mchemicaldata$mz)==TRUE),]
						
						#final_isp_annot_res<-mchemicaldata
			
						#bool_check<-0
											
						mzid_cur<-paste(mchemicaldata$mz,mchemicaldata$time,sep="_")
					
						cor_mz<-global_cor[which(mzid%in%mzid_cur),which(mzid%in%mzid_cur)]
						
						cor_mz<-abs(cor_mz)
						
						#print("length")
						#print(length(cor_mz))
						
						if(length(cor_mz)>1){
									check_cor<-sapply(1:dim(cor_mz)[1],function(k){
										
										count_mz<-length(which(cor_mz[k,-k]>=corthresh))
										return(count_mz)
									})
						topquant_cor<-max(cor_mz[upper.tri(cor_mz)])
						
						#print(check_cor)
						
						check_cor2<-check_cor/length(check_cor)
								#at least score of 2
								if(length(which(check_cor>0)==TRUE)>0)
								{
							
																hub_mz<-which(check_cor==max(check_cor)[1])
																
																hub_mz<-hub_mz[1]
														
																sub_cor_mz<-cor_mz[hub_mz,which(cor_mz[hub_mz,]>=corthresh)]
																
																#topquant_cor<-quantile(sub_cor_mz,0.9,na.rm=TRUE)[1]
																
																 if(is.na(topquant_cor)==TRUE){
                                                                                                                                	topquant_cor=0
                                                                                                                                }
                                                                                                                                
                                                               #  topquant_cor<-(topquant_cor)^6
                                                                 
																
																
														
																
																
																mchemicaldata<-mchemicaldata[which(cor_mz[hub_mz,]>=corthresh),]
																
																
																
																if(nrow(mchemicaldata)<2){
																next
																}
																
																mchemicaldata<-as.data.frame(mchemicaldata)
												
																
					
						
												
		
		if(length(which(is.na(mchemicaldata$time))==TRUE)>0){
		mchemicaldata<-mchemicaldata[-which(is.na(mchemicaldata$time)==TRUE),]
		}	 
		
		conf_level<-get_confidence_stage2(curdata=mchemicaldata,adduct_weights=adduct_weights)
		if(diff_rt<=max_diff_rt)
		{
			
			dup_mz_check<-which(duplicated(mchemicaldata$Adduct)==TRUE)
																							
																							if(length(dup_mz_check)>0){
																							mchemicaldata<-mchemicaldata[-c(dup_mz_check),]
																							}
																							if(dim(mchemicaldata)[1]>1){
																									#chemical_score<-3				
				
																									k_power<-1					
																									mchemicaldata$time<-as.numeric(as.character(mchemicaldata$time))					
																						#	#print("here")
																							
																						#	#print(topquant_cor)
																						#	#print(diff_rt)
																									
																			chemical_score<-100*length(which(mchemicaldata$Adduct%in%adduct_weights[,1]))+1*(topquant_cor)*(1/((diff_rt*0.1)+1)^k_power)*length(unique(mchemicaldata$Adduct))						
																								
																								}
																							names(chemical_score)<-chemicalid[1]
																							
									if(length(which(mchemicaldata$Adduct%in%adduct_weights[,1]))>0){
									
									fname<-paste(chemicalid,"score.txt",sep="_")

									#mchemicaldata<-cbind(chemical_score,mchemicaldata)
									
									#write.table(mchemicaldata,file=fname,sep="\t",row.names=FALSE)
	
									#     return(list("chemical_score"=chemical_score,"filtdata"=mchemicaldata))                       
                                                                        }
									
									#return(list("chemical_score"=chemical_score,"filtdata"=mchemicaldata))
			}

								}else{
									
									chemical_score<-0
												
									}
														
							}
							
							
							mchemicaldata<-mchemicaldata[order(mchemicaldata$Adduct),]
	
							
							
							#if(FALSE)
							{
							
							#print(mchemicaldata)
							
							#print("score is")
							#print(chemical_score)
							#print("topquant_cor")
							#print(topquant_cor)
							#print(diff_rt)
							#print(k_power)
							#print(conf_level)
							
							#print(length(mchemicaldata$Adduct))
							}
					}
					
			#######add code for only correlation criteria here
					

	
	return(list("chemical_score"=chemical_score,"filtdata"=mchemicaldata))
	}	
	
}
