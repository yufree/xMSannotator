multilevelannotationstep2 <-
function(outloc1,list_number){
    
    #,adduct_weights=NA,max.time.diff=NA,filter.by=c("M+H"),max_isp=100,numnodes=2,
    #MplusH.abundance.ratio.check=FALSE,mass_defect_window=0.01,mass_defect_mode="pos"){
	
	setwd(outloc1)
	
	load("step1_results.Rda")
	load("global_cor.Rda")
       unlink("allmatches_with_isotopes.txt")
	
	
	
	outloc<-outloc1
	if(is.na(max.rt.diff)==FALSE){
	
		max_diff_rt<-max.rt.diff
	}

	if(list_number>length(chemids_split)){
		list_number<-length(chemids_split)
		# stop(paste("Invalid list number. Must be less than or equal to ",length(chemids_split),sep=""))
		return(0)	
	}	
	#rm(dataA)
	#rm(levelA_res1)

	
	
		if(list_number>num_sets){
			list_number<-length(chemids_split)
			#stop(paste("Invalid list number. Must be less than or equal to ",num_sets,sep=""))
			return(0)		
        }


outloc1<-paste(outloc,"/stage2/",sep="")
suppressWarnings(dir.create(outloc1))
setwd(outloc1)			
if(is.na(adduct_weights)==TRUE){
 data(adduct_weights)

		adduct_weights1<-matrix(nrow=2,ncol=2,0)
		adduct_weights1[1,]<-c("M+H",1)
		adduct_weights1[2,]<-c("M-H",1)
		adduct_weights<-as.data.frame(adduct_weights1)
		colnames(adduct_weights)<-c("Adduct","Weight")
}
		
		#mchemdata<-m2
		cnames<-colnames(mchemdata)
		cnames[2]<-"time"
		colnames(mchemdata)<-as.character(cnames)
		#mchemdata<-unique(mchemdata)
		mchemdata$mz<-as.numeric(as.character(mchemdata$mz))
		mchemdata$time<-as.numeric(as.character(mchemdata$time))

	
		
			chem_score<-lapply(chemids_split[[list_number]],function(j){
			chemid<-chemids[j]
			chemscoremat<-{}
			 curmchemdata<-mchemdata[which(mchemdata$chemical_ID==chemid),]
			curmchemdata$mz<-as.numeric(as.character(curmchemdata$mz))
			curmchemdata$time<-as.numeric(as.character(curmchemdata$time))
			curmchemdata<-as.data.frame(curmchemdata)
			curmchemdata$Module_RTclust<-gsub(curmchemdata$Module_RTclust,pattern="_[0-9]*",replacement="")
			isop_res_md$Module_RTclust<-gsub(isop_res_md$Module_RTclust,pattern="_[0-9]*",replacement="")
			isp_masses_mz_data<-{}
                      isp_masses_mz_data<-lapply(1:length(curmchemdata$mz),function(m)
		       {
                                                        isp_group<-as.character(curmchemdata$ISgroup[m])
                                                        module_rt_group<-as.character(curmchemdata$Module_RTclust[m])
							module_rt_group<-gsub(module_rt_group,pattern="_[0-9]*",replacement="")
                                                        query_md<-curmchemdata$mz[m]-round(curmchemdata$mz[m])
                                                        put_isp_masses_curmz_data<-isop_res_md[which(abs((isop_res_md$MD)-(query_md))<mass_defect_window & isop_res_md$Module_RTclust==module_rt_group),]
                                                        put_isp_masses_curmz_data<-as.data.frame(put_isp_masses_curmz_data)
							return(put_isp_masses_curmz_data)
                        })
			isp_masses_mz_data<-ldply(isp_masses_mz_data,rbind)
			cnames<-colnames(isp_masses_mz_data)
 cnames[5]<-"AvgIntensity"
colnames(isp_masses_mz_data)<-cnames
			isp_masses_mz_data<-as.data.frame(isp_masses_mz_data)
	if(is.na(mass_defect_mode)==TRUE){
		mass_defect_mode="pos"
	}

				isp_masses_mz_data$mz<-as.numeric(as.character(isp_masses_mz_data$mz))
			isp_masses_mz_data$time<-as.numeric(as.character(isp_masses_mz_data$time))
			isp_masses_mz_data<-as.data.frame(isp_masses_mz_data)

	
        chem_score<-get_chemscorev1.6.71(chemicalid=chemid,mchemicaldata=curmchemdata,corthresh=corthresh,global_cor=global_cor,mzid,max_diff_rt=max_diff_rt,
        level_module_isop_annot=isp_masses_mz_data,adduct_table=adduct_table,adduct_weights=adduct_weights,filter.by,max_isp=max_isp,
	MplusH.abundance.ratio.check=MplusH.abundance.ratio.check,mass_defect_window=mass_defect_window,mass_defect_mode=mass_defect_mode,outlocorig=outloc)
	
	
	

	
        if(length(chem_score)>0){
                if(chem_score$chemical_score>=(-100)){

                        chem_score$filtdata<-chem_score$filtdata[order(chem_score$filtdata$mz),]

			if(nrow(chem_score$filtdata)>0){
                cur_chem_score<-chem_score$chemical_score
				chemscoremat<-cbind(cur_chem_score,chem_score$filtdata)
			}else{
				#print(chem_score)
			}
                        rm(chem_score)
			chemscoremat<-na.omit(chemscoremat)
                        chemscoremat<-as.data.frame(chemscoremat)
                }
        }else{
                rm(chem_score)
        }
       		rm("curmchemdata","isp_masses_mz_data","mzid_cur","chemid")
	
	suppressWarnings(rm(hmdbCompMZ))
    suppressWarnings(rm(hmdbAllinf))
	
	#chemscoremat<-chemscoremat[,c(1:7,9,11,14)]
	#chemscoremat$MatchCategory[which(chemscoremat$MatchCategory=="Multiple")]<-"M"
	#chemscoremat$MatchCategory[which(chemscoremat$MatchCategory=="Unique")]<-"U"
	#chemscoremat<-as.matrix(chemscoremat)

	return(chemscoremat)
			}) #,mc.cores=num_nodes,mc.preschedule=FALSE)
			
			#stopCluster(cl)
			cur_fname<-paste("chem_score",list_number,"_a.Rda",sep="")	
			#save(chem_score,file=cur_fname)
	
			#mc.cores=num_nodes,mc.preschedule=FALSE)
			#stopCluster(cl)
				
			chem_score2<-chem_score[which(chem_score!="NULL")]
			
			rm(chem_score)
			
			curchemscoremat <- ldply(chem_score2, rbind)
			rm(chem_score2)
			
			#save(chem_score2,file="chem_scoreA.Rda")
			cur_fname<-paste("chem_score",list_number,".Rda",sep="")
			save(curchemscoremat,file=cur_fname)
			#chemscoremat<-rbind(chemscoremat,curchemscoremat)
			
			Sys.sleep(1)

			rm("curchemscoremat","mchemdata","chemids","adduct_table","global_cor","mzid","max_diff_rt","isop_res_md","corthresh","level_module_isop_annot",
"chemids_split","corthresh","max.mz.diff","outloc","num_sets","db_name","num_nodes","num_sets","adduct_weights","filter.by")

			
	rm(list=ls())
	

	
}
