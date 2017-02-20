getscore <-
function(k,mchemdata,chemids,corthresh,global_cor,mzid,max_diff_rt,level_module_isop_annot,adduct_table,adduct_weights,mass_defect_window=0.01,mass_defect_mode="pos"){

        chemscoremat<-{}
        chemid<-chemids[k]
        curmchemdata<-mchemdata[which(mchemdata$chemical_ID==chemid),]
	rm(mchemdata)
	rm(chemids)
        curmchemdata<-as.data.frame(curmchemdata)
	
	#level_module_isop_annot<-level_module_isop_annot[which(level_module_isop_annot$mz

        chem_score<-get_chemscorev1.6.71(chemicalid=chemid,mchemicaldata=curmchemdata,corthresh=corthresh,global_cor,mzid,max_diff_rt=max_diff_rt,
        level_module_isop_annot=level_module_isop_annot,adduct_table=adduct_table,adduct_weights=adduct_weights,mass_defect_window=mass_defect_window,mass_defect_mode=mass_defect_mode)

        if(length(chem_score)>0){
                if(chem_score$chemical_score>=(-100)){

                     
                        
                        chem_score$filtdata<-chem_score$filtdata[order(chem_score$filtdata$mz),]

                        chemscoremat<-cbind(chem_score$chemical_score,chem_score$filtdata)
                        rm(chem_score)
                        
			chemscoremat<-na.omit(chemscoremat)
                        chemscoremat<-as.data.frame(chemscoremat)

                }
        }else{
                rm(chem_score)
        }
       
	return(chemscoremat)
}
