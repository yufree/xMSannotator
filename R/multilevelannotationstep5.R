multilevelannotationstep5 <-
function(outloc,max.mz.diff=5,max.rt.diff=30,adduct_weights=NA,filter.by=NA,min_ions_perchem=1,boostIDs=NA,max_isp=5,dbAllinf=NA,db_name="HMDB",chemscoremat=NA,num_nodes=2){
    
    setwd(outloc)
    max_diff_rt=max.rt.diff
   
if(is.na(chemscoremat)==TRUE){
    curated_res<-read.csv("Stage4.csv")
}else{

	curated_res<-chemscoremat
	rm(chemscoremat)
}
    curated_res<-as.data.frame(curated_res)
   curated_res$mz<-as.numeric(curated_res$mz)
    
    scorethresh<-0
    
  
   
    
    data(adduct_table)
    if(is.na(adduct_weights)==TRUE){
        data(adduct_weights)
        adduct_weights1<-matrix(nrow=2,ncol=2,0)
        adduct_weights1[1,]<-c("M+H",1)
        adduct_weights1[2,]<-c("M-H",1)
        adduct_weights<-as.data.frame(adduct_weights1)
        colnames(adduct_weights)<-c("Adduct","Weight")
    }
    adduct_table<-adduct_table[order(adduct_table$Adduct),]
   
	
	 
    t1<-table(curated_res$mz,curated_res$chemical_ID)
    t2<-apply(t1,1,sum)
    t2<-apply(t1,1,sum)
    
    multi_mz<-names(t2[which(t2>1)])
    
    curated_res$MatchCategory<-gsub(as.character(curated_res$MatchCategory),pattern="Multiple",replacement="Unique")
    curated_res$MatchCategory[which(curated_res$mz%in%multi_mz)]<-"Multiple"
    

    
    
    curated_res<-curated_res[order(curated_res$Confidence,curated_res$chemical_ID,curated_res$score,curated_res$Adduct,decreasing=TRUE),]
    
    

   

	#curated_res<-curated_res[order(curated_res$Confidence,decreasing=TRUE),]
    


    t1<-table(curated_res$mz,curated_res$chemical_ID)
    
    s1<-apply(t1,1,sum)
    
    rm(t1)
    
    mzunique<-colnames(which(s1==1))
    
    dunique<-curated_res[which(curated_res$mz%in%mzunique),]
    
    s2<-s1[which(s1>1)]
    mzdup<-names(s2)
    mzdup<-as.data.frame(mzdup)
    mzdup<-mzdup[,1]
    
    bad_ind<-{}
    #for(mind in 1:length(mzdup)){
    cl<-makeSOCKcluster(num_nodes)
	#print(format(Sys.time(), "%a %b %d %X %Y"))

   bad_ind<-foreach(mind=1:length(mzdup), .combine=rbind) %dopar%
  {
        mznum<-mzdup[mind]
        dmultsub<-curated_res[which(curated_res$mz%in%mznum),]
	dgood_add<-which(dmultsub$Adduct%in%adduct_weights[,1]) #=="M+H")
	if(length(dgood_add)>0){
	dmultsub$score[dgood_add]<-(dmultsub$score[dgood_add])*100
	}
        com_ind<-which(curated_res$mz%in%mznum)
        good_ind<-which(dmultsub$score==max(dmultsub$score,na.rm=TRUE))
        
        for(com_indval in 1:length(com_ind)){
            scoreval<-{}
            if(com_indval%in%good_ind==FALSE){
                #print(com_indval)
                dmat_com<-curated_res[which(curated_res$chemical_ID%in%dmultsub$chemical_ID[com_indval]),]
                scoreval<-((dim(dmat_com)[1])-1)*dmat_com$score[1]/(dim(dmat_com)[1])
                scorevec<-c(rep(scoreval,length(which(curated_res$chemical_ID%in%dmultsub$chemical_ID[com_indval]))))
                if(length(scorevec)<length(which(curated_res$chemical_ID%in%dmultsub$chemical_ID[com_indval]))){
                    break;
                }
                curated_res$score[which(curated_res$chemical_ID%in%dmultsub$chemical_ID[com_indval])]<-scorevec
            }

        }
        com_ind<-com_ind[-good_ind]
        
        #bad_ind<-c(bad_ind,com_ind)
        return(com_ind)
    }
   # print(format(Sys.time(), "%a %b %d %X %Y"))

      stopCluster(cl)
    
    if(length(bad_ind)>0){
        curated_res_unique<-curated_res[-c(bad_ind),]
    }else{
        curated_res_unique<-curated_res
    }
    
    good_ind<-which(curated_res_unique$score>=scorethresh)
    
    curated_res_unique_highconf<-{}
    if(length(good_ind)>0){
        curated_res_unique_highconf<-curated_res_unique[good_ind,]
    }
    
   curated_res<-curated_res_unique_highconf
 t2<-table(curated_res$mz)
    
    same1<-which(t2==1)
    
    uniquemz<-names(same1)
    
    curated_res$MatchCategory=rep("Multiple",dim(curated_res)[1])
    
    curated_res$MatchCategory[which(curated_res$mz%in%uniquemz)]<-"Unique"
    
    #write.table(curated_res,file="Stage5.txt",sep="\t",row.names=FALSE)

   write.csv(curated_res,file="Stage5.csv",row.names=FALSE)

   curated_res$mz<-sprintf("%.5f",as.numeric(curated_res$mz))
    
    
    
    
    
    chemIDs<-curated_res$chemical_ID
    htmllink<-curated_res$chemical_ID
    
    link_text<-chemIDs[1]
    t2<-gregexpr(pattern="HMDB",perl=FALSE,text=link_text)
    
    #if(length(t2)>1){
    if(db_name=="HMDB"){
            htmllink<-paste("<a href=http://www.hmdb.ca/metabolites/",chemIDs,">",chemIDs,"</a>",sep="")
    }else{
        
        if(db_name=="KEGG"){
            htmllink<-paste("<a href=http://www.genome.jp/dbget-bin/www_bget?",chemIDs,">",chemIDs,"</a>",sep="")
        }else{
            
            if(db_name=="LipidMaps"){
                htmllink<-paste("<a href=http://www.lipidmaps.org/data/LMSDRecord.php?LMID=",chemIDs,">",chemIDs,"</a>",sep="")
                
                            }else{
                
                if(db_name=="T3DB"){
                    
                    
                    
                    htmllink<-paste("<a href=http://www.t3db.ca/toxins/",chemIDs,">",chemIDs,"</a>",sep="")
                    
                }
                
            }
            
            
        }
        
        
    }
    
    fname=paste("Stage5_annotation_results",sep="")
    unlink(fname)
    
    curated_res$chemical_ID<-htmllink
    
#    HTMLInitFile(filename=fname,Title="Stage 5 annotation results", outdir=outloc)
#    fname=paste(outloc,"/Stage5_annotation_results.html",sep="")
 #   HTML(curated_res,file=fname,Border=1,innerBorder=1,useCSS=TRUE)
 #   HTMLEndFile(file=fname)
   
   #if(FALSE)
{ 
    outloc2<-paste(outloc,"/stage2/",sep="") 

    unlink(outloc2,force=TRUE,recursive=TRUE)
    
    outloc2<-paste(outloc,"/stage2",sep="")
    
    unlink(outloc2,force=TRUE,recursive=TRUE)
    
    file.remove(dir(outloc2, full.names=TRUE))

    outloc2<-paste(outloc,"/stage3/",sep="") 
   unlink(outloc2,force=TRUE,recursive=TRUE)
   file.remove(dir(outloc2, full.names=TRUE))
   outloc2<-paste(outloc,"/stage4/",sep="")
   file.remove(dir(outloc2, full.names=TRUE))
    unlink(outloc2,force=TRUE,recursive=TRUE)
    outloc2<-paste(outloc,"/stage5/",sep="")
   file.remove(dir(outloc2, full.names=TRUE))
   unlink(outloc2,force=TRUE,recursive=TRUE)

    suppressWarnings(unlink("*.Rda"))
    

    try(unlink("step1_results.Rda"),silent=TRUE)
    try(unlink("plot.pdf"),silent=TRUE)
    try(unlink("Rplots.pdf"),silent=TRUE)
    try(unlink("Rplots.pdf"),silent=TRUE)
}

 curated_res$chemical_ID<-chemIDs

curated_res<-as.data.frame(curated_res)

curated_res<-curated_res[order(curated_res$Confidence,decreasing=TRUE),]
   
print("Stage 5 confidence level distribution for unique chemical/metabolite IDs")
print(table(curated_res$Confidence[-which(duplicated(curated_res$chemical_ID)==TRUE)]))


print("Stage 5 confidence level distribution for unique chemical/metabolite formulas")
print(table(curated_res$Confidence[-which(duplicated(curated_res$Formula)==TRUE)]))


    return(curated_res)
    
}
