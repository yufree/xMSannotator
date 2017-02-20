ChemSpider.annotation <-
function(dataA,max.mz.diff=10, queryadductlist=c("M+H"), xMSannotator.outloc, numnodes=1,datasources=c("KEGG"),tokenstr=NA,maxhits=30,syssleep=5)
{
    
    
    data(adduct_table)
    adduct_table<-as.data.frame(adduct_table)
    
    data_a<-as.data.frame(dataA)
    
    #print("Using the 1st column as \"mz\" for annotation.")
    
    if(is.na(tokenstr)==TRUE){
        
        stop("Please specify a valid ChemSpider security token. \nCreate a ChemSpider account to obtain a token: \n www.chemspider.com/Register.aspx")
        
    }
    mzlist<-data_a[,1]
    
    #=c("PubChem","MassBank", "EPA DSSTox","EPA Toxcast","NIST Chemistry WebBook","KEGG",
    #"Human Metabolome Database", "ChEMBL", "ChEBI", "NIAID","Pesticide Common Names","SMPDB Small Molecule Pathway Database",
    #"MeSH","LipidMAPS","ChemBank","BioCyc")
    
    dir.create(xMSannotator.outloc,showWarnings=FALSE)
    setwd(xMSannotator.outloc)
    
    #adduct_table<-read.table("/Users/karanuppal/Documents/Emory/JonesLab/Projects/xMSannotator/adduct_table.txt",sep="\t",header=TRUE)
    #adduct_table<-adduct_table[c(which(adduct_table[,6]=="S"),which(adduct_table[,6]=="Acetonitrile")),]
    
    adduct_names<-as.character(adduct_table[,1])
    adductlist<-adduct_table[,4]
    mult_charge<-adduct_table[,3]
    num_mol<-adduct_table[,2]
    names(adductlist)<-as.character(adduct_names)
    names(mult_charge)<-as.character(adduct_names)
    names(num_mol)<-as.character(adduct_names)
    alladducts<-adduct_names
    
    if(queryadductlist[1]=="positive")
    {
        queryadductlist<-adduct_names[which(adduct_table[,5]=="positive")]
        
    }else{
        if(queryadductlist[1]=="negative")
        {
            
            queryadductlist<-adduct_names[which(adduct_table[,5]=="negative")]
            
        }else{
            if(queryadductlist[1]=="all"){
                
                
                queryadductlist<-alladducts
                
                
            }else{
                if(length(which(queryadductlist%in%alladducts==FALSE))>0){
                    
                    errormsg<-paste("Adduct should be one of:",sep="")
                    for(i in alladducts){errormsg<-paste(errormsg, i,sep=" ; ")}
                    stop(errormsg, "\n\nUsage: feat.batch.annotation.KEGG(dataA,max.mz.diff=10, queryadductlist=c(\"M+H\", \"M+Na\"), xMSannotator.outloc, numnodes=1)",
                    "\n\n OR  feat.batch.annotation.KEGG(dataA,max.mz.diff=10, queryadductlist=c(\"positive\"), xMSannotator.outloc, numnodes=1)",
                    "\n\n OR  feat.batch.annotation.KEGG(dataA,max.mz.diff=10, queryadductlist=c(\"negative\"), xMSannotator.outloc, numnodes=1)",
                    "\n\n OR  feat.batch.annotation.KEGG(dataA,max.mz.diff=10, queryadductlist=c(\"all\"), xMSannotator.outloc, numnodes=1)"
                    )
                }
                
            }
        }
    }
    
    parentres={}
    #cnames<-c("","Adduct","Query.m/z", "Search mass \n tolerance range (+/-)", "Metlin", "Compound.Name", "CASID", "KEGG.Compound.ID", "KEGG.Pathway.name", "KEGG.Pathway.ID", "HMDB.ID", "PubChem.Substance.ID", "PubChem.Compound.ID", "ChEBI.ID", "LIPID.MAPS.ID")
    #cnames<-c("","Adduct","Query.m/z", "Search mass \n tolerance range (+/-)", "ChemspiderID", "Link.to.ChemSpider","CommonName", "Molecular.Formula","SMILES","InChI", "InChIKey",
    #		"AverageMass","MolecularWeight","MonoisotopicMass", "NominalMass", "ALogP", "XLogP")
    
    cnames<-c("","Adduct","Query.m/z", "Search mass \n tolerance range (+/-)", "ChemspiderID", "CommonName", "Molecular.Formula","SMILES","InChI", "InChIKey",
    "AverageMass","MolecularWeight","MonoisotopicMass", "NominalMass", "ALogP", "XLogP", "Structure","KEGG.Compound.ID", "HMDB", "LipidMAPS", "PubChem",
    "MassBank", "BioCyc", "SMPDB","EPA.DSSTox","EPA.Toxcast", "Pesticide.Common.Names","ChEMBL","ChEBI","NIST.Chem.WebBook",
    "WikiPathways", "DrugBank", "Comparative Toxicogenomics Database",
    "ACToR: Aggregated Computational Toxicology Resource")
    for(adnum in 1:length(queryadductlist))
    {
        Sys.sleep(syssleep)
        
        adductname=queryadductlist[adnum]
        adductmass=adductlist[as.character(adductname)]
        
        cl<-makeSOCKcluster(numnodes)
        
        clusterEvalQ(cl, library(XML))
        clusterEvalQ(cl, library(R2HTML))
        clusterEvalQ(cl, library(RCurl))
        clusterEvalQ(cl, library(SSOAP))
        clusterEvalQ(cl, "processWSDL")
        clusterEvalQ(cl, library(png))
        clusterEvalQ(cl, "chspider.batch.annotation.child")
        
        #print(paste("Query adduct: ",adductname,sep=""))
        
        mz.annot.res<-new("list")
        min_mz<-min(mzlist)
        max_mz<-max(mzlist)
        
        #mz_group<-ceiling(max_mz/min_mz)
        mz_group<-10
        #mz_group<-round(max_mz/10)
        #length(mzlist)
        num_mz<-1
        for(mind in seq(1,length(mzlist),mz_group)){
            
            stopind<-mind+mz_group-1
            if(stopind>length(mzlist)){
                
                stopind<-length(mzlist)
            }
            s1<-mzlist[mind:stopind]
            s1<-unique(s1)
            num_mz<-num_mz+length(s1)
            if(num_mz%%50>0){
                Sys.sleep(syssleep)
            }else{
                Sys.sleep(syssleep)
            }
            
            if(length(s1)>2){
                
                repeat{
                    cur.annot.res<-parLapply(cl,s1,chspider.batch.annotation.child,max.mz.diff=max.mz.diff,adductname=adductname,datasources,tokenstr,maxhits,syssleep,adduct_table=adduct_table)
                    if(is(cur.annot.res,"try-error")){
                        Sys.sleep(10)
                        cur.annot.res<-parLapply(cl,s1,chspider.batch.annotation.child,max.mz.diff=max.mz.diff,adductname=adductname,datasources,tokenstr,maxhits,syssleep,adduct_table=adduct_table)
                        
                    }else{
                        break
                    }
                    
                }
                
                mz.annot.res<-c(mz.annot.res,cur.annot.res)
                
            }else{
                for(i in 1:length(s1)){
                    #print("mz is")
                    #print(s1[i])
                    rescur<-chspider.batch.annotation.child(mz.val=s1[i],max.mz.diff=max.mz.diff,adductname=adductname,datasources,tokenstr,maxhits,syssleep,adduct_table=adduct_table)
                    #print(length(rescur))
                    if(length(rescur)>0){
                        rescur<-as.matrix(rescur)
                        #print(dim(rescur))
                        if(dim(rescur)[2]==1){
                            rescur<-t(rescur)
                            rescur<-as.data.frame(rescur)
                        }
                        rescur<-as.data.frame(rescur)
                        #print(dim(rescur))
                        mz.annot.res[[length(mz.annot.res)+1L]]<-rescur
                    }
                }
            }
            if(mind%%10>0){
                Sys.sleep(syssleep/2)
            }else{
                Sys.sleep(1)
                stopCluster(cl)
                cl<-makeSOCKcluster(numnodes)
                
                
                clusterEvalQ(cl, library(XML))
                clusterEvalQ(cl, library(RCurl))
                clusterEvalQ(cl, library(SSOAP))
                clusterEvalQ(cl, library(png))
                clusterEvalQ(cl, "processWSDL")
                clusterEvalQ(cl, "feat.batch.annotation.child")
            }
            
            fname<-paste("cur_res.Rda",sep="")
            #save(mz.annot.res,file=fname)
        }
        
        stopCluster(cl)
        
        res={}
        
        #print(length(mz.annot.res))
        if(length(mz.annot.res)>0){
            #print(adductname)
            for(mzl in 1:length(mz.annot.res))
            {
                #print(dim(mz.annot.res[[mzl]]))
                res=rbind(res,mz.annot.res[[mzl]])
                
            }
        }
        #print(mz.annot.res)
        res<-unique(res)
        
        #print("length of res")
        #print(length(res))
        
        fname2<-paste(adductname,"res.Rda",sep="")
        #save(res,file=fname2)
        
        text_res<-{}
        
        if(length(res)>0){
            
            adductname1=c(rep(adductname,dim(res)[1]))
            
            temp_res<-cbind(adductname1,res)
            
            temp_res<-as.matrix(temp_res)
            
            
            
            if(dim(temp_res)[2]==1){
                
                temp_res<-t(temp_res)
                temp_res<-as.data.frame(temp_res)
                
            }
            
            bad_rows<-which(temp_res[,2]=="1")
            
            if(length(bad_rows)>0){
                temp_res<-temp_res[-bad_rows,]
                temp_res<-as.matrix(temp_res)
                
                #temp_res<-t(temp_res)
                if(dim(temp_res)[2]==1){
                    
                    temp_res<-t(temp_res)
                    
                    
                }
                
                
            }
            #temp_res<-as.data.frame(temp_res)
            colnames(temp_res)=NULL
            
            
            #text_resindex<-c(1,2,5,6,7,8,11,10,13,15,17,19,21)
            #text_resindex<-c(1,2,3,4,6:16)
            text_resindex<-c(1,2,3,6:9,10:17,seq(18,50,2))
            text_resindex<-text_resindex+1
            text_res<-temp_res[,c(1,text_resindex)]
            
            text_res<-as.matrix(text_res)
            
            
            
            if(dim(text_res)[2]==1){
                
                text_res<-t(text_res)
                
                
            }
            text_res<-as.data.frame(text_res)
            bad_rows<-which(text_res[,2]=="1")
            
            
            if(length(bad_rows)>0){
                text_res<-text_res[-bad_rows,]
                text_res<-as.matrix(text_res)
                text_res<-t(text_res)
            }
            text_res<-as.data.frame(text_res)
            
            sernum=seq(1,dim(text_res)[1])
            text_res<-cbind(sernum,text_res)
            colnames(text_res)=cnames
            
            
            parentres=rbind(parentres,temp_res)
            rm(temp_res)
            colnames(parentres)=NULL
            
        }
        
        fname=paste(xMSannotator.outloc,"/ChemSpider_annotation_results_",queryadductlist[adnum],".txt",sep="")
        write.table(text_res,file=fname,sep="\t",row.names=FALSE)
        
        
        
        Sys.sleep(2)
        
    }
    
    html_res<-{}
    
    
    
    if(length(parentres)>0){
        
        
        res<-parentres[order(parentres[,2]),]
        
        #html_resindex<-c(1,2,5,4,6:7,9,11:12,14,16,18,20,22)
        res<-as.matrix(res)
        
        if(dim(res)[2]==1){res<-t(res)}
        
        
        
        #html_resindex<-c(1,2,5,6:7,9,11:12,14,16,18,20,22)
        
        
        #text_resindex<-c(1,2,3,6:9,10:17,seq(18,50,2))
        html_resindex<-c(1,2,5:16,17,seq(19,51,2))
        html_resindex<-html_resindex+1
        html_res<-res[,c(1,html_resindex)]
        
        html_res<-as.matrix(html_res)
        if(dim(html_res)[2]==1){html_res<-t(html_res)}
        
        sernum=seq(1,dim(html_res)[1])
        
        #html_res<-cbind(sernum,html_res)
        html_res<-as.data.frame(html_res)
        #"ChemspiderID",
        #cnames<-c("","Adduct","Query.m/z", "Search mass \n tolerance range (+/-)", "Metlin.match.ID", "Metlin.match.mass", "Metlin.compound.name", "CASID", "KEGG.Compound.ID", "KEGG.Pathway.name", "KEGG.Pathway.ID", "HMDB.ID", "PubChem.Substance.ID", "PubChem.Compound.ID", "ChEBI.ID", "LIPID.MAPS.ID")
        
        cnames<-c("Adduct","Query.m/z", "Search mass \n tolerance range (+/-)", "Link.to.ChemSpider","CommonName", "Molecular.Formula","SMILES","InChI", "InChIKey",
        "AverageMass","MolecularWeight","MonoisotopicMass", "NominalMass", "ALogP", "XLogP","Structure","KEGG.Compound.ID", "HMDB", "LipidMAPS", "PubChem",
        "MassBank", "BioCyc", "SMPDB","EPA.DSSTox","EPA.Toxcast", "Pesticide.Common.Names","ChEMBL","ChEBI","NIST.Chem.WebBook",
        "WikiPathways", "DrugBank", "Comparative Toxicogenomics Database",
        "ACToR: Aggregated Computational Toxicology Resource")
        
        
        colnames(html_res)<-cnames
        rownames(html_res)<-sernum
        
        cnames2<-cnames[c(1:3,4,5:6,12,16,7:9,17:33)]
        html_res<-html_res[,c(1:3,4,5:6,12,16,7:9,17:33)]
        colnames(html_res)<-cnames2
        
        fname=paste("ChemSpider_annotation_results",sep="")
        unlink(fname)
        #fname=paste(xMSannotator.outloc,"/KEGG_annotation_results.html",sep="")
        HTMLInitFile(filename=fname,Title="ChemSpider annotation results", outdir=xMSannotator.outloc)
        fname=paste(xMSannotator.outloc,"/ChemSpider_annotation_results.html",sep="")
        HTML(html_res,file=fname,Border=1,innerBorder=1,useCSS=TRUE)
        HTMLEndFile(file=fname)
        
        if(length(queryadductlist)>1){
            
            text_res<-{}
            text_resindex<-c(1,2,3,6:16,17,seq(18,50,2))
            text_resindex<-text_resindex+1
            text_res<-res[,c(1,text_resindex)]
            
            text_res<-as.matrix(text_res)
            if(dim(text_res)[2]==1){text_res<-t(text_res)}
            
            if(length(text_res)>0){
                sernum=seq(1,dim(text_res)[1])
            }else{
                sernum={}
            }
            #text_res<-cbind(sernum,text_res)
            #"Link.to.ChemSpider",
            cnames<-c("Adduct","Query.m/z", "Search mass \n tolerance range (+/-)", "ChemspiderID", "CommonName", "Molecular.Formula","SMILES","InChI", "InChIKey",
            "AverageMass","MolecularWeight","MonoisotopicMass", "NominalMass", "ALogP", "XLogP", "Structure","KEGG.Compound.ID", "HMDB", "LipidMAPS", "PubChem",
            "MassBank", "BioCyc", "SMPDB","EPA.DSSTox","EPA.Toxcast", "Pesticide.Common.Names","ChEMBL","ChEBI","NIST.Chem.WebBook",
            "WikiPathways", "DrugBank", "Comparative Toxicogenomics Database",
            "ACToR: Aggregated Computational Toxicology Resource")
            text_res<-as.data.frame(text_res)
            colnames(text_res)=cnames
            #text_res<-text_res[,c(1:3,4,6,13,17,7:9,18:29)]
            #c("","Adduct","Query.m/z", "Search mass \n tolerance range (+/-)", "Metlin.match.ID", "Metlin.match.mass", "Metlin.compound.name", "CASID", "KEGG.Compound.ID", "KEGG.Pathway.ID", "KEGG.Pathway.name", "HMDB.ID", "PubChem.Substance.ID", "PubChem.Compound.ID", "ChEBI.ID", "LIPID.MAPS.ID")
            
            fname=paste(xMSannotator.outloc,"/ChemSpider_annotation_results_alladducts.txt",sep="")
            write.table(text_res,file=fname,sep="\t",row.names=FALSE)
        }
        
        
    }
    return(list("text.res"=text_res,"html.res"=html_res))
}
