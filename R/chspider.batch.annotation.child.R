chspider.batch.annotation.child <-
function(mz.val,max.mz.diff, adductname,datasources=c("KEGG"),tokenstr,maxhits=30,syssleep=10,adduct_table){
    
    #data(adduct_table)
    annot_res={}
    
    Sys.sleep(syssleep);
    
    cs<-processWSDL("http://www.chemspider.com/MassSpecAPI.asmx?WSDL")
    o<-genSOAPClientInterface(,cs)
    
    cs2<-processWSDL("http://www.chemspider.com/Search.asmx?WSDL")
    o2<-genSOAPClientInterface(,cs2)
    
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
    
    
    #print(mz.val)
    #convert to neutral mass
    #mz=mz.val-adductmass
    
    adductmass=adductlist[as.character(adductname)]
    adduct_charge=mult_charge[as.character(adductname)]
    adductnmol=num_mol[as.character(adductname)]
    #mz=(mz.val-adductmass)*adductcharge
    #mz=(exact_mass/adductcharge)+adductmass
    
    #mz=(mz.val*adductcharge)-(adductmass)
    #mz=#((nmol*M)/charge)+adductMass
    
    # mz=((nmol*M)+adductMass)/charge
    #M=((mz*charge)-adductMass)/nmol
    
    
    #M=((mz.val-adductmass)*adductcharge)/(adductnmol)
    
    #reverse
    # mono_mass=((mz.val-adductmass)*adduct_charge)/(adductnmol)
    
    mono_mass=((mz.val*adduct_charge)-(adductmass))/(adductnmol)
    
    mz=mono_mass
    
    #mz=(exact_mass/adductcharge)+adductmass
    delta_ppm=(max.mz.diff)*(mz/1000000)
    min_mz=round((mz-delta_ppm),5)
    max_mz=round((mz+delta_ppm),5)
    
    
    
    #csids<-o@functions$SearchByMass2(mass=mz,range=delta_ppm)
    
    #csids<-o@functions$SearchByMassAsync(mass=mz,range=delta_ppm,dbs="KEGG",token=tokenstr)
    # csids<-o@functions$SearchByMassAsync(mass=mz,range=delta_ppm,dbs=c("PubChem","MassBank", "EPA DSSTox","EPA Toxcast","NIST Chemistry WebBook","KEGG",
    #"Human Metabolome Database", "ChEMBL", "ChEBI", "NIAID","Pesticide Common Names","SMPDB Small Molecule Pathway Database",
    #"MeSH","LipidMAPS","ChemBank","BioCyc"),token=tokenstr)
    #ttest_res=try(t.test(x,y,paired=TRUE),silent=TRUE)
    csids<-try(o@functions$SearchByMassAsync(mass=mz,range=delta_ppm,dbs=datasources,token=tokenstr),silent=TRUE)
    
    if(is(csids,"try-error")){
        
        annotres<-c(mz.val,rep("NC",50))
        write.table(mz.val,file="chemspider_bad_mzs.txt",sep="\t",append=TRUE,row.names=FALSE)
        
    }else{
        csids<-try(o2@functions$GetAsyncSearchResult(rid=csids,token=tokenstr),silent=TRUE)
        
        
        if(is(csids,"try-error")){
            
            annotres<-c(mz.val,rep("NC",50))
            
            write.table(mz.val,file="chemspider_bad_mzs.txt",sep="\t",append=TRUE,row.names=FALSE)
            
        }else{
            
            #WikiPathways, DrugBank, Comparative Toxicogenomics Database, ACToR: Aggregated Computational Toxicology Resource
            
            external_sources<-c("PubChem","MassBank", "EPA DSSTox","EPA Toxcast","NIST Chemistry WebBook","KEGG",
            "Human Metabolome Database", "ChEMBL", "ChEBI", "NIAID","Pesticide Common Names","SMPDB Small Molecule Pathway Database",
            "MeSH","LipidMAPS","ChemBank","BioCyc","WikiPathways", "DrugBank", "Comparative Toxicogenomics Database",
            "ACToR: Aggregated Computational Toxicology Resource")
            
            write.table(mz.val,file="mzval.txt",sep="\t",row.names=FALSE)
            
            if(length(csids)>0){
                
                
                countid<-0
                
                if(length(csids)<maxhits){
                    maxhits=length(csids)
                }
                for(id_ind in 1:maxhits)
                {
                    Sys.sleep(syssleep)
                    csid<-csids[id_ind]
                    csid1<-o2@functions$CSID2ExtRefs(CSID=as.integer(csid),datasources=external_sources,token=tokenstr)
                    info<-o@functions$GetExtendedCompoundInfo(CSID=csid,token=tokenstr)
                    urllink1<-paste("www.chemspider.com/Chemical-Structure.",csid,".html",sep="")
                    urllink2<-paste("<a href=http://www.chemspider.com/Chemical-Structure.",csid,".html>",csid,"</a>",sep="")
                    CommonName="-"
                    MF="-"
                    SMILES="-"
                    InChI="-"
                    InChIKey="-"
                    AverageMass="-"
                    MolecularWeight="-"
                    MonoisotopicMass="-"
                    NominalMass="-"
                    ALogP="-"
                    XLogP="-"
                    
                    slotnames<-slotNames(class(info))
                    
                    for(i in slotnames){
                        if(i=="CommonName"){
                            if(length(info@CommonName)>0){
                                CommonName=info@CommonName}
                        }else{
                            if(i=="MF"){
                                if(length(info@MF)>0){
                                    MF=info@MF
                                }
                            }else{
                                if(i=="SMILES"){
                                    if(length(info@SMILES)>0){
                                        SMILES=info@SMILES
                                    }
                                }else{
                                    if(i=="InChI"){
                                        if(length(info@InChI)>0){
                                            InChI=info@InChI
                                        }
                                    }else{
                                        if(i=="InChIKey"){
                                            if(length(info@InChIKey)>0){
                                                InChIKey=info@InChIKey
                                            }
                                        }else{
                                            if(i=="AverageMass"){
                                                if(length(info@AverageMass)>0){
                                                    AverageMass=info@AverageMass
                                                }
                                            }else{
                                                if(i=="MolecularWeight"){
                                                    if(length(info@MolecularWeight)>0){
                                                        MolecularWeight=info@MolecularWeight
                                                    }
                                                }else{
                                                    if(i=="MonoisotopicMass"){
                                                        if(length(info@MonoisotopicMass)>0){
                                                            MonoisotopicMass=info@MonoisotopicMass
                                                        }
                                                    }else{
                                                        if(i=="NominalMass"){
                                                            if(length(info@NominalMass)>0){
                                                                NominalMass=info@NominalMass
                                                            }
                                                        }else{
                                                            if(i=="ALogP"){
                                                                if(length(info@ALogP)>0){
                                                                    ALogP=info@ALogP
                                                                }
                                                            }else{
                                                                if(i=="XLogP"){
                                                                    if(length(info@XLogP)>0){
                                                                        XLogP=info@XLogP
                                                                    }
                                                                }
                                                                
                                                            }
                                                            
                                                        }
                                                        
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                    
                    
                    csid1NAME<-o2@functions$GetCompoundThumbnail(id=csid,token=tokenstr)
                    
                    dec<-base64Decode(csid1NAME,"raw")
                    pngoutloc<-"pngfiles"
                    dir.create(pngoutloc)
                    png_image<-readPNG(dec)
                    fname<-paste(pngoutloc,"/",csid,".png",sep="")
                    
                    
                    writePNG(png_image,target=fname)
                    
                    #img_url<-paste("<img src=",fname," alt=", csid," width=100 height=100>",sep="")
                    img_url<-paste("<img src=\"",fname,"\" alt=", csid," width=100 height=100>",sep="")
                    
                    #if(FALSE)
                    {
                        
                        
                        ID.PubChem="-"
                        ID.MassBank="-"
                        ID.EPA.DSSTox="-"
                        ID.EPA.Toxcast="-"
                        ID.NIST.Chem.WebBook="-"
                        ID.KEGG="-"
                        ID.HMDB="-"
                        ID.ChEMBL="-"
                        ID.ChEBI="-"
                        ID.Pesticide.Common.Names="-"
                        ID.SMPDB="-"
                        ID.MeSH="-"
                        ID.LipidMAPS="-"
                        ID.ChemBank="-"
                        ID.BioCyc="-"
                        ID.WikiPathways="-"
                        ID.DrugBank="-"
                        ID.CTD="-"
                        ID.ACToR="-"
                        
                        
                        URL.PubChem="-"
                        URL.MassBank="-"
                        URL.EPA.DSSTox="-"
                        URL.EPA.Toxcast="-"
                        URL.NIST.Chem.WebBook="-"
                        URL.KEGG="-"
                        URL.HMDB="-"
                        URL.ChEMBL="-"
                        URL.ChEBI="-"
                        URL.Pesticide.Common.Names="-"
                        URL.SMPDB="-"
                        URL.MeSH="-"
                        URL.LipidMAPS="-"
                        URL.ChemBank="-"
                        URL.BioCyc="-"
                        URL.WikiPathways="-"
                        URL.DrugBank="-"
                        URL.CTD="-"
                        URL.ACToR="-"
                        
                        if(length(csid1)>0){
                            for(dbindex in 1:length(csid1)){
                                if(length(csid1[dbindex]$ExtRef@ext_id)>0 && length(csid1[dbindex]$ExtRef@ds_name)>0){
                                    if(csid1[dbindex]$ExtRef@ds_name=="PubChem"){
                                        ID.PubChem=csid1[dbindex]$ExtRef@ext_id
                                        URL.PubChem<-paste("<a href=",csid1[dbindex]$ExtRef@ext_url,">",ID.PubChem,"</a>",sep="")
                                    }else{
                                        
                                        if(csid1[dbindex]$ExtRef@ds_name=="MassBank"){
                                            
                                            ID.MassBank=csid1[dbindex]$ExtRef@ext_id
                                            #URL.MassBank<-paste("<a href=",csid1[dbindex]$ExtRef@ext_url,">",ID.MassBank,"</a>",sep="")
                                            #http://massbank.normandata.eu/MassBank/jsp/Dispatcher.jsp?type=disp&id=PB000410&site=25
                                            URL.MassBank<-csid1[dbindex]$ExtRef@ext_id
                                            #paste("<a href=http://massbank.normandata.eu/MassBank/jsp/Dispatcher.jsp?type=disp&id=",ID.MassBank,"&site=25>",ID.MassBank,"</a>",sep="")
                                            
                                        }else{
                                            
                                            if(csid1[dbindex]$ExtRef@ds_name=="EPA DSSTox"){
                                                ID.EPA.DSSTox=csid1[dbindex]$ExtRef@ext_id
                                                URL.EPA.DSSTox<-paste("<a href=",csid1[dbindex]$ExtRef@ext_url,">",ID.EPA.DSSTox,"</a>",sep="")
                                            }else{
                                                
                                                if(csid1[dbindex]$ExtRef@ds_name=="EPA Toxcast"){
                                                    ID.EPA.Toxcast=csid1[dbindex]$ExtRef@ext_id
                                                    URL.EPA.Toxcast<-paste("<a href=",csid1[dbindex]$ExtRef@ext_url,">",ID.EPA.Toxcast,"</a>",sep="")
                                                }else{
                                                    
                                                    if(csid1[dbindex]$ExtRef@ds_name=="NIST Chemistry WebBook"){
                                                        ID.NIST.Chem.WebBook=csid1[dbindex]$ExtRef@ext_id
                                                        URL.NIST.Chem.WebBook<-paste("<a href=",csid1[dbindex]$ExtRef@ext_url,">",ID.NIST.Chem.WebBook,"</a>",sep="")
                                                    }else{
                                                        
                                                        if(csid1[dbindex]$ExtRef@ds_name=="KEGG"){
                                                            ID.KEGG=csid1[dbindex]$ExtRef@ext_id
                                                            #URL.KEGG=csid1[dbindex]$ExtRef@ext_url
                                                            URL.KEGG<-paste("<a href=",csid1[dbindex]$ExtRef@ext_url,">",ID.KEGG,"</a>",sep="")
                                                            
                                                        }else{
                                                            
                                                            if(csid1[dbindex]$ExtRef@ds_name=="Human Metabolome Database"){
                                                                ID.HMDB=csid1[dbindex]$ExtRef@ext_id
                                                                URL.HMDB<-paste("<a href=",csid1[dbindex]$ExtRef@ext_url,">",ID.HMDB,"</a>",sep="")
                                                            }else{
                                                                
                                                                if(csid1[dbindex]$ExtRef@ds_name=="ChEMBL"){
                                                                    ID.ChEMBL=csid1[dbindex]$ExtRef@ext_id
                                                                    URL.ChEMBL<-paste("<a href=",csid1[dbindex]$ExtRef@ext_url,">",ID.ChEMBL,"</a>",sep="")
                                                                }else{
                                                                    
                                                                    if(csid1[dbindex]$ExtRef@ds_name=="ChEBI"){
                                                                        ID.ChEBI=csid1[dbindex]$ExtRef@ext_id
                                                                        URL.ChEBI<-paste("<a href=",csid1[dbindex]$ExtRef@ext_url,">",ID.ChEBI,"</a>",sep="")
                                                                    }else{
                                                                        
                                                                        if(csid1[dbindex]$ExtRef@ds_name=="Pesticide Common Names"){
                                                                            ID.Pesticide.Common.Names=csid1[dbindex]$ExtRef@ext_id
                                                                            URL.Pesticide.Common.Names<-paste("<a href=",csid1[dbindex]$ExtRef@ext_url,">",ID.Pesticide.Common.Names,"</a>",sep="")
                                                                        }else{
                                                                            
                                                                            if(csid1[dbindex]$ExtRef@ds_name=="SMPDB Small Molecule Pathway Database"){
                                                                                ID.SMPDB=csid1[dbindex]$ExtRef@ext_id
                                                                                URL.SMPDB<-paste("<a href=",csid1[dbindex]$ExtRef@ext_url,">",ID.SMPDB,"</a>",sep="")
                                                                            }else{
                                                                                
                                                                                if(csid1[dbindex]$ExtRef@ds_name=="MeSH"){
                                                                                    ID.MeSH=csid1[dbindex]$ExtRef@ext_id
                                                                                    URL.MeSH<-paste("<a href=",csid1[dbindex]$ExtRef@ext_url,">",ID.MeSH,"</a>",sep="")
                                                                                }else{
                                                                                    
                                                                                    if(csid1[dbindex]$ExtRef@ds_name=="LipidMAPS"){
                                                                                        ID.LipidMAPS=csid1[dbindex]$ExtRef@ext_id
                                                                                        URL.LipidMAPS<-paste("<a href=",csid1[dbindex]$ExtRef@ext_url,">",ID.LipidMAPS,"</a>",sep="")
                                                                                    }else{
                                                                                        
                                                                                        if(csid1[dbindex]$ExtRef@ds_name=="ChemBank"){
                                                                                            ID.ChemBank=csid1[dbindex]$ExtRef@ext_id
                                                                                            URL.ChemBank<-paste("<a href=",csid1[dbindex]$ExtRef@ext_url,">",ID.ChemBank,"</a>",sep="")
                                                                                        }else{
                                                                                            
                                                                                            if(csid1[dbindex]$ExtRef@ds_name=="BioCyc"){
                                                                                                ID.BioCyc=csid1[dbindex]$ExtRef@ext_id
                                                                                                URL.BioCyc<-paste("<a href=",csid1[dbindex]$ExtRef@ext_url,">",ID.BioCyc,"</a>",sep="")
                                                                                            }else{
                                                                                                
                                                                                                if(csid1[dbindex]$ExtRef@ds_name=="WikiPathways"){
                                                                                                    ID.WikiPathways=csid1[dbindex]$ExtRef@ext_id
                                                                                                    URL.WikiPathways<-paste("<a href=",csid1[dbindex]$ExtRef@ext_url,">",ID.WikiPathways,"</a>",sep="")
                                                                                                }else{
                                                                                                    
                                                                                                    if(csid1[dbindex]$ExtRef@ds_name=="DrugBank"){
                                                                                                        ID.DrugBank=csid1[dbindex]$ExtRef@ext_id
                                                                                                        #URL.DrugBank<-paste("<a href=",csid1[dbindex]$ExtRef@ext_url,">",ID.DrugBank,"</a>",sep="")
                                                                                                        URL.DrugBank<-paste("<a href=http://www.drugbank.ca/drugs/", ID.DrugBank,">", ID.DrugBank,"</a>",sep="")
                                                                                                    }else{
                                                                                                        
                                                                                                        if(csid1[dbindex]$ExtRef@ds_name=="Comparative Toxicogenomics Database"){
                                                                                                            ID.CTD=csid1[dbindex]$ExtRef@ext_id
                                                                                                            URL.CTD<-paste("<a href=",csid1[dbindex]$ExtRef@ext_url,">",ID.CTD,"</a>",sep="")
                                                                                                        }else{
                                                                                                            
                                                                                                            if(csid1[dbindex]$ExtRef@ds_name=="ACToR: Aggregated Computational Toxicology Resource"){
                                                                                                                if(length(csid1[dbindex]$ExtRef@ext_url)>0){
                                                                                                                    ID.ACToR=csid1[dbindex]$ExtRef@ext_id
                                                                                                                    #print(ID.ACToR)
                                                                                                                    #print(csid1[dbindex]$ExtRef)
                                                                                                                    URL.ACToR<-paste("<a href=",csid1[dbindex]$ExtRef@ext_url,">",ID.ACToR,"</a>",sep="")
                                                                                                                }
                                                                                                            }
                                                                                                        }
                                                                                                    }
                                                                                                }
                                                                                            }
                                                                                        }
                                                                                    }
                                                                                }
                                                                            }
                                                                        }
                                                                    }
                                                                }
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }}
                        }
                        annot_res<-rbind(annot_res,c(mz.val,delta_ppm,info@CSID,urllink1,urllink2,CommonName,MF,SMILES,InChI,InChIKey,AverageMass,MolecularWeight,MonoisotopicMass,
                        NominalMass,ALogP,XLogP, img_url, ID.KEGG,URL.KEGG, ID.HMDB, URL.HMDB, ID.LipidMAPS, URL.LipidMAPS,
                        ID.PubChem,URL.PubChem,ID.MassBank,URL.MassBank,ID.BioCyc,URL.BioCyc, ID.SMPDB, URL.SMPDB,
                        ID.EPA.DSSTox,URL.EPA.DSSTox,ID.EPA.Toxcast,
                        URL.EPA.Toxcast,ID.Pesticide.Common.Names,URL.Pesticide.Common.Names,ID.ChEMBL,URL.ChEMBL,ID.ChEBI,
                        URL.ChEBI,ID.NIST.Chem.WebBook,URL.NIST.Chem.WebBook,ID.WikiPathways,URL.WikiPathways,
                        ID.DrugBank,URL.DrugBank,ID.CTD,URL.CTD,ID.ACToR,URL.ACToR))
                        countid<-countid+1
                        if(countid%%50>0){
                            Sys.sleep(syssleep)			
                        }
                    }
                }
                
            }
            
            
            write.table(annot_res,file="chspider_annot_res.txt",sep="\t",row.names=FALSE,append=TRUE)	
        }#else{
        #	annot_res<-c(mz.val,rep("-",50)) #"-","-","-","-","-","-","-","-","-","-","-","-","-","-","-")
        
    }
    Sys.sleep(syssleep*2)
    
    
    return(annot_res)
    
}
