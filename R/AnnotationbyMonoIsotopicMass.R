AnnotationbyMonoIsotopicMass <-
function(inputmassmat,dataA,queryadductlist=c("M+H"),adduct_table,max.mz.diff=10){
	
#inputmassmat<-read.table("~/Documents/Emory/JonesLab/Projects/DOD/PAHs_list_dpj.txt",sep="\t",header=TRUE)

#d2<-unique(d1[,1])

inputmassmat<-inputmassmat[-which(duplicated(inputmassmat[,2])==TRUE),]

mz_search_list<-{}
for(m in 1:dim(inputmassmat)[1]){
	
	mz_search_list<-rbind(mz_search_list,get_mz_by_monoisotopicmass(monoisotopicmass=inputmassmat[m,1],name=as.character(inputmassmat[m,2]),queryadductlist = c("positive")))
	
	
}
	
ghilic<-getVenn(dataA=dataA,name_a="Experimental data",name_b="Search list",dataB=mz_search_list_1,mz.thresh=max.mz.diff,time.thresh=NA,
xMSanalyzer.outloc=outloc,alignment.tool=NA)


mhilic<-merge(mz_search_list,ghilic$common,by.x="mz",by.y="mz.data.B")

mhilic_2<-merge(mhilic,hilic_data,by.x="mz.data.A",by.y="mz")

mhilic_3<-mhilic_2[order(mhilic_2$Name),]

mhilic_3<-mhilic_3[,-c(3,8,9)]

cnames<-colnames(mhilic_3)
cnames[1]<-"mz"
cnames[2]<-"theoretical.mz"
colnames(mhilic_3)<-as.character(cnames)

#write.table(mhilic_3,file="PAH_Hilic_pos_mode_matches.txt",sep="\t",row.names=FALSE)

return(mhilic_3)

}
