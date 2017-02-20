find.Unique.mzs <-
function(dataA, dataB, mz.thresh=10, time.thresh=NA, alignment.tool=NA)
{

	
        data_a<-as.data.frame(dataA)
        data_b<-as.data.frame(dataB)
	data_a<-unique(data_a)
	rm(dataA)
	rm(dataB)
        
        data_b<-unique(data_b)
        
        com_mz_num=1
        unique_mz={}
        ppm_v={}
        rt_v={}

        commat={}

	col.names.dataA=colnames(data_a)
	col.names.dataB=colnames(data_b)

	if(is.na(alignment.tool)==FALSE){
	 if(alignment.tool=="apLCMS")
        {
              sample.col.start=5
        }
        else
        {
              if(alignment.tool=="XCMS")
              {
                    sample.col.start=9
                    col.names.dataA[1]="mz"
                    col.names.dataA[2]="time"
		    col.names.dataB[1]="mz"
                    col.names.dataB[2]="time"
                    colnames(data_a)=col.names.dataA
                    colnames(data_b)=col.names.dataB
              }
	      
	}}else{
                    #stop(paste("Invalid value for alignment.tool. Please use either \"apLCMS\" or \"XCMS\"", sep=""))
		   
		    col.names.dataA[1]="mz"
		   
		    col.names.dataB[1]="mz"
		    if(is.na(time.thresh)==FALSE){
		     col.names.dataA[2]="time"
		     col.names.dataB[2]="time"
		      #print("Using the 1st columns as \"mz\" and 2nd columsn as \"retention time\"")
		     }else{
		      #print("Using the 1st columns as \"mz\"")
		     }
		    colnames(data_a)=col.names.dataA
                    colnames(data_b)=col.names.dataB
	}

       #data_a<-data_a[order(data_a$mz),]
       #data_b<-data_b[order(data_b$mz),]
       data_a<-as.data.frame(data_a)
	data_b<-as.data.frame(data_b)
	colnames(data_a)=col.names.dataA
	colnames(data_b)=col.names.dataB
	
        #create header for the matrix with common features
	if(is.na(time.thresh)==FALSE){
	mznames=c("index.A","mz.data.A", "time.data.A", "index.B","mz.data.B","time.data.B", "time.difference") 
        }else{
	mznames=c("index.A","mz.data.A", "index.B","mz.data.B") 
	}
	
	overlap_res<-find.Overlapping.mzs(dataA=data_a, dataB=data_b, mz.thresh, time.thresh, alignment.tool)
	
	

	if(length(overlap_res$index.A)>0){
	uniqueA<-data_a[-c(overlap_res$index.A),]
	}else{
		uniqueA<-data_a
	}
	
	if(length(overlap_res$index.B)>0){
	uniqueB<-data_b[-c(overlap_res$index.B),]
	}else{
		uniqueB<-data_b
	}
        return(list("uniqueA"=uniqueA,"uniqueB"=uniqueB))
}
