HMDB.Annotation <-
function(dataA,max.mz.diff=10,num_nodes=2,queryadductlist=c("M+2H","M+H+NH4","M+ACN+2H","M+2ACN+2H","M+H","M+NH4","M+Na","M+ACN+H","M+ACN+Na","M+2ACN+H","2M+H","2M+Na","2M+ACN+H"),
gradienttype="Acetonitrile",mode="pos",outloc){
    
    res<-simpleAnnotation(dataA,max.mz.diff=max.mz.diff,num_nodes=num_nodes,queryadductlist=queryadductlist,
    gradienttype=gradienttype,mode=mode,outloc,db_name="HMDB")
    
    return(res)
}
