docortest <-
function(n,r){
	#n<-length(x)
t<-r*sqrt((n-2)/(1-r^2))
pvalue<-2*pt(-abs(t),df=n-2)
	
	return(pvalue)
	
}
