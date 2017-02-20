do_minmax <-
function(data_m, newmin, newmax)
	{
		data_m=apply(data_m,2,function(x){
		minx=min(x,na.rm=TRUE)
		maxx=max(x,na.rm=TRUE)
		if(minx!=maxx)
		{
			(((x-min(x,na.rm=TRUE))/(max(x,na.rm=TRUE)-min(x,na.rm=TRUE)))*(newmax-newmin))+newmin
		}else
		{
			(x-min(x,na.rm=TRUE))/(max(x,na.rm=TRUE)-min(x,na.rm=TRUE)+1)
		}
		})
		return(data_m)

	}
