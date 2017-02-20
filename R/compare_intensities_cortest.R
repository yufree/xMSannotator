compare_intensities_cortest <-
function(other_feats,y,merge.eval.pvalue){
											cortest_sum=apply(other_feats,1,function(x)
				                                                        {
				                                                        	
																		x<-as.matrix(x)
				                                                        	y<-as.matrix(y)
												yind<-which(y==0)
												xind<-which(x==0)
												naind<-c(yind,xind)
												
				                                                                if(dim(x)[1]>dim(y)[1])
												{
													x<-t(x)
												}
												
												if(dim(y)[1]>dim(x)[1])
												{
													y<-t(y)
												}
												if(length(naind)>0){
												 x<-x[-naind]
												 y<-y[-naind]
				                                                                 
												 }
												 #if(max(x)!=0 & max(y)!=0)
												
													
													if(length(x)>2 & length(y)>2)
													{
																cortest_res=try(cor.test(as.numeric(x),as.numeric(y),method="pearson"),silent=TRUE)
																if (is(cortest_res, "try-error")){
																	cortest_pval<-1
																	cortest_est<-0
																}else{
																cortest_pval=cortest_res$p.value
																cortest_est=cortest_res$estimate
																if(cortest_est<0.1){
																	cortest_pval=1
																}
																}
													}else{
															cortest_pval<-1
													}
													
											      return(cortest_pval)
				                                                        })
return(cortest_sum)
}
