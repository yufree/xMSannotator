multilevelannotationstep4 <- function(outloc, max_diff_rt, 
    adduct_weights = NA, filter.by = NA, min_ions_perchem = 1, 
    boostIDs = NA, max_isp = 5, dbAllinf = NA) {
    
    setwd(outloc)
    
    outloc2 <- paste(outloc, "/stage3/", sep = "")
    setwd(outloc)
    
    # chemscoremat_highconf<-read.table('Stage3.txt',sep='\t',header=TRUE)
    # #read.table('Stage3.txt',sep='\t',header=TRUE)
    
    chemscoremat_highconf <- read.csv("Stage3.csv")
    
    chemscoremat_highconf <- as.data.frame(chemscoremat_highconf)
    chemscoremat_highconf$mz <- as.numeric(chemscoremat_highconf$mz)
    
    if (FALSE) {
        chemscoremat_highconf$Name <- gsub(chemscoremat_highconf$Name, 
            pattern = "[\\\"']", replacement = "")
        
        chemscoremat_highconf$Formula <- gsub(chemscoremat_highconf$Formula, 
            pattern = "[\\\"']", replacement = "")
        
        chemscoremat_highconf$Adduct <- gsub(chemscoremat_highconf$Adduct, 
            pattern = "[\\\"']", replacement = "")
    }
    cnames <- colnames(chemscoremat_highconf)
    
    cnames <- gsub(cnames, pattern = ".x", replacement = "")
    
    colnames(chemscoremat_highconf) <- cnames
    
    chemids <- chemscoremat_highconf$chemical_ID
    
    chemids <- unique(chemids)
    
    chemscoremat_conf_levels <- rep("High", length(chemids))
    
    # print('boostIDs') print(boostIDs)
    
    data(adduct_table)
    if (is.na(adduct_weights) == TRUE) {
        data(adduct_weights)
        adduct_weights1 <- matrix(nrow = 2, ncol = 2, 
            0)
        adduct_weights1[1, ] <- c("M+H", 1)
        adduct_weights1[2, ] <- c("M-H", 1)
        adduct_weights <- as.data.frame(adduct_weights1)
        colnames(adduct_weights) <- c("Adduct", "Weight")
    }
    adduct_table <- adduct_table[order(adduct_table$Adduct), 
        ]
    chemscoremat_conf_levels1 <- lapply(1:length(chemids), 
        function(c) {
            
            cur_chemid <- chemids[c]
            
            curdata <- chemscoremat_highconf[which(chemscoremat_highconf$chemical_ID == 
                cur_chemid), ]
            
            bool_check = 1
            
            
            curdata <- curdata[order(curdata$Adduct), 
                ]
            
            if ((is.na(filter.by) == FALSE) && (bool_check == 
                1)) {
                
                check_adduct <- which(curdata$Adduct %in% 
                  filter.by)
                if (length(check_adduct) > 0) {
                  bool_check = 1
                } else {
                  bool_check = 0
                }
                
            }
            
            if (bool_check == 1) {
                
                # print(cur_chemid) print(curdata)
                final_res <- get_confidence_stage4(curdata, 
                  max_diff_rt, adduct_weights = adduct_weights, 
                  filter.by = filter.by, max_isp = max_isp, 
                  min_ions_perchem = min_ions_perchem)
                
                Confidence <- 0
                # print(final_res)
                if (final_res != "None") {
                  if (is.na(final_res[1, 1]) == FALSE) {
                    
                    
                    Confidence <- as.numeric(as.character(final_res[, 
                      1]))
                    
                    curdata <- final_res  #[,-c(1)]
                    
                    # print(final_res)
                    if (Confidence < 2) {
                      
                      if (length(which(curdata$Adduct %in% 
                        adduct_weights[which(as.numeric(adduct_weights[, 
                          2]) > 0), 1])) > 0) {
                        
                        if (curdata$score > 10) {
                          
                          mnum <- max(as.numeric(as.character(adduct_weights[which(adduct_weights[, 
                            1] %in% curdata$Adduct), 
                            2])))[1]
                          curdata <- curdata[which(curdata$Adduct %in% 
                            adduct_weights[which(as.numeric(as.character(adduct_weights[, 
                              2])) >= mnum), 1]), 
                            ]
                          Confidence <- 2
                        }
                      }
                    }
                    
                  }
                }
            } else {
                Confidence <- 0
                if (length(which(curdata$Adduct %in% 
                  adduct_weights[, 1])) > 0) {
                  
                  if (curdata$score >= 10) {
                    
                    # curdata<-curdata[which(curdata$Adduct%in%adduct_weights[,1]),]
                    mnum <- max(as.numeric(as.character(adduct_weights[which(adduct_weights[, 
                      1] %in% curdata$Adduct), 2])))[1]
                    # curdata<-curdata[which(curdata$Adduct%in%adduct_weights[which(as.numeric(as.character(adduct_weights[,2]))>=mnum),1]),]
                    if (length(which(curdata$Adduct %in% 
                      filter.by)) > 0) {
                      curdata <- curdata[which(curdata$Adduct %in% 
                        filter.by), ]
                      Confidence <- 2
                    }
                    
                  }
                }
                
            }
            
            if (nrow(curdata) > 1) {
                if (curdata$score < 10) {
                  if (length(unique(curdata$Adduct)) < 
                    2) {
                    
                    Confidence <- 0
                  } else {
                    if (FALSE) {
                      if (Confidence < 2) {
                        
                        if (length(which(curdata$Adduct %in% 
                          adduct_weights[which(adduct_weights[, 
                            2] > 1), 1])) > 0) {
                          
                          if (curdata$score > 10) {
                            # Confidence<-2
                            
                            mnum <- max(as.numeric(as.character(adduct_weights[which(adduct_weights[, 
                              1] %in% curdata$Adduct), 
                              2])))[1]
                            curdata <- curdata[which(curdata$Adduct %in% 
                              adduct_weights[which(as.numeric(as.character(adduct_weights[, 
                                2])) >= mnum), 1]), 
                              ]
                            Confidence <- 2
                            
                          }
                        }
                      }
                    }
                  }
                }
            }
            if (cur_chemid %in% boostIDs) {
                
                Confidence <- 4
            }
            
            curdata <- cbind(Confidence, curdata)
            curdata <- as.data.frame(curdata)
            
            
            return(curdata)
        })
    # chemids<-c('HMDB00277','HMDB00222','HMDB00043')
    
    # stopCluster(cl)
    
    # save(list=c('chemscoremat_conf_levels1','chemids'),file='stage4conf_levels.Rda')
    
    # save(chemscoremat_conf_levels1,file='stage4conf_levels.Rda')
    
    
    chemscoremat_conf_levels <- ldply(chemscoremat_conf_levels1, 
        rbind)  #unlist(chemscoremat_conf_levels)
    
    chemscoremat_highconf <- as.data.frame(chemscoremat_highconf)
    
    
    chemscoremat_conf_levels <- as.data.frame(chemscoremat_conf_levels)
    
    
    # chemconf_levels<-cbind(chemscoremat_conf_levels,chemids)
    # chemconf_levels<-as.data.frame(chemconf_levels)
    # colnames(chemconf_levels)<-c('Confidence','chemids')
    
    
    # write.table(chemconf_levels,file='confidence_levels_chemicals.txt',sep='\t',row.names=FALSE)
    
    curated_res <- chemscoremat_conf_levels  # merge(chemconf_levels,chemscoremat_highconf,by.x='chemids',by.y='chemical_ID')
    cnames <- colnames(curated_res)
    cnames[1] <- "Confidence"
    colnames(curated_res) <- as.character(cnames)
    
    curated_res_isp_check <- gregexpr(text = curated_res$Adduct, 
        pattern = "(_\\[(\\+|\\-)[0-9]*\\])")
    
    isp_ind <- which(curated_res_isp_check > 0)
    
    if (length(isp_ind) > 0) {
        
        formula_vec <- curated_res$Formula  #aregexpr(text=curated_res$Formula,pattern='(_\\[(\\+|\\-)[0-9]*\\])')
        formula_vec <- gsub(x = formula_vec, pattern = "[a-z|+|-|0-9]*_", 
            replacement = "", ignore.case = T)
        # print(formula_vec) print(formula_vec[isp_ind])
        temp_adducts <- as.character(paste("M", "_", 
            formula_vec[isp_ind], sep = ""))
        curated_res$Adduct <- as.character(curated_res$Adduct)
        # print(curated_res$Adduct) print(temp_adducts)
        # print(curated_res$Adduct[isp_ind])
        curated_res$Adduct[isp_ind] <- as.character(temp_adducts)
        # print(curated_res$Adduct[isp_ind])
    }
    
    
    # write.table(curated_res,file='confidence_levels_chemicals.txt',sep='\t',row.names=FALSE)
    
    cnames <- colnames(curated_res)
    
    # curated_res<-cbind(curated_res[,3],curated_res[,2],curated_res[,-c(2:3)])
    
    
    # cnames[1]<-'chemical_ID'
    # cnames[2]<-'Score_category'
    # colnames(curated_res)<-as.character(cnames)
    
    # curated_res<-curated_res[order(curated_res[,2],curated_res[,1],curated_res[,3],decreasing=TRUE),]
    
    # curated_res<-curated_res[order(curated_res$Confidence,curated_res$chemical_ID,curated_res$score,curated_res$Adduct,decreasing=TRUE),]
    
    # curated_res<-curated_res[order(curated_res$Confidence,curated_res$chemical_ID,curated_res$score,curated_res$mean_int_vec,decreasing=TRUE),]
    
    
    
    # outloc3<-paste(outloc,'/stage4/',sep='')
    
    outloc3 <- outloc
    dir.create(outloc3)
    setwd(outloc3)
    
    # print(curated_res[1:2,])
    
    # curated_res<-curated_res[,-c(15)]
    # print(curated_res[1:2,])
    
    
    if (FALSE) {
        d2 <- read.table("Stage2.txt", sep = "\t", 
            header = TRUE)
        
        d3pres <- d2[-which(d2$mz %in% curated_res$mz), 
            ]
        d3pres <- d3pres[-which(d3pres$chemical_ID %in% 
            curated_res$chemical_ID), ]
        d3pres <- d3pres[which(d3pres$Adduct %in% 
            adduct_weights[, 1]), ]
        
        conf_vec <- rep(1, length(d3pres[, 1]))
        score_vec <- rep(1, length(d3pres[, 1]))
        d6pres <- cbind(conf_vec, d3pres[, 5], score_vec, 
            d3pres[, 11], d3pres[, 1:4], d3pres[, 
                6:10], d3pres[, 13:14])
        
        print(curated_res[1:2, 1:15])
        print(d6pres[1:2, 1:15])
        print(dim(d6pres))
        print(dim(curated_res))
        
        curated_res <- curated_res[, c(1:15)]
        
        colnames(d6pres) <- colnames(curated_res)
        curated_res <- rbind(curated_res, d6pres)
    }
    curated_res <- as.data.frame(curated_res)
    
    # curated_res<-curated_res[order(curated_res$Confidence,curated_res$chemical_ID,curated_res$score,curated_res$mean_int_vec,decreasing=TRUE),]
    
    # curated_res<-curated_res[order(curated_res$Confidence,curated_res$score,curated_res$mean_int_vec,decreasing=TRUE),]
    # curated_res<-curated_res[order(curated_res$Confidence,curated_res$score,decreasing=TRUE),]
    t2 <- table(curated_res$mz)
    
    same1 <- which(t2 == 1)
    
    uniquemz <- names(same1)
    
    curated_res$MatchCategory = rep("Multiple", dim(curated_res)[1])
    
    curated_res$MatchCategory[which(curated_res$mz %in% 
        uniquemz)] <- "Unique"
    
    # curated_res<-curated_res[,-c(3,12:14)]
    
    
    cnames1 <- colnames(curated_res)
    
    rmindex1 <- which(cnames1 == "score_level")
    if (length(rmindex1) > 0) {
        curated_res <- curated_res[, -c(rmindex1)]
    }
    
    
    
    curated_res <- as.data.frame(curated_res)
    
    curated_res$mz <- as.numeric(as.character(curated_res$mz))
    
    # curated_res_temp<-curated_res[,c(5,8)]
    curated_res$theoretical.mz <- as.numeric(as.character(curated_res$theoretical.mz))
    
    curated_res_temp <- curated_res[, c("mz", "theoretical.mz")]
    
    
    curated_res_temp <- apply(curated_res_temp, 1, 
        as.numeric)
    curated_res_temp <- t(curated_res_temp)
    curated_res_temp <- as.data.frame(curated_res_temp)
    
    delta_ppm <- apply(curated_res_temp, 1, function(x) {
        
        ppmerror = 10^6 * abs(x[2] - x[1])/(x[2])
        return(ppmerror)
    })
    delta_ppm <- round(delta_ppm, 2)
    
    
    
    curated_res <- cbind(curated_res[, 1:8], delta_ppm, 
        curated_res[, 9:dim(curated_res)[2]])
    
    curated_res <- curated_res[order(curated_res$Confidence, 
        decreasing = TRUE), ]
    
    # write.table(curated_res,file='Stage4.txt',sep='\t',row.names=FALSE)
    
    write.csv(curated_res, file = "Stage4.csv", row.names = FALSE)
    
    
    
    if (FALSE) {
        fname = paste("Stage4_annotation_results", 
            sep = "")
        unlink(fname)
        
        HTMLInitFile(filename = fname, Title = "Stage 4 annotation results", 
            outdir = outloc)
        fname = paste(outloc, "/Stage4_annotation_results.html", 
            sep = "")
        HTML(curated_res, file = fname, Border = 1, 
            innerBorder = 1, useCSS = TRUE)
        HTMLEndFile(file = fname)
    }
    
    # print(head(curated_res))
    
    curated_res <- as.data.frame(curated_res)
    
    print("Stage 4 confidence level distribution for unique chemical/metabolite IDs")
    print(table(curated_res$Confidence[-which(duplicated(curated_res$chemical_ID) == 
        TRUE)]))
    
    
    # print(table(curated_res$Confidence))
    return(curated_res)
    
    
}
