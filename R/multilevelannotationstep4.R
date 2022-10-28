multilevelannotationstep4 <- function(outloc,
                                      max.mz.diff = 5,
                                      max.rt.diff = 30,
                                      adduct_weights = NA,
                                      filter.by = NA,
                                      min_ions_perchem = 1,
                                      boostIDs = NA,
                                      max_isp = 5,
                                      dbAllinf = NA,
                                      num_nodes = 2) {
  setwd(outloc)
  
  
  max_diff_rt = max.rt.diff
  
  # chemscoremat_highconf<-read.table('Stage3.txt',sep='\t',header=TRUE)
  # #read.table('Stage3.txt',sep='\t',header=TRUE)
  
  chemscoremat_highconf <- read.csv("Stage3.csv")
  
  chemscoremat_highconf <- as.data.frame(chemscoremat_highconf)
  chemscoremat_highconf$mz <- as.numeric(chemscoremat_highconf$mz)
  
  cnames <- colnames(chemscoremat_highconf)
  
  cnames <- gsub(cnames, pattern = ".x", replacement = "")
  
  colnames(chemscoremat_highconf) <- cnames
  
  chemids <- chemscoremat_highconf$chemical_ID
  
  chemids <- unique(chemids)
  
  
  
  chemscoremat_conf_levels <- rep("High", length(chemids))
  
  
  data(adduct_table)
  if (is.na(adduct_weights)[1] == TRUE) {
    data(adduct_weights)
    adduct_weights1 <- matrix(nrow = 2, ncol = 2, 0)
    adduct_weights1[1, ] <- c("M+H", 1)
    adduct_weights1[2, ] <- c("M-H", 1)
    adduct_weights <- as.data.frame(adduct_weights1)
    colnames(adduct_weights) <- c("Adduct", "Weight")
  }
  adduct_table <- adduct_table[order(adduct_table$Adduct), ]
  
  
  chemscoremat_conf_levels <- {
    
  }
  
  # cl<-makeSOCKcluster(num_nodes)
  
  # if(length(chemids)>1000)
  cur_chemid <- chemids
  curdata <-
    chemscoremat_highconf[which(chemscoremat_highconf$chemical_ID ==
                                  cur_chemid), ]
  
  bool_check = 1
  
  
  curdata <- curdata[order(curdata$Adduct), ]
  tttt <- (is.na(filter.by) == FALSE) && (bool_check ==
                                            1)
  if (tttt[1]) {
    check_adduct <- which(curdata$Adduct %in%
                            filter.by)
    if (length(check_adduct) > 0) {
      bool_check = 1
    } else {
      bool_check = 0
    }
    
    
    if (bool_check == 1) {
      # print(cur_chemid) print(curdata)
      final_res <- get_confidence_stage4(
        curdata,
        max_diff_rt,
        adduct_weights = adduct_weights,
        filter.by = filter.by,
        max_isp = max_isp,
        min_ions_perchem = min_ions_perchem
      )
      
      Confidence <- 0
      # print(final_res)
      if ((final_res[1] != "None")[1]) {
        if (is.na(final_res[1, 1])[1] == FALSE) {
          Confidence <- as.numeric(as.character(final_res[,
                                                          1]))
          
          curdata <- final_res  #[,-c(1)]
          rm(final_res)
          if (Confidence[1] < 2) {
            if (length(which(curdata$Adduct %in%
                             adduct_weights[which(as.numeric(adduct_weights[,
                                                                            2]) > 0), 1])) > 0) {
              if (curdata$score[1] > 10) {
                mnum <-
                  max(as.numeric(as.character(adduct_weights[which(adduct_weights[,
                                                                                  1] %in% curdata$Adduct), 2])))[1]
                curdata <- curdata[which(curdata$Adduct %in%
                                           adduct_weights[which(as.numeric(as.character(adduct_weights[,
                                                                                                       2])) >= mnum), 1]), ]
                Confidence <- 2
              }
            }
          }
          
        }
      }
    } else {
      Confidence <- 0
      if (length(which(curdata$Adduct %in% adduct_weights[,
                                                          1])) > 0) {
        if (curdata$score[1] >= 10) {
          # curdata<-curdata[which(curdata$Adduct%in%adduct_weights[,1]),]
          mnum <-
            max(as.numeric(as.character(adduct_weights[which(adduct_weights[,
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
      if (curdata$score[1] < 10) {
        if (length(unique(curdata$Adduct)) < 2) {
          Confidence <- 0
        }
      }
    }
    curdata <- cbind(Confidence, curdata)
    curdata <- as.data.frame(curdata)
    
    curdata <- curdata[, c("Confidence", "chemical_ID")]
    chemscoremat_conf_levels <- unique(curdata)
    # chemids<-c('HMDB00277','HMDB00222','HMDB00043')
    
    # stopCluster(cl)
    # save(list=c('chemscoremat_conf_levels1','chemids'),file='stage4conf_levels.Rda')
    
    chemscoremat_conf_levels <-
      as.data.frame(chemscoremat_conf_levels)
    # save(chemscoremat_conf_levels2,file='stage4conf_levelsB.Rda')
    # print(mem_used())
    
    # stop('done saving')
    chemscoremat_conf_levels <-
      chemscoremat_conf_levels  #ldply(chemscoremat_conf_levels2,rbind) #unlist(chemscoremat_conf_levels)
    # rm(chemscoremat_conf_levels2)
    # chemscoremat_conf_levels<-rbind(chemscoremat_conf_levels,chemscoremat_conf_levels_temp)
  }
  
  # chemscoremat_highconf<-as.data.frame(chemscoremat_highconf)
  
  chemscoremat_highconf <- unique(chemscoremat_highconf)
  chemscoremat_conf_levels <-
    as.data.frame(chemscoremat_conf_levels)  #[,c(1,3:13,15:16)])
  # print(dim(chemscoremat_conf_levels))
  # save(chemscoremat_conf_levels,file='chemscoremat_conf_levels.Rda')
  # save(chemscoremat_highconf,file='chemscoremat_highconf.Rda')
  
  # chemconf_levels<-cbind(chemscoremat_conf_levels,chemids)
  # chemconf_levels<-as.data.frame(chemconf_levels)
  # colnames(chemconf_levels)<-c('Confidence','chemids')
  # write.table(chemconf_levels,file='confidence_levels_chemicals.txt',sep='\t',row.names=FALSE)
  
  # curated_res<-cbind(chemscoremat_conf_levels[,1],chemscoremat_highconf)
  curated_res <-
    merge(chemscoremat_conf_levels, chemscoremat_highconf,
          by = "chemical_ID")
  
  cnames <- colnames(curated_res)
  # cnames[3]<-'score' #'Confidence'
  colnames(curated_res) <- as.character(cnames)
  
  rm(chemscoremat_highconf)
  
  curated_res_isp_check <- gregexpr(text = curated_res$Adduct,
                                    pattern = "(_\\[(\\+|\\-)[0-9]*\\])")
  
  isp_ind <- which(curated_res_isp_check > 0)
  
  
  # write.table(curated_res,file='confidence_levels_chemicals.txt',sep='\t',row.names=FALSE)
  
  cnames <- colnames(curated_res)
  
  
  outloc3 <- outloc
  suppressWarnings(dir.create(outloc3))
  setwd(outloc3)
  
  curated_res <- as.data.frame(curated_res)
  
  
  
  cnames1 <- colnames(curated_res)
  
  
  
  
  curated_res <- as.data.frame(curated_res)
  
  curated_res$mz <- as.numeric(as.character(curated_res$mz))
  
  
  curated_res$theoretical.mz <-
    as.numeric(as.character(curated_res$theoretical.mz))
  
  curated_res_temp <- curated_res[, c("mz", "theoretical.mz")]
  
  
  curated_res_temp <- apply(curated_res_temp, 1, as.numeric)
  curated_res_temp <- t(curated_res_temp)
  curated_res_temp <- as.data.frame(curated_res_temp)
  
  delta_ppm <- apply(curated_res_temp, 1, function(x) {
    ppmerror = 10 ^ 6 * abs(x[2] - x[1]) / (x[2])
    return(ppmerror)
  })
  delta_ppm <- round(delta_ppm, 2)
  
  
  
  curated_res <-
    cbind(curated_res[, 1:8], delta_ppm, curated_res[,
                                                     9:dim(curated_res)[2]])
  
  curated_res <- curated_res[order(curated_res$Confidence,
                                   decreasing = TRUE), ]
  
  if (is.na(boostIDs)[1] == FALSE) {
    cnames_boost <- colnames(boostIDs)
    
    if (length(cnames_boost) > 1) {
      curated_res_mzrt <- curated_res[, c("mz", "time")]
      validated_mzrt <- boostIDs[, c("mz", "time")]
      
      ghilicpos <- getVenn(
        curated_res_mzrt,
        name_a = "exp",
        validated_mzrt,
        name_b = "boost",
        mz.thresh = max.mz.diff,
        time.thresh = max_diff_rt,
        alignment.tool = NA,
        xMSanalyzer.outloc = getwd(),
        use.unique.mz = FALSE,
        plotvenn = FALSE
      )
      
      save(ghilicpos, file = "ghilicpos.Rda")
      
      g1 <- ghilicpos$common
      rm(ghilicpos)
      
      if (is.na(max_diff_rt)[1] == FALSE) {
        # g1<-g1[order(g1$index.B,g1$time.difference),]
        
        t1 <- table(g1$index.B)
        
        ind_names <- names(t1)
        
        parent_bad_ind <- {
          
        }
        
        
        
        # dup_index_g1<-parent_bad_ind
        # #which(duplicated(g1$index.B)==TRUE)
        # if(length(dup_index_g1)>0){ g1<-g1[-dup_index_g1,] }
        
      }
      
      t1 <-
        table(curated_res$Confidence, curated_res$chemical_ID)
      cnames <- colnames(t1)
      cnames <- cnames[which(cnames %in% boostIDs$ID)]
      
      good_ind_1 <- {
        
      }
      
      for (ind2 in 1:dim(g1)[1]) {
        temp_ind1 <- g1$index.A[ind2]
        temp_ind2 <- g1$index.B[ind2]
        
        ttt <- curated_res$chemical_ID[temp_ind1] %in%
          boostIDs$ID[temp_ind2]
        if (ttt[1]) {
          good_ind_1 <- c(good_ind_1, g1$index.A[ind2])
        }
      }
      
      overlap_mz_time_id <-
        good_ind_1  #which(bool_vec3==2) #mz_time_index #which(curated_res$chemical_ID%in%good_ids)
      
      curated_res$Confidence[overlap_mz_time_id] <- 4
      curated_res$score[overlap_mz_time_id] <-
        curated_res$score[overlap_mz_time_id] *
        100
      t1 <- table(curated_res$Confidence[overlap_mz_time_id],
                  curated_res$chemical_ID[overlap_mz_time_id])
      
      cnames1 <- colnames(t1)
      cnames2 <- cnames1[which(t1 > 0)]
      good_ind <- {
        
      }  #which(curated_res$chemical_ID%in%cnames2)
      if (length(good_ind) > 0) {
        curated_res$Confidence[good_ind] <- 4
        curated_res$score[good_ind] <-
          curated_res$score[good_ind] *
          100
      }
      
    } else {
      good_ind <- which(curated_res$chemical_ID %in%
                          boostIDs)
      if (length(good_ind) > 0) {
        curated_res$Confidence[good_ind] <- 4
        curated_res$score[good_ind] <-
          curated_res$score[good_ind] *
          100
      }
    }
    
  }
  t2 <- table(curated_res$mz)
  
  same1 <- which(t2 == 1)
  
  uniquemz <- names(same1)
  
  curated_res$MatchCategory = rep("Multiple", dim(curated_res)[1])
  
  curated_res$MatchCategory[which(curated_res$mz %in% uniquemz)] <-
    "Unique"
  
  # write.table(curated_res,file='Stage4.txt',sep='\t',row.names=FALSE)
  
  write.csv(curated_res, file = "Stage4.csv", row.names = FALSE)
  
  
  # print(head(curated_res))
  
  curated_res <- as.data.frame(curated_res)
  curated_res <- curated_res[order(curated_res$Confidence,
                                   decreasing = TRUE), ]
  
  print("Stage 4 confidence level distribution for unique chemical/metabolite IDs")
  print(table(curated_res$Confidence[-which(duplicated(curated_res$chemical_ID) ==
                                              TRUE)]))
  
  print("Stage 4 confidence level distribution for unique chemical/metabolite formulas")
  print(table(curated_res$Confidence[-which(duplicated(curated_res$Formula) ==
                                              TRUE)]))
  
  
  
  # print(table(curated_res$Confidence))
  return(curated_res)
  
  
}
