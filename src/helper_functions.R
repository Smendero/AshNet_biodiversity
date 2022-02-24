## Helper functions for data analysis ##


Hclust_dimension_reduc <- function(num_dat, cor_thresh, outfile=NULL, method="pearson"){
        ## For use when data is too highly dimensional to run PCA. Finds groups as determined by heirarchical clustering of a correlation matrix correlated at a given threshold
        # outputs a dataframe with new group members and the member with the highest average correlation to the other variables
        cor_mat <- cor(num_dat, method=method)
        
        cor_dist <- dist(cor_mat)
        
        cor_clust <- hclust(cor_dist)
        
        
        ## go down "leaves" of grouped variables 
        ## until pulled single variable that is highly correlated with related variables (.9 or higher).
        finalgroups <- c()
        
        group_num <- 0
        
        captured <- c()
        
        for(i in unique(cor_clust$height)[order(unique(cor_clust$height), decreasing=TRUE)]){
                cor_groups <- cutree(cor_clust, h = i)
                cor_groups <- cor_groups[!names(cor_groups) %in% captured]
                if(length(cor_groups)<=1){
                        break
                }
                group_cor_chk <- data.frame(group = c(), maxavcor = c(), rep_var=c(),  members=c())
                for(g in unique(cor_groups)){
                        temporary_groups <- cor_mat[row.names(cor_mat) %in% names(cor_groups[cor_groups==g]), 
                                                            colnames(cor_mat) %in% names(cor_groups[cor_groups==g])]
                        temporary_groups[temporary_groups==1] <- NA
                        
                        
                        if(length(temporary_groups) > 1){
                                group_cor_chk <- rbind(group_cor_chk,
                                                       data.frame(group=g,
                                                                  maxavcor=max(rowMeans(abs(temporary_groups), na.rm = T)),
                                                                  rep_var = ifelse(!is.null(nrow(temporary_groups[rowMeans(abs(temporary_groups), na.rm=T) == max(rowMeans(abs(temporary_groups), na.rm = T)), ])),
                                                                                   rownames(temporary_groups)[rowMeans(abs(temporary_groups), na.rm=T) == max(rowMeans(abs(temporary_groups), na.rm = T))][1], 
                                                                                   rownames(temporary_groups)[rowMeans(abs(temporary_groups), na.rm=T) == max(rowMeans(abs(temporary_groups), na.rm = T))]), 
                                                                  members = colnames(temporary_groups), 
                                                                  correlation_value = if(!is.null(nrow(temporary_groups[rowMeans(abs(temporary_groups), na.rm=T) == max(rowMeans(abs(temporary_groups), na.rm = T)), ]))){
                                                                          temporary_groups[rowMeans(abs(temporary_groups), na.rm=T) == max(rowMeans(abs(temporary_groups), na.rm = T)), ][1, ]
                                                                  }else{
                                                                          temporary_groups[rowMeans(abs(temporary_groups), na.rm=T) == max(rowMeans(abs(temporary_groups), na.rm = T)), ]  
                                                                  }
                                                       ))
                        }else{
                                group_cor_chk <- rbind(group_cor_chk,
                                                       data.frame(group=g,
                                                                  maxavcor=1,
                                                                  rep_var = names(cor_groups)[cor_groups==g],
                                                                  members = names(cor_groups)[cor_groups==g], 
                                                                  correlation_value = NA))
                        }
                        
                }
                group_cor_chk <- group_cor_chk[group_cor_chk$maxavcor > cor_thresh, ]
                captured <- c(captured, names(cor_groups)[cor_groups %in% group_cor_chk$group])
                
                group_cor_chk$group <- as.numeric(as.factor(group_cor_chk$group)) + group_num
                
                # use factorization to rename group numbers
                finalgroups <- rbind(finalgroups, group_cor_chk)
                
                if(nrow(finalgroups)>0){group_num <- max(finalgroups$group)}
        }
        
        finalgroups$correlation_value[is.na(finalgroups$correlation_value)] <- 1
        
        ## write out the finalgroups to a CSV for record
        if(!is.null(outfile)){write.csv(finalgroups, outfile)}
        return(finalgroups)
        
}

attr(Hclust_dimension_reduc, "help") <- "Finds groups of correlated variables using cor and hclust and selects the variable that has the highest average correlation (above a provided threshold) to serve as a representative variable for the group."
