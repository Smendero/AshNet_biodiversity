### Microbial Community Helper functions ###

#### NMDS ggplot code function ####

## create datasets for ggplotting

Ord_extract_sites <- function(Ord, characteristics_df=NULL){
      ## Ord: the Ordination object ## Can be either NMDS or RDA
      ## characteristics_df: the characteristics for the sites that went into the creation of the Ordination. 
      ## These must match the order of the data used to create the Ordination
      if("metaMDS" %in% attr(Ord, "class")){
            return(data.frame(Ord$points,
                              characteristics_df))}
      if("rda" %in% attr(Ord, "class")){
            Ord_plot <- plot(Ord)
            return(data.frame(Ord_plot$sites,
                              characteristics_df))
      }
      
      
}


Ord_extract_spec <- function(Ord, spec_characteristics_df=NULL){
      ## Ord: the Ordination object ## Can be either NMDS or RDA
      ## spec_characteristics_df: the characteristics for the species that went into the creation of the Ordination. 
      ## These must match the order of the columns of data used to create the Ordination
      if("metaMDS" %in% attr(Ord, "class")){
            return(data.frame(Ord$species,
                              spec_characteristics_df))}
      if("rda" %in% attr(Ord, "class")){
            Ord_plot <- plot(Ord)
            return(data.frame(Ord_plot$species,
                              spec_characteristics_df))
      }
      
      
}

Ordination_hulls <- function(Ord_ggsites, groups=NULL){
      ## Ord_ggsites: the ggplot compatible data from the Ord_extract_sites function
      ## groups: a character vector containing the characteristics to use for creating hulls (will calculate hulls for all combinations of variables)  
      ## order data
      if(sum(groups %in% colnames(Ord_ggsites)) < length(Ord_ggsites)) {
         errorCondition(message=paste(str_c(groups[!groups %in% colnames(Ord_ggsites)], collapse = ", "), "not in columns of Ord_ggsites"))
      }
   
      ordervec <- do.call(order, as.list(Ord_ggsites[colnames(Ord_ggsites) %in% groups]))
      
      Ord_ggsites <- Ord_ggsites[ordervec,]
      
      ##combine variables into group
      Ord_ggsites$hullgrp <- apply(Ord_ggsites, 1, function(x){str_c(x[names(x) %in% groups], collapse="")})
      
      ## create find_hull function
      find_hull <- function(df) df[chull(df[,1], df[,2]), ]
      
      ## Find hull
      hull <- plyr::ddply(Ord_ggsites, "hullgrp", find_hull) ## hullgrp is grouping variable

}

Ordination_load_extract <- function(Ord){
      ## Ord: an rda class object
      if("rda" %in% attr(Ord, "class")){
            Ord_plot <- plot(Ord)
            return(data.frame(Ord_plot$biplot, Names=row.names(Ord_plot$biplot)))
      }else{
            error("Only takes RDA")
      }
}

Ordination_axis_var <- function(Ord){
      ## Ord: an rda class object
      if("rda" %in% attr(Ord, "class")){
            return(Ord$CCA$eig/Ord$tot.chi*100)
      }else{
            error("Only takes RDA")
      }
}

pairwise_adonis <- function(Orddf, grouping_columns, characteristics_df, method="bray"){
      ## Orddf: a species abundance matrix
      ## grouping_columns: column combinations over which to perform pairwise comparisons
      ## Characteristics_df: contains variables of interest to run adonis, needs to be in same order as species matrix
      ## Does not work for random variables in formula, only currently works for single grouping
      
      ## grab variables and pull list of unique pairings from data
      var_cols <- grouping_columns
      
      if(length(var_cols)==1){
         characteristics_df$unique_sets <- characteristics_df[,var_cols]
         
      }else{
         characteristics_df$unique_sets <- apply(characteristics_df[, colnames(characteristics_df) %in% var_cols], 1, function(x){str_c(x, collapse = "-")})
      }
      
      unique_set_values <- unique(characteristics_df$unique_sets)
      
      pairwise_results <- data.frame()
      
      ##Get unique combinations
      for(i in 1:(length(unique_set_values)-1)){
            for(ii in (i+1):length(unique_set_values)){
                  subset = unique_set_values[c(i, ii)]
            
                  pair_df <- Orddf[characteristics_df$unique_sets %in% subset, ]
                  test_set <- characteristics_df[characteristics_df$unique_sets %in% subset, ]
                  
                  ## run new adonis
                  test_pair = str_c(subset, collapse = ":")
                  
                  new_formula = formula(paste("pair_df", "~", "unique_sets"))
                  pairwise_results <- rbind(pairwise_results, 
                                               data.frame(test_pair=test_pair, data.frame(adonis(new_formula, test_set, method=method)$aov.tab)[1,]))
            }
      }
      return(pairwise_results)
}


pairwise_betadisper <- function(Orddf, grouping_columns, characteristics_df, method="bray"){
   ## Orddf: a species abundance matrix
   ## grouping_columns: column combinations over which to perform pairwise comparisons
   ## Characteristics_df: contains variables of interest to run adonis, needs to be in same order as species matrix
   ## Does not work for random variables in formula, only currently works for single grouping
   
   ## grab variables and pull list of unique pairings from data
   var_cols <- grouping_columns

   if(length(var_cols)==1){
      characteristics_df$unique_sets <- characteristics_df[,var_cols]
      
   }else{
      characteristics_df$unique_sets <- apply(characteristics_df[, colnames(characteristics_df) %in% var_cols], 1, function(x){str_c(x, collapse = "-")})
      
   }
   
   unique_set_values <- unique(characteristics_df$unique_sets)
   
   pairwise_results <- data.frame()
   
   ##Get unique combinations
   for(i in 1:(length(unique_set_values)-1)){
      for(ii in (i+1):length(unique_set_values)){
         subset = unique_set_values[c(i, ii)]
         
         pair_df <- Orddf[characteristics_df$unique_sets %in% subset, ]
         test_set <- characteristics_df$unique_sets[characteristics_df$unique_sets %in% subset]
         
         ## run new adonis
         test_pair = str_c(subset, collapse = ":")
         
         pairwise_results <- rbind(pairwise_results, 
                                   data.frame(test_pair=test_pair, anova(betadisper(vegdist(pair_df, method=method), test_set))[1,]))
      }
   }
   return(pairwise_results)
}

