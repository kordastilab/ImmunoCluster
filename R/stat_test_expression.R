#' @import broom
#' @import reshape2
#'
#' @export
#'
stat_test_expression <- function(
  indata,
  assay = 'scaled',
  grouping = 'group', # The condensing function
  feature = 'condition', # The contrast
  clusterAssign = 'cell_annotation', # The clustering
  test = "wilcox"
){

  plotobj <- data.frame(Cluster = indata@metadata[,clusterAssign],
                        Grouping = indata@metadata[,grouping],
                        as.data.frame(t(as.matrix(assay(indata)))))

  data.median = plotobj %>% group_by( Grouping, Cluster) %>%
    summarize_all(list(median))

  ## Convert the datatframe using Melt
  expr_median_sample_cluster_melt = melt(data.median,
                                         id.vars = c("Grouping", "Cluster"), value.name = "median_expression",
                                         variable.name = "antigen")

  ## Rearange so the rows represent clusters and markers
  expr_median_sample_cluster <- dcast(expr_median_sample_cluster_melt,
                                      Cluster + antigen ~ Grouping, value.var = "median_expression")
  rownames(expr_median_sample_cluster) <- paste0(expr_median_sample_cluster$Cluster,
                                                 "_", expr_median_sample_cluster$antigen)

  # Summarize sampleID to metadata feature
  feature.summary = data.frame(Grouping = indata@metadata[,grouping], Feature = indata@metadata[,feature]) %>%
    distinct()

  ## Eliminate cases with zero expression in all samples
  keep0 <- rowSums(as.data.frame(expr_median_sample_cluster[, as.character(feature.summary$Grouping)]), na.rm = TRUE) > 0
  expr_median_sample_cluster <- expr_median_sample_cluster[keep0, as.character(feature.summary$Grouping)]

  ## Add condition info
  mm = match(colnames(expr_median_sample_cluster), feature.summary$Grouping)


  # Need to recover the original matrix
  grouping_vector = factor(feature.summary$Feature[mm])

  pval_df = NULL

  for(i in 1:length(rownames(expr_median_sample_cluster))){

    obs_test = c()

    ## Find the condition positions in the grouping vector
    feature_variables =  unique(grouping_vector)

    # Count the number fo NAs for each condition
    for(k in 1:length(feature_variables)){

      vals = expr_median_sample_cluster[i,]

      # count the number of not NAs
      idx = which(grouping_vector == feature_variables[k])
      num_obs = sum(!is.na(vals[idx]))

      # If too few obs then add FALSE to vector
      if(num_obs > 2){
        obs_test = c(obs_test, TRUE)
      } else {obs_test = c(obs_test, FALSE)}


    }

    # If obs test contains a FALSE abandon the stat test and return an NA pval
    if(FALSE %in% obs_test){

      # pval is NA as too few observations in at least one condition
      pval = data.frame(p.value = NA)

    } else { # run the stat test

    if(test == "wilcox"){
     pval = tidy(pairwise.wilcox.test(as.numeric(expr_median_sample_cluster[i,]),
                                      g = grouping_vector,
                                      paired = F, correct = F, exact = F))
    } else if(test == "t_test"){
      pval = tidy(pairwise.t.test(as.numeric(expr_median_sample_cluster[i,]),
                                       g = grouping_vector,
                                       paired = F, alternative = "two.sided", pool.sd = F, p.adjust.method = "none"))
    } else{
      print("Error: Test not found, please select available test 'wilcox'")
      return()
    }
    }

    ord = order(grouping_vector)

    vals = cbind(t(expr_median_sample_cluster[i,][ord]), condition = grouping_vector[ord])

    means = as.data.frame(vals) %>% group_by(condition) %>% summarize_all(list(~ mean(., na.rm = TRUE)))

    log_fc = log(means[1,2]/means[2,2])

    out_data = data.frame(mean_2 = means[1,2],
                          mean_2 = means[2,2], logfc = log_fc[1,1])

   colnames(out_data) = c('mean_1', 'mean_2', 'logfc')

   out_data = cbind(out_data, "p_val" = pval$p.value)

    # Add mean1, mean2 and LogFC
    pval_df = rbind(pval_df, out_data)

  }

  expression_df = cbind(expr_median_sample_cluster, pval_df)
  expression_df$padj = p.adjust(expression_df$p_val, method = "BH")
  expression_df = cbind("cluster" = rownames(expression_df), expression_df)
  rownames(expression_df) = NULL

  return(expression_df)


}
