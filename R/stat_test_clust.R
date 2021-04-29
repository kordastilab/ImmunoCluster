#' @rdname stat_test_clust
#'
#' @title Perform statistical test on cluster abundances as a proportion of total cells within a sample across conditions.
#' Generates p adjusted values with "BH" p.adjust method.
#'
#' @param sct a \code{\link[SingleCellExperiment]{SingleCellExperiment}}.
#' @param group the sample_id to generate the calculate the cluster abundance proportions from.
#' @param clustering the clustering to the cluster abundance proportions from.
#' @param feature the metadata condition to test the cluster abundance proportions between.
#' @param test either 'wilcox' calls pairwise.wilcox.test or 't_test' calls pairwise.t.test.
#' @param var_equal sets the pool.sd parameter to apply welch's correction as a parameter of pairwise.t.test
#'
#' @returns a dataframe containing the aubundance proportions per group of each cluster
#' along with pval and BH adjusted pvalue.
#'
#' @import broom
#' @export
#'
stat_test_clust <- function(
  sct,
  group = 'group',
  clustering = 'cell_annotation',
  feature = 'condition',
  test = "wilcox",
  var_equal = TRUE
){

  # count by cluster and sample specified in group
  cell_counts_tbl = table(sct@metadata[clustering][,1], sct@metadata[,group])

  # turn the counts into proportion table
  cell_proportion_tbl = t(t(cell_counts_tbl) / colSums(cell_counts_tbl)) * 100

  # Convert to dataframe
  cell_proportion_df = as.data.frame.matrix(cell_proportion_tbl)

  # Generate cell proportions dataframe
  prop_abundance_df = melt(data.frame(cluster = rownames(cell_proportion_df), cell_proportion_df),
               id.vars = "cluster", value.name = "proportion", variable.name = "sample_id")

  prop_abundance_df$cluster = factor(prop_abundance_df$cluster)
  prop_abundance_df$sample_id = gsub("X", "", prop_abundance_df$sample_id)

  # Summarize sampleID to metadata feature
  feature.summary = data.frame(sample_id = sct@metadata[,group], feature = sct@metadata[,feature]) %>%
    distinct()

  # Add condition info
  match_idx = match(prop_abundance_df$sample_id, feature.summary$sample_id)

  # Need to recover the original matrix
  prop_abundance_df$condition = factor(feature.summary$feature[match_idx])

  # Perform statistical testing
  pop_split_df = split(prop_abundance_df, prop_abundance_df$cluster)

  pval_df = data.frame(cluster = NULL, pval = NULL, padj = NULL)

  for(i in 1:length(pop_split_df)){

    population = unique(pop_split_df[[i]]$cluster)

    if(test == "wilcox"){
      # data = pop_split_df[[i]],
      p_val = tidy(pairwise.wilcox.test(x = pop_split_df[[i]]$proportion,
                                   g = pop_split_df[[i]]$condition,
                                  paired = F, correct = F, exact = F))
    } else if(test == "t_test"){
      # data = pop_split_df[[i]],
      p_val = tidy(pairwise.t.test(x = pop_split_df[[i]]$proportion,
                                        g = pop_split_df[[i]]$condition,
                                        paired = F, alternative = "two.sided", pool.sd = var_equal, p.adjust.method = "none"))

    } else{
      print("Error: Test not found, please select available test 'wilcox' or 't_test'")
      return()
    }

    ord = order(pop_split_df[[i]]$condition)

    vals = as.data.frame(t(pop_split_df[[i]]$proportion[ord]))
    colnames(vals) = pop_split_df[[i]]$sample_id[ord]

    means = pop_split_df[[i]][,c("proportion","condition")] %>% group_by(condition) %>% summarize_all(list(mean))

    log_fc = log(means[1,2]/means[2,2])

    # Add mean1, mean2 and LogFC
    pval_df = rbind(pval_df, data.frame(cluster = population, vals, mean_2 = means[1,2],
                                        mean_2 = means[2,2], logfc = log_fc[1,1], pval = p_val$p.value))


  }

  pval_df$padj = p.adjust(pval_df$pval, method = "BH")

  return(pval_df)
}
