#' @import broom
#' @export
#'
stat_test_clust <- function(
  sct,
  group = 'group',
  clustering = 'cell_annotation',
  feature = 'condition',
  p_val = 'padj',
  threshold = 0.1,
  test = "wilcox"
){

  ## Create the props table
  counts_table = table(sct@metadata[clustering][,1], sct@metadata[,group])

  props_table = t(t(counts_table) / colSums(counts_table)) * 100
  counts = as.data.frame.matrix(counts_table)
  props = as.data.frame.matrix(props_table)

  ggdf = melt(data.frame(cluster = rownames(props), props),
               id.vars = "cluster", value.name = "proportion", variable.name = "sample_id")
  ggdf$cluster = factor(ggdf$cluster)
  ggdf$sample_id = gsub("X", "", ggdf$sample_id)

  # Summarize sampleID to metadata feature
  feature.summary = data.frame(sample_id = sct@metadata[,group], feature = sct@metadata[,feature]) %>%
    distinct()

  ## Add condition info
  mm = match(ggdf$sample_id, feature.summary$sample_id)


  # Need to recover the original matrix
  ggdf$condition = factor(feature.summary$feature[mm])

  # Perform statistical testing
  pop_split_df = split(ggdf, ggdf$cluster)

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
                                        paired = F, alternative = "two.sided", pool.sd = F, p.adjust.method = "none"))

    } else{
      print("Error: Test not found, please select available test 'wilcox'")
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
