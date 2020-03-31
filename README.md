immunoCluster
================
James Opzoomer, Kevin Blighe, Jessica Timms
2020-03-31

## 1\. Introduction to immunoCluster

The ImmunoCluster package provides a complete toolkit for carry out
immune profiling from High-dimensional mass (CyTOF) and flow cytometry.
The package features standardized data infrastructure, using the
SingleCellExperiment class object, interactive visualization tools and
convenient implementations for popular algorithms designed for a
non-specialist. Provided below is a simple walthrough exploring some of
the packages functionality.

## 2\. Main functionalities of immunoCluster

**NOTE: THIS PACKAGE IS STILL UNDER DEVELOPMENT AND SO SOME OF THE
FUNCTIONALITY IS NOT FULLY TESTED**

# Installation

## 1\. Install the package from github

Install the development version of the package:

``` r
  # Install devtools for github installation if not present
  require(devtools)
  
  # Install the ImmunoCluster package
  devtools::install_github("kordastilab/ImmunoCluster")
  
  # May need to run this line if an R build version error while installing appears
  # Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS="true")
```

## 2\. Load the package and the dependancies into the R session

``` r
  library(immunoCluster)

  # Load dependencies
  require(flowCore)
  library(stringr)
  library(cowplot)
  library(RColorBrewer)
  library(Biobase)
  library(ConsensusClusterPlus)
  library(FlowSOM)
  library(limma)
  library(ggrepel)
  library(reshape2)
  library(Rtsne)
  library(dplyr)
  library(pheatmap)
```

# Walkthrough: Immune profiling PBMC from simple experimental setup

## 1\. Introduction

In this walkthrough we will start with FCS data files to examine changes
in immune populations from PBMC. The samples were taken from patients
following hematopoetic reconsitution in leukemia patients following BMT.
The dataset represents 15 individuals, sampled at multiple time points
after BMT, for a total of 28 samples. Of these patients, a small subset
suffered from graft versus host disease (GvHD, n = 3), whereas most
other patients did not experience such complications.

The data was originally reported in [Comprehensive Immune Monitoring of
Clinical Trials to Advance Human
Immunotherapy](https://www.sciencedirect.com/science/article/pii/S2211124719308228).

The raw data is available
[here](http://flowrepository.org/id/FR-FCM-Z244). We have performed
gating cleanup on a subset of the dataset which can be dowloaded from
zenodo. \[Insert zenodo button\]

## 2\. Import data

``` r
# Insert zenodo link
```

``` r
# Insert zenodo RDS link to replace local
sample_metadata_df = readRDS(file = "../github_readme/rds_objects/sample_metadata.rds")

head(sample_metadata_df)
```

    ## # A tibble: 6 x 5
    ##   file_name           condition patient_id sample_id day_id
    ##   <chr>               <fct>     <chr>      <chr>     <fct> 
    ## 1 BMT01_D30_clean.fcs None      P1         S1        D30   
    ## 2 BMT01_D90_clean.fcs None      P1         S2        D90   
    ## 3 BMT02_D30_clean.fcs None      P2         S3        D30   
    ## 4 BMT02_D90_clean.fcs None      P2         S4        D90   
    ## 5 BMT08_D30_clean.fcs None      P8         S15       D30   
    ## 6 BMT09_D30_clean.fcs None      P9         S16       D30

``` r
# Insert zenodo RDS link to replace local
panel_metadata_df = readRDS(file = "../github_readme/rds_objects/panel_metadata.rds")

head(panel_metadata_df)
```

    ## # A tibble: 6 x 4
    ##   metal     antigen  clustering activation
    ##   <chr>     <chr>         <dbl>      <dbl>
    ## 1 BCKG190Di blank_00          0          0
    ## 2 Ba138Di   blank_01          0          0
    ## 3 Bi209Di   CD16              1          0
    ## 4 Ce140Di   beads             0          0
    ## 5 Ce142Di   beads_1           0          0
    ## 6 Center    blank_8           0          0

``` r
# Create marker grouping designations
clustering_markers = panel_metadata_df$antigen[panel_metadata_df$clustering ==1]
activation_markers = panel_metadata_df$antigen[panel_metadata_df$activation == 1]
col_names = panel_metadata_df$antigen[panel_metadata_df$clustering == 1 | panel_metadata_df$activation == 1 ]

# Select columns to dicard from fcs file
discard = setdiff(panel_metadata_df$antigen, c(clustering_markers, activation_markers))
queries = str_detect(panel_metadata_df$antigen, paste(discard, collapse = "|"))

if(length(discard) > 0){colsDiscard = panel_metadata_df$metal[queries]
}else{colsDiscard = ''}
```

``` r
filelist = list.files(
      path = "../github_readme/fcs_data",
    pattern = "*.fcs|*.FCS",
    full.names = TRUE)

metadata = data.frame(
      file = filelist,
      group = sample_metadata_df$sample_id,
      condition = sample_metadata_df$condition,
      patient_id = sample_metadata_df$patient_id,
      day_id = sample_metadata_df$day_id,
      row.names = filelist,
      stringsAsFactors = FALSE)

# Re-jig the accept discard columns to get a complete limited dataset.
require(SingleCellExperiment)
  sce_gvhd = processFCS(
    files = filelist,
    metadata = metadata,
    transformation = TRUE,
    filter = FALSE,
    transFun = function (x) asinh(x),
    downsampleVar = NULL, 
    downsample = NULL, # Downsample to n cells per downsample_grouping
    downsample_grouping = NULL, # Downdample to downsample cells by the specified metadata slot
    newColnames = col_names,
    colsDiscard = colsDiscard)

sce_gvhd
```

Aternately you can download the imported RDS directly to skip the import
step and save time:

``` r
# Insert zenodo RDS link
sce_gvhd = readRDS(file = "../github_readme/rds_objects/sce_gvhd.rds")

sce_gvhd
```

    ## class: SingleCellExperiment 
    ## dim: 32 912813 
    ## metadata(5): file group condition patient_id day_id
    ## assays(1): scaled
    ## rownames(32): CD16 CD152 ... HLA_DR CD56
    ## rowData names(0):
    ## colnames(912813): cell1 cell2 ... cell912812 cell912813
    ## colData names(0):
    ## reducedDimNames(0):
    ## spikeNames(0):
    ## altExpNames(0):

## 2\. Data Exploration

``` r
mds1 = mdsplot(sce_gvhd, 
               feature = "condition", 
               colkey = c(None = 'royalblue', GvHD = 'red2'))

mds2 = mdsplot(sce_gvhd, 
               feature = "day_id", 
               colkey = c(D30  = 'darkorange1', D90 = 'darkgreen'))

plot_grid(mds1, mds2,
    labels = c('A','B'),
    ncol = 2, align = "l", label_size = 20)
```

![](README_files/figure-gfm/mds_plot-1.png)<!-- -->

``` r
feature = medianHeatmap(sce_gvhd, grouping = "group", feature = "condition", feature_cols = c('royalblue', 'red2'))
```

![](README_files/figure-gfm/sample%20heatmap-1.png)<!-- -->

## 3\. Dimensionality reduction and clustering

Run the UMAP on selected lineage markers, downampling to 2000 cells per
sample. tSNE is also supported.

``` r
require(umap)
# Run UMAP and store in sce object
sce_gvhd = performUMAP(sce_gvhd, downsample = 1000, useMarkers = clustering_markers)

# Run tSNE and store in sce object
sce_gvhd = performTSNE(sce_gvhd, downsample = 1000, useMarkers = clustering_markers)
```

``` r
expression_markers =  c('CD3', 'CD4', 'CD8a', 'CD11b', 'CD19', 'CD56')

exp_plot_umap = markerExpression(sce_gvhd,
    markers = expression_markers,
    reducedDim = 'UMAP',
    title = 'Marker exprssion in UMAP space',
    nrow = 1, ncol = 6,
    pointSize = 0.05,
    legendKeyHeight = 1.0,
    legendLabSize = 14,
    stripLabSize = 20,
    axisLabSize = 18,
    titleLabSize = 24,
    captionLabSize = 22)

exp_plot_tsne = markerExpression(sce_gvhd,
    markers = expression_markers,
    reducedDim = 'Rtsne',
    title = 'Marker exprssion in tSNE space',
    nrow = 1, ncol = 6,
    pointSize = 0.05,
    legendKeyHeight = 1.0,
    legendLabSize = 14,
    stripLabSize = 20,
    axisLabSize = 18,
    titleLabSize = 24,
    captionLabSize = 22)

plot_grid(exp_plot_umap, exp_plot_tsne,
    labels = c('A','B'),
    nrow = 2, align = "l", label_size = 24)
```

![](README_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

Run consensus method of FLowSOM and ConsensusClusterPlus on all cells
using a selection of markers, with a final desired cluster number of uo
to 10.

``` r
sce_gvhd = runFlowSOM(sce_gvhd, k = 10, markers = clustering_markers)
```

    ## Building SOM

    ## Mapping data to SOM

``` r
fsom_k10 = metadataPlot(sce_gvhd,
    colby = 'flowsom_cc_k10',
    title = 'flowSOM k=10',
    legendPosition = 'top',
    legendLabSize = 16,
    axisLabSize = 16,
    titleLabSize = 16,
    subtitleLabSize = 16,
    captionLabSize = 16)

fsom_clusters = unique(sce_gvhd@metadata$flowsom_cc_k10)

fsom_abundance = plotAbundance(sce_gvhd,
              clusters = fsom_clusters,
              clusterAssign = 'flowsom_cc_k10',
              feature = NULL,
              legendPosition = 'none', 
              stripLabSize = 10,
              axisLabSize = 14,
              titleLabSize = 12,
              subtitleLabSize = 10,
              captionLabSize = 18)


plot_grid(fsom_k10, fsom_abundance,
    labels = c('A','B'),
    ncol = 2, align = "l", label_size = 20)
```

![](README_files/figure-gfm/viz%20flowSOM-1.png)<!-- -->

``` r
# With feature
feature = medianHeatmap(sce_gvhd, grouping = "flowsom_cc_k10", feature = NULL, scale_01 = T)
```

![](README_files/figure-gfm/diagnostic%20heatmap-1.png)<!-- -->

``` r
library(scran)
library(tibble)

cytof_markers = findMarkers(assay(sce_gvhd), sce_gvhd@metadata$flowsom_cc_k10, test.type="wilcox", direction="up",lfc=1.5, pval.type="some")

cluster_markers = NULL
top_n = 3
for (i in 1:length(cytof_markers)) {
  
  deg_markers = cytof_markers[[i]][1:top_n,1:2]
  
  deg_markers = rownames_to_column(as.data.frame(deg_markers), var = "marker")
  
  markers = data.frame(cluster = rep(i, top_n), deg_markers)
  
  cluster_markers = rbind(cluster_markers, markers)
  
}

# Top3 DEG markers per cluster
cluster_markers
```

    ##    cluster marker       p.value           FDR
    ## 1        1   CD14  0.000000e+00  0.000000e+00
    ## 2        1   CD33  0.000000e+00  0.000000e+00
    ## 3        1  CD11b  0.000000e+00  0.000000e+00
    ## 4        2  CD11c  0.000000e+00  0.000000e+00
    ## 5        2   CD16  0.000000e+00  0.000000e+00
    ## 6        2  CD123  0.000000e+00  0.000000e+00
    ## 7        3   CD14  0.000000e+00  0.000000e+00
    ## 8        3   CD33  0.000000e+00  0.000000e+00
    ## 9        3    CD3  0.000000e+00  0.000000e+00
    ## 10       4    CD4  0.000000e+00  0.000000e+00
    ## 11       4    CD3  0.000000e+00  0.000000e+00
    ## 12       4  CD127  0.000000e+00  0.000000e+00
    ## 13       5   CD19  0.000000e+00  0.000000e+00
    ## 14       5 HLA_DR 2.572184e-315 4.115495e-314
    ## 15       5   CD16  1.000000e+00  1.000000e+00
    ## 16       6    CD3  0.000000e+00  0.000000e+00
    ## 17       6   Tbet  0.000000e+00  0.000000e+00
    ## 18       6   CD16  4.417942e-25  4.712472e-24
    ## 19       7   CD8a  0.000000e+00  0.000000e+00
    ## 20       7    CD3  0.000000e+00  0.000000e+00
    ## 21       7   Tbet  0.000000e+00  0.000000e+00
    ## 22       8  CD123  0.000000e+00  0.000000e+00
    ## 23       8  FceRI  1.793673e-42  2.869876e-41
    ## 24       8   CD16  1.000000e+00  1.000000e+00
    ## 25       9   CD25  0.000000e+00  0.000000e+00
    ## 26       9    CD4  0.000000e+00  0.000000e+00
    ## 27       9    CD3  0.000000e+00  0.000000e+00
    ## 28      10   CD56  0.000000e+00  0.000000e+00
    ## 29      10   Tbet  0.000000e+00  0.000000e+00
    ## 30      10   CD16  1.486227e-60  1.585309e-59

``` r
# Assignments
original_cluster_ident = as.character(rep(1:10))
new_cluster_ident = c("Monocytes", "cDC", "NKT Cells", "CD4+ T Cells", "B Cells","gd T Cells", "CD8+ T Cells", "Basophils", "CD4+ Treg", "NK")

# Set the new annotations
sce_gvhd = setClusterIdents(sce_gvhd, orig.ident = original_cluster_ident , new.ident = new_cluster_ident, clustering = 'flowsom_cc_k10')

# Vizualise in UMAP space
cell_annotation = metadataPlot(sce_gvhd,
    colby = 'cell_annotation',
    title = 'Annotated clusters',
    legendLabSize = 8,
    axisLabSize = 20,
    titleLabSize = 20,
    subtitleLabSize = 16,
    captionLabSize = 16)

annotated_clusters = unique(sce_gvhd@metadata$flowsom_cc_k10)

annotated_abundance = plotAbundance(sce_gvhd,
              clusters = annotated_clusters,
              clusterAssign = 'cell_annotation',
              feature = NULL,
              legendPosition = 'none', 
              stripLabSize = 10,
              axisLabSize = 14,
              titleLabSize = 12,
              subtitleLabSize = 10,
              captionLabSize = 18)

plot_grid(cell_annotation, annotated_abundance, 
    labels = c('A','B'),
    ncol = 2, nrow = 1, align = "l", label_size = 20)
```

![](README_files/figure-gfm/merge%20clusters-1.png)<!-- -->

``` r
bar_plot = plotAbundance(sce_gvhd,
              graph_type = 'bar',
              clusters = annotated_clusters,
              clusterAssign = 'cell_annotation',
              feature = 'condition',
              legendLabSize = 8,
              stripLabSize = 22,
              axisLabSize = 22,
              titleLabSize = 22,
              subtitleLabSize = 18,
              captionLabSize = 18)

bar_plot
```

![](README_files/figure-gfm/bar%20plot-1.png)<!-- -->

``` r
# With feature
feature = medianHeatmap(sce_gvhd, grouping = "cell_annotation", feature = NULL, scale_01 = T, heat_bar = "Greys")
```

![](README_files/figure-gfm/final%20heatmap-1.png)<!-- -->

``` r
cell_annotation_dif = metadataPlot(sce_gvhd,
    colby = 'condition',
    title = 'UMAP labelled by GvHD status',
    colkey = c(None = 'royalblue', GvHD = 'red2'),
    legendPosition = 'top',
    legendLabSize = 16,
    axisLabSize = 20,
    titleLabSize = 20,
    subtitleLabSize = 16,
    captionLabSize = 16)

cell_annotation_dif = cell_annotation_dif + facet_wrap(~condition)

abundance_dif = plotAbundance(sce_gvhd,
              clusters = annotated_clusters,
              clusterAssign = 'cell_annotation',
              feature = 'condition',
              colkey = c(None = 'royalblue', GvHD = 'red2'),
              legendPosition = 'none', 
              stripLabSize = 10,
              axisLabSize = 14,
              titleLabSize = 12,
              subtitleLabSize = 10,
              captionLabSize = 18)

plot_grid(cell_annotation_dif, abundance_dif,
    labels = c('A','B'),
    ncol = 2, align = "l", label_size = 20)
```

![](README_files/figure-gfm/cell%20annotations-1.png)<!-- -->

``` r
 library(broom)
# This need the individual values
# This also need the right contrasts
porportions = diffClust(sce_gvhd, 
                        clustering = 'cell_annotation')

porportions[,c(1,14:17)]
```

    ##         cluster proportion proportion.1       logfc       pval
    ## 1       B Cells  0.7571028    6.4246483 -2.13839808 0.01040562
    ## 2     Basophils  0.4749399    0.5833272 -0.20556014 0.52183939
    ## 3  CD4+ T Cells  8.4219643   12.6510706 -0.40689875 0.26233168
    ## 4     CD4+ Treg  1.1673934    0.7770172  0.40706623 0.87278012
    ## 5  CD8+ T Cells 18.5345116   18.8898577 -0.01899066 0.87278012
    ## 6           cDC  0.8772951    2.0589283 -0.85309741 0.07816909
    ## 7    gd T Cells  1.7899011    2.4753240 -0.32421095 0.74877404
    ## 8     Monocytes 55.0714116   42.1360672  0.26772666 0.42333964
    ## 9            NK 12.2719019   13.5053902 -0.09577663 0.52183939
    ## 10    NKT Cells  0.6335783    0.4983693  0.24004216 0.74877404

``` r
# Need to catch errors - cluster not there causes crash
marker_test = diffExpression(sce_gvhd,
                assay = 'scaled',
                grouping = 'group', # The condensing function
                feature = 'condition', # The contrast
                clusterAssign = 'cell_annotation')

sig_markers = marker_test[which(marker_test$p.value < 0.05),]

sig_markers[c(1,12:14),c(13:19)]
```

    ##                        mean_1     mean_2     logfc group1 group2     p.value
    ## cDC_PDL_1           0.2864016 0.01616282  2.874682   None   GvHD 0.015567638
    ## CD8+ T Cells_CD45RA 0.4558150 1.56237437 -1.231875   None   GvHD 0.003947752
    ## CD8+ T Cells_CD27   2.6619634 0.55987973  1.559097   None   GvHD 0.037372988
    ## CD8+ T Cells_PD_1   1.5150474 0.31405018  1.573649   None   GvHD 0.044951243
    ##                          padj
    ## cDC_PDL_1           0.5929411
    ## CD8+ T Cells_CD45RA 0.5929411
    ## CD8+ T Cells_CD27   0.5929411
    ## CD8+ T Cells_PD_1   0.5929411

``` r
# Difference in activation markers per cluster
markerExpressionPerSample(sce_gvhd,
    caption = '',
    clusters = c("CD8+ T Cells", "cDC"),
    feature = 'condition',
    clusterAssign = 'cell_annotation',
    markers = activation_markers,
    colkey = c('royalblue', 'red2'),
    title = '',
    stripLabSize = 12,
    axisLabSize = 10,
    titleLabSize = 18,
    subtitleLabSize = 10,
    captionLabSize = 10)
```

![](README_files/figure-gfm/activation%20markers-1.png)<!-- -->

# Session info

``` r
sessionInfo()
```

    ## R version 3.6.2 (2019-12-12)
    ## Platform: x86_64-apple-darwin15.6.0 (64-bit)
    ## Running under: macOS High Sierra 10.13.6
    ## 
    ## Matrix products: default
    ## BLAS:   /Library/Frameworks/R.framework/Versions/3.6/Resources/lib/libRblas.0.dylib
    ## LAPACK: /Library/Frameworks/R.framework/Versions/3.6/Resources/lib/libRlapack.dylib
    ## 
    ## locale:
    ## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
    ## 
    ## attached base packages:
    ## [1] parallel  stats4    stats     graphics  grDevices utils     datasets 
    ## [8] methods   base     
    ## 
    ## other attached packages:
    ##  [1] broom_0.5.5                 tibble_2.1.3               
    ##  [3] scran_1.14.6                umap_0.2.5.0               
    ##  [5] pheatmap_1.0.12             Rtsne_0.15                 
    ##  [7] reshape2_1.4.3              limma_3.42.2               
    ##  [9] FlowSOM_1.18.0              igraph_1.2.5               
    ## [11] ConsensusClusterPlus_1.50.0 RColorBrewer_1.1-2         
    ## [13] cowplot_1.0.0               stringr_1.4.0              
    ## [15] flowCore_1.52.1             immunoCluster_0.1.0        
    ## [17] dplyr_0.8.5                 ggrepel_0.8.2              
    ## [19] ggplot2_3.3.0               SingleCellExperiment_1.8.0 
    ## [21] SummarizedExperiment_1.16.1 DelayedArray_0.12.2        
    ## [23] BiocParallel_1.20.1         matrixStats_0.56.0         
    ## [25] Biobase_2.46.0              GenomicRanges_1.38.0       
    ## [27] GenomeInfoDb_1.22.1         IRanges_2.20.2             
    ## [29] S4Vectors_0.24.3            BiocGenerics_0.32.0        
    ## [31] kableExtra_1.1.0            knitr_1.28                 
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] backports_1.1.5          plyr_1.8.6               splines_3.6.2           
    ##   [4] fda_2.4.8.1              scater_1.14.6            digest_0.6.25           
    ##   [7] htmltools_0.4.0          viridis_0.5.1            fansi_0.4.1             
    ##  [10] magrittr_1.5             CytoML_1.12.1            cluster_2.1.0           
    ##  [13] ks_1.11.7                readr_1.3.1              RcppParallel_5.0.0      
    ##  [16] R.utils_2.9.2            askpass_1.1              flowWorkspace_3.34.1    
    ##  [19] jpeg_0.1-8.1             colorspace_1.4-1         rvest_0.3.5             
    ##  [22] rrcov_1.5-2              xfun_0.12                crayon_1.3.4            
    ##  [25] RCurl_1.98-1.1           jsonlite_1.6.1           hexbin_1.28.1           
    ##  [28] graph_1.64.0             glue_1.3.2               flowClust_3.24.0        
    ##  [31] gtable_0.3.0             zlibbioc_1.32.0          XVector_0.26.0          
    ##  [34] webshot_0.5.2            ggcyto_1.14.1            BiocSingular_1.2.2      
    ##  [37] IDPmisc_1.1.20           Rgraphviz_2.30.0         DEoptimR_1.0-8          
    ##  [40] scales_1.1.0             mvtnorm_1.1-0            edgeR_3.28.1            
    ##  [43] Rcpp_1.0.4               viridisLite_0.3.0        clue_0.3-57             
    ##  [46] reticulate_1.13          dqrng_0.2.1              openCyto_1.24.0         
    ##  [49] rsvd_1.0.3               mclust_5.4.5             tsne_0.1-3              
    ##  [52] httr_1.4.1               ellipsis_0.3.0           pkgconfig_2.0.3         
    ##  [55] XML_3.99-0.3             R.methodsS3_1.8.0        farver_2.0.3            
    ##  [58] flowViz_1.50.0           locfit_1.5-9.1           utf8_1.1.4              
    ##  [61] flowStats_3.44.0         tidyselect_1.0.0         labeling_0.3            
    ##  [64] rlang_0.4.5              munsell_0.5.0            tools_3.6.2             
    ##  [67] cli_2.0.2                generics_0.0.2           evaluate_0.14           
    ##  [70] yaml_2.2.1               robustbase_0.93-6        purrr_0.3.3             
    ##  [73] nlme_3.1-145             RBGL_1.62.1              R.oo_1.23.0             
    ##  [76] xml2_1.2.2               compiler_3.6.2           rstudioapi_0.11         
    ##  [79] beeswarm_0.2.3           png_0.1-7                statmod_1.4.34          
    ##  [82] pcaPP_1.9-73             stringi_1.4.6            RSpectra_0.16-0         
    ##  [85] lattice_0.20-40          Matrix_1.2-18            vctrs_0.2.4             
    ##  [88] pillar_1.4.3             lifecycle_0.2.0          BiocNeighbors_1.4.2     
    ##  [91] data.table_1.12.8        bitops_1.0-6             irlba_2.3.3             
    ##  [94] corpcor_1.6.9            R6_2.4.1                 latticeExtra_0.6-29     
    ##  [97] KernSmooth_2.23-16       gridExtra_2.3            vipor_0.4.5             
    ## [100] MASS_7.3-51.5            gtools_3.8.1             assertthat_0.2.1        
    ## [103] openssl_1.4.1            withr_2.1.2              mnormt_1.5-6            
    ## [106] GenomeInfoDbData_1.2.2   hms_0.5.3                ncdfFlow_2.32.0         
    ## [109] grid_3.6.2               tidyr_1.0.2              rmarkdown_2.1           
    ## [112] DelayedMatrixStats_1.8.0 base64enc_0.1-3          ggbeeswarm_0.6.0        
    ## [115] ellipse_0.4.1

# Contact

For any queries relating to software:

  - James W Opzoomer (<james.opzoomer@kcl.ac.uk>)
