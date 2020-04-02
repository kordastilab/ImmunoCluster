immunoCluster
================
James Opzoomer, Kevin Blighe, Jessica Timms
2020-04-02

**NOTE: THIS PACKAGE IS STILL UNDER DEVELOPMENT AND SO SOME OF THE
FUNCTIONALITY IS NOT FULLY TESTED**

## 1\. Introduction to immunoCluster

The immunoCluster package provides a complete toolkit to carry out
immune profiling from both liquid and imaging high-dimensional mass
(CyTOF) and flow cytometry. The package features standardized data
infrastructure, making use of the SingleCellExperiment class object,
interactive visualization tools and convenient implementations for
popular dimenionality reduction and clustering algorithms designed for a
non-specialist. To learn using immunoCluster we have provided
walkthrough below exploring some of the packages basic functionality.

## 2\. Main functionalities of immunoCluster

The basic pipeline of immunoCluster analysis involves:

  - Creating file and marker labels to generate the experimental design

  - Importing the .fcs data into the SingelCellExperiment object

  - Exploratory analysis of the data

  - Dimensionality reduction

  - Clustering and biological interpretation of the clusters

  - Differential abundace and differential expression analysis

# Installation

## 1\. Install the package from github

Install the development version of the package, installation will
usually take a few minutes, depending on the number of dependencies that
have been installed on your pc. Install immunoCluster from github here:

``` r
  # Install devtools for github installation if not present
  require(devtools)
  
  # Install the ImmunoCluster package
  devtools::install_github("kordastilab/ImmunoCluster")
  
  # May need to run this line if an R build version error while installing appears
  # Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS="true")
```

## 2\. Load the package and the dependancies into the R session

Load necessary packages:

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

# Walkthrough: Immune profiling PBMC from a two condition experimental design

## 1\. Introduction

In this walkthrough we will learn to perform immune profiling and
examine the from changes in immune populations from human PBMC. We will
start with pre-gated .fcs data files. The samples were taken from
patients following hematopoetic reconsitution in leukemia patients
following BMT. The dataset represents 15 individuals, sampled at
multiple time points after BMT, for a total of 28 samples. Of these
patients, a small subset suffered from graft versus host disease (GvHD,
n = 3), whereas most other patients did not experience such
complications.

We have performed gating cleanup on a subset of the dataset which can be
dowloaded from zenodo, along with all the other data required for the
walkthough. \[Insert zenodo button\].

The data was originally reported in [Comprehensive Immune Monitoring of
Clinical Trials to Advance Human
Immunotherapy](https://www.sciencedirect.com/science/article/pii/S2211124719308228).
The full raw data (not used here, but which was used in the
immunoCluster paper) is available
[here](http://flowrepository.org/id/FR-FCM-Z244).

## 2\. Download and import data

You can download all the .fcs files from zenodo here here, or
alternatively you can download the SingleCellExperiment object directly
into memory and skip the import process just below.

``` r
# Insert zenodo link
```

The experimental design informs the what metadata we will need to
incorporate into the SingleCellExperiment object, here we have created
one table (from an excel file), that contains five columns, and one row
for each .fcs file. We will include the .fcs filename, the GvHD
“condition” of the patient, patient\_id, a unique sample\_id, and the
day\_id representing the day following BMT that the sample was taken.
(Link to excel file).

``` r
# Excel file can be read in using this line
#

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

We also include a table that represents a row for every channel in the
original fcs file. The columns represent the original name if the
channel, the name we want to rename it in antigen. We then optionally
include columns of binary 1 and 0s that allow us to group markers into
relevant sets, here for clustering or differential expression analysis.
Any marker that is not included in either group will be discarded. This
table must be an exact match to the order that the markers are
represnted in the raw .fcs files. (Link to excel file).

``` r
# Excel file can be read in using this line
#

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

We next load the .fcs data and the metadata into the
SingleCellExperiment object. Here the number of rows in the metadata
must match the number of .fcs files to be analysed. We have the option
to specify the scaling transformation and the cofactor. We have used
asinh with cofactor 5. It is possible to downsaple from each file or by
a specified metadata grouping, although here we have not downsampled.
The columns are also renamed as specified in the panel\_metadata table
and columns not needed are discarded.

``` r
# Collect all .fcs files from a directory
filelist = list.files(
      path = "../github_readme/fcs_data",
    pattern = "*.fcs|*.FCS",
    full.names = TRUE)

# Select which metadata to include
metadata = data.frame(
      file = filelist, 
      group = sample_metadata_df$sample_id,
      condition = sample_metadata_df$condition,
      patient_id = sample_metadata_df$patient_id,
      day_id = sample_metadata_df$day_id,
      row.names = filelist,
      stringsAsFactors = FALSE)

# Create SCE object with CyTOF data inside
require(SingleCellExperiment)
  sce_gvhd = processFCS(
    files = filelist,
    metadata = metadata,
    transformation = TRUE, # Transform data
    filter = FALSE, 
    transFun = function (x) asinh(x), # Using asinh function (cofactor 5 is used by default)
    downsampleVar = NULL, 
    downsample = NULL, # Downsample to n cells per downsample_grouping
    downsample_grouping = NULL, # Downsample to downsample cells by the specified metadata slot
    newColnames = col_names, # Rename columns from panel metadata 
    colsDiscard = colsDiscard) # Discard columns not selected for downstream analysis

sce_gvhd
```

Using the code chunk below you can download the imported r data object
with the metatdat incorporated directly into memory to skip the import
steps and save time:

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

We provide several functions to examine globla patterns within the data,
they can be used to get an idea of the influence of conditions/treatment
or batch effects within the data. The mdsplot() function will generate a
multi-dimensional scaling (MDS) plot on median expresion values for each
channel. MDS plots represent an unsupervised way to indicate
similarities between samples at a global level berfore more in depth
analsysi. In the example here we can see that there is no observabe
high-level difference between the two GvHD condition samples, however
there does appear to be a weak separation between samples collected at
D30 and D90 post BMT. Here and throughout the pipeline, colkey can be
used to specify plotting colors for graphing conditions.

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

Further to the mds plots,the medianHeatmap() function will generate a
heatmap with median marker expression across all (default) or a
selection (with markers = ). The heatmap is clustered over rows and
columns. The heatmap can be generated on the asinh tranformed values
(default) or 0-\>1 scaled (with scale\_01 = T). The heatmap should
provide insight as to which markers are highly expressed and wether
thier marker expression profile is assosciated with a metadata
condition. Highly expressed markers will have a greater influence on
clustering and this will help understand clustering results generated
downstream.

``` r
feature = medianHeatmap(sce_gvhd, grouping = "group", feature = "condition", feature_cols = c('royalblue', 'red2'))
```

![](README_files/figure-gfm/sample%20heatmap-1.png)<!-- -->

## 3\. Dimensionality reduction and clustering

immunoCluster supports multiple different dimensionality reduction
algorithms and has UMAP and tSNE functionality implemented below with
the ability to downsample before running the algorithms. We run both
UMAP and tSNE on selected lineage markers, downsampling to 1000 cells
per sample to reduce runtime. The dimensionality reductions are stored
in the reducedDimanmes slot and any other dimensionality reduction, like
PCA, can also be stored in the SCE object in parallel.

``` r
require(umap)
# Run UMAP and store in sce object
sce_gvhd = performUMAP(sce_gvhd, downsample = 1000, useMarkers = clustering_markers)

# Run tSNE and store in sce object
sce_gvhd = performTSNE(sce_gvhd, downsample = 1000, useMarkers = clustering_markers)

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
    ## reducedDimNames(2): UMAP Rtsne
    ## spikeNames(0):
    ## altExpNames(0):

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

``` r
library(diffcyt)
```

``` r
# SubsetSCE function
```

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
    ##  [1] diffcyt_1.6.4               broom_0.5.5                
    ##  [3] tibble_3.0.0                scran_1.14.6               
    ##  [5] umap_0.2.5.0                pheatmap_1.0.12            
    ##  [7] Rtsne_0.15                  reshape2_1.4.3             
    ##  [9] limma_3.42.2                FlowSOM_1.18.0             
    ## [11] igraph_1.2.5                ConsensusClusterPlus_1.50.0
    ## [13] RColorBrewer_1.1-2          cowplot_1.0.0              
    ## [15] stringr_1.4.0               flowCore_1.52.1            
    ## [17] immunoCluster_0.1.0         dplyr_0.8.5                
    ## [19] ggrepel_0.8.2               ggplot2_3.3.0              
    ## [21] SingleCellExperiment_1.8.0  SummarizedExperiment_1.16.1
    ## [23] DelayedArray_0.12.2         BiocParallel_1.20.1        
    ## [25] matrixStats_0.56.0          Biobase_2.46.0             
    ## [27] GenomicRanges_1.38.0        GenomeInfoDb_1.22.1        
    ## [29] IRanges_2.20.2              S4Vectors_0.24.3           
    ## [31] BiocGenerics_0.32.0         kableExtra_1.1.0           
    ## [33] knitr_1.28                 
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] utf8_1.1.4               reticulate_1.13          R.utils_2.9.2           
    ##   [4] ks_1.11.7                tidyselect_1.0.0         lme4_1.1-21             
    ##   [7] grid_3.6.2               munsell_0.5.0            codetools_0.2-16        
    ##  [10] statmod_1.4.34           withr_2.1.2              colorspace_1.4-1        
    ##  [13] flowViz_1.50.0           rstudioapi_0.11          flowClust_3.24.0        
    ##  [16] robustbase_0.93-6        openCyto_1.24.0          labeling_0.3            
    ##  [19] GenomeInfoDbData_1.2.2   mnormt_1.5-6             farver_2.0.3            
    ##  [22] flowWorkspace_3.34.1     TH.data_1.0-10           vctrs_0.2.4             
    ##  [25] generics_0.0.2           xfun_0.12                R6_2.4.1                
    ##  [28] ggbeeswarm_0.6.0         clue_0.3-57              rsvd_1.0.3              
    ##  [31] locfit_1.5-9.4           bitops_1.0-6             assertthat_0.2.1        
    ##  [34] scales_1.1.0             multcomp_1.4-12          beeswarm_0.2.3          
    ##  [37] gtable_0.3.0             sandwich_2.5-1           rlang_0.4.5             
    ##  [40] GlobalOptions_0.1.1      splines_3.6.2            hexbin_1.28.1           
    ##  [43] yaml_2.2.1               backports_1.1.5          IDPmisc_1.1.20          
    ##  [46] RBGL_1.62.1              tools_3.6.2              ellipsis_0.3.0          
    ##  [49] Rcpp_1.0.4               plyr_1.8.6               base64enc_0.1-3         
    ##  [52] zlibbioc_1.32.0          purrr_0.3.3              RCurl_1.98-1.1          
    ##  [55] openssl_1.4.1            GetoptLong_0.1.8         viridis_0.5.1           
    ##  [58] zoo_1.8-7                cluster_2.1.0            fda_2.4.8.1             
    ##  [61] magrittr_1.5             ncdfFlow_2.32.0          data.table_1.12.8       
    ##  [64] RSpectra_0.16-0          circlize_0.4.8           mvtnorm_1.1-0           
    ##  [67] hms_0.5.3                evaluate_0.14            XML_3.99-0.3            
    ##  [70] jpeg_0.1-8.1             mclust_5.4.5             gridExtra_2.3           
    ##  [73] shape_1.4.4              ggcyto_1.14.1            compiler_3.6.2          
    ##  [76] scater_1.14.6            ellipse_0.4.1            flowStats_3.44.0        
    ##  [79] KernSmooth_2.23-16       crayon_1.3.4             minqa_1.2.4             
    ##  [82] R.oo_1.23.0              htmltools_0.4.0          corpcor_1.6.9           
    ##  [85] pcaPP_1.9-73             tidyr_1.0.2              rrcov_1.5-2             
    ##  [88] RcppParallel_5.0.0       ComplexHeatmap_2.2.0     MASS_7.3-51.5           
    ##  [91] boot_1.3-24              Matrix_1.2-18            readr_1.3.1             
    ##  [94] cli_2.0.2                R.methodsS3_1.8.0        pkgconfig_2.0.3         
    ##  [97] xml2_1.2.5               vipor_0.4.5              dqrng_0.2.1             
    ## [100] webshot_0.5.2            XVector_0.26.0           rvest_0.3.5             
    ## [103] digest_0.6.25            tsne_0.1-3               graph_1.64.0            
    ## [106] rmarkdown_2.1            edgeR_3.28.1             DelayedMatrixStats_1.8.0
    ## [109] gtools_3.8.2             rjson_0.2.20             nloptr_1.2.2.1          
    ## [112] lifecycle_0.2.0          nlme_3.1-145             jsonlite_1.6.1          
    ## [115] BiocNeighbors_1.4.2      viridisLite_0.3.0        askpass_1.1             
    ## [118] fansi_0.4.1              pillar_1.4.3             lattice_0.20-40         
    ## [121] survival_3.1-11          httr_1.4.1               DEoptimR_1.0-8          
    ## [124] glue_1.3.2               png_0.1-7                Rgraphviz_2.30.0        
    ## [127] stringi_1.4.6            BiocSingular_1.2.2       CytoML_1.12.1           
    ## [130] latticeExtra_0.6-29      irlba_2.3.3

# Contact

For any queries relating to software:

  - James W Opzoomer (<james.opzoomer@kcl.ac.uk>)
