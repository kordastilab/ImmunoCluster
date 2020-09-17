immunoCluster
================
James Opzoomer, Kevin Blighe, Jessica Timms
2020-09-17

## 1\. Introduction to immunoCluster

The immunoCluster package usese the [scDataViz bioconductor
package’s](https://bioconductor.org/packages/devel/bioc/vignettes/scDataviz/inst/doc/scDataviz.html)
vizualization tools and adaptation of the
[SingleCellExperiment](https://www.bioconductor.org/packages/release/bioc/vignettes/SingleCellExperiment/inst/doc/intro.html)
data structure to build simple and flexible cytometry analysis workflows
like those outlined in Nowicka et al. (2017) [CyTOF workflow:
differential discovery in high-throughput high-dimensional cytometry
datasets](https://f1000research.com/articles/6-748).

The package provides a broad toolkit to carry out immune profiling from
both liquid and imaging high-dimensional mass (CyTOF) and flow cytometry
data. immunoCluster features standardised data infrastructure, making
use of the SingleCellExperiment class object, scDataViz interactive
visualization tools along with convenient implementations for popular
dimensionality reduction and clustering algorithms designed for a
non-specialist. To learn using immunoCluster we have provided a
walkthrough below exploring some of the package’s basic functionality on
a published dataset.

An in depth report of immunoCluster’s use to create liquid and imaging
mass cytometry analysis pipelines, as well as its use for fluorescent
flow cytometry analysis is outlined in our recent preprint
[here](https://www.biorxiv.org/content/10.1101/2020.09.09.289033v1).

## 2\. Main functionalities of immunoCluster

The basic pipeline of immunoCluster analysis involves:

  - Creating file and marker labels to generate the experimental design

  - Importing .fcs data into the SingleCellExperiment object

  - Exploratory analysis of the data

  - Dimensionality reduction

  - Clustering and biological interpretation of the clusters

  - Differential abundance and differential expression analysis

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

## 2\. Load required packages

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

In this walkthrough we will learn to perform immune profiling with
immunoCluster to examine changes in immune populations from human PBMC
following bone marrow transplantation (BMT).

We will start with pre-gated .fcs data files. The samples were taken
from patients following hematopoetic reconstitution in leukemia patients
following BMT. The dataset represents 15 individuals, sampled at two
time points (Day 30 and 90) after BMT, for a total of 28 samples. Of
these patients, a small subset suffered from graft versus host disease
(GvHD, n = 3), whereas most other patients did not experience such
complications.

We have performed gating cleanup on a subset of the dataset which can be
dowloaded from zenodo, along with all the other data required for the
walkthough.

The data was originally reported in [Comprehensive Immune Monitoring of
Clinical Trials to Advance Human
Immunotherapy](https://www.sciencedirect.com/science/article/pii/S2211124719308228).
The full raw data (not used here, but which was used in the
immunoCluster paper) is available
[here](http://flowrepository.org/id/FR-FCM-Z244).

All the data required can be downloaded from Zenodo.
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3801882.svg)](https://doi.org/10.5281/zenodo.3801882)

## 2\. Download and import data

You can download all the .fcs files from zenodo
[here](https://doi.org/10.5281/zenodo.3801882), or alternatively you can
download the SingleCellExperiment object directly into memory and skip
the import process just below. The fcs files can be downloaded from this
[download\_link](https://zenodo.org/record/3801882/files/fcs_data.zip?download=1).

The experimental design informs the what metadata we will need to
incorporate into the SingleCellExperiment (SCE) object, here we have
created one table (from an excel file), that contains five columns, and
one row for each .fcs file. We will include the .fcs filename, the GvHD
“condition” of the patient, patient\_id, a unique sample\_id, and the
day\_id representing the day following BMT that the sample was taken.
[Download\_link](https://zenodo.org/record/3801882/files/sample_metadata.xlsx?download=1)

``` r
# Download sample metadata from Zenodo
sample_metadata_df = readRDS(url("https://zenodo.org/record/3801882/files/sample_metadata.rds"))

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

We also create a table that represents the marker panel information.
This table contains a row for every channel in the original fcs file.
The columns represent the original name of the channel and the channel
name we want to rename it to in antigen. We then optionally include
columns of binary 1 and 0s that allow us to group markers into relevant
sets, here we have two sets, one group of markers for clustering and
another for downstream expression analysis on our cell populations of
interest. Any marker that is not included in either group will be
discarded when we build the SCE object. This table must be an exact
match to the order that the markers are represnted in the raw .fcs files
or the channel renaming will be incorrect or fail.
[Download\_link](https://zenodo.org/record/3801882/files/panel_metadata.xlsx?download=1).

``` r
# Download panel metadata from Zenodo
panel_metadata_df = readRDS(url("https://zenodo.org/record/3801882/files/panel_metadata.rds"))

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

We next load the .fcs data and the metadata into the SCE object. Here
the number of rows in the metadata must match the number of .fcs files
to be analysed. We have the option to specify the scaling transformation
and the cofactor. We have used asinh with cofactor 5. It is possible to
downsaple from each file or by a specified metadata grouping, although
here we have not downsampled in this case. The columns are also renamed
as specified in the panel\_metadata table and columns not needed are
discarded.

``` r
# Collect all .fcs files from a directory
filelist = list.files(
      path = "fcs_data/",
    pattern = "*.fcs|*.FCS",
    full.names = TRUE)

# Select which metadata to include from the sample metadata table
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
    filter = FALSE, # Do not perform noise filtering
    transformation = TRUE, # Transform data
    transFun = function (x) asinh(x), # Using asinh function (cofactor 5 is used by default)
    downsampleVar = NULL, # No downsampling by variance
    downsample = NULL, # Downsample to n cells per downsample_grouping 
    downsample_grouping = NULL, # if downsample = n, downsample n cells by the specified metadata grouping.
    newColnames = col_names, # Rename columns from panel metadata 
    colsDiscard = colsDiscard) # Discard columns not selected for downstream analysis

sce_gvhd
```

Using the code chunk below you can download the imported r data object
with the metatdata incorporated, directly into memory to skip the import
steps to save time:

``` r
# Download complete ImmunoCluster SCE object from zenodo 
sce_gvhd = readRDS(url("https://zenodo.org/record/3801882/files/sce_gvhd.rds"))

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

We provide several functions to examine global patterns within the data,
they can be used to get an idea of the influence of conditions/treatment
or batch effects within the data. The mdsplot() function will generate a
multi-dimensional scaling (MDS) plot on median expression values for
each channel. MDS plots represent an unsupervised way to indicate
similarities between samples at a global level before more in depth
analsysis. In the example here we can see that there is no observable
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

<img src="README_files/figure-gfm/mds_plot-1.png" style="display: block; margin: auto;" />

Further to the mds plots,the medianHeatmap() function will generate a
heatmap with median marker expression across all (default) or a
selection (with markers = ). The heatmap is clustered over rows and
columns. The heatmap can be generated on the asinh tranformed values
(default) or 0-min to 1-max scaled (with scale\_01 = T). The heatmap
should provide insight as to which markers are highly expressed and
whether their marker expression profile is associated with a metadata
condition. Highly expressed markers will have a greater influence on
clustering and this will help understand clustering results generated
downstream.

``` r
feature = medianHeatmap(sce_gvhd, grouping = "group", feature = "condition", feature_cols = c('royalblue', 'red2'))
```

<img src="README_files/figure-gfm/sample heatmap-1.png" style="display: block; margin: auto;" />

## 3\. Dimensionality reduction and clustering

immunoCluster supports multiple different dimensionality reduction
algorithms and has UMAP and tSNE functionality implemented below with
the ability to downsample before running the algorithms. We run both
UMAP and tSNE on selected lineage markers, downsampling to 1000 randomly
selected cells per sample to reduce runtime. The dimensionality
reductions are stored in the reducedDimNames slot and any other
dimensionality reduction, like PCA, can also be stored in the SCE object
in parallel.

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

The markerExpression() function will allow us to view the UMAP
dimensionality reductions, overlaying marker expression of our major
immune supopulation lineage markers as an exporatory analysis of the
resultant UMAP. The dimensionality reduction data to use is specified
with the reducedDim parameter.The markers and several ggplot style
plotting parameters can be specified resulting in tiled plots of the
UMAP with marker expression as colour.

``` r
expression_markers =  c('CD3', 'CD4', 'CD8a', 'CD11b', 'CD19', 'CD56')

exp_plot_umap = markerExpression(sce_gvhd,
    markers = expression_markers,
    reducedDim = 'UMAP',
    title = 'UMAP',
    nrow = 1, ncol = 6,
    pointSize = 0.05,
    legendKeyHeight = 1.0,
    legendLabSize = 14,
    stripLabSize = 20,
    axisLabSize = 18,
    titleLabSize = 20,
    captionLabSize = 22)

exp_plot_tsne = markerExpression(sce_gvhd,
    markers = expression_markers,
    reducedDim = 'Rtsne',
    title = 'tSNE',
    nrow = 1, ncol = 6,
    pointSize = 0.05,
    legendKeyHeight = 1.0,
    legendLabSize = 14,
    stripLabSize = 20,
    axisLabSize = 18,
    titleLabSize = 20,
    captionLabSize = 22)

plot_grid(exp_plot_umap, exp_plot_tsne,
    labels = c('A','B'),
    nrow = 2, align = "l", label_size = 24)
```

<img src="README_files/figure-gfm/unnamed-chunk-2-1.png" style="display: block; margin: auto;" />

immunoCluster provides several wrapper functions to perform unsupervised
clustering on your data. Currently users can implement either
Rphenograph or FlowSOM clustering. Here we have used an ensemble method
of FlowSOM and ConsensusCLusterPlus, which allows us to perform a fast,
high resolution clustering on our dataset. We perform the clustering on
all cells using a selection of lineage markers, with a final desired
cluster number of up to 10, specified by the parameter k. FlowSOM will
cluster cells into 100 SOM codes defined by the (dimensions of som\_x
and som\_y) and these will be clustered into 2-k clusters and all of
these clustering results are saved as metadata in the SCE
object.

``` r
sce_gvhd = runFlowSOM(sce_gvhd, k = 10, markers = clustering_markers, som_x = 10, som_y = 10)
```

    ## Building SOM

    ## Mapping data to SOM

The clustering identities are stored in the metadata dataframe. The per
cell membership of the 100 SOM codes is stored under som\_codes and the
respective clusterings are stored under flowsom\_cc\_k(k=n). We can
vizualise the clustering or any other metatdata on the dimnesionality
plot using the function metadataPlot() (which maps to UMAP by default
but can specify of dimension reductions that are stored like tSNE). We
can also vizualise the abundance of each cluster as a proportion of
total cells per sample using the plotAbundance() function.

``` r
# Plot UMAP with flowsom k=10 clustering overlay
fsom_k10 = metadataPlot(sce_gvhd,
    colby = 'flowsom_cc_k10',
    title = '',
    reducedDim = 'UMAP',
    legendPosition = 'right',
    legendLabSize = 8,
    axisLabSize = 10,
    titleLabSize = 10,
    subtitleLabSize = 1,
    captionLabSize = 16)

# Specify all clusters to plot abundance for
fsom_clusters = unique(sce_gvhd@metadata$flowsom_cc_k10)

# Create boxplot of flowsom k=10 abundance 
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
    ncol = 2, align = "l", label_size = 20, rel_widths = c(1.4, 1))
```

<img src="README_files/figure-gfm/viz flowSOM-1.png" style="display: block; margin: auto;" />

## 4\. Biological annotation of the clustering

The medianHeatmap() function can also be used to display the median
marker expression values of each flowSOM metacluster, by specifying the
metadata clustering in the grouping parameter. Here the median values
are 0-1 scaled to emphasize the differences between clusters (using
scale\_01 = T).

``` r
# With feature
feature = medianHeatmap(sce_gvhd, grouping = "flowsom_cc_k10", feature = NULL, scale_01 = T)
```

<img src="README_files/figure-gfm/diagnostic heatmap-1.png" style="display: block; margin: auto;" />

The heatmap can inform the supposed biological identity of the clusters
generated by flowSOM and consesusClusterPlus. Another way to understand
which markers are enriched per cluster is by using the finMarkers()
function from the scran package. Here we perform a wilcox test to
understand which markers have a LogFC \> 1.5 enrichment between one
cluster and several others.

``` r
library(scran)
library(tibble)

# Generate list of upregulated markers
cytof_markers = findMarkers(assay(sce_gvhd), sce_gvhd@metadata$flowsom_cc_k10, test.type="wilcox", direction="up",lfc=1.5, pval.type="some")

cluster_markers = NULL
top_n = 3 # The number of top upregulated markers to display
for (i in 1:length(cytof_markers)) {
  
  deg_markers = cytof_markers[[i]][1:top_n,1:2]
  
  deg_markers = rownames_to_column(as.data.frame(deg_markers), var = "marker")
  
  markers = data.frame(cluster = rep(i, top_n), deg_markers)
  
  cluster_markers = rbind(cluster_markers, markers)
  
}

# Show top 3 DEG markers per cluster
head(cluster_markers)
```

    ##   cluster marker p.value FDR
    ## 1       1   CD14       0   0
    ## 2       1   CD33       0   0
    ## 3       1  CD11b       0   0
    ## 4       2  CD11c       0   0
    ## 5       2   CD16       0   0
    ## 6       2  CD123       0   0

From this information we can assign identities to the clusters using the
setClusterIdents() function. We feed a vector of original clusters and
the new clusters we wish to assign these original cluster assignments
to. The new cluster assignments will be stored in the metadata as
cell\_annotation. We can then assign this cell\_annotation to a new
metadata slot if we wish to store and reassign using setClusterIdents()
again.

``` r
# Assignments
original_cluster_ident = as.character(rep(1:10))
new_cluster_ident = c("Monocytes", "cDC", "NKT Cells", "CD4+ T Cells", "B Cells","gd T Cells", "CD8+ T Cells", "Basophils", "CD4+ Treg", "NK")

# Set the new annotations
sce_gvhd = setClusterIdents(sce_gvhd, orig.ident = original_cluster_ident , new.ident = new_cluster_ident, clustering = 'flowsom_cc_k10')

# Define population colours
cols = brewer.pal(10, "Set3")

# Assign to populations for plotting
col_key = c("Monocytes" = cols[1], "cDC" = cols[2], "NKT Cells" = cols[3], "CD4+ T Cells" = cols[4], "B Cells" = cols[5], "gd T Cells" = cols[6], "CD8+ T Cells" = cols[7], "Basophils" = cols[8], "CD4+ Treg" = cols[9], "NK" = cols[10])

# Vizualise in UMAP space
cell_annotation = metadataPlot(sce_gvhd,
    colby = 'cell_annotation',
    colkey = col_key, # Add colours for plotting 
    title = '',
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
              colkey = col_key,
              legendPosition = 'none', 
              stripLabSize = 10,
              axisLabSize = 14,
              titleLabSize = 12,
              subtitleLabSize = 10,
              captionLabSize = 18)

plot_grid(cell_annotation, annotated_abundance, 
    labels = c('A','B'), align = "l", label_size = 20, rel_widths = c(1.5,1), rel_heights = c(1.3,1))
```

<img src="README_files/figure-gfm/merge clusters-1.png" style="display: block; margin: auto;" />

The plotAbundance() function additionally allows to creat bar plots but
specifying graph\_type = “bar”. Specfying a metadata feature to split
the data by will allow to another level of comparison if wanted, first
by plotting abundance by sample as a stacked bar and then arranging the
samples by GvHD condition.

``` r
bar_plot = plotAbundance(sce_gvhd,
              graph_type = 'bar',
              clusters = annotated_clusters,
              clusterAssign = 'cell_annotation',
              feature = 'condition',
              colkey = col_key,
              legendLabSize = 7,
              stripLabSize = 22,
              axisLabSize = 22,
              titleLabSize = 22,
              subtitleLabSize = 18,
              captionLabSize = 18)

bar_plot
```

<img src="README_files/figure-gfm/bar plot-1.png" style="display: block; margin: auto;" />

Once we have defined the biological identity of our clusters we can
create a final heatmap to demonstrate that the marker expression
patterns are consistent with our labels. Here we use medianHeatmap to
create 0-1 normalised heatmap per cluster. We can also modify the the
heat bar from the default greens to a black-white gradient (using
heat\_bar = “greys”).

``` r
# With feature
feature = medianHeatmap(sce_gvhd, 
                        grouping = "cell_annotation",
                        markers = clustering_markers,
                        feature = NULL, 
                        scale_01 = T, 
                        heat_bar = "Greys")
```

<img src="README_files/figure-gfm/final heatmap-1.png" style="display: block; margin: auto;" />

## 5\. Differential abundance and expression analysis

We can now visualise the data by disease status to better understand the
differences that may help us define biomarkers to predict GvHD status
post BMT. We can also use the metadataPlot() function to split the UMAP
by GvHD status by defining colkey and splitting the UMAP plot by
condition with facet\_wrap(~condition). Also by adding a feature
parameter to the plotAbundance function we can split the boxplot by GvHD
status allowing us to investigate changes changes in subset abundance by
disease status or indeed, any other metadata feature slot.

``` r
cell_annotation_dif = metadataPlot(sce_gvhd,
    colby = 'condition',
    title = '',
    colkey = c(None = 'royalblue', GvHD = 'red2'),
    legendPosition = 'right',
    pointSize = 0.3,
    legendLabSize = 8,
    axisLabSize = 14,
    titleLabSize = 16,
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
    ncol = 2, align = "l", label_size = 20, rel_widths = c(1.6,1))
```

<img src="README_files/figure-gfm/cell annotations-1.png" style="display: block; margin: auto;" />

Finally we can perform a wilcoxon rank sum test to see if we can detect
statistically significant changes between our two experimental
conditions. diffClust() will perform the wilcoxon test on our data and
produce a results table along with Log fold change and the proportion of
each sample that each cluster represents. We can thes use diffPlot() to
generate a volcano plot, which shows us that there is a significant
difference between the abundance of B cells between our two conditions.
We can see that in our dataset it appears that B Cells are markedly
reduced in those patients that go on to develop GvHD.

``` r
 library(broom)
# This need the individual values
# This also need the right contrasts
diffClust_out = diffClust(sce_gvhd, 
                        clustering = 'cell_annotation', 
                        feature = 'condition')

head(diffClust_out[,c(1,14:18)])
```

    ##        cluster proportion proportion.1       logfc       pval      padj
    ## 1      B Cells  0.7571028    6.4246483 -2.13839808 0.01040562 0.1040562
    ## 2    Basophils  0.4749399    0.5833272 -0.20556014 0.52183939 0.8697323
    ## 3 CD4+ T Cells  8.4219643   12.6510706 -0.40689875 0.26233168 0.8697323
    ## 4    CD4+ Treg  1.1673934    0.7770172  0.40706623 0.87278012 0.8727801
    ## 5 CD8+ T Cells 18.5345116   18.8898577 -0.01899066 0.87278012 0.8727801
    ## 6          cDC  0.8772951    2.0589283 -0.85309741 0.07816909 0.3908454

``` r
gg_volcano = diffPlot(diffClust_out, p_val = "padj", threshold = 0.12)

abundance_bcell = plotAbundance(sce_gvhd,
              clusters = c("B Cells", "Basophils"),
              clusterAssign = 'cell_annotation',
              feature = 'condition',
              colkey = c(None = 'royalblue', GvHD = 'red2'),
              legendPosition = 'right', 
              stripLabSize = 10,
              axisLabSize = 14,
              titleLabSize = 12,
              subtitleLabSize = 10,
              captionLabSize = 18)

plot_grid(gg_volcano, abundance_bcell,
    labels = c('A','B'),
    ncol = 2, align = "l", label_size = 20, rel_widths = c(1,1.4))
```

<img src="README_files/figure-gfm/abundance test-1.png" style="display: block; margin: auto;" />

We can perform the same process of applying a wilcoxon test for
expression levels of markers, across all markers on all populations by
using the diffExpression() function. We can plot the results a volcano
plot as before using diffPlot(). Beside this can vizualise this cluster
specific expression of certain makers with markerExpressionPerSample()
to look at the median marker expression per sample on these clusters (or
all clusters\!).

``` r
# Need to catch errors - cluster not there causes crash
marker_test = diffExpression(sce_gvhd,
                assay = 'scaled',
                grouping = 'group', # The condensing function
                feature = 'condition', # The contrast
                clusterAssign = 'cell_annotation')

sig_markers = marker_test[which(marker_test$p_val < 0.05),]

head(sig_markers)
```

    ##                 cluster        S1         S2         S3         S4       S15
    ## 43            cDC_PDL_1 0.0000000 0.03833666 0.03454712 0.02409312 0.0000000
    ## 74     NKT Cells_CD45RA 2.9245811 2.82914412 1.35815724 1.62859842 1.3404085
    ## 88        NKT Cells_CD3 4.6261090 4.66563326 4.56496722 4.57661471 4.5837929
    ## 106 CD4+ T Cells_CD45RA 1.6022981 0.87746351 0.39035619 0.46861896 0.1388336
    ## 112   CD4+ T Cells_PD_1 0.1060011 0.30094336 0.17684411 0.35370009 0.7500233
    ## 133        B Cells_CCR7 1.4479768 1.90229903 1.99799022 2.32094205 1.0633981
    ##           S16        S23        S24       S25        S26       S27        S28
    ## 43  0.0000000 0.07893154 0.02029057 0.2797268 0.76921816 0.4672593 0.10298348
    ## 74  1.8112743 1.34732314 1.48795015 0.5681961         NA 0.8301662 0.74637372
    ## 88  4.8404246 4.57702149 4.48036593 4.2798380         NA 4.4999407 3.86990581
    ## 106 0.2393415 0.17488804 0.21700140 0.1083829 0.05085623 0.3204029 0.06000536
    ## 112 0.2040139 3.22044827 2.55991190 1.7497636 0.29420018 1.4976574 0.52886947
    ## 133 0.5621361 0.83526030 1.41891478 1.0311894 0.04672376 0.2211306 0.00000000
    ##        mean_1     mean_2       logfc      p_val      padj
    ## 43  0.2864016 0.01616282  2.87468184 0.01556764 0.5929411
    ## 74  0.9960019 1.98202728 -0.68812634 0.02845974 0.5929411
    ## 88  4.3414144 4.64292362 -0.06714407 0.01762209 0.5929411
    ## 106 0.1552561 0.61948531 -1.38381275 0.02497468 0.5929411
    ## 112 1.6418085 0.31525431  1.65017398 0.02497468 0.5929411
    ## 133 0.5922031 1.54912372 -0.96159499 0.02497468 0.5929411

``` r
# Create volcano plot
gg_volcano = diffPlot(marker_test, p_val = "p_val", threshold = 0.05)

# Difference in activation markers per cluster
pd_expression = markerExpressionPerSample(sce_gvhd,
    caption = '',
    clusters = c("CD8+ T Cells", "cDC"),
    feature = 'condition',
    clusterAssign = 'cell_annotation',
    markers = activation_markers[3:5],
    colkey = c('royalblue', 'red2'),
    title = '',
    stripLabSize = 12,
    axisLabSize = 10,
    titleLabSize = 18,
    subtitleLabSize = 10,
    captionLabSize = 10)


plot_grid(gg_volcano, pd_expression,
    labels = c('A','B'),
    ncol = 2, align = "l", label_size = 20, rel_widths = c(1,1.4))
```

<img src="README_files/figure-gfm/marker test-1.png" style="display: block; margin: auto;" />

# Subsetting the SCE object

The SCE object can be subsetted for further clustering, for instance, in
a situation when increased cluster resolution is desired or if in a
situiation where you want to diplay the proportion of a certain cluster
as a percentage of another parent populations, for instance CD4+ Tregs
as a proportion of all CD4+ T cells. Subsetting can be performed using
the subsetSCE() function. The function can subset cells on any specified
metadata slot using dplyr syntax conditional statements. An couple of
examples are presented below:

``` r
# Subset by patient_ID
sce_subset = subsetSCE(sce_gvhd, patient_id %in% c("P1", "P9", "P15"))

table(sce_subset@metadata$patient_id)
```

    ## 
    ##     P1    P15     P9 
    ## 156422 131263  98571

``` r
# Subset by clustering idenity
sce_cd4 = subsetSCE(sce_gvhd, cell_annotation %in% c("CD4+ T Cells", "Basophils"))

table(sce_cd4@metadata$cell_annotation)
```

    ## 
    ##    Basophils CD4+ T Cells 
    ##         5241       101908

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
    ##  [1] broom_0.7.0                 tibble_3.0.3               
    ##  [3] scran_1.14.6                umap_0.2.6.0               
    ##  [5] pheatmap_1.0.12             Rtsne_0.15                 
    ##  [7] reshape2_1.4.4              limma_3.42.2               
    ##  [9] FlowSOM_1.18.0              igraph_1.2.5               
    ## [11] ConsensusClusterPlus_1.50.0 RColorBrewer_1.1-2         
    ## [13] cowplot_1.0.0               stringr_1.4.0              
    ## [15] flowCore_1.52.1             immunoCluster_0.1.0        
    ## [17] dplyr_1.0.1                 ggrepel_0.8.2              
    ## [19] ggplot2_3.3.2               SingleCellExperiment_1.8.0 
    ## [21] SummarizedExperiment_1.16.1 DelayedArray_0.12.3        
    ## [23] BiocParallel_1.20.1         matrixStats_0.56.0         
    ## [25] Biobase_2.46.0              GenomicRanges_1.38.0       
    ## [27] GenomeInfoDb_1.22.1         IRanges_2.20.2             
    ## [29] S4Vectors_0.24.4            BiocGenerics_0.32.0        
    ## [31] kableExtra_1.1.0            knitr_1.28                 
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] backports_1.1.8          plyr_1.8.6               splines_3.6.2           
    ##   [4] fda_5.1.5.1              scater_1.14.6            digest_0.6.25           
    ##   [7] htmltools_0.4.0          viridis_0.5.1            fansi_0.4.1             
    ##  [10] magrittr_1.5             CytoML_1.12.1            cluster_2.1.0           
    ##  [13] ks_1.11.7                readr_1.3.1              RcppParallel_5.0.2      
    ##  [16] R.utils_2.9.2            askpass_1.1              flowWorkspace_3.34.1    
    ##  [19] jpeg_0.1-8.1             colorspace_1.4-1         rvest_0.3.5             
    ##  [22] rrcov_1.5-5              xfun_0.13                crayon_1.3.4            
    ##  [25] RCurl_1.98-1.2           jsonlite_1.7.0           hexbin_1.28.1           
    ##  [28] graph_1.64.0             glue_1.4.1               flowClust_3.24.0        
    ##  [31] gtable_0.3.0             zlibbioc_1.32.0          XVector_0.26.0          
    ##  [34] webshot_0.5.2            ggcyto_1.14.1            BiocSingular_1.2.2      
    ##  [37] IDPmisc_1.1.20           Rgraphviz_2.30.0         DEoptimR_1.0-8          
    ##  [40] scales_1.1.1             mvtnorm_1.1-1            edgeR_3.28.1            
    ##  [43] Rcpp_1.0.5               viridisLite_0.3.0        clue_0.3-57             
    ##  [46] dqrng_0.2.1              reticulate_1.16          openCyto_1.24.0         
    ##  [49] rsvd_1.0.3               mclust_5.4.6             tsne_0.1-3              
    ##  [52] httr_1.4.1               ellipsis_0.3.1           pkgconfig_2.0.3         
    ##  [55] XML_3.99-0.3             R.methodsS3_1.8.0        farver_2.0.3            
    ##  [58] flowViz_1.50.0           locfit_1.5-9.4           utf8_1.1.4              
    ##  [61] flowStats_3.44.0         tidyselect_1.1.0         labeling_0.3            
    ##  [64] rlang_0.4.7              munsell_0.5.0            tools_3.6.2             
    ##  [67] cli_2.0.2                generics_0.0.2           evaluate_0.14           
    ##  [70] yaml_2.2.1               robustbase_0.93-6        purrr_0.3.4             
    ##  [73] RBGL_1.62.1              R.oo_1.23.0              xml2_1.3.2              
    ##  [76] compiler_3.6.2           rstudioapi_0.11          beeswarm_0.2.3          
    ##  [79] png_0.1-7                statmod_1.4.34           pcaPP_1.9-73            
    ##  [82] stringi_1.4.6            RSpectra_0.16-0          lattice_0.20-41         
    ##  [85] Matrix_1.2-18            vctrs_0.3.2              pillar_1.4.6            
    ##  [88] lifecycle_0.2.0          BiocNeighbors_1.4.2      data.table_1.12.8       
    ##  [91] bitops_1.0-6             irlba_2.3.3              corpcor_1.6.9           
    ##  [94] R6_2.4.1                 latticeExtra_0.6-29      KernSmooth_2.23-17      
    ##  [97] gridExtra_2.3            vipor_0.4.5              MASS_7.3-51.6           
    ## [100] gtools_3.8.2             assertthat_0.2.1         openssl_1.4.2           
    ## [103] withr_2.2.0              mnormt_1.5-7             GenomeInfoDbData_1.2.2  
    ## [106] hms_0.5.3                ncdfFlow_2.32.0          grid_3.6.2              
    ## [109] tidyr_1.1.1              DelayedMatrixStats_1.8.0 rmarkdown_2.3           
    ## [112] base64enc_0.1-3          ggbeeswarm_0.6.0         ellipse_0.4.2

# Contact

For any queries relating to software:

  - James W Opzoomer (<james.opzoomer@kcl.ac.uk>)
