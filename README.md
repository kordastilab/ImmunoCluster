immunoCluster
================
James Opzoomer, Kevin Blighe, Jessica Timms
2020-03-30

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
```

## 2\. Load the package and the dependancies into the R session

``` r
  library(immunnoCluster)

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

Aternately you can download the imported RDS directly to skip the import
step and save time:

``` r
# Insert zenodo RDS link
```

# Contact

For any queries relating to software:

  - James W Opzoomer (<james.opzoomer@kcl.ac.uk>)
