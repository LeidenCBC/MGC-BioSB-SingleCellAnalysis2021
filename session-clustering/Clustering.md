Clustering
================

Created by: Ahmed Mahfouz

Edited by: Mohammed Charrout, Lieke Michielsen

# Overview

In this tutorial we will look at different approaches to cluster
scRNA-seq datasets in order to characterize the different subgroups of
cells. Using unsupervised clustering, we will try to identify groups of
cells based on the similarities of the transcriptomes without any prior
knowledge of the labels.

Load required packages:

``` r
suppressMessages(require(Seurat))
```

## Datasets

Again, we will continue with dataset the you have preprocessed and
visualized in the previous practicals. Let’s start by loading the data
again.

During the previous practical, we have already selected highly variable
genes. This step is also to decide which genes to use when clustering
the cells. Single cell RNA-seq can profile a huge number of genes in a
lot of cells. But most of the genes are not expressed enough to provide
a meaningful signal and are often driven by technical noise. Including
them could potentially add some unwanted signal that would blur the
biological variation. Moreover gene filtering can also speed up the
computational time for downstream analysis.

``` r
pbmc <- readRDS('../session-dimensionalityreduction/pbmc3k.rds')
```

## Clustering

### Hierarchical clustering

``` r
# Get scaled counts from the Seurat object
scaled_pbmc <- pbmc@assays$RNA@scale.data

# Calculate Distances (default: Euclidean distance)
distance_euclidean <- dist(t(scaled_pbmc))

#Perform hierarchical clustering using ward linkage
ward_hclust_euclidean <- hclust(distance_euclidean,method = "ward.D2")
plot(ward_hclust_euclidean, main = "dist = eucledian, Ward linkage", labels=FALSE)
```

![](Clustering_files/figure-gfm/hierarchical_eucledian_ward-1.png)<!-- -->

Now cut the dendrogram to generate 10 clusters and plot the cluster
labels and the previously given celltype labels on the t-SNE plot. For
now, we just pick 10, but you can of course vary this number to see how
it influences your results.

``` r
#Cutting the cluster tree to make 10 groups
cluster_hclust <- cutree(ward_hclust_euclidean,k = 10)
pbmc@meta.data$cluster_hclust <- factor(cluster_hclust)

p1 <- DimPlot(pbmc, reduction="tsne", group.by = "cluster_hclust")
p2 <- DimPlot(pbmc, reduction="tsne", group.by = "celltype")

p1+p2
```

![](Clustering_files/figure-gfm/hierarchical_eucledian_ward_pcaplot-1.png)<!-- -->

Now let’s try a different distance measure. A commonly used distance
measure is 1 - correlation.

``` r
# Calculate Distances (1 - correlation)
C <- cor(scaled_pbmc)

# Run clustering based on the correlations, where the distance will 
# be 1-correlation, e.g. higher distance with lower correlation.
distance_corr <- as.dist(1-C) 
    
#Perform hierarchical clustering using ward linkage
ward_hclust_corr <- hclust(distance_corr,method="ward.D2")
plot(ward_hclust_corr, main = "dist = 1-corr, Ward linkage", labels=FALSE)
```

![](Clustering_files/figure-gfm/hierarchical_corr_ward-1.png)<!-- -->

Again, let’s cut the dendrogram to generate 10 clusters and plot the
cluster labels on the t-SNE plot.

``` r
#Cutting the cluster tree to make 10 groups
cluster_hclust <- cutree(ward_hclust_corr,k = 10)
pbmc@meta.data$cluster_hclust <- factor(cluster_hclust)

p1 <- DimPlot(pbmc, reduction="tsne", group.by = "cluster_hclust")
p2 <- DimPlot(pbmc, reduction="tsne", group.by = "celltype")

p1+p2
```

![](Clustering_files/figure-gfm/hierarchical_corr_ward_pcaplot-1.png)<!-- -->

Instead of changing the distance metric, we can change the linkage
method. Instead of using Ward’s method, let’s use complete linkage.

``` r
#Perform hierarchical clustering using complete linkage & euclidean distance
comp_hclust_eucledian <- hclust(distance_euclidean,method = "complete")
plot(comp_hclust_eucledian, main = "dist = euclidean, complete linkage", labels=FALSE)
```

![](Clustering_files/figure-gfm/hierarchical_eucledian_complete-1.png)<!-- -->

Once more, let’s cut the dendrogram to generate 10 clusters and plot the
cluster labels on the t-SNE plot.

``` r
#Cutting the cluster tree to make 10 groups
cluster_hclust <- cutree(comp_hclust_eucledian,k = 10)
pbmc@meta.data$cluster_hclust <- factor(cluster_hclust)

p1 <- DimPlot(pbmc, reduction="tsne", group.by = "cluster_hclust")
p2 <- DimPlot(pbmc, reduction="tsne", group.by = "celltype")

p1+p2
```

![](Clustering_files/figure-gfm/hierarchical_eucledian_complete_pcaplot-1.png)<!-- -->
As you can see, these linkage methods and distances cluster the data
differently. If you want, there are even more distance measures and
linkage methods to play around with.

### K-means

Next, we will try the k-means algorithm on the scaled data.

``` r
pbmc_kmeans <- kmeans(x = t(scaled_pbmc), centers = 10)
pbmc@meta.data$cluster_kmeans <- factor(pbmc_kmeans$cluster)

p1 <- DimPlot(pbmc, reduction="tsne", group.by = "cluster_kmeans")
p2 <- DimPlot(pbmc, reduction="tsne", group.by = "celltype")

p1+p2
```

![](Clustering_files/figure-gfm/kmeans-1.png)<!-- -->

### Graph based clustering

The clustering algorithm of Seurat itself is based on graph based
clustering. The output of the clustering, will be saved automatically in
the metadata as ‘seurat\_clusters’. As explained in the lecture, the
resolution parameter is related to the number of clusters. You can play
around with this parameters to see how it influences the results.

``` r
pbmc <- FindNeighbors(pbmc, dims = 1:10, verbose = FALSE)
pbmc <- FindClusters(pbmc, resolution = 0.25, verbose = FALSE)

p1 <- DimPlot(pbmc, reduction="tsne", group.by = "seurat_clusters")
p2 <- DimPlot(pbmc, reduction="tsne", group.by = "celltype")

p1+p2
```

![](Clustering_files/figure-gfm/graph_clust-1.png)<!-- -->

## Visualizing marker genes and annotating the cells

Once, you are satisfied with the clusters, these can be annotated by
visualizing known marker genes or by looking at differentially expressed
genes. In a later practical, you will learn how to select these, for now
we will just focus on known marker genes. A commonly used approach is
that the data is annotated in a hierarchical fashion. First the data is
annotated at a low resolution (e.g. only 2-3 cell types) and afterwards
each cluster is subsetted from the data, clustered and annotated again.
This process can continue until you’re satisfied with the resolution.

``` r
pbmc <- FindNeighbors(pbmc, dims = 1:10, verbose = FALSE)
pbmc <- FindClusters(pbmc, resolution = 0.01, verbose = FALSE)

p1 <- DimPlot(pbmc, reduction="tsne", group.by = "seurat_clusters")
p2 <- DimPlot(pbmc, reduction="tsne", group.by = "celltype")

p1+p2
```

![](Clustering_files/figure-gfm/graph_clust_lowres-1.png)<!-- -->

So now that we have clustered the data at a low resolution, we can
visualize some marker genes: CD19 (B cells), CD3D (T cells), CD14
(Monocytes), NKG7 (NK cells).

``` r
FeaturePlot(pbmc, reduction='tsne', features=c('CD19', 'CD3D', 'CD14', 'NKG7'))
```

![](Clustering_files/figure-gfm/featplot-1.png)<!-- -->

For a new, more complex dataset, you will probably need to visualize
more genes before you can label a cluster. For now, we will assume that
cluster 0 are NK and T cells, cluster 1 are Monocytes and cluster 2 are
B cells. In the code below, you will assign these labels to your
cluster.

``` r
new.cluster.ids <- c("NK and T cells", "Monocytes", "B cells")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
DimPlot(pbmc, reduction = "tsne", label = TRUE) + NoLegend()
```

![](Clustering_files/figure-gfm/clusnames-1.png)<!-- -->

If you want to cluster the cells at a higher resolution, you could for
instance subset the data now and repeat these steps. For now, we will
just save the object for the next practicals.

``` r
saveRDS(pbmc, file = "pbmc3k.rds")
```

### Session info

``` r
sessionInfo()
```

    ## R version 4.0.5 (2021-03-31)
    ## Platform: x86_64-w64-mingw32/x64 (64-bit)
    ## Running under: Windows 10 x64 (build 18363)
    ## 
    ## Matrix products: default
    ## 
    ## locale:
    ## [1] LC_COLLATE=Dutch_Netherlands.1252  LC_CTYPE=Dutch_Netherlands.1252   
    ## [3] LC_MONETARY=Dutch_Netherlands.1252 LC_NUMERIC=C                      
    ## [5] LC_TIME=Dutch_Netherlands.1252    
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ## [1] SeuratObject_4.0.2 Seurat_4.0.4      
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] nlme_3.1-152          spatstat.sparse_2.0-0 matrixStats_0.61.0   
    ##   [4] RcppAnnoy_0.0.19      RColorBrewer_1.1-2    httr_1.4.2           
    ##   [7] sctransform_0.3.2     tools_4.0.5           utf8_1.2.2           
    ##  [10] R6_2.5.1              irlba_2.3.3           rpart_4.1-15         
    ##  [13] KernSmooth_2.23-18    uwot_0.1.10           mgcv_1.8-34          
    ##  [16] lazyeval_0.2.2        colorspace_2.0-2      tidyselect_1.1.1     
    ##  [19] gridExtra_2.3         compiler_4.0.5        plotly_4.9.4.1       
    ##  [22] labeling_0.4.2        scales_1.1.1          lmtest_0.9-38        
    ##  [25] spatstat.data_2.1-0   ggridges_0.5.3        pbapply_1.5-0        
    ##  [28] goftest_1.2-2         stringr_1.4.0         digest_0.6.28        
    ##  [31] spatstat.utils_2.2-0  rmarkdown_2.11        pkgconfig_2.0.3      
    ##  [34] htmltools_0.5.2       parallelly_1.28.1     highr_0.9            
    ##  [37] fastmap_1.1.0         htmlwidgets_1.5.4     rlang_0.4.11         
    ##  [40] shiny_1.7.0           farver_2.1.0          generics_0.1.0       
    ##  [43] zoo_1.8-9             jsonlite_1.7.2        ica_1.0-2            
    ##  [46] dplyr_1.0.7           magrittr_2.0.1        patchwork_1.1.1      
    ##  [49] Matrix_1.3-4          Rcpp_1.0.7            munsell_0.5.0        
    ##  [52] fansi_0.5.0           abind_1.4-5           reticulate_1.22      
    ##  [55] lifecycle_1.0.1       stringi_1.7.4         yaml_2.2.1           
    ##  [58] MASS_7.3-54           Rtsne_0.15            plyr_1.8.6           
    ##  [61] grid_4.0.5            parallel_4.0.5        listenv_0.8.0        
    ##  [64] promises_1.2.0.1      ggrepel_0.9.1         crayon_1.4.1         
    ##  [67] deldir_0.2-10         miniUI_0.1.1.1        lattice_0.20-41      
    ##  [70] cowplot_1.1.1         splines_4.0.5         tensor_1.5           
    ##  [73] knitr_1.36            pillar_1.6.3          igraph_1.2.6         
    ##  [76] spatstat.geom_2.2-2   future.apply_1.8.1    reshape2_1.4.4       
    ##  [79] codetools_0.2-18      leiden_0.3.9          glue_1.4.2           
    ##  [82] evaluate_0.14         data.table_1.14.2     png_0.1-7            
    ##  [85] vctrs_0.3.8           httpuv_1.6.3          polyclip_1.10-0      
    ##  [88] gtable_0.3.0          RANN_2.6.1            purrr_0.3.4          
    ##  [91] spatstat.core_2.3-0   tidyr_1.1.4           scattermore_0.7      
    ##  [94] future_1.22.1         ggplot2_3.3.5         xfun_0.26            
    ##  [97] mime_0.12             xtable_1.8-4          later_1.3.0          
    ## [100] survival_3.2-10       viridisLite_0.4.0     tibble_3.1.4         
    ## [103] cluster_2.1.1         globals_0.14.0        fitdistrplus_1.1-6   
    ## [106] ellipsis_0.3.2        ROCR_1.0-11
