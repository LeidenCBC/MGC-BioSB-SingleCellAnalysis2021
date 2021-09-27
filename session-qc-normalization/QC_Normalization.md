Quality Control and Normalization
================

Created by: Ahmed Mahfouz  
Edited by: Mohammed Charrout, Lieke Michielsen

# Overview

In this practical, we will walk through a pipeline to analyze single
cell RNA-sequencing (scRNA-seq) data. Starting from a count matrix, we
will cover the following steps of the analysis:

1.  Quality control
2.  Normalization
3.  Feature selection

## Datasets

For this tutorial we will use 3 different PBMC datasets from the 10x
Genomics website
(<https://support.10xgenomics.com/single-cell-gene-expression/datasets>).

-   1k PBMCs using 10x v2 chemistry
-   1k PBMCs using 10x v3 chemistry
-   1k PBMCs using 10x v3 chemistry in combination with cell surface
    proteins, but disregarding the protein data and only looking at gene
    expression.

The datasets are available in this repository.

Load required packages:

``` r
library(Seurat)
library(scater)
library(scran)
library(Matrix)
```

# Read the data and create a Seurat object

Here, we use the function `Read10X_h5` of the Seurat package to read in
the expression matrices. R stores these matrices as sparse matrix
objects, which are essentially memory-efficient tables of values. In
this case the values represent the RNA counts in each cell.

``` r
v3.1k <- Read10X_h5("pbmc_1k_v3_filtered_feature_bc_matrix.h5")
v2.1k <- Read10X_h5("pbmc_1k_v2_filtered_feature_bc_matrix.h5")
p3.1k <- Read10X_h5("pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5")
```

    ## Genome matrix has multiple modalities, returning a list of matrices for this genome

``` r
# select only gene expression data from the CITE-seq data.
p3.1k <- p3.1k$`Gene Expression`
```

Rather than working directly with matrices, Seurat works with custom
objects that wrap around them. These Seurat objects also conveniently
contain tables of metadata for the cells and features, which avoids the
clutter of managing them as separate objects. As we will see later,
normalized expression values are stored in a separate matrix within the
Seurat object, which allows us to play around with different
normalization strategies without manually keeping a backup of the
original values. In addition to RNA counts, we are able to store
additional data types (called *assays*) within the Seurat object, such
as protein measurements measured by CITE-seq, though we will stick to
the default RNA assay here.

First, create Seurat objects for each of the datasets, and then merge
into one large Seurat object. We will use the cell metadata to keep
track of which dataset the cell originated from.

``` r
sdata.v2.1k <- CreateSeuratObject(v2.1k, project = "v2.1k")
sdata.v3.1k <- CreateSeuratObject(v3.1k, project = "v3.1k")
sdata.p3.1k <- CreateSeuratObject(p3.1k, project = "p3.1k")
# Merge into one single Seurat object. 
# Prefix cell ids with dataset name (`all.cell.ids`) just in case you have 
# overlapping barcodes between the datasets.
alldata <- merge(sdata.v2.1k, c(sdata.v3.1k, sdata.p3.1k), add.cell.ids=c("v2.1k","v3.1k","p3.1k"))
# Also add in a metadata column that indicates v2 vs v3 chemistry.
chemistry <- rep("v3", ncol(alldata))
chemistry[Idents(alldata) == "v2.1k"] <- "v2"
alldata <- AddMetaData(alldata, chemistry, col.name = "Chemistry")
alldata
```

    ## An object of class Seurat 
    ## 33538 features across 2931 samples within 1 assay 
    ## Active assay: RNA (33538 features, 0 variable features)

The metadata of the Seurat object, which itself is a data frame, can be
accessed using the slot operator (`@`) like so `alldata@meta.data`.
Alternatively one can call the object with double empty square brackets:
`alldata[[]]`. Another slot to be aware of is `alldata@active.ident`, or
alternatively `Idents(alldata)`, which stores a column of the metadata
that should be used to identify groups of cells. The value of the
identities are by default chosen to be whatever is passed to the
`project` parameter in the `CreateSeuratObject` call, and is stored in
the `orig.ident` column of the metadata object. We are free to change
the column that represent the cell identities but for this tutorial (and
in the general case) we keep it as is.

Let’s check number of cells from each sample using the idents.

``` r
table(Idents(alldata))
```

    ## 
    ## p3.1k v2.1k v3.1k 
    ##   713   996  1222

## 1. Quality control

On object creation, Seurat automatically calculates some QC-stats such
as the number of UMIs and features per cell. This information is stored
in the columns `nCount_RNA` and `nFeature_RNA` of the metadata.

``` r
head(alldata@meta.data)
```

    ##                          orig.ident nCount_RNA nFeature_RNA Chemistry
    ## v2.1k_AAACCTGAGCGCTCCA-1      v2.1k       6631         2029        v2
    ## v2.1k_AAACCTGGTGATAAAC-1      v2.1k       2196          881        v2
    ## v2.1k_AAACGGGGTTTGTGTG-1      v2.1k       2700          791        v2
    ## v2.1k_AAAGATGAGTACTTGC-1      v2.1k       3551         1183        v2
    ## v2.1k_AAAGCAAGTCTCTTAT-1      v2.1k       3080         1333        v2
    ## v2.1k_AAAGCAATCCACGAAT-1      v2.1k       5769         1556        v2

Note that the `_RNA` suffix is due to the aforementioned potential to
hold multiple assays. The default assay is named `RNA`, accessible by
the assays slot `alldata@assays$RNA`, which is by default set to be the
standard active assay (see `alldata@active.assay`). Effectively this
means that any calls that are done on the Seurat object are applied on
the `RNA` assay data.

### Calculate mitochondrial proportion

We will manually calculate the proportion of mitochondrial reads and add
it to the metadata table. Mitochondrial genes start with a `MT-` prefix.

``` r
percent.mito <- PercentageFeatureSet(alldata, pattern = "^MT-")
alldata <- AddMetaData(alldata, percent.mito, col.name = "percent.mito")
```

#### Calculate ribosomal proportion

In the same manner we will calculate the proportion of the counts that
come from ribosomal proteins, identified by the `RPS` and `RPL`
prefixes.

``` r
percent.ribo <- PercentageFeatureSet(alldata, pattern = "^RP[SL]")
alldata <- AddMetaData(alldata, percent.ribo, col.name = "percent.ribo")
```

Now have another look at the metadata table.

``` r
head(alldata@meta.data)
```

    ##                          orig.ident nCount_RNA nFeature_RNA Chemistry
    ## v2.1k_AAACCTGAGCGCTCCA-1      v2.1k       6631         2029        v2
    ## v2.1k_AAACCTGGTGATAAAC-1      v2.1k       2196          881        v2
    ## v2.1k_AAACGGGGTTTGTGTG-1      v2.1k       2700          791        v2
    ## v2.1k_AAAGATGAGTACTTGC-1      v2.1k       3551         1183        v2
    ## v2.1k_AAAGCAAGTCTCTTAT-1      v2.1k       3080         1333        v2
    ## v2.1k_AAAGCAATCCACGAAT-1      v2.1k       5769         1556        v2
    ##                          percent.mito percent.ribo
    ## v2.1k_AAACCTGAGCGCTCCA-1     5.172674     25.84829
    ## v2.1k_AAACCTGGTGATAAAC-1     4.143898     20.81056
    ## v2.1k_AAACGGGGTTTGTGTG-1     3.296296     51.55556
    ## v2.1k_AAAGATGAGTACTTGC-1     5.885666     29.25936
    ## v2.1k_AAAGCAAGTCTCTTAT-1     2.987013     17.53247
    ## v2.1k_AAAGCAATCCACGAAT-1     2.010747     45.69249

#### Plot QC

Now we can plot some of the QC-features as violin plots. Note that
Seurat by default will generate a violin plot per identity class.

``` r
VlnPlot(alldata, features = c("nFeature_RNA", "nCount_RNA", "percent.mito", "percent.ribo"), 
        ncol = 2, pt.size = 0.1) + NoLegend()
```

![](QC_Normalization_files/figure-gfm/plot_qc-1.png)<!-- -->

As you can see, the v2 chemistry gives lower gene detection, but higher
detection of ribosomal proteins. As the ribosomal proteins are highly
expressed they will make up a larger proportion of the transcriptional
landscape when fewer of the lowly expressed genes are detected.

We can also plot the different QC-measures as scatter plots.

``` r
p1 <- FeatureScatter(alldata, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + NoLegend()
p2 <- FeatureScatter(alldata, feature1 = "nFeature_RNA", feature2 = "percent.mito") + NoLegend()
p3 <- FeatureScatter(alldata, feature1="percent.ribo", feature2="nFeature_RNA")
p1 + p2 + p3
```

![](QC_Normalization_files/figure-gfm/plot_scatter1-1.png)<!-- -->

We can also subset the data to only plot one sample.

``` r
FeatureScatter(alldata, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", 
               cells = WhichCells(alldata, expression = orig.ident == "v3.1k") )
```

![](QC_Normalization_files/figure-gfm/plot_scatter2-1.png)<!-- -->

### Filtering

#### Mitochondrial filtering

We have quite a lot of cells with high proportion of mitochondrial
reads. It could be wise to remove those cells, if we have enough cells
left after filtering. Another option would be to either remove all
mitochondrial reads from the dataset and hope that the remaining genes
still have enough biological signal. A third option would be to just
regress out the `percent.mito` variable during scaling.

In this case we have as much as 99.7% mitochondrial reads in some of the
cells, so it is quite unlikely that there is much cell type signature
left in those.

By eyeballing the plots we can make reasonable decisions on where to
draw the cutoff. In this case, the bulk of the cells are below 25%
mitochondrial reads and that will be used as a cutoff.

``` r
# Select cells with percent.mito < 25
idx <- which(alldata$percent.mito < 25)
selected <- WhichCells(alldata, cells = idx)
length(selected)
```

    ## [1] 2703

``` r
# and subset the object to only keep those cells.
data.filt <- subset(alldata, cells = selected)
# plot violins for new data
VlnPlot(data.filt, features = "percent.mito")
```

![](QC_Normalization_files/figure-gfm/mito.filt-1.png)<!-- -->

As you can see, there is still quite a lot of variation in percent mito,
so it will have to be dealt with in the data analysis step.

#### Gene detection filtering

Extremely high number of detected genes could indicate doublets.
However, depending on the cell type composition in your sample, you may
have cells with higher number of genes (and also higher counts) from one
cell type.

In our datasets, we observe a clear difference between the v2 vs v3 10x
chemistry with regards to gene detection, so it may not be fair to apply
the same cutoffs to all of them.

Also, in the protein assay data there is a lot of cells with few
detected genes giving a bimodal distribution. This type of distribution
is not seen in the other 2 datasets. Considering that they are all pbmc
datasets it makes sense to regard this distribution as low quality
libraries.

Filter the cells with high gene detection (putative doublets) with
cutoffs 4100 for v3 chemistry and 2000 for v2.

``` r
# Start with cells with many genes detected.
high.det.v3 <- WhichCells(data.filt, expression = nFeature_RNA > 4100)
high.det.v2 <- WhichCells(data.filt, expression = nFeature_RNA > 2000 & orig.ident == "v2.1k")
# Remove these cells.
data.filt <- subset(data.filt, cells=setdiff(WhichCells(data.filt),c(high.det.v2,high.det.v3)))
# Check number of cells.
ncol(data.filt)
```

    ## [1] 2631

Filter the cells with low gene detection (low quality libraries) with
less than 1000 genes for v2 and &lt; 500 for v2.

``` r
#start with cells with many genes detected.
low.det.v3 <- WhichCells(data.filt, expression = nFeature_RNA < 1000 & orig.ident != "v2.1k")
low.det.v2 <- WhichCells(data.filt, expression = nFeature_RNA < 500 & orig.ident == "v2.1k")
# remove these cells
data.filt <- subset(data.filt, cells=setdiff(WhichCells(data.filt),c(low.det.v2,low.det.v3)))
# check number of cells
ncol(data.filt)
```

    ## [1] 2531

#### Plot QC-stats again

Lets plot the same qc-stats another time.

``` r
VlnPlot(data.filt, features = c("nFeature_RNA", "nCount_RNA", "percent.mito", "percent.ribo"), 
        ncol = 2, pt.size = 0.1) + NoLegend()
```

![](QC_Normalization_files/figure-gfm/vln.plot2-1.png)<!-- -->

``` r
# and check the number of cells per sample before and after filtering
table(Idents(alldata))
```

    ## 
    ## p3.1k v2.1k v3.1k 
    ##   713   996  1222

``` r
table(Idents(data.filt))
```

    ## 
    ## p3.1k v2.1k v3.1k 
    ##   526   933  1072

### Calculate cell-cycle scores

Seurat has a function for calculating cell cycle scores based on a list
of know S-phase and G2/M-phase genes.

``` r
data.filt <- CellCycleScoring(
  object = data.filt,
  g2m.features = cc.genes$g2m.genes,
  s.features = cc.genes$s.genes
)
```

    ## Warning: The following features are not present in the object: MLF1IP, not
    ## searching for symbol synonyms

    ## Warning: The following features are not present in the object: FAM64A, HN1, not
    ## searching for symbol synonyms

``` r
VlnPlot(data.filt, features = c("S.Score","G2M.Score"))
```

![](QC_Normalization_files/figure-gfm/cc-1.png)<!-- -->

In this case it looks like we only have a few cycling cells in the
datasets.