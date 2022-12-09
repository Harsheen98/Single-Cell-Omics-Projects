### Loading libraries
library(DropletUtils)
library(ggplot2)
library(scater)
library(Seurat)
library(dplyr)
library(patchwork)

############# LOADING THE GENE EXPRESSION MATRIX #######################

df_75688 <- read.table('Data/GSE75688_GEO_processed_Breast_Cancer_raw_TPM_matrix.txt')  # loading the data
df_75688 <- t(df_75688)                                                                 # transposing
class(df_75688)                                                                         # checking the class, it is matrix
df_75688 <- as.data.frame(df_75688)                                                     # changing it to dataframe
class(df_75688)                                                                         
colnames(df_75688) <- df_75688[2,]                                                      # making the gene name column as the header
df_75688 <- df_75688[-c(1:3),]                                                          # removing the first three rows ; gene id, gene name(made into header), gene type
colnames(df_75688)[colnames(df_75688) == 'gene_name'] <- "!Sample_title"                # Changing the column for ease of merging
dim(df_75688)

############# LOADING THE METADATA #####################

df_75688_meta <- read.table('GSE75688_series_matrix.txt')                               # loading the metadata file
df_75688_meta <- t(df_75688_meta)                                                       # transposing the file
df_75688_meta <- as.data.frame(df_75688_meta)                                           # converting the matrix datatype into dataframe
colnames(df_75688_meta) <- df_75688_meta[1,]                                            # making the first row as the header
df_75688_meta <- df_75688_meta[-1,]                                                     # deleting the first row since it is now header
dim(df_75688_meta)

############# SINCE THE TWO DATAFRAMES ARE NOW THE SAME DIMENSIONS, WE ARE GOING TO MERGE BY THE SAMPLE TITLE COLUMN ##############

df_merged <- merge(df_75688, df_75688_meta, by='!Sample_title')
dim(df_merged)
write.table(df_merged, 'temp_.tsv', sep='\t', row.names = FALSE)                         # save the file

# filter the dataframe using terminal to remove samples corresponding to pooled rna seq, tumor rna seq, er+, her2+, double positive breast cancer
# load the filtered df 
df_filtered <- read.table('temp4_.tsv', header=T, sep = '\t')            
df_filtered2 <- df_filtered[-c(57822:57918)]                                             # removed the ercc spike in columns and sample char, sample type columns  
dim(df_filtered2)

############ PREPROCESSING AND QC ###################

df_filtered2 <- t(df_filtered2)
colnames(df_filtered2) <- df_filtered2[1,]
df_filtered2 <- df_filtered2[-1,]

# 1. Removing cells with zero total umi count / empty barcodes
# 1.1 first use package dropletutils to rank the barcodes
df_filtered2 <- as.data.frame(df_filtered2)                                              # converting matrix into dataframe
str(df_filtered2)                                                                        # values in the dataframe are characters
df_filtered2[] <- sapply(df_filtered2,as.numeric)                                        # convert into numeric to use barcoderanks function
str(df_filtered2)        

barcode_ranks <- barcodeRanks(df_filtered2)
class(barcode_ranks)                                                                    # before plotting with ggplot we need to convert into dataframe
barcode_ranks_dum <- as.data.frame(barcode_ranks)
class(barcode_ranks_dum)
# making barcode ranks plot 
br_plot <- ggplot(barcode_ranks_dum, aes(x=rank, y=total)) + geom_point(alpha=0.5)
br_plot <- br_plot + geom_hline(aes(yintercept=694682, linetype='knee'), color='red')                 # add knee point ;  save it into a variable so that we can add another h-line
br_plot <- br_plot + geom_hline(aes(yintercept=177674, linetype='inflection'))                        # add inflection point
br_plot + xlab("Rank") + ylab("Total UMI Count")

# 1.2 Use emptyDrops function to remove the empty droplets
# for the empty drops function we need to decide upon the value of ambient gene expression
# for this I will be plotting a density plot for the total column in the barcode ranks DFrame
tot_umi_plt <- ggplot(barcode_ranks_dum, aes(x=total)) + geom_density()
tot_umi_plt + ylab('Density') + xlab('Total UMI Count') + geom_vline(xintercept=625000, linetype='dashed', color='red')  # selected 625000 

empty_drops_df <- emptyDrops(df_filtered2, lower = 625000)
empty_drops_df_dum <- as.data.frame(empty_drops_df)                                                   # also keep a dataframe version

# we need to now drop the index/samples with NA value from our gene expression matrix
index_dump <- which(is.na(empty_drops_df_dum), arr.ind=TRUE)                                          # this returns us the samples with NA value once for all columns, i.e., samples with NA in Total, then in LogProb, PValue, Limited, FDR                             
index_dump <- as.data.frame(index_dump)                                                               # slice for only the samples in the first column, because after that it is repeating
index_dump <- index_dump[c(1:21),]
sample_drop <- index_dump$row
df_filtered3 <- df_filtered2[-c(sample_drop)]                                                         # dataframe with the empty droplets (samples with ambient gene expression dropped)
dim(df_filtered3)

# 2. Dropping low quality cells
# 2.1 filtering cells with more than 5% mitochondrial gene proportion
is.mito =  grep("MT.", rownames(df_filtered3))                                                        # selecting mitochondrial genes(rows)
per.cell <- perCellQCMetrics(df_filtered3, subsets=list(Mito=is.mito))
per.cell
per.cell_df <- as.data.frame(per.cell)                                                                # converting into df for plotting
mito_prop_plot <- ggplot(per.cell_df, aes(x=subsets_Mito_percent)) + geom_histogram(binwidth = 1,     # plotting the mitochondrial proportion in each cell
                                                  color = 'black', fill='white')
mito_prop_plot + ylab('Number of cells') + xlab('Mitochondrial prop. (%)')
# perform filtering
qc.stats <- perCellQCFilters(per.cell, sub.fields="subsets_Mito_percent")                             # filtering the cells with more than 5% mt gene percentage
discard.index <- which(qc.stats$discard, arr.ind = FALSE, useNames = TRUE)                            # to pick out the indices that have TRUE value in the discard column
df_filtered4 <- df_filtered3[-c(discard.index)]
dim(df_filtered4)
# after filtering
is.mito2 <- grep("MT.", rownames(df_filtered4))
per.cell2 <- perCellQCMetrics(df_filtered4, subsets=list(Mito=is.mito))
per.cell2
per.cell2_df <- as.data.frame(per.cell2)
mito_prop_plot2 <- ggplot(per.cell2_df, aes(x=subsets_Mito_percent)) + geom_histogram(binwidth = 1,    
                                                                                     color = 'black', fill='white')
mito_prop_plot2 + ylab('Number of cells') + xlab('Mitochondrial prop. (%)')

# making plot of top 20 most expressed genes
a <- as.data.frame(t(df_filtered4))
umi_counts_per_gene <- perCellQCMetrics(a)
umi_counts_per_gene_df <- as.data.frame(umi_counts_per_gene)
top_genes <- umi_counts_per_gene_df[order(umi_counts_per_gene_df$sum, decreasing = TRUE),]
top_20_genes <- top_genes[c(1:20),]
top_20_genes <- tibble::rownames_to_column(top_20_genes, "gene")                                      # converting the index into column with gene name for ease of plotting
ggplot(top_20_genes, aes(x=sum, y=reorder(gene,sum))) + geom_bar(stat="identity") 
+ labs(x="Total UMI counts", y = "Top expressed genes")

# save the Final Quality controlled dataset into csv
GSE_75688_QC <- write.csv(df_filtered4 , sep = '\t', header=T)

# 3. Normalization of the Expression Matrix
seurat.obj <- CreateSeuratObject(df_filtered4)
seurat.obj
str(seurat.obj)
seurat.obj <- NormalizeData(seurat.obj)
str(seurat.obj)

############### PRINCIPAL COMPONENT ANALYSIS (PCA) #######################

# 1. Feature Selection ; Identify highly variable genes
seurat.obj <- FindVariableFeatures(seurat.obj)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(seurat.obj), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(seurat.obj)
LabelPoints(plot = plot1, points = top10, repel = TRUE)

# 2. Scaling the data before PCA
all.genes <- rownames(seurat.obj)
seurat.obj <- ScaleData(seurat.obj, features = all.genes)

# 3. Perform Linear dimensionality reduction using PCA
seurat.obj <- RunPCA(seurat.obj)

# Visualise PCA results
seurat.obj[['pca']]
print(seurat.obj[['pca']], dims = 1:5, nfeatures = 5)                                     # this is to visualise only the first 5 pcs
DimHeatmap(seurat.obj, dims= 1:9, cells=184, balanced=TRUE)
head(Embeddings(seurat.obj, reduction='pca')[,1:5])                                       # cell embeddings
head(Loadings(seurat.obj, reduction = "pca")[,1:5])                                       # feature loadings

# visualising the loadings
VizDimLoadings(seurat.obj, dims = 1:2, reduction = "pca") 
VizDimLoadings(seurat.obj, dims = 3:4, reduction = "pca")                                 # not visualising all 4 PC loadings in one go since the output plot is getting too cluttered

# determine the dimensionality of the data by selecting the PCs which carry approx 70% variance of the data
ElbowPlot(seurat.obj)                                                                     # selecting 15 PCs for neighborhood graph construction

# visualising clusters using UMAP
seurat.obj <- FindNeighbors(seurat.obj, reduction = 'pca', dims = 1:15)                    
seurat.obj <- FindClusters(seurat.obj)
seurat.obj <- runUMAP(seurat.obj, dims=1:15)
p1 <- DimPlot(seurat.obj, reduction = 'umap', group.by = 'obj.ident')

################# LOADING THE GENE EXPRESSION MATRIX FOR THE GSE_118389 DATASET ################

# this dataset is already quality controlled as is mentioned in the Paper, so we will the QC
df_118389 <- read.table('Data/GSE118389_tpm_rsem.txt', sep = '\t', header = T)
class(df_118389)
dim(df_118389)

# loading the df into the seurat object
seurat.obj.2 <- CreateSeuratObject(df_118389)
seurat.obj.2
str(seurat.obj.2)

# Normalization
seurat.obj.2 <- NormalizeData(seurat.obj.2)

############### PRINCIPAL COMPONENT ANALYSIS ###################

# 1. Feature Selection
seurat.obj.2 <- FindVariableFeatures(seurat.obj.2)

# Identify the 10 most important genes
top10.2 <- head(VariableFeatures(seurat.obj.2), 10)

# plot variable features with or without labels
plot2 <- VariableFeaturePlot(seurat.obj.2)
LabelPoints(plot=plot2, points=top10.2, repel=T)

# 2. Scaling the data
all.genes.2 <- rownames(seurat.obj.2)
seurat.obj.2 <- ScaleData(seurat.obj.2, features = all.genes.2)

# 3. Perform linear dimensionality reduction using PCA
seurat.obj.2 <- RunPCA(seurat.obj.2)

# visualising PCA results
seurat.obj.2[['pca']]
print(seurat.obj.2[['pca']], dims = 1:5, nfeatures = 5)
DimHeatmap(seurat.obj.2, dims=1:9, cells=1534, balanced=T)
head(Embeddings(seurat.obj.2, reduction = 'pca')[,1:5])
head(Loadings(seurat.obj.2, reduction = 'pca')[,1:5])

# visualising the loadings
VizDimLoadings(seurat.obj.2, dims = 1:2, reduction = "pca") 
VizDimLoadings(seurat.obj.2, dims = 3:4, reduction = "pca")                                 

# determine the dimensionality of the data by selecting the PCs which carry approx 70% variance of the data
ElbowPlot(seurat.obj.2)    

# Running UMAP

#################### INTEGRATING THE DATASETS ###############################

# 1.Select integration features
features <- SelectIntegrationFeatures(object.list = list(seurat.obj, seurat.obj.2))

# 2.Find integration anchors (CCA)
anchors <- FindIntegrationAnchors(object.list = list(seurat.obj,seurat.obj.2),
                                  anchor.features = features,
                                  dims = 1:9, k.filter = 9, k.score = 9)
anchors 

# 3.Integrate data
seurat.integrated <- IntegrateData(anchorset=anchors, k.weight = 10)
seurat.integrated
str(seurat.integrated)

# Scale data, run PCA 
all.genes.3 <- rownames(seurat.integrated)
seurat.integrated <- ScaleData(seurat.integrated, features = all.genes.3)
seurat.integrated <- RunPCA(seurat.integrated)
ElbowPlot(seurat.integrated)                                                     # we will again be selecting 15 PCs for neighborhood graph construction

########################### CLUSTERING ##########################

seurat.integrated <- FindNeighbors(seurat.integrated, dims=1:15)
seurat.integrated <- FindClusters(seurat.integrated, resolution = 0.013)

# 1.visualising the clusters using UMAP
seurat.integrated <- RunUMAP(seurat.integrated, dims=1:15)
DimPlot(seurat.integrated, reduction = 'umap', label = T)

# 2. Cluster annotation via identification of Marker Genes
# 2.1 find markers for every cluster compared to all remaining cells, report only the positive ones
markers <- FindAllMarkers(seurat.integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 1)      # keeping only those genes which have logfc >= 1
markers %>% group_by(cluster) %>% slice_max(n = 2, order_by = avg_log2FC)
markers
write.csv(markers, "marker_gene_4.csv")                                                                   # saving in csv file for easier access

# 2.2 manual annotation according to the marker genes
new.cluster.ids <- c("Tumor", "Stromal", "Myeloid", "Bcell/Tcell")
names(new.cluster.ids) <- levels(seurat.integrated)
seurat.integrated <- RenameIdents(seurat.integrated, new.cluster.ids)
DimPlot(seurat.integrated, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(seurat.integrated, features = top10$gene) + NoLegend()

#### Reference
# Jiang K, Dong M, Li C, Sheng J. Unraveling Heterogeneity of Tumor Cells and 
# Microenvironment and Its Clinical Implications for Triple Negative Breast Cancer. 
# Front Oncol. 2021 Mar 29;11:557477. 
# doi: 10.3389/fonc.2021.557477. PMID: 33854958; PMCID: PMC8040954.

