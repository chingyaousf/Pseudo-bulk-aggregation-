#clean up your environment
rm(list=ls())

# load library
library(Seurat)
library(SeuratData)
library(SeuratDisk)

#Load Seurat Object:
h5_friendly.h5seurat <- LoadH5Seurat("C:/Users/Administrator/Desktop/single_cell_analysis_aggregation/GSE162631_merge_tumor_h5.h5seurat") 

# Create the directory if it doesn't exist
output_folder <- "C:/Users/Administrator/Desktop/single_cell_analysis_aggregation/GSE162631_merge_tumor_h5.h5seurat_pseudo_bulk"
if (!dir.exists(output_folder)) {
  dir.create(output_folder, recursive = TRUE)
}

# Save the metadata to a text file in the specified path with column names
write.table(h5_friendly.h5seurat@meta.data, file.path(output_folder, "GSE162631_merge_tumor_h5.h5seurat_meta.txt"), sep="\t", row.names = TRUE, col.names = TRUE)

# Display Head of Metadata:
head(h5_friendly.h5seurat@meta.data) 

# Aggregate Expression Data
agg = AggregateExpression(h5_friendly.h5seurat, return.seurat = T, group.by = c('orig.ident', 'seurat_clusters_gabby_annotation_revised_v4'), normalization.method = "LogNormalize",scale.factor = 10000)

head(agg)

# generating pseudo bulk level meta file
cnames = colnames(agg@assays$RNA$data)
cnames

meta1 = do.call(rbind, strsplit(cnames, "_"))
meta1

meta_text = cbind(cnames, meta1)
meta_text

# Specify File Paths:
file <- file.path(output_folder, "GSE162631_merge_tumor_expression_matrix.txt")
meta_file <- file.path(output_folder, "GSE162631_merge_tumor_information.txt")

# Write Expression Matrix to File:
write.table(2^ agg@assays$RNA$data, file = file, sep="\t", row.names = TRUE,col.names=NA)

# Write Meta Information to File:
write.table(meta_text, file = meta_file, sep="\t", row.names=F) 

