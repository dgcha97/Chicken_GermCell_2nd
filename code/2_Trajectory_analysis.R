library(Seurat)
library(leidenbase)
library(monocle3)

dr_male  <- readRDS("Male_2.5-8.rds")
dr_female  <- readRDS("Female_2.5-8.rds")


####Male####
#create cds object for monocle3
cell_metadata = dr_male@meta.data
gene_metadata = data.frame(gene_short_name = rownames(dr_male), row.names = rownames(dr_male))

cds <- new_cell_data_set(dr_male@assays$RNA@data,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_metadata)

cds <- preprocess_cds(cds, "PCA", num_dim = 7, norm_method = "none", use_genes = VariableFeatures(dr_male))
monocle3::plot_pc_variance_explained(cds)

cds <- reduce_dimension(cds, reduction_method = "UMAP")
plot_cells(cds, color_cells_by = "seurat_clusters", cell_size = 1, group_label_size = 5)
plot_cells(cds, color_cells_by = "condition", cell_size = 1, group_label_size = 5)

cds <- cluster_cells(cds, resolution = 0.05)

cds <- learn_graph(cds)

#set cell with maximum Pou5f3 expr as earliest principal node
get_earliest_principal_node <- function(cds, max.gene="Pou5f3"){
  cell_ids <- which.max(as.vector(cds@assays@data$counts[max.gene,]))
  
  closest_vertex <-
    cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
    igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
                                                              (which.max(table(closest_vertex[cell_ids,]))))]
  
  root_pr_nodes
}

cds <- order_cells(cds)

plot_cells(cds,
           color_cells_by = "pseudotime",
           cell_size = 1,
           label_groups_by_cluster = FALSE,
           label_leaves = FALSE,
           label_branch_points = FALSE)

cds_male <- cds


####Female####
#create cds object for monocle3
cell_metadata = dr_female@meta.data
gene_metadata = data.frame(gene_short_name = rownames(dr_female), row.names = rownames(dr_female))

cds <- new_cell_data_set(dr_female@assays$RNA@data,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_metadata)

cds <- preprocess_cds(cds, "PCA", num_dim = 7, norm_method = "none", use_genes = VariableFeatures(dr_female))
monocle3::plot_pc_variance_explained(cds)

cds <- reduce_dimension(cds, reduction_method = "UMAP")
plot_cells(cds, color_cells_by = "seurat_clusters", cell_size = 1, group_label_size = 5)
plot_cells(cds, color_cells_by = "condition", cell_size = 1, group_label_size = 5)

cds <- cluster_cells(cds, resolution = 0.01)
cds <- learn_graph(cds)


#set cell with maximum Pou5f3 expr as earliest principal node
get_earliest_principal_node <- function(cds, max.gene="Pou5f3"){
  cell_ids <- which.max(as.vector(cds@assays@data$counts[max.gene, cds@colData$condition %in% c("E6_Female", "E8_Female")]))
  
  closest_vertex <-
    cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
    igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
                                                              (which.max(table(closest_vertex[cell_ids,]))))]
  
  root_pr_nodes
}

cds <- order_cells(cds,
                   root_pr_nodes=get_earliest_principal_node(cds)
)

plot_cells(cds,
           color_cells_by = "pseudotime",
           cell_size = 1,
           label_groups_by_cluster = FALSE,
           label_leaves = FALSE,
           label_branch_points = FALSE)

cds_female <- cds



#add coordinate of monocle-calculated UMAP to Seurat object
umap_female <- cds_female@int_colData$reducedDims$UMAP
umap_female <- umap_female[colnames(dr_female),]
colnames(umap_female) <- c("1", "2")
dr_female[["MonoUMAP"]] <- CreateDimReducObject(embeddings = as.matrix(umap_female),
                                                               key = "MonoUMAP_", assay = DefaultAssay(dr_female))


umap_male <- cds_male@int_colData$reducedDims$UMAP
umap_male <- umap_male[colnames(dr_male),]
colnames(umap_male) <- c("1", "2")
dr_male[["MonoUMAP"]] <- CreateDimReducObject(embeddings = as.matrix(umap_male),
                                                             key = "MonoUMAP_", assay = DefaultAssay(dr_male))


dr_male$pseudotime <- pseudotime(cds_male)
dr_female$pseudotime <- pseudotime(cds_female)

dr_male$pseudotime[dr_male$pseudotime == Inf] <- NA


#save objects
saveRDS(cds_male, "Monocle3_Male_E2.5-E8.rds")
saveRDS(cds_female, "Monocle3_Female_E2.5-E8.rds")
