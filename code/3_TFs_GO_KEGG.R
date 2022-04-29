library(biomaRt)
library(limma)
library(Seurat)
library(org.Gg.eg.db)
library(topGO)


source("E:/Device_D/R_source/Seurat_v3_VizSource.R")
source("E:/Device_D/R_source/utils.R")

setwd("E:/Device_D/SNU_Chicken/DR/")

seurat_male <- readRDS("data/rds/Male_2.5-8_Deva.rds")
seurat_female <- readRDS("data/rds/Female_2.5-8_Deva.rds")


######GOBP#######
#load genes
tfs <- list(m2.5 = read.csv("data/GO_TFs/Up_Male_E2.5.csv", header = F)$V1,
            f2.5 = read.csv("data/GO_TFs/Up_Female_E2.5.csv", header = F)$V1,
            m68 = read.csv("data/GO_TFs/Up_Male_E6-E8.csv", header = F)$V1,
            f68 = read.csv("data/GO_TFs/Up_Female_E6-E8.csv", header = F)$V1)


#Male
backgroundGenes = rownames(seurat_male)
for(j in c("m2.5", "m68")){
  targetGenes = tfs[[j]]
  allGene = factor(as.integer(backgroundGenes %in% targetGenes))
  names(allGene) = backgroundGenes
  library(plyr)
  onts= c("BP")
  tab = as.list(onts)
  names(tab) = onts
  for (i in 1:length(onts)) {
    tgd = new("topGOdata", ontology = onts[i], allGenes=allGene, nodeSize=10,
              annot=annFUN.org, mapping="org.Gg.eg.db", ID = "symbol")
    resultTopGO.elim = runTest(tgd, algorithm="elim", statistic = "Fisher")
    tab[[i]] = GenTable(tgd, Fisher.elim = resultTopGO.elim,
                        orderBy="Fisher.classic", topNodes=300, numChar = 1000L)
  }
  
  topGOResults = rbind.fill(tab)
  topGOResults$Fisher.elim <- as.numeric(topGOResults$Fisher.elim)
  topGOResults <- topGOResults[topGOResults$Fisher.elim < 0.05, ]
  topGOResults[["'-log10p"]] <- -log10(topGOResults$Fisher.elim)
  
  #GenesInTerm
  allGO = genesInTerm(tgd)
  
  go_ids <- topGOResults$GO.ID
  go_genelist <- list()
  for(i in c(1:nrow(topGOResults))){
    topGOResults[i, "Genes"] <- paste(unlist(allGO[go_ids[i]])[unlist(allGO[go_ids[i]]) %in% targetGenes], collapse = "; ")
  }
  
  write.csv(topGOResults, file=paste("data/GO_TFs/GO_Results/GOBP_Male_Upreg_", j, ".csv", sep=""))
}

#Female
backgroundGenes = rownames(seurat_female)
for(j in c("f2.5", "f68")){
  targetGenes = tfs[[j]]
  allGene = factor(as.integer(backgroundGenes %in% targetGenes))
  names(allGene) = backgroundGenes
  library(plyr)
  onts= c("BP")
  tab = as.list(onts)
  names(tab) = onts
  for (i in 1:length(onts)) {
    tgd = new("topGOdata", ontology = onts[i], allGenes=allGene, nodeSize=10,
              annot=annFUN.org, mapping="org.Gg.eg.db", ID = "symbol")
    resultTopGO.elim = runTest(tgd, algorithm="elim", statistic = "Fisher")
    tab[[i]] = GenTable(tgd, Fisher.elim = resultTopGO.elim,
                        orderBy="Fisher.classic", topNodes=300, numChar = 1000L)
  }
  
  topGOResults = rbind.fill(tab)
  topGOResults$Fisher.elim <- as.numeric(topGOResults$Fisher.elim)
  topGOResults <- topGOResults[topGOResults$Fisher.elim < 0.05, ]
  topGOResults[["'-log10p"]] <- -log10(topGOResults$Fisher.elim)
  
  
  #GenesInTerm
  allGO = genesInTerm(tgd)
  
  go_ids <- topGOResults$GO.ID
  go_genelist <- list()
  for(i in c(1:nrow(topGOResults))){
    topGOResults[i, "Genes"] <- paste(unlist(allGO[go_ids[i]])[unlist(allGO[go_ids[i]]) %in% targetGenes], collapse = "; ")
  }
  
  write.csv(topGOResults, file=paste("data/GO_TFs/GO_Results/GOBP_Female_Upreg_", j, ".csv", sep=""))
}


#####KEGG#####
#load genes
tfs <- list(m2.5 = read.csv("data/GO_TFs/Up_Male_E2.5.csv", header = F)$V1,
            f2.5 = read.csv("data/GO_TFs/Up_Female_E2.5.csv", header = F)$V1,
            m68 = read.csv("data/GO_TFs/Up_Male_E6-E8.csv", header = F)$V1,
            f68 = read.csv("data/GO_TFs/Up_Female_E6-E8.csv", header = F)$V1)

ensembl <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset="ggallus_gene_ensembl", host="www.ensembl.org")
ensemblGenes <- getBM(attributes=c('ensembl_gene_id', 'external_gene_name',  'entrezgene_id', 'gene_biotype'), mart=ensembl)
ensemblGenes[ensemblGenes$external_gene_name == "", "external_gene_name"] <- ensemblGenes[ensemblGenes$external_gene_name == "", "ensembl_gene_id"]

tfs_entrez <- list()
for(i in names(tfs)){
   tfs_entrez[[i]] <- ensemblGenes[ensemblGenes$external_gene_name %in% tfs[[i]],]$entrezgene_id
}
#omitted: m2.5-5, m68-6,39,40, f68-5
tfs_entrez$m2.5[5] <- 396205 #NFYB

tfs_entrez$m68[6] <- 769693 #DMRT1
tfs_entrez$m68[39] <- 107054265 #IRX6
tfs_entrez$m68[40] <- 107054273 #IRX3

tfs_entrez$f68[5] <- 769693 #DMRT1

for(i in names(tfs_entrez)){
  tfs_entrez[[i]] <- as.character(tfs_entrez[[i]])
}


df_gse <- list()
for(i in names(tfs_entrez)){
  df <- limma::kegga(tfs_entrez[[i]], species = "Gg", FDR = 0.05)
  df <- df[df$P.DE != 1,]
  df <- df[order(df$P.DE), ]
  colnames(df) <- c("Pathway", "number_of_genes", "number_of_genes_from_query", "p_adj")
  df[["-log10p"]] <- -log10(df[["p_adj"]])
  
  paths <- substr(rownames(df), 6, 20L)
  .libPaths(Sys.getenv("R_LIBS_USER"))
  library(KEGGREST)
  
  gene_vec <- c()
  for(j in c(1:length(paths))){
    genes_path <- keggGet(paths[j])[[1]]$GENE[c(TRUE, FALSE)]
    genes_path <- genes_path[genes_path %in% tfs_entrez[[i]]]
    genes_symbol <- ensemblGenes[ensemblGenes$entrezgene_id %in% genes_path,]$external_gene_name
    
    gene_vec <- c(gene_vec, paste(genes_symbol, collapse = "; "))
  }

  df$Genes <- gene_vec
  write.csv(df, paste0("data/GO_TFs/GO_Results/KEGG/KEGG_", i, ".csv"))
}


paths <- substr(rownames(df), 6, 20L)
genes_path <- keggGet(paths[1])[[1]]$GENE[c(TRUE, FALSE)]