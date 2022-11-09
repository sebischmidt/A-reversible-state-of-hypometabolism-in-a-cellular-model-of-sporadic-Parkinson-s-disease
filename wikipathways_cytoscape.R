######Load packages######
library(biomaRt)
library(org.Hs.eg.db)
library(clusterProfiler)
library(DOSE)
library(devtools)
library(enrichplot)
library(devEMF)
library(extrafont)
library(rWikiPathways)
library(ggplot2)
library(DOSE, dplyr, stringr)
library(RCy3)
library(GSEABase)
library(tidyr)
library(RColorBrewer)
library(readxl)
cytoscapePing()  #this will tell you if you're able to successfully connect to Cytoscape or not
installApp('WikiPathways') 
installApp('CyTargetLinker') 
installApp('stringApp') 

######Load data######
setwd('')


DEPs <- read.delim("proteome_diff_res_DESeq2.csv", sep = ",")
DEGs <- read_xlsx("41467_2022_32229_MOESM8_ESM_bulk-like.xlsx")
metabolites <- read.delim('non-targeted-metabolomics_diff_res.csv.csv', sep = ",")

######Prepare data for merge######

genes <- as.character(DEPs[,2])
mart = useDataset("hsapiens_gene_ensembl", useEnsembl(biomart="ensembl", version=98))
G_list <- getBM(filters= "hgnc_symbol", attributes= c("ensembl_gene_id", "hgnc_symbol", "entrezgene_id"),values=genes,mart=mart)
DEPs$ENTREZID = G_list[match(genes, G_list$hgnc_symbol),]$entrezgene_id
DEPs$ensemblid = G_list[match(genes, G_list$hgnc_symbol),]$ensembl_gene_id
colnames(DEPs)[2] <- "SYMBOL"
colnames(DEPs)[4] <- "log2fc"

genes <- as.character(DEGs[,2])
G_list <- getBM(filters= "hgnc_symbol", attributes= c("ensembl_gene_id", "hgnc_symbol", "entrezgene_id"),values=genes,mart=mart)
DEGs$ENTREZID = G_list[match(genes, G_list$hgnc_symbol),]$entrezgene_id
DEGs$ensemblid = G_list[match(genes, G_list$hgnc_symbol),]$ensembl_gene_id
colnames(DEGs)[2] <- "SYMBOL"
colnames(DEPs)[4] <- "log2fc"
DEGs$log2fc <- log2(DEGs$log2fc)

metabolites$fc <- log2(metabolites$fc)
colnames(metabolites)[16] <- "log2fc"

######Plot figures in Cytoscape - metabolites + DEPs######
##TCA
RCy3::commandsRun('wikipathways import-as-pathway id=WP78') 
loadTableData(DEPs, data.key.column = "ENTREZID", table.key.column = "XrefId")
loadTableData(metabolites, data.key.column = "HMDB", table.key.column = "XrefId")
data.values = c(-6,-0.3,-0.2,-0.1,-0.05,0,0.05,0.1,0.2,0.3,6)
node.colors <- c(rev(brewer.pal(length(data.values), "RdBu")))
setNodeColorMapping("log2fc", data.values, node.colors, default.color = "#FFFFFF", style.name = "WikiPathways")
exportImage(filename="proteome_metabolites_tca.pdf", type='PDF')

##OXPHOS
RCy3::commandsRun('wikipathways import-as-pathway id=WP111') 
loadTableData(DEPs, data.key.column = "ENTREZID", table.key.column = "XrefId")
loadTableData(DEPs, data.key.column = "ensemblid", table.key.column = "XrefId")
loadTableData(metabolites, data.key.column = "PUBCHEM", table.key.column = "XrefId")
loadTableData(metabolites, data.key.column = "HMDB", table.key.column = "XrefId")
loadTableData(metabolites, data.key.column = "KEGG", table.key.column = "XrefId")
data.values = c(-6,-0.3,-0.2,-0.1,-0.05,0,0.05,0.1,0.2,0.3,6)
node.colors <- c(rev(brewer.pal(length(data.values), "RdBu")))
setNodeColorMapping("log2fc", data.values, node.colors, default.color = "#FFFFFF", style.name = "WikiPathways")
exportImage(filename="proteome_metabolites_oxphos.pdf", type='PDF')

######Plot figures in Cytoscape - metabolites + DEGs######
##TCA
RCy3::commandsRun('wikipathways import-as-pathway id=WP78') 
loadTableData(DEGs, data.key.column = "ENTREZID", table.key.column = "XrefId")
loadTableData(metabolites, data.key.column = "HMDB", table.key.column = "XrefId")
data.values = c(-6,-0.3,-0.2,-0.1,-0.05,0,0.05,0.1,0.2,0.3,6)
node.colors <- c(rev(brewer.pal(length(data.values), "RdBu")))
setNodeColorMapping("log2fc", data.values, node.colors, default.color = "#FFFFFF", style.name = "WikiPathways")
exportImage(filename="proteome_metabolites_tca.pdf", type='PDF')

##OXPHOS
RCy3::commandsRun('wikipathways import-as-pathway id=WP111') 
loadTableData(DEGs, data.key.column = "ENTREZID", table.key.column = "XrefId")
loadTableData(DEGs, data.key.column = "ensemblid", table.key.column = "XrefId")
loadTableData(metabolites, data.key.column = "PUBCHEM", table.key.column = "XrefId")
loadTableData(metabolites, data.key.column = "HMDB", table.key.column = "XrefId")
loadTableData(metabolites, data.key.column = "KEGG", table.key.column = "XrefId")
data.values = c(-6,-0.3,-0.2,-0.1,-0.05,0,0.05,0.1,0.2,0.3,6)
node.colors <- c(rev(brewer.pal(length(data.values), "RdBu")))
setNodeColorMapping("log2fc", data.values, node.colors, default.color = "#FFFFFF", style.name = "WikiPathways")
exportImage(filename="proteome_metabolites_oxphos.pdf", type='PDF')
