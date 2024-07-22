install.packages("readr")
install.packages("gdata")
install.packages("tidyr")
install.packages("devtools")
install_bitbucket("ibi_group/disgenet2r")

BiocManager::install("clusterProfiler")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("DESeq2")
BiocManager::install("Biostrings")
BiocManager::install("EnhancedVolcano")
BiocManager::install("ReactomePA")
BiocManager::install("reactome.db")
BiocManager::install("disgenet2r ")




library(disgenet2r )
library(readr)
library(DESeq2)
library(dplyr)
library(tidyr)
library(EnhancedVolcano)
library(pheatmap)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(DOSE)
library(gdata)
library(ReactomePA)
library(ggplot2)
library(reactome.db)
library(devtools)






#reading the samples names into a varible
getwd()
setwd("E:/Courses/RNA-seq/Bakry/normal")
list_of_normal_samples<-file.path(list.files())[1:5]

#reading the first file of normal sample in a varible 
#to us it as a counter for the other files
counter_norma_sample<-read.delim(list_of_normal_samples[1] ,sep ="\t",header = T,skip = 1,row.names = 1)

#preprocessing the file to remove the first rows that is side info not count data
counter_norma_sample <- counter_norma_sample[-(1:4), ]

#save it for annot later
gene_symbols <- counter_norma_sample %>%
 dplyr:: select(c(gene_name,gene_type))

#preprocess
counter_norma_sample <- counter_norma_sample %>%
  dplyr::select(-c(gene_name,gene_type)) %>%
  dplyr:: select(unstranded)#chosed the unstranded as the used platform is illumina

#reading & preprocess solid normal samples 
for (sample_index in 2:5) {
  combine_with_counter<-read.delim(list_of_normal_samples[sample_index]  ,sep ="\t",header = T,skip = 1,row.names = 1)
  combine_with_counter <- combine_with_counter[-(1:4), ]
  combine_with_counter <- combine_with_counter %>%
    #delete coulmns
    dplyr:: select(-c(gene_name,gene_type)) %>%
    dplyr::select(unstranded)
  counter_norma_sample <- cbind(counter_norma_sample, combine_with_counter)
}

#setting path & reading data
setwd("E:/Courses/RNA-seq/Bakry/tumor")
list_of_tumor_samples<-file.path(list.files())[1:5]
counter_tumor_sample<-read.delim(list_of_tumor_samples[1] ,sep ="\t",header = T,skip = 1,row.names = 1)
counter_tumor_sample <- counter_tumor_sample[-(1:4), ]

#preprocess
counter_tumor_sample <- counter_tumor_sample %>%
  dplyr::select(-c(gene_name,gene_type)) %>%
  dplyr::select(unstranded)

#reading the rest of tumor samples and combining them toghether 
for (sample_index in 2:5) {
  combine_with_counter<-read.delim(list_of_tumor_samples[sample_index] ,sep ="\t",header = T,skip = 1,row.names = 1)
  combine_with_counter <- combine_with_counter[-(1:4), ]
  combine_with_counter <- combine_with_counter %>%
    #delete coulmns
    dplyr::select(-c(gene_name,gene_type)) %>%
    dplyr::select(unstranded)
  counter_tumor_sample <- cbind(counter_tumor_sample, combine_with_counter)
}

# tumor and normal combined 
pre_deseq_data<-cbind(counter_tumor_sample,counter_norma_sample)

# Prepare metadata
metadata <- data.frame(condition = rep(c("tumor", "solid_normal"), each = 5))
metadata$condition <- factor(metadata$condition)  # Convert to factor for DESeq2

# Rename columns to match metadata
list_sample <- c(list_of_normal_samples, list_of_tumor_samples)

#making each coulmn unique
for (i in 1:10) {colnames(pre_deseq_data)[i] <- paste0(list_sample[i])}

# Ensure the row names of metadata match the column names of your data
rownames(metadata) <- colnames(pre_deseq_data)

#preperaing the object to be used in deseq
dds <- DESeqDataSetFromMatrix(countData = pre_deseq_data, colData = metadata, design = ~  condition)

# Normalization , Handling variance
dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)

#run deseq2
dds <- DESeq(dds)

# rlog transformation 
rld <- rlog(dds)

#dds results
res <- results(dds,contrast = c("condition","tumor","solid_normal"))

# LFC shrink 
resNorm <- lfcShrink(dds = dds, res = res, type = "normal", coef = 2)

# Plotting PCA plot ,the size factors,dispersion estimates,
#MA-Plot using the results from the DESeq2
pdf("./outputs/DEseq following Plots")
plotPCA(rld)+ ggtitle("PCA plot of rlog transformation")
barplot(sizeFactors(dds), main="Size Factors", ylab="Size Factor", xlab="Sample Index", col='blue')
plotDispEsts(dds, main="Dispersion Estimates")
plotMA(res, main="MA-Plot before lfsh", ylim=c(-5, 5))
plotMA(resNorm ,main="MA-Plot after shrink", ylim=c(-5, 5))
dev.off()

#clear workplace 
keep(gene_symbols,dds,res,resNorm,rld,metadata,sure = T)
gc()

#adding the ENGS as a coulmn
resdf<-as.data.frame(resNorm)# 'ENSG' are row names, make them a column
resdf$ENSG <- row.names(resdf)
gene_symbols$ENSG <- row.names(gene_symbols)

# Merging the two data frames on the 'ENSG' column
combined_data <- left_join(resdf, gene_symbols, by = "ENSG")

#volcanoplot
volcano<-EnhancedVolcano(combined_data, lab = combined_data$gene_name, pCutoff = 0.05,
               FCcutoff = 2,
               x = "log2FoldChange", y = "padj")

# Extract the normalized counts for pheatmap
normalized_counts <- counts(dds, normalized = TRUE)
top_genes <- res[order(res$pvalue),][1:50,] 
top_genes_matrix <- normalized_counts[rownames(top_genes), ]

# use the rlog for better visualization 
top_genes_matrix_transformed <- assay(rld)[rownames(top_genes), ]

# Generate the heatmap for the top genes 
pheat_map<-pheatmap(top_genes_matrix_transformed,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_rownames = TRUE,
         show_colnames = FALSE, # Turn off if you have too many samples
         color = colorRampPalette(c("blue", "white", "red"))(255), # Adjust the color palette if needed
         annotation_col = metadata # This should be your metadata data frame
)

# Save the heatmap, volcanoplot to a file
ggsave(plot = pheat_map,"heatmap.pdf", device = "pdf", path = getwd(), width = 10, height = 8) # Adjust dimensions as needed
ggsave(plot = volcano,"volcanoplot.pdf", device = "pdf", path = getwd(), width = 10, height = 8) # Adjust dimensions as needed

#extract the DE genes with the gene symbols and save it
DE_genes <- combined_data %>%
  filter(padj < 0.05& abs(log2FoldChange)>2) 
  #dplyr::select(gene_name)
write.csv(as.data.frame(DE_genes), "DE_genes.csv")

#get gene symbols
significant_genes<-DE_genes$gene_name

#Preprocess for GO
entrez_ids <- bitr(significant_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
entrez_ids <- na.omit(entrez_ids) #remove NA

#GO
go_results <- enrichGO(gene          = entrez_ids$ENTREZID,
                       OrgDb         = org.Hs.eg.db,
                       keyType       = 'ENTREZID',
                       ont           = "BP", # Biological Process
                       pAdjustMethod = "BH", # Benjamini-Hochberg adjustment
                       pvalueCutoff  = 0.05,
                       qvalueCutoff  = 0.2)

#KEGG Pathway Enrichment Analysis
kegg_results <- enrichKEGG(gene = entrez_ids$ENTREZID,
                           organism = 'hsa',
                           keyType = 'kegg',
                           pvalueCutoff = 0.05,
                           qvalueCutoff = 0.2)

#disgenet
disgenet_results <- enrichDGN(entrez_ids$ENTREZID, pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.2)

#saving the results
write.csv(as.data.frame(kegg_results@result), "kegg_results.csv")
write.csv(as.data.frame(go_results@result), "go_results.csv")
write.csv(as.data.frame(disgenet_results@result), "disgenet_results.csv")

#Enrichment_Plots.pdf
pdf("Enrichment_Plots.pdf")
dotplot(go_results) + ggtitle("GO Enrichment Results") + theme_minimal()
dotplot(kegg_results) + ggtitle("KEGG Pathway Results") + theme_minimal()
dotplot(disgenet_results) + ggtitle("DisGeNET Enrichment Results") + theme_minimal()
barplot(disgenet_results, showCategory = 20) + ggtitle("Top 20 Enriched DisGeNET Terms") + theme_minimal()
barplot(go_results, showCategory=15) + ggplot2::theme_minimal() + ggtitle("Top 20 GO ")
dotplot(kegg_results, showCategory=20) + ggplot2::theme_minimal()
dev.off()