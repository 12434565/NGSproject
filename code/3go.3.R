install.packages("BiocManager")
BiocManager::install("clusterProfiler")
BiocManager::install("org.Hs.eg.db")

############################################
# GO + KEGG Enrichment Analysis
############################################

# 📦 加载包
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
setwd("/wellstein/user_homes/il1380/project/6go2/")
############################################
# 1. 读取 DESeq2 结果
############################################

res_df <- read.csv("../5clustering/edgeR_DE_results.csv", row.names = 1)

# 去掉NA
res_df <- res_df[!is.na(res_df$FDR), ]

############################################
# 2. 筛选显著基因（和你volcano一致）
############################################

threshold_p <- 0.01
threshold_fc <- 2   # log2FC >1 = FC >2

deg <- res_df[
  res_df$FDR < threshold_p & abs(res_df$logFC) > threshold_fc,
]

cat("Number of DEGs:", nrow(deg), "\n")

############################################
# 3. 提取 gene symbol
############################################

genes <- rownames(deg)

############################################
# 4. 转换为 ENTREZID（关键🔥）
############################################

gene_df <- bitr(
  genes,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

# 去重复
gene_df <- gene_df[!duplicated(gene_df$ENTREZID), ]

############################################
# 5. GO enrichment（BP）
############################################

ego <- enrichGO(
  gene = gene_df$ENTREZID,
  OrgDb = org.Hs.eg.db,
  ont = "BP",   # Biological Process
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)
ego_mf <- enrichGO(
  gene = gene_df$ENTREZID,
  OrgDb = org.Hs.eg.db,
  ont = "MF",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)
ego_cc <- enrichGO(
  gene = gene_df$ENTREZID,
  OrgDb = org.Hs.eg.db,
  ont = "CC", 
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)

############################################
# 6. GO 可视化
############################################

pdf("GO_dotplot.pdf", width=8, height=8)
dotplot(ego, showCategory=20)
dev.off()

pdf("GO_barplot.pdf", width=8, height=8)
barplot(ego, showCategory=20)
dev.off()

pdf("GO_MF_dotplot.pdf", width=8, height=8)
print(
  dotplot(ego_mf, showCategory=20) +
    ggtitle("GO Molecular Function")
)
dev.off()
pdf("GO_CC_dotplot.pdf", width=8, height=8)
print(
  dotplot(ego_cc, showCategory=20) +
    ggtitle("GO Cellular Component")
)
dev.off()
write.csv(as.data.frame(ego), "GO_BP_results.csv")
write.csv(as.data.frame(ego_mf), "GO_MF_results.csv")
write.csv(as.data.frame(ego_cc), "GO_CC_results.csv")
############################################
# 7. KEGG enrichment（加分🔥）
############################################
deg_sig <- res_df[
  res_df$FDR < 0.05 &
  abs(res_df$logFC) > 1,
]
gene_df <- bitr(
  rownames(deg_sig),
  fromType="SYMBOL",
  toType="ENTREZID",
  OrgDb=org.Hs.eg.db
)
# library(clusterProfiler)
library(KEGG.db)
# kegg_ids <- bitr_kegg(
#   gene_df$ENTREZID,
#   fromType = "ncbi-geneid",
#   toType = "kegg",
#   organism = "hsa"
# )
# kk <- enrichKEGG(
#   gene = gene_df$ENTREZID,
#   organism = "hsa",
#   use_internal_data = TRUE
# )
# kegg_gene <- paste0("hsa:", gene_df$ENTREZID) 
# kk <- enrichKEGG(
#   gene      = kegg_gene,
#   organism  = "hsa",
#   keyType   = "ncbi-geneid"   # 🔥 必加
# )

# pdf("KEGG_dotplot.pdf", width=8, height=6)
# dotplot(kk, showCategory=20)
# dev.off()

library(msigdbr)
library(dplyr)
# > unique(msigdbr()$gs_subcollection)
#  [1] "MIR:MIR_LEGACY"  "TFT:TFT_LEGACY"  "CGP"             "TFT:GTRD"       
#  [5] ""                "VAX"             "CP:BIOCARTA"     "CGN"            
#  [9] "3CA"             "GO:BP"           "GO:CC"           "IMMUNESIGDB"    
# [13] "GO:MF"           "HPO"             "CP:KEGG_LEGACY"  "CP:KEGG_MEDICUS"
# [17] "MIR:MIRDB"       "CM"              "CP:PID"          "CP:REACTOME"    
# [21] "CP"              "CP:WIKIPATHWAYS"
m_df <- msigdbr(
  species = "Homo sapiens",
  collection = "C2",
  subcollection = "CP:KEGG_LEGACY"
)

kegg_list <- split(
  x = as.character(m_df$ncbi_gene),
  f = m_df$gs_name
)
############################################
# 8. 保存结果表
############################################

write.csv(as.data.frame(ego), "GO_results.csv")
write.csv(as.data.frame(kk), "KEGG_results.csv")

deg_sig <- res_df[
  res_df$FDR < 0.05 &
  abs(res_df$logFC) > 1,
]
gene_map <- bitr(
  rownames(deg_sig),
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)
deg <- as.character(unique(gene_map$ENTREZID))
############################################
# 4️⃣ 背景基因（很重要）
############################################
bg_map <- bitr(
  rownames(res_df),
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)
bg <- as.character(unique(bg_map$ENTREZID))
############################################
# 5️⃣ KEGG enrichment（Fisher test）
############################################
kegg_res <- data.frame()
for (pw in names(kegg_list)) {
  geneset <- kegg_list[[pw]]
  overlap <- length(intersect(deg, geneset))
  if (overlap < 3) next
  a <- overlap
  b <- length(deg) - a
  c <- length(geneset) - a
  d <- length(bg) - a - b - c
  if (d < 0) next
  mat <- matrix(c(a, b, c, d), nrow = 2)
  pval <- fisher.test(mat)$p.value
  kegg_res <- rbind(
    kegg_res,
    data.frame(
      pathway = pw,
      overlap = a,
      pvalue = pval
    )
  )
}
############################################
# 6️⃣ 多重校正 + 排序
############################################
kegg_res$padj <- p.adjust(kegg_res$pvalue, method = "BH")
kegg_res <- kegg_res[order(kegg_res$padj), ]
############################################
# 7️⃣ 保存结果
############################################
write.csv(kegg_res, "KEGG_results.csv", row.names = FALSE)
############################################
# 8️⃣ 画 KEGG dotplot
############################################
top <- head(kegg_res, 20)
pdf("KEGG_dotplot.pdf", width=10, height=8)
ggplot(top, aes(x=reorder(pathway, overlap), y=overlap)) +
  geom_point(aes(size = -log10(padj), color = padj)) +
  coord_flip() +
  theme_minimal() +
  ggtitle("KEGG Pathway Enrichment")
dev.off()
############################################
# DONE

############################################

print("KEGG enrichment complete!")

###
edge_df <- data.frame()

for (pw in names(kegg_list)) {

  geneset <- kegg_list[[pw]]
  overlap_genes <- intersect(deg, geneset)

  if (length(overlap_genes) < 3) next

  temp <- data.frame(
    pathway = pw,
    gene = overlap_genes
  )

  edge_df <- rbind(edge_df, temp)
}

write.csv(edge_df, "cytoscape_edges.csv", row.names = FALSE, quote = FALSE)
# gene 节点
node_df <- merge(
  deg_sig,
  gene_map,
  by.x = "row.names",   # SYMBOL
  by.y = "SYMBOL"
)
colnames(node_df)[1] <- "SYMBOL"
gene_nodes <- data.frame(
  id = as.character(node_df$ENTREZID),
  label = node_df$SYMBOL,
  log2FC = node_df$logFC,
  padj = node_df$FDR,
  type = "gene"
)
gene_nodes <- gene_nodes[

  !is.na(gene_nodes$id) & gene_nodes$id != "",

]
# pathway 节点
pathway_nodes <- data.frame(
  id = unique(edge_df$pathway),
  label = unique(edge_df$pathway),
  log2FC = NA,
  padj = NA,
  type = "pathway"
)

# 合并
node_all <- rbind(gene_nodes, pathway_nodes)

write.csv(node_all, "cytoscape_nodes.csv",
          row.names = FALSE,
          quote = FALSE)