############################################
# GO enrichment (SYMBOL version)
############################################

library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)
library(ggplot2)

# 读取 DESeq2 或 edgeR 结果
res_df <- read.csv("../5clustering/edgeR_DE_results.csv", row.names = 1)

# 去NA
res_df <- res_df[!is.na(res_df$FDR), ]
# 筛选 DEG
deg <- res_df[
  res_df$PValue < 0.05 &
  abs(res_df$logFC) > 1,
]

genes <- rownames(deg)

############################################
# GO enrichment（🔥核心）
############################################

ego <- enrichGO(
  gene = genes,
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  keyType = "SYMBOL"
)

############################################
# 可视化
############################################

pdf("1GO_dotplot.2.pdf")
dotplot(ego, showCategory = 20)
dev.off()

pdf("1GO_barplot.2.pdf")
barplot(ego, showCategory = 20)
dev.off()

############################################
# 网络图（加分🔥）
############################################

fold_change <- deg$logFC
names(fold_change) <- rownames(deg)

pdf("2GO_cnetplot.pdf")
cnetplot(ego, foldChange = fold_change, showCategory = 5)
dev.off()

pdf("2GO_heatplot.pdf")
heatplot(ego, foldChange = fold_change, showCategory = 10)
dev.off()

############################################
# emapplot（高级🔥）
############################################

library(enrichplot)
sim <- pairwise_termsim(ego)

pdf("GO_emapplot.pdf")
emapplot(sim, showCategory = 5)
dev.off()

############################################
# 保存结果
############################################

write.csv(as.data.frame(ego), "GO_results.csv")