library(pheatmap)
library(reshape2)
library(viridis)
color_scheme <- colorRampPalette(c("white", "#582630"))(20)
factor_r2 <- read.csv('/data/scPAFA_paper_formal/revise/GSE164522/plot/20240722_factor_r2.csv',row.names = 1)
factor_r2 <- factor_r2[factor_r2$Factor %in% c('Factor4','Factor5','Factor6'),]

wide_data <- dcast(factor_r2, Factor ~ View, value.var = "R2")
row.names(wide_data) <- wide_data$Factor
wide_data$Factor <- NULL

pheatmap(wide_data)

factor_plot <- wide_data['Factor4',]
factor_plot  <- data.frame(t(factor_plot))
row <- as.numeric(factor_plot[,1])
# 获取排序后的索引
sorted_indices <- order(-row)[1:10]
pdf('/data/scPAFA_paper_formal/revise/GSE164522/plot/20240723_factor4_r2.pdf',width = 3,height = 6)
pheatmap(factor_plot[sorted_indices,,drop = FALSE],cluster_cols = F,cluster_rows = F,color = color_scheme)
dev.off()


factor_plot <- wide_data['Factor5',]
factor_plot  <- data.frame(t(factor_plot))
row <- as.numeric(factor_plot[,1])
# 获取排序后的索引
sorted_indices <- order(-row)[1:10]
pdf('/data/scPAFA_paper_formal/revise/GSE164522/plot/20240723_factor5_r2.pdf',width = 3,height = 6)
pheatmap(factor_plot[sorted_indices,,drop = FALSE],cluster_cols = F,cluster_rows = F,color = color_scheme)
dev.off()


factor_plot <- wide_data['Factor6',]
factor_plot  <- data.frame(t(factor_plot))
row <- as.numeric(factor_plot[,1])
# 获取排序后的索引
sorted_indices <- order(-row)[1:10]
pdf('/data/scPAFA_paper_formal/revise/GSE164522/plot/20240723_factor6_r2.pdf',width = 3,height = 6)
pheatmap(factor_plot[sorted_indices,,drop = FALSE],cluster_cols = F,cluster_rows = F,color = color_scheme)
dev.off()
