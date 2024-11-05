library(pheatmap)
library(RColorBrewer)

scpa_out_500 <- read.csv('/data/scPAFA_paper_formal/scPA_and_scPAFA/CRC_aviv_dataset/20240725_qval_matrix_CRC_500cells.csv',header = T,row.names = 1)
# 获取每列最大5个值的行名
row_names_list <- vector()
for (col in colnames(scpa_out_500)) {
  top5 <- order(scpa_out_500[[col]], decreasing = TRUE)[1:3]
  row_names_list <- c(row_names_list, row.names(scpa_out_500)[top5])
}
# 获取这些行名的并集
unique_rows <- unique(row_names_list)
unique_rows_scPA_500 <- unique_rows
scpa_out_500_plot <- scpa_out_500[unique_rows, ]
scpa_out_500_plot <- scpa_out_500_plot[order(rowSums(scpa_out_500_plot),decreasing = T), ]

colnames(scpa_out_500_plot)

scpa_out_500_plot <- scpa_out_500_plot[,c("Tumor_cE01_.StemTA.like._qval","Tumor_cE02_.StemTA.likeImmature_Goblet._qval","Tumor_cE03_.StemTA.like_prolif._qval","Tumor_cE05_.Enterocyte_2._qval","Tumor_cE06_.Immature_Goblet._qval")]

breaks <- seq(0, max(scpa_out_500_plot), length.out = 101)
pdf('/data/scPAFA_paper_formal/scPA_and_scPAFA/final_heatmap/20240806_CRC_scPA_500.pdf',width = 6,height = 6)
pheatmap(scpa_out_500_plot,breaks = breaks,cluster_cols = FALSE,cluster_rows = F,color = colorRampPalette(brewer.pal(9, "BuPu"))(100),border_color = NA)
dev.off()

scpa_out_2000 <- read.csv('/data/scPAFA_paper_formal/scPA_and_scPAFA/CRC_aviv_dataset/20240725_qval_matrix_CRC_2000cells.csv',header = T,row.names = 1)
# 获取每列最大5个值的行名
row_names_list <- vector()
for (col in colnames(scpa_out_2000)) {
  top5 <- order(scpa_out_2000[[col]], decreasing = TRUE)[1:3]
  row_names_list <- c(row_names_list, row.names(scpa_out_2000)[top5])
}
# 获取这些行名的并集
unique_rows <- unique(row_names_list)
unique_rows_scPA_2000 <- unique_rows
scpa_out_2000_plot <- scpa_out_2000[unique_rows, ]
scpa_out_2000_plot[scpa_out_2000_plot == Inf] = max(max(scpa_out_2000_plot[scpa_out_2000_plot != Inf]))
scpa_out_2000_plot <- scpa_out_2000_plot[order(rowSums(scpa_out_2000_plot),decreasing = T), ]

scpa_out_2000_plot <- scpa_out_2000_plot[,c("Tumor_cE01_.StemTA.like._qval","Tumor_cE02_.StemTA.likeImmature_Goblet._qval","Tumor_cE03_.StemTA.like_prolif._qval","Tumor_cE05_.Enterocyte_2._qval","Tumor_cE06_.Immature_Goblet._qval")]

breaks <- seq(0, max(scpa_out_2000_plot), length.out = 101)
pdf('/data/scPAFA_paper_formal/scPA_and_scPAFA/final_heatmap/20240806_CRC_scPA_2000.pdf',width = 6,height = 6)
pheatmap(scpa_out_2000_plot,breaks = breaks,cluster_cols = FALSE,cluster_rows = F,color = colorRampPalette(brewer.pal(9, "BuPu"))(100),border_color = NA)
dev.off()

#draw_scPAFA_pathway_heatmap
scPAFA_weight <- read.csv('/data/scPAFA_paper_formal/CRC/CRC_data_reverse/20240725_scPAFA_celltype_pathway_weight.csv',row.names = 1)
row_names_list <- vector()
for (col in colnames(scPAFA_weight)) {
  top5 <- order(scPAFA_weight[[col]], decreasing = TRUE)[1:3]
  row_names_list <- c(row_names_list, row.names(scPAFA_weight)[top5])
}

unique_rows <- unique(row_names_list)
unique_rows_scPAFA <- unique_rows

scPAFA_weight <- scPAFA_weight[unique_rows,]
scPAFA_weight[is.na(scPAFA_weight)] = 0
scPAFA_weight <- scPAFA_weight[order(rowSums(scPAFA_weight),decreasing = T), ]
breaks <- seq(0, max(scPAFA_weight), length.out = 101)
pdf('/data/scPAFA_paper_formal/scPA_and_scPAFA/final_heatmap/20240806_CRC_scPAFA.pdf',width = 6.5,height = 5)
pheatmap(scPAFA_weight,breaks = breaks,cluster_cols = FALSE,cluster_rows = F,color = colorRampPalette(brewer.pal(9, "BuPu"))(100),border_color = NA)
dev.off()

intersect(unique_rows_scPA_500 ,unique_rows_scPA_2000)
intersect(unique_rows_scPA_500 ,unique_rows_scPAFA)
intersect(unique_rows_scPA_2000 ,unique_rows_scPAFA)

