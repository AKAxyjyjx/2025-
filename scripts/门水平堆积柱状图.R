
# 导入数据并整理
bacteria <- read.delim("temperary_data/phylum_rel.txt", row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
# 每一行求和，并添加在最后一列
bacteria$sum <- rowSums(bacteria)
# 按照sum列排序
bacteria <- bacteria[order(bacteria$sum, decreasing = TRUE),]
# 导出表格
#write.table(bacteria, "C:状_fungus_phylum.txt", sep = '\t', quote = FALSE)
#这一部分是拖到excel里手操哦，整理导出的文件，如果分类多于10个的话可以把多的一部分加和在一起，命名为Others
# ============================================
# 9个最高丰度的 + 1个 "Others" = 10行
# ============================================
# 取分组平均值
samples <- unique(gsub("\\d+", "", colnames(bacteria))) # 提取样品组名（去掉所有数字）
bacteria_avg <- data.frame(row.names = rownames(bacteria))

for (sample in samples) {
  # 找出当前组的所有列
  sample_cols <- grep(paste0("^", sample, "\\d*$"), colnames(bacteria), value = TRUE)
  
  # 计算当前组的行平均值
  if(length(sample_cols) > 1) {
    bacteria_avg [[sample]] <- rowMeans(bacteria[, sample_cols, drop = FALSE])
  } else {
    bacteria_avg [[sample]] <- bacteria[, sample_cols] # 如果只有一列就直接取值
  }
}

# 2. 把行名变成一列
bacteria_avg <- cbind(Phylum = rownames(bacteria_avg), bacteria_avg)
rownames(bacteria_avg) <- NULL

# 3. 现在数据是
#   Phylum           EL1    EL2    EL3 ...
# 1 Proteobacteria  0.5    0.3    0.2
# 2 Firmicutes      0.3    0.4    0.3
# ...


# 5. 取前9行
top9 <- bacteria_avg[1:9, ]

others <- bacteria_avg[10:nrow(bacteria_avg), setdiff(colnames(bacteria_avg), c("Phylum", "sum")), drop = FALSE]
others_sum <- colSums(others)

others_row <- as.data.frame(t(c(Phylum = "Others", others_sum)))
colnames(others_row) <- c("Phylum", setdiff(colnames(bacteria_avg), c("Phylum", "sum")))


# 7. 合并（正好10行）

bacteria_top10 <- rbind(bacteria_avg[1:9, c("Phylum", setdiff(colnames(bacteria_avg), c("Phylum", "sum")))], others_row)
bacteria_top10 <- bacteria_top10[, !names(bacteria_top10) %in% "sum"]

# 导入数据并整理
#bacteria <- read.delim("C:积柱状_fungus_phylum.txt", row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
# 加一列行名，便于后续的长宽转换
#bacteria_df <- cbind(Phylum = row.names(bacteria_top10), bacteria_top10)
# 将bacteria转换为长格式
#bacteria_plot <- bacteria_df %>% gather(key = "Sample", value = "Abundance", -Phylum)
bacteria_plot <- bacteria_top10 %>% gather(key = "Sample", value = "Abundance", -Phylum)
# 开始绘图
# 自定义颜色
#install.packages("MetBrewer")
#library(MetBrewer) #还是一样，已经在main里加载过了，放在这里说明一下
mycol <- met.brewer("Hokusai1",10)
mycol
#install.packages("wesanderson")
#library(wesanderson)
#install.packages("ggpubr")
#library(ggpubr)

# 将 Abundance 转换为数值型
bacteria_plot$Abundance <- as.numeric(bacteria_plot$Abundance)





ggbarplot(bacteria_plot, x = "Sample", y="Abundance", color="black", fill="Phylum",
          legend="right", 
          font.main = c(14,"bold", "black"), font.x = c(12, "bold"), font.y=c(12,"bold")) + 
  theme_bw() +
  rotate_x_text() + 
  # 更改颜色
  scale_fill_manual(values = c(rev(mycol), "gray")) +
  #facet_grid(~ group, scales = "free_x", space='free') + # 分面设置
  labs(x = "Sample", y = "Relative Abundance") + 
  theme(text = element_text(family = "serif"),
    axis.text.x=element_text(angle = 45, hjust = 1, vjust = 1, size = 16, face = "bold", colour = "black"),  # 设置x轴刻度，倾斜45°
        axis.text.y = element_text(size = 16, face = "bold", colour = "black"),   # 设置y刻度
        axis.title.x = element_blank(),    # 隐藏x轴标题
        axis.title.y = element_text(size = 20, face = "bold"),   # 设置y轴标题
        legend.title = element_text(size = 16, face = "bold"),   # 设置图例标题
        legend.text = element_text(size = 12, face = "bold"),    # 设置图例文本字体
        panel.border = element_rect(colour = "black", fill = NA, size = 2),   # 设置面板边框
        panel.grid.major = element_line(colour = "gray90", size = 0.75, linetype = "solid"),   # 设置主要网格线
        panel.grid.minor = element_line(colour = "gray95", size = 0.5, linetype = "solid"),   # 设置次要网格线
  )




# 将图片以.pdf的格式导出
ggsave("output/水平堆积柱状_bacteria_phylum.pdf", width = 8, height = 5, units = "in")


