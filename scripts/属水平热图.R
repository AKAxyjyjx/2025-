

# 导入数据，属相对丰度3不同分类水平物种数量统计/genus_rel
asv <- read.table("temperary_data/genus_rel.txt", header = T, row.names = 1, sep = "\t", stringsAsFactors = F)
# 每一行求和，并添加在最后一列
asv$sum <- rowSums(asv)
# 根据sum列的大小，对数据进行排序
asv_1 <- asv[order(asv$sum, decreasing = T),]
asv_2 <- asv_1

# 将第一行unclassified 和第7行unidentified合并为一行，并将合并后的值相加，行名改为unclassified 
#asv_2["unclassified",] <- colSums(asv_2[1:2,])
# 删去第7行
#asv_2 <- asv_2[-1,]

# 选取前15行，另存在asv_3中，删去最后一列，作为后续画图使用的数据
asv_3 <- asv_2[1:15,]
asv_3 <- asv_3[,-ncol(asv_3)]
#write.table(asv_3, "C:/Users/15428/Desktop/家庭数据/ASV结果汇总/bac_genus_paixu.txt", sep = "\t", quote = F)
#整理数据，取平均值
# 取分组平均值
samples <- unique(gsub("\\d+", "", colnames(asv_3))) # 提取样品组名（去掉所有数字）
genus_avg <- data.frame(row.names = rownames(asv_3))

for (sample in samples) {
  # 找出当前组的所有列
  sample_cols <- grep(paste0("^", sample, "\\d*$"), colnames(asv_3), value = TRUE)
  
  # 计算当前组的行平均值
  if(length(sample_cols) > 1) {
    genus_avg [[sample]] <- rowMeans(asv_3[, sample_cols, drop = FALSE])
  } else {
    genus_avg [[sample]] <- asv_3[, sample_cols] # 如果只有一列就直接取值
  }
}
# 取平均值
#samples <- unique(gsub("\\d", "", colnames(asv_6w_5))) # 提取样品名称（去掉数字）
#asv_6w_plot <- data.frame(row.names = rownames(asv_6w_5))
#for (sample in samples) {
# 计算每个样品的平均值
#asv_6w_plot[[sample]] <- rowMeans(asv_6w_5[grep(sample, colnames(asv_6w_5))])
#}

#asv_plot <- read.table("C:/Users/15428/Desktop/家庭数据/ASV结果汇总/bac_genus_paixu.txt", header = T, row.names = 1, sep = "\t", stringsAsFactors = F)
# 将bacteria_df按顺序排列
#####排序没问题就不需要这一步
#new_column_order <- c( )
#asv_plot <- asv_plot[, new_column_order]

# 根据asv_6w_plot绘制热图,首先将行名变为一列
genus_avg$Genus <- rownames(genus_avg)
# 将数据从宽格式转换为长格式，适合ggplot
asv_plot_long <- reshape2::melt(genus_avg, id.vars = "Genus")

# 重新设置 Genus 的因子水平为数据中的顺序
asv_plot_long$Genus <- factor(asv_plot_long$Genus, levels = rev(unique(asv_plot_long$Genus)))

write.table(asv_plot_long, "C:/Users/15428/Desktop/高通量测序学习资料/fungus_asv_plot_long.txt", sep = "\t", quote = F)

# 确保数据中有值以生成图例,看相对丰度的范围，根据生成的数值来判断颜色断点及图例刻度
print(summary(asv_plot_long$value))

# 定义颜色断点
breaks <- c(0, 0.25, 0.5) # 你想要颜色改变的相对丰度值，可根据上步生成的数值来定义
colors <- c("#2a4f28", "#faf2ab", "#7d2421") # 对应于断点的颜色

# 图例的刻度通常是由ggplot2自动生成的，定义图例的刻度和标签
legend_breaks <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5)#这个值也是可以改动的，根据breaks中的值
legend_labels <- c("0.0", "0.1", "0.2", "0.3", "0.4", "0.5")


# 绘制热图
p2 <- ggplot(asv_plot_long, aes(x = variable, y = Genus, fill = value)) +
  geom_tile(width = 0.9, height = 0.9) + # 创建热图的方块
  scale_fill_gradientn(colors = colors, values = rescale(breaks)) +
  theme_minimal() + # 使用简约主题
  theme(text = element_text(family = "serif"),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 16, face = "bold", color = "black"),    # 旋转X轴的文本以更好的显示
        axis.text.y = element_text(size = 16, face = "bold", color = "black"),   # 可以根据需要调整Y轴文本的大小
        axis.ticks = element_blank(),   # 移除刻度线
        panel.grid.major = element_blank(),   # 移除主要网格线
        legend.position = "left") + 
  scale_y_discrete(position = "right") +   # 将Y轴标签放在右侧
  labs(x = NULL, y = NULL, fill = "Relative\nAbundance") # 添加和修改图例标题


p2
# 保存图片
ggsave("output/bac属水平热图.pdf", p2, width = 9, height = 8, dpi = 300)
