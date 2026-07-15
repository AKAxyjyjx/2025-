

# 基于Bray-Curtis距离的PCoA分析
# 载入数据，先载入要测的相对丰度表,同Alpha多样性
otu <- read.table("temperary_data/relative_table_only.csv", header = TRUE, row.names = 1, sep = "\t")
# 载入分组信息
group <- read.table("temperary_data/sample_group_info.csv", header = TRUE, sep = ",")

# 计算Bray-Curtis距离，dist表示该对象为一个距离矩阵
distance_bray <- vegdist(t(otu), method="bray",diag=T, upper=T)
class(distance_bray)

# 把dist转换成matrix
distance_bray_matrix <- as.matrix(distance_bray)
dim(distance_bray_matrix)

# 再把matrix转换成dataframe
distance_bray_frame <- as.data.frame(distance_bray_matrix)
dim(distance_bray_frame)

# 基于Bray-Curtis距离做PCoA分析
pcoa_bray <- cmdscale(distance_bray_frame, k = 2,eig=T)  #默认维度为k = 2
pcoa_bray$points  #查看每个样本在坐标的位置

# 获取各维度的占比与解释率的dataframe
eig_pcoa_bray <- summary(eigenvals(pcoa_bray))
head(eig_pcoa_bray)

# 主坐标分析（PCoA）的结果，其中包含了特征值（Eigenvalue）、解释的方差比例（Proportion Explained）以及累积方差比例（Cumulative Proportion）
# 特征值（Eigenvalue）代表了每个主坐标轴的方差量，即该轴上数据分布的广度，一个较高的特征值意味着该维度能够捕获更多的数据变异性。
# 解释的方差比例（Proportion Explained）表示每个主坐标轴解释了总变异性的多少比例，有助于我们理解哪些轴对群落差异贡献最大。
# 累积方差比例（Cumulative Proportion）通常显示了从第一个轴到当前轴为止的总方差解释比例的累积，通常这应该是递增的数值，反映了随着维度的增加，累积解释的总变异性的增长。

# 设置各维度的名字，从PCoA1开始
axis_pcoa_bray <- paste0("PCo", 1:ncol(eig_pcoa_bray))
axis_pcoa_bray

# 生成各轴的解释度数据，合并维度名字和数据
eig_pcoa_bray <- data.frame(Axis = axis_pcoa_bray, t(eig_pcoa_bray)[, -3])
head(eig_pcoa_bray)

# 获取前两个维度的解释率并保留两位小数
pco1_bray = round(eig_pcoa_bray[1, 3] * 100, 2)
pco2_bray = round(eig_pcoa_bray[2, 3] * 100, 2)

# 设置画图时的x轴和y轴的标题，衔接逗号中的内容
xlab_pcoa_bray = paste0("PCo1 (",pco1_bray,"%)")
ylab_pcoa_bray = paste0("PCo2 (",pco2_bray,"%)")

# 获取各样本在前两个维度的坐标
pcoa_bray_points = as.data.frame(pcoa_bray$points)[,c(1:2)]
head(pcoa_bray_points)

# 合并样本坐标和对应的分组信息
pcoa_bray_points = data.frame(pcoa_bray_points, 
                              group)
head(pcoa_bray_points)

# 进行Adonis分析（PERMANOVA），导入的数据为距离矩阵
adonis_bray <- adonis2(distance_bray ~ Group, data = group, permutations = 999)
print(adonis_bray)
class(adonis_bray)
# Df (Degrees of freedom)：自由度，Management变量的自由度为11（表示Management变量的类别数减去1），剩余的自由度（Residual）为24（表示总样本数减去Management类别数）。
# SumOfSqs (Sum of Squares)：平方和，表示由Management变量解释的多样性（或距离）变异量为8.7370，剩余未解释的变异量为2.4796，总变异量（Total）为11.2165。
# R2：决定系数，表示Management变量解释的总变异比例为77.894%（0.77894），剩余未解释的变异比例为22.106%（0.22106）。R2值越接近1，说明模型解释的变异比例越高。
# F：F统计量，用于衡量Management变量解释的变异与剩余变异之间的比例大小，值为7.6878。F值较大通常意味着解释变量对样本间距离的影响较为显著。
# Pr(>F)：P值，用于测试Management变量对样本间距离是否有显著影响，值为0.001。在这里，***表示P值非常小（小于0.001），意味着Management变量显著影响了样本间的距离，或者说，不同Management分组间的群落组成存在显著差异。

# 提取Adonis分析的P值和R2值，存储在r_squared和p_value中
r_square <- adonis_bray$R2
p_value <- adonis_bray$`Pr(>F)`[1]
if (p_value < 0.001) {
  p_text <- "P < 0.001"
} else {
  p_text <- sprintf("P = %.3f", p_value)
}

# 添加椭圆，注意，这些方法提供的椭圆只是近似表示，不应该被解释为严格的统计置信区间，三个平行加不了
p_pcoa_bray+tat_ellipse(aes(fill = Management, color = Management), geom = 'polygon', level = 0.95, alpha = 0.1,
show.legend = FALSE)
#————————————————————————————————————————————————————————————————————————————————————————————————————————————#
#————————————————————————————————————————————————————————————————————————————————————————————————————————————#
# 基于Bray-Curtis距离做NMDS分析
NMDS_bray <- metaMDS(distance_bray_frame, k = 2,eig=T)  #默认维度为k = 2
NMDS_bray$points  #查看每个样本在坐标的位置
NMDS_bray$stress 
# NMDS不能给出每个轴的解释度，但需要标明stress
#对于NMDS二维分析，通常认为stress<0.2时有一定的解释意义；
#当stress<0.1时，可认为是一个好的排序；当 stress<0.05时，则具有很好的代表性。

# 获取各样本在前两个维度的坐标
NMDS_bray_points = as.data.frame(NMDS_bray$points)[,c(1:2)]
head(NMDS_bray_points)
#合并样本坐标和对应的分组信息
NMDS_bray_points = data.frame(NMDS_bray_points, 
                              group)
head(NMDS_bray_points)
#绘图
#通常需要标记stress信息，不标记轴的权重信息。



# 按顺序指定颜色
mycol <- met.brewer("Hokusai1", 9)  # 确保颜色数量与分组数量一致
names(mycol) <- c("EL","EM","EH", "ML", "MM", "MH", "LL", "LM", "LH")

# 确保 Management 列是因子，并按照指定顺序设置水平
NMDS_bray_points$Group <- factor(NMDS_bray_points$Group, levels = c("EL","EM","EH", "ML", "MM", "MH", "LL", "LM", "LH"))

# 首先计算统一的绘图范围
x_combined_range <- range(c(pcoa_bray_points$V1, NMDS_bray_points$MDS1))
y_combined_range <- range(c(pcoa_bray_points$V2, NMDS_bray_points$MDS2))

# 扩展10%的边距
x_expanded <- x_combined_range + c(-0.1, 0.1)*diff(x_combined_range)
y_expanded <- y_combined_range + c(-0.1, 0.1)*diff(y_combined_range)






# 1. 修改NMDS图
p_NMDS_bray <- ggplot(NMDS_bray_points, aes(MDS1, MDS2)) +
  geom_point(size = 6, aes(shape = Group, color = Group), alpha = 0.7) +
  scale_shape_manual(values = c(15, 16, 17, 18, 15, 16, 17, 18, 15, 16, 17, 18)) +
  scale_color_manual(values = mycol) +
  labs(
    x = "NMDS1", 
    y = "NMDS2",
    title = "NMDS"
  ) +
  geom_hline(yintercept = 0, linetype = 4) +
  geom_vline(xintercept = 0, linetype = 4) +
  coord_cartesian(xlim = x_expanded, ylim = y_expanded) + # 使用统一的范围
  theme_bw() +
  theme(
    text = element_text(family = "serif", face = "bold"),
    aspect.ratio = 1,
    plot.title = element_text(size = 16, hjust = 0.5, colour = "black"),
    axis.title = element_text(size = 20, colour = "black"),
    axis.text = element_text(size = 16, colour = "black"),
    legend.title = element_text(size = 16, colour = "black"),
    legend.text = element_text(size = 16, colour = "black"),
    legend.position = "right",
    panel.border = element_rect(colour = "black", fill = NA, size = 2),
    panel.grid.major = element_line(colour = "grey90", size = 0.75),
    panel.grid.minor = element_line(colour = "grey95", size = 0.5),
    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm") # 统一边距
  ) +
  annotate(
    "text", 
    label = "stress = 0.127", 
    x = min(x_expanded) + 0.4*diff(x_expanded),
    y = max(y_expanded) - 0.05*diff(y_expanded),
    size = 5, 
    family = "serif",
    fontface = "bold"
  )

##如果需要标签附着
#library(ggrepel)  # 用于避免标签重叠 #在main里已经装了，这里就是说明一下

p_NMDS_bray <- ggplot(NMDS_bray_points, aes(MDS1, MDS2)) +
  geom_point(size = 6, aes(shape = Group, color = Group), alpha = 0.7) +
  # 添加样本标签
  geom_text_repel(aes(label = rownames(NMDS_bray_points)),# 假设行名是样本ID
  size = 4,
  family = "serif",
  fontface = "bold",
  box.padding = 0.5,
  max.overlaps = 100
) +
  scale_shape_manual(values = c(15, 16, 17, 18, 15, 16, 17, 18, 15, 16, 17, 18)) +
  scale_color_manual(values = mycol) +
  labs(
    x = "NMDS1", 
    y = "NMDS2",
    title = "NMDS"
  ) +
  geom_hline(yintercept = 0, linetype = 4) +
  geom_vline(xintercept = 0, linetype = 4) +
  coord_cartesian(xlim = x_expanded, ylim = y_expanded) + # 使用统一的范围
  theme_bw() +
  theme(
    text = element_text(family = "serif", face = "bold"),
    aspect.ratio = 1,
    plot.title = element_text(size = 16, hjust = 0.5, colour = "black"),
    axis.title = element_text(size = 20, colour = "black"),
    axis.text = element_text(size = 16, colour = "black"),
    legend.title = element_text(size = 16, colour = "black"),
    legend.text = element_text(size = 16, colour = "black"),
    legend.position = "right",
    panel.border = element_rect(colour = "black", fill = NA, size = 2),
    panel.grid.major = element_line(colour = "grey90", size = 0.75),
    panel.grid.minor = element_line(colour = "grey95", size = 0.5),
    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm") # 统一边距
  ) +
  annotate(
    "text", 
    label = "stress = 0.127", 
    x = min(x_expanded) + 0.4*diff(x_expanded),
    y = max(y_expanded) - 0.05*diff(y_expanded),
    size = 5, 
    family = "serif",
    fontface = "bold"
  )

  # 其余主题设置保持不变...
print(p_NMDS_bray)






# 2. 修改PCoA图使用相同参数
# 指定配色方案
mycol <- met.brewer("Hokusai1",9)
mycol

# 按顺序指定颜色
mycol <- met.brewer("Hokusai1", 9)  # 确保颜色数量与分组数量一致
names(mycol) <- c("EL", "EM", "EH", "ML", "MM", "MH", "LL", "LM", "LH")

# 确保 Management 列是因子，并按照指定顺序设置水平
pcoa_bray_points$Group <- factor(pcoa_bray_points$Group, levels = c("EL","EM","EH", "ML", "MM", "MH", "LL", "LM", "LH"))


p_pcoa_bray <- ggplot(pcoa_bray_points, aes(V1, V2)) +
  geom_point(size = 6, aes(shape = Group, color = Group), alpha = 0.7) +
  scale_shape_manual(values = c(15, 16, 17, 18, 15, 16, 17, 18, 15, 16, 17, 18)) +
  scale_color_manual(values = c(rev(mycol), "gray")) +
  labs(x = xlab_pcoa_bray, y = ylab_pcoa_bray, title = "PCoA") +
  geom_hline(yintercept = 0, linetype = 4) +
  geom_vline(xintercept = 0, linetype = 4) +
  coord_cartesian(xlim = x_expanded, ylim = y_expanded) + # 使用相同范围
  theme_bw() +
  theme(
    text = element_text(family = "serif", face = "bold"),
    aspect.ratio = 1,
    axis.title = element_text(size = 20, colour = "black"),
    axis.text = element_text(size = 16, colour = "black"),
    plot.title = element_text(size = 16, hjust = 0.5, colour = "black"),
    legend.title = element_text(size = 16, colour = "black"),
    legend.text = element_text(size = 16, colour = "black"),
    legend.position = "right", # 与NMDS一致
    panel.border = element_rect(colour = "black", fill = NA, size = 2),
    panel.grid.major = element_line(colour = "grey90", size = 0.75),
    panel.grid.minor = element_line(colour = "grey95", size = 0.5),
    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm") # 统一边距
  ) +
  annotate("text", 
           x = min(x_expanded) + 0.45*diff(x_expanded),
           y = max(y_expanded) - 0.15*diff(y_expanded),
           label = sprintf("%s\nR² = %.3f", p_text, r_square[1]),
           size = 5,
           family = "serif",
           fontface = "bold")


##同样的如果需要标签附着
p_pcoa_bray <- ggplot(pcoa_bray_points, aes(V1, V2)) +
  geom_point(size = 6, aes(shape = Group, color = Group), alpha = 0.7) +
  # 添加样本标签（使用ggrepel避免重叠）
  geom_text_repel(
    aes(label = rownames(pcoa_bray_points)),  # 假设行名是样本ID
    size = 4,
    family = "serif",
    fontface = "bold",
    box.padding = 0.5,  # 调整标签与点的间距
    max.overlaps = 100  # 允许更多重叠检查
  ) +
  scale_shape_manual(values = c(15, 16, 17, 18, 15, 16, 17, 18, 15, 16, 17, 18)) +
  scale_color_manual(values = c(rev(mycol), "gray")) +
  labs(x = xlab_pcoa_bray, y = ylab_pcoa_bray, title = "PCoA") + scale_color_manual(values = c(rev(mycol), "gray")) +
  labs(x = xlab_pcoa_bray, y = ylab_pcoa_bray, title = "PCoA") +
  geom_hline(yintercept = 0, linetype = 4) +
  geom_vline(xintercept = 0, linetype = 4) +
  coord_cartesian(xlim = x_expanded, ylim = y_expanded) + # 使用相同范围
  theme_bw() +
  theme(
    text = element_text(family = "serif", face = "bold"),
    aspect.ratio = 1,
    axis.title = element_text(size = 20, colour = "black"),
    axis.text = element_text(size = 16, colour = "black"),
    plot.title = element_text(size = 16, hjust = 0.5, colour = "black"),
    legend.title = element_text(size = 16, colour = "black"),
    legend.text = element_text(size = 16, colour = "black"),
    legend.position = "right", # 与NMDS一致
    panel.border = element_rect(colour = "black", fill = NA, size = 2),
    panel.grid.major = element_line(colour = "grey90", size = 0.75),
    panel.grid.minor = element_line(colour = "grey95", size = 0.5),
    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm") # 统一边距
  ) +
  annotate("text", 
           x = min(x_expanded) + 0.45*diff(x_expanded),
           y = max(y_expanded) - 0.15*diff(y_expanded),
           label = sprintf("%s\nR² = %.3f", p_text, r_square[1]),
           size = 5,
           family = "serif",
           fontface = "bold")

  # 其余主题设置保持不变...
#为了生成多个不同版本的图片，我这里不会注释掉其中一个版本，由于相互覆盖，最后会是带标签的。
print(p_pcoa_bray)


combined <- (p_pcoa_bray + theme(legend.position = "none")) | 
  (p_NMDS_bray + theme(legend.position = "right"))
print(combined)

ggsave("output/细菌beta多样性.pdf", width = 28, height = 14, units = "in")

