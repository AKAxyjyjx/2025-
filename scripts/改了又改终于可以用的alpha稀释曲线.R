

raw <- read.delim("asv_table_only.csv", sep = ',', stringsAsFactors = FALSE, check.names = FALSE)
# 提取样品名称，保留前两个字符
for (i in 1:ncol(raw)) {
samples <- unique(substr(colnames(raw), 1,2))

}
asv <- data.frame(row.names = rownames(raw))

for (sample in samples) {
  # 计算每个样品的平均值
  asv[[sample]] <- rowMeans(raw[,grep(sample, colnames(raw))])
}

# 将asv中的数值转换为integers
asv <- as.data.frame(lapply(asv, as.integer))

otu <- t(asv)


# 假设 raw 是ASV表（行=ASV，列=样本）统计样本深度
sample_depths <- colSums(raw)
print(sample_depths)



#定义函数
#计算多种Alpha多样性指数，结果返回至向量
alpha_index <- function(x, method = 'richness', tree = NULL, base = exp(1)) {
  if (method == 'richness') result <- rowSums(x > 0)    #丰富度指数
  else if (method == 'chao1') result <- estimateR(x)[2, ]    #Chao1 指数
  else if (method == 'ace') result <- estimateR(x)[4, ]    #ACE 指数
  else if (method == 'shannon') result <- diversity(x, index = 'shannon', base = base)    #Shannon 指数
  else if (method == 'simpson') result <- diversity(x, index = 'simpson')    #Gini-Simpson 指数
  else if (method == 'pielou') result <- diversity(x, index = 'shannon', base = base) / log(estimateR(x)[1, ], base)    #Pielou 均匀度
  else if (method == 'gc') result <- 1 - rowSums(x == 1) / rowSums(x)    #goods_coverage
  else if (method == 'pd' & !is.null(tree)) {    #PD_whole_tree
    pd <- pd(x, tree, include.root = FALSE)
    result <- pd[ ,1]
    names(result) <- rownames(pd)
  }
  result
}

#根据抽样步长，统计每个稀释梯度下的Alpha多样性指数，结果返回至列表
alpha_curves <- function(x, step, method = 'richness', rare = NULL, tree = NULL, base = exp(1)) {
  x_nrow <- nrow(x)
  if (is.null(rare)) rare <- rowSums(x) else rare <- rep(rare, x_nrow)
  alpha_rare <- list()
  
  for (i in 1:x_nrow) {
    step_num <- seq(0, rare[i], step)
    if (max(step_num) < rare[i]) step_num <- c(step_num, rare[i])
    
    alpha_rare_i <- NULL
    for (step_num_n in step_num) alpha_rare_i <- c(alpha_rare_i, alpha_index(x = rrarefy(x[i, ], step_num_n), method = method, tree = tree, base = base))
    names(alpha_rare_i) <- step_num
    alpha_rare <- c(alpha_rare, list(alpha_rare_i))
  }
  
  names(alpha_rare) <- rownames(x)
  alpha_rare
}

#如果基于一次抽样的结果可能存在误差，希望在同深度下多次抽样并统计均值和标准差，并最终绘制成带有误差棒的曲线图
##多计算几次以获取均值 ± 标准差，然后再展示出也是一个不错的选择
#重复抽样 5 次
plot_richness <- data.frame()
## 6w深度的步长设置为2000，20w深度的设置为6000
for (n in 1:5) {
  richness_curves <- alpha_curves(otu, step = 1000, method = 'richness')
  
  for (i in names(richness_curves)) {
    richness_curves_i <- (richness_curves[[i]])
    richness_curves_i <- data.frame(rare = names(richness_curves_i), alpha = richness_curves_i, sample = i, stringsAsFactors = FALSE)
    plot_richness <- rbind(plot_richness, richness_curves_i)
  }
}

#计算均值 ± 标准差（doBy 包中的 summaryBy() 函数）
#install.packages("doBy")

plot_richness_stat <- summaryBy(alpha~sample+rare, plot_richness, FUN = c(mean, sd))
plot_richness_stat$rare <- as.numeric(plot_richness_stat$rare)
plot_richness_stat[which(plot_richness_stat$rare == 0),'alpha.sd'] <- NA


mycol <- met.brewer("Hokusai1",9)#几个处理设置几个颜色
mycol
# 定义处理顺序
treatment_order <- c("EL","EM", "EH", "ML", "MM", "MH", "LL", "LM", "LH")

# 将sample列转换为因子，并设置顺序
plot_richness_stat$sample <- factor(plot_richness_stat$sample, levels = treatment_order)

# 确保rare列是数值类型
plot_richness_stat$rare <- as.numeric(plot_richness_stat$rare)




# 计算每个样本的饱和点（95%阈值）
sat_points <- plot_richness_stat %>%
  group_by(sample) %>%
  mutate(max_alpha = max(alpha.mean)) %>%
  filter(alpha.mean >= 0.95 * max_alpha) %>%
  summarise(sat_depth = min(rare))

# 打印结果
print(sat_points)




ggplot(plot_richness_stat, aes(rare, alpha.mean, color = sample)) +
  ##x 轴和 y 轴的数据分别来自 plot_richness_stat 数据框中的 rare 和 alpha.mean 列
  geom_line() +
  geom_point() +
  geom_errorbar(aes(ymin = alpha.mean - alpha.sd, ymax = alpha.mean + alpha.sd), width = 500) +
  labs(x = 'Number of sequences', y = 'Richness', color = 'Treatment') +
  theme(panel.grid = element_blank(), 
        panel.background = element_rect(fill = 'transparent', color = 'black'), 
        legend.key = element_rect(fill = 'transparent')) +
  geom_vline(xintercept = 29000, linetype = 2) +   # 饱和点需要另外计算
  scale_x_continuous(breaks = seq(0, 50000, 25000), labels = as.character(seq(0, 50000,25000))) + 
  scale_y_continuous(breaks = seq(0, 2500, 2000), labels = as.character(seq(0, 2500, 2000)), limits = c(0, 2500)) +  
  scale_color_manual(values = c(rev(mycol)),
                     breaks = treatment_order  # 关键修复：按指定顺序排列图例
  ) +
  theme(
    # 全局字体设置为 serif##
    text = element_text(family = "serif"),
    axis.text.x = element_text(size = 16, face = "bold", color = "black"),
    axis.text.y = element_text(size = 16, face = "bold", color = "black"),
    axis.title.x = element_text(size = 20, face = "bold"), 
    axis.title.y = element_text(size = 20, face = "bold"),
    legend.title = element_text(face = "bold", size = 14),   # 修改图例标题的字体
    legend.text = element_text(size = 14, face = "bold", color = "black"),    # 修改图例文本的字体
    # axis.line = element_line(size = 1),
    panel.border = element_rect(colour = "black", fill = NA, size = 2),  # 加粗面板边框
    panel.grid.major = element_line(colour = "grey90", size = 0.75, linetype = "solid"),  # 设置主要网格线
    panel.grid.minor = element_line(colour = "grey95", size = 0.5, linetype = "solid"))  # 设置次要网格线
# 导出图片
ggsave("output/alpha稀释曲线.pdf", width = 8, height = 6, dpi = 300)
