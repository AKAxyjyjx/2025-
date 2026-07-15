#R包在main里加载了，这里就删掉了

# Alpha多样性是反映样本丰富度和均匀度的综合指标
## 计算alpha多样性 所用数据为相对丰度中所得的相对丰度表格，名字叫relative-XXX.csv的。
asv1 <- read.table("asv_table_only.csv", header=T, sep=",",stringsAsFactors = FALSE) 
# 读取分组文件group。
Group <- read.table("sample_group_info.csv", header = TRUE, sep = ",")
asv1_t <- t(asv1)

# 检查 diversity 来自哪个包
find("diversity")


# 创建函数一步计算alpha多样性
calculate_alpha <- function(otu){
  data_richness <- rowSums(otu>0)
  data_shannon <- vegan::diversity(otu, index = 'shannon', base = exp(1))
  data_pielou <- data_shannon/log(data_richness, base = exp(1))
  data_simpson <- vegan::diversity(otu, index = 'simpson', base = exp(1))
  alpha_matrix <- cbind(data_pielou, data_richness, data_shannon, data_simpson)
  alpha_df <- as.data.frame(alpha_matrix)
  return(alpha_df)
}
# 利用函数，创建绘图所需数据框并重命名其中的列
plot_df1 <- calculate_alpha(asv1_t)%>%
  cbind(Group$Group)%>%
  rename_with(~"Group", 5)

# 逐一绘制预览alpha多样性箱线图
# 先设置颜色
mycol <- met.brewer("Hokusai1",9)

#宽表转化为长表
df_long <- pivot_longer(plot_df1, cols = -Group, names_to = "type", values_to = "alpha_index")

# 将整个ANOVA分析和显著性标记过程打包成一个函数anova_sig()便于后续绘图，函数中还添加了分组最大值和sd的计算
anova_sig <- function(df, alpha_diversity, group){
  anova <- aov(alpha_diversity~group, data = plot_df1)
  pair_comparison <- TukeyHSD(anova)
  pair_comparison_df <- pair_comparison$group
  pair_comparison_df <- as.data.frame(pair_comparison_df)
  group_mean <- aggregate(x = alpha_diversity, by = list(group), FUN = mean)%>%
    rename_with(~c("Group", "mean_val"), 1:2)
  group_max <- aggregate(x = alpha_diversity, by = list(group), FUN = max)%>%
    rename_with(~c("Group", "max"), 1:2)
  group_sd <- aggregate(x = alpha_diversity, by = list(group), FUN = sd)%>%
    rename_with(~c("Group", "sd"), 1:2)
  ntr <- nrow(group_mean)
  mat <- matrix(1, ncol = ntr, nrow = ntr)
  p <- pair_comparison_df$`p adj`
  k <- 0
  for (i in 1:(ntr - 1)) {
    for (j in (i + 1):ntr) {
      k <- k + 1
      mat[i, j] <- p[k]
      mat[j, i] <- p[k]
    }
  }
  treatments <- as.vector(group_mean$Group)
  means <- as.vector(group_mean$mean_val)
  alpha <- 0.05
  pvalue <- mat
  output <- orderPvalue(treatments, means, alpha, pvalue, console = TRUE)
  output$Group <- rownames(output)
  output <- left_join(output, group_max, by = "Group")
  output <- left_join(output, group_sd, by = "Group")
  # 确保每个处理组只返回一行数据
  output <- output %>% 
    distinct(Group, .keep_all = TRUE)
  
  return(output)
}


# 丰富度的ANOVA检验及结果
data_richness <- plot_df1$data_richness
Group = plot_df1$Group
richness_out <- anova_sig(plot_df1, data_richness, Group)
richness_out$type <- "data_richness" #添加一列alpha多样性类别

# 均匀度的ANOVA检验及结果
data_pielou <- plot_df1$data_pielou
Group = plot_df1$Group
pielou_out <- anova_sig(plot_df1, data_pielou, Group)
pielou_out$type <- "data_pielou"

# 香农指数的的ANOVA检验及结果
data_shannon <- plot_df1$data_shannon
Group = plot_df1$Group
shannon_out <- anova_sig(plot_df1, data_shannon, Group)
shannon_out$type <- "data_shannon"

#合并三者结果，但是后面合并四者结果的有了，这个先留着吧
#alpha_out <- rbind(pielou_out, shannon_out, richness_out)%>%rename_with(~"marker", 2)
#将长表与差异分析结果合并
#df_long_all <- left_join(df_long, alpha_out, by = c("type", "Group"))

############之前的代码最后一个图片缺少显著性标识是因为未进行计算，这是进行计算并与之前的进行合并
# Simpson 多样性的 ANOVA 检验及结果
data_simpson <- plot_df1$data_simpson
Group = plot_df1$Group
simpson_out <- anova_sig(plot_df1, data_simpson, Group)
simpson_out$type <- "data_simpson"  # 添加一列 alpha 多样性类别

# 合并四者结果
alpha_out <- rbind(pielou_out, shannon_out, richness_out, simpson_out) %>%
  rename_with(~"marker", 2)

# 将长表与差异分析结果合并
df_long_all <- left_join(df_long, alpha_out, by = c("type", "Group"))

# 检查并处理 marker 列的缺失值
df_long_all <- df_long_all %>%
  mutate(marker = ifelse(is.na(marker), "", marker))
#




# 将 Group 列转换为因子，并设置水平顺序
df_long_all$Group <- factor(df_long_all$Group, levels = c("EL","EM","EH","ML","MM","MH","LL","LM","LH"))

# 按顺序指定颜色
mycol <- met.brewer("Hokusai1", 9)  # 确保颜色数量与分组数量一致
names(mycol) <- c("EL","EM","EH","ML","MM","MH","LL","LM","LH")

# 拆分 marker 列中的多个字母
df_long_all <- df_long_all %>%
mutate(marker = strsplit(as.character(marker), ", ")) %>%  # 按逗号拆分
unnest(marker) %>%  # 将拆分的字母展开为多行
group_by(Group, type) %>%
mutate(text_y = max + sd + (row_number() - 1) * 0.1) %>%  # 调整文本的垂直位置
ungroup()


# 修改后续的数据处理步骤，不再拆分标记
df_long_all <- left_join(df_long, alpha_out, by = c("type", "Group")) %>%
  # 处理缺失的标记
  mutate(marker = ifelse(is.na(marker), "", marker)) %>%
  # 为每个Group-type组合计算一个text_y位置
  group_by(type, Group) %>%
  mutate(text_y = max + sd + 0.1) %>%  # 统一调整文本位置
  ungroup()

# 绘制分面箱线图
p <- ggplot(df_long_all, aes(Group, alpha_index)) +
  scale_fill_manual(values = mycol) +  # 按指定顺序分配颜色
  geom_boxplot(aes(fill = Group), lwd = 0.8, fatten = 1.5) +  # 调整箱体线条和中位线粗细
  theme_bw() +
  labs(x = "Treatment", y = "Alpha diversity") +
  scale_x_discrete(limits = c("EL","EM","EH","ML","MM","MH","LL","LM","LH")) +
  geom_text(
    aes(x = Group, y = text_y, label = marker),  # 使用调整后的 text_y
    size = 5.5, 
    fontface = "bold", 
    family = "serif",  # 设置差异标识为 serif 字体
    vjust = 0  # 调整文本的垂直对齐方式
  ) +
  facet_wrap(. ~ type, scales = "free_y") +
  theme(
    text = element_text(family = "serif"),  # 全局字体设置
    strip.text = element_text(size = 14, face = "bold", color = "black"), # 分面字体调整
    legend.title = element_text(size = 14, face = "bold", colour = "black"),
    legend.text = element_text(size = 14, face = "bold", colour = "black"),  # 右侧标签字体加粗
    axis.text.x = element_text(size = 16, face = "bold", colour = "black"),   # x轴标签正常横向显示
    axis.text.y = element_text(size = 16, face = "bold", colour = "black"),   # 加粗并放大y轴标签
    axis.title.x = element_text(size = 20, face = "bold", colour = "black"),  # 加粗并放大x轴标题
    axis.title.y = element_text(size = 20, face = "bold", colour = "black"),  # 加粗并放大y轴标题
    panel.border = element_rect(colour = "black", fill = NA, size = 2),  # 加粗面板边框
    panel.grid.major = element_line(colour = "grey90", size = 0.75, linetype = "solid"),  # 设置主要网格线
    panel.grid.minor = element_line(colour = "grey95", size = 0.5, linetype = "solid"),  # 设置次要网格线
  )

# 显示图形
print(p)




# 将图片以.pdf的格式导出,图像堆积时记得调整大小
ggsave("output/Alpha多样性", width = 28, height = 14, units = "in")
