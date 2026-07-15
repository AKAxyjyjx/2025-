
# 读取otu表！可使用去除叶绿体、线粒体等的相对丰度或者绝对丰度的ASV表格(2中结果)，细菌和真菌放在一个文件夹里面，注意各有多少列
otu <- read.table("asv_table_only.csv", header = T,  sep = ',')
env <- read.table("env_data.csv", header = TRUE,
                  sep = ",")


otu_beta <- as.data.frame(t(otu))
# 方法2：Sample 设为行名
rownames(env) <- env$Sample
env <- env[, -1]
# 进行Mantel检验
mantel = mantel_test(
  spec = otu_beta, env = env,
  #spec_select = list(bacteria =1:13556,fungi = 13557:19053), #根据细菌真菌ASV数量进行修改
  spec_dist =  dist_func(.FUN = "vegdist", method = "bray"), # 样本距离使用的vegdist()计算，可以选择适合自己数据的距离指数。
  env_dist = dist_func(.FUN = "vegdist", method = "euclidean"),
  mantel_fun = 'mantel', # mantel.partial：partial mantel test则需要设置env_ctrl
  na_omit=TRUE,
)
mantel # 查看环境因子与细菌、真菌群落的mantel相关的r和p值。
mantel$P_adj_BH <- p.adjust(mantel$p, method = 'BH')#矫正P
# 导出mantel
#write.xlsx(mantel,"C:/Users/15428/Desktop/高通量测序学习资料/12Mantel分析_去除TC-05.27_mantel.xlsx")

#设置相关系数、显著性标签
mantel = mutate(mantel, 
                r = cut(r, right = TRUE,# 表示分割区间形式为(a1,a2],前开后闭。
                        breaks = c(-Inf, 0.4, 0.6, 0.8, Inf),
                        labels = c('< 0.4', '0.4 - 0.6', '0.6 - 0.8', '>= 0.8')),
                p = cut(p, right = FALSE,# 表示分割区间形式为[a1,a2)。
                        breaks = c(-Inf, 0.001, 0.01, 0.05, Inf), 
                        labels = c('***', '< 0.01', '0.01 - 0.05', '>= 0.05'),))
mantel = data.frame(mantel,stringsAsFactors = FALSE)

# 开始绘图

p <- 
  # 计算与绘制土壤因子之间相关性热图
  correlate(env, method = "pearson") %>% 
  qcorrplot(type = "upper",  # 上三角
            diag = FALSE,    # 去除自相关
            grid_col = "black", 
            grid_size = 0.25) +
  geom_square() +
  # 添加r值与显著性标记
  geom_mark(
    sep = '\n', 
    size = 4,
    sig_level = c(0.05, 0.01, 0.001),
    sig_thres = 0.05,
    only_mark = TRUE,
    colour = "white"
  ) +
  scale_fill_gradientn(
    colours = rev(RColorBrewer::brewer.pal(11, "RdBu")),
    limits = c(-1, 1),
    breaks = seq(-1, 1, 0.5)
  ) +
  # 绘制连线
  geom_couple(aes(color = p, size = r), 
              data = mantel, 
              label.size = 4,
              drop = TRUE,
              label.colour = "black",
              label.fontface = "bold",  # 标签加粗
              nudge_x = 1,
              curvature = 0.1) +
  scale_size_manual(values = c(0.2, 0.5, 0.8, 1.2)) +
  scale_color_manual(values = c("#D95F02", "#1B9E77", "#A2A2A2")) +
  guides(
    size = guide_legend(title = "Mantel's r",
                        override.aes = list(colour = "grey35"),
                        order = 1),
    colour = guide_legend(title = "Mantel's p",
                          override.aes = list(size = 3),
                          order = 2),
    fill = guide_colorbar(title = "Pearson's r", order = 3)
  ) +
  theme(
    # 全局字体设置
    text = element_text(family = "serif", face = "bold"),
    # 坐标轴文字
    axis.text = element_text(size = 14, colour = "black"),
    axis.text.x.top = element_text(angle = 45, vjust = 0),
    # 图例
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12),
    legend.key = element_blank()
  )



p

# 保存图片
ggsave("Cd-Mantel-test.pdf", width = 9, height = 8, units = "in")
