
# 加载数据，包括ASV表（2中去除线粒体、叶绿体的表格）、样本分组信息以及环境因子数据
data <- read.table(file="asv_table_only.csv", header=T,check.names=T ,row.names=1)
group <- read.table(file="sample_group_info.csv", header=T,check.names=T,row.names=NULL )
rownames(group) <- group$sample
env <- read.table(file="env_data.csv",header=T,check.names=T, fileEncoding = "UTF-8", sep = ",", row.names = 1)
sampledata <- t(data)


# 标准化数据，添加分类颜色
# hellinger转换
sampledata <- decostand(sampledata,method = "hellinger")
#定义分组的填充颜色

mycol <- met.brewer("Hokusai1",9)#几个处理设置几个颜色
mycol
col <- mycol
#col <- c("#E7B800","#43CD80","#2E9FDF")

# 进行RDA分析之前，先进行DCA分析,根据DCA结果中的Axis lengths的数值来进行判断:
# 如果其中最大的数值大于4，则应选择CCA，如果最大的数值小于3，则选择RDA，如果最大的数值在3-4之间，则两种分析方法都可以。
dca = decorana(veg = sampledata)
dca

# RDA分析    
rda = rda(sampledata, env, scale = TRUE)
# RDA得分
rda.sample=data.frame(rda$CCA$u[,1:2]) #样本得分
rda.sample$group = group$Group#添加分组
#write.csv(rda.sample,file="C:/Users/15428/Desktop/家庭数据/RDA_fungus_meihuo_result.csv")
rda.env=data.frame(rda$CCA$ biplot[,1:2]) #环境因子得分
rda1 =round(rda$CCA$eig[1]/sum(rda$CCA$eig)*100,2) # 第一轴解释量
rda2 =round(rda$CCA$eig[2]/sum(rda$CCA$eig)*100,2) # 第二轴解释量
#置换检验 
envfit <- envfit(rda, env, permutations  = 999)
r <- as.matrix(envfit$vectors$r)
p <- as.matrix(envfit$vectors$pvals)
env.p <- cbind(r,p)
colnames(env.p) <- c("r2","p-value")
KK <- as.data.frame(env.p)
KK$p.adj = p.adjust(KK$`p-value`, method = 'BH')
KK

#1. 将分组转换为因子并指定顺序
rda.sample$group <- factor(rda.sample$group, 
                           levels = c("EL","EM","EH","ML","MM","MH","LL","LM","LH"))

# 2. 定义颜色并命名（确保与分组顺序一致）
mycol <- met.brewer("Hokusai1", 9)
names(mycol) <- levels(rda.sample$group)  # 直接使用因子的levels保证顺序匹配



# 绘图（所有字体设置为serif+加粗）
ggplot(rda.sample, aes(RDA1, RDA2)) +
  # 样本点
  geom_point(aes(fill = group, color = group), size = 3, shape = 21) + 
  scale_color_manual(values = col) +
  scale_fill_manual(values = col) +
  # 坐标轴标签（带百分比）
  xlab(paste("RDA1 ( ", rda1, "%", " )", sep = "")) + 
  ylab(paste("RDA2 ( ", rda2, "%", " )", sep = "")) +
  # 分组椭圆
  geom_mark_ellipse(aes(fill = group), alpha = 0.1, tol = 0.6, expand = unit(0, "mm")) + 
  coord_cartesian(clip = "off") +
  # 环境因子箭头
  geom_segment(
    data = rda.env, 
    aes(x = 0, y = 0, xend = rda.env[,1], yend = rda.env[,2]),
    arrow = arrow(length = unit(0.35, "cm"), type = "closed", angle = 22.5),
    linetype = 1, colour = "black", size = 0.6
  ) + 
  # 环境因子标签
  geom_text_repel(
    data = rda.env,
    aes(x = rda.env[,1], y = rda.env[,2], label = rownames(rda.env)),
    size = 6, segment.colour = "black",
    family = "serif", fontface = "bold"  # 新增：标签字体设置
  ) +
  # 参考线
  geom_vline(aes(xintercept = 0), linetype = "dotted") +
  geom_hline(aes(yintercept = 0), linetype = "dotted") +
  # 主题设置（统一字体）
  theme_bw() +
  theme(
    # 全局字体
    text = element_text(family = "serif", face = "bold"),
    # 坐标轴标题
    axis.title = element_text(size = 20, colour = "black"),
    # 坐标轴刻度
    axis.text = element_text(size = 16, color = "black"),
    # 图例文字
    legend.text = element_text(size = 14),
    # 图标题
    plot.title = element_text(size = 20, hjust = 0.5),
    panel.grid = element_blank(),
    legend.position = "right"
  ) +
  # 标题（字体已通过theme统一设置）
  labs(title = "Bacteria")


# 保存图片
ggsave("C:/Users/15428/Desktop/家庭数据/ASV结果汇总/RDA_fungus_meihuo_Zn.pdf",width=8,height=8,dpi=300)












