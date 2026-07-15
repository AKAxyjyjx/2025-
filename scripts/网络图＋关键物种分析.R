

# 读取OTU/ASV相对丰度/绝对含量的原始数据表格
otu <- read.table("temperary_data/asv_table_only.csv",sep=",", header=T)
group <- read.table("temperary_data/sample_group_info.csv",sep=",", header=T)
# ============================================
# 从ASV表中提取单个组的数据
# ============================================
# 获取 LH 组的样本名
samples_LH <- group$Sample[group$Group == "LH"]
# 如果 otu 行名是数字，用 group 的行名或索引来匹配
group_LH <- otu[, group$Group == "LH"]

otu1 <- group_LH
# 如果你的 ASV 名称是 ASV1, ASV2, ...
rownames(otu1) <- paste0("ASV", 1:nrow(otu1))

# 查看
head(rownames(otu1))
# 将丰度值大于0的值替换为1，便于计算不同属的覆盖度
otu1[otu1>0] <- 1
#  例如只保留在3/5个及以上样本中出现的属,保留多少的属，下面对应更改！
otu2 <- otu1[which(rowSums(otu1) >=5), ]

# 1. 从 otu 原始数据恢复
if (exists("otu1") && nrow(otu1) == nrow(otu2)) {
  rownames(otu2) <- rownames(otu1)
}
# 加上前面的筛选，这里剩下的属是样品覆盖度为3/5，导出筛选的数据
write.table(otu2, file="temperary_data/bacteria_LH_sample.txt", quote=F, sep="\t", row.names=T, col.names=NA)

# 读取矩阵，行为sample，列为genus，就是刚才导出的表格
otu3 <- otu2

# adjust为校准r值，alpha为显著性水平，这里使用FDR校准，显著性水平为0.05
# 这一步耗时最久，与筛选后的样本量有关
occor <- corr.test(t(otu3),method="spearman",adjust="fdr",alpha=0.05)
# 取相关性矩阵R值
occor.r <-  occor$r
# 取相关性矩阵p值
occor.p <- occor$p
# p 值校正，这里使用 BH 法校正 p 值
p <- p.adjust(occor.p, method = 'BH') 
# 确定物种间存在相互作用关系的阈值，将相关性R矩阵内不符合的数据转换为0
occor.r[occor.p>0.05|abs(occor.r)<0.6] = 0
diag(occor.r) <- 0
# 将occor.r保存
write.table(occor.r, file="temperary_data/bacteria_LH_相关性计算结果.txt", quote=F, sep="\t", row.names=T, col.names=NA)
# 根据上述筛选的 r 值和 p 值保留数据
z <- occor.r * occor.p
# 将相关矩阵中对角线中的值（代表了自相关）转为 0
diag(z) <- 0
head(z)[1:6,1:6]
z
dim(z)
# 得到邻接矩阵格式的网络文件（微生物属的相关系数矩阵）
write.table(z, file="temperary_data/bacteria_LH_绘图矩阵.txt", quote=F, sep="\t", row.names=T, col.names=NA)

##获得网络
# 将邻接矩阵转化为 igraph 网络的邻接列表
# 构建含权的无向网络，权重代表了微生物属间丰度的 spearman 相关系数 
g <- graph.adjacency(z, weighted = TRUE, mode = 'undirected')
g
#自相关也可以通过该式去除
g <- simplify(g)
#孤立节点的删除（删除度为 0 的节点）
g <- delete.vertices(g, names(degree(g)[degree(g) == 0]))

###在此处插入分类信息
# 读取分类数据
tax <- read.table("temperary_data/taxonomy_table.csv",
                  header=TRUE, row.names=1, sep=",", stringsAsFactors=FALSE)

# 确保分类数据与网络节点匹配
tax <- tax[V(g)$name, ]

# 方法1：使用set_vertex_attr函数（推荐）
for (level in c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")) {
  g <- set_vertex_attr(g, name=level, value=tax[,level])
}
# 方法2：使用vertex_attr<-语法
#vertex_attr(g) <- as.list(tax)  # 一次性添加所有分类列

# 验证添加结果
vertex_attr_names(g)  # 查看所有顶点属性
head(vertex_attr(g))  # 查看前几个节点的属性
# ------------------------------ #

#该模式下，边权重代表了相关系数
#由于权重通常为正值，因此最好取个绝对值，相关系数重新复制一列
E(g)$correlation <- E(g)$weight
E(g)$weight <- abs(E(g)$weight)
E(g)$cor[E(g)$correlation>0] <- 1
E(g)$cor[E(g)$correlation<0] <- -1

# 导入数据并整合在一起，使用到phyloseq包
#查看网络图
g
plot(g)
#graphml 格式，可使用 gephi 软件打开并进行可视化编辑
write_graph(g, 'output/bacteria_LH_网络图.graphml', format = 'graphml')
#write_graph(g, 'network.graphml', format = 'graphml')






# 关键物种分析

# 确保矩阵对称（处理模拟数据的非对称问题）
z <- (z + t(z)) / 2  # 强制对称
diag(z) <- 0         # 对角线设为0
# 检查是否有 NA
sum(is.na(z))
z[is.na(z)] <- 0  # NA 设为 0
#以上代码都是为了弥补模拟数据带来的问题的，真实数据需要去掉这些

# 关键物种分析
z[abs(z)>0]=1
z

adjacency_unweight <- z
#这是一个微生物互作网络，数值“1”表示微生物 OTU 之间存在互作，“0”表示无互作
head(adjacency_unweight)[1:6]    #邻接矩阵类型的网络文件

#邻接矩阵 -> igraph 的邻接列表，获得非含权的无向网络
igraph <- graph_from_adjacency_matrix(as.matrix(adjacency_unweight), mode = 'undirected', weighted = NULL, diag = FALSE)
igraph    #igraph 的邻接列表

#计算节点度
V(igraph)$degree <- degree(igraph)

#模块划分，详情 ?cluster_fast_greedy，有多种模型
set.seed(123)
V(igraph)$modularity <- membership(cluster_fast_greedy(igraph))

#输出各节点（微生物 OTU）名称、节点度、及其所划分的模块的列表
nodes_list <- data.frame(
  nodes_id = V(igraph)$name, 
  degree = V(igraph)$degree, 
  modularity = V(igraph)$modularity
)
head(nodes_list)    #节点列表，包含节点名称、节点度、及其所划分的模块
write.table(nodes_list, 'temperary_data/bacteria_LH_nodes_list.txt', sep = '\t', row.names = FALSE, quote = FALSE)

##计算模块内连通度（Zi）和模块间连通度（Pi）
source("script/zi_pi.r")

#上述的邻接矩阵类型的网络文件
adjacency_unweight 

#节点属性列表，包含节点所划分的模块
nodes_list <- read.delim('temperary_data/bacteria_LH_nodes_list.txt', row.names = 1, sep = '\t', check.names = FALSE)

#两个文件的节点顺序要一致
nodes_list <- nodes_list[rownames(adjacency_unweight), ]

#计算模块内连通度（Zi）和模块间连通度（Pi）
#指定邻接矩阵、节点列表、节点列表中节点度和模块度的列名称
zi_pi <- zi.pi(nodes_list, adjacency_unweight, degree = 'degree', modularity_class = 'modularity')
head(zi_pi)

##可再根据阈值对节点划分为 4 种类型，并作图展示其分布
zi_pi <- na.omit(zi_pi)   #NA 值最好去掉，不要当 0 处理
zi_pi[which(zi_pi$within_module_connectivities < 2.5 & zi_pi$among_module_connectivities < 0.62),'type'] <- 'Peripherals'
zi_pi[which(zi_pi$within_module_connectivities < 2.5 & zi_pi$among_module_connectivities > 0.62),'type'] <- 'Connectors'
zi_pi[which(zi_pi$within_module_connectivities > 2.5 & zi_pi$among_module_connectivities < 0.62),'type'] <- 'Module hubs'
zi_pi[which(zi_pi$within_module_connectivities > 2.5 & zi_pi$among_module_connectivities > 0.62),'type'] <- 'Network hubs'
write.csv(zi_pi,"bacteria_LH_zipi结果.csv")

ggplot(zi_pi, aes(among_module_connectivities, within_module_connectivities)) +
  geom_point(aes(color = type), alpha = 0.8, size = 6, shape = 16) +
  scale_y_continuous(breaks = c(-4, -2, 0, 2, 4), limits = c(-4, 4)) +
  scale_color_manual(
    values = c("#8491B4FF", "#91D1C2FF", "#F39B7FFF", "#4DBBD5FF"), 
    limits = c('Peripherals', 'Connectors', 'Module hubs', 'Network hubs')
  ) +
  theme_bw() +
  theme(
    text = element_text(family = "serif", face = "bold"),  # 全局字体设置为 serif 并加粗
    axis.text = element_text(size = 16, color = "black"), # 设置坐标轴字体大小和颜色
    axis.title = element_text(size = 20, color = "black"), # 设置坐标轴标题大小和颜色
    axis.ticks = element_line(linetype = "solid", linewidth = 1, color = "black"), # 设置刻度线粗细和颜色
    panel.border = element_rect(linetype = "solid", color = "black", linewidth = 2), # 设置图形边框
    panel.grid.major = element_line(linetype = "solid", color = "gray90", linewidth = 0.75), # 设置主要网格线
    panel.grid.minor = element_line(linetype = "solid", color = "gray95", linewidth = 0.5), # 设置次要网格线
    legend.text = element_text(size = 16, color = "black"),  # 设置图例文本
    legend.title = element_text(size = 16, color = "black")  # 设置图例标题
  ) +
  labs(x = 'Among-module connectivities', y = 'Within-module connectivities', color = '') +
  geom_vline(xintercept = 0.62, linetype = 2, size = 1) +
  geom_hline(yintercept = 2.5, linetype = 2, size = 1)


# 将图片以.pdf的格式导出
ggsave("output/bacteria_LH_关键物种分析.pdf", width = 8, height = 6, units = "in")

