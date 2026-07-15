# ============================================
# 水香薷-微生物-重金属联合分析主程序
# 作者：陈逸阳
# 最后更新：2026-07-09
# ============================================

#  R 包
library(tidyverse)
library(ggplot2)
library(vegan)
#install.packages("agricolae")
library(agricolae)
library(ggpubr)
library(MetBrewer)
#install.packages("picante")
library(picante)
library(ape)
#install.packages("ggforce")
library(ggforce)
library(dplyr)
#install.packages("ggrepel")
library(ggrepel)
#install.packages("wesanderson")
library(wesanderson)
library(pheatmap)
library(scales)
#install.packages("reshape2")
library(reshape2)
#install.packages("devtools")
library(devtools)
install_github('Hy4m/linkET', force = TRUE)   # 安装包
#devtools::install_github("Hy4m/linkET", force = TRUE)
library(linkET)
##install.packages("openxlsx")
#library(openxlsx)
#install.packages("igraph")
library(igraph)
#install.packages("psych")
library(psych)
library(patchwork)
library(doBy)
library(FD)

# 读取数据
# 原始数据因保密协议不公开，此处为模拟数据结构
source("scripts/模拟数据生成用代码.R")  # 如果你把数据读取单独放的话
source("scripts/相对丰度计算.R")  
source("scripts/统计合并作业.R")  


# Alpha多样性
source("scripts/Alpha多样性.R")

# Alpha稀释曲线
source("scripts/改了又改终于可以用的alpha稀释曲线.R")

# Beta多样性
source("scripts/Beta多样性.R")

#在生成数据之前先在丰度基础上做不同分类水平物种统计
source("scripts/不同分类水平物种统计（丰度基础上）.R")

# 门水平堆积柱状图
source("scripts/门水平堆积柱状图.R")

# 属水平热图
source("scripts/属水平热图.R")


# 每个重金属元素的Mantel test
source("scripts/Mantel分析.R")

# 每个重金属元素的关键物种网络图
source("scripts/网络图+关键物种分析.R")

# RDA
source("scripts/RDA分析.R")

# ============================================
# ~咱也只做了这么多desuwa
# ============================================