
#relative2为注释信息（无误）+相对丰度的结果
# 导入相对丰度表
asv_raw <- read.table("relative_table_tax.txt", header = T, row.names = 1, sep = "\t")

# 先处理门水平
# 提取phylum
phylum <- asv_raw[, c(3, 9:53)]
# 将phylum按照相同的Phylum进行合并
phylum1 <- aggregate(. ~ Phylum, data = phylum, sum)
phylum1 <- aggregate(phylum[, -1], by = list(Phylum = phylum$Phylum), FUN = sum)
# 输出处理后的表格
write.table(phylum1, "output/phylum_rel.txt", sep = "\t", quote = F, row.names = F)

# 再处理纲水平
#class <- asv_raw[, c(3, 8:67)]
# 将class按照相同的Class进行合并
#class1 <- aggregate(. ~ Class, data = class, sum)
#write.table(class1, "C:/Users/user/Desktop/小论文数据整理与分析/真菌/3不同分类水平物种数量统计/class_rel.txt", sep = "\t", quote = F, row.names = F)

#目
#order <- asv_raw[, c(4, 8:67)]
#family <- asv_raw[, c(5, 8:67)]
#family1 <- aggregate(. ~ Family, data = family, sum)
#write.table(family1, "C:/Users/user/Desktop/小论文数据整理与分析/真菌/3不同分类水平物种数量统计/family_rel.txt", sep = "\t", quote = F, row.names = F)

#属
genus <- asv_raw[, c(7, 9:53)]
genus1 <- aggregate(. ~ Genus, data = genus, sum)
genus1 <- aggregate(genus[, -1], by = list(Genus = genus$Genus), FUN = sum)
write.table(genus1, "output/genus_rel.txt", sep = "\t", quote = F, row.names = F)


