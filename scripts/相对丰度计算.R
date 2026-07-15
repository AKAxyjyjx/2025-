
asv4 <- read.table("temperary_data/asv_table_only.csv", header=T, sep=",", stringsAsFactors = FALSE) 
asv4 <- rbind(asv4, colSums(asv4))
for (i in 1:ncol(asv4)) {
  asv4[, i] <- asv4[, i] / asv4[nrow(asv4), i]
}

# 去除最后一行
asv4 <- asv4[-nrow(asv4), ]
# 检验
sum <- colSums(asv4)
sum

# 将相对丰度表导出，导出后记得打开表格手动调整
write.table(asv4, "temperary_data/relative_table_only.csv", sep = "\t", quote = F)

#生成的relative表格第一行也有错位，需要调整
