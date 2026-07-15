
# 读取
asv_table <- read.csv("relative_table_only.csv")
taxonomy_table <- read.csv("taxonomy_table.csv")

#asv_table <-asv4
#taxonomy_table 
#asv4_rel <- as.data.frame(asv4)
#asv4_df <- cbind(ASV_ID = rownames(asv4_rel), asv4_rel)
#rownames(asv4_df) <- NULL
# 合并
merged <- left_join(taxonomy_table, asv_table, by = "ASV_ID")
#merged <- left_join(taxonomy_table, asv4_df, by = "ASV_ID")
# 查看
head(merged[, 1:10])

write.table(merged, "relative_table_tax.csv", sep = "\t", quote = F)