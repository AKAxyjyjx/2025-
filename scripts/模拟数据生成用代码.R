#生成ASV表

generate_asv_table <- function(
    n_samples = 45,
    n_asvs = 5002,
    group_vector = NULL,
    group_names = NULL,
    samples_per_group = NULL,
    seed = 2026,
    include_taxonomy = TRUE,
    include_sequences = TRUE
) {
  
  set.seed(seed)
  
  # ============================================
  # 1. 处理分组信息
  # ============================================
  
  if (!is.null(group_vector)) {
    if (length(group_vector) != n_samples) {
      stop("group_vector长度必须等于n_samples")
    }
    groups <- factor(group_vector)
    n_groups <- length(unique(groups))
    samples_per_group <- table(groups)
    group_levels <- levels(groups)
    
  } else if (!is.null(samples_per_group) && !is.null(group_names)) {
    if (length(samples_per_group) != length(group_names)) {
      stop("samples_per_group长度必须等于group_names长度")
    }
    n_groups <- length(group_names)
    groups <- rep(group_names, times = samples_per_group)
    if (length(groups) != n_samples) {
      stop("总样本数不等于n_samples")
    }
    groups <- factor(groups)
    group_levels <- levels(groups)
    
  } else if (!is.null(group_names)) {
    n_groups <- length(group_names)
    samples_per_group <- rep(floor(n_samples / n_groups), n_groups)
    remainder <- n_samples - sum(samples_per_group)
    if (remainder > 0) {
      samples_per_group[1:remainder] <- samples_per_group[1:remainder] + 1
    }
    groups <- rep(group_names, times = samples_per_group)
    groups <- factor(groups)
    group_levels <- levels(groups)
    
  } else {
    n_groups <- 3
    group_names <- c("Control", "Treatment", "Other")
    samples_per_group <- rep(floor(n_samples / n_groups), n_groups)
    remainder <- n_samples - sum(samples_per_group)
    if (remainder > 0) {
      samples_per_group[1:remainder] <- samples_per_group[1:remainder] + 1
    }
    groups <- rep(group_names, times = samples_per_group)
    groups <- factor(groups)
    group_levels <- levels(groups)
    cat("使用默认分组：", paste(group_names, collapse = ", "), "\n")
  }
  
  n_groups <- length(unique(groups))
  group_levels <- levels(groups)
  
  cat("分组信息：\n")
  print(table(groups))
  
  # ============================================
  # 2. 生成样本名（关键修改在这里！）
  # ============================================
  
  sample_names <- character(n_samples)
  for (g in 1:n_groups) {
    group_indices <- which(groups == group_levels[g])
    for (i in seq_along(group_indices)) {
      sample_names[group_indices[i]] <- paste0(group_levels[g], i)
    }
  }
  
  asv_names <- paste0("ASVID", 1:n_asvs)
  
  cat("样本名示例：", paste(head(sample_names), collapse = ", "), "\n")
  
  # ============================================
  # 3. 生成ASV计数表
  # ============================================
  
  cat("\n生成ASV表...\n")
  
  asv_table <- matrix(0, nrow = n_samples, ncol = n_asvs)
  
  n_core <- round(n_asvs * 0.03)
  n_medium <- round(n_asvs * 0.10)
  n_rare <- n_asvs - n_core - n_medium
  
  asv_roles <- rep("rare", n_asvs)
  asv_roles[1:n_core] <- "core"
  asv_roles[(n_core + 1):(n_core + n_medium)] <- "medium"
  asv_roles <- sample(asv_roles)
  
  group_specific_n <- round(n_asvs * 0.02)
  group_specific_asvs <- list()
  all_asvs <- 1:n_asvs
  
  for (g in 1:n_groups) {
    available <- which(asv_roles != "core")
    if (length(available) > group_specific_n) {
      group_specific_asvs[[g]] <- sample(available, group_specific_n)
    } else {
      group_specific_asvs[[g]] <- sample(all_asvs, group_specific_n)
    }
  }
  
  for (i in 1:n_samples) {
    group_idx <- which(group_levels == groups[i])
    depth <- rnorm(1, mean = 15000, sd = 4000)
    depth <- max(depth, 3000)
    
    for (j in 1:n_asvs) {
      if (asv_roles[j] == "core") {
        mu <- runif(1, 50, 200)
        prob <- 0.95
      } else if (j %in% group_specific_asvs[[group_idx]]) {
        mu <- runif(1, 30, 100)
        prob <- 0.80
      } else if (asv_roles[j] == "medium") {
        mu <- runif(1, 10, 50)
        prob <- 0.50
      } else {
        mu <- runif(1, 1, 10)
        prob <- 0.20
      }
      
      if (runif(1) < prob) {
        asv_table[i, j] <- rnbinom(1, size = 0.5, mu = mu)
      }
    }
  }
  
  rownames(asv_table) <- sample_names
  colnames(asv_table) <- asv_names
  
  # ============================================
  # 4. 生成分类信息
  # ============================================
  
  cat("生成分类信息...\n")
  
  if (include_taxonomy) {
    phyla_list <- c(
      "Proteobacteria", "Firmicutes", "Bacteroidetes", 
      "Actinobacteria", "Acidobacteria", "Chloroflexi",
      "Cyanobacteria", "Planctomycetes", "Verrucomicrobia",
      "Gemmatimonadetes", "Nitrospirae"
    )
    
    genera_list <- c(
      "Escherichia", "Bacillus", "Lactobacillus", "Streptococcus",
      "Staphylococcus", "Pseudomonas", "Bacteroides", "Clostridium",
      "Acinetobacter", "Klebsiella", "Enterococcus", "Listeria",
      "Corynebacterium", "Mycobacterium", "Streptomyces", "Salmonella",
      "Shigella", "Vibrio", "Yersinia", "Legionella",
      "Campylobacter", "Helicobacter", "Neisseria", "Haemophilus",
      "Prevotella", "Ruminococcus", "Faecalibacterium", "Bifidobacterium"
    )
    
    taxonomy <- data.frame(
      ASV_ID = asv_names,
      Kingdom = rep("Bacteria", n_asvs),
      Phylum = sample(phyla_list, n_asvs, replace = TRUE),
      Class = character(n_asvs),
      Order = character(n_asvs),
      Family = character(n_asvs),
      Genus = character(n_asvs),
      Species = character(n_asvs),
      stringsAsFactors = FALSE
    )
    
    for (i in 1:n_asvs) {
      phylum <- taxonomy$Phylum[i]
      class_options <- switch(phylum,
                              "Proteobacteria" = c("Alphaproteobacteria", "Betaproteobacteria", 
                                                   "Gammaproteobacteria", "Deltaproteobacteria"),
                              "Firmicutes" = c("Bacilli", "Clostridia", "Erysipelotrichia"),
                              "Bacteroidetes" = c("Bacteroidia", "Flavobacteriia", "Sphingobacteriia"),
                              "Actinobacteria" = c("Actinobacteria", "Coriobacteriia"),
                              "Acidobacteria" = c("Acidobacteriia", "Blastocatellia"),
                              "Cyanobacteria" = c("Cyanophyceae", "Oscillatoriophycideae"),
                              "Chloroflexi" = c("Chloroflexia", "Anaerolineae"),
                              "Planctomycetes" = c("Planctomycetacia", "Phycisphaerae"),
                              "Verrucomicrobia" = c("Verrucomicrobiae", "Opitutae"),
                              "Gemmatimonadetes" = "Gemmatimonadetes",
                              "Nitrospirae" = "Nitrospira"
      )
      taxonomy$Class[i] <- sample(class_options, 1)
      taxonomy$Genus[i] <- sample(genera_list, 1)
      taxonomy$Species[i] <- paste0(taxonomy$Genus[i], "_sp_", sample(100:999, 1))
      taxonomy$Family[i] <- paste0(taxonomy$Genus[i], "aceae")
      taxonomy$Order[i] <- paste0(taxonomy$Genus[i], "ales")
    }
  } else {
    taxonomy <- NULL
  }
  
  # ============================================
  # 5. 生成ASV序列
  # ============================================
  
  if (include_sequences) {
    cat("生成ASV序列...\n")
    sequences <- character(n_asvs)
    for (i in 1:n_asvs) {
      seq_length <- sample(180:250, 1)
      bases <- sample(c("A", "T", "C", "G"), seq_length, replace = TRUE)
      sequences[i] <- paste(bases, collapse = "")
    }
    names(sequences) <- asv_names
  } else {
    sequences <- NULL
  }
  
  # ============================================
  # 6. 创建分组信息
  # ============================================
  
  group_info <- data.frame(
    Sample = sample_names,
    Group = as.character(groups),
    stringsAsFactors = FALSE
  )
  rownames(group_info) <- sample_names
  
  # ============================================
  # 7. 统计信息
  # ============================================
  
  statistics <- list(
    n_samples = n_samples,
    n_asvs = n_asvs,
    n_groups = n_groups,
    group_names = group_levels,
    samples_per_group = as.vector(table(groups)),
    sparsity = sum(asv_table == 0) / (n_samples * n_asvs),
    mean_depth = mean(rowSums(asv_table)),
    sd_depth = sd(rowSums(asv_table)),
    seed = seed,
    generated_date = Sys.time()
  )
  
  cat("✓ 生成完成！\n")
  cat("样本数:", n_samples, "\n")
  cat("ASV数:", n_asvs, "\n")
  cat("分组数:", n_groups, "\n")
  cat("稀疏度:", round(statistics$sparsity * 100, 2), "%\n")
  
  return(list(
    asv_table = asv_table,
    taxonomy = taxonomy,
    sequences = sequences,
    group_info = group_info,
    statistics = statistics
  ))
}




data <- generate_asv_table(
  n_samples = 45,
  n_asvs = 5002,
  group_names = c("EL", "EM", "EH", 'ML', 'MM', 'MH', 'LL', 'LM', 'LH'),
  samples_per_group = c(5, 5, 5, 5, 5, 5, 5, 5, 5)
)






save_asv_data <- function(
    data,                          # generate_asv_table() 返回的数据
    output_dir = "temperary_data/"   # 输出目录
) 
  
  # 1. 创建输出目录
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # ============================================
  # 2. 保存分类信息
  # ============================================
  
  
  if (!is.null(data$taxonomy)) {
  taxonomy_table <- data$taxonomy
  tax_cols <- c("ASV_ID", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  tax_cols_exist <- tax_cols[tax_cols %in% names(taxonomy_table)]
  taxonomy_table <- taxonomy_table[, tax_cols_exist, drop = FALSE]
  }
  
  
  #write.csv(taxonomy_table, 
            file = paste0(output_dir, "taxonomy_table.csv"), 
            row.names = FALSE)
  write.csv(taxonomy_table, 
            "temperary_data/taxonomy_table.csv", 
            row.names = FALSE)
  # ============================================
  # 保存不带分类的纯ASV表（样本为行）
  # ============================================
  
  
  asv_only <- t(as.data.frame(data$asv_table))
  
  
  write.csv(asv_only, 
            file = paste0(output_dir, "asv_table_only.csv"), 
            row.names = FALSE)
  
  # ============================================
  # 4. 保存分组信息
  # ============================================
  
  
  # 每个样本的分组
  df <- data$group_info[, c(2, 1)]
  #write.csv(df, 
            file = paste0(output_dir, "sample_group_info.csv"), 
            row.names = FALSE)
  
  
  write.csv(df, 
            "temperary_data/sample_group_info.csv", 
            row.names = FALSE)

  
#生成环境信息
  # ============================================
  # 生成环境因素表
  # ============================================
  
generate_environment_data <- function(
    sample_names = paste0("EL", 1:3),  # 样本名
    seed = 2026
  ) {
    set.seed(seed)
    
    n <- length(sample_names)
    
    # 生成环境因子数据
    data <- data.frame(
      Sample = sample_names,
      
      # 植物组织 Cd 含量 (mg/kg)
      `Leaf-Cd` = round(runif(n, 0, 1.2), 5),
      `Stem-Cd` = round(runif(n, 0.5, 3.5), 5),
      
      # Cd 形态分级 (mg/kg)
      Residual = round(runif(n, 0, 0.15), 5),
      Oxidizable = round(runif(n, 0, 0.08), 5),
      Reducible = round(runif(n, 0, 0.08), 5),
      `Acid-extractable` = round(runif(n, 0, 0.08), 5),
      
      # 土壤理化性质
      Eh = round(runif(n, 20, 430), 2),
      pH = round(runif(n, 5, 7), 2)
    )
    
    return(data)
  }
  
# ============================================
# 使用示例
# ============================================

  # 示例3：生成所有样本（45个）
all_samples <- c(
    paste0("EL", 1:5),
    paste0("EM", 1:5),
    paste0("EH", 1:5),
    paste0("ML", 1:5),
    paste0("MM", 1:5),
    paste0("MH", 1:5),
    paste0("LL", 1:5),
    paste0("LM", 1:5),
    paste0("LH", 1:5)
  )
env_data <- generate_environment_data(sample_names = all_samples, seed = 123)
head(env_data)
  
  # 保存
write.csv(env_data, "temperary_data/env_data.csv", row.names = FALSE)
