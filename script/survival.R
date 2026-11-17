

setwd("C:\\Users\\17860\\Desktop\\FRG\\28.External")

# 加载必要的包
library(survival)
library(survminer)

# 读取基因表达矩阵文件
gene_matrix <- read.table("geneMatrix.txt", header = TRUE, row.names = 1, sep = "\t")

# 读取生存时间和生存状态文件
survival_data <- read.table("time.txt", header = TRUE, sep = "\t")

# 检查基因表达矩阵和生存数据的前几行
head(gene_matrix)
head(survival_data)

# 确保生存数据中的样本ID列名与基因表达矩阵的行名匹配
# 例如，如果生存数据中的样本ID列名为 "id"，我们需要重命名它以匹配基因表达矩阵的行名
colnames(survival_data)[colnames(survival_data) == "id"] <- "Sample_ID"

# 设置样本ID为行名
rownames(survival_data) <- survival_data$Sample_ID

# 筛选预后基因的表达数据
prognostic_genes <- c("XDH", "MMACHC", "BCL2", "PTK2", "PIK3CG", "ATF4", "LCK","IL18","RUNX2")
filtered_gene_matrix <- gene_matrix[rownames(gene_matrix) %in% prognostic_genes, ]

# 检查筛选后的基因表达矩阵
head(filtered_gene_matrix)

# 转置基因表达矩阵以便样本作为行
filtered_gene_matrix_t <- t(filtered_gene_matrix)

# 打印转置后的基因表达矩阵的前几行
head(filtered_gene_matrix_t)

# 合并基因表达矩阵和生存数据，确保样本匹配
merged_data <- merge(survival_data, filtered_gene_matrix_t, by = "row.names", all = TRUE)
colnames(merged_data)[1] <- "Sample"

# 检查合并后的数据框
head(merged_data)

# 转换时间和状态列为数值型数据
merged_data$time <- as.numeric(merged_data$time)
merged_data$status <- as.numeric(merged_data$status)

# 检查时间和状态列是否有缺失值或非数值数据
sum(is.na(merged_data$time))
sum(is.na(merged_data$status))

# 去除包含缺失值的样本
merged_data <- merged_data[complete.cases(merged_data$time, merged_data$status), ]

# 计算预后打分
merged_data$PI <- (0.296028305971066 * merged_data$XDH) + 
  (0.26340117648018 * merged_data$MMACHC) + 
  (0.190492580140298 * merged_data$BCL2) -  
  (0.172031701467542 * merged_data$PTK2) - 
  (0.155098091497332 * merged_data$PIK3CG) + 
  (0.316484640419813 * merged_data$ATF4) - 
  (0.222412916132292 * merged_data$LCK) -
  (0.250971333082422 * merged_data$IL18) - 
  (0.284812710331733 * merged_data$RUNX2)
# 检查预后打分是否有缺失值
sum(is.na(merged_data$PI))

# 去除包含缺失值的样本
merged_data <- merged_data[complete.cases(merged_data$PI), ]

# 按预后打分中位数将样本分为高低风险组
median_PI <- median(merged_data$PI, na.rm = TRUE)
merged_data$risk_group <- ifelse(merged_data$PI >= median_PI, "high", "low")

# 创建生存对象
surv_object <- Surv(time = merged_data$time, event = merged_data$status)

# 检查生存对象和风险组变量是否有缺失值
sum(is.na(surv_object))
sum(is.na(merged_data$risk_group))

# 确保生存对象和风险组变量没有缺失值
merged_data <- merged_data[complete.cases(merged_data$time, merged_data$status, merged_data$risk_group), ]

# 再次创建生存对象
surv_object <- Surv(time = merged_data$time, event = merged_data$status)

# 生存分析
fit <- survfit(surv_object ~ risk_group, data = merged_data)

# 绘制生存曲线
ggsurvplot(fit, 
           data = merged_data, 
           risk.table = TRUE, 
           pval = TRUE, 
           conf.int = TRUE,
           xlab = "Time",
           ylab = "Survival Probability",
           title = "Survival Analysis by Risk Group")

# 保存图像
ggsave("survival_plot.png")
