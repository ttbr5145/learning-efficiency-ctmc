###############################################################################
# 01_recon_data.R —— 标定数据重建
# 原文未公开原始问卷，故以原文报告的统计特征（519 份、5 类学生、效率分布、
# 表3 数值、S480 案例等）为锚点重建数据；随机性由固定种子控制，可完全复现。
###############################################################################
set.seed(20261121)
# 仓库根目录：优先环境变量 PROJ_ROOT，其次当前工作目录（应在仓库根），最后回退本机默认路径
root <- Sys.getenv("PROJ_ROOT", unset = getwd())
if (!file.exists(file.path(root, "代码", "dea_fun.R"))) root <- "D:/AI/随机过程大作业2"
setwd(root)
source(file.path(root, "代码", "dea_fun.R"))

n <- 519
item_names <- c("interest", "revision", "examtimes", "onspot", "mental", "rest",
                "tutoring", "mocktest", "reciting", "listening", "speaking",
                "reading", "extracur", "abroad")
item_cn <- c("英语学习兴趣", "课后复习重视程度", "参加考试次数", "考场临场发挥",
             "考前心理状态", "休息充分程度", "课后辅导次数", "模拟考试次数",
             "朗读背诵时间", "听力训练时间", "口语练习时间", "阅读时间",
             "课外积累", "出国意愿")

# ---- 1. 五类学生的潜在画像（0~100 标度）---------------------------------
prof <- rbind(
  C1 = c(22, 16, 34, 28, 34, 52, 14, 24, 26, 22, 18, 22, 16, 22),  # 被动型
  C2 = c(46, 32, 82, 74, 66, 46, 40, 78, 38, 38, 32, 38, 30, 38),  # 潜力型
  C3 = c(56, 60, 84, 62, 58, 48, 58, 76, 62, 54, 48, 58, 52, 44),  # 积累型
  C4 = c(68, 70, 86, 72, 66, 54, 66, 84, 70, 66, 56, 68, 62, 56),  # 积累型
  C5 = c(88, 88, 70, 86, 86, 74, 50, 62, 88, 84, 78, 86, 84, 68)   # 均衡型
)
sizes <- c(80, 95, 96, 105, 143)
stopifnot(sum(sizes) == n)
cluster_true <- rep(1:5, times = sizes)

# ---- 2. 生成三种题型的原始作答 ------------------------------------------
# 题型 A：5 级有序选项（态度类）；题型 B：0~12 次计数（频次类）；
# 题型 C：0~20 小时的连续时长（时间类）
type <- c("A", "A", "B", "A", "A", "A", "B", "B", "C", "C", "C", "C", "A", "A")
lat <- matrix(NA_real_, n, length(item_names))
for (j in seq_along(item_names)) {
  mu <- (prof[, j] - 50) / 11
  lat[, j] <- mu[cluster_true] + rnorm(n, 0, 0.80)
}
ordinal <- raw_counts <- hours <- matrix(NA_real_, n, length(item_names))
thr <- c(-1.2816, -0.5244, 0.5244, 1.2816)
for (j in seq_along(item_names)) {
  if (type[j] == "A") {
    ordinal[, j] <- 6L - as.integer(cut(lat[, j], c(-Inf, thr, Inf), labels = FALSE))
  } else if (type[j] == "B") {
    raw_counts[, j] <- pmin(round(pmax(2 + 2.2 * lat[, j] + rnorm(n, 0, 1), 0)), 12)
  } else {
    hours[, j] <- pmin(pmax(4 + 3.2 * lat[, j] + rnorm(n, 0, 2), 0), 20)
  }
}

# ---- 3. 复合量化：累积百分比法 -------------------------------------------
# 有序选项：s_k = 100 * (优于或等于 k 档的累计人数比例)，最优档为 100；
# 计数/时长：按经验分布函数（即累积百分比）映射至 0~100，并作对数压缩以
# 抑制长尾。三种题型统一到同一标度，供聚类与 DEA 使用。
item_score <- matrix(NA_real_, n, length(item_names))
quant_example <- NULL
for (j in seq_along(item_names)) {
  if (type[j] == "A") {
    tab <- as.numeric(table(factor(ordinal[, j], levels = 1:5)) / n)
    # 累积百分比法（分位中点口径）：越优的选项分越高
    s <- numeric(5); acc <- 0
    for (k in 1:5) {
      s[k] <- 100 * (1 - (acc + 0.5 * tab[k]))
      acc <- acc + tab[k]
    }
    item_score[, j] <- s[ordinal[, j]]
    if (j == 1) quant_example <- data.frame(level = 1:5, share = round(tab * 100, 2),
                                            score = round(s, 2))
  } else if (type[j] == "B") {
    v <- raw_counts[, j]
    item_score[, j] <- 100 * (rank(v, ties.method = "average") - 0.5) / n
  } else {
    v <- log1p(hours[, j])
    item_score[, j] <- 100 * (rank(v, ties.method = "average") - 0.5) / n
  }
  if (var(item_score[, j]) < 10) item_score[, j] <- 100 * (rank(lat[, j]) - 0.5) / n
}
colnames(item_score) <- item_names

# ---- 4. 投入-产出变量 ----------------------------------------------------
x1 <- rowMeans(item_score[, c("reciting", "listening", "speaking", "reading")])
x2 <- rowMeans(item_score[, c("examtimes", "mocktest")])
x3 <- rowMeans(item_score[, c("revision", "tutoring", "extracur")])
x4 <- rowMeans(item_score[, c("mental", "rest")])

base <- 45 + 0.52 * (0.45 * x1 + 0.26 * x2 + 0.18 * x3 + 0.11 * x4)
c_eff <- c(-11, -4, 0, 4, 9)[cluster_true]
# 期末成绩 = 投入可解释部分 + 个体“转化能力”残差（个体异质性）
e_i <- rnorm(n, c_eff, 11)
score <- pmin(pmax(base + e_i, 30), 99)

# ---- 5. 校准个案 S480：类内 DEA 效率命中 0.6771 ---------------------------
idx480 <- 480
stopifnot(cluster_true[idx480] == 5)
X <- cbind(x1, x2, x3, x4); Y <- matrix(score, ncol = 1)
S5 <- which(cluster_true == 5)
theta_of <- function(e480) {
  sc <- score; sc[idx480] <- base[idx480] + e480
  dea_ccr(X, matrix(sc, ncol = 1), idx480, S = S5)$theta
}
lo <- -40; hi <- 30
for (k in 1:50) {
  mid <- (lo + hi) / 2
  if (theta_of(mid) > 0.6771) hi <- mid else lo <- mid
}
e_i[idx480] <- (lo + hi) / 2
score[idx480] <- base[idx480] + e_i[idx480]
cat("校准后 S480 成绩 =", round(score[idx480], 2), "\n")

dat <- data.frame(id = 1:n, cluster_true = cluster_true, item_score,
                  x1 = x1, x2 = x2, x3 = x3, x4 = x4,
                  score = round(score, 2), e = e_i)
write.csv(dat, file.path(root, "数据", "learning_data.csv"), row.names = FALSE)
write.csv(quant_example, file.path(root, "数据", "quant_example.csv"), row.names = FALSE)

cat("n =", n, " 各类：", paste(tabulate(cluster_true), collapse = " / "), "\n")
cat("成绩均值", round(mean(score), 2), " 范围", round(range(score), 2),
    " S480 =", round(score[idx480], 2), "\n")

