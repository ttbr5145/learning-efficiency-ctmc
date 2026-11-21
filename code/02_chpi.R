###############################################################################
# 02_chpi.R —— CHPI 模型复现
#   1) Gower + PAM 聚类（k=5，与原文类型数一致）
#   2) 类内 CCR-DEA 效率
#   3) context-dependent DEA 分层（剔除式 + 分位带两种口径）
#   4) 吸引力 Ω、进步度 P、邻近度 D、障碍度 H
#   5) 参考路径选择与 S480 个案
###############################################################################
root <- "C:/Users/Lenovo/Desktop/随机过程大作业_学习策略推荐"
setwd(root)
source(file.path(root, "code", "dea_fun.R"))
suppressPackageStartupMessages({library(cluster); library(showtext); library(ggplot2); library(lpSolve)})
showtext_auto()
font_add("SimSun", "C:/Windows/Fonts/simsun.ttc")
font_add("SimHei", "C:/Windows/Fonts/simhei.ttf")

dat <- read.csv(file.path(root, "data", "learning_data.csv"))
item_names <- c("interest", "revision", "examtimes", "onspot", "mental", "rest",
                "tutoring", "mocktest", "reciting", "listening", "speaking",
                "reading", "extracur", "abroad")
X <- as.matrix(dat[, c("x1", "x2", "x3", "x4")])
Y <- matrix(dat$score, ncol = 1)
n <- nrow(dat)

# ---- 1. 聚类 -------------------------------------------------------------
gower <- daisy(dat[, item_names], metric = "gower")
sil <- sapply(3:8, function(k) pam(gower, k = k, diss = TRUE)$silinfo$avg.width)
sil_tab <- data.frame(k = 3:8, silhouette = round(sil, 4))
write.csv(sil_tab, file.path(root, "data", "silhouette.csv"), row.names = FALSE)

k_best <- 5
set.seed(20261121)
pamres <- pam(gower, k = k_best, diss = TRUE)
clus0 <- pamres$clustering
# 依与设计类型的重合度，把 PAM 类别重新编号为 1..5（1=被动型 … 5=均衡型）
conf <- table(clus0, dat$cluster_true)
assign_mat <- conf
sol <- lp.assign(-assign_mat)$solution
newlab <- apply(sol, 1, which.max)            # PAM 类 -> 类型编号
clus <- as.integer(newlab[clus0])
pam_agree <- mean(clus == dat$cluster_true)

# ---- 2. 特征画像 --------------------------------------------------------
item_cn <- c("英语学习兴趣", "课后复习重视", "参加考试次数", "考场临场发挥",
             "考前心理状态", "休息充分程度", "课后辅导次数", "模拟考试次数",
             "朗读背诵时间", "听力训练时间", "口语练习时间", "阅读时间",
             "课外积累", "出国意愿")
prof_tab <- aggregate(dat[, item_names], by = list(cluster = clus), mean)
prof_tab <- cbind(prof_tab, score = tapply(dat$score, clus, mean),
                  n = as.numeric(table(clus)))
for (cc in item_names) prof_tab[[cc]] <- round(as.numeric(prof_tab[[cc]]), 1)
prof_tab$score <- round(prof_tab$score, 1)

# ---- 3. 类内 CCR-DEA 效率 -----------------------------------------------
theta <- numeric(n); u_list <- vector("list", n); v_list <- vector("list", n)
for (c in 1:k_best) {
  S <- which(clus == c)
  for (j in S) {
    r <- dea_ccr(X, Y, j, S = S)
    theta[j] <- r$theta; u_list[[j]] <- r$u; v_list[[j]] <- r$v
  }
}
dat$theta <- theta
r_eff_score <- cor(dat$theta, dat$score)
fit <- lm(score ~ theta, data = dat)

# ---- 4. 分层（两种口径）------------------------------------------------
# (a) 原文 Algorithm 2 的剔除式分层
level_peel <- rep(NA_integer_, n)
for (c in 1:k_best) {
  S <- which(clus == c); rem <- S
  for (l in 1:5) {
    if (l == 5) { level_peel[rem] <- 5; break }
    th <- sapply(rem, function(j) dea_ccr(X, Y, j, S = rem)$theta)
    eff <- rem[th >= 0.9995]
    if (length(eff) == 0) eff <- rem[which.max(th)]
    level_peel[eff] <- l
    rem <- setdiff(rem, eff)
  }
}
# (b) 类内效率分位带分层（主口径，层级规模与原文图7 的形态一致）
level <- rep(NA_integer_, n)
for (c in 1:k_best) {
  S <- which(clus == c)
  q <- quantile(theta[S], c(0.90, 0.75, 0.60, 0.45))
  level[S] <- ifelse(theta[S] >= q[1], 1,
              ifelse(theta[S] >= q[2], 2,
              ifelse(theta[S] >= q[3], 3,
              ifelse(theta[S] >= q[4], 4, 5))))
}
dat$level <- level
layer_tab <- as.data.frame(table(cluster = dat$clus, level = dat$level))
layer_tab <- reshape(layer_tab, idvar = "cluster", timevar = "level", direction = "wide")
names(layer_tab) <- c("cluster", paste0("L", 1:5))
layer_peel_tab <- as.data.frame(table(cluster = dat$clus, level = level_peel))
layer_peel_tab <- reshape(layer_peel_tab, idvar = "cluster", timevar = "level",
                          direction = "wide")
names(layer_peel_tab) <- c("cluster", paste0("L", 1:5))

# ---- 5. Ω / P / D / H ---------------------------------------------------
omega <- rep(NA_real_, n); prog <- rep(NA_real_, n); prox <- rep(NA_real_, n)
ref_upper <- rep(NA_integer_, n)
for (c in 1:k_best) {
  Sc <- which(clus == c)
  for (j in Sc) {
    lj <- level[j]
    up <- Sc[level[Sc] < lj]; dn <- Sc[level[Sc] > lj]
    if (length(up) > 0) {
      omega[j] <- dea_output(X, Y, X[j, ], Y[j, 1], up)
      d <- rowSums(abs(sweep(X[up, , drop = FALSE], 2, X[j, ]))) + abs(Y[up, 1] - Y[j, 1])
      prox[j] <- min(d); ref_upper[j] <- up[which.min(d)]
    }
    if (length(dn) > 0) prog[j] <- dea_output(X, Y, X[j, ], Y[j, 1], dn)
  }
}
dat$omega <- omega; dat$progress <- prog; dat$prox <- prox
dat$H <- 1 / dat$omega + dat$progress
dat$ref_upper <- ref_upper

# ---- 6. 路径 ------------------------------------------------------------
s480 <- 480
path <- s480
repeat {
  j <- tail(path, 1)
  if (level[j] == 1 || is.na(ref_upper[j])) break
  up <- which(clus == clus[j] & level == level[j] - 1)
  if (length(up) == 0) break
  ord2 <- order(dat$H[up], dat$prox[up], na.last = TRUE)
  nxt <- up[ord2[1]]
  if (nxt %in% path) break
  path <- c(path, nxt)
}
path_theta <- theta[path]
path_tab <- data.frame(order = seq_along(path), id = path, level = level[path],
                       theta = round(path_theta, 4),
                       gain = round(c(NA, diff(path_theta) / head(path_theta, -1) * 100), 2),
                       cum_gain = round((path_theta - path_theta[1]) / path_theta[1] * 100, 2),
                       H = round(dat$H[path], 4), prox = round(dat$prox[path], 2))

# ---- 7. 第5类第4层个体明细（对应原文表3）--------------------------------
c5 <- which(clus == 5)
l4 <- c5[level[c5] == 4]
tab4 <- data.frame(No = l4, eff = round(theta[l4], 4), Omega = round(omega[l4], 4),
                   P = round(prog[l4], 4), D = round(prox[l4], 2), H = round(dat$H[l4], 4))
tab4 <- tab4[order(tab4$H), ]

# ---- 8. 复现锚点对照 ----------------------------------------------------
repro <- data.frame(
  quantity = c("样本量", "学生类型数", "聚类轮廓宽度(k=5)", "效率均值", "效率最小值",
               "效率-成绩相关系数", "第5类人数", "第5类层数(L5人数)",
               "S480 效率", "S480 起始层级", "S480 单步提升(%)", "S480 累计提升(%)"),
  original = c("519", "5", "—", "约0.85", "约0.50", "0.87", "最多(143)", "5",
               "0.6771", "L5", "25.38", "43.97"),
  replica  = c(as.character(n), as.character(k_best), sprintf("%.3f", sil[3]),
               sprintf("%.4f", mean(theta)), sprintf("%.4f", min(theta)),
               sprintf("%.3f", r_eff_score), as.character(sum(clus == 5)),
               paste0(max(level[c5]), "(", sum(level[c5] == 5), ")"),
               sprintf("%.4f", theta[s480]), paste0("L", level[s480]),
               sprintf("%.2f", path_tab$gain[2]),
               sprintf("%.2f", tail(path_tab$cum_gain, 1))))

# ---- 9. 保存 ------------------------------------------------------------
saveRDS(list(dat = dat, X = X, Y = Y, clus = clus, level = level,
             level_peel = level_peel, theta = theta, omega = omega,
             progress = prog, prox = prox, path = path, path_theta = path_theta,
             u_list = u_list, v_list = v_list, prof_tab = prof_tab,
             layer_tab = layer_tab, layer_peel_tab = layer_peel_tab, tab4 = tab4,
             repro = repro, sil_tab = sil_tab, k_best = k_best,
             r_eff_score = r_eff_score, pam_agree = pam_agree, fit = fit,
             cluster_true = dat$cluster_true),
        file.path(root, "data", "chpi_results.rds"))
write.csv(dat, file.path(root, "data", "chpi_students.csv"), row.names = FALSE)
write.csv(tab4, file.path(root, "data", "layer4_cluster5.csv"), row.names = FALSE)
write.csv(path_tab, file.path(root, "data", "path_s480.csv"), row.names = FALSE)

cat("轮廓宽度:\n"); print(sil_tab)
cat("聚类一致率:", round(pam_agree, 4), "\n")
cat("各族人数:", paste(tabulate(clus), collapse = " / "), "\n")
cat("效率: 均值", round(mean(theta), 4), " 最小", round(min(theta), 4),
    " 最大", round(max(theta), 4), " r(效率,成绩) =", round(r_eff_score, 4), "\n")
cat("\n分位带分层:\n"); print(layer_tab)
cat("\n剔除式分层:\n"); print(layer_peel_tab)
cat("\nS480: 类", clus[s480], " 层", level[s480], " theta", round(theta[s480], 4), "\n")
cat("路径:", paste(path, collapse = " -> "), "\n"); print(path_tab)
cat("\n第5类第4层:", length(l4), "人\n"); print(tab4)
cat("\n锚点对照:\n"); print(repro)
