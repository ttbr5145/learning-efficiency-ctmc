###############################################################################
# 04_figures.R —— 论文插图（PNG, 300 dpi）
###############################################################################
# 仓库根目录：优先环境变量 PROJ_ROOT，其次当前工作目录（应在仓库根），最后回退本机默认路径
root <- Sys.getenv("PROJ_ROOT", unset = getwd())
if (!file.exists(file.path(root, "代码", "dea_fun.R"))) root <- "D:/AI/随机过程大作业2"
setwd(root)
suppressPackageStartupMessages({library(showtext); library(Matrix)})
showtext_auto(); showtext_opts(dpi = 300)
font_add("SimSun", "C:/Windows/Fonts/simsun.ttc")
font_add("SimHei", "C:/Windows/Fonts/simhei.ttf")

res <- readRDS(file.path(root, "数据", "chpi_results.rds"))
sp <- readRDS(file.path(root, "数据", "sp_results.rds"))
dat <- res$dat; clus <- dat$clus; theta <- dat$theta; level <- dat$level

openfig <- function(name, w, h) {
  png(file.path(root, "论文", "fig", paste0(name, ".png")), width = w, height = h,
      units = "in", res = 300, type = "cairo", family = "SimSun")
  par(family = "SimSun", las = 1)
}

# ---- 1. 雷达图 -----------------------------------------------------------
items <- c("interest", "revision", "examtimes", "onspot", "mental", "rest",
           "mocktest", "reciting", "listening", "speaking", "reading", "extracur")
lab_cn <- c("兴趣", "复习", "考试次数", "临场", "心理", "休息",
            "模拟考", "背诵", "听力", "口语", "阅读", "课外")
M <- t(sapply(1:5, function(c) colMeans(dat[clus == c, items])))
openfig("fig_radar", 6.2, 5.2)
par(mar = c(0.5, 0.5, 0.5, 0.5))
p <- ncol(M); ang <- seq(0, 2 * pi, length.out = p + 1)[-(p + 1)]
plot(0, 0, type = "n", xlim = c(-1.85, 1.35), ylim = c(-1.3, 1.3), axes = FALSE, xlab = "", ylab = "")
for (r in seq(0, 1, 0.2))
  lines(r * cos(seq(0, 2 * pi, 0.01)), r * sin(seq(0, 2 * pi, 0.01)), col = "grey85")
for (i in 1:p) lines(c(0, 1.15 * cos(ang[i])), c(0, 1.15 * sin(ang[i])), col = "grey85")
text(1.33 * cos(ang), 1.33 * sin(ang), lab_cn, cex = 0.85)
cols <- c("#4C72B0", "#DD8452", "#55A868", "#C44E52", "#8172B2")
for (c in 1:nrow(M)) polygon((M[c, ] / 100) * cos(ang), (M[c, ] / 100) * sin(ang),
                            border = cols[c], lwd = 1.6)
legend(-1.82, 1.18, paste0("第", 1:5, "类"), col = cols, lwd = 1.6, bty = "n", cex = 0.85)
dev.off()

# ---- 2. 效率—成绩散点 ---------------------------------------------------
openfig("fig_eff_score", 7.2, 4.2)
par(mar = c(4.2, 4.2, 1.2, 0.8), mgp = c(2.6, 0.8, 0))
plot(theta, dat$score, pch = 19, cex = 0.5, col = "#4C72B0AA",
     xlab = "学习效率（CCR-DEA）", ylab = "期末英语成绩")
abline(lm(score ~ theta, data = dat), col = "#C44E52", lwd = 2)
text(0.45, 95, sprintf("r = %.3f", res$r_eff_score), cex = 1, pos = 4)
dev.off()

# ---- 3. 分层散点（第 5 类与第 1 类）------------------------------------
openfig("fig_layer", 7.2, 3.8)
par(mfrow = c(1, 2), mar = c(4.0, 4.0, 1.8, 0.6), mgp = c(2.4, 0.7, 0))
cols <- c("#4C72B0", "#DD8452", "#55A868", "#C44E52", "#8172B2")
for (c in c(1, 5)) {
  sub <- which(clus == c)
  plot(theta[sub], level[sub] + runif(length(sub), -0.15, 0.15), pch = 19, cex = 0.65,
       col = cols[level[sub]], xlim = c(0.3, 1.02), ylim = c(0.6, 5.4),
       xlab = "学习效率", ylab = "效率层级", main = paste0("第", c, "类（n = ", length(sub), "）"))
  abline(h = 1:5, col = "grey88", lty = 2)
}
dev.off()

# ---- 4. CTMC 样本路径 ---------------------------------------------------
openfig("fig_traj", 7.2, 3.8)
par(mar = c(4.0, 4.0, 1.2, 0.8), mgp = c(2.4, 0.7, 0))
plot(NA, xlim = c(0, 40), ylim = c(1.1, 5.2), xlab = "时间（周）", ylab = "效率状态")
cols8 <- c("#B03A2E", "#1F618D", "#1E8449", "#D68910", "#6C3483", "#117A65")
for (b in unique(sp$traj$rep)) {
  d <- sp$traj[sp$traj$rep == b, ]
  lines(d$t, d$state, type = "s", col = cols8[b], lwd = 1.3)
}
abline(h = 1, col = "grey55", lty = 2)
text(39, 1.3, "高效态", cex = 0.85, col = "grey30", pos = 2)
dev.off()

# ---- 5. 首达时间分布 ---------------------------------------------------
openfig("fig_hitting", 7.2, 3.6)
set.seed(11); cc <- sp$c_scale; r480 <- sp$rho[480]
lam <- cc * r480 / (1 + r480); mu <- cc / (1 + r480)
hits <- numeric(3000)
for (b in 1:3000) {
  k <- 5; t <- 0
  while (k > 1) {
    up <- lam; dn <- if (k < 5) mu else 0
    t <- t + rexp(1, up + dn)
    k <- if (runif(1) < up / (up + dn)) k - 1 else k + 1
  }
  hits[b] <- t
}
par(mar = c(4.0, 4.0, 1.2, 0.8), mgp = c(2.4, 0.7, 0))
hist(hits, breaks = 40, freq = FALSE, col = "#9ec9e2", border = "white",
     xlab = "首达时间 T₁（周）", ylab = "密度", main = "")
abline(v = sp$dat$ET[480], col = "#C44E52", lwd = 2, lty = 2)
legend("topright", bty = "n", cex = 0.85,
       legend = c(sprintf("解析期望 %.1f 周", sp$dat$ET[480]),
                  sprintf("模拟均值 %.1f 周", mean(hits))))
dev.off()

# ---- 6. 平稳分布 -------------------------------------------------------
openfig("fig_stationary", 7.2, 3.6)
par(mar = c(4.0, 4.0, 1.2, 0.8), mgp = c(2.4, 0.7, 0))
barplot(rbind(sp$obs_tab / sum(sp$obs_tab), sp$exp_tab / sum(sp$exp_tab)),
        beside = TRUE, names.arg = paste0("L", 1:5), col = c("#4C72B0", "#DD8452"),
        xlab = "效率层级", ylab = "比例", border = NA, las = 1, cex.axis = 0.9)
legend("topleft", c("观测比例", "模型期望"), fill = c("#4C72B0", "#DD8452"),
       border = NA, bty = "n", cex = 0.85)
dev.off()

# ---- 7. 行为杠杆效应 ---------------------------------------------------
dlt <- sp$delta; se <- sp$se[2:6]; nm <- sp$cov_cn
openfig("fig_strategy", 7.2, 3.6)
par(mar = c(4.2, 5.2, 1.2, 0.8), mgp = c(2.5, 0.8, 0))
ord <- order(dlt)
plot(dlt[ord], 1:5, xlim = c(-0.55, 0.55), ylim = c(0.6, 5.4), yaxt = "n", pch = 19, cex = 1.2,
     col = ifelse(dlt[ord] > 0, "#C44E52", "#4C72B0"),
     xlab = expression(delta[j]), ylab = "")
axis(2, at = 1:5, labels = nm[ord], las = 1, cex.axis = 0.95, tick = FALSE)
segments(dlt[ord] - 1.96 * se[ord], 1:5, dlt[ord] + 1.96 * se[ord], 1:5, col = "grey45", lwd = 2)
abline(v = 0, lty = 2, col = "grey70")
dev.off()

# ---- 7b. 个案各动作的跃层概率增量 -------------------------------------
openfig("fig_action", 7.2, 3.4)
par(mar = c(4.2, 5.2, 1.2, 0.8), mgp = c(2.5, 0.8, 0))

# 个案 S480：单独实施每个动作的跃层概率变化
Q_of <- function(r, cc) {
  lam <- cc * r / (1 + r); mu <- cc / (1 + r)
  Q <- matrix(0, 5, 5)
  for (k in 1:5) { if (k > 1) Q[k, k - 1] <- lam; if (k < 5) Q[k, k + 1] <- mu }
  diag(Q) <- -rowSums(Q); Q
}
DZ <- matrix(unlist(sp$act$dz), nrow = 5, byrow = TRUE)
base_p <- sp$dat$p_up16[480]
dp <- sapply(1:5, function(j) {
  rnew <- sp$rho[480] * exp(sum(sp$delta * DZ[j, ]))
  P <- as.matrix(expm(Q_of(rnew, sp$c_scale) * 16))
  sum(P[5, 1:4]) - base_p
})
barplot(rev(dp), horiz = TRUE, names.arg = rev(c("复习", "模拟考", "听说读写", "作息", "兴趣拓展")),
        col = rev(ifelse(dp > 0, "#55A868", "#C44E52")), border = NA, las = 1,
        xlab = "跃层概率增量", cex.names = 0.9, cex.axis = 0.9)
abline(v = 0, col = "grey40")
dev.off()

# ---- 8. 三方案期望效率轨迹 ---------------------------------------------
openfig("fig_gain", 7.2, 4.2)
par(mar = c(4.0, 4.0, 1.2, 0.8), mgp = c(2.4, 0.7, 0))
tb <- sp$theta_bar; tt <- seq(0, 32, length.out = 65)
trajE <- function(r) sapply(tt, function(t) sum(as.matrix(expm(Q_of(r, sp$c_scale) * t))[5, ] * tb))
plot(tt, trajE(sp$rho[480]), type = "l", lwd = 2.2, col = "#4C72B0", ylim = c(0.66, 0.98),
     xlab = "时间（周）", ylab = expression(E*"["*theta[t]*"]"))
lines(tt, trajE(sp$rho_ref), lwd = 2.2, col = "#DD8452")
lines(tt, trajE(sp$dat$optRho[480]), lwd = 2.2, col = "#55A868")
abline(h = tb[1], lty = 3, col = "grey55")
legend("bottomright", c("现状", "模仿参照个体 S477", "模型最优策略"),
       col = c("#4C72B0", "#DD8452", "#55A868"), lwd = 2.2, bty = "n", cex = 0.85)
dev.off()

cat("figures done\n")

