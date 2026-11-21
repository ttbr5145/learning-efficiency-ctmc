###############################################################################
# 03_stochastic.R —— 以随机过程替换 context-dependent DEA
###############################################################################
# 仓库根目录：优先环境变量 PROJ_ROOT，其次当前工作目录（应在仓库根），最后回退本机默认路径
root <- Sys.getenv("PROJ_ROOT", unset = getwd())
if (!file.exists(file.path(root, "code", "dea_fun.R"))) root <- "D:/AI/随机过程大作业2"
setwd(root)
suppressPackageStartupMessages({library(Matrix); library(numDeriv)})
res <- readRDS(file.path(root, "data", "chpi_results.rds"))
dat <- res$dat
n <- nrow(dat)
s <- dat$level; theta <- dat$theta; clus <- dat$clus
K <- 5

# ---- A. 协变量 -----------------------------------------------------------
z_raw <- cbind(
  复习巩固 = rowMeans(dat[, c("revision", "tutoring")]),
  应试训练 = rowMeans(dat[, c("examtimes", "mocktest")]),
  输入训练 = rowMeans(dat[, c("reciting", "listening", "speaking", "reading")]),
  心理作息 = rowMeans(dat[, c("mental", "rest")]),
  兴趣拓展 = rowMeans(dat[, c("interest", "extracur", "abroad")])
)
z <- scale(z_raw); colnames(z) <- colnames(z_raw)
cov_cn <- colnames(z_raw)

# ---- B. 生灭过程的参数估计 ----------------------------------------------
nll <- function(par) {
  a <- par[1]; d <- par[2:6]
  lp <- pmin(pmax(as.numeric(a + z %*% d), -30), 30)
  lnum <- outer(lp, 0:(K - 1), function(x, k) -k * x)
  -sum(lnum[cbind(seq_len(n), s)] - log(rowSums(exp(lnum))))
}
opt <- optim(c(0, rep(0, 5)), nll, method = "BFGS", hessian = TRUE,
             control = list(maxit = 3000, reltol = 1e-13))
alpha <- opt$par[1]; delta <- opt$par[2:6]
se <- sqrt(diag(solve(opt$hessian)))
rho <- exp(pmin(pmax(as.numeric(alpha + z %*% delta), -30), 30))
lrho <- log(rho)
p_hat <- exp(outer(lrho, 0:(K - 1), function(x, k) -k * x))
p_hat <- p_hat / rowSums(p_hat)
acc <- mean(max.col(p_hat) == s)
obs_tab <- as.numeric(table(factor(s, levels = 1:K)))
exp_tab <- as.numeric(colSums(p_hat))
chi <- sum((obs_tab - exp_tab)^2 / exp_tab)
p_chi <- pchisq(chi, df = K - 2, lower.tail = FALSE)

# ---- C. 首达时间与时间标定 ----------------------------------------------
hit_time <- function(rho_i, cc) {
  lam <- cc * rho_i / (1 + rho_i); mu <- cc / (1 + rho_i)
  A <- matrix(0, 4, 4); b <- rep(-1, 4)
  for (idx in 1:4) {
    k <- idx + 1
    A[idx, idx] <- -(lam + mu)
    if (k - 1 >= 2) A[idx, idx - 1] <- lam
    if (k + 1 <= 5) A[idx, idx + 1] <- mu
  }
  A[4, 4] <- -lam
  c(0, solve(A, b))
}
Q_of <- function(rho_i, cc) {
  lam <- cc * rho_i / (1 + rho_i); mu <- cc / (1 + rho_i)
  Q <- matrix(0, K, K)
  for (k in 1:K) {
    if (k > 1) Q[k, k - 1] <- lam
    if (k < K) Q[k, k + 1] <- mu
  }
  diag(Q) <- -rowSums(Q)
  Q
}
u_unit <- t(sapply(rho, hit_time, cc = 1))
c_scale <- median(u_unit[clus == 5, 5]) / 16
E_T <- u_unit / c_scale
dat$rho <- rho; dat$ET <- E_T[cbind(1:n, s)]

P16 <- array(NA_real_, c(K, K, n))
for (i in 1:n) P16[, , i] <- as.matrix(expm(Q_of(rho[i], c_scale) * 16))
p_up16 <- sapply(1:n, function(i) if (s[i] == 1) 0 else sum(P16[s[i], 1:(s[i] - 1), i]))
dat$p_up16 <- p_up16
share_drift_down <- mean(rho < 1)

# ---- 状态—效率映射（把转移概率折算成效率）-------------------------------
theta_bar <- sapply(1:K, function(k) mean(theta[clus == 5 & s == k]))
Etheta16 <- sapply(1:n, function(i) sum(P16[s[i], , i] * theta_bar))
dat$Etheta16 <- Etheta16
dat$gain16 <- (Etheta16 - theta) / theta

# ---- D. 蒙特卡洛校验 ----------------------------------------------------
mc_ctmc <- function(rho_i, cc, k0 = 5, tmax = 400, npath = 1000) {
  hits <- numeric(npath)
  lam <- cc * rho_i / (1 + rho_i); mu <- cc / (1 + rho_i)
  for (b in 1:npath) {
    k <- k0; t <- 0
    if (k == 1) { hits[b] <- 0; next }
    while (t < tmax) {
      up <- if (k > 1) lam else 0
      dn <- if (k < K) mu else 0
      tot <- up + dn
      if (tot <= 0) break
      t <- t + rexp(1, tot)
      k <- if (runif(1) < up / tot) k - 1 else k + 1
      if (k == 1) { hits[b] <- t; break }
    }
  }
  hits[hits == 0] <- NA
  mean(hits, na.rm = TRUE)
}
set.seed(20261121)
mc_ids <- c(480, 477, 456, 494, 293)
mc_rows <- data.frame(
  id = mc_ids, state = s[mc_ids], rho = round(rho[mc_ids], 3),
  解析解 = round(dat$ET[mc_ids], 2),
  蒙特卡洛 = round(mapply(mc_ctmc, rho[mc_ids], MoreArgs = list(cc = c_scale),
                        k0 = s[mc_ids]), 2)
)
mc_rows$相对误差 <- round(abs(mc_rows$蒙特卡洛 - mc_rows$解析解) / mc_rows$解析解 * 100, 1)

sim_path <- function(rho_i, cc, tmax = 60) {
  k <- 5; t <- 0; tt <- 0; kk <- 5
  lam <- cc * rho_i / (1 + rho_i); mu <- cc / (1 + rho_i)
  while (t < tmax) {
    up <- if (k > 1) lam else 0
    dn <- if (k < K) mu else 0
    tot <- up + dn
    if (tot <= 0) break
    t <- t + rexp(1, tot)
    tt <- c(tt, t); kk <- c(kk, k)
    k <- if (runif(1) < up / tot) k - 1 else k + 1
    if (k == 1) { tt <- c(tt, t); kk <- c(kk, 1); break }
  }
  data.frame(t = tt, state = kk)
}
set.seed(7)
traj <- do.call(rbind, lapply(1:6, function(b) { p <- sim_path(rho[480], c_scale); p$rep <- b; p }))

# ---- E. OU 扩散过程 -----------------------------------------------------
theta480 <- theta[480]
kappa <- -log(1 - 0.2538 * theta480 / (1 - theta480)) / 16
ou_ET <- function(theta0, theta_star, sig, kap, xbar = 0.95, m = 3000) {
  x <- seq(0, xbar, length.out = m); dx <- x[2] - x[1]
  ls <- -(kap / sig^2) * (x - theta_star)^2; ls <- ls - max(ls)
  sf <- exp(ls)
  cum <- c(0, cumsum((sf[-m] + sf[-1]) / 2 * dx))
  du <- -2 * cum / (sig^2 * sf)
  u <- -rev(cumsum(rev((du[1:(m-1)] + du[2:m]) / 2 * dx))); u <- c(u, 0)
  approx(x, u, xout = min(theta0, xbar - 1e-6))$y
}
sig_c <- sapply(1:K, function(c) sqrt(2 * kappa) * sd(theta[clus == c]))
sig_i <- sig_c[clus]
mu_c <- sapply(1:K, function(c) mean(theta[clus == c]))
ET_ou <- mapply(function(th, sg) ou_ET(th, 1, sg, kappa), theta, sig_i)
ET_ou_no <- mapply(function(th, sg, mu) ou_ET(th, mu, sg, kappa), theta, sig_i, mu_c[clus])
cat("OU: sig_c =", paste(round(sig_c, 4), collapse = ", "), "\n")
E_theta16_ou <- 1 - (1 - theta) * exp(-kappa * 16)
gain16_ou <- (E_theta16_ou - theta) / theta
T43 <- -log(1 - 0.4397 * theta480 / (1 - theta480)) / kappa

# ---- F. 预算约束下的策略优化 -------------------------------------------
act <- data.frame(
  name = c("提高复习频次", "增加模拟考与真题", "增加听说读写训练",
           "改善作息与考前调适", "提升兴趣与课外拓展"),
  hours = c(2.0, 2.0, 3.0, 1.0, 1.5),
  dz = I(list(c(1,0,0,0,0), c(0,1,0,0,0), c(0,0,1,0,0), c(0,0,0,1,0), c(0,0,0,0,1)))
)
DZ <- matrix(unlist(act$dz), nrow = 5, byrow = TRUE)
budget <- 3
combos <- as.matrix(expand.grid(rep(list(0:1), 5)))
combos <- combos[as.numeric(combos %*% act$hours) <= budget, , drop = FALSE]
dlogrho_all <- combos %*% DZ %*% delta
best_for <- function(i) {
  newrho <- as.numeric(rho[i] * exp(dlogrho_all))
  newP <- sapply(newrho, function(r) if (s[i] == 1) 0 else {
    P <- as.matrix(expm(Q_of(r, c_scale) * 16)); sum(P[s[i], 1:(s[i] - 1)]) })
  newET <- sapply(newrho, function(r) hit_time(r, c_scale)[s[i]])
  o <- order(-newP, newET)[1]
  list(k = combos[o, ], p = newP[o], ET = newET[o], rho = newrho[o])
}
best_all <- lapply(1:n, best_for)
dat$optP <- sapply(best_all, `[[`, "p")
dat$optET <- sapply(best_all, `[[`, "ET")
dat$optRho <- sapply(best_all, `[[`, "rho")
opt_actions <- t(sapply(best_all, `[[`, "k"))
colnames(opt_actions) <- act$name
act_share <- colMeans(opt_actions)

ref <- res$path[2]
dz_ref <- as.numeric(z[ref, ] - z[480, ])
rho_ref <- rho[480] * exp(sum(delta * dz_ref))
P_ref <- {P <- as.matrix(expm(Q_of(rho_ref, c_scale) * 16)); sum(P[5, 1:4])}
ET_ref <- hit_time(rho_ref, c_scale)[5]
hours_ref <- sum(abs(dz_ref)) * 4
Etheta16_ref <- sum({P <- as.matrix(expm(Q_of(rho_ref, c_scale) * 16)); P[5, ]} * theta_bar)
Etheta16_opt <- sum({P <- as.matrix(expm(Q_of(dat$optRho[480], c_scale) * 16)); P[5, ]} * theta_bar)

# ---- G. 导出表格 -------------------------------------------------------
write_tex <- function(df, file, caption, label, align = NULL, digits = 3, size = "\\small") {
  if (is.null(align)) align <- paste0("l", paste(rep("c", ncol(df) - 1), collapse = ""))
  con <- file(file, "w", encoding = "UTF-8")
  writeLines(c("\\begin{table}[htbp]", "\\centering", size,
               paste0("\\caption{", caption, "}"), paste0("\\label{", label, "}"),
               paste0("\\begin{tabular}{", align, "}"), "\\toprule",
               paste(paste(colnames(df), collapse = " & "), "\\\\"), "\\midrule"), con)
  for (i in seq_len(nrow(df))) {
    writeLines(paste(paste(sapply(df[i, ], function(v)
      if (is.numeric(v)) formatC(v, format = "f", digits = digits) else as.character(v)), collapse = " & "), "\\\\"), con)
  }
  writeLines(c("\\bottomrule", "\\end{tabular}", "\\end{table}"), con)
  close(con)
}

par_tab <- data.frame(变量 = c("常数项 $\\alpha$", cov_cn),
                      估计 = c(alpha, delta), 标准误 = c(se[1], se[2:6]),
                      "$e^{\\hat{\\delta}}$" = c(exp(alpha), exp(delta)))
write_tex(par_tab, file.path(root, "tex", "tab_ctmc_par.tex"),
          "生灭过程参数估计（$\\rho_i=\\exp(\\alpha+\\delta^{\\top}z_i)$）", "tab:ctmcpar",
          align = "lccc", digits = 3)

f5 <- which(clus == 5)
l4 <- f5[s[f5] == 4]
stu_tab <- data.frame(个体 = paste0("S", l4), 层级 = s[l4], 效率 = theta[l4],
                      当前 = rho[l4], 首达 = dat$ET[l4], 优化 = dat$optRho[l4],
                      优化首达 = dat$optET[l4], 跃层概率 = dat$p_up16[l4],
                      优化概率 = dat$optP[l4])
colnames(stu_tab) <- c("个体", "层级", "效率", "$\\rho$", "$E[T_1]$", "$\\rho'$",
                       "$E[T_1]'$", "$p$", "$p'$")
stu_tab <- stu_tab[order(stu_tab[, 5]), ]
write_tex(stu_tab, file.path(root, "tex", "tab_time.tex"),
          "第 5 类第 4 层个体的期望首达时间与策略优化效果", "tab:time",
          align = "lcccccccc", digits = 2)

s480_tab <- data.frame(
  方案 = c("现状（不干预）", "原文式：模仿参照个体 S477", "模型最优（预算 3 学时/周）"),
  学时 = c(0, round(hours_ref, 1), 3),
  rho = c(rho[480], rho_ref, dat$optRho[480]),
  ET = c(dat$ET[480], ET_ref, dat$optET[480]),
  p = c(dat$p_up16[480], P_ref, dat$optP[480]),
  Eth = c(dat$Etheta16[480], Etheta16_ref, Etheta16_opt)
)
s480_tab$增益 <- (s480_tab$Eth - theta480) / theta480 * 100
colnames(s480_tab) <- c("方案", "额外学时/周", "$\\rho$", "$E[T_1]$（周）",
                        "16 周跃层概率", "16 周期望效率", "相对提升(\\%)")
write_tex(s480_tab, file.path(root, "tex", "tab_s480.tex"),
          "S480 个案：原文参考路径与随机过程策略的对照", "tab:s480",
          align = "lcccccc", digits = 3)

gain_tab <- data.frame(
  时间尺度 = c("一学期（16 周）", "两学期（32 周）", "反解原文累计增益 43.97\\%"),
  期望效率 = c(1 - (1 - theta480) * exp(-kappa * 16),
              1 - (1 - theta480) * exp(-kappa * 32),
              1 - (1 - theta480) * exp(-kappa * T43)),
  相对提升 = c(gain16_ou[480] * 100,
              (1 - theta480) * (1 - exp(-kappa * 32)) / theta480 * 100, 43.97),
  所需周数 = c(16, 32, T43))
colnames(gain_tab) <- c("时间尺度", "期望效率", "相对提升(\\%)", "所需周数")
write_tex(gain_tab, file.path(root, "tex", "tab_gain.tex"),
          "OU 扩散模型给出的效率提升时间尺度", "tab:gain", align = "lccc", digits = 2)

# ---- 保存 -------------------------------------------------------------
saveRDS(list(dat = dat, z = z, z_raw = z_raw, alpha = alpha, delta = delta, se = se,
             rho = rho, c_scale = c_scale, E_T = E_T, P16 = P16, p_up16 = p_up16,
             acc = acc, obs_tab = obs_tab, exp_tab = exp_tab, chi = chi, p_chi = p_chi,
             mc_rows = mc_rows, traj = traj, kappa = kappa, sig_c = sig_c, sig_i = sig_i,
             ET_ou = ET_ou, ET_ou_no = ET_ou_no, gain16_ou = gain16_ou, T43 = T43,
             E_theta16_ou = E_theta16_ou, theta_bar = theta_bar,
             opt_actions = opt_actions, act_share = act_share, act = act,
             stu_tab = stu_tab, s480_tab = s480_tab, gain_tab = gain_tab,
             par_tab = par_tab, share_drift_down = share_drift_down,
             ref = ref, rho_ref = rho_ref, P_ref = P_ref, ET_ref = ET_ref,
             hours_ref = hours_ref, cov_cn = cov_cn, mu_c = mu_c),
        file.path(root, "data", "sp_results.rds"))
write.csv(dat, file.path(root, "data", "sp_students.csv"), row.names = FALSE)
write.csv(mc_rows, file.path(root, "data", "mc_check.csv"), row.names = FALSE)

cat("alpha =", round(alpha, 3), "  delta =", paste(round(delta, 3), collapse = ", "), "\n")
cat("se =", paste(round(se, 3), collapse = ", "), "\n")
cat("rho 范围:", round(range(rho), 3), " rho<1 占比:", round(share_drift_down, 3), "\n")
cat("分类准确率:", round(acc, 3), " 卡方:", round(chi, 2), " p =", round(p_chi, 4), "\n")
cat("观测:", paste(obs_tab, collapse = "/"), " 期望:", paste(round(exp_tab, 1), collapse = "/"), "\n")
cat("c =", round(c_scale, 3), " E[T] 中位(第5类):", round(median(dat$ET[f5]), 2), "\n")
print(mc_rows)
cat("kappa =", round(kappa, 4), " T43 =", round(T43, 1), "\n")
cat("ET_ou 范围:", paste(round(range(ET_ou, na.rm = TRUE), 1), collapse = " ~ "),
    " 无干预:", paste(round(range(ET_ou_no, na.rm = TRUE), 1), collapse = " ~ "), "\n")
cat("动作比例:", paste(round(act_share, 3), collapse = ", "), "\n")
print(s480_tab)
cat("theta_bar:", paste(round(theta_bar, 3), collapse = ", "), "\n")





