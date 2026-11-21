# =============================================================================
#  06_textbook.R
#  与教材（《随机过程》，Ch3 Poisson / Ch5 Markov 链 / Ch6 鞅 / Ch7 Brown
#  运动 / Ch8 随机积分）逐章对接的补充计算：
#    (1) 嵌入链与平稳分布、谱隙与混合时间（Ch5.3 极限定理）
#    (2) Poisson 一致化：跳事件计数、上行事件计数（Ch3.1–3.2）
#    (3) 调和函数与鞅停时定理互验首达期望（Ch6.2）
#    (4) Brown 运动反射原理、最大值与反正弦律（Ch7.5）
#    (5) Itô 公式给出的 OU 闭式解、方差与平稳分布（Ch8.3–8.4）
#  输出：数据/textbook_checks.txt、数据/textbook_table.csv
# =============================================================================

suppressPackageStartupMessages({
  library(Matrix)
})

# 仓库根目录：优先环境变量 PROJ_ROOT，其次当前工作目录（应在仓库根），最后回退本机默认路径
root <- Sys.getenv("PROJ_ROOT", unset = getwd())
if (!file.exists(file.path(root, "代码", "dea_fun.R"))) root <- "D:/AI/随机过程大作业2"
setwd(root)

sp <- readRDS(file.path(root, "数据", "sp_results.rds"))
d  <- sp$dat
cc <- sp$c_scale
K  <- 5

rho_c <- tapply(d$rho, d$cluster_true, median)   # 类内代表性行强度比（中位数）

# ---- 生成矩阵 ---------------------------------------------------------------
make_Q <- function(rho, c = cc, K = 5) {
  lam <- c * rho / (1 + rho)   # 上行强度 lambda_k
  mu  <- c / (1 + rho)        # 下行强度 mu_k
  Q <- matrix(0, K, K)
  for (k in 1:K) {
    if (k > 1) Q[k, k - 1] <- lam
    if (k < K) Q[k, k + 1] <- mu
    Q[k, k] <- -sum(Q[k, ])
  }
  Q
}

# ---- (1) 谱隙与混合时间 ------------------------------------------------------
spec_tab <- data.frame()
for (cl in 1:5) {
  Q  <- make_Q(rho_c[cl])
  ev <- eigen(Q)$values
  ev <- sort(Re(ev), decreasing = TRUE)          # 0 = lambda_1 > lambda_2 > ...
  gap <- -ev[2]
  t_mix <- log(1 / 0.01) / gap                  # 到 1% 总变差距离
  half  <- log(2) / gap                         # 半衰期
  # 嵌入链（跳链）：P^J = I + Q/c，因 lambda+mu=c 行和为 0 的对角恰为 0
  PJ <- diag(K) + Q / cc
  # 跳链平稳分布：解 pi^J PJ = pi^J
  A <- t(PJ) - diag(K); A[K, ] <- 1
  b <- c(rep(0, K - 1), 1)
  piJ <- as.numeric(solve(A, b))
  # 连续时间链平稳分布 pi_k ∝ rho^{-(k-1)}
  r <- 1 / rho_c[cl]
  pik <- r^(0:(K - 1)); pik <- pik / sum(pik)
  # 16 周转移矩阵与平稳分布的总变差距离（从最不利的起点出发）
  P16 <- as.matrix(expm(Q * 16))
  tv  <- max(apply(P16, 1, function(row) 0.5 * sum(abs(row - pik))))
  spec_tab <- rbind(spec_tab, data.frame(
    类 = cl, rho = rho_c[cl], lambda2 = -gap, 半衰期周 = half,
    混合周 = t_mix, 跳链_pi1 = piJ[1], 连续_pi1 = pik[1],
    十六周总变差 = tv,
    跳链与连续最大差 = max(abs(piJ - pik))
  ))
}

# ---- (2) Poisson 一致化 -----------------------------------------------------
Q480 <- make_Q(d$rho[d$id == 480])
T16  <- 16
# 用 Poisson 加权和重建 P(16)：sum_n e^{-cT}(cT)^n/n! (P^J)^n
PJ480 <- diag(K) + Q480 / cc
P_uni <- matrix(0, K, K)
term  <- diag(K)
w     <- exp(-cc * T16)
P_uni <- w * term
for (n in 1:400) {
  term <- term %*% PJ480
  w <- w * (cc * T16) / n
  P_uni <- P_uni + w * term
}
P_exp <- as.matrix(expm(Q480 * T16))
uni_err <- max(abs(P_uni - P_exp))
lam480 <- cc * d$rho[d$id == 480] / (1 + d$rho[d$id == 480])
mu480  <- cc / (1 + d$rho[d$id == 480])
ev480  <- sort(Re(eigen(Q480)$values), decreasing = TRUE)

# ---- (3) 调和函数与鞅停时定理 -----------------------------------------------------
# h 满足 Qh = 0，边界 h(1)=0, h(5)=1：h(k) = P_k(先到 1 层还是先到 5 层) 中的
# 吸收概率；而 m(S_t) + t 是鞅，停时定理给出 E[T_1] = m(S_0)，可与式(28) 互验
harm <- function(Q) {
  A <- Q[2:4, 2:4]
  rhs <- -Q[2:4, c(1, 5)] %*% c(0, 1)
  h <- c(0, as.numeric(solve(A, rhs)), 1)
  h
}
h480 <- harm(Q480)
m480 <- sp$E_T[480, ]
cat("S480 调和函数 h(k):", round(h480, 4), "\n")
cat("S480 首达期望 m_k:", round(m480, 4), "\n")
et_mc <- sp$mc_rows[["蒙特卡洛"]][sp$mc_rows$state == 5]   # Gillespie 均值
cat("停时定理校验 |m(5)-Gillespie 均值| =", abs(m480[5] - et_mc), "\n")
cat("先到 1 层的概率 g = 1-h（k=2,3,4）:", round(1 - h480[2:4], 4), "\n")

# ---- (4) Brown 运动的反射原理 -------------------------------------------------
sig480 <- sp$sig_c[d$cluster_true[d$id == 480]]
kap    <- sp$kappa
x0     <- d$theta[d$id == 480]
b      <- 0.95
dd     <- b - x0
T      <- 16
mu0    <- kap * (b - x0)                        # 取局部漂移 κ(θ*-θ0)，θ*=b
P_drift <- pnorm((mu0 * T - dd) / (sig480 * sqrt(T))) +
           exp(2 * mu0 * dd / sig480^2) *
           pnorm((-mu0 * T - dd) / (sig480 * sqrt(T)))
P_pure  <- 2 * (1 - pnorm(dd / (sig480 * sqrt(T))))   # κ=0，纯波动
# 漂移到 1（θ*=1）时的同一概率
mu1 <- kap * (1 - x0)
P_to1 <- pnorm((mu1 * T - dd) / (sig480 * sqrt(T))) +
         exp(2 * mu1 * dd / sig480^2) *
         pnorm((-mu1 * T - dd) / (sig480 * sqrt(T)))
# 反正弦律：漂移为 0 时，最大值时刻在 [0,T] 上的位置 ~ Beta(1/2,1/2)
arc_q <- qbeta(c(0.25, 0.5, 0.75), 0.5, 0.5)
arc_25 <- (2 / pi) * asin(sqrt(0.25))

# ---- (5) Itô 闭式解：方差与平稳分布 ------------------------------------------
var_tab <- data.frame()
for (cl in 1:5) {
  s <- sp$sig_c[cl]
  sd16  <- s * sqrt((1 - exp(-2 * kap * 16)) / (2 * kap))
  sdinf <- s / sqrt(2 * kap)
  sdobs <- sd(d$theta[d$cluster_true == cl])
  var_tab <- rbind(var_tab, data.frame(
    类 = cl, sigma = s, sd_16周 = sd16, 平稳标准差 = sdinf, 类内实测标准差 = sdobs
  ))
}

# ---- (4b) OU 的蒙特卡洛校验（反射原理不能跨模型照搬） ----------------------
ou_mc <- function(theta_star, M = 2000, dt = 0.005, T = 16, seed = 20261122) {
  set.seed(seed)
  kap <- 0.0475; sig <- 0.02793; x0 <- 0.6771; b <- 0.95
  X <- matrix(x0, M, 1); hit <- rep(FALSE, M)
  for (i in 1:(T / dt)) {
    X <- X + kap * (theta_star - X) * dt + sig * rnorm(M, 0, sqrt(dt))
    hit <- hit | (X >= b)
  }
  mean(hit)
}
P_mc_095 <- ou_mc(0.95)
P_mc_1   <- ou_mc(1)
sd_eff   <- sig480 * sqrt((1 - exp(-2 * kap * 16)) / (2 * kap))

# ---- 输出 -------------------------------------------------------------------
sink(file.path(root, "数据", "textbook_checks.txt"))
cat("=== 1. 谱隙与混合时间（Ch5.3） ===\n")
print(spec_tab, digits = 5)
cat("\n=== 2. Poisson 一致化（Ch3） ===\n")
cat("c =", cc, "周^-1；16 周内跳事件期望次数 =", round(cc * T16, 2), "\n")
cat("P(N=0) =", exp(-cc * T16), "\n")
cat("一致化重建 P(16) 与 expm(Q*16) 最大绝对误差 =", uni_err, "\n")
cat("S480: lambda =", round(lam480, 4), " mu =", round(mu480, 4),
    " 上行事件 16 周期望 =", round(lam480 * T16, 2), "\n")
cat("S480 上行等待均值(周) =", round(1 / lam480, 3),
    " 中位 =", round(log(2) / lam480, 3), "\n")
cat("S480 Q 的特征值:", round(ev480, 4), "\n")
cat("\n=== 3. 鞅与停时定理（Ch6.2） ===\n")
cat("S480 调和函数 h(k):", round(h480, 4), "\n")
cat("S480 m_k:", round(m480, 4), "\n")
cat("E[T_1|5] (来自 m):", round(m480[5], 4), " (来自 Gillespie 均值):", et_mc, "\n")
cat("\n=== 4. Brown 运动反射原理与反正弦律（Ch7.5） ===\n")
cat("S480: x0 =", round(x0, 4), " sigma =", round(sig480, 4), " kappa =", kap, "\n")
cat("16 周内触达 0.95 的概率（含漂移 μ=κ(0.95-x0)）:", round(P_drift, 4), "\n")
cat("16 周内触达 0.95 的概率（纯波动，κ=0）:", round(P_pure, 4), "\n")
cat("16 周内触达 0.95 的概率（漂移指向 1）:", round(P_to1, 4), "\n")
cat("反正弦律分位数（0.25/0.5/0.75）:", round(arc_q, 4), "；P(argmax≤T/4) =", round(arc_25, 4), "\n")
cat("OU 蒙特卡洛：θ*=0.95 时 P(T≤16) =", round(P_mc_095, 4),
    "；θ*=1 时 =", round(P_mc_1, 4), "\n")
cat("16 周有效标准差 =", round(sd_eff, 4), "；BM 的 σ√T =", round(sig480 * 4, 4), "\n")
cat("\n=== 5. Itô 闭式解与方差（Ch8） ===\n")
print(var_tab, digits = 4)
sink()

write.csv(spec_tab, file.path(root, "数据", "textbook_table.csv"), row.names = FALSE, fileEncoding = "UTF-8")
cat("textbook checks done\n")

