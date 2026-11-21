###############################################################################
# dea_fun.R —— 数据包络分析（DEA）工具函数（基于 lpSolve 线性规划）
###############################################################################
suppressPackageStartupMessages(library(lpSolve))

# CCR 投入导向包络模型：min theta
#   s.t. sum_j lambda_j x_ij - theta x_i0 <= 0
#        -sum_j lambda_j y_rj <= -y_r0
# 返回 theta、lambda、对偶乘子 (u,v)
dea_ccr <- function(X, Y, idx, S = NULL, tol = 1e-9) {
  X <- as.matrix(X); Y <- as.matrix(Y)
  if (is.null(S)) S <- seq_len(nrow(X))
  Xs <- X[S, , drop = FALSE]; Ys <- Y[S, , drop = FALSE]
  k <- length(S)
  x0 <- X[idx, ]; y0 <- Y[idx, ]
  m <- ncol(X); s <- ncol(Y)
  # 变量：lambda_1..lambda_k, theta
  obj <- c(rep(0, k), 1)
  con <- rbind(
    cbind(t(Xs), -x0),          # 投入约束 <= 0
    cbind(-t(Ys), rep(0, s))    # 产出约束 <= -y0
  )
  dir <- c(rep("<=", m), rep("<=", s))
  rhs <- c(rep(0, m), -y0)
  res <- lp("min", obj, con, dir, rhs, all.int = FALSE)
  if (res$status != 0) return(list(theta = NA, lambda = NULL, u = NULL, v = NULL, refs = NULL))
  lam <- res$solution[1:k]
  theta <- res$solution[k + 1]
  # 对偶乘子：解乘子型 LP（规模更小，变量 m+s 个，约束 k+1 个）
  obj2 <- c(y0, rep(0, m))                       # max u'y0
  con2 <- rbind(
    c(rep(0, s), x0),                            # v'x0 = 1
    cbind(Ys, -Xs)                               # u'y_j - v'x_j <= 0
  )
  dir2 <- c("=", rep("<=", k))
  rhs2 <- c(1, rep(0, k))
  res2 <- lp("max", obj2, con2, dir2, rhs2, all.int = FALSE)
  u <- if (res2$status == 0) res2$solution[1:s] else rep(NA, s)
  v <- if (res2$status == 0) res2$solution[(s + 1):(s + m)] else rep(NA, m)
  refs <- S[which(lam > tol)]
  list(theta = theta, lambda = lam, u = u, v = v, refs = refs, S = S)
}

# 产出导向模型（用于 context-dependent DEA 的吸引力/进步度）
#   max phi  s.t. sum_j lambda_j x_ij <= x_i0 ; sum_j lambda_j y_rj >= phi y_r0
dea_output <- function(X, Y, x0, y0, S) {
  Xs <- X[S, , drop = FALSE]; Ys <- Y[S, , drop = FALSE]
  k <- length(S); m <- ncol(X); s <- ncol(Y)
  obj <- c(rep(0, k), 1)
  con <- rbind(
    cbind(t(Xs), rep(0, m)),
    cbind(-t(Ys), y0)
  )
  dir <- c(rep("<=", m), rep("<=", s))
  rhs <- c(x0, rep(0, s))
  res <- lp("max", obj, con, dir, rhs, all.int = FALSE)
  if (res$status != 0) return(NA_real_)
  res$solution[k + 1]
}

# 使用某组最优乘子 (u,v) 与“参照个体的投入”计算反事实效率
#   theta_new = u'y0 / v'x_ref
counterfactual_eff <- function(u, v, y0, x_ref) {
  as.numeric(sum(u * y0) / sum(v * x_ref))
}
