###############################################################################
# 07_export_web.R —— 把论文与随机过程模型的关键数值导出为网页用 JSON
#   输出：数据/web_data.json
###############################################################################
root <- Sys.getenv("PROJ_ROOT", unset = getwd())
if (!file.exists(file.path(root, "\u4ee3\u7801", "dea_fun.R"))) root <- "D:/AI/\u968f\u673a\u8fc7\u7a0b\u5927\u4f5c\u4e1a2"
setwd(root)

suppressPackageStartupMessages({ library(jsonlite); library(Matrix) })

sp  <- readRDS(file.path(root, "\u6570\u636e", "sp_results.rds"))
ch  <- read.csv(file.path(root, "\u6570\u636e", "chpi_students.csv"))
tb  <- read.csv(file.path(root, "\u6570\u636e", "textbook_table.csv"))

alpha <- sp$alpha; delta <- sp$delta; c_scale <- sp$c_scale
Z <- sp$z                      # 519 x 5 标准化协变量
rho <- sp$rho                  # 个体强度比
lev <- sp$dat$level            # 层级 1..5（1 = 最高效）
clus <- sp$dat$clus
stopifnot(nrow(Z) == length(rho))

# ---- 情景：兴趣拓展每小时对应的标准化协变量增量（取自 03 的动作网格）------
act <- sp$act
dz_per_hour <- as.numeric(act$dz[[5]]) / act$hours[5]

rho_opt <- exp(alpha + as.numeric(Z %*% delta) + as.numeric(Z %*% dz_per_hour) * 0 +
               as.numeric(act$dz[[5]] %*% delta))

# ---- 分层平均强度比与率 -----------------------------------------------------
k_by_lev <- lapply(1:5, function(k) which(lev == k))
rho_k     <- sapply(k_by_lev, function(i) mean(rho[i]))
rho_k_opt <- sapply(k_by_lev, function(i) mean(rho_opt[i]))
n_k       <- sapply(k_by_lev, length)
n_clus    <- as.numeric(table(clus))

# ---- 首达时间（对起始层级取均值）-------------------------------------------
ET_by_start <- colMeans(sp$E_T[, 1:5])

# ---- 16 周转移矩阵（按个体平均）---------------------------------------------
P16_mean <- apply(sp$P16, c(1, 2), mean)

# ---- 复现拟合：观测 vs 期望转移计数 ----------------------------------------
fit <- data.frame(obs = as.numeric(sp$obs_tab), exp = as.numeric(sp$exp_tab))

# ---- rho 直方图 -------------------------------------------------------------
brks <- seq(min(rho), max(rho), length.out = 25)
cnt  <- hist(rho, breaks = brks, plot = FALSE)$counts

# ---- 组装 JSON --------------------------------------------------------------
# ---- 表格：改成 ASCII 键名，便于网页端直接读取 -----------------------------
s480 <- sp$s480_tab; names(s480) <- c("plan", "hours", "rho", "ET", "p16", "theta16", "gain")
gain <- sp$gain_tab; names(gain)  <- c("scale", "theta", "gain", "weeks")
par  <- sp$par_tab;  names(par)   <- c("term", "est", "se", "exp_delta")
tb2  <- tb;          names(tb2)   <- c("cls", "rho", "lambda2", "halflife", "mixing",
                                       "pi1_jump", "pi1_cont", "tv16", "maxdiff")

# ---- 个体级数据：网页端在浏览器里直接重算 rho 与转移矩阵 -------------------
students <- data.frame(
  id    = sp$dat$id,
  layer = as.integer(lev),
  clus  = as.integer(clus),
  rho   = round(rho, 5),
  theta = round(sp$dat$theta, 4),
  score = round(sp$dat$score, 1),
  z1 = round(Z[, 1], 4), z2 = round(Z[, 2], 4), z3 = round(Z[, 3], 4),
  z4 = round(Z[, 4], 4), z5 = round(Z[, 5], 4)
)
actions <- data.frame(name = act$name, hours = act$hours)
actions$dz <- I(lapply(act$dz, function(v) round(as.numeric(v), 4)))

out <- list(
  meta = list(
    n = length(rho), layers = n_k, classes = n_clus,
    theta_mean = mean(ch$theta), eff_score_r = cor(sp$dat$theta, sp$dat$score),
    alpha = alpha, c = c_scale,
    rho_median = median(rho), share_lt1 = mean(rho < 1),
    rho_opt_mean = mean(rho_opt),
    ET_first = ET_by_start[1], T43 = sp$T43,
    kappa = sp$kappa, chi = sp$chi, p_chi = sp$p_chi,
    p16_layer5 = mean(sp$p_up16[lev == 5]),
    hours_opt = act$hours[5], act_name = act$name[5]
  ),
  covary = list(
    names = sp$cov_cn,
    delta = as.numeric(delta),
    zbar = t(sapply(k_by_lev, function(i) colMeans(Z[i, , drop = FALSE]))),
    dz_per_hour = dz_per_hour
  ),
  actions = actions,
  students = students,
  rho = list(values = rho, hist_breaks = brks, hist_counts = cnt,
             layer_mean = rho_k, layer_mean_opt = rho_k_opt),
  chain = list(P16 = P16_mean, ET_by_start = as.numeric(ET_by_start),
               theta_bar = as.numeric(sp$theta_bar)),
  fit  = fit,
  strategy = list(s480 = s480, gain = gain, par = par),
  textbook = tb2,
  layer_sizes = n_k
)

json <- toJSON(out, auto_unbox = TRUE, digits = 6, pretty = TRUE)
writeLines(json, file.path(root, "\u6570\u636e", "web_data.json"), useBytes = TRUE)

cat("=== \u5bfc\u51fa\u5b8c\u6210 ===\n")
cat("n =", length(rho), " | \u5c42\u7ea7\u4eba\u6570:", n_k, "\n")
cat("rho \u4e2d\u4f4d\u6570 =", round(median(rho), 4), " | rho<1 \u5360\u6bd4 =", round(mean(rho < 1) * 100, 1), "%\n")
cat("\u5c42\u5747 rho:", round(rho_k, 4), " | \u6700\u4f18:", round(rho_k_opt, 4), "\n")
cat("\u6700\u4f18\u60c5\u666f\u5e73\u5747 rho =", round(mean(rho_opt), 4), " (\u8bba\u6587 0.9541)\n")
cat("\u9996\u8fbe\u65f6\u95f4(\u6309\u8d77\u59cb\u5c42):", round(ET_by_start, 3), "\n")
cat("\u7b2c5\u5c42 16 \u5468\u8dc3\u5c42\u6982\u7387\u5747\u503c =", round(mean(sp$p_up16[lev == 5]), 4), "\n")
cat("\u5199\u51fa:", file.path(root, "\u6570\u636e", "web_data.json"), "\n")
