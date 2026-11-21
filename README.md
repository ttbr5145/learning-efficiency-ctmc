# 学习效率层级演化的随机过程建模与策略推荐

以 519 名学生的英语学习问卷为例：把学习效率聚成五级层级，先复现原文的 context-dependent DEA「分层—吸引力—进步度—策略路径」证据链，再把层级当作状态空间，改用**带协变量的生灭过程 + OU 扩散**，用同一个对象回答"层级怎么演化、多久能升一级、每周多出的学时该投到哪里"。

东北林业大学《随机过程》课程大作业。作者：吴桐、赵睿、张文琦、吴梦豪（理学院 数学与应用数学 2024 级）；授课教师：艾晓辉。论文分**个人版**（`论文/paper_solo.pdf`）与**团队版**（`论文/paper_team.pdf`）两份，正文、图表、公式完全一致，差别只在封面署名与页眉署名。篇幅 49 页、9 幅插图、21 张表、57 个编号公式。

除了 PDF，还有一份**可交互的网页沙盘**：<https://xunsuiwuchu.site/learning-efficiency/>，层级分布、动作剂量、个体跃层、转移矩阵都由 519 条个体数据在浏览器里现场重算。

原文献：Zhao L T, et al. A recommendation system for effective learning strategies: An integrated approach using context-dependent DEA. *Expert Systems with Applications*, 2023, 211: 118535.

---

## 一句话结论

生灭过程的强度比 $\rho_i=\exp(\alpha+\bm\delta^\top\bm z_i)$ 把"谁在往上走、多快走上去"压成 6 个可估参数：$\rho$ 中位数 0.682、98.1% 的个体 $\rho<1$，说明层级整体向下漂移；但只要有 $+1.5$ h/周 投到「兴趣与课外拓展」，$\rho$ 就升到 0.954、期望效率提升 30.39%，而原文推荐的"模仿 S477、每周多投 4.1 h"反而把 $\rho$ 压到 0.669。

## 方法要点

| 步骤 | 做法 | 章节 |
| --- | --- | --- |
| 破解 | 把原文五环证据链拆成可检验的数学表述：复合量化 → Gower+PAM 聚类 → 类内 CCR-DEA → context-dependent 分层 → 三准则与参考路径；逐环标出原文未披露的自由度与内部不一致（含 $H=1/\Omega+P$ 与表中数值不自洽一处） | 三（一） |
| 复现 | 原始问卷不可得，用四条锚点（$n=519$、5 类、S480 效率 0.6771、单步 25.38%／累计 43.97%）逆向重建 519 份样本，用 R 复现 CHPI：PAM-Gower 聚类（$k=3,\dots,8$ 轮廓系数全报、与设计标签一致率 93.6%）、CCR-DEA（lpSolve 求解，包络式与乘子式互验）、分位带与剔除式两种分层口径、三准则与 S480 四步路径 | 三（二） |
| 换方法 | 五级层级作状态空间，建立带协变量生灭过程（$\lambda_i+\mu_i=c$ 恒定）；BFGS 极大似然 + 观测信息阵标准误；解析求首达时间 $m_k$，鞅停时定理严格化、Gillespie 模拟校验；一致化重建 $\bm P(t)$ 与谱隙；OU 扩散补连续效率层与时间尺度，用反射原理、反正弦律与 Itô 公式给闭式结果；最后在学时预算下做动作子集优化 | 三（三）（四） |

模型结构：$\rho_i=\exp(\alpha+\bm\delta^\top\bm z_i)$，$\alpha=-0.382$、$\bm\delta=(-0.135,-0.075,0.003,-0.272,0.305)$，$c=2.54$ 周$^{-1}$；五类规模 80/95/96/105/143，$\theta$ 均值 0.7935、与成绩相关 $r=0.675$。

## 主要结果

| 问题 | 结果 |
| --- | --- |
| 复现是否成立 | 结构性结论（五类、分层、路径方向）成立；数值性结论在锚点约束下成立（S480 16 周末已不在第 5 层概率 0.642、期望提升 25.53%，原文 25.38%） |
| 升学最快多久 | 首达最高层期望 14.48 周（Gillespie 14.50 周）；OU 反解累计 43.97% 需 53.7 周 |
| 只靠自然波动能否跃层 | 16 周内触及 0.95 的机会 1.5%（纯 Brown）；OU 下同目标 10.5%，两者不可混用 |
| 会不会先跌回去 | 从第 4 层出发先跌回第 5 层的概率 86% |
| 学时投在哪最划算 | 「兴趣与课外拓展 +1.5 h/周」（$\rho\to0.954$，提升 30.39%）；原文式「模仿 S477 +4.1 h/周」使 $\rho$ 降至 0.669 |

## 仓库结构

```
论文/
  paper_solo.tex           个人版入口（封面与页眉署「吴桐」）
  paper_team.tex           团队版入口（封面与页眉署四人）
  main.tex                 两版共用正文（版式、宏定义、章节 \input 顺序）
  paper_solo.pdf           个人版成品（49 页）
  paper_team.pdf           团队版成品（49 页）
  tex/                     章节源文件，由 main.tex 依次 \input
    00_cover_solo.tex      个人版封面
    00_cover_team.tex      团队版封面（授课教师：艾晓辉）
    01_abstract.tex        摘要与关键词
    02_problem.tex         一、问题陈述
    03_analysis.tex        二、问题分析
    04_parse.tex           三（一）原文解析
    05_repro.tex           三（二）CHPI 复现
    06_stochastic.tex      三（三）随机过程重构
    07_compare.tex         三（四）新旧方法系统对比与结论
    08_refs.tex            参考文献（GB/T 7714，15 条）
    09_appendix.tex        附录：符号表、六段关键推导、代码与数据说明
  fig/                     9 幅插图（PNG，300 dpi）
代码/
  dea_fun.R                DEA 线性规划工具函数（CCR、产出导向、反事实效率）
  01_recon_data.R          标定数据重建（set.seed(20261121)，二分法校准 S480）
  02_chpi.R                量化、PAM-Gower、CCR-DEA、分层、三准则与路径
  03_stochastic.R          生灭过程 MLE、首达时间、Gillespie、OU、策略优化
  04_figures.R             生成 9 幅插图并写入 论文/fig/
  05_digest.R              导出全部结果到 数据/digest.txt（核对用）
  06_textbook.R            谱隙、一致化误差、调和函数、反射原理、OU 蒙特卡洛、Itô 方差表
  07_export_web.R          把个体级结果导出为 数据/web_data.json（网页沙盘的数据源）
  build_web.py             把 JSON 内联进 网页/template.html，生成单文件 网页/index.html
数据/
  learning_data.csv        重建的 519×14 作答数据
  chpi_students.csv        量化得分、效率、层级、三准则
  sp_students.csv          协变量、ρ、平稳分布、首达时间
  chpi_results.rds / sp_results.rds   完整结果对象
  mc_check.csv             Gillespie 校验
  silhouette.csv           不同 k 的轮廓宽度
  path_s480.csv            S480 四步路径
  textbook_table.csv       谱隙与混合时间表
  web_data.json            网页用数据（meta/个体/系数/表格）
  digest.txt / chpi_log.txt / sp_log.txt / textbook_checks.txt   计算日志与全部数值
网页/
  template.html            沙盘模板（占位符 @@JSON@@）
  index.html               构建产物：单文件网页，双击可离线打开，也是部署到网站的那一份
查重/
  随机过程大作业_查重报告.pdf
  随机过程大作业_AIGC检测报告.pdf
```

## 复现

环境：R 4.6.0、TeX Live 2025（XeLaTeX）。R 依赖 `cluster`、`lpSolve`、`Matrix`、`numDeriv`、`ggplot2`、`showtext`、`sysfonts`。

**1. 编译论文**（在 `论文/` 目录下执行，两版各连续编译两遍，第一遍生成目录与交叉引用）

```bash
cd 论文
xelatex paper_solo.tex && xelatex paper_solo.tex
xelatex paper_team.tex && xelatex paper_team.tex
```

本机 `xelatex` 不在 PATH，用绝对路径：

```
& "D:\texlive\2025\bin\windows\xelatex.exe" -interaction=nonstopmode paper_team.tex
```

或一条命令：`& "D:\texlive\2025\bin\windows\latexmk.exe" -xelatex paper_team.tex`。插图路径由 `\graphicspath{{fig/}}` 指向 `论文/fig/`；两版共用 `main.tex`，版本差异由入口文件里的 `\def\soloversion` / `\def\teamversion` 切换封面与页眉。

**2. 重跑计算**（在仓库根目录执行）

```
& "D:\Program Files\R\R-4.6.0\bin\x64\Rscript.exe" 代码\01_recon_data.R
& "D:\Program Files\R\R-4.6.0\bin\x64\Rscript.exe" 代码\02_chpi.R
& "D:\Program Files\R\R-4.6.0\bin\x64\Rscript.exe" 代码\03_stochastic.R
& "D:\Program Files\R\R-4.6.0\bin\x64\Rscript.exe" 代码\04_figures.R
& "D:\Program Files\R\R-4.6.0\bin\x64\Rscript.exe" 代码\05_digest.R
& "D:\Program Files\R\R-4.6.0\bin\x64\Rscript.exe" 代码\06_textbook.R
```

若已把 `D:\Program Files\R\R-4.6.0\bin\x64` 加入 PATH，简写为 `Rscript 代码\01_recon_data.R` 等。

脚本开头的 `root` 自动定位仓库根：优先环境变量 `PROJ_ROOT`，其次当前工作目录，都不命中时回退到本机默认 `D:/AI/随机过程大作业2`。克隆到别的机器后 `cd` 到仓库根执行即可，也可先设 `$env:PROJ_ROOT = "<你的仓库路径>"`。

所有随机成分由固定种子控制（`set.seed(20261121)`；Gillespie 与 OU 蒙特卡洛各另有种子 20261122），重跑结果完全一致。注意重跑 `01` 会覆盖 `数据/` 下的数据文件；论文表格数字即来自这些文件，只想看数字不必重跑。

生成插图时用 `png(..., family="SimSun")` 并设 `las=1`（y 轴刻度横排），因此图在 Windows 上直接可用；Linux 无 SimSun 时把 `family` 换成 `"Noto Sans CJK SC"` 之类的可用中文字体即可。

**3. 重建网页**（在仓库根目录执行）

```
& "D:\Program Files\R\R-4.6.0\bin\x64\Rscript.exe" 代码\07_export_web.R
python 代码\build_web.py
```

产物是 `网页/index.html`：一个不含外部依赖的单文件，浏览器直接打开即可交互，`python -m http.server` 或任意静态托管都能部署。

## 数据与口径说明

1. **数据是「标定重建」而非「原始复原」**：原文未公开问卷、题项编码与量表，任何声称逐条复原的做法都不可信。本文以原文公开事实为锚点，构造统计结构同构的数据，并在论文中把「哪些是锚点、哪些是假设」逐条写清（三（二）1 与表 4、表 5）。因此复现结论分两级：结构性结论为「复现成立」，数值性结论为「在锚点约束下成立」。
2. **口径必须一致**：文中「16 周跃层概率 0.642」一律指「16 周末已不在第 5 层」的净位移概率，而非「期间从未跃层」（后者约 $5\times10^{-8}$）；「1.5%」只适用于纯 Brown 运动，OU 下同一目标的机会是 10.5%（Euler–Maruyama 模拟），两者不可混用。

## 团队

| 项目 | 内容 |
| --- | --- |
| 学院 | 理学院 |
| 专业班级 | 数学与应用数学 2024-01 / 2024-02 班 |
| 授课教师 | 艾晓辉 |
| 团队成员 | 吴桐、赵睿、张文琦、吴梦豪（2024 级） |
| 论文版本 | 个人版 `论文/paper_solo.pdf`、团队版 `论文/paper_team.pdf` |

## 许可证

按内容类型双许可，条款与分区说明见 [`LICENSING.md`](LICENSING.md)：

- **代码**（`代码/` 下全部 R 与 Python 脚本）：**MIT**，全文见 [`LICENSE`](LICENSE)，版权归**吴桐**个人所有。
- **论文、插图、网页与说明文字**：**CC BY-NC 4.0**，署名「吴桐、赵睿、张文琦、吴梦豪」，非商业使用、相同方式共享。
- 第三方原文（Zhao et al., ESWA 211 (2023) 118535）与「标定重建」数据不在本许可范围内。

## 参考文献

1. Charnes A, Cooper W W, Rhodes E. Measuring the efficiency of decision making units. *EJOR*, 1978.
2. Seiford L M, Zhu J. Context-dependent DEA—Measuring attractiveness and progress. *Omega*, 2003.
3. Zhao L T, et al. A recommendation system for effective learning strategies. *ESWA*, 2023.
4. Kalbfleisch J D, Lawless J F. The analysis of panel data under a Markov assumption. *JASA*, 1985.
5. Kendall D G. On the generalized birth-and-death process. *Ann. Math. Stat.*, 1948.
6. Gillespie D T. Exact stochastic simulation of coupled chemical reactions. *J. Phys. Chem.*, 1977.
7. Uhlenbeck G E, Ornstein L S. On the theory of the Brownian motion. *Phys. Rev.*, 1930.
8. 魏权龄. 数据包络分析. 科学出版社, 2004.

---

数据集与原文版权归原作者所有，仅用于课程学习。
