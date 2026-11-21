# 《随机过程》大作业 · 学习效率层级演化的随机过程建模与策略推荐

——对一篇 context-dependent DEA 论文的解析、复现与方法替换

作者：赵睿（理学院 数学与应用数学 2024 级 01 班，学号 2024213044）
授课教师：艾晓辉
篇幅：49 页（封面、摘要、目录、四节正文、参考文献、附录），9 幅插图、21 张表、57 个编号公式。

仓库地址：<https://github.com/ttbr5145/learning-efficiency-ctmc>

```bash
git clone https://github.com/ttbr5145/learning-efficiency-ctmc.git
```

---

## 一、怎么编译（TeXworks）

1. 打开 TeXworks，用「文件 → 打开」选择本文件夹下的 `main.tex`。
2. 在工具栏的**编译引擎下拉框**里选 **XeLaTeX**（不要用 pdfLaTeX，中文会报错）。
3. 点绿色三角编译。**请连续编译两遍**：第一遍生成目录与交叉引用，第二遍才把页码和「式 (n)」「表 n」「图 n」的编号对齐。
4. 生成的 `main.pdf` 就在同一目录，可直接阅读。

命令行等价做法（在本文件夹内打开 PowerShell）：

```
xelatex main.tex
xelatex main.tex
```

或一条命令搞定：

```
latexmk -xelatex main.tex
```

若本机 `xelatex` 不在 PATH，用 TeX Live 的绝对路径（本机实测可用）：

```
& "D:\texlive\2025\bin\windows\xelatex.exe" -interaction=nonstopmode main.tex
```

编译产生的 `main.aux`、`main.log`、`main.out`、`main.toc`、`main.synctex.gz` 属中间产物，已在 `.gitignore` 中忽略，不会进入仓库。

> 依赖宏包均为 TeX Live 自带：ctex、xeCJK、geometry、amsmath、amssymb、amsthm、bm、booktabs、tabularx、multirow、array、longtable、graphicx、float、caption、enumitem、xcolor、hyperref。无需额外安装。

---

## 二、论文做了什么（三部分）

| 部分 | 内容 | 对应章节 |
|---|---|---|
| 第一步 破解 | 把原文（Zhao et al., *Expert Systems with Applications* 211 (2023) 118535）的五环证据链拆成可检验的数学表述：复合量化 → Gower+PAM 聚类 → 类内 CCR-DEA → context-dependent 分层 → 三准则与参考路径；逐环标出原文未披露的自由度与内部不一致（含 $H=1/\Omega+P$ 与表中数值不自洽一处） | 三（一） |
| 第二步 复现 | 在原始问卷不可得的条件下，用四条锚点（$n=519$、5 类、S480 效率 0.6771、单步 25.38%／累计 43.97%）逆向重建 519 份样本，用 R 完整复现 CHPI：PAM-Gower 聚类（轮廓系数 $k=3,\dots,8$ 全部报告、与设计标签一致率 93.6%）、CCR-DEA（lpSolve 求解，包络式与乘子式互验）、分位带与剔除式两种分层口径、三准则与 S480 四步路径 | 三（二） |
| 第三步 换方法 | 弃用 context-dependent DEA，改用随机过程：把五级层级作为状态空间，建立带协变量的生灭过程（强度比 $\rho_i=\exp(\alpha+\bm{\delta}^{\top}\bm{z}_i)$，$\lambda_i+\mu_i=c$ 恒定），BFGS 极大似然 + 观测信息阵标准误；解析求首达时间 $m_k$，并用鞅的停时定理给出严格化、Gillespie 模拟校验；用一致化重建 $\bm{P}(t)$ 与谱隙；用 Ornstein–Uhlenbeck 扩散补出连续效率层与时间尺度，以反射原理、反正弦律与 Itô 公式给出闭式结果；最后在学时预算下做动作子集优化，并与原文逐项对比 | 三（三）（四） |

关键数字（全部由 `code/` 下的 R 脚本算出，可在 `data/` 里核对）：类规模 80/95/96/105/143；$\theta$ 均值 0.7935，与成绩相关 $r=0.675$；$\alpha=-0.382$，$\delta=(-0.135,-0.075,0.003,-0.272,0.305)$；$\rho$ 中位数 0.682，98.1% 的个体 $\rho<1$；$c=2.54$ 周$^{-1}$；S480 的 16 周末已不在第 5 层的概率 0.642、期望效率提升 25.53%（原文 25.38%）；首达最高层期望 14.48 周（Gillespie 14.50 周）；OU 反解累计 43.97% 需 53.7 周；仅靠自然波动 16 周内触及 0.95 的机会 1.5%（OU 下 10.5%）；从第 4 层出发先跌回第 5 层的概率 86%。最优策略为「兴趣与课外拓展 +1.5 h/周」（$\rho$ 升至 0.954，提升 30.39%），而原文式「模仿 S477、每周多投入 4.1 h」反而使 $\rho$ 降至 0.669。

---

## 三、目录结构

仓库根目录：<https://github.com/ttbr5145/learning-efficiency-ctmc>（下文所有路径均相对仓库根）

```
learning-efficiency-ctmc/
├── main.tex                 论文主文件（XeLaTeX 编译入口，版式与宏定义）
├── main.pdf                 已编译成品（49 页）
├── README.md                本说明
├── tex/                     各章节源文件，由 main.tex 依次 \input
│   ├── 00_cover.tex         封面（按《随机过程》大作业模板，授课教师：艾晓辉）
│   ├── 01_abstract.tex      摘要与关键词
│   ├── 02_problem.tex       一、问题陈述
│   ├── 03_analysis.tex      二、问题分析
│   ├── 04_parse.tex         三（一）原文解析
│   ├── 05_repro.tex         三（二）CHPI 复现
│   ├── 06_stochastic.tex    三（三）随机过程重构
│   ├── 07_compare.tex       三（四）新旧方法系统对比与结论
│   ├── 08_refs.tex          参考文献（GB/T 7714，15 条）
│   └── 09_appendix.tex      附录：符号表、六段关键推导、代码与数据说明
├── code/                    R 脚本（按编号顺序执行即可完整复现）
│   ├── dea_fun.R            DEA 线性规划工具函数（CCR、产出导向、反事实效率）
│   ├── 01_recon_data.R      标定数据重建（set.seed(20261121)，二分法校准 S480）
│   ├── 02_chpi.R            量化、PAM-Gower、CCR-DEA、分层、三准则与路径
│   ├── 03_stochastic.R      生灭过程 MLE、首达时间、Gillespie、OU、策略优化
│   ├── 04_figures.R         生成 9 幅插图（PNG，300 dpi，SimSun 字体）
│   ├── 05_digest.R          导出全部结果到 data/digest.txt（核对用）
│   └── 06_textbook.R        谱隙、一致化误差、调和函数、反射原理、OU 蒙特卡洛、Itô 方差表
├── data/                    中间数据与结果（CSV / RDS / 日志）
│   ├── learning_data.csv    重建的 519×14 作答数据
│   ├── chpi_students.csv    量化得分、效率、层级、三准则
│   ├── sp_students.csv      协变量、ρ、平稳分布、首达时间
│   ├── chpi_results.rds / sp_results.rds   完整结果对象
│   ├── mc_check.csv         Gillespie 校验
│   ├── silhouette.csv       不同 k 的轮廓宽度
│   ├── path_s480.csv        S480 四步路径
│   ├── textbook_table.csv   谱隙与混合时间表
│   └── digest.txt / chpi_log.txt / sp_log.txt / textbook_checks.txt   计算日志与全部数值
├── figures/                 9 幅插图（PNG，300 dpi）
├── .gitignore               忽略编译中间产物与 R 临时文件
├── 免费_PDF打印版_查重报告_[main] (1).pdf      查重报告（公开保留）
└── 免费_PDF打印版_AIGC检测报告_[main].pdf      AIGC 检测报告（公开保留）
```

## 四、怎么复现（R）

R 版本：R 4.6.0。依赖 `cluster`、`lpSolve`、`Matrix`、`numDeriv`、`ggplot2`、`showtext`、`sysfonts`。

在本仓库根目录依次执行（本机 `Rscript` 未在 PATH，故先给出绝对路径版本）：

```
& "D:\Program Files\R\R-4.6.0\bin\x64\Rscript.exe" code\01_recon_data.R
& "D:\Program Files\R\R-4.6.0\bin\x64\Rscript.exe" code\02_chpi.R
& "D:\Program Files\R\R-4.6.0\bin\x64\Rscript.exe" code\03_stochastic.R
& "D:\Program Files\R\R-4.6.0\bin\x64\Rscript.exe" code\04_figures.R
& "D:\Program Files\R\R-4.6.0\bin\x64\Rscript.exe" code\05_digest.R
& "D:\Program Files\R\R-4.6.0\bin\x64\Rscript.exe" code\06_textbook.R
```

若已把 `D:\Program Files\R\R-4.6.0\bin\x64` 加入 PATH，则简写为 `Rscript code\01_recon_data.R` 等。

脚本开头的 `root` 会自动定位仓库根目录：优先读环境变量 `PROJ_ROOT`，其次用当前工作目录（即仓库根），两者都不命中时回退到本机默认路径 `D:/AI/随机过程大作业2`。所以克隆到别的机器后，只要`cd` 到仓库根目录再执行即可；也可先设 `$env:PROJ_ROOT = "<你的仓库路径>"`。

每一步都会把结果写入 `data/` 并打印关键数值。所有随机成分都由固定种子控制（`set.seed(20261121)`；Gillespie 与 OU 蒙特卡洛各另有独立种子 20261122），重跑结果完全一致。注意：重跑 01 会覆盖 `data/` 下的数据文件；论文里的表格数字即来自这些文件，若只想看数字不必重跑。

`06_textbook.R` 用于核对正文中每一处公式与模拟的对接触发的数字：一致化重建 $\bm{P}(16)$ 与 `expm` 的最大误差（$3.8\times10^{-15}$）、五类谱隙与混合时间、S480 的调和函数与「先到第 1 层」概率、纯 Brown 运动的反射原理、OU 的 Euler–Maruyama 蒙特卡洛（2000 条路径、步长 0.005 周）、Itô 公式给出的 16 周与平稳标准差。

生成插图时用 `png(..., family="SimSun")` 并设 `las=1`（y 轴刻度横排），因此图在 Windows 上直接可用；Linux 无 SimSun 时把 `family` 换成 `"Noto Sans CJK SC"` 之类的可用中文字体即可。

---

## 五、两点说明

1. **数据是「标定重建」而非「原始复原」**：原文没有公开问卷、题项编码与量表，任何声称逐条复原的做法都不可信。本文用原文公开的事实作锚点，构造一份统计结构同构的数据，并在论文中把「哪些是锚点、哪些是假设」逐条写清楚（第三节（二）1 与表 4、表 5）。因此复现结论分两级：结构性结论（五类、分层、路径方向）为「复现成立」，数值性结论为「在锚点约束下成立」。
2. **口径必须一致**：文中「16 周跃层概率 0.642」一律指「16 周末已不在第 5 层」的净位移概率，而非「期间从未跃层」（后者约 $5\times10^{-8}$）；「1.5%」只适用于纯 Brown 运动，OU 下同一目标的机会是 10.5%（Euler–Maruyama 模拟），两者不可混用——这一条同时是正文的方法论提醒。
