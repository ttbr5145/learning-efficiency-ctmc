# 授权说明 / Licensing

本仓库按**内容类型**分别授权。`LICENSE` 文件只放代码用的 MIT 全文（便于 GitHub 自动识别），
其余内容的分区与条款写在本文件里。

| 内容 | 路径 | 许可证 | 版权人 |
| --- | --- | --- | --- |
| 代码 | `代码/`（R 脚本、Python 脚本、构建脚本） | [MIT](./LICENSE) | 吴桐 |
| 论文、插图、网页、文字 | `论文/`、`网页/`、`README.md`、`数据/web_data.json` 的字段说明 | [CC BY-NC 4.0](https://creativecommons.org/licenses/by-nc/4.0/legalcode.zh-hans) | 吴桐、赵睿、张文琦、吴梦豪 |
| 第三方原文与数据 | 引文 Zhao et al. (2023)；原始问卷 | 版权归原作者，仅用于课程学习 | — |

---

## 1. 代码：MIT License（版权归吴桐个人所有）

`代码/` 目录下全部 R 与 Python 脚本（含 `01_recon_data.R` … `07_export_web.R`、`build_web.py`）
以 **MIT License** 授权，全文见仓库根目录的 [`LICENSE`](./LICENSE)。

```
Copyright (c) 2026 吴桐
```

这是一份宽松许可：可以自由使用、修改、合并、发布、再许可乃至商业使用，
只需在副本或实质部分中保留上述版权声明与许可声明。软件按“原样”提供，不含任何担保。

论文、插图与网页等**非代码内容**不适用 MIT，其署名仍是四位作者（见第 2 节）。

## 2. 论文、插图、网页与文字：CC BY-NC 4.0

本仓库中的论文（`论文/paper_solo.pdf`、`论文/paper_team.pdf` 及其 LaTeX 源码）、
插图（`论文/fig/`）、网页（`网页/`）与说明文字，采用
**Creative Commons Attribution-NonCommercial 4.0 International（CC BY-NC 4.0）** 许可。

你可以自由复制、转载、改编并在非商业场景下使用，但须遵守三条：

1. **署名**：注明作者「吴桐、赵睿、张文琦、吴梦豪」与本仓库地址；
2. **非商业**：不得用于商业目的；
3. **相同方式共享**：改编后须以相同许可发布，并说明是否作出修改。

许可全文（中文）：<https://creativecommons.org/licenses/by-nc/4.0/legalcode.zh-hans>

## 3. 第三方内容

- 被解析的原文：Zhao L T, et al. *A recommendation system for effective learning
  strategies: An integrated approach using context-dependent DEA*. Expert Systems
  with Applications, 2023, 211: 118535。版权归 Elsevier 及原作者。
- 原始问卷数据未公开；本仓库 `数据/learning_data.csv` 是依公开锚点**标定重建**的
  同构数据，不代表任何个体的真实作答。

上述内容不在本仓库许可范围内，引用请遵循原出版方的条款。

## 4. 免责

本项目为东北林业大学《随机过程》课程作业，模型结论仅用于教学与方法演示，
不构成对任何个人学习能力、学业成绩的评价或预测依据。
