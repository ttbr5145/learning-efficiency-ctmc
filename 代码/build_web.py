#!/usr/bin/env python3
"""build_web.py —— 把 数据/web_data.json 内联进 网页/template.html 生成单文件网页。

网页是动态沙盘：所有曲线、分布、直方图、热图都由 JSON 里的个体数据在浏览器里重算，
不内联任何论文插图。用法（在仓库根目录）：

    python 代码/build_web.py

输出：网页/index.html（单文件，双击即可离线打开，也是部署到网站的那一份）
"""
import json
import pathlib
import sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
TPL = ROOT / "网页" / "template.html"
OUT = ROOT / "网页" / "index.html"
DATA = ROOT / "数据" / "web_data.json"


def main() -> int:
    if not TPL.exists():
        print(f"缺少模板：{TPL}", file=sys.stderr)
        return 1
    if not DATA.exists():
        print(f"缺少数据：{DATA}（先运行 代码/07_export_web.R）", file=sys.stderr)
        return 1

    payload = DATA.read_text(encoding="utf-8")
    json.loads(payload)  # 提前校验，避免把坏 JSON 塞进页面

    html = TPL.read_text(encoding="utf-8")
    if "@@JSON@@" not in html:
        print("模板里没有 @@JSON@@ 占位符", file=sys.stderr)
        return 1
    html = html.replace("@@JSON@@", payload)

    OUT.write_text(html, encoding="utf-8")
    print(f"写出 {OUT.relative_to(ROOT)}（{OUT.stat().st_size/1024:.1f} KB）")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
