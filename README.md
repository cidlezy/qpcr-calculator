# qPCR Calculator

计算QPCR的简易小工具：输入 Ct 值，自动计算 ΔCt、ΔΔCt 和相对表达量 (2^-ΔΔCt)，并可直接出图。

纯静态页面，所有计算在浏览器本地完成，不上传任何数据。直接打开 `index.html` 即可使用。

在线版：<https://cidlezy.github.io/qpcr-calculator/>

## 功能

- 每行输入样品名、内参 Ct、目的基因 Ct，实时计算结果（目的基因名和数量可自选）
- 每行可选对照样品，以该行的 ΔCt 作为基准
- 支持在 Ct 输入框中输入多个技术重复（逗号分隔），自动计算 SEM 并绘制带误差线的柱形图
- 自动生成"按样品分组"和"按基因分组"两张柱形图，风格仿照 Auto-qPCR 的 matplotlib 输出
- 横轴标签和图标题可编辑，图表可保存为 PNG / SVG
- **用 Python 出图**：在页面内直接运行本仓库的 `qpcr_plot.py`（浏览器端 Pyodide + matplotlib），得到与本地脚本完全一致的图；结果可下载 300 dpi PNG 和矢量 SVG，数据可导出 CSV
- 结果表可一键复制为制表符分隔格式，直接粘贴到 Excel
- 支持打印或导出 PDF

## 本地跑 Python 脚本（可选）

```bash
python qpcr_plot.py --csv qpcr_data.csv
```

CSV 格式为 `sample,gene,ct` 三列，可用页面上的"导出数据 CSV"按钮直接生成。网页里的 Python 环境不含中文字体，图内文字请用英文或数字。

## 计算公式

```
ΔCt = Ct(目的基因) - Ct(内参基因)
ΔΔCt = ΔCt(处理组) - ΔCt(对照组)
Fold Change = 2^(-ΔΔCt)
```

前提是目的基因和内参基因的扩增效率都接近 100%（效率 E ≈ 2）。
