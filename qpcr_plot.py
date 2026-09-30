# -*- coding: utf-8 -*-
"""
qPCR 相对定量计算 + 柱形图脚本
================================
输入 Ct 值，自动完成:
    1. 每个样品每个基因的 Ct 均值 / SEM（技术重复用列表表示）
    2. dCt  = Ct(目的基因) - Ct(内参)
    3. ddCt = dCt(样品) - dCt(对照组)
    4. RQ   = 2 ^ (-ddCt)          <- 最终倍数结果
    5. 误差线: 把技术重复的 SEM 通过误差传播算到 RQ 上
    6. 输出柱形图（PNG / SVG / PDF），样式接近 Auto-qPCR

用法一（最简单）: 直接改下面 DATA 区域的数值，然后运行
    python qpcr_plot.py

用法二（命令行）: 用长表 CSV(sample,gene,ct) 跑
    python qpcr_plot.py --csv data.csv --ref "Actin" --control "0H"

依赖: numpy, matplotlib
"""

from __future__ import annotations

import argparse
import csv
import math
import os
import sys
from typing import Dict, Iterable, List, Sequence, Tuple

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


# ======================================================================
# 1.  DATA —— 在这里填你的 Ct 值
#    结构: { "样品名": {"基因名": [重复1, 重复2, ...], ...}, ... }
#    只输入一个值也可以（那样不画误差线）；多个技术重复就多写几个数。
#    内参基因（如 Actin）也要作为其中一个"基因"写进去。
# ======================================================================

REFERENCE_GENE = "Actin"        # 内参基因名（必须出现在下面的基因里）
CONTROL_SAMPLE = "0H"           # 对照样品名（作为 ddCt 的基准，RQ = 1）

DATA: Dict[str, Dict[str, List[float]]] = {
    "0H":  {"CNNM1": [22.48], "Actin": [12.88]},
    "12H": {"CNNM1": [22.35], "Actin": [13.81]},
    "24H": {"CNNM1": [22.18], "Actin": [12.99]},
}

# 画图外观设置 ------------------------------------------------------
CHART_TITLE = "Crotonic acid treatment (A549)"   # 图标题，不需要就设为 ""
X_LABEL = ""                                     # 横轴标题
Y_LABEL = "Relative CNNM1 mRNA level (RQ)"       # 纵轴标题
SHOW_GENE_CHART = None                           # True/False/None(自动)
OUTPUT_PREFIX = "qpcr_result"                    # 输出文件名前缀
OUTPUT_DIR = os.path.dirname(os.path.abspath(__file__))


# ======================================================================
# 2.  计算部分
# ======================================================================

def sem(values: Sequence[float]) -> float:
    """样本标准误 SEM（样本标准差 / sqrt(n)），不足 2 个点返回 0。"""
    arr = np.asarray(list(values), dtype=float)
    arr = arr[~np.isnan(arr)]
    if arr.size < 2:
        return 0.0
    return float(arr.std(ddof=1) / math.sqrt(arr.size))


def mean(values: Sequence[float]) -> float:
    arr = np.asarray(list(values), dtype=float)
    arr = arr[~np.isnan(arr)]
    return float(arr.mean()) if arr.size else float("nan")


class QpcrResult:
    """一个样品、一个基因的计算结果。"""

    __slots__ = ("sample", "gene", "ct", "ct_sem", "dct", "dct_sem",
                 "ddct", "ddct_sem", "rq", "rq_err_low", "rq_err_up")

    def __init__(self, sample: str, gene: str):
        self.sample = sample
        self.gene = gene
        self.ct = float("nan")
        self.ct_sem = 0.0
        self.dct = float("nan")
        self.dct_sem = 0.0
        self.ddct = float("nan")
        self.ddct_sem = 0.0
        self.rq = float("nan")
        self.rq_err_low = 0.0
        self.rq_err_up = 0.0


def compute(data: Dict[str, Dict[str, List[float]]],
            reference: str,
            control: str) -> Tuple[List[str], List[str], List[QpcrResult]]:
    """核心计算，返回 (样品名列表, 目的基因列表, 结果列表)。"""
    samples = list(data.keys())
    if control not in data:
        raise KeyError(f"对照样品 {control!r} 不在数据里，现有样品: {samples}")

    genes = list(data[control].keys())
    if reference not in genes:
        raise KeyError(f"内参基因 {reference!r} 不在基因列表里: {genes}")
    targets = [g for g in genes if g != reference]
    if not targets:
        raise ValueError("除了内参基因以外没有别的目的基因，无法计算。")

    results: Dict[Tuple[str, str], QpcrResult] = {}

    # 第一步: Ct 均值 / SEM / dCt
    for s in samples:
        if reference not in data[s]:
            raise KeyError(f"样品 {s!r} 缺少内参基因 {reference!r}")
        ref_ct = mean(data[s][reference])
        ref_sem = sem(data[s][reference])
        for g in targets:
            r = QpcrResult(s, g)
            r.ct = mean(data[s][g])
            r.ct_sem = sem(data[s][g])
            r.dct = r.ct - ref_ct
            r.dct_sem = math.hypot(r.ct_sem, ref_sem)
            results[(s, g)] = r

    # 第二步: 对照组平均 dCt
    ctrl_dct: Dict[str, float] = {}
    ctrl_dct_sem: Dict[str, float] = {}
    for g in targets:
        ctrl_dct[g] = results[(control, g)].dct
        ctrl_dct_sem[g] = results[(control, g)].dct_sem

    # 第三步: ddCt 与 RQ
    for (s, g), r in results.items():
        r.ddct = r.dct - ctrl_dct[g]
        r.ddct_sem = math.hypot(r.dct_sem, ctrl_dct_sem[g])
        r.rq = 2.0 ** (-r.ddct)
        # RQ 的 1xSEM 区间: RQ * 2^(+/-ddCt_sem) - RQ
        r.rq_err_up = max(r.rq * 2.0 ** r.ddct_sem - r.rq, 0.0)
        r.rq_err_low = max(min(r.rq - r.rq * 2.0 ** (-r.ddct_sem), r.rq), 0.0)

    ordered = [results[(s, g)] for s in samples for g in targets]
    return samples, targets, ordered


def print_table(samples: List[str], targets: List[str], results: List[QpcrResult]) -> None:
    lookup = {(r.sample, r.gene): r for r in results}
    print("\n=== qPCR 计算结果 ===")
    header = f"{'样品':<16}{'基因':<14}{'Ct':>8}{'dCt':>9}{'ddCt':>9}{'RQ(2^-ddCt)':>14}"
    print(header)
    print("-" * len(header))
    for s in samples:
        for g in targets:
            r = lookup[(s, g)]
            print(f"{s:<16}{g:<14}{r.ct:>8.2f}{r.dct:>9.2f}{r.ddct:>9.2f}{r.rq:>14.3f}")
    print()


# ======================================================================
# 3.  画图部分（样式参考 Auto-qPCR 的 matplotlib 输出）
# ======================================================================

PALETTE = ["#4e79a7", "#f28e2b", "#e15759", "#76b7b2", "#59a14f",
           "#edc948", "#b07aa1", "#ff9da7", "#9c755f", "#bab0ac"]


def setup_chinese_font() -> None:
    """让图中中文正常显示，找不到中文字体就退回默认。"""
    from matplotlib import font_manager
    available = {f.name for f in font_manager.fontManager.ttflist}
    for name in ["Microsoft YaHei", "SimHei", "Noto Sans CJK SC",
                 "Source Han Sans SC", "PingFang SC", "Arial Unicode MS"]:
        if name in available:
            plt.rcParams["font.sans-serif"] = [name]
            break
    plt.rcParams["axes.unicode_minus"] = False


def draw_bar_chart(ax: plt.Axes,
                   series: List[dict],
                   group_names: Sequence[str],
                   y_label: str,
                   x_label: str = "",
                   title: str = "") -> None:
    """在 ax 上画分组柱形图。

    series: [{"name": 图例名, "color": 颜色, "values": [QpcrResult, ...]}, ...]
    每个 value 与 group_names 一一对应。
    """
    n_series = len(series)
    n_groups = len(group_names)
    x = np.arange(n_groups)
    total_width = 0.8 if n_series > 1 else 0.55
    bar_w = total_width / max(n_series, 1)

    for si, s in enumerate(series):
        offset = (si - (n_series - 1) / 2.0) * bar_w
        heights = np.array([v.rq for v in s["values"]], dtype=float)
        err_low = np.array([v.rq_err_low for v in s["values"]], dtype=float)
        err_up = np.array([v.rq_err_up for v in s["values"]], dtype=float)
        xpos = x + offset

        ax.bar(xpos, heights, width=bar_w * 0.92, color=s["color"],
               edgecolor="white", linewidth=0.6, label=s["name"], zorder=3)

        if np.any(err_up > 0) or np.any(err_low > 0):
            ax.errorbar(xpos, heights, yerr=np.vstack([err_low, err_up]),
                        fmt="none", ecolor="#333333", elinewidth=1.2,
                        capsize=3.5, capthick=1.2, zorder=4)

    # 参考线 RQ = 1
    ax.axhline(1.0, color="#999999", linewidth=1.0,
               linestyle=(0, (4, 3)), zorder=1)

    ax.set_xticks(x)
    ax.set_xticklabels([str(g) for g in group_names])
    ax.set_ylabel(y_label)
    if x_label:
        ax.set_xlabel(x_label)
    if title:
        ax.set_title(title, fontsize=13, pad=12)

    # 只要左轴和下轴，去掉上、右框线
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_linewidth(1.1)
        ax.spines[side].set_color("#333333")
    ax.tick_params(direction="out", length=4, width=1.0, colors="#333333")

    peak = float(np.nanmax([v.rq + v.rq_err_up
                            for s in series for v in s["values"]]))
    ax.set_ylim(0, max(1.05, peak * 1.18))
    ax.grid(axis="y", color="#dddddd", linewidth=0.7, zorder=0)
    ax.set_axisbelow(True)

    if n_series > 1:
        ax.legend(frameon=False, loc="upper right", fontsize=10)


def make_figures(samples: List[str],
                 targets: List[str],
                 results: List[QpcrResult],
                 args) -> List[Tuple[str, plt.Figure]]:
    """生成一到两张图: 按样品分组（每个目的基因一根柱），以及按基因分组。"""
    lookup = {(r.sample, r.gene): r for r in results}
    figures: List[Tuple[str, plt.Figure]] = []

    series_by_gene = [
        {"name": g,
         "color": PALETTE[gi % len(PALETTE)],
         "values": [lookup[(s, g)] for s in samples]}
        for gi, g in enumerate(targets)
    ]

    fig1, ax1 = plt.subplots(figsize=(6.4, 4.6))
    draw_bar_chart(ax1, series_by_gene, samples,
                   y_label=args.ylabel, x_label=args.xlabel, title=args.title)
    fig1.tight_layout()
    figures.append(("by_sample", fig1))

    show_gene_chart = args.gene_chart
    if show_gene_chart is None:
        show_gene_chart = len(targets) > 1
    if show_gene_chart and len(targets) > 1:
        series_by_sample = [
            {"name": s,
             "color": PALETTE[ri % len(PALETTE)],
             "values": [lookup[(s, g)] for g in targets]}
            for ri, s in enumerate(samples)
        ]
        fig2, ax2 = plt.subplots(figsize=(6.4, 4.6))
        draw_bar_chart(ax2, series_by_sample, targets,
                       y_label=args.ylabel, x_label=args.xlabel,
                       title=(args.title + " (by gene)") if args.title else "")
        fig2.tight_layout()
        figures.append(("by_gene", fig2))

    return figures


# ======================================================================
# 4.  读 CSV（可选）
# ======================================================================

def load_csv(path: str) -> Dict[str, Dict[str, List[float]]]:
    """长表 CSV，表头形如: sample,gene,ct  （一行一个技术重复）。"""
    out: Dict[str, Dict[str, List[float]]] = {}
    with open(path, "r", encoding="utf-8-sig", newline="") as fh:
        reader = csv.DictReader(fh)
        cols = {c.strip().lower(): c for c in (reader.fieldnames or [])}
        need = ["sample", "gene", "ct"]
        if not all(k in cols for k in need):
            raise ValueError(f"CSV 需要 sample/gene/ct 三列，实际表头: {reader.fieldnames}")
        for row in reader:
            s = (row[cols["sample"]] or "").strip()
            g = (row[cols["gene"]] or "").strip()
            raw = (row[cols["ct"]] or "").strip()
            if not s or not g or not raw:
                continue
            out.setdefault(s, {}).setdefault(g, []).append(float(raw))
    return out


# ======================================================================
# 5.  主程序
# ======================================================================

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="qPCR 2^-ddCt 计算与柱形图（默认使用脚本内 DATA）",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--csv", help="长表 CSV 文件 (sample,gene,ct)；提供后忽略脚本内 DATA")
    p.add_argument("--ref", default=REFERENCE_GENE, help="内参基因名")
    p.add_argument("--control", default=CONTROL_SAMPLE, help="对照样品名")
    p.add_argument("--title", default=CHART_TITLE, help="图标题，留空则不显示")
    p.add_argument("--xlabel", default=X_LABEL, help="横轴标题")
    p.add_argument("--ylabel", default=Y_LABEL, help="纵轴标题")
    p.add_argument("--out", default=OUTPUT_PREFIX, help="输出文件名前缀")
    p.add_argument("--outdir", default=OUTPUT_DIR, help="输出目录")
    p.add_argument("--formats", default="png,svg", help="输出格式，逗号分隔: png,svg,pdf")
    p.add_argument("--dpi", type=int, default=600, help="位图分辨率")
    p.add_argument("--gene-chart", dest="gene_chart", action="store_true",
                   default=SHOW_GENE_CHART, help="额外输出按基因分组的图")
    p.add_argument("--no-gene-chart", dest="gene_chart", action="store_false",
                   help="不输出按基因分组的图")
    return p


def main(argv: Iterable[str] | None = None) -> int:
    args = build_parser().parse_args(list(argv) if argv is not None else None)

    data = load_csv(args.csv) if args.csv else DATA

    samples, targets, results = compute(data, reference=args.ref, control=args.control)
    print_table(samples, targets, results)

    setup_chinese_font()
    figures = make_figures(samples, targets, results, args)

    os.makedirs(args.outdir, exist_ok=True)
    formats = [f.strip().lower() for f in args.formats.split(",") if f.strip()]
    written = []
    for tag, fig in figures:
        for ext in formats:
            path = os.path.join(args.outdir, f"{args.out}_{tag}.{ext}")
            fig.savefig(path, dpi=args.dpi, bbox_inches="tight", facecolor="white")
            written.append(path)
        plt.close(fig)

    print("已输出图片:")
    for path in written:
        print("  " + path)
    return 0


if __name__ == "__main__":
    sys.exit(main())
