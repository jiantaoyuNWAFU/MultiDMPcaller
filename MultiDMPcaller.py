import sys
import argparse
import os
import re
import pandas as pd
import numpy as np
from scipy.stats import fisher_exact
import glob
import time
from collections import defaultdict
from pathlib import Path
from typing import Union, Tuple, Optional, List
import matplotlib
import matplotlib.pyplot as plt
from scipy.special import betainc
from sklearn.mixture import GaussianMixture
import bisect
import subprocess
import shutil
import contextlib
from concurrent.futures import ProcessPoolExecutor, as_completed
from collections import namedtuple
import math
import matplotlib.ticker as ticker
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib import colors as mcolors
from scipy.ndimage import gaussian_filter1d

# DMP的Q值显著性阈值配置
DMP_QVALUE_THRESHOLDS = {
    'CpG': 0.05,
    'CHH': 0.045,
    'CHG': 0.04,
}

# DMR的Q值显著性阈值
DMR_QVALUE_THRESHOLD = 0.05

# 最终DMP/DMR跨replicate组合投票阈值；默认2/3，兼容原代码行为
VOTE_THRESHOLD = 2 / 3

# 自动估计 DMP/DMR 投票阈值。默认关闭，保证不传新参数时完全沿用 --vote-threshold。
AUTO_DMP_VOTE_THRESHOLD = False
AUTO_DMR_VOTE_THRESHOLD = False
AUTO_VOTE_THRESHOLD_REPORT_ONLY = False
AUTO_DMP_VOTE_THRESHOLDS = {'CpG': None, 'CHH': None, 'CHG': None}
AUTO_DMR_VOTE_THRESHOLDS = {'CpG': None, 'CHH': None, 'CHG': None}

# DMR candidate generation engine: 'python' keeps the original implementation; 'cpp' uses dmr_step1 + dmr_step2_dynamic.
DMR_ENGINE = "python"

# Adaptive q-value threshold for two-step FDR contexts only.
# Disabled by default to keep the original behavior unchanged.
AUTO_QVALUE_TWOSTEP = False
AUTO_QVALUE_REPORT_ONLY = False
AUTO_QVALUE_P_CUTOFF = 0.05
AUTO_QVALUE_MIN_CANDIDATES = 10
AUTO_QVALUE_USE_SMOOTH = False
AUTO_QVALUE_SMOOTH_SIGMA = 4.0

# Final-DMP low-difference strict voting post-filter.
# Disabled by default to keep legacy behavior unchanged. When enabled, the
# provisional final DMPs are first called by the usual q-value + vote rule;
# candidates whose boundary abs(MethDiff) is <= DMP_LOWDIFF_CUTOFF must then
# satisfy a stricter vote requirement: ceil((base_required + total_pairs) / 2).
DMP_LOWDIFF_STRICT_VOTE = False
DMP_LOWDIFF_CUTOFF = 0.3
DMP_LOWDIFF_STRICT_VOTE_REPORT_ONLY = False

def get_dmp_threshold(methylation_type):
    """获取指定甲基化类型的DMP Q值阈值"""
    return DMP_QVALUE_THRESHOLDS.get(methylation_type, 0.05)


def calc_meth_diff(m1, u1, m2, u2):
    """计算两个样本之间的绝对甲基化水平差异，范围为0-1。"""
    total1 = m1 + u1
    total2 = m2 + u2
    ratio1 = m1 / total1 if total1 > 0 else 0
    ratio2 = m2 / total2 if total2 > 0 else 0
    return abs(ratio1 - ratio2)


# =========================================================
# Mianjifa auto-methdiff module (embedded, no GUI)
# Adapted from mianjifa_automethdiff.py. It only:
#   1) estimates one global abs(methdiff) threshold from pairwise DMP distributions;
#   2) saves m×n MethDiff distribution plots for each methylation context.
# It does not modify pairwise Fisher/FDR logic.
# =========================================================
AUTO_METHDIFF_X_MIN = -1.0
AUTO_METHDIFF_X_MAX = 1.0
AUTO_METHDIFF_BIN_WIDTH = 0.01
AUTO_METHDIFF_N_BINS = int(round((AUTO_METHDIFF_X_MAX - AUTO_METHDIFF_X_MIN) / AUTO_METHDIFF_BIN_WIDTH))
AUTO_METHDIFF_BINS = np.linspace(AUTO_METHDIFF_X_MIN, AUTO_METHDIFF_X_MAX, AUTO_METHDIFF_N_BINS + 1)
AUTO_METHDIFF_COLOR_NORMAL = "#4C72B0"
AUTO_METHDIFF_COLOR_CUT = "#DD8452"
AUTO_METHDIFF_COLOR_ZERO_LINE = "black"
AUTO_METHDIFF_COLOR_THRESHOLD = "red"
AUTO_METHDIFF_RAW_CACHE = {}


def _mianjifa_safe_name(name: str) -> str:
    """Generate a filesystem-safe filename component."""
    return re.sub(r'[\\/:*?"<>|]+', "_", str(name))


def _mianjifa_normalize_chr(chr_value):
    """Normalize chromosome labels: chr1/Chr1/1 -> 1."""
    chr_value = str(chr_value).strip()
    chr_value = re.sub(r"^chr", "", chr_value, flags=re.IGNORECASE)
    return chr_value


def _mianjifa_normalize_context(context_value):
    """Normalize methylation contexts, allowing CG and CpG to match."""
    context_value = str(context_value).strip()
    if context_value.upper() == "CG" or context_value.lower() == "cpg":
        return "CpG"
    context_upper = context_value.upper()
    if context_upper in {"CHG", "CHH"}:
        return context_upper
    return context_value


def _mianjifa_extract_chr_from_name(filename: str):
    """Extract chromosome information from names such as DMP_*_Chr1.txt."""
    m = re.search(r"chr([A-Za-z0-9]+)", filename, flags=re.IGNORECASE)
    if m:
        return _mianjifa_normalize_chr(m.group(1))
    return None


def _mianjifa_parse_output_pair(output_dir_name: str):
    """Parse output_wt1_mut2-like folders into group/replicate labels."""
    name = output_dir_name
    if not name.lower().startswith("output_"):
        raise ValueError(f"不是 output_ 开头的文件夹名：{name}")
    body = name[len("output_"):]
    m = re.match(r"^(.+?)(\d+)[_-](.+?)(\d+)$", body)
    if not m:
        raise ValueError(f"无法解析 output 文件夹名：{name}，期望格式类似 output_wt1_mut1")
    group1 = m.group(1)
    rep1 = int(m.group(2))
    group2 = m.group(3)
    rep2 = int(m.group(4))
    label = f"{group1}{rep1} - {group2}{rep2}"
    sort_key = (group1, rep1, group2, rep2)
    return {
        "group1": group1,
        "rep1": rep1,
        "group2": group2,
        "rep2": rep2,
        "label": label,
        "sort_key": sort_key,
    }


def _mianjifa_find_case_insensitive_dir(parent: Path, dirname: str) -> Path:
    exact = parent / dirname
    if exact.exists() and exact.is_dir():
        return exact
    for p in parent.iterdir():
        if p.is_dir() and p.name.lower() == dirname.lower():
            return p
    raise FileNotFoundError(f"{parent} 下找不到原始组文件夹：{dirname}")


def _mianjifa_resolve_group_dir(selected_dir: Path, group_name: str, mut_dir: str, wt_dir: str) -> Path:
    """Resolve output-folder group labels to actual input directories."""
    key = str(group_name).lower().strip("_- ")
    if key in {"wt", "wild", "wildtype", "control", "ctrl", "dcntrol", "dcontrol", "dcontrols"}:
        return Path(wt_dir)
    if key in {"mut", "mutation", "mutant", "case", "treat", "treated", "treatment", "dcases", "dcase"}:
        return Path(mut_dir)

    # Fallback: try a real folder under the run root, preserving original standalone behavior.
    return _mianjifa_find_case_insensitive_dir(selected_dir, group_name)


def _mianjifa_choose_grid(n: int):
    """Choose subplot layout for n comparisons."""
    if n <= 0:
        return 1, 1
    ncols = math.ceil(math.sqrt(n))
    nrows = math.ceil(n / ncols)
    return nrows, ncols


def _mianjifa_read_dmp_sites(dmp_file: Path) -> pd.DataFrame:
    """Read DMP site positions from a DMP* file; first column is position or chr:pos."""
    chr_from_name = _mianjifa_extract_chr_from_name(dmp_file.name)
    records = []
    with open(dmp_file, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split()
            if not parts:
                continue
            first = parts[0]
            if ":" in first:
                try:
                    chr_part, pos_part = first.split(":", 1)
                    chr_value = _mianjifa_normalize_chr(chr_part)
                    pos = int(float(pos_part))
                    records.append((chr_value, pos))
                    continue
                except Exception:
                    continue
            try:
                pos = int(float(first))
            except ValueError:
                continue
            records.append((chr_from_name, pos))
    if not records:
        return pd.DataFrame(columns=["chr", "pos"])
    df = pd.DataFrame(records, columns=["chr", "pos"])
    return df.drop_duplicates(subset=["chr", "pos"])


def _mianjifa_load_dmp_sites_from_methylation_type_dir(methylation_type_dir: Path) -> pd.DataFrame:
    dmp_files = sorted([p for p in methylation_type_dir.glob("DMP*") if p.is_file()])
    if not dmp_files:
        raise FileNotFoundError(f"{methylation_type_dir} 中没有找到 DMP* 文件")
    all_sites = []
    for dmp_file in dmp_files:
        sites = _mianjifa_read_dmp_sites(dmp_file)
        if not sites.empty:
            all_sites.append(sites)
    if not all_sites:
        return pd.DataFrame(columns=["chr", "pos"])
    result = pd.concat(all_sites, ignore_index=True)
    return result.drop_duplicates(subset=["chr", "pos"])



def _mianjifa_load_q_significant_sites_from_fdr(methylation_type_dir: Path):
    """Load q-significant sites from the complete pairwise FDR table.

    This deliberately ignores MethDiff so that mianjifa estimates its threshold
    from q-significant sites that have not already been truncated by a user or
    previous auto-methdiff threshold.
    """
    fdr_files = sorted(
        p for p in methylation_type_dir.glob("FDR_corrected_results_*.txt")
        if p.is_file()
    )
    if not fdr_files:
        raise FileNotFoundError(
            f"{methylation_type_dir} 中没有找到 FDR_corrected_results_*.txt"
        )
    if len(fdr_files) > 1:
        raise RuntimeError(
            f"{methylation_type_dir} 中找到多个FDR文件，无法唯一确定："
            f"{[p.name for p in fdr_files]}"
        )

    fdr_file = fdr_files[0]
    try:
        df = pd.read_csv(fdr_file, sep="\t")
        if len(df.columns) < 6:
            df = pd.read_csv(fdr_file, sep=r"\s+", engine="python")
    except Exception as exc:
        raise RuntimeError(f"读取FDR文件失败 {fdr_file}: {exc}") from exc

    required = {"Chromosome", "Position", "Qvalue"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(
            f"{fdr_file} 缺少必要列 {sorted(missing)}；实际列={list(df.columns)}"
        )

    qvalue = pd.to_numeric(df["Qvalue"], errors="coerce")
    fallback_q = get_dmp_threshold(methylation_type_dir.name)
    if "Qvalue_Threshold_Used" in df.columns:
        threshold = pd.to_numeric(
            df["Qvalue_Threshold_Used"], errors="coerce"
        ).fillna(float(fallback_q))
        threshold_source = "Qvalue_Threshold_Used"
    else:
        threshold = pd.Series(float(fallback_q), index=df.index)
        threshold_source = "fixed_context_threshold"

    sig_mask = qvalue.notna() & threshold.notna() & (qvalue <= threshold)
    sig_df = df.loc[sig_mask, ["Chromosome", "Position"]].copy()
    sig_df["chr"] = sig_df["Chromosome"].map(_mianjifa_normalize_chr)
    sig_df["pos"] = pd.to_numeric(sig_df["Position"], errors="coerce")
    sig_df = sig_df.dropna(subset=["pos"])
    sig_df["pos"] = sig_df["pos"].astype(np.int64)
    sites = sig_df[["chr", "pos"]].drop_duplicates().reset_index(drop=True)

    stats = {
        "candidate_source": "FDR_q_significant_without_methdiff_filter",
        "fdr_file": fdr_file.name,
        "fdr_total_rows": int(len(df)),
        "q_significant_site_count": int(len(sites)),
        "q_threshold_source": threshold_source,
        "q_threshold_min": float(threshold.min()) if len(threshold) else np.nan,
        "q_threshold_max": float(threshold.max()) if len(threshold) else np.nan,
    }
    return sites, stats


def _mianjifa_find_replicate_file(group_dir: Path, group_name: str, rep_id: int) -> Path:
    """Find raw methylation input file for a replicate in a group directory."""
    candidates = sorted([p for p in group_dir.glob(f"{rep_id}-*.txt") if p.is_file()])
    if not candidates:
        raise FileNotFoundError(f"{group_dir} 下找不到 {rep_id}-*.txt")

    exact_name = f"{rep_id}-{group_name}.txt".lower()
    for p in candidates:
        if p.name.lower() == exact_name:
            return p

    contains_group = [p for p in candidates if group_name.lower() in p.name.lower()]
    if contains_group:
        return sorted(contains_group, key=lambda x: len(x.name))[0]

    bad_keywords = [
        "bothmeunme", "diffchromo", "norepeated", "norepeat", "repeated", "repeat",
        "result", "dmp", "fet", "pvalue", "qvalue",
    ]
    filtered = []
    for p in candidates:
        lower = p.name.lower()
        if any(k in lower for k in bad_keywords):
            continue
        filtered.append(p)
    if filtered:
        return sorted(filtered, key=lambda x: len(x.name))[0]
    return sorted(candidates, key=lambda x: len(x.name))[0]


def _mianjifa_read_raw_methylation_file(file_path: Path, meth_type: str) -> pd.DataFrame:
    """Read raw 5-column methylation file and return chr,pos,context,meth_rate for one context."""
    normalized_meth_type = _mianjifa_normalize_context(meth_type)
    cache_key = (str(file_path.resolve()), str(normalized_meth_type).lower())
    if cache_key in AUTO_METHDIFF_RAW_CACHE:
        return AUTO_METHDIFF_RAW_CACHE[cache_key].copy()

    df = pd.read_csv(
        file_path,
        sep=r"\s+",
        header=None,
        names=["chr", "pos", "meth", "unmeth", "context"],
        usecols=[0, 1, 2, 3, 4],
        dtype={"chr": "string", "pos": "int64", "meth": "float64", "unmeth": "float64", "context": "string"},
        engine="python",
    )
    df["chr"] = df["chr"].map(_mianjifa_normalize_chr)
    df["context"] = df["context"].map(_mianjifa_normalize_context)
    df = df[df["context"].str.lower() == normalized_meth_type.lower()].copy()
    if df.empty:
        result = pd.DataFrame(columns=["chr", "pos", "context", "meth_rate"])
        AUTO_METHDIFF_RAW_CACHE[cache_key] = result.copy()
        return result

    df = df.groupby(["chr", "pos", "context"], as_index=False)[["meth", "unmeth"]].sum()
    total = df["meth"] + df["unmeth"]
    df = df[total > 0].copy()
    df["meth_rate"] = df["meth"] / (df["meth"] + df["unmeth"])
    result = df[["chr", "pos", "context", "meth_rate"]].copy()
    AUTO_METHDIFF_RAW_CACHE[cache_key] = result.copy()
    return result


def _mianjifa_filter_raw_by_dmp_sites(raw_df: pd.DataFrame, dmp_sites: pd.DataFrame) -> pd.DataFrame:
    if raw_df.empty or dmp_sites.empty:
        return raw_df.iloc[0:0].copy()
    parts = []
    with_chr = dmp_sites[dmp_sites["chr"].notna()].copy()
    without_chr = dmp_sites[dmp_sites["chr"].isna()].copy()
    if not with_chr.empty:
        with_chr["chr"] = with_chr["chr"].map(_mianjifa_normalize_chr)
        tmp = pd.merge(raw_df, with_chr[["chr", "pos"]].drop_duplicates(), on=["chr", "pos"], how="inner")
        parts.append(tmp)
    if not without_chr.empty:
        allowed_pos = set(without_chr["pos"].astype(np.int64).tolist())
        tmp = raw_df[raw_df["pos"].isin(allowed_pos)].copy()
        parts.append(tmp)
    if not parts:
        return raw_df.iloc[0:0].copy()
    out = pd.concat(parts, ignore_index=True)
    return out.drop_duplicates(subset=["chr", "pos", "context"])


def _mianjifa_calculate_methdiff_from_raw(
        selected_dir: Path,
        mut_dir: str,
        wt_dir: str,
        group1: str,
        rep1: int,
        group2: str,
        rep2: int,
        meth_type: str,
        dmp_sites: pd.DataFrame,
):
    """Calculate signed MethDiff = group1_rep1_meth_rate - group2_rep2_meth_rate."""
    group1_dir = _mianjifa_resolve_group_dir(selected_dir, group1, mut_dir=mut_dir, wt_dir=wt_dir)
    group2_dir = _mianjifa_resolve_group_dir(selected_dir, group2, mut_dir=mut_dir, wt_dir=wt_dir)

    group1_file = _mianjifa_find_replicate_file(group1_dir, group1, rep1)
    group2_file = _mianjifa_find_replicate_file(group2_dir, group2, rep2)

    df1 = _mianjifa_read_raw_methylation_file(group1_file, meth_type)
    df2 = _mianjifa_read_raw_methylation_file(group2_file, meth_type)
    df1 = _mianjifa_filter_raw_by_dmp_sites(df1, dmp_sites)
    df2 = _mianjifa_filter_raw_by_dmp_sites(df2, dmp_sites)

    merged = pd.merge(
        df1,
        df2,
        on=["chr", "pos", "context"],
        how="inner",
        suffixes=(f"_{group1}{rep1}", f"_{group2}{rep2}"),
    )

    stats = {
        "group1_file": group1_file.name,
        "group2_file": group2_file.name,
        "group1_dmp_rows": len(df1),
        "group2_dmp_rows": len(df2),
        "common_dmp_rows": len(merged),
    }
    if merged.empty:
        return np.array([], dtype=float), stats

    rate1_col = f"meth_rate_{group1}{rep1}"
    rate2_col = f"meth_rate_{group2}{rep2}"
    merged["MethDiff"] = merged[rate1_col] - merged[rate2_col]
    return merged["MethDiff"].to_numpy(dtype=float), stats


def _mianjifa_compute_area_threshold_from_hist(counts, bin_edges, side: str, cut_fraction: float):
    """Compute threshold by cutting a fraction of histogram area outward from zero."""
    widths = np.diff(bin_edges)
    if side == "left":
        idx = np.where(bin_edges[1:] <= 0)[0]
        idx = idx[::-1]
    elif side == "right":
        idx = np.where(bin_edges[:-1] >= 0)[0]
    else:
        raise ValueError("side must be 'left' or 'right'")
    if len(idx) == 0:
        return None, 0.0, 0.0
    total_area = float(np.sum(counts[idx] * widths[idx]))
    if total_area <= 0:
        return None, 0.0, 0.0
    target_area = total_area * float(cut_fraction)
    cumulative_area = 0.0
    for i in idx:
        left = bin_edges[i]
        right = bin_edges[i + 1]
        height = counts[i]
        width = widths[i]
        area = height * width
        if height <= 0 or area <= 0:
            continue
        if cumulative_area + area < target_area:
            cumulative_area += area
            continue
        remain_area = target_area - cumulative_area
        cut_width = remain_area / height
        if side == "left":
            threshold = right - cut_width
        else:
            threshold = left + cut_width
        threshold = max(left, min(right, threshold))
        return threshold, total_area, target_area
    if side == "left":
        return bin_edges[idx[-1]], total_area, target_area
    return bin_edges[idx[-1] + 1], total_area, target_area


def _mianjifa_draw_hist_with_area_cut(ax, counts, bin_edges, left_threshold, right_threshold):
    """Draw histogram and highlight the central area cut from zero."""
    widths = np.diff(bin_edges)
    for i, h in enumerate(counts):
        if h <= 0:
            continue
        left = bin_edges[i]
        right = bin_edges[i + 1]
        width = widths[i]
        if right <= 0 and left_threshold is not None:
            if right <= left_threshold:
                ax.bar(left, h, width=width, align="edge", color=AUTO_METHDIFF_COLOR_NORMAL, edgecolor="black", linewidth=0.3)
            elif left >= left_threshold:
                ax.bar(left, h, width=width, align="edge", color=AUTO_METHDIFF_COLOR_CUT, edgecolor="black", linewidth=0.3)
            else:
                ax.bar(left, h, width=left_threshold - left, align="edge", color=AUTO_METHDIFF_COLOR_NORMAL, edgecolor="black", linewidth=0.3)
                ax.bar(left_threshold, h, width=right - left_threshold, align="edge", color=AUTO_METHDIFF_COLOR_CUT, edgecolor="black", linewidth=0.3)
        elif left >= 0 and right_threshold is not None:
            if right <= right_threshold:
                ax.bar(left, h, width=width, align="edge", color=AUTO_METHDIFF_COLOR_CUT, edgecolor="black", linewidth=0.3)
            elif left >= right_threshold:
                ax.bar(left, h, width=width, align="edge", color=AUTO_METHDIFF_COLOR_NORMAL, edgecolor="black", linewidth=0.3)
            else:
                ax.bar(left, h, width=right_threshold - left, align="edge", color=AUTO_METHDIFF_COLOR_CUT, edgecolor="black", linewidth=0.3)
                ax.bar(right_threshold, h, width=right - right_threshold, align="edge", color=AUTO_METHDIFF_COLOR_NORMAL, edgecolor="black", linewidth=0.3)
        else:
            ax.bar(left, h, width=width, align="edge", color=AUTO_METHDIFF_COLOR_NORMAL, edgecolor="black", linewidth=0.3)


def _mianjifa_plot_one_comparison_on_ax(ax, methdiff_values, title, cut_fraction):
    """Plot one comparison's signed MethDiff histogram."""
    values = np.asarray(methdiff_values, dtype=float)
    values = values[np.isfinite(values)]
    values = values[(values >= AUTO_METHDIFF_X_MIN) & (values <= AUTO_METHDIFF_X_MAX)]
    if len(values) == 0:
        ax.text(0.5, 0.5, "No common DMP sites", ha="center", va="center", transform=ax.transAxes)
        ax.set_title(title, fontsize=10)
        ax.set_xlim(AUTO_METHDIFF_X_MIN, AUTO_METHDIFF_X_MAX)
        ax.set_xlabel("MethDiff")
        ax.set_ylabel("Frequency")
        return {
            "n": 0,
            "left_threshold": None,
            "right_threshold": None,
            "left_total_area": 0,
            "right_total_area": 0,
            "left_target_area": 0,
            "right_target_area": 0,
        }

    counts, bin_edges = np.histogram(values, bins=AUTO_METHDIFF_BINS)
    left_threshold, left_total_area, left_target_area = _mianjifa_compute_area_threshold_from_hist(
        counts, bin_edges, side="left", cut_fraction=cut_fraction
    )
    right_threshold, right_total_area, right_target_area = _mianjifa_compute_area_threshold_from_hist(
        counts, bin_edges, side="right", cut_fraction=cut_fraction
    )
    _mianjifa_draw_hist_with_area_cut(ax, counts, bin_edges, left_threshold, right_threshold)
    ax.axvline(0, color=AUTO_METHDIFF_COLOR_ZERO_LINE, linestyle="-", linewidth=1.2)
    if left_threshold is not None:
        ax.axvline(left_threshold, color=AUTO_METHDIFF_COLOR_THRESHOLD, linestyle="--", linewidth=1.5)
    if right_threshold is not None:
        ax.axvline(right_threshold, color=AUTO_METHDIFF_COLOR_THRESHOLD, linestyle="--", linewidth=1.5)

    threshold_text = []
    if left_threshold is not None:
        threshold_text.append(f"L={left_threshold:.4f}")
    if right_threshold is not None:
        threshold_text.append(f"R={right_threshold:.4f}")
    if threshold_text:
        title = title + "\n" + ", ".join(threshold_text)
    ax.set_title(title, fontsize=10)
    ax.set_xlabel("MethDiff")
    ax.set_ylabel("Frequency")
    ax.set_xlim(AUTO_METHDIFF_X_MIN, AUTO_METHDIFF_X_MAX)
    ax.grid(axis="y", alpha=0.3)
    return {
        "n": len(values),
        "left_threshold": left_threshold,
        "right_threshold": right_threshold,
        "left_total_area": left_total_area,
        "right_total_area": right_total_area,
        "left_target_area": left_target_area,
        "right_target_area": right_target_area,
    }


def _mianjifa_find_output_dirs(selected_dir: Path):
    """Find output_*_* folders that can be parsed as pairwise comparisons."""
    output_dirs = []
    for p in selected_dir.iterdir():
        if not p.is_dir() or not p.name.lower().startswith("output_"):
            continue
        try:
            info = _mianjifa_parse_output_pair(p.name)
        except Exception:
            continue
        output_dirs.append((info["sort_key"], info, p))
    return sorted(output_dirs, key=lambda x: x[0])


def _mianjifa_collect_methylation_types(output_dir_records):
    """Collect context folders containing complete pairwise FDR tables."""
    meth_types = set()
    for _, _, output_dir in output_dir_records:
        for sub in output_dir.iterdir():
            if not sub.is_dir():
                continue
            has_fdr = any(
                p.is_file()
                for p in sub.glob("FDR_corrected_results_*.txt")
            )
            if has_fdr:
                meth_types.add(sub.name)
    preferred = ["CpG", "CHG", "CHH"]
    return sorted(
        meth_types,
        key=lambda x: (preferred.index(x) if x in preferred else 99, x),
    )


def _mianjifa_plot_one_methylation_type(
        meth_type: str,
        output_dir_records,
        selected_dir: Path,
        mut_dir: str,
        wt_dir: str,
        output_dir: Path,
        cut_fraction: float):
    """Plot signed raw MethDiff for q-significant, MethDiff-unfiltered sites."""
    comparisons = []
    for sort_key, info, out_dir in output_dir_records:
        meth_dir = out_dir / meth_type
        if not meth_dir.exists() or not meth_dir.is_dir():
            comparisons.append({
                "sort_key": sort_key,
                "info": info,
                "output_dir": out_dir,
                "values": np.array([]),
                "stats": {},
                "error": f"missing {meth_type} folder",
            })
            continue
        try:
            qsig_sites, source_stats = _mianjifa_load_q_significant_sites_from_fdr(
                meth_dir
            )
            values, raw_stats = _mianjifa_calculate_methdiff_from_raw(
                selected_dir=selected_dir,
                mut_dir=mut_dir,
                wt_dir=wt_dir,
                group1=info["group1"],
                rep1=info["rep1"],
                group2=info["group2"],
                rep2=info["rep2"],
                meth_type=meth_type,
                dmp_sites=qsig_sites,
            )
            stats = {**source_stats, **raw_stats}
            stats["area_input_site_count"] = len(qsig_sites)
            comparisons.append({
                "sort_key": sort_key,
                "info": info,
                "output_dir": out_dir,
                "values": values,
                "stats": stats,
                "error": "",
            })
        except Exception as e:
            comparisons.append({
                "sort_key": sort_key,
                "info": info,
                "output_dir": out_dir,
                "values": np.array([]),
                "stats": {},
                "error": str(e),
            })

    comparisons = sorted(comparisons, key=lambda x: x["sort_key"])
    n = len(comparisons)
    nrows, ncols = _mianjifa_choose_grid(n)
    fig_width = max(6 * ncols, 10)
    fig_height = max(5 * nrows, 6)
    fig, axes = plt.subplots(nrows, ncols, figsize=(fig_width, fig_height))
    axes_flat = np.array(axes).reshape(-1)
    summary_records = []

    for ax, comp in zip(axes_flat, comparisons):
        info = comp["info"]
        title = f"{info['label']}\n{comp['output_dir'].name}"
        if comp["error"]:
            ax.text(
                0.5, 0.5, comp["error"], ha="center", va="center",
                transform=ax.transAxes, fontsize=8,
            )
            ax.set_title(title, fontsize=10)
            ax.set_xlim(AUTO_METHDIFF_X_MIN, AUTO_METHDIFF_X_MAX)
            ax.set_xlabel("MethDiff")
            ax.set_ylabel("Frequency")
            plot_stats = {
                "n": 0,
                "left_threshold": None,
                "right_threshold": None,
                "left_total_area": 0,
                "right_total_area": 0,
                "left_target_area": 0,
                "right_target_area": 0,
            }
        else:
            plot_stats = _mianjifa_plot_one_comparison_on_ax(
                ax, comp["values"], title, cut_fraction=cut_fraction
            )

        left_abs = np.nan
        right_abs = np.nan
        if plot_stats.get("left_threshold") is not None:
            left_abs = abs(float(plot_stats["left_threshold"]))
        if plot_stats.get("right_threshold") is not None:
            right_abs = abs(float(plot_stats["right_threshold"]))
        side_values = [
            v for v in [left_abs, right_abs] if np.isfinite(v)
        ]
        pair_abs_threshold = (
            float(np.mean(side_values)) if side_values else np.nan
        )

        record = {
            "methylation_type": meth_type,
            "output_folder": comp["output_dir"].name,
            "comparison": info["label"],
            "group1": info["group1"],
            "rep1": info["rep1"],
            "group2": info["group2"],
            "rep2": info["rep2"],
            "methdiff_definition": (
                f"{info['group1']}{info['rep1']} - "
                f"{info['group2']}{info['rep2']}"
            ),
            "error": comp["error"],
            "left_abs_threshold": left_abs,
            "right_abs_threshold": right_abs,
            "pair_abs_threshold_mean_lr": pair_abs_threshold,
            **comp.get("stats", {}),
            **plot_stats,
        }
        summary_records.append(record)

    for ax in axes_flat[len(comparisons):]:
        ax.axis("off")

    fig.suptitle(
        f"{selected_dir.name} | {meth_type} | {n} output comparisons | "
        f"q-significant sites without MethDiff prefilter | raw signed MethDiff | "
        f"cut {cut_fraction * 100:.2f}% histogram area from 0",
        fontsize=15,
    )
    legend_items = [
        Patch(
            facecolor=AUTO_METHDIFF_COLOR_NORMAL,
            edgecolor="black",
            label="Retained area",
        ),
        Patch(
            facecolor=AUTO_METHDIFF_COLOR_CUT,
            edgecolor="black",
            label=f"Cut {cut_fraction * 100:.2f}% area from 0",
        ),
        Line2D(
            [0], [0], color=AUTO_METHDIFF_COLOR_ZERO_LINE,
            lw=1.5, linestyle="-", label="x = 0",
        ),
        Line2D(
            [0], [0], color=AUTO_METHDIFF_COLOR_THRESHOLD,
            lw=1.5, linestyle="--", label="Cut threshold",
        ),
    ]
    fig.legend(
        handles=legend_items, loc="upper center", ncol=4,
        frameon=True, bbox_to_anchor=(0.5, 0.97),
    )
    plt.tight_layout(rect=[0, 0, 1, 0.92])
    output_png = output_dir / (
        f"{_mianjifa_safe_name(selected_dir.name)}_"
        f"{_mianjifa_safe_name(meth_type)}_"
        f"{n}outputs_Qsig_rawMethDiff_area{cut_fraction * 100:.2f}pct.png"
    )
    plt.savefig(output_png, dpi=300)
    plt.close()
    print(f"已保存 auto-methdiff 图片：{output_png}")
    return summary_records


def _mianjifa_aggregate_threshold(summary_df: pd.DataFrame, aggregate: str, fallback: float):
    """Aggregate left/right plot thresholds into one global abs(methdiff) threshold."""
    candidates = []
    for col in ["left_abs_threshold", "right_abs_threshold"]:
        if col in summary_df.columns:
            values = pd.to_numeric(summary_df[col], errors="coerce").replace([np.inf, -np.inf], np.nan).dropna()
            candidates.extend([float(v) for v in values if float(v) > 0])
    if not candidates:
        return float(fallback), 0, "fallback_no_valid_side_thresholds"
    arr = np.asarray(candidates, dtype=float)
    aggregate = str(aggregate).lower()
    if aggregate == "mean":
        threshold = float(np.mean(arr))
    elif aggregate == "max":
        threshold = float(np.max(arr))
    elif aggregate == "min":
        threshold = float(np.min(arr))
    else:
        threshold = float(np.median(arr))
        aggregate = "median"
    threshold = max(0.0, min(1.0, threshold))
    return threshold, int(len(arr)), f"ok_{aggregate}"


def estimate_mianjifa_auto_methdiff_threshold(
        mut_dir,
        wt_dir,
        work_dir=".",
        cut_fraction=0.05,
        fallback=0.3,
        aggregate="median",
        output_dir="and_output/auto_methdiff_thresholds",
        report_only=False):
    """Estimate one global abs threshold from q-significant sites before any MethDiff filtering."""
    selected_dir = Path(work_dir)
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    AUTO_METHDIFF_RAW_CACHE.clear()

    output_dir_records = _mianjifa_find_output_dirs(selected_dir)
    if not output_dir_records:
        summary_df = pd.DataFrame([{
            "status": "fallback_no_output_dirs",
            "auto_methdiff_threshold": np.nan,
            "used_methdiff_threshold": float(fallback),
            "fallback_methdiff": float(fallback),
            "cut_fraction": float(cut_fraction),
            "aggregate": aggregate,
            "report_only": bool(report_only),
        }])
        out_file = output_path / "mianjifa_auto_methdiff_threshold_summary.tsv"
        summary_df.to_csv(out_file, sep="\t", index=False)
        print(f"自动 methdiff 阈值估计失败：未找到 output_*_* 文件夹；回退到 {fallback}")
        print(f"自动 methdiff 诊断表保存至: {out_file}")
        return float(fallback), summary_df

    meth_types = _mianjifa_collect_methylation_types(output_dir_records)
    if not meth_types:
        summary_df = pd.DataFrame([{
            "status": "fallback_no_methylation_type_dirs_with_fdr",
            "auto_methdiff_threshold": np.nan,
            "used_methdiff_threshold": float(fallback),
            "fallback_methdiff": float(fallback),
            "cut_fraction": float(cut_fraction),
            "aggregate": aggregate,
            "report_only": bool(report_only),
        }])
        out_file = output_path / "mianjifa_auto_methdiff_threshold_summary.tsv"
        summary_df.to_csv(out_file, sep="\t", index=False)
        print(f"自动 methdiff 阈值估计失败：未找到包含完整FDR表的甲基化类型目录；回退到 {fallback}")
        print(f"自动 methdiff 诊断表保存至: {out_file}")
        return float(fallback), summary_df

    print("\n估计 auto-methdiff 阈值：mianjifa q-significant raw MethDiff area-cut 方法")
    print(f"  work_dir = {selected_dir}")
    print(f"  mut_dir = {mut_dir}")
    print(f"  wt_dir = {wt_dir}")
    print(f"  cut_fraction = {cut_fraction}")
    print(f"  aggregate = {aggregate}")

    all_summary = []
    for meth_type in meth_types:
        print(f"\n正在处理 auto-methdiff 甲基化类型：{meth_type}")
        records = _mianjifa_plot_one_methylation_type(
            meth_type=meth_type,
            output_dir_records=output_dir_records,
            selected_dir=selected_dir,
            mut_dir=mut_dir,
            wt_dir=wt_dir,
            output_dir=output_path,
            cut_fraction=float(cut_fraction),
        )
        all_summary.extend(records)

    summary_df = pd.DataFrame(all_summary)
    threshold, n_threshold_values, status = _mianjifa_aggregate_threshold(
        summary_df=summary_df,
        aggregate=aggregate,
        fallback=float(fallback),
    )

    summary_df["auto_methdiff_global_status"] = status
    summary_df["auto_methdiff_global_threshold"] = threshold if status.startswith("ok") else np.nan
    summary_df["used_methdiff_threshold"] = threshold
    summary_df["fallback_methdiff"] = float(fallback)
    summary_df["cut_fraction"] = float(cut_fraction)
    summary_df["aggregate"] = aggregate
    summary_df["n_threshold_values_for_aggregate"] = n_threshold_values
    summary_df["report_only"] = bool(report_only)

    summary_tsv = output_path / "mianjifa_auto_methdiff_threshold_summary.tsv"
    summary_csv = output_path / "mianjifa_auto_methdiff_threshold_summary.csv"
    summary_df.to_csv(summary_tsv, sep="\t", index=False)
    summary_df.to_csv(summary_csv, index=False, encoding="utf-8-sig")

    aggregate_df = pd.DataFrame([{
        "method": "mianjifa_Qsig_rawMethDiff_area_cut",
        "status": status,
        "auto_methdiff_threshold": threshold if status.startswith("ok") else np.nan,
        "used_methdiff_threshold": threshold,
        "fallback_methdiff": float(fallback),
        "cut_fraction": float(cut_fraction),
        "aggregate": aggregate,
        "n_threshold_values_for_aggregate": n_threshold_values,
        "report_only": bool(report_only),
    }])
    aggregate_tsv = output_path / "mianjifa_auto_methdiff_threshold_aggregate.tsv"
    aggregate_df.to_csv(aggregate_tsv, sep="\t", index=False)

    print(f"自动 methdiff 诊断表保存至: {summary_tsv}")
    print(f"自动 methdiff 聚合阈值表保存至: {aggregate_tsv}")
    print(f"mianjifa auto-methdiff global abs threshold = {threshold:.6g} ({status})")
    return threshold, summary_df



def estimate_auto_qvalue_threshold_twostep(
        pvalues,
        qvalues,
        pvalue_cutoff=0.05,
        fallback_q=0.05,
        min_candidates=10,
        use_smooth=False,
        smooth_sigma=4.0):
    """Estimate an adaptive q-value threshold for two-step FDR results.

    This function is intended ONLY for contexts where the workflow has already:
      1) kept sites with raw p-value <= pvalue_cutoff; and
      2) recalculated q-values within that pre-filtered subset.

    Rule:
      - clean finite p/q pairs;
      - keep p <= pvalue_cutoff;
      - sort by increasing p-value;
      - compute diff = qvalue - pvalue;
      - choose the row with maximum diff;
      - return its qvalue as the adaptive threshold.

    If the estimation is not possible, return fallback_q and a diagnostic dict.
    The function does not apply an upper q-value cap because, in the intended
    two-step setting, q-values should naturally remain near the pre-filter range.
    """
    fallback_q = float(fallback_q)
    p = np.asarray(pvalues, dtype=float)
    q = np.asarray(qvalues, dtype=float)

    valid = np.isfinite(p) & np.isfinite(q)
    p = p[valid]
    q = q[valid]

    if len(p) == 0:
        return fallback_q, {
            "auto_q_status": "fallback_no_valid_pq",
            "auto_q_threshold": fallback_q,
            "auto_q_raw_threshold": np.nan,
            "auto_q_pvalue_at_max": np.nan,
            "auto_q_diff_at_max": np.nan,
            "auto_q_rank_at_max": 0,
            "auto_q_n_candidates": 0,
            "auto_q_pvalue_cutoff": float(pvalue_cutoff),
            "auto_q_fallback_q": fallback_q,
            "auto_q_used_smooth": bool(use_smooth),
        }

    candidate_mask = p <= float(pvalue_cutoff)
    p_sub = p[candidate_mask]
    q_sub = q[candidate_mask]

    if len(p_sub) < int(min_candidates):
        return fallback_q, {
            "auto_q_status": "fallback_too_few_candidates",
            "auto_q_threshold": fallback_q,
            "auto_q_raw_threshold": np.nan,
            "auto_q_pvalue_at_max": np.nan,
            "auto_q_diff_at_max": np.nan,
            "auto_q_rank_at_max": 0,
            "auto_q_n_candidates": int(len(p_sub)),
            "auto_q_pvalue_cutoff": float(pvalue_cutoff),
            "auto_q_fallback_q": fallback_q,
            "auto_q_used_smooth": bool(use_smooth),
        }

    order = np.argsort(p_sub, kind="mergesort")
    p_sub = p_sub[order]
    q_sub = q_sub[order]

    diff = q_sub - p_sub
    diff_for_max = diff

    if use_smooth:
        try:
            from scipy.ndimage import gaussian_filter1d
            diff_for_max = gaussian_filter1d(diff, sigma=float(smooth_sigma))
        except Exception:
            diff_for_max = diff
            use_smooth = False

    if not np.any(np.isfinite(diff_for_max)):
        return fallback_q, {
            "auto_q_status": "fallback_no_finite_diff",
            "auto_q_threshold": fallback_q,
            "auto_q_raw_threshold": np.nan,
            "auto_q_pvalue_at_max": np.nan,
            "auto_q_diff_at_max": np.nan,
            "auto_q_rank_at_max": 0,
            "auto_q_n_candidates": int(len(p_sub)),
            "auto_q_pvalue_cutoff": float(pvalue_cutoff),
            "auto_q_fallback_q": fallback_q,
            "auto_q_used_smooth": bool(use_smooth),
        }

    max_idx = int(np.nanargmax(diff_for_max))
    best_q = float(q_sub[max_idx])

    info = {
        "auto_q_status": "ok",
        "auto_q_threshold": best_q,
        "auto_q_raw_threshold": best_q,
        "auto_q_pvalue_at_max": float(p_sub[max_idx]),
        "auto_q_diff_at_max": float(diff[max_idx]),
        "auto_q_rank_at_max": int(max_idx + 1),
        "auto_q_n_candidates": int(len(p_sub)),
        "auto_q_pvalue_cutoff": float(pvalue_cutoff),
        "auto_q_fallback_q": fallback_q,
        "auto_q_used_smooth": bool(use_smooth),
    }

    return best_q, info


def save_auto_qvalue_report(output_dir, replicate_x, replicate_y, methylation_type, info, threshold_used, report_only=False):
    """Save the adaptive q-value diagnostic table for one pairwise comparison/context."""
    if info is None:
        return
    os.makedirs(output_dir, exist_ok=True)
    record = dict(info)
    record.update({
        "methylation_type": methylation_type,
        "replicate_x_mut": int(replicate_x),
        "replicate_y_wt": int(replicate_y),
        "threshold_used_for_dmp_calling": float(threshold_used),
        "report_only": bool(report_only),
    })
    out_file = os.path.join(output_dir, f"auto_qvalue_threshold_{methylation_type}_wt_replicate{replicate_y}_vs_mut_replicate{replicate_x}.tsv")
    pd.DataFrame([record]).to_csv(out_file, sep="\t", index=False)

    print(f"    自动q-value阈值诊断表保存至: {out_file}")

# =========================================================
# Auto q-value plotting module
# adapted from auto_qvalue_and_plot.py
# =========================================================

AUTO_Q_PLOT_FILE_EXTENSIONS = (".txt", ".tsv", ".csv")
AUTO_Q_PLOT_LINE_WIDTH = 1.5
AUTO_Q_PLOT_DASH_PATTERN = (3, 1.5)
AUTO_Q_PLOT_FONT_FAMILY = "Arial"
AUTO_Q_PLOT_LABEL_SIZE = 10
AUTO_Q_PLOT_TITLE_SIZE = 12
AUTO_Q_PLOT_TICK_SIZE = 9

AUTO_Q_PLOT_DEFAULT_COLOR_PAIRS = {
    "wt2_mut1": ("#ffc400", "#ffaa00"),
    "wt2_mut2": ("#a6cee3", "#1f78b4"),
    "wt1_mut2": ("#7bc96f", "#33a02c"),
    "wt1_mut1": ("#ff9999", "#e31a1c"),

    "wt1_mut3": ("#b2df8a", "#006d2c"),
    "wt2_mut3": ("#cab2d6", "#6a3d9a"),
    "wt3_mut1": ("#fdbf6f", "#ff7f00"),
    "wt3_mut2": ("#fb9a99", "#e31a1c"),
    "wt3_mut3": ("#bdbdbd", "#525252"),
}

AUTO_Q_PLOT_PREFERRED_ORDER_KEYS = [
    "wt2_mut1", "wt2_mut2", "wt1_mut2", "wt1_mut1",
    "wt1_mut1", "wt1_mut2", "wt1_mut3",
    "wt2_mut1", "wt2_mut2", "wt2_mut3",
    "wt3_mut1", "wt3_mut2", "wt3_mut3",
]


def _autoq_natural_key(s):
    return [
        int(x) if x.isdigit() else x.lower()
        for x in re.split(r"(\d+)", str(s))
    ]


def _autoq_normalize_key(text):
    key = str(text).lower()
    key = re.sub(r"(?i)_?replicate", "", key)
    key = re.sub(r"(?i)_?rep", "", key)
    key = key.replace("_vs_", "_")
    key = re.sub(r"[^a-z0-9]+", "_", key)
    key = re.sub(r"_+", "_", key).strip("_")
    return key


def _autoq_make_display_name(comparison_name):
    name = str(comparison_name)
    name = re.sub(r"(?i)_?replicate", "_rep", name)
    name = re.sub(r"(?i)(?:[_-]?cpg)$", "", name)
    name = re.sub(r"_+", "_", name).strip("_")
    return name


def _autoq_blend_with_white(color, alpha=0.55):
    rgb = np.array(mcolors.to_rgb(color))
    white = np.array([1.0, 1.0, 1.0])
    blended = rgb * (1 - alpha) + white * alpha
    return mcolors.to_hex(blended)


def _autoq_build_color_pairs(ordered_keys):
    color_pairs = {}

    for key in ordered_keys:
        if key in AUTO_Q_PLOT_DEFAULT_COLOR_PAIRS:
            color_pairs[key] = AUTO_Q_PLOT_DEFAULT_COLOR_PAIRS[key]

    cmap = plt.get_cmap("tab20")
    auto_i = 0

    for key in ordered_keys:
        if key in color_pairs:
            continue

        dark = mcolors.to_hex(cmap(auto_i % cmap.N))
        light = _autoq_blend_with_white(dark, alpha=0.55)
        color_pairs[key] = (light, dark)
        auto_i += 1

    return color_pairs


def _autoq_infer_subplot_layout(n):
    if n <= 0:
        return 1, 1

    ncols = int(np.ceil(np.sqrt(n)))
    nrows = int(np.ceil(n / ncols))

    if n in (5, 6):
        nrows, ncols = 2, 3

    return nrows, ncols


def _autoq_setup_sci_style():
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": [AUTO_Q_PLOT_FONT_FAMILY, "Helvetica", "DejaVu Sans"],
        "font.size": AUTO_Q_PLOT_LABEL_SIZE,
        "axes.linewidth": 0.8,
        "axes.edgecolor": "black",
        "xtick.direction": "out",
        "ytick.direction": "out",
        "xtick.labelsize": AUTO_Q_PLOT_TICK_SIZE,
        "ytick.labelsize": AUTO_Q_PLOT_TICK_SIZE,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "legend.frameon": False,
        "mathtext.default": "regular",
    })


def _autoq_format_rank_axis(ax):
    ax.xaxis.set_major_formatter(ticker.ScalarFormatter(useMathText=True))
    ax.ticklabel_format(style="sci", axis="x", scilimits=(0, 0))


def _autoq_smooth(y):
    y = np.asarray(y, dtype=float)
    if len(y) < 3:
        return y
    return gaussian_filter1d(y, sigma=float(AUTO_QVALUE_SMOOTH_SIGMA))


def _autoq_read_pq_file(filepath):
    try:
        table = pd.read_csv(filepath, sep=None, engine="python")
        lower_cols = {str(c).strip().lower(): c for c in table.columns}

        p_candidates = [
            "pvalue", "p_value", "p.value", "p-value", "p", "pval", "p_val"
        ]
        q_candidates = [
            "qvalue", "q_value", "q.value", "q-value",
            "adj_p", "adj.p", "adjusted_p", "adjusted.p",
            "padj", "p_adj", "fdr"
        ]

        p_col = next((lower_cols[c] for c in p_candidates if c in lower_cols), None)
        q_col = next((lower_cols[c] for c in q_candidates if c in lower_cols), None)

        if p_col is not None and q_col is not None:
            df = pd.DataFrame({
                "p": pd.to_numeric(table[p_col], errors="coerce"),
                "adj_p": pd.to_numeric(table[q_col], errors="coerce"),
            })
        else:
            df = table.iloc[:, :2].copy()
            df.columns = ["p", "adj_p"]
            df["p"] = pd.to_numeric(df["p"], errors="coerce")
            df["adj_p"] = pd.to_numeric(df["adj_p"], errors="coerce")

    except Exception as e:
        print(f"  [auto-q plot] 读取失败 {filepath}: {e}")
        return None, None, None

    df = df.dropna(subset=["p", "adj_p"])
    if df.empty:
        return None, None, None

    df = df.sort_values("p").reset_index(drop=True)
    df.index = np.arange(1, len(df) + 1)

    return (
        df.index.values.astype(float),
        df["p"].values.astype(float),
        df["adj_p"].values.astype(float),
    )


def _autoq_keep_pvalue_below_cutoff(p, adj, cutoff):
    p = np.asarray(p, dtype=float)
    adj = np.asarray(adj, dtype=float)

    mask = p < float(cutoff)
    if not np.any(mask):
        return None, None, None

    p_sub = p[mask]
    adj_sub = adj[mask]
    rank_sub = np.arange(1, len(p_sub) + 1, dtype=float)

    return rank_sub, p_sub, adj_sub


def _autoq_compute_max_diff_point(rank, p, adj):
    rank = np.asarray(rank, dtype=float)
    p = np.asarray(p, dtype=float)
    adj = np.asarray(adj, dtype=float)

    p_s = _autoq_smooth(p)
    adj_s = _autoq_smooth(adj)

    raw_diff = adj - p
    smooth_diff = _autoq_smooth(raw_diff)

    diff_for_max = smooth_diff if AUTO_QVALUE_USE_SMOOTH else raw_diff

    if not np.any(np.isfinite(diff_for_max)):
        return None

    max_idx = int(np.nanargmax(diff_for_max))

    return {
        "rank": rank,
        "p": p,
        "adj": adj,
        "p_s": p_s,
        "adj_s": adj_s,
        "raw_diff": raw_diff,
        "smooth_diff": smooth_diff,
        "max_idx": max_idx,
        "max_rank": rank[max_idx],
        "max_p": p[max_idx],
        "max_adj_p": adj[max_idx],
        "max_diff_raw": raw_diff[max_idx],
        "max_diff_smooth": smooth_diff[max_idx],
        "max_diff_used": diff_for_max[max_idx],
        "max_method": "smooth_diff" if AUTO_QVALUE_USE_SMOOTH else "raw_diff",
    }


def _autoq_collect_context_fdr_files(work_dir, m, n, methylation_type):
    input_files = {}
    display_name = {}

    for mut_idx in range(1, int(m) + 1):
        for wt_idx in range(1, int(n) + 1):
            fpath = os.path.join(
                work_dir,
                f"output_wt{wt_idx}_mut{mut_idx}",
                methylation_type,
                f"FDR_corrected_results_wt_replicate{wt_idx}_vs_mut_replicate{mut_idx}.txt"
            )

            if not os.path.isfile(fpath):
                continue

            comparison = f"wt_replicate{wt_idx}_vs_mut_replicate{mut_idx}"
            key = _autoq_normalize_key(comparison)
            label = _autoq_make_display_name(comparison)

            input_files[key] = fpath
            display_name[key] = label

    def sort_key(k):
        if k in AUTO_Q_PLOT_PREFERRED_ORDER_KEYS:
            return (0, AUTO_Q_PLOT_PREFERRED_ORDER_KEYS.index(k), _autoq_natural_key(display_name.get(k, k)))
        return (1, 9999, _autoq_natural_key(display_name.get(k, k)))

    ordered_keys = sorted(input_files.keys(), key=sort_key)
    return input_files, ordered_keys, display_name


def _autoq_save_max_diff_table(diff_cache, ordered_keys, display_name, save_dir):
    rows = []

    for key in ordered_keys:
        if key not in diff_cache:
            continue

        info = diff_cache[key]
        rows.append({
            "sample": display_name[key],
            "max_method": info["max_method"],
            "max_rank": info["max_rank"],
            "p_value_at_max": info["max_p"],
            "adjusted_p_value_at_max": info["max_adj_p"],
            "max_diff_raw_adj_minus_p": info["max_diff_raw"],
            "smooth_diff_at_max": info["max_diff_smooth"],
        })

    if not rows:
        print("  [auto-q plot] 没有 max-diff 结果可保存")
        return None

    out_df = pd.DataFrame(rows)
    out_csv = os.path.join(save_dir, f"max_pvalue_adjustedP_diff_points_{len(rows)}comparisons.csv")
    out_df.to_csv(out_csv, index=False)

    print(f"  [auto-q plot] max-diff 表保存至: {out_csv}")
    return out_csv


def _autoq_plot_fdr_panels(diff_cache, ordered_keys, display_name, color_pairs, save_dir, p_cutoff):
    n = len(ordered_keys)
    nrows, ncols = _autoq_infer_subplot_layout(n)

    fig_width = max(8.5, 3.8 * ncols)
    fig_height = max(6.5, 3.3 * nrows)

    fig, axes = plt.subplots(nrows, ncols, figsize=(fig_width, fig_height), sharey=True)
    axes = np.atleast_1d(axes).ravel()

    for i, key in enumerate(ordered_keys):
        ax = axes[i]

        if key not in diff_cache:
            ax.axis("off")
            continue

        c_light, c_dark = color_pairs[key]
        info = diff_cache[key]

        ax.plot(info["rank"], info["p_s"], color=c_light, lw=AUTO_Q_PLOT_LINE_WIDTH, zorder=10)
        ax.plot(
            info["rank"], info["adj_s"],
            color=c_dark,
            lw=AUTO_Q_PLOT_LINE_WIDTH,
            linestyle="-",
            dashes=AUTO_Q_PLOT_DASH_PATTERN,
            dash_capstyle="butt",
            zorder=11,
        )

        y_mark = info["adj_s"][info["max_idx"]]
        ax.plot(
            info["max_rank"], y_mark,
            marker="o",
            mfc="white",
            mec=c_dark,
            mew=1.4,
            markersize=6,
            zorder=20,
        )
        ax.axvline(info["max_rank"], color=c_dark, lw=0.8, ls=":", alpha=0.45, zorder=5)

        ax.set_title(display_name[key], fontsize=AUTO_Q_PLOT_TITLE_SIZE, fontweight="bold", pad=4)
        ax.set_ylim(0, float(p_cutoff))
        ax.set_yticks(np.linspace(0, float(p_cutoff), 6))
        ax.grid(True, axis="y", ls=":", lw=0.6, color="#bbbbbb", alpha=0.5)
        _autoq_format_rank_axis(ax)
        ax.axhline(0.049, color="black", ls=":", lw=1.2, zorder=5)

        row, _ = divmod(i, ncols)
        if row == nrows - 1:
            ax.set_xlabel("sites ranked by increasing p-value", linespacing=1.0)
        else:
            ax.set_xticklabels([])

    for j in range(n, len(axes)):
        axes[j].axis("off")

    fig.supylabel("p-value / adjusted p-value", fontsize=AUTO_Q_PLOT_LABEL_SIZE)

    legend_handles = [
        Line2D([0], [0], color="black", lw=1.5, label="p-value"),
        Line2D([0], [0], color="black", lw=1.5, linestyle="--",
               dashes=AUTO_Q_PLOT_DASH_PATTERN, label="adjusted p-value"),
        Line2D([0], [0], marker="o", color="w", mfc="white", mec="black",
               mew=1.4, markersize=6, label="max diff point"),
    ]

    fig.legend(
        handles=legend_handles,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.99),
        ncol=3,
        fontsize=AUTO_Q_PLOT_LABEL_SIZE,
        frameon=False,
        labelspacing=0.6,
        handletextpad=0.5,
    )

    plt.subplots_adjust(
        left=0.08,
        bottom=0.10,
        right=0.98,
        top=0.90,
        wspace=0.18,
        hspace=0.32,
    )

    save_path = os.path.join(
        save_dir,
        f"FDR_{n}_comparisons_p_lt_{p_cutoff}_panels_with_max_diff.png"
    )
    fig.savefig(save_path, dpi=600, bbox_inches="tight")
    fig.savefig(save_path.replace(".png", ".pdf"), bbox_inches="tight")
    plt.close(fig)

    print(f"  [auto-q plot] FDR panel 图保存至: {save_path}")
    return save_path


def _autoq_plot_diff_panels(diff_cache, ordered_keys, display_name, color_pairs, save_dir):
    n = len(ordered_keys)
    nrows, ncols = _autoq_infer_subplot_layout(n)

    fig_width = max(8.5, 3.8 * ncols)
    fig_height = max(6.5, 3.3 * nrows)

    fig, axes = plt.subplots(nrows, ncols, figsize=(fig_width, fig_height), sharey=False)
    axes = np.atleast_1d(axes).ravel()

    for i, key in enumerate(ordered_keys):
        ax = axes[i]

        if key not in diff_cache:
            ax.axis("off")
            continue

        info = diff_cache[key]
        _, c_dark = color_pairs[key]

        ax.plot(info["rank"], info["smooth_diff"], color=c_dark, lw=AUTO_Q_PLOT_LINE_WIDTH, zorder=10)

        y_mark = info["smooth_diff"][info["max_idx"]]
        ax.plot(
            info["max_rank"], y_mark,
            marker="o",
            mfc="white",
            mec=c_dark,
            mew=1.4,
            markersize=6.5,
            zorder=20,
        )
        ax.axvline(info["max_rank"], color=c_dark, lw=0.8, ls=":", alpha=0.45, zorder=5)

        ax.set_title(display_name[key], fontsize=AUTO_Q_PLOT_TITLE_SIZE, fontweight="bold")
        ax.set_xlabel("sites ranked by increasing p-value", fontsize=AUTO_Q_PLOT_LABEL_SIZE)
        ax.set_ylabel("adjusted p-value - p-value", fontsize=AUTO_Q_PLOT_LABEL_SIZE)
        ax.grid(True, axis="y", ls=":", lw=0.6, color="#bbbbbb", alpha=0.5)
        _autoq_format_rank_axis(ax)

    for j in range(n, len(axes)):
        axes[j].axis("off")

    legend_handles = [
        Line2D(
            [0], [0],
            marker="o",
            color="w",
            mfc="white",
            mec="black",
            mew=1.4,
            markersize=6.5,
            label="max diff point",
        )
    ]

    fig.legend(
        handles=legend_handles,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.98),
        fontsize=AUTO_Q_PLOT_LABEL_SIZE,
        frameon=False,
        labelspacing=0.6,
        handletextpad=0.5,
    )

    fig.suptitle(
        "Difference between adjusted p-value and p-value",
        fontsize=AUTO_Q_PLOT_TITLE_SIZE + 1,
        fontweight="bold",
        y=1.02,
    )

    fig.tight_layout(rect=[0, 0, 1, 0.95])

    out_png = os.path.join(save_dir, f"Pvalue_AdjustedP_Difference_{n}panels.png")
    fig.savefig(out_png, dpi=600, bbox_inches="tight")
    fig.savefig(out_png.replace(".png", ".pdf"), bbox_inches="tight")
    plt.close(fig)

    print(f"  [auto-q plot] diff panel 图保存至: {out_png}")
    return out_png


def plot_auto_qvalue_panels_for_context(work_dir, m, n, methylation_type, unfilter_mtypes):
    """
    对一个 methylation context 汇总所有 pair 的 FDR_corrected_results_*.txt，
    生成与 auto_qvalue_and_plot.py 视觉一致的两张图：
      1. p-value / adjusted p-value 曲线 + max diff point
      2. adjusted p-value - p-value 差值曲线 + max diff point
    """
    if not (AUTO_QVALUE_TWOSTEP or AUTO_QVALUE_REPORT_ONLY):
        return

    # auto-qvalue 只用于两步法 context；非两步法 context 不画，避免误解。
    if methylation_type in unfilter_mtypes:
        print(f"  [auto-q plot] {methylation_type} 不是两步法 context，跳过 auto-qvalue 绘图")
        return

    input_files, ordered_keys, display_name = _autoq_collect_context_fdr_files(
        work_dir=work_dir,
        m=m,
        n=n,
        methylation_type=methylation_type
    )

    if not ordered_keys:
        print(f"  [auto-q plot] 未找到 {methylation_type} 的 FDR_corrected_results 文件，跳过")
        return

    save_dir = os.path.join(work_dir, "and_output", "auto_qvalue_plots", methylation_type)
    os.makedirs(save_dir, exist_ok=True)

    _autoq_setup_sci_style()
    color_pairs = _autoq_build_color_pairs(ordered_keys)
    diff_cache = {}

    for key in ordered_keys:
        rank, p, adj = _autoq_read_pq_file(input_files[key])
        if rank is None:
            print(f"  [auto-q plot] 读取失败，跳过: {input_files[key]}")
            continue

        rank_sub, p_sub, adj_sub = _autoq_keep_pvalue_below_cutoff(
            p,
            adj,
            cutoff=AUTO_QVALUE_P_CUTOFF
        )
        if rank_sub is None:
            print(f"  [auto-q plot] {display_name[key]} 无 p < {AUTO_QVALUE_P_CUTOFF} 位点，跳过")
            continue

        info = _autoq_compute_max_diff_point(rank_sub, p_sub, adj_sub)
        if info is None:
            print(f"  [auto-q plot] {display_name[key]} max-diff 计算失败，跳过")
            continue

        diff_cache[key] = info

    if not diff_cache:
        print(f"  [auto-q plot] {methylation_type} 没有可绘图数据")
        return

    print(f"  [auto-q plot] {methylation_type}: 检测到 {len(diff_cache)} 个 pairwise FDR 文件")
    _autoq_save_max_diff_table(diff_cache, ordered_keys, display_name, save_dir)
    _autoq_plot_fdr_panels(diff_cache, ordered_keys, display_name, color_pairs, save_dir, AUTO_QVALUE_P_CUTOFF)
    _autoq_plot_diff_panels(diff_cache, ordered_keys, display_name, color_pairs, save_dir)


def plot_all_auto_qvalue_panels(work_dir, m, n, unfilter_mtypes):
    if not (AUTO_QVALUE_TWOSTEP or AUTO_QVALUE_REPORT_ONLY):
        return

    print("\n生成 auto q-value 阈值诊断图...")
    for mtype in ['CpG', 'CHH', 'CHG']:
        plot_auto_qvalue_panels_for_context(
            work_dir=work_dir,
            m=m,
            n=n,
            methylation_type=mtype,
            unfilter_mtypes=unfilter_mtypes
        )


def _auto_vote_fallback_required_count(total_groups):
    """Fallback integer threshold corresponding to the configured VOTE_THRESHOLD."""
    return int(np.floor(float(VOTE_THRESHOLD) * int(total_groups) + 0.5))


def _auto_vote_group_columns(wide_df):
    return [c for c in wide_df.columns if isinstance(c, int) or (isinstance(c, str) and str(c).isdigit())]


def _estimate_vote_required_count_from_counts(counts, total_groups, allow_truncated_adjust=False):
    """Estimate an integer vote threshold from a support-count distribution.

    This is the embedded version of the standalone vote-threshold logic. It keeps
    the same main idea: fit a two-component GMM on support counts and choose the
    valley between the two components; if the valley cannot be found, use the
    Otsu-style between-class variance fallback. For DMR, the historical module
    additionally adjusts truncated distributions by lowering the threshold by 1.
    """
    total_groups = int(total_groups)
    fallback = _auto_vote_fallback_required_count(total_groups)

    counts = counts.reindex(range(1, total_groups + 1), fill_value=0)
    nonzero_counts = counts[counts > 0]
    if len(nonzero_counts) == 0:
        return None, {
            "status": "fallback_no_candidates",
            "method": "fallback",
            "auto_required_count": np.nan,
            "fallback_required_count": fallback,
            "valley_x": np.nan,
            "mean_low": np.nan,
            "mean_high": np.nan,
            "n_candidates": 0,
        }

    # 单峰时 GMM 没有可靠谷底，使用 2/3/当前 vote-threshold 对应的兜底整数阈值。
    if len(nonzero_counts) == 1:
        single_value = int(nonzero_counts.index[0])
        return fallback, {
            "status": "fallback_single_support_count",
            "method": "fallback_vote_threshold",
            "auto_required_count": fallback,
            "fallback_required_count": fallback,
            "valley_x": np.nan,
            "mean_low": float(single_value),
            "mean_high": float(single_value),
            "n_candidates": int(counts.sum()),
        }

    x = np.repeat(np.arange(1, total_groups + 1), counts.values)
    if len(x) == 0:
        return None, {
            "status": "fallback_no_candidates",
            "method": "fallback",
            "auto_required_count": np.nan,
            "fallback_required_count": fallback,
            "valley_x": np.nan,
            "mean_low": np.nan,
            "mean_high": np.nan,
            "n_candidates": 0,
        }

    gmm = GaussianMixture(n_components=2, random_state=42, max_iter=500)
    gmm.fit(x.reshape(-1, 1))
    means = np.sort(gmm.means_.flatten())

    x_plot = np.linspace(0.5, total_groups + 0.5, 500)
    pdf = np.exp(gmm.score_samples(x_plot.reshape(-1, 1)))

    search = (x_plot > means[0]) & (x_plot < means[1])
    if np.any(search):
        valley_idx = int(np.argmin(pdf[search]))
        valley_x = float(x_plot[search][valley_idx])
        best_t = int(np.ceil(valley_x))
        method = "two_component_gmm_valley"
    else:
        total = counts.sum()
        indices = np.arange(1, total_groups + 1)
        sum_val = np.sum(indices * counts)
        sumc = np.cumsum(counts)
        weight0 = sumc / total
        weight1 = 1 - weight0
        mean0 = np.cumsum(indices * counts) / (sumc + 1e-10)
        mean1 = (sum_val - np.cumsum(indices * counts)) / (total - sumc + 1e-10)
        var_between = weight0 * weight1 * (mean0 - mean1) ** 2
        best_t = int(np.argmax(var_between) + 1)
        valley_x = float(best_t - 0.5)
        method = "otsu_fallback"

    status = "ok"
    if allow_truncated_adjust:
        half_limit = total_groups // 2
        low_range_zero = all(counts.loc[i] == 0 for i in range(1, half_limit + 1))
        if low_range_zero:
            original_t = best_t
            best_t = max(1, best_t - 1)
            status = "ok_truncated_adjusted"
            method = f"{method}_truncated_adjust_{original_t}_to_{best_t}"

    best_t = max(1, min(int(best_t), total_groups))
    return best_t, {
        "status": status,
        "method": method,
        "auto_required_count": int(best_t),
        "fallback_required_count": fallback,
        "valley_x": valley_x,
        "mean_low": float(means[0]),
        "mean_high": float(means[1]),
        "n_candidates": int(counts.sum()),
    }


def _plot_vote_support_distribution(counts, total_groups, required_count, info, out_file, title, ylabel):
    """Save a support-count distribution plot. Plotting failures should not stop analysis."""
    try:
        counts = counts.reindex(range(1, total_groups + 1), fill_value=0)
        fig, ax1 = plt.subplots(figsize=(16, 9))
        max_count = max(counts.values) if len(counts.values) else 0
        ax1.bar(range(1, total_groups + 1), counts.values, alpha=0.6, width=0.6,
                label='Distribution of support counts')
        ax1.set_xlabel(f"Support count (1-{total_groups})", fontsize=14)
        ax1.set_ylabel(ylabel, fontsize=12)
        ax1.set_xticks(range(1, total_groups + 1))
        ax1.set_ylim(0, max_count * 1.1 if max_count > 0 else 1)
        ax1.grid(alpha=0.2, ls=':')

        mean_low = info.get("mean_low", np.nan)
        mean_high = info.get("mean_high", np.nan)
        valley_x = info.get("valley_x", np.nan)
        if np.isfinite(mean_low):
            ax1.axvline(mean_low, lw=1.5, alpha=0.8, label=f'Component 1 mean: {mean_low:.2f}')
        if np.isfinite(mean_high) and mean_high != mean_low:
            ax1.axvline(mean_high, lw=1.5, alpha=0.8, label=f'Component 2 mean: {mean_high:.2f}')
        if np.isfinite(valley_x):
            ax1.axvline(valley_x, ls='--', alpha=0.7, label=f'Theoretical boundary: {valley_x:.2f}')
        if required_count is not None:
            ax1.axvline(required_count, lw=2.5, label=f'Recommended threshold: t={required_count}')
        ax1.legend(loc='upper right', frameon=True)
        plt.title(title, fontsize=14)
        plt.tight_layout()
        plt.savefig(out_file, dpi=300, bbox_inches='tight')
        plt.close()
    except Exception as e:
        print(f"  警告：自动投票阈值分布图保存失败 {out_file}: {e}")
        try:
            plt.close()
        except Exception:
            pass


def _write_auto_vote_summary(records, output_dir, filename):
    os.makedirs(output_dir, exist_ok=True)
    out_file = os.path.join(output_dir, filename)
    pd.DataFrame(records).to_csv(out_file, sep='\t', index=False)
    print(f"自动投票阈值汇总表保存至: {out_file}")


def compute_dmp_vote_thresholds(
        m: int,
        n: int,
        chromosomes: Optional[Union[List[str], str]] = None,
        base_dir: Optional[str] = None,
        output_dir: Optional[str] = None,
        meth_diff_threshold: float = 0.0,
) -> dict:
    """Estimate DMP required_count from consistently filtered FDR tables.

    Each pair contributes one support vote only when both conditions hold:
      1) Qvalue <= that row/pair's Qvalue_Threshold_Used (or context fallback);
      2) abs(MethDiff) >= meth_diff_threshold.

    Therefore, when mianjifa auto-methdiff is enabled, auto-vote is calculated
    from the sites remaining after the newly estimated MethDiff threshold, not
    from stale pairwise DMP files generated before that threshold was known.
    """
    if base_dir is None:
        base_path = Path(__file__).parent
    else:
        base_path = Path(base_dir)

    plot_path = (
        Path(output_dir)
        if output_dir is not None
        else base_path / "and_output" / "auto_vote_thresholds"
    )
    plot_path.mkdir(parents=True, exist_ok=True)

    group_label = f"{m}×{n}"
    total_groups = int(m) * int(n)
    meth_diff_threshold = float(meth_diff_threshold)
    print(
        f"分析 DMP 自动投票阈值：{m} 行 × {n} 列比较组"
        f"（共 {total_groups} 组）"
    )
    print(
        "  数据源 = 完整FDR表；pair支持条件 = "
        f"q通过 且 abs(MethDiff)>={meth_diff_threshold:.6g}"
    )

    group_dirs = []
    for mut_idx in range(1, int(m) + 1):
        for wt_idx in range(1, int(n) + 1):
            dir_path = base_path / f"output_wt{wt_idx}_mut{mut_idx}"
            if not dir_path.is_dir():
                print(f"  ⚠ 警告：目录不存在 {dir_path}")
            group_dirs.append((wt_idx, mut_idx, dir_path))

    if chromosomes is None or chromosomes == "all":
        chrom_filter = None
    else:
        chrom_filter = {
            _mianjifa_normalize_chr(c).lower() for c in chromosomes
        }

    contexts = ["CpG", "CHH", "CHG"]
    thresholds = {}
    records = []
    pair_filter_records = []

    for ctx in contexts:
        print(f"\n▶ 处理 DMP {ctx} 上下文 ...")
        all_records = []

        for group_id, (wt_idx, mut_idx, grp_dir) in enumerate(group_dirs):
            ctx_dir = grp_dir / ctx
            if not ctx_dir.is_dir():
                continue

            fdr_files = sorted(
                p for p in ctx_dir.glob("FDR_corrected_results_*.txt")
                if p.is_file()
            )
            if not fdr_files:
                print(f"  ⚠ {ctx_dir} 中未找到完整FDR文件")
                continue
            if len(fdr_files) > 1:
                print(
                    f"  ⚠ {ctx_dir} 中有多个FDR文件，使用第一个："
                    f"{fdr_files[0].name}"
                )

            fdr_file = fdr_files[0]
            try:
                df = pd.read_csv(fdr_file, sep="\t")
                if len(df.columns) < 6:
                    df = pd.read_csv(
                        fdr_file, sep=r"\s+", engine="python"
                    )
            except Exception as exc:
                print(f"  ⚠ 读取FDR文件失败 {fdr_file}: {exc}")
                continue

            required = {
                "Chromosome", "Position", "Pvalue", "Qvalue", "MethDiff"
            }
            missing = required - set(df.columns)
            if missing:
                print(
                    f"  ⚠ {fdr_file} 缺少列 {sorted(missing)}，跳过"
                )
                continue

            qvalue = pd.to_numeric(df["Qvalue"], errors="coerce")
            methdiff = pd.to_numeric(
                df["MethDiff"], errors="coerce"
            ).abs()
            fallback_q = float(get_dmp_threshold(ctx))
            if "Qvalue_Threshold_Used" in df.columns:
                qthreshold = pd.to_numeric(
                    df["Qvalue_Threshold_Used"], errors="coerce"
                ).fillna(fallback_q)
                qthreshold_source = "Qvalue_Threshold_Used"
            else:
                qthreshold = pd.Series(fallback_q, index=df.index)
                qthreshold_source = "fixed_context_threshold"

            q_pass = qvalue.notna() & qthreshold.notna() & (qvalue <= qthreshold)
            methdiff_pass = methdiff.notna() & (
                methdiff >= meth_diff_threshold
            )
            support_mask = q_pass & methdiff_pass

            selected = df.loc[
                support_mask,
                ["Chromosome", "Position", "Pvalue"],
            ].copy()
            selected["chromosome"] = selected["Chromosome"].map(
                lambda value: "chr" + _mianjifa_normalize_chr(value)
            )
            selected["position"] = pd.to_numeric(
                selected["Position"], errors="coerce"
            )
            selected["pvalue"] = pd.to_numeric(
                selected["Pvalue"], errors="coerce"
            )
            selected = selected.dropna(
                subset=["position", "pvalue"]
            )
            selected["position"] = selected["position"].astype(int)

            if chrom_filter is not None:
                selected = selected[
                    selected["chromosome"].map(
                        lambda x: _mianjifa_normalize_chr(x).lower()
                    ).isin(chrom_filter)
                ]

            selected["group_id"] = group_id
            if not selected.empty:
                all_records.append(
                    selected[[
                        "chromosome", "position", "group_id", "pvalue"
                    ]]
                )

            pair_filter_records.append({
                "context": ctx,
                "wt_replicate": wt_idx,
                "mut_replicate": mut_idx,
                "group_id": group_id,
                "fdr_file": str(fdr_file),
                "total_rows": int(len(df)),
                "q_pass_count": int(q_pass.sum()),
                "methdiff_pass_count": int(methdiff_pass.sum()),
                "pair_support_count": int(support_mask.sum()),
                "meth_diff_threshold_used": meth_diff_threshold,
                "q_threshold_source": qthreshold_source,
                "q_threshold_min": (
                    float(qthreshold.min()) if len(qthreshold) else np.nan
                ),
                "q_threshold_max": (
                    float(qthreshold.max()) if len(qthreshold) else np.nan
                ),
            })

        fallback = _auto_vote_fallback_required_count(total_groups)

        if not all_records:
            print(
                f"  未找到任何满足q与MethDiff条件的 {ctx} 位点，"
                "使用 --vote-threshold 兜底"
            )
            thresholds[ctx] = None
            record = {
                "target": "DMP",
                "context": ctx,
                "total_groups": total_groups,
                "fallback_required_count": fallback,
                "auto_required_count": np.nan,
                "used_required_count": fallback,
                "report_only": bool(AUTO_VOTE_THRESHOLD_REPORT_ONLY),
                "n_candidates": 0,
                "status": "fallback_no_candidates",
                "method": "fallback",
                "data_source": "FDR_q_and_methdiff_filtered",
                "meth_diff_threshold_used": meth_diff_threshold,
                "valley_x": np.nan,
                "mean_low": np.nan,
                "mean_high": np.nan,
            }
            for i in range(1, total_groups + 1):
                record[f"count_{i}"] = 0
            records.append(record)
            continue

        long_df = pd.concat(all_records, ignore_index=True)
        long_df = long_df.drop_duplicates(
            subset=["chromosome", "position", "group_id"],
            keep="first",
        )

        support_long_file = (
            plot_path /
            f"DMP_{ctx}_pair_support_after_methdiff_filter.tsv"
        )
        long_df.to_csv(support_long_file, sep="\t", index=False)

        wide_df = long_df.pivot(
            index=["chromosome", "position"],
            columns="group_id",
            values="pvalue",
        )
        wide_df.reset_index(inplace=True)
        wide_df.columns.name = None

        group_cols = [
            c for c in wide_df.columns
            if isinstance(c, (int, np.integer))
            or (isinstance(c, str) and c.isdigit())
        ]

        if len(group_cols) == 0:
            print(f"  {ctx} 未形成有效group列，使用比例阈值兜底")
            thresholds[ctx] = None
            record = {
                "target": "DMP",
                "context": ctx,
                "total_groups": total_groups,
                "fallback_required_count": fallback,
                "auto_required_count": np.nan,
                "used_required_count": fallback,
                "report_only": bool(AUTO_VOTE_THRESHOLD_REPORT_ONLY),
                "n_candidates": 0,
                "status": "fallback_no_group_columns",
                "method": "fallback",
                "data_source": "FDR_q_and_methdiff_filtered",
                "meth_diff_threshold_used": meth_diff_threshold,
                "valley_x": np.nan,
                "mean_low": np.nan,
                "mean_high": np.nan,
            }
            for i in range(1, total_groups + 1):
                record[f"count_{i}"] = 0
            records.append(record)
            continue

        wide_df["vote_counts"] = wide_df[group_cols].notna().sum(axis=1)
        N = len(group_cols)
        counts = wide_df["vote_counts"].value_counts().reindex(
            range(1, N + 1), fill_value=0
        )
        x = np.repeat(range(1, N + 1), counts.values)

        if len(x) < 2:
            best_t = fallback
            valley_x = np.nan
            means = np.array([np.nan, np.nan])
            method = "fallback_too_few_candidates"
            status = "fallback_too_few_candidates"
            x_plot = None
            pdf = None
        else:
            gmm = GaussianMixture(
                n_components=2, random_state=42, max_iter=500
            )
            gmm.fit(x.reshape(-1, 1))
            means = np.sort(gmm.means_.flatten())

            x_plot = np.linspace(0.5, N + 0.5, 3000)
            pdf = np.exp(gmm.score_samples(x_plot.reshape(-1, 1)))
            search = (x_plot > means[0]) & (x_plot < means[1])
            if np.any(search):
                valley_idx = np.argmin(pdf[search])
                valley_x = float(x_plot[search][valley_idx])
                best_t = int(np.ceil(valley_x))
                method = "two_component_gmm_valley"
                status = "ok"
            else:
                total = counts.sum()
                indices = np.arange(1, N + 1)
                sum_val = np.sum(indices * counts)
                sumc = np.cumsum(counts)
                weight0 = sumc / total
                weight1 = 1 - weight0
                mean0 = np.cumsum(indices * counts) / (sumc + 1e-10)
                mean1 = (
                    sum_val - np.cumsum(indices * counts)
                ) / (total - sumc + 1e-10)
                var_between = weight0 * weight1 * (mean0 - mean1) ** 2
                best_t = int(np.argmax(var_between) + 1)
                valley_x = float(best_t - 0.5)
                method = "otsu_fallback"
                status = "ok"
            best_t = max(1, min(int(best_t), total_groups))

        img_file = (
            plot_path /
            f"DMP_{ctx}_vote_support_distribution_{total_groups}_groups.jpeg"
        )
        try:
            fig, ax1 = plt.subplots(figsize=(16, 9))
            max_count = max(counts.values) if len(counts.values) else 0
            ax1.bar(
                range(1, N + 1), counts.values,
                alpha=0.6, color="skyblue", width=0.6,
                label="Distribution of support counts",
            )
            ax1.set_xlabel(f"Support count (1-{total_groups})", fontsize=14)
            ax1.set_ylabel("Number of DMPs", fontsize=12)
            ax1.set_xticks(range(1, N + 1))
            ax1.set_ylim(0, max_count * 1.1 if max_count > 0 else 1)
            ax1.grid(alpha=0.2, ls=":")

            if np.isfinite(means[0]):
                ax1.axvline(means[0], color="darkblue", lw=1.5, alpha=0.8)
                ax1.text(
                    means[0], max_count * 0.85 if max_count > 0 else 0.85,
                    f"Noise Component: {means[0]:.2f}", ha="center",
                    color="darkblue", fontsize=11,
                    bbox=dict(facecolor="white", alpha=0.7),
                )
            if np.isfinite(means[1]):
                ax1.axvline(means[1], color="darkred", lw=1.5, alpha=0.8)
                ax1.text(
                    means[1], max_count * 0.85 if max_count > 0 else 0.85,
                    f"Signal Component: {means[1]:.2f}", ha="center",
                    color="darkred", fontsize=11,
                    bbox=dict(facecolor="white", alpha=0.7),
                )
            if np.isfinite(valley_x):
                ax1.axvline(
                    valley_x, color="purple", ls="--", alpha=0.7,
                    label=f"Theoretical boundary: {valley_x:.2f}",
                )
            ax1.axvline(
                best_t, color="green", lw=2.5,
                label=f"Recommended threshold: t={best_t}",
            )

            if x_plot is not None and pdf is not None:
                ax2 = ax1.twinx()
                ax2.plot(x_plot, pdf, "r-", lw=2, label="GMM-fitted density")
                ax2.set_ylabel("Probability density", fontsize=12)
                lines1, lab1 = ax1.get_legend_handles_labels()
                lines2, lab2 = ax2.get_legend_handles_labels()
                ax1.legend(
                    lines1 + lines2, lab1 + lab2,
                    loc="upper right", frameon=True,
                )
            else:
                ax1.legend(loc="upper right", frameon=True)

            plt.title(
                f"DMP ({ctx}) support counts in {group_label}; "
                f"MethDiff >= {meth_diff_threshold:.4g}",
                fontsize=14,
            )
            plt.tight_layout()
            plt.savefig(img_file, dpi=300, bbox_inches="tight")
            plt.close()
        except Exception as e:
            print(f"  ⚠ {ctx} 自动投票阈值图片保存失败 {img_file}: {e}")
            try:
                plt.close()
            except Exception:
                pass

        thresholds[ctx] = int(best_t)
        used = (
            fallback
            if AUTO_VOTE_THRESHOLD_REPORT_ONLY or best_t is None
            else int(best_t)
        )
        counts_full = counts.reindex(
            range(1, total_groups + 1), fill_value=0
        )
        record = {
            "target": "DMP",
            "context": ctx,
            "total_groups": total_groups,
            "fallback_required_count": fallback,
            "auto_required_count": (
                int(best_t) if best_t is not None else np.nan
            ),
            "used_required_count": used,
            "report_only": bool(AUTO_VOTE_THRESHOLD_REPORT_ONLY),
            "n_candidates": int(counts.sum()),
            "status": status,
            "method": method,
            "data_source": "FDR_q_and_methdiff_filtered",
            "meth_diff_threshold_used": meth_diff_threshold,
            "valley_x": valley_x,
            "mean_low": float(means[0]) if np.isfinite(means[0]) else np.nan,
            "mean_high": float(means[1]) if np.isfinite(means[1]) else np.nan,
        }
        for i in range(1, total_groups + 1):
            record[f"count_{i}"] = int(counts_full.loc[i])
        records.append(record)
        print(
            f"  ✓ {ctx}: 推荐阈值={best_t}；"
            f"MethDiff过滤阈值={meth_diff_threshold:.6g}；"
            f"图片={img_file}"
        )

    _write_auto_vote_summary(
        records, str(plot_path), "DMP_vote_threshold_summary.tsv"
    )
    pd.DataFrame(pair_filter_records).to_csv(
        plot_path / "DMP_pair_filtering_summary.tsv",
        sep="\t", index=False,
    )
    print(
        "DMP pair过滤汇总表保存至: "
        f"{plot_path / 'DMP_pair_filtering_summary.tsv'}"
    )
    return thresholds


def compute_dmr_vote_thresholds(
        m: int,
        n: int,
        chromosomes: Optional[Union[List[str], str]] = None,
        base_dir: Optional[str] = None,
        output_dir: Optional[str] = None
) -> dict:
    """
    自动读取 m×n 个比较组 and_output/dmr_analysis_wt<wt>_mut<mut>/<context>/dmr_fisher_significant_*.txt，
    筛选 qvalue <= DMR_QVALUE_THRESHOLD 的显著 DMR，合并后计算每个 DMR 的支持次数，
    使用双峰 GMM 拟合分布并自动选择 DMR 投票阈值。
    同时保存 DMR 支持次数分布图和 DMR_vote_threshold_summary.tsv。

    特殊处理：
      1. 单峰分布：阈值使用当前 --vote-threshold 对应的 fallback required_count；
      2. 截断型分布：低支持次数段全无数据时，将 GMM 阈值减 1。
    """
    if base_dir is None:
        base_path = Path(__file__).parent
    else:
        base_path = Path(base_dir)

    # 兼容两种情况：
    # 1) base_dir 是运行根目录，DMR 文件在 base_dir/and_output/dmr_analysis_...
    # 2) base_dir 已经是 and_output，DMR 文件在 base_dir/dmr_analysis_...
    and_output_path = base_path / 'and_output'
    if not and_output_path.is_dir():
        and_output_path = base_path

    plot_path = Path(output_dir) if output_dir is not None else and_output_path / 'auto_vote_thresholds'
    plot_path.mkdir(parents=True, exist_ok=True)

    group_label = f"{m}×{n}"
    total_groups = int(m) * int(n)
    print(f"分析 DMR 自动投票阈值：{m} 行 × {n} 列 比较组（共 {total_groups} 组）")

    group_dirs = []
    for mut_idx in range(1, int(m) + 1):
        for wt_idx in range(1, int(n) + 1):
            dir_path = and_output_path / f"dmr_analysis_wt{wt_idx}_mut{mut_idx}"
            if not dir_path.is_dir():
                print(f"  ⚠ 警告：目录不存在 {dir_path}")
            group_dirs.append(dir_path)

    if chromosomes is None or chromosomes == "all":
        chrom_filter = None
    else:
        chrom_filter = {str(c).lower() for c in chromosomes}

    contexts = ["CpG", "CHH", "CHG"]
    thresholds = {}
    records = []

    for ctx in contexts:
        print(f"\n▶ 处理 DMR {ctx} 上下文 ...")
        all_records = []
        fallback = _auto_vote_fallback_required_count(total_groups)

        for group_id, grp_dir in enumerate(group_dirs):
            ctx_dir = grp_dir / ctx
            if not ctx_dir.is_dir():
                continue

            dmr_files = [
                f for f in ctx_dir.iterdir()
                if f.is_file() and re.match(
                    r'dmr_fisher_significant_chr[0-9A-Za-z]+\.txt$',
                    f.name,
                    re.IGNORECASE
                )
            ]

            for file in dmr_files:
                chrom_match = re.search(r'chr[0-9A-Za-z]+', file.name, re.IGNORECASE)
                if not chrom_match:
                    continue

                chrom = chrom_match.group(0)
                if chrom_filter is not None and chrom.lower() not in chrom_filter:
                    continue

                try:
                    df = pd.read_csv(file, sep='\t')
                    if df.empty:
                        continue

                    if 'qvalue' in df.columns:
                        df = df[df['qvalue'] <= DMR_QVALUE_THRESHOLD].copy()

                    if df.empty:
                        continue

                    df['chromosome'] = chrom
                    df['group_id'] = group_id
                    all_records.append(df[['chromosome', 'DMR_start', 'DMR_end', 'group_id']])

                except Exception as e:
                    print(f"  ⚠ 警告：读取 DMR 文件失败 {file}: {e}")
                    continue

        if not all_records:
            print(f"  未找到任何显著 {ctx} DMR，使用 --vote-threshold 兜底")
            thresholds[ctx] = None

            record = {
                'target': 'DMR',
                'context': ctx,
                'total_groups': total_groups,
                'fallback_required_count': fallback,
                'auto_required_count': np.nan,
                'used_required_count': fallback,
                'report_only': bool(AUTO_VOTE_THRESHOLD_REPORT_ONLY),
                'n_candidates': 0,
                'status': 'fallback_no_candidates',
                'method': 'fallback',
                'valley_x': np.nan,
                'mean_low': np.nan,
                'mean_high': np.nan,
            }
            for i in range(1, total_groups + 1):
                record[f'count_{i}'] = 0
            records.append(record)
            continue

        long_df = pd.concat(all_records, ignore_index=True)
        long_df = long_df.drop_duplicates(
            subset=['chromosome', 'DMR_start', 'DMR_end', 'group_id'],
            keep='first'
        )
        long_df['present'] = 1

        wide_df = long_df.pivot(
            index=['chromosome', 'DMR_start', 'DMR_end'],
            columns='group_id',
            values='present'
        )
        wide_df.reset_index(inplace=True)
        wide_df.columns.name = None

        group_cols = [
            c for c in wide_df.columns
            if isinstance(c, (int, np.integer)) or (isinstance(c, str) and c.isdigit())
        ]

        if len(group_cols) == 0:
            print(f"  {ctx} 未形成有效 group 列，使用 --vote-threshold 兜底")
            thresholds[ctx] = None

            record = {
                'target': 'DMR',
                'context': ctx,
                'total_groups': total_groups,
                'fallback_required_count': fallback,
                'auto_required_count': np.nan,
                'used_required_count': fallback,
                'report_only': bool(AUTO_VOTE_THRESHOLD_REPORT_ONLY),
                'n_candidates': 0,
                'status': 'fallback_no_group_columns',
                'method': 'fallback',
                'valley_x': np.nan,
                'mean_low': np.nan,
                'mean_high': np.nan,
            }
            for i in range(1, total_groups + 1):
                record[f'count_{i}'] = 0
            records.append(record)
            continue

        wide_df['vote_counts'] = wide_df[group_cols].notna().sum(axis=1)
        counts = wide_df['vote_counts'].value_counts().reindex(
            range(1, total_groups + 1),
            fill_value=0
        )

        nonzero_counts = counts[counts > 0]

        # ========== 特殊处理：单峰分布 ==========
        if len(nonzero_counts) == 1:
            single_value = int(nonzero_counts.index[0])
            best_t = int(fallback)
            valley_x = np.nan
            means = np.array([float(single_value), float(single_value)])
            status = 'fallback_single_support_count'
            method = 'fallback_vote_threshold'

            print(
                f"  ⚠ 检测到单峰分布：所有 DMR 的支持次数均为 {single_value}，"
                f"采用 fallback 阈值 t={best_t}"
            )

            thresholds[ctx] = best_t
            used = fallback if AUTO_VOTE_THRESHOLD_REPORT_ONLY or best_t is None else best_t

            record = {
                'target': 'DMR',
                'context': ctx,
                'total_groups': total_groups,
                'fallback_required_count': fallback,
                'auto_required_count': best_t,
                'used_required_count': used,
                'report_only': bool(AUTO_VOTE_THRESHOLD_REPORT_ONLY),
                'n_candidates': int(counts.sum()),
                'status': status,
                'method': method,
                'valley_x': valley_x,
                'mean_low': float(means[0]),
                'mean_high': float(means[1]),
            }
            for i in range(1, total_groups + 1):
                record[f'count_{i}'] = int(counts.loc[i])
            records.append(record)
            continue

        # ========== 正常双峰拟合 ==========
        x = np.repeat(range(1, total_groups + 1), counts.values)

        if len(x) < 2:
            best_t = int(fallback)
            valley_x = np.nan
            means = np.array([np.nan, np.nan])
            status = 'fallback_too_few_candidates'
            method = 'fallback_vote_threshold'
            x_plot = None
            pdf = None
        else:
            gmm = GaussianMixture(n_components=2, random_state=42, max_iter=500)
            gmm.fit(x.reshape(-1, 1))
            means = np.sort(gmm.means_.flatten())

            x_plot = np.linspace(0.5, total_groups + 0.5, 3000)
            pdf = np.exp(gmm.score_samples(x_plot.reshape(-1, 1)))

            search = (x_plot > means[0]) & (x_plot < means[1])
            if np.any(search):
                valley_idx = np.argmin(pdf[search])
                valley_x = float(x_plot[search][valley_idx])
                best_t = int(np.ceil(valley_x))
                status = 'ok'
                method = 'two_component_gmm_valley'
            else:
                total = counts.sum()
                indices = np.arange(1, total_groups + 1)
                sum_val = np.sum(indices * counts)
                sumc = np.cumsum(counts)
                weight0 = sumc / total
                weight1 = 1 - weight0
                mean0 = np.cumsum(indices * counts) / (sumc + 1e-10)
                mean1 = (sum_val - np.cumsum(indices * counts)) / (total - sumc + 1e-10)
                var_between = weight0 * weight1 * (mean0 - mean1) ** 2
                best_t = int(np.argmax(var_between) + 1)
                valley_x = float(best_t - 0.5)
                status = 'ok'
                method = 'otsu_fallback'

            # ========== 特殊处理：截断型分布 ==========
            half_limit = total_groups // 2
            low_range_zero = all(counts.loc[i] == 0 for i in range(1, half_limit + 1))

            if low_range_zero:
                original_t = int(best_t)
                best_t = max(1, int(best_t) - 1)
                status = 'ok_truncated_adjusted'
                method = f'{method}_truncated_adjust_{original_t}_to_{best_t}'
                print(
                    f"  ⚠ 检测到截断型分布：支持次数 1 到 {half_limit} 均无 DMR，"
                    f"将阈值从 {original_t} 调整为 {best_t}"
                )

            best_t = max(1, min(int(best_t), total_groups))

        # ========== 绘图 ==========
        img_file = plot_path / f"DMR_{ctx}_vote_support_distribution_{total_groups}_groups.jpeg"

        try:
            fig, ax1 = plt.subplots(figsize=(16, 9))
            max_count = max(counts.values) if len(counts.values) else 0

            ax1.bar(
                range(1, total_groups + 1),
                counts.values,
                alpha=0.6,
                color='skyblue',
                width=0.6,
                label='Distribution of support counts'
            )
            ax1.set_xlabel(f"Support count (1-{total_groups})", fontsize=14)
            ax1.set_ylabel("Number of DMRs", fontsize=12)
            ax1.set_xticks(range(1, total_groups + 1))
            ax1.set_ylim(0, max_count * 1.1 if max_count > 0 else 1)
            ax1.grid(alpha=0.2, ls=':')

            if np.isfinite(means[0]):
                ax1.axvline(means[0], color='darkblue', lw=1.5, alpha=0.8)
                ax1.text(
                    means[0],
                    max_count * 0.85 if max_count > 0 else 0.85,
                    f'Noise Component: {means[0]:.2f}',
                    ha='center',
                    color='darkblue',
                    fontsize=11,
                    bbox=dict(facecolor='white', alpha=0.7)
                )

            if np.isfinite(means[1]):
                ax1.axvline(means[1], color='darkred', lw=1.5, alpha=0.8)
                ax1.text(
                    means[1],
                    max_count * 0.85 if max_count > 0 else 0.85,
                    f'Signal Component: {means[1]:.2f}',
                    ha='center',
                    color='darkred',
                    fontsize=11,
                    bbox=dict(facecolor='white', alpha=0.7)
                )

            if np.isfinite(valley_x):
                ax1.axvline(
                    valley_x,
                    color='purple',
                    ls='--',
                    alpha=0.7,
                    label=f'Theoretical boundary (minimum point): {valley_x:.2f}'
                )

            ax1.axvline(
                best_t,
                color='green',
                lw=2.5,
                label=f'Recommended threshold: t={best_t}'
            )

            if x_plot is not None and pdf is not None:
                ax2 = ax1.twinx()
                ax2.plot(x_plot, pdf, 'r-', lw=2, label='GMM-fitted density')
                ax2.set_ylabel("Probability density", fontsize=12)

                lines1, lab1 = ax1.get_legend_handles_labels()
                lines2, lab2 = ax2.get_legend_handles_labels()
                ax1.legend(lines1 + lines2, lab1 + lab2, loc='upper right', frameon=True)
            else:
                ax1.legend(loc='upper right', frameon=True)

            plt.title(
                f"Distribution fitting of DMR ({ctx}) support counts in {group_label} comparisons",
                fontsize=14
            )
            plt.tight_layout()
            plt.savefig(img_file, dpi=300, bbox_inches='tight')
            plt.close()

        except Exception as e:
            print(f"  ⚠ 警告：{ctx} 自动投票阈值图片保存失败 {img_file}: {e}")
            try:
                plt.close()
            except Exception:
                pass

        thresholds[ctx] = int(best_t)
        used = fallback if AUTO_VOTE_THRESHOLD_REPORT_ONLY or best_t is None else int(best_t)

        record = {
            'target': 'DMR',
            'context': ctx,
            'total_groups': total_groups,
            'fallback_required_count': fallback,
            'auto_required_count': int(best_t) if best_t is not None else np.nan,
            'used_required_count': used,
            'report_only': bool(AUTO_VOTE_THRESHOLD_REPORT_ONLY),
            'n_candidates': int(counts.sum()),
            'status': status,
            'method': method,
            'valley_x': valley_x,
            'mean_low': float(means[0]) if np.isfinite(means[0]) else np.nan,
            'mean_high': float(means[1]) if np.isfinite(means[1]) else np.nan,
        }

        for i in range(1, total_groups + 1):
            record[f'count_{i}'] = int(counts.loc[i])

        records.append(record)

        print(f"  ✓ {ctx}: 推荐阈值 = {best_t}，图片已保存至 {img_file}")

    _write_auto_vote_summary(records, str(plot_path), 'DMR_vote_threshold_summary.tsv')
    return thresholds


def update_auto_dmp_vote_thresholds(m, n, work_dir='.', meth_diff_threshold=0.0):
    global AUTO_DMP_VOTE_THRESHOLDS
    if not (AUTO_DMP_VOTE_THRESHOLD or AUTO_VOTE_THRESHOLD_REPORT_ONLY):
        AUTO_DMP_VOTE_THRESHOLDS = {'CpG': None, 'CHH': None, 'CHG': None}
        print("DMP自动投票阈值未启用，使用 --vote-threshold")
        return AUTO_DMP_VOTE_THRESHOLDS
    AUTO_DMP_VOTE_THRESHOLDS = compute_dmp_vote_thresholds(
        m=m, n=n, chromosomes='all', base_dir=work_dir,
        output_dir=os.path.join(work_dir, 'and_output', 'auto_vote_thresholds'),
        meth_diff_threshold=meth_diff_threshold
    )
    print(f"DMP自动投票阈值: {AUTO_DMP_VOTE_THRESHOLDS}")
    return AUTO_DMP_VOTE_THRESHOLDS


def update_auto_dmr_vote_thresholds(m, n, work_dir='.'):
    global AUTO_DMR_VOTE_THRESHOLDS
    if not (AUTO_DMR_VOTE_THRESHOLD or AUTO_VOTE_THRESHOLD_REPORT_ONLY):
        AUTO_DMR_VOTE_THRESHOLDS = {'CpG': None, 'CHH': None, 'CHG': None}
        print("DMR自动投票阈值未启用，使用 --vote-threshold")
        return AUTO_DMR_VOTE_THRESHOLDS
    AUTO_DMR_VOTE_THRESHOLDS = compute_dmr_vote_thresholds(
        m=m, n=n, chromosomes='all', base_dir=work_dir,
        output_dir=os.path.join(work_dir, 'and_output', 'auto_vote_thresholds')
    )
    print(f"DMR自动投票阈值: {AUTO_DMR_VOTE_THRESHOLDS}")
    return AUTO_DMR_VOTE_THRESHOLDS


DmrRecord = namedtuple('DmrRecord', ['exp_methy', 'exp_unmethy', 'wild_methy', 'wild_unmethy', 'qvalue', 'direction'])


def _apply_worker_config(config):
    """Apply CLI-derived global settings inside a worker process."""
    global DMR_QVALUE_THRESHOLD, VOTE_THRESHOLD, DMR_ENGINE
    global AUTO_DMP_VOTE_THRESHOLD, AUTO_DMR_VOTE_THRESHOLD, AUTO_VOTE_THRESHOLD_REPORT_ONLY
    global AUTO_DMP_VOTE_THRESHOLDS, AUTO_DMR_VOTE_THRESHOLDS
    global AUTO_QVALUE_TWOSTEP, AUTO_QVALUE_REPORT_ONLY
    global AUTO_QVALUE_P_CUTOFF, AUTO_QVALUE_MIN_CANDIDATES
    global AUTO_QVALUE_USE_SMOOTH, AUTO_QVALUE_SMOOTH_SIGMA
    global DMP_LOWDIFF_STRICT_VOTE, DMP_LOWDIFF_CUTOFF, DMP_LOWDIFF_STRICT_VOTE_REPORT_ONLY
    if not config:
        return
    DMP_QVALUE_THRESHOLDS.update(config.get("dmp_qvalue_thresholds", {}))
    DMR_QVALUE_THRESHOLD = config.get("dmr_qvalue_threshold", DMR_QVALUE_THRESHOLD)
    VOTE_THRESHOLD = config.get("vote_threshold", VOTE_THRESHOLD)
    AUTO_DMP_VOTE_THRESHOLD = config.get("auto_dmp_vote_threshold", AUTO_DMP_VOTE_THRESHOLD)
    AUTO_DMR_VOTE_THRESHOLD = config.get("auto_dmr_vote_threshold", AUTO_DMR_VOTE_THRESHOLD)
    AUTO_VOTE_THRESHOLD_REPORT_ONLY = config.get("auto_vote_threshold_report_only", AUTO_VOTE_THRESHOLD_REPORT_ONLY)
    AUTO_DMP_VOTE_THRESHOLDS.update(config.get("auto_dmp_vote_thresholds", {}))
    AUTO_DMR_VOTE_THRESHOLDS.update(config.get("auto_dmr_vote_thresholds", {}))
    DMR_ENGINE = config.get("dmr_engine", DMR_ENGINE)
    AUTO_QVALUE_TWOSTEP = config.get("auto_qvalue_twostep", AUTO_QVALUE_TWOSTEP)
    AUTO_QVALUE_REPORT_ONLY = config.get("auto_qvalue_report_only", AUTO_QVALUE_REPORT_ONLY)
    AUTO_QVALUE_P_CUTOFF = config.get("auto_qvalue_p_cutoff", AUTO_QVALUE_P_CUTOFF)
    AUTO_QVALUE_MIN_CANDIDATES = config.get("auto_qvalue_min_candidates", AUTO_QVALUE_MIN_CANDIDATES)
    AUTO_QVALUE_USE_SMOOTH = config.get("auto_qvalue_use_smooth", AUTO_QVALUE_USE_SMOOTH)
    AUTO_QVALUE_SMOOTH_SIGMA = config.get("auto_qvalue_smooth_sigma", AUTO_QVALUE_SMOOTH_SIGMA)
    DMP_LOWDIFF_STRICT_VOTE = config.get("dmp_lowdiff_strict_vote", DMP_LOWDIFF_STRICT_VOTE)
    DMP_LOWDIFF_CUTOFF = config.get("dmp_lowdiff_cutoff", DMP_LOWDIFF_CUTOFF)
    DMP_LOWDIFF_STRICT_VOTE_REPORT_ONLY = config.get("dmp_lowdiff_strict_vote_report_only", DMP_LOWDIFF_STRICT_VOTE_REPORT_ONLY)


def _run_replicate_pair_worker(payload):
    """Worker for one replicate-pair task."""
    _apply_worker_config(payload.get("config"))
    replicate_x = payload["replicate_x"]
    replicate_y = payload["replicate_y"]
    log_file = payload.get("log_file")
    start = time.time()

    def _run():
        success_count, test_count = process_replicate_pair(
            replicate_x=replicate_x,
            replicate_y=replicate_y,
            files1=payload["files1"],
            files2=payload["files2"],
            dir1=payload["dir1"],
            dir2=payload["dir2"],
            dir1_name=payload["dir1_name"],
            dir2_name=payload["dir2_name"],
            unfilter_mtypes=payload["unfilter_mtypes"],
            work_dir=payload.get("work_dir", "."),
            meth_diff_threshold=payload.get("meth_diff_threshold", 0.0),
            skip_dmr=payload.get("skip_dmr", False),
            skip_window=payload.get("skip_window", False),
        )
        return success_count, test_count

    if log_file:
        os.makedirs(os.path.dirname(log_file), exist_ok=True)
        with open(log_file, "w") as lf, contextlib.redirect_stdout(lf), contextlib.redirect_stderr(lf):
            print(f"[START] wt{replicate_y}_mut{replicate_x} {time.ctime()}")
            result = _run()
            print(f"[END] wt{replicate_y}_mut{replicate_x} {time.ctime()}")
    else:
        result = _run()

    elapsed = time.time() - start
    return {
        "replicate_x": replicate_x,
        "replicate_y": replicate_y,
        "success_count": int(result[0]),
        "test_count": int(result[1]),
        "elapsed": elapsed,
        "log_file": log_file,
    }


def _run_common_sites_to_dmr_worker(payload):
    """Worker for common-sites-to-DMR candidate generation of one methylation context."""
    _apply_worker_config(payload.get("config"))
    mtype = payload["methylation_type"]
    log_file = payload.get("log_file")
    start = time.time()

    def _run():
        return process_common_sites_to_dmr(methylation_type=mtype, work_dir=payload.get("work_dir", "."))

    if log_file:
        os.makedirs(os.path.dirname(log_file), exist_ok=True)
        with open(log_file, "w") as lf, contextlib.redirect_stdout(lf), contextlib.redirect_stderr(lf):
            print(f"[START] common DMR candidate {mtype} {time.ctime()}")
            result = _run()
            print(f"[END] common DMR candidate {mtype} {time.ctime()}")
    else:
        result = _run()

    return {
        "methylation_type": mtype,
        "ok": result is not None,
        "n_chromosomes": len(result) if result else 0,
        "elapsed": time.time() - start,
        "log_file": log_file,
    }


def _run_summarize_dmr_worker(payload):
    """Worker for common DMR reads summation/Fisher/FDR for one pair × context task."""
    _apply_worker_config(payload.get("config"))
    log_file = payload.get("log_file")
    start = time.time()
    replicate_x = payload["replicate_x"]
    replicate_y = payload["replicate_y"]
    mtype = payload["methylation_type"]

    def _run():
        summarize_dmr_methylation(
            methy_dir=payload["methy_output_dir"],
            replicate_x=replicate_x,
            replicate_y=replicate_y,
            file1_path=payload["file1_path"],
            file2_path=payload["file2_path"],
            methylation_type=mtype,
            custom_dmr_dir=payload["custom_dmr_dir"],
        )

    if log_file:
        os.makedirs(os.path.dirname(log_file), exist_ok=True)
        with open(log_file, "w") as lf, contextlib.redirect_stdout(lf), contextlib.redirect_stderr(lf):
            print(f"[START] summarize DMR wt{replicate_y}_mut{replicate_x} {mtype} {time.ctime()}")
            _run()
            print(f"[END] summarize DMR wt{replicate_y}_mut{replicate_x} {mtype} {time.ctime()}")
    else:
        _run()

    return {
        "replicate_x": replicate_x,
        "replicate_y": replicate_y,
        "methylation_type": mtype,
        "methy_output_dir": payload["methy_output_dir"],
        "elapsed": time.time() - start,
        "log_file": log_file,
    }



def _run_newtoboth_worker(payload):
    """Worker for converting one raw methylation file to bothMeUnme matrix files."""
    filepath = payload["filepath"]
    output_dir = payload["output_dir"]
    num = payload["num"]
    chr_series = payload["chr_series"]
    label = payload.get("label", os.path.basename(filepath))
    log_file = payload.get("log_file")
    start = time.time()

    def _run():
        print(f"处理文件{filepath}")
        single_newtoboth(filepath, output_dir, num, chr_series)

    if log_file:
        os.makedirs(os.path.dirname(log_file), exist_ok=True)
        with open(log_file, "w") as lf, contextlib.redirect_stdout(lf), contextlib.redirect_stderr(lf):
            print(f"[START] newtoboth {label} {time.ctime()}")
            _run()
            print(f"[END] newtoboth {label} {time.ctime()}")
    else:
        _run()

    return {
        "label": label,
        "filepath": filepath,
        "output_dir": output_dir,
        "num": int(num),
        "elapsed": time.time() - start,
        "log_file": log_file,
    }

def _build_parallel_config():
    """Collect current global parameters to pass into worker processes."""
    return {
        "dmp_qvalue_thresholds": dict(DMP_QVALUE_THRESHOLDS),
        "dmr_qvalue_threshold": DMR_QVALUE_THRESHOLD,
        "vote_threshold": VOTE_THRESHOLD,
        "auto_dmp_vote_threshold": AUTO_DMP_VOTE_THRESHOLD,
        "auto_dmr_vote_threshold": AUTO_DMR_VOTE_THRESHOLD,
        "auto_vote_threshold_report_only": AUTO_VOTE_THRESHOLD_REPORT_ONLY,
        "auto_dmp_vote_thresholds": dict(AUTO_DMP_VOTE_THRESHOLDS),
        "auto_dmr_vote_thresholds": dict(AUTO_DMR_VOTE_THRESHOLDS),
        "dmr_engine": DMR_ENGINE,
        "auto_qvalue_twostep": AUTO_QVALUE_TWOSTEP,
        "auto_qvalue_report_only": AUTO_QVALUE_REPORT_ONLY,
        "auto_qvalue_p_cutoff": AUTO_QVALUE_P_CUTOFF,
        "auto_qvalue_min_candidates": AUTO_QVALUE_MIN_CANDIDATES,
        "auto_qvalue_use_smooth": AUTO_QVALUE_USE_SMOOTH,
        "auto_qvalue_smooth_sigma": AUTO_QVALUE_SMOOTH_SIGMA,
        "dmp_lowdiff_strict_vote": DMP_LOWDIFF_STRICT_VOTE,
        "dmp_lowdiff_cutoff": DMP_LOWDIFF_CUTOFF,
        "dmp_lowdiff_strict_vote_report_only": DMP_LOWDIFF_STRICT_VOTE_REPORT_ONLY,
    }



def process_common_sites_dmr_and_summarize(dir1, dir2, m, n, methylation_types=['CpG', 'CHH', 'CHG'], work_dir=".", threads=1):
    """
    完整流程：处理 common sites DMR + 循环求和分析 + 汇总

    并行策略：
      1. common significant sites -> common DMR candidate：按 methylation type 并行
      2. common DMR 在各 replicate pair 中 reads 汇总/Fisher/FDR：在 summarize_all_dmr_methylation 中并行
    """
    print("第三阶段：处理共同显著位点 DMR")

    threads = max(1, int(threads))
    config = _build_parallel_config()

    if threads <= 1 or len(methylation_types) <= 1:
        for mtype in methylation_types:
            dmr_results = process_common_sites_to_dmr(methylation_type=mtype, work_dir=work_dir)
            if dmr_results:
                print(f"\n{mtype} 类型 DMR 分析完成，共 {len(dmr_results)} 个染色体有DMR结果")
    else:
        max_workers = min(threads, len(methylation_types))
        log_dir = os.path.join(work_dir, "parallel_logs", "common_dmr_candidates")
        os.makedirs(log_dir, exist_ok=True)
        tasks = []
        for mtype in methylation_types:
            tasks.append({
                "methylation_type": mtype,
                "work_dir": work_dir,
                "config": config,
                "log_file": os.path.join(log_dir, f"common_dmr_candidate_{mtype}.log"),
            })

        print(f"启用并行 common DMR candidate 处理: workers={max_workers}, tasks={len(tasks)}")
        with ProcessPoolExecutor(max_workers=max_workers) as executor:
            future_to_task = {executor.submit(_run_common_sites_to_dmr_worker, task): task for task in tasks}
            for idx, future in enumerate(as_completed(future_to_task), 1):
                task = future_to_task[future]
                mtype = task["methylation_type"]
                try:
                    res = future.result()
                    print(
                        f"[DONE] {idx}/{len(tasks)} common DMR {mtype}: "
                        f"{res['n_chromosomes']} chromosomes, {res['elapsed']:.2f}s, log={res['log_file']}"
                    )
                except Exception as e:
                    print(f"[FAILED] common DMR {mtype}: {e}")
                    raise

    # 第二步：对所有 output_x_y 的 DMR 进行求和分析并汇总
    summarize_all_dmr_methylation(dir1, dir2, m, n, methylation_types, work_dir=work_dir, threads=threads)

def generate_final_significant_dmr(dmr_data_dict, methylation_type, output_dir,m,n,dir1,dir2):
    """
    生成最终的显著 DMR 文件（基于贝叶斯判定）

    参数：
        dmr_data_dict: 染色体->DMR区域(start,end)->数据列表(一系列每个位点的统计数据) 的映射字典
        methylation_type: 甲基化类型
        output_dir: 输出目录
        threshold: 贝叶斯判定阈值（默认 2/3）
    """
    print(f"\n处理 {methylation_type} 的最终 DMR 汇总...")

    if not dmr_data_dict:
        print(f"  {methylation_type} 无 DMR 数据")
        return

    auto_vote_threshold = AUTO_DMR_VOTE_THRESHOLDS.get(methylation_type)
    if AUTO_VOTE_THRESHOLD_REPORT_ONLY:
        auto_vote_threshold = None
    if auto_vote_threshold is not None:
        print(f"  使用 {methylation_type} DMR自动投票阈值: {auto_vote_threshold}")
    else:
        print(f"  {methylation_type} DMR自动投票阈值不可用或仅报告，使用比例阈值 {VOTE_THRESHOLD}")

    final_dmr_list = []

    # 遍历所有染色体
    for chr_num in sorted(dmr_data_dict.keys()):
        print(f"  处理染色体 {chr_num}...")
        chr_dmr_dict = dmr_data_dict[chr_num]  # 获取当前染色体中的该映射：DMR区域(start,end)->数据列表(一系列每个位点的统计数据)
                                                                    # 此处数据列表中有m*n个元素，代表每次检验中该区域的统计信息

        # 遍历该染色体的所有 DMR 区域
        for dmr_key, data_list in chr_dmr_dict.items():
            start, end = dmr_key

            # 统计显著次数
            sig_count = sum(1 for item in data_list if item.qvalue <= DMR_QVALUE_THRESHOLD)
            total_count = len(data_list)

            # 判定
            is_significant = bayes_deciding(sig_count, total_count - sig_count, auto_vote_threshold=auto_vote_threshold)

            if not is_significant:
                continue

            # 计算平均值
            avg_exp_m = np.mean([item.exp_methy for item in data_list])
            avg_exp_u = np.mean([item.exp_unmethy for item in data_list])
            avg_wild_m = np.mean([item.wild_methy for item in data_list])
            avg_wild_u = np.mean([item.wild_unmethy for item in data_list])

            filtered_values = [item.qvalue for item in data_list if item.qvalue <= DMR_QVALUE_THRESHOLD]
            avg_qvalue = np.mean(filtered_values) if filtered_values else 1

            # 投票决定 direction（多数投票）
            direction_votes = [item.direction for item in data_list]
            direction = 1 if sum(direction_votes) >= len(direction_votes) / 2 else 0

            # 计算概率
            prob = sig_count / total_count

            dmr_record = {
                'Chromosome': chr_num,
                'Methylation_Type': methylation_type,
                'DMR_start': start,
                'DMR_end': end,
                'Length': end - start + 1,
                'Avg_exp_methy': avg_exp_m,
                'Avg_exp_unmethy': avg_exp_u,
                'Avg_wild_methy': avg_wild_m,
                'Avg_wild_unmethy': avg_wild_u,
                'Sig_Avg_qvalue': avg_qvalue,
                'Direction': direction,
                'Sig_count': sig_count,
                'Total_count': total_count,
                'Sig_probability': prob
            }
            # 添加所有replicate的qvalue（按顺序）
            idx = 0
            for x in range(1, m + 1):
                for y in range(1, n + 1):
                    col_name = f'qvalue_{os.path.basename(dir2.rstrip("/"))}{y}_{os.path.basename(dir1.rstrip("/"))}{x}'
                    if idx < len(data_list):
                        dmr_record[col_name] = data_list[idx].qvalue
                    else:
                        dmr_record[col_name] = 1.0
                    idx += 1
            final_dmr_list.append(dmr_record)

    if not final_dmr_list:
        print(f"  {methylation_type} 无显著 DMR")
        return

    # 转换为 DataFrame 并排序
    final_df = pd.DataFrame(final_dmr_list)
    final_df = final_df.sort_values(['Chromosome', 'DMR_start'])
    column_order = [
        'Chromosome', 'Methylation_Type', 'DMR_start', 'DMR_end', 'Length',
        'Direction', 'Sig_count', 'Total_count', 'Sig_probability',
        'Avg_exp_methy', 'Avg_exp_unmethy', 'Avg_wild_methy', 'Avg_wild_unmethy',
        'Sig_Avg_qvalue'
    ]
    replicate_columns = sorted(
        [col for col in final_df.columns if col.startswith('qvalue_')],
        key=lambda x: tuple(map(int, re.findall(r'\d+', x)))
    )
    column_order = column_order + replicate_columns
    final_df = final_df[column_order]
    # 保存结果
    output_file = os.path.join(output_dir, f"{methylation_type}-final_significant_regions_DMRs.txt")
    final_df.to_csv(output_file, sep='\t', index=False)

    print(f"  {methylation_type} 最终显著 DMR: {len(final_df)} 个")
    print(f"  保存至: {output_file}")

    # 统计信息
    hyper_count = len(final_df[final_df['Direction'] == 1])
    hypo_count = len(final_df[final_df['Direction'] == 0])
    print(f"    - 高甲基化: {hyper_count} ({hyper_count / len(final_df) * 100:.1f}%)")
    print(f"    - 低甲基化: {hypo_count} ({hypo_count / len(final_df) * 100:.1f}%)")

    return final_df

def collect_dmr_results(methy_dir, methylation_type, all_dmr_results):
    """
    收集单次 DMR 分析的结果

    参数：
        methy_dir: 甲基化类型目录（如 ./and_output/CpG/）
        methylation_type: 甲基化类型
        all_dmr_results: 汇总字典，输入的时候只有三个甲基化类型键，后面内容是空的
    """
    # 查找所有 dmr_fisher_Chr*.txt 文件
    fisher_files = glob.glob(os.path.join(methy_dir, "dmr_fisher_Chr*.txt"))

    for fisher_file in fisher_files:
        # 提取染色体号
        match = re.search(r'Chr(\d+)\.txt$', fisher_file)
        if not match:
            continue
        chr_num = int(match.group(1))

        # 读取文件
        try:
            df = pd.read_csv(fisher_file, sep='\t')

            if df.empty:
                continue

            # 初始化该染色体的字典
            if chr_num not in all_dmr_results[methylation_type]:
                all_dmr_results[methylation_type][chr_num] = defaultdict(list)

            # 收集每个 DMR 区域的数据
            for _, row in df.iterrows():
                dmr_key = (int(row['DMR_start']), int(row['DMR_end']))

                # 存储 (exp_m, exp_u, wild_m, wild_u, qvalue, direction)
                dmr_data = DmrRecord(
                    exp_methy=int(row['exp_methy_sum']),
                    exp_unmethy=int(row['exp_unmethy_sum']),
                    wild_methy=int(row['wild_methy_sum']),
                    wild_unmethy=int(row['wild_unmethy_sum']),
                    qvalue=float(row['qvalue']) if not pd.isna(row['qvalue']) else 1.0,
                    direction=int(row['direction'])
                )

                all_dmr_results[methylation_type][chr_num][dmr_key].append(dmr_data)

        except Exception as e:
            print(f"    警告：读取 {fisher_file} 失败: {e}")
            continue

def summarize_all_dmr_methylation(dir1, dir2, m, n, methylation_types=['CpG', 'CHH', 'CHG'], work_dir=".", threads=1):
    """
    对所有 output_x_y 目录，使用 common DMR 区域进行甲基化读段求和分析。
    并行版中，worker 只写自己的 dmr_analysis_wt*_mut*/context 目录；
    汇总 all_dmr_results 的步骤仍由主进程顺序完成，避免共享字典并发写入。
    """
    print("第四阶段：Common DMR 在各组合中的甲基化读段求和分析")

    and_output_dir = os.path.join(work_dir, "and_output")
    os.makedirs(and_output_dir, exist_ok=True)

    all_dmr_results = {mtype: {} for mtype in methylation_types}
    threads = max(1, int(threads))
    config = _build_parallel_config()

    tasks = []
    for replicate_x in range(1, m + 1):
        for replicate_y in range(1, n + 1):
            file1_path_template = os.path.join(dir1, f"{replicate_x}-bothMeUnme_diffChromo_NOREPEATED_methy_sites_{{}}.txt")
            file2_path_template = os.path.join(dir2, f"{replicate_y}-bothMeUnme_diffChromo_NOREPEATED_methy_sites_{{}}.txt")

            for mtype in methylation_types:
                f1 = file1_path_template.format(mtype)
                f2 = file2_path_template.format(mtype)
                if not os.path.exists(f1) or not os.path.exists(f2):
                    print(f"  跳过 wt{replicate_y}_mut{replicate_x} {mtype}：both 文件不存在")
                    continue

                methy_output_dir = os.path.join(and_output_dir, f"dmr_analysis_wt{replicate_y}_mut{replicate_x}", mtype)
                os.makedirs(methy_output_dir, exist_ok=True)

                tasks.append({
                    "replicate_x": replicate_x,
                    "replicate_y": replicate_y,
                    "methylation_type": mtype,
                    "file1_path": f1,
                    "file2_path": f2,
                    "methy_output_dir": methy_output_dir,
                    "custom_dmr_dir": and_output_dir,
                    "config": config,
                })

    completed = []

    if threads <= 1 or len(tasks) <= 1:
        for idx, task in enumerate(tasks, 1):
            print(f"\n处理组合 ({dir1}{task['replicate_x']}, {dir2}{task['replicate_y']}) {task['methylation_type']}...")

            try:
                summarize_dmr_methylation(
                    methy_dir=task["methy_output_dir"],
                    replicate_x=task["replicate_x"],
                    replicate_y=task["replicate_y"],
                    file1_path=task["file1_path"],
                    file2_path=task["file2_path"],
                    methylation_type=task["methylation_type"],
                    custom_dmr_dir=and_output_dir,
                )
                completed.append(task)
            except Exception as e:
                print(f"  错误：处理 {task['methylation_type']} 失败: {e}")
                import traceback
                traceback.print_exc()
                continue
    else:
        max_workers = min(threads, len(tasks))
        log_dir = os.path.join(work_dir, "parallel_logs", "common_dmr_summarize")
        os.makedirs(log_dir, exist_ok=True)
        for task in tasks:
            task["log_file"] = os.path.join(
                log_dir,
                f"dmr_summary_wt{task['replicate_y']}_mut{task['replicate_x']}_{task['methylation_type']}.log"
            )

        print(f"启用并行 common DMR reads/Fisher/FDR 处理: workers={max_workers}, tasks={len(tasks)}")
        with ProcessPoolExecutor(max_workers=max_workers) as executor:
            future_to_task = {executor.submit(_run_summarize_dmr_worker, task): task for task in tasks}
            for idx, future in enumerate(as_completed(future_to_task), 1):
                task = future_to_task[future]
                label = f"wt{task['replicate_y']}_mut{task['replicate_x']} {task['methylation_type']}"
                try:
                    res = future.result()
                    completed.append(task)
                    print(f"[DONE] {idx}/{len(tasks)} {label}: {res['elapsed']:.2f}s, log={res.get('log_file')}")
                except Exception as e:
                    print(f"[FAILED] {label}: {e}")
                    raise

    # 主进程统一收集结果，避免并发修改 all_dmr_results。
    # 注意：并行任务完成顺序是不确定的；如果直接按完成顺序 collect，
    # data_list 中的 qvalue 会被追加成随机顺序，导致最终 DMR 文件的 qvalue_* 列错位。
    # 因此这里必须恢复原串行版顺序：replicate_x -> replicate_y -> methylation_type。
    completed = sorted(
        completed,
        key=lambda t: (int(t["replicate_x"]), int(t["replicate_y"]), str(t["methylation_type"]))
    )

    print("收集 common DMR 在各组合中的 Fisher/FDR 结果...")
    for task in completed:
        collect_dmr_results(task["methy_output_dir"], task["methylation_type"], all_dmr_results)

    if AUTO_DMR_VOTE_THRESHOLD or AUTO_VOTE_THRESHOLD_REPORT_ONLY:
        update_auto_dmr_vote_thresholds(m=m, n=n, work_dir=work_dir)

    print("第五阶段：汇总 DMR 结果并进行贝叶斯判定")

    for mtype in methylation_types:
        generate_final_significant_dmr(
            all_dmr_results[mtype],
            methylation_type=mtype,
            output_dir=and_output_dir,
            m=m,
            n=n,
            dir1=dir1,
            dir2=dir2
        )

def process_common_sites_to_dmr(methylation_type='CpG',work_dir="."):
    """
    读取 final_significant_sites_DMPs.txt 文件，按染色体分组处理DMR

    参数：

    methylation_type: 甲基化类型（CpG, CHH, CHG）
    返回：
        dmr_results: 字典 {chr_num: dmr_list_file_path}
    """
    print(f"\n开始处理 {methylation_type} 的共同显著位点 DMR 分析...")

    and_output_dir = os.path.join(work_dir, "and_output")
    os.makedirs(and_output_dir, exist_ok=True)

    # 1. 读取 final_significant_sites_DMPs 文件
    common_file = os.path.join(and_output_dir, f"{methylation_type}-final_significant_sites_DMPs.txt")
    if not os.path.exists(common_file):
        print(f"错误：文件不存在 {common_file}")
        return None

    try:
        df = pd.read_csv(common_file, sep='\t')
        print(f"成功读取文件，共 {len(df)} 个位点")
    except Exception as e:
        print(f"读取文件失败: {e}")
        return None

    if df.empty:
        print(f"警告：{methylation_type} 文件为空")
        return None

    # 2. 按染色体分组
    chromosomes = sorted(df['Chromosome'].unique(), key=natural_sort_key)
    print(f"找到 {len(chromosomes)} 个染色体: {chromosomes}")

    dmr_results = {}
    total_chromosomes = len(chromosomes)

    # 3. 对每个染色体进行处理
    for chr_num in chromosomes:
        print(f"\n  处理染色体 {chr_num} ({chr_num}/{total_chromosomes})...")

        # 筛选当前染色体的数据
        chr_df = df[df['Chromosome'] == chr_num].copy()

        # 提取需要的三列：Position, Sig_Mean_Qvalue, Methylation_Change
        dmp_data = chr_df[['Position', 'Sig_Mean_Qvalue', 'Methylation_Change']].copy()

        # 确保数据类型正确
        dmp_data['Position'] = dmp_data['Position'].astype(int)
        dmp_data['Sig_Mean_Qvalue'] = dmp_data['Sig_Mean_Qvalue'].astype(float)
        dmp_data['Methylation_Change'] = dmp_data['Methylation_Change'].astype(int)

        # 按位点号排序
        dmp_data = dmp_data.sort_values('Position').reset_index(drop=True)

        print(f"    染色体 {chr_num} 共 {len(dmp_data)} 个位点")

        if len(dmp_data) == 0:
            print(f"    跳过：染色体 {chr_num} 无有效位点")
            continue

        # 4. 创建临时 DMP 文件（格式：pos qvalue change）
        temp_dmp_file = os.path.join(and_output_dir,
                                     f"DMP_common_{methylation_type}_Chr{chr_num}.txt")

        # 写入 DMP 格式文件（第一行是 "first line"，后续是 pos qvalue change）
        with open(temp_dmp_file, 'w') as f:
            f.write("first line\n")
            for _, row in dmp_data.iterrows():
                f.write(f"{int(float(row['Position']))} {float(row['Sig_Mean_Qvalue'])} {int(float(row['Methylation_Change']))}\n")

        print(f"    创建 DMP 文件: {os.path.basename(temp_dmp_file)}")

        # 5. 调用 DMR 分析函数
        # 注意：run_dmr_pipeline_on_dmp_file 需要 chromoNo 参数
        # 我们传入总染色体数
        try:
            dmr_list_file = run_dmr_pipeline_on_dmp_file_auto(
                dmp_file=temp_dmp_file,
                chromoNo=total_chromosomes
            )

            if dmr_list_file:
                dmr_results[chr_num] = dmr_list_file
                print(f"     染色体 {chr_num} DMR 分析完成")
            else:
                print(f"    染色体 {chr_num} DMR 分析失败（可能无有效DMR）")

        except Exception as e:
            print(f"     染色体 {chr_num} 处理出错: {e}")
            import traceback
            traceback.print_exc()
            continue

    print(f"\n{methylation_type} 共同显著位点 DMR 分析完成！")
    print(f"成功处理 {len(dmr_results)}/{len(chromosomes)} 个染色体")

    return dmr_results

def summarize_dmr_methylation(methy_dir, replicate_x, replicate_y, file1_path, file2_path, methylation_type='CpG', custom_dmr_dir=None):
    """
    对 DMR 区域进行甲基化读段求和，并计算 Fisher p 值、FDR q 值和甲基化方向。
    所有输出文件将保存在 methy_dir 目录下（例如：./output_1_1/CpG/）。
    注意！！！custom_dmr_dir: 自定义 DMR 文件所在目录（如果为 None，则使用 methy_dir），这样能区分每个output_x_y自己生成的dmr和common_sites生成的dmr
    """
    # 这里是只处理一次检验的
    print(f"    开始 DMR 甲基化读段求和分析...")

    n_chromosomes = get_column_count(file1_path)
    if n_chromosomes is None:
        print("    无法获取染色体数量，跳过 DMR 求和")
        return

    chromosomes = [f'Chr{i}' for i in range(1, n_chromosomes + 1)]

    # 收集所有染色体的 DMR 数据（不含 p 值）
    all_dmr_data = []  # (chrom, start, end, exp_m, exp_u, wild_m, wild_u)

    for idx, chrom in enumerate(chromosomes):
        chrom_num = idx + 1
        if custom_dmr_dir is not None:
            # 读取 common DMR 文件
            dmr_file = os.path.join(custom_dmr_dir, f"DMR_list_DMP_common_{methylation_type}_Chr{chrom_num}.txt")
        else:
            # 读取 output目录 的 DMR 文件
            dmr_file = os.path.join(methy_dir, f"DMR_list_DMP_wt_replicate{replicate_y}_mut_replicate{replicate_x}_Chr{chrom_num}.txt")
        if not os.path.exists(dmr_file):
            print(f"      跳过 {chrom}，DMR 文件不存在")
            continue

        # 读取 DMR 区域
        dmr_list = []
        with open(dmr_file, 'r') as f:
            for line in f:
                parts = line.strip().split()
                if len(parts) >= 2:
                    try:
                        start = int(parts[0])
                        end = int(parts[1])
                        dmr_list.append((start, end))
                    except ValueError:
                        continue

        if not dmr_list:
            print(f"      {chrom} 无有效 DMR 区域")
            continue

        col_start = idx * 3

        # 读取实验组数据
        exp_sites, exp_methy_dict, exp_unmethy_dict = [], {}, {}
        try:
            with open(file1_path, 'r') as f:
                for line in f:
                    parts = line.strip().split()
                    if len(parts) < col_start + 3:
                        continue
                    try:
                        site = int(float(parts[col_start]))
                        methy = int(float(parts[col_start + 1]))
                        unmethy = int(float(parts[col_start + 2]))
                        if site > 0:
                            exp_sites.append(site)
                            exp_methy_dict[site] = methy
                            exp_unmethy_dict[site] = unmethy
                    except:
                        continue
            exp_sites.sort()
        except Exception as e:
            print(f"      读取实验组失败 ({chrom}): {e}")
            continue

        # 读取对照组数据
        wild_sites, wild_methy_dict, wild_unmethy_dict = [], {}, {}
        try:
            with open(file2_path, 'r') as f:
                for line in f:
                    parts = line.strip().split()
                    if len(parts) < col_start + 3:
                        continue
                    try:
                        site = int(float(parts[col_start]))
                        methy = int(float(parts[col_start + 1]))
                        unmethy = int(float(parts[col_start + 2]))
                        if site > 0:
                            wild_sites.append(site)
                            wild_methy_dict[site] = methy
                            wild_unmethy_dict[site] = unmethy
                    except:
                        continue
            wild_sites.sort()
        except Exception as e:
            print(f"      读取对照组失败 ({chrom}): {e}")
            continue

        # 汇总每个 DMR 区域的 reads
        for start, end in dmr_list:
            l1 = bisect.bisect_left(exp_sites, start)
            r1 = bisect.bisect_right(exp_sites, end)
            exp_m_sum = sum(exp_methy_dict[site] for site in exp_sites[l1:r1])
            exp_u_sum = sum(exp_unmethy_dict[site] for site in exp_sites[l1:r1])

            l2 = bisect.bisect_left(wild_sites, start)
            r2 = bisect.bisect_right(wild_sites, end)
            wild_m_sum = sum(wild_methy_dict[site] for site in wild_sites[l2:r2])
            wild_u_sum = sum(wild_unmethy_dict[site] for site in wild_sites[l2:r2])

            all_dmr_data.append((chrom, start, end, exp_m_sum, exp_u_sum, wild_m_sum, wild_u_sum))

    if not all_dmr_data:
        print("    无任何有效 DMR 数据，跳过后续分析")
        return

    # === 第一步：输出 dmr_summary_{chrom}.txt ===
    chrom_summary_dict = defaultdict(list)
    for item in all_dmr_data:
        chrom, start, end, exp_m, exp_u, wild_m, wild_u = item
        chrom_summary_dict[chrom].append((start, end, exp_m, exp_u, wild_m, wild_u))

    for chrom in chromosomes:
        if chrom in chrom_summary_dict:
            summary_file = os.path.join(methy_dir, f"dmr_summary_{chrom}.txt")
            with open(summary_file, 'w') as f:
                f.write("DMR_start\tDMR_end\texp_methy_sum\texp_unmethy_sum\twild_methy_sum\twild_unmethy_sum\n")
                for row in chrom_summary_dict[chrom]:
                    f.write(f"{row[0]}\t{row[1]}\t{row[2]}\t{row[3]}\t{row[4]}\t{row[5]}\n")
            print(f"      {chrom} DMR 汇总完成 → {summary_file}")

    # === 第二步：为每个 DMR 计算 p 值 ===
    dmr_with_pvals = []
    for (chrom, start, end, exp_m, exp_u, wild_m, wild_u) in all_dmr_data:
        if (exp_m + exp_u < 4) or (wild_m + wild_u < 4):
            pval = np.nan
        else:
            try:
                _, pval = fisher_exact([[exp_m, exp_u], [wild_m, wild_u]], alternative='two-sided')
                pval = float(pval)
            except:
                pval = np.nan
        dmr_with_pvals.append((chrom, start, end, exp_m, exp_u, wild_m, wild_u, pval))

    # === 第三步：全局 FDR 校正（跨染色体） ===
    pvalues = np.array([item[7] for item in dmr_with_pvals])
    qvalues = calculate_qvalues(pvalues, pi=1.0)

    # === 第四步：按染色体组织结果，添加 direction ===
    chrom_data_dict = defaultdict(list)
    for i, (chrom, start, end, exp_m, exp_u, wild_m, wild_u, pval) in enumerate(dmr_with_pvals):
        qval = qvalues[i]

        # 计算甲基化方向
        exp_total = exp_m + exp_u
        wild_total = wild_m + wild_u
        if exp_total > 0 and wild_total > 0:
            exp_rate = exp_m / exp_total
            wild_rate = wild_m / wild_total
            direction = 1 if exp_rate > wild_rate else 0
        else:
            direction = 0  # 无法判断时设为 0

        chrom_data_dict[chrom].append((start, end, exp_m, exp_u, wild_m, wild_u, pval, qval, direction))

    # === 第五步：输出完整 Fisher 结果 + 显著子集 ===
    for chrom in chromosomes:
        if chrom not in chrom_data_dict:
            continue

        # 完整结果
        fisher_file = os.path.join(methy_dir, f"dmr_fisher_{chrom}.txt")
        with open(fisher_file, 'w') as f:
            f.write("DMR_start\tDMR_end\texp_methy_sum\texp_unmethy_sum\twild_methy_sum\twild_unmethy_sum\tpvalue\tqvalue\tdirection\n")
            for row in chrom_data_dict[chrom]:
                start, end, exp_m, exp_u, wild_m, wild_u, pval, qval, direction = row
                p_out = pval if not np.isnan(pval) else 'nan'
                q_out = qval if not np.isnan(qval) else 'nan'
                f.write(f"{start}\t{end}\t{exp_m}\t{exp_u}\t{wild_m}\t{wild_u}\t{p_out:.6g}\t{q_out:.6g}\t{direction}\n")
        print(f"      {chrom} Fisher + FDR + direction 完成 → {fisher_file}")

        # 显著结果（使用可配置的 DMR q-value 阈值）
        sig_rows = [
            row for row in chrom_data_dict[chrom]
            if not np.isnan(row[7]) and row[7] <= DMR_QVALUE_THRESHOLD
        ]
        if sig_rows:
            sig_file = os.path.join(methy_dir, f"dmr_fisher_significant_{chrom}.txt")
            with open(sig_file, 'w') as f:
                f.write("DMR_start\tDMR_end\texp_methy_sum\texp_unmethy_sum\twild_methy_sum\twild_unmethy_sum\tpvalue\tqvalue\tdirection\n")
                for row in sig_rows:
                    start, end, exp_m, exp_u, wild_m, wild_u, pval, qval, direction = row
                    f.write(f"{start}\t{end}\t{exp_m}\t{exp_u}\t{wild_m}\t{wild_u}\t{pval:.6g}\t{qval:.6g}\t{direction}\n")
            print(f"        → 显著 DMR (q<={DMR_QVALUE_THRESHOLD}): {sig_file}")
        else:
            print(f"        → {chrom} 无显著 DMR (q<={DMR_QVALUE_THRESHOLD})")

    print("    DMR 甲基化分析完成。")

def run_dmr_pipeline_on_dmp_file(dmp_file: str, chromoNo: int = 10):
    """
    从 DMP 文件生成 DMR 区域 , 生成dmr_list文件
    - methylation_matrix_file: 用于动态获取染色体数量（如 1-bothMeUnme_...txt）
    - 所有输出文件保存在 dmp_file 所在目录。
    """
    # 总的核心思想如下：
    # 用滑动窗口找到 DMP 密集的区域
    # 通过"跳跃合并"连接相邻的密集区域
    # 最终筛选出足够长、足够显著的 DMR
    sWinN = 1000  # 滑动窗口大小（1000 bp）
    M0 = 4  # 窗口内最少 DMP 数量（至少 4 个）
    M1 = 10  # 最终 DMR 内最少 DMP 数量（至少 10 个）
    M2 = 10  # 跳跃步长，保留原始分割行为，

    # 为安全起见，确保 chromoNo >= 6（因为用到 arrayMethy1_script1[5]）
    chromoNo = max(chromoNo, 6)

    class PositionNoNode:
        def __init__(self, pos=0, end=0, pV=0.0, ratio=0.0, num=0, num2=0, numCom=0, markR=0, DMR_S=0, DMR_E=0):
            self.pos = pos  # 窗口起始位置
            self.end = end  # 窗口结束位置
            self.pV = pV
            self.ratio = ratio
            self.num = num  # 高甲基化数量
            self.num2 = num2  # 低甲基化数量
            self.numCom = numCom
            self.markR = markR
            self.DMR_S = DMR_S
            self.DMR_E = DMR_E
            self.posV = []
            self.logPV = []
            self.meUnV = []

    output_dir = os.path.dirname(dmp_file)
    base_name = os.path.basename(dmp_file)

    arrayMethy1 = [[] for _ in range(chromoNo)]  # 创建染色体个数个子列表

    # === 读取 DMP 位点 ===
    try:
        with open(dmp_file, 'r') as fin1:
            lines = fin1.readlines()
            if not lines or len(lines) < 2:
                return None
            for line in lines[1:]:
                parts = line.strip().split()
                if len(parts) < 3:
                    continue
                try:
                    pos = int(parts[0])
                    pValue = float(parts[1])
                    change = int(parts[2])
                    tmpNode = PositionNoNode()
                    tmpNode.pos = pos
                    tmpNode.pV = pValue
                    tmpNode.num = change  # 以上几行建立了当前dmp行所对应的位点信息
                    arrayMethy1[0].append(tmpNode)
                except (ValueError, IndexError):
                    continue
    except Exception as e:
        print(f"    DMP 文件读取失败 {dmp_file}: {e}")
        return None

    if not arrayMethy1[0]:
        return None

    arrayMethy1[0].sort(key=lambda x: x.pos)
    firstP = arrayMethy1[0][0].pos
    lastP = arrayMethy1[0][-1].pos

    # === 滑动窗口 ===
    ratioI = 0.1
    cirN0 = (firstP - 1) // sWinN
    cir2 = 0
    while True:
        startP = cirN0 * sWinN + int(sWinN * ratioI * cir2) + 1
        endP = cirN0 * sWinN + int(sWinN * (ratioI * cir2 + 1))
        if startP > lastP:
            break

        numS1 = numS2 = 0
        for node in arrayMethy1[0]:
            if node.pos > endP:
                break
            if startP <= node.pos <= endP:
                if node.num == 1:
                    numS1 += 1
                else:
                    numS2 += 1

        tmpNode = PositionNoNode()
        tmpNode.pos = startP
        tmpNode.end = endP
        tmpNode.num = numS1
        tmpNode.num2 = numS2
        arrayMethy1[1].append(tmpNode)
        cir2 += 1

    # === 构建滑动窗口列表 ===
    arrayMethy1[2] = []
    for node in arrayMethy1[1]:
        new_node = PositionNoNode()
        new_node.pos = node.pos
        new_node.end = node.end
        new_node.num = node.num
        new_node.num2 = node.num2
        new_node.numCom = node.num + node.num2
        new_node.markR = 0
        arrayMethy1[2].append(new_node)

    # === 标准化输出 ===
    if arrayMethy1[2]:
        maxCom = max(node.numCom for node in arrayMethy1[2])
        maxCom = max(maxCom, 1)
        std_file = os.path.join(output_dir, f"noTitle_allDMCs_new_Standardized_slidingW_{base_name}")
        with open(std_file, 'w') as f:
            for node in arrayMethy1[2]:
                if node.end <= lastP:
                    std_val = node.numCom / maxCom
                    f.write(f"{node.pos:<20}{node.end:<20}{node.numCom:<20}{std_val:.6f}\n")

    # === 跳跃合并（使用你提供的代码）===
    dmr_out_file = os.path.join(output_dir, f"DMR_{base_name}")
    boundary_file = os.path.join(output_dir, f"DMR_boundaries_{base_name}")

    arrayMethy1[3] = []
    for i, it02 in enumerate(arrayMethy1[2]):
        it03 = it02
        it04 = it02
        num1 = it02.pos
        num2 = it02.end
        num3 = it02.numCom

        num01 = num1
        num02 = num2
        num03 = num3
        num04 = 0

        flag1 = False
        flag2 = False
        flag3 = False
        flag4 = False
        flag5 = False
        flag6 = False
        flag7 = False

        # 向左扩展
        while num03 >= M0 and num01 >= firstP:
            if not flag3:
                flag1 = True
                flag3 = True
                num07 = num02
                num09 = num01
                num10 = num02

            num04 += num03
            num06 = num01
            it03.markR = 1

            idx = arrayMethy1[2].index(it03)
            if idx > 0:
                for k in range(M2):
                    idx -= 1
                    if idx <= 0:
                        flag5 = True
                        break
                    it03 = arrayMethy1[2][idx]
                if flag5:
                    break
                else:
                    num01 = it03.pos
                    num03 = it03.numCom
            else:
                break

        # 向右扩展
        if flag1:
            idx = arrayMethy1[2].index(it04)
            for k in range(M2):
                idx += 1
                if idx >= len(arrayMethy1[2]) - 1:
                    flag6 = True
                    break
                it04 = arrayMethy1[2][idx]

            if not flag6:
                while idx < len(arrayMethy1[2]):
                    it04 = arrayMethy1[2][idx]
                    num02 = it04.end
                    num03 = it04.numCom
                    # M0:4
                    if num03 >= M0 and num02 <= lastP:
                        if not flag4:
                            flag2 = True
                            flag4 = True

                        num07 = num02
                        num04 += num03
                        it04.markR = 1

                        for k in range(M2):
                            idx += 1
                            if idx >= len(arrayMethy1[2]) - 1:
                                flag7 = True
                                break
                            it04 = arrayMethy1[2][idx]

                        if flag7:
                            break
                    else:
                        break

        if flag1 or flag2:
            tmpNode = PositionNoNode()
            tmpNode.pos = num06
            tmpNode.end = num07
            tmpNode.numCom = num04
            tmpNode.DMR_S = num09
            tmpNode.DMR_E = num10
            arrayMethy1[3].append(tmpNode)

    print("=" * 60)
    print(f"识别到 {len(arrayMethy1[3])} 个DMR区域")

    # 输出DMR结果和边界文件
    with open(dmr_out_file, 'w') as cout05:
        with open(boundary_file, 'w') as bound_out:
            for node in arrayMethy1[3]:
                cout05.write(f"{node.pos:<20}{node.end:<20}{node.numCom:<20}[{node.DMR_S} {node.DMR_E}]\n")
                bound_out.write(f"{node.DMR_S} {node.DMR_E}\n")

    print(f"已生成边界文件: {boundary_file}")

    # === 合并重叠边界（动态 chromoL）===
    chromoL = lastP + 100000  # 动态基因组长度
    chromoArray = [0] * chromoL
    try:
        with open(boundary_file, 'r') as f_in:
            for line in f_in:
                parts = line.strip().split()
                if len(parts) >= 2:
                    try:
                        s, e = int(parts[0]), int(parts[1])
                        for idx in range(s, min(e + 1, chromoL)):
                            chromoArray[idx] = 1
                    except:
                        continue
    except:
        return boundary_file

    arrayMethy1_script1 = [[] for _ in range(chromoNo)]
    i = 0
    while i < chromoL:
        start = None
        for j in range(i, chromoL):
            if chromoArray[j] == 1:
                start = j
                break
        if start is None:
            break
        end = start
        while end + 1 < chromoL and chromoArray[end + 1] == 1:
            end += 1
        node = PositionNoNode()
        node.DMR_S = start
        node.DMR_E = end
        arrayMethy1_script1[0].append(node)
        i = end + 1

    no_overlap_file = os.path.join(output_dir, f"boundaries_noOverlapping_{base_name}")
    with open(no_overlap_file, 'w') as f_no:
        for node in arrayMethy1_script1[0]:
            f_no.write(f"{node.DMR_S}    {node.DMR_E}\n")

    arrayMethy1_script1[5] = arrayMethy1[0]
    final_dmr_list = []
    for region in arrayMethy1_script1[0]:
        s, e = region.DMR_S, region.DMR_E
        hyper = hypo = 0
        for site in arrayMethy1_script1[5]:
            if site.pos < s:
                continue
            if site.pos > e:
                break
            if site.num == 1:
                hyper += 1
            else:
                hypo += 1

        total = hyper + hypo
        length = e - s + 1
        if total >= M1 and length >= sWinN:
            node = PositionNoNode()
            node.pos = s
            node.end = e
            node.num = hyper
            node.num2 = hypo
            node.numCom = total
            final_dmr_list.append(node)

    final_file = os.path.join(output_dir, f"DMR_list_{base_name}")
    final_dmr_list.sort(key=lambda x: x.pos)
    with open(final_file, 'w') as f_final:
        for node in final_dmr_list:
            f_final.write(f"{node.pos:<20}{node.end:<20}{node.num:<20}{node.num2:<20}{node.numCom:<20}\n")

    print(f"    DMR 分析完成: {base_name} → {final_file}")
    return final_file


def resolve_cpp_dmr_binary(binary_name: str) -> str:
    """
    Locate a bundled C++ DMR executable.

    Search order:
    1. directory containing this Python script
    2. current working directory
    3. PATH
    """
    script_dir = os.path.dirname(os.path.abspath(__file__))
    candidates = [
        os.path.join(script_dir, binary_name),
        os.path.join(os.getcwd(), binary_name),
    ]

    for candidate in candidates:
        if os.path.isfile(candidate) and os.access(candidate, os.X_OK):
            return os.path.abspath(candidate)

    found = shutil.which(binary_name)
    if found:
        return os.path.abspath(found)

    raise FileNotFoundError(
        f"Cannot find executable '{binary_name}'. "
        f"Please put it beside {os.path.basename(__file__)} or in the current directory."
    )


def _copy_file_if_needed(src: str, dst: str):
    """Copy src to dst unless they are already the same absolute path."""
    if os.path.abspath(src) != os.path.abspath(dst):
        shutil.copyfile(src, dst)


def run_dmr_pipeline_on_dmp_file_cpp(dmp_file: str, chromoNo: int = 10):
    """
    C++ replacement for DMP -> DMR candidate/list generation.

    It calls:
      1. dmr_step1 <DMP_file>
      2. dmr_step2_dynamic <DMR_boundaries_file> <DMP_file>

    The C++ step2 naturally writes names such as:
      DMR_list_DMR_boundaries_<base_name>

    This wrapper adapts them back to the Python-expected names:
      DMR_list_<base_name>
      boundaries_noOverlapping_<base_name>
    """
    dmp_abs = os.path.abspath(dmp_file)
    output_dir = os.path.dirname(dmp_abs) or "."
    base_name = os.path.basename(dmp_abs)

    if not os.path.exists(dmp_abs):
        print(f"    C++ DMR 输入文件不存在: {dmp_abs}")
        return None

    # Match the Python behavior for empty / header-only DMP files.
    try:
        with open(dmp_abs, "r") as f:
            lines = f.readlines()
        valid_data_lines = [line for line in lines[1:] if line.strip()]
        if len(lines) < 2 or not valid_data_lines:
            return None
    except Exception as e:
        print(f"    C++ DMR 读取 DMP 文件失败 {dmp_abs}: {e}")
        return None

    try:
        step1_bin = resolve_cpp_dmr_binary("dmr_step1")
        step2_bin = resolve_cpp_dmr_binary("dmr_step2_dynamic")
    except Exception as e:
        print(f"    C++ DMR executable 查找失败: {e}")
        raise

    # dmr_step1 historically writes DMR_<base_name>, not DMR_boundaries_<base_name>.
    # Some previous bridge tests also used DMR_boundaries_<base_name>, so we support both,
    # but prefer the true step1 output DMR_<base_name>.
    boundary_candidates = [
        f"DMR_{base_name}",
        f"DMR_boundaries_{base_name}",
        f"DMR_noOverlap_{base_name}",
        f"DMR_questionNoOverlap_{base_name}",
    ]

    boundary_base = None
    boundary_file = None
    for cand in boundary_candidates:
        cand_path = os.path.join(output_dir, cand)
        if os.path.exists(cand_path):
            boundary_base = cand
            boundary_file = cand_path
            break

    # The raw step2 outputs are named from argv[1], i.e. from boundary_base.
    raw_noover_base = None
    raw_list_base = None
    raw_noover_file = None
    raw_list_file = None

    expected_noover_file = os.path.join(output_dir, f"boundaries_noOverlapping_{base_name}")
    expected_list_file = os.path.join(output_dir, f"DMR_list_{base_name}")

    log1 = os.path.join(output_dir, f"cpp_dmr_step1_{base_name}.log")
    log2 = os.path.join(output_dir, f"cpp_dmr_step2_{base_name}.log")

    print(f"    使用 C++ DMR engine 处理: {base_name}")

    # Run from output_dir so that C++ output filenames remain clean and match the original naming convention.
    with open(log1, "w") as lf1:
        subprocess.run(
            [step1_bin, base_name],
            cwd=output_dir,
            stdout=lf1,
            stderr=subprocess.STDOUT,
            check=True
        )

    # Re-scan after step1, because these files are created by dmr_step1.
    boundary_base = None
    boundary_file = None
    for cand in boundary_candidates:
        cand_path = os.path.join(output_dir, cand)
        if os.path.exists(cand_path):
            boundary_base = cand
            boundary_file = cand_path
            break

    if boundary_file is None:
        existing = sorted(
            f for f in os.listdir(output_dir)
            if ("DMR" in f or "boundar" in f or "Boundary" in f)
        )
        raise RuntimeError(
            "C++ step1 did not generate a usable boundary/input file. "
            f"Checked candidates: {boundary_candidates}. Existing DMR-like files: {existing}"
        )

    # Build a clean two-column boundary file for step2.
    # dmr_step1 writes DMR_<base_name> with columns like:
    #   expanded_start expanded_end numCom [DMR_S DMR_E]
    # The original Python pipeline feeds DMR_S/DMR_E boundaries to step2, not necessarily
    # the expanded_start/expanded_end columns. Extracting the bracket avoids 1-bp
    # off-by-one differences such as 1512299 vs 1512300.
    clean_boundary_base = f"DMR_boundaries_{base_name}"
    clean_boundary_file = os.path.join(output_dir, clean_boundary_base)

    extracted = 0
    with open(boundary_file, "r") as src, open(clean_boundary_file, "w") as dst:
        for line in src:
            line = line.strip()
            if not line:
                continue

            m = re.search(r"\[\s*(\d+)\s+(\d+)\s*\]", line)
            if m:
                dst.write(f"{m.group(1)} {m.group(2)}\n")
                extracted += 1
                continue

            # Fallback for already-clean two-column files.
            parts = line.split()
            if len(parts) >= 2 and parts[0].isdigit() and parts[1].isdigit():
                dst.write(f"{parts[0]} {parts[1]}\n")
                extracted += 1

    if extracted == 0:
        raise RuntimeError(f"Failed to extract clean boundaries from {boundary_file}")

    boundary_base = clean_boundary_base
    boundary_file = clean_boundary_file

    raw_noover_base = f"boundaries_noOverlapping_{boundary_base}"
    raw_list_base = f"DMR_list_{boundary_base}"
    raw_noover_file = os.path.join(output_dir, raw_noover_base)
    raw_list_file = os.path.join(output_dir, raw_list_base)

    with open(log2, "w") as lf2:
        subprocess.run(
            [step2_bin, boundary_base, base_name],
            cwd=output_dir,
            stdout=lf2,
            stderr=subprocess.STDOUT,
            check=True
        )

    if not os.path.exists(raw_list_file):
        existing = sorted(
            f for f in os.listdir(output_dir)
            if f.startswith("DMR_list") or f.startswith("boundaries_noOverlapping")
        )
        raise RuntimeError(
            f"C++ step2 did not generate expected DMR list file: {raw_list_file}. "
            f"Existing step2-like files: {existing}"
        )

    # Adapt C++ natural output names to names expected by summarize_dmr_methylation().
    _copy_file_if_needed(raw_list_file, expected_list_file)
    if os.path.exists(raw_noover_file):
        _copy_file_if_needed(raw_noover_file, expected_noover_file)

    print(f"    C++ DMR 分析完成: {base_name} → {expected_list_file}")
    return expected_list_file


def run_dmr_pipeline_on_dmp_file_auto(dmp_file: str, chromoNo: int = 10):
    """
    Dispatch DMP -> DMR candidate/list generation according to --dmr-engine.
    """
    if DMR_ENGINE == "cpp":
        return run_dmr_pipeline_on_dmp_file_cpp(dmp_file, chromoNo=chromoNo)
    return run_dmr_pipeline_on_dmp_file(dmp_file, chromoNo=chromoNo)


def process_chr_in_one_file(df):
    """将输入的单个文件的前缀改为chr并将该文件所有的染色体信息返回"""

    # 该函数使得以chr开头的染色体号仍以chr开头，若不以chr开头的会加上前缀chr作为新的染色体号
    def to_chr(chr_str):
        chr_str = str(chr_str).strip()
        if chr_str.lower().startswith('chr'):
            return f"chr{chr_str[3:]}"
        else:
            return f"chr{chr_str}"

    # 将所有染色体号改为chr开头
    df['染色体号'] = df['染色体号'].apply(to_chr)
    # 将所有染色体号的不重复集合返回
    return df['染色体号'].unique()  # 返回numpy.array

def natural_sort_key(chr_name):
    """
    自定义染色体排序函数
    数字染色体按数值排序，字母染色体(X,Y,M等)排在最后按字母排序
    """
    chr_name = str(chr_name).lower()
    # 去除 'chr' 前缀
    if chr_name.startswith('chr'):
        suffix = chr_name[3:]
    else:
        suffix = chr_name

    # 尝试转换为整数
    try:
        # 如果是数字，返回 (0, 数字值, '')
        return (0, int(suffix), '')
    except ValueError:
        # 如果是字母(如X,Y,M)，返回 (1, 0, 字母)
        return (1, 0, suffix)

def _normalize_chr_series(s):
    s = s.astype("string").str.strip()
    has_chr = s.str.lower().str.startswith("chr").fillna(False)
    return pd.Series(
        np.where(has_chr, "chr" + s.str.slice(3), "chr" + s),
        index=s.index
    )

def scan_all_files_for_chr_mapping_fast(m, n, dir1, dir2, chunksize=2_000_000):
    all_chromosomes = set()

    paths = []
    for i in range(1, m + 1):
        paths.append(os.path.join(dir1, f"{i}-{os.path.basename(dir1)}.txt"))
    for j in range(1, n + 1):
        paths.append(os.path.join(dir2, f"{j}-{os.path.basename(dir2)}.txt"))

    for filepath in paths:
        if not os.path.exists(filepath):
            print(f"警告: 文件不存在 {filepath}")
            continue

        try:
            for chunk in pd.read_csv(
                filepath,
                sep=r"\s+",
                header=None,
                usecols=[0],
                names=["染色体号"],
                dtype={"染色体号": "string"},
                chunksize=chunksize,
            ):
                chrs = _normalize_chr_series(chunk["染色体号"])
                all_chromosomes.update(chrs.dropna().unique().tolist())

            print(f"扫描染色体完成: {filepath}")

        except Exception as e:
            print(f"读取文件 {filepath} 时出错: {e}")

    unique_chrs = sorted(all_chromosomes, key=natural_sort_key)
    chr_series = pd.Series(range(len(unique_chrs)), index=unique_chrs)

    print(f"统一染色体映射: {chr_series}")
    return chr_series

def scan_all_files_for_chr_mapping(m, n, dir1, dir2):
    """扫描m+n个、即所有文件，利用上一个函数收集所有染色体信息，以生成统一映射"""

    all_chromosomes = set()  # 利用集合不重复的特性存储所有染色体号

    # 扫描第一个目录的m个文件，名称格式为 i-dir1.txt,如3-msv.txt
    for i in range(1, m + 1):
        filepath = os.path.join(dir1, f"{i}-{os.path.basename(dir1)}.txt")
        if not os.path.exists(filepath):
            print(f"警告: 文件不存在 {filepath}")
            continue

        try:
            df = pd.read_csv(filepath,
                             sep=r'\s+',
                             header=None,
                             names=['染色体号', '位点号', '甲基化读段数', '非甲基化读段数', '甲基化类型'],
                             dtype=str)
            chromosomes = process_chr_in_one_file(df)  # 获取到该文件中的所有染色体号的不重复集合
            all_chromosomes.update(chromosomes)  # 将该文件的所有不重复的染色体号放到all_chromosomes这个集合里
            # 注：update接受任何可迭代对象并将其中元素添加到集合中
            print(f"文件 {filepath} 包含染色体: {sorted(chromosomes, key=natural_sort_key)}")
        except Exception as e:
            print(f"读取文件 {filepath} 时出错: {e}")

    # 扫描第二个目录的文件
    for j in range(1, n + 1):
        filepath = os.path.join(dir2, f"{j}-{os.path.basename(dir2)}.txt")
        if not os.path.exists(filepath):
            print(f"警告: 文件不存在 {filepath}")
            continue

        try:
            df = pd.read_csv(filepath,
                             sep=r'\s+',
                             header=None,
                             names=['染色体号', '位点号', '甲基化读段数', '非甲基化读段数', '甲基化类型'],
                             dtype=str)
            chromosomes = process_chr_in_one_file(df)  # 同上
            all_chromosomes.update(chromosomes)  # 同上
            print(f"文件 {filepath} 包含染色体: {sorted(chromosomes, key=natural_sort_key)}")
        except Exception as e:
            print(f"读取文件 {filepath} 时出错: {e}")

    # 创建统一的染色体映射，
    # 注意：此时all_chromosomes这个集合中包含两个基因型文件夹中所有文件的不重复染色体号
    unique_chrs = sorted(all_chromosomes, key=natural_sort_key)
    # 建立一个两个基因型目录中所有 不重复的 且排序了的 染色体号->数值的Series，或者说数组
    chr_series = pd.Series(range(len(unique_chrs)), index=unique_chrs)

    print(f"统一染色体映射: {chr_series}")
    return chr_series

# def single_newtoboth(filepath1, output_dir, num1, chr_series):
    '''此处参数：filepath1即如1-wt.txt，是处理的新格式文件，
    output_dir是输出结果both文件到的目录，一般是filepath1所在的目录，
    num1是当前处理第几个新格式文件，即当前正处理num1-基因型.txt，
    chr_series是所有 不重复的 且排序了的 染色体号->数值的Series'''

    df = pd.read_csv(filepath1,
                     sep=r'\s+',
                     header=None,
                     names=['染色体号', '位点号', '甲基化读段数', '非甲基化读段数', '甲基化类型'],
                     dtype=str)

    def to_chr(chr_str):
        chr_str = str(chr_str).strip()
        if chr_str.lower().startswith('chr'):
            return f"chr{chr_str[3:]}"
        else:
            return f"chr{chr_str}"

    # 之前构建chr_series时没有修改源文件，所以此时读进来仍可能有不符合chr开头这个格式的染色体号，所以得处理
    df['染色体号'] = df['染色体号'].apply(to_chr)
    df['染色体号'] = df['染色体号'].map(chr_series)  # 将每个染色体号修改为其对应的数值，如chr1->0,chr2->1
    # 以下三行，将位点号、两个读段数转换为数值
    df['位点号'] = pd.to_numeric(df['位点号'])
    df['甲基化读段数'] = pd.to_numeric(df['甲基化读段数'])
    df['非甲基化读段数'] = pd.to_numeric(df['非甲基化读段数'])
    # 将df中的染色体号换为对应的chr_series中的索引，位点号和两个读段数换成数值类型，此外把甲基化类型CG转换为CpG
    df['甲基化类型'] = df['甲基化类型'].str.replace('CG', 'CpG')
    # 对三类数据进行分组，后续分别操作
    data_groups = df.groupby('甲基化类型')
    # 根据分组的键、对应的子数据表进行遍历
    for methy_type, data_ind in data_groups:
        if data_ind.empty:
            continue
        # 获取当前甲基化类型实际存在的染色体号并排序，
        # 是当前甲基化类型有的染色体号，注意：！！！是转换后从零开始的下标！！！
        actual_chrs = sorted(data_ind['染色体号'].dropna().unique())
        # 同时获取总染色体数，无论该甲基化类型是否具有对应的染色体信息
        chr_count = len(chr_series)
        chr_data_dict = {}
        # mlen用于记录当前甲基化类型中所有 染色体数据条数 中的最大值
        mlen = 0

        for chr_num in actual_chrs:  # 遍历该甲基化类型的每个染色体号（从零开始的数值）
            # 将当前染色体号的所有数据按照位点号进行升序排序，存储在chr_data中
            chr_data = data_ind[data_ind['染色体号'] == chr_num].sort_values('位点号').reset_index(drop=True)
            # 将三个关心的数据 位点号-甲基化读段数-非甲基化读段数 以染色体号为键存储在chr_data_dict字典中
            # 每个键对应的值是该三个属性组成的numpy数组
            chr_data_dict[chr_num] = chr_data[['位点号', '甲基化读段数', '非甲基化读段数']].values
            # 更新可能存在的更大的染色体数据条数到mlen中
            mlen = max(mlen, len(chr_data_dict[chr_num]))

        # 创建输出矩阵，行数为最大的染色体数据条数，列数为两个目录中所有存在的不重复染色体数*3，先用0填充
        output_matrix = np.zeros((mlen, chr_count * 3), dtype=np.int32)

        # 遍历输出矩阵中的所有行
        for i in range(mlen):
            # 遍历当前甲基化类型有的染色体号（从零开始编号的数值）
            for chr_num in actual_chrs:  # 使用连续索引
                col_start = chr_num * 3  # 基于连续索引计算列位置
                if i < len(chr_data_dict[chr_num]):
                    output_matrix[i, col_start:col_start + 3] = chr_data_dict[chr_num][i]

        # 输出（用pandas导出）
        output_df = pd.DataFrame(output_matrix)  # 将输出矩阵转换为DataFrame方便输出
        output_file = f"{num1}-bothMeUnme_diffChromo_NOREPEATED_methy_sites_{methy_type}.txt"
        output_path = os.path.join(output_dir, output_file)
        output_df.to_csv(output_path, sep='\t', header=False, index=False)

def single_newtoboth(filepath1, output_dir, num1, chr_series):
    """
    加速版 newtoboth 单文件转换。
    保持原输出格式：
    每个 methylation type 输出一个 bothMeUnme_diffChromo_NOREPEATED 文件；
    每个染色体占 3 列：Position, methylated reads, unmethylated reads。
    """

    col_names = [
        "染色体号",
        "位点号",
        "甲基化读段数",
        "非甲基化读段数",
        "甲基化类型",
    ]

    df = pd.read_csv(
        filepath1,
        sep=r"\s+",
        header=None,
        names=col_names,
        dtype={
            "染色体号": "string",
            "位点号": "int64",
            "甲基化读段数": "int32",
            "非甲基化读段数": "int32",
            "甲基化类型": "string",
        },
    )

    chr_map = chr_series.to_dict()

    df["染色体号"] = _normalize_chr_series(df["染色体号"]).map(chr_map)
    df = df.dropna(subset=["染色体号"]).copy()
    df["染色体号"] = df["染色体号"].astype(np.int32)

    # 只把真正的 CG 改成 CpG，避免无意替换其他字符串
    df["甲基化类型"] = df["甲基化类型"].replace({
        "CG": "CpG",
        "cg": "CpG",
        "cG": "CpG",
        "Cg": "CpG",
    })

    chr_count = len(chr_series)
    out_dtype = np.int32

    for methy_type, data_ind in df.groupby("甲基化类型", sort=False):
        if data_ind.empty:
            continue

        # 一次排序，避免每个染色体单独 sort
        data_ind = data_ind.sort_values(
            ["染色体号", "位点号"],
            kind="mergesort",
        )

        chr_blocks = []
        mlen = 0

        for chr_num, chr_df in data_ind.groupby("染色体号", sort=True):
            arr = chr_df[
                ["位点号", "甲基化读段数", "非甲基化读段数"]
            ].to_numpy(dtype=out_dtype, copy=True)

            chr_num = int(chr_num)
            chr_blocks.append((chr_num, arr))
            if arr.shape[0] > mlen:
                mlen = arr.shape[0]

        if mlen == 0:
            continue

        output_matrix = np.zeros((mlen, chr_count * 3), dtype=out_dtype)

        # 关键加速点：按染色体整块赋值，而不是逐行逐染色体赋值
        for chr_num, arr in chr_blocks:
            col_start = chr_num * 3
            output_matrix[:arr.shape[0], col_start:col_start + 3] = arr

        output_file = (
            f"{num1}-bothMeUnme_diffChromo_NOREPEATED_"
            f"methy_sites_{methy_type}.txt"
        )
        output_path = os.path.join(output_dir, output_file)

        pd.DataFrame(output_matrix).to_csv(
            output_path,
            sep="\t",
            header=False,
            index=False,
        )

def get_chr_name(chr_num, chr_series):
    """
    根据染色体序号获取染色体名称
    参数：
        chr_num: 第几个染色体（从1开始）
        chr_series: 染色体映射 Series
    返回：
        染色体名称字符串（如 'chr1', 'chrX'）
    """
    chr_num = int(chr_num)
    if chr_series is not None:
        try:
            # chr_num 是从1开始的，所以要减1
            if 0 <= chr_num - 1 < len(chr_series):
                return chr_series.index[chr_num - 1]
        except Exception as e:
            print(f"警告：获取染色体名称失败: {e}")

    # 如果出错或没有映射，返回数字编号
    return f"chr{chr_num}"

def newtoboth(m, n, dir1, dir2, threads=1, work_dir="."):
    """Convert raw methylation files to bothMeUnme matrix files.

    Parallel design:
      1. Keep chromosome mapping construction serial to guarantee one shared chr_series.
      2. If threads > 1, convert the m+n independent raw files concurrently.
         Each worker writes only its own replicate-numbered output files, so filenames do not collide.
    """
    # 没必要检查目录是否存在了，因为main函数主逻辑检查过了
    # 获取两个基因型目录中所有 不重复的 且排序了的 染色体号->数值的Series
    chr_series = scan_all_files_for_chr_mapping_fast(m, n, dir1, dir2)
    print(f"映射关系: {chr_series}")

    threads = max(1, int(threads))
    tasks = []

    # 循环m+n次，分别处理两个基因型的文件。任务顺序保持原版：dir1 全部在前，dir2 全部在后。
    for i in range(1, m + 1):
        filepath = os.path.join(dir1, f"{i}-{os.path.basename(dir1)}.txt")
        if not os.path.exists(filepath):
            print(f"警告: 文件不存在 {filepath}")
            continue
        tasks.append({
            "filepath": filepath,
            "output_dir": dir1,
            "num": i,
            "chr_series": chr_series,
            "label": f"{os.path.basename(dir1)}{i}",
        })

    for j in range(1, n + 1):
        filepath = os.path.join(dir2, f"{j}-{os.path.basename(dir2)}.txt")
        if not os.path.exists(filepath):
            print(f"警告: 文件不存在 {filepath}")
            continue
        tasks.append({
            "filepath": filepath,
            "output_dir": dir2,
            "num": j,
            "chr_series": chr_series,
            "label": f"{os.path.basename(dir2)}{j}",
        })

    if threads <= 1 or len(tasks) <= 1:
        for task in tasks:
            print(f"处理文件{task['filepath']}")
            single_newtoboth(task["filepath"], task["output_dir"], task["num"], chr_series)
    else:
        max_workers = min(threads, len(tasks))
        log_dir = os.path.join(work_dir, "parallel_logs", "newtoboth")
        os.makedirs(log_dir, exist_ok=True)
        for task in tasks:
            safe_label = sanitize_filename(task["label"])
            task["log_file"] = os.path.join(log_dir, f"newtoboth_{safe_label}.log")

        print(f"启用并行 newtoboth 文件转换: workers={max_workers}, tasks={len(tasks)}")
        with ProcessPoolExecutor(max_workers=max_workers) as executor:
            future_to_task = {executor.submit(_run_newtoboth_worker, task): task for task in tasks}
            for idx, future in enumerate(as_completed(future_to_task), 1):
                task = future_to_task[future]
                try:
                    res = future.result()
                    print(
                        f"[DONE] newtoboth {idx}/{len(tasks)} {res['label']}: "
                        f"{res['elapsed']:.2f}s, log={res['log_file']}"
                    )
                except Exception as e:
                    print(f"[FAILED] newtoboth {idx}/{len(tasks)} {task['label']}: {e}")
                    raise

    return chr_series

def sanitize_filename(name):
    """清理文件名中的不允许在文件名中出现的特殊字符"""
    return re.sub(r'[\\/*?:"<>|]', "", name)  # 将name中的特殊字符替换为空，以使其无害化


def get_column_count(file_path):
    """获取文件的列数并返回列数除以3的结果"""
    try:
        with open(file_path, 'r') as file:
            first_line = file.readline().strip()
            column_count = len(first_line.split())
            return column_count // 3
    except Exception as e:
        print(f"读取文件时发生错误: {e}")
        return None


def parse_filename(filename):
    """解析单个文件名，使用正则表达式提取文件编号和甲基化类型"""
    # 用捕获组（这里有括号的部分）捕获对应的文件序号和甲基化类型
    pattern = r'^(\d+)-bothMeUnme_diffChromo_NOREPEATED_methy_sites_(.+)\.txt$'
    match = re.match(pattern, filename)  # 将特定的文件名和此正则表达式进行匹配
    if match:  # 匹配到后，用match.group(i)获取对应第i个捕获组的内容
        file_id = int(match.group(1))
        methylation_type = match.group(2)
        return file_id, methylation_type
    return None


def scan_sample_files_by_replicates(sample_dir, max_replicates):
    """在sample_dir目录中，搜寻出序号 1~max_replicates-both...甲基化类型.txt 文件"""
    files_by_replicates = {}  # 创建字典用于根据前端的序号存储对应的文件名
    methylation_types = ['CpG', 'CHH', 'CHG']
    # 若不存在对应目录，则返回空字典，这步合并后可省略，因为main函数中一开始就判断了
    # if not os.path.exists(sample_dir):
    #     return files_by_replicates

    # 遍历该文件夹中的所有内容（包括文件和子目录）
    for filename in os.listdir(sample_dir):
        if filename.endswith('.txt'):  # 若找到txt文件（只有可能是both文件或新格式文件）
            parsed = parse_filename(filename)  # 解析该文件获取可能存在的 序号,甲基化类型
            if parsed:  # 如果对于该文件解析到了 序号，甲基化类型
                file_id, methylation_type = parsed  # （元组，逗号是元组的标志，括号只是为了为避免歧义的产物）
                if file_id <= max_replicates and methylation_type in methylation_types:  # 序号和甲基化类型都合法
                    if file_id not in files_by_replicates:
                        files_by_replicates[file_id] = {}  # 使得files_by_replicates这个字典成为 file_id->{子字典}
                    files_by_replicates[file_id][methylation_type] = filename  # 使得子字典中格式为 甲基化类型->filename
    return files_by_replicates  # 就是返回一个链路 序号->甲基化类型->文件名 ，可以通过序号和甲基化类型获取到对应文件名


def process_methylation_type_with_collection(file1_path, file2_path, methylation_type, output_dir,
                                             dir1_name, dir2_name, chr_num, meth_diff_threshold=0.0):
    """处理甲基化类型和染色体的Fisher检验 处理1次检验的一个染色体
    参数解释：file_path-当前处理的both文件的相对路径
    methylation_type-当前处理的甲基化类型
    output_dir-输出文件夹./output_x_y/
    dir_name-输入的两个基因型目录路径的最后一个部分
    chr_num-处理的第几个染色体，而不是染色体号"""

    print(f"      处理甲基化类型 {methylation_type}，染色体 {chr_num}...")

    # 计算当前染色体的数据列
    mOrder = 3 * (chr_num - 1)  # 这里的chr_num是第几个染色体，而不是染色体号！！！
    usecols = [mOrder, mOrder + 1, mOrder + 2]

    try:
        # 分块读取第一个文件的所有数据到字典
        data1_dict = {}
        for chunk in pd.read_csv(file1_path, sep=r'\s+', header=None, usecols=usecols, chunksize=100000):
            chunk.columns = ['pos', 'methy', 'unmethy']
            # 用zip每次获取对应列的一个元素，共三个元素放到(pos,methy,unmethy)元组中
            for pos, methy, unmethy in zip(chunk['pos'].values, chunk['methy'].values, chunk['unmethy'].values):
                # 按照如下格式将数据放入对应字典中
                data1_dict[pos] = (methy, unmethy)
        # 至此第一个文件的所有位点数据都录入了data1_dict字典中

        # 分块读取第二个文件并查找共同位点，这是因为如果同时将两个文件的所有数据读到内存里，可能会占用太大的内存
        # 而第一个文件必须全部载入内存，因为我们需要对它进行快速的随机查找以确认某个位点是否在两个文件中都有，都有的话就得检验
        # 这里reader由于有chunksize这个参数，所以read_csv的返回值是一个迭代器，遍历它每次可以每次最多返回100000行数据
        reader2 = pd.read_csv(file2_path, sep=r'\s+', header=None, usecols=usecols, chunksize=100000)

    except Exception as e:
        print(f"        加载数据失败: {e}")
        return False

    # 创建输出文件夹:output_x_y/甲基化类型/
    methylation_dir = os.path.join(output_dir, methylation_type)
    os.makedirs(methylation_dir, exist_ok=True)

    # 创建输出文件路径
    stats_filename = os.path.join(methylation_dir,
                                  f"FET_results_{methylation_type}_{dir1_name}_{dir2_name}_Chr{chr_num}.txt")
    all_output = os.path.join(methylation_dir, f"all_simple_Chr{chr_num}.txt")
    sig_output = os.path.join(methylation_dir, f"sig_simple_Chr{chr_num}.txt")
    combine_output = os.path.join(methylation_dir, f"combineResult_Chr{chr_num}.txt")

    # 设置参数
    M0, M2, Th_pValue = 2, 4, 0.05

    try:
        # 创建四个列表用于存储各类数据以便后续输出
        all_results, sig_results, fet_results, combine_results = [], [], [], []

        # 处理第二个文件的每个数据块
        for chunk2 in reader2:
            chunk2.columns = ['pos', 'methy', 'unmethy']
            # 对于第二个文件的每个数据块，遍历每一行三个值，存储到pos,m2,u2中
            for pos, m2, u2 in zip(chunk2['pos'].values, chunk2['methy'].values, chunk2['unmethy'].values):
                # 若该位点在第一个文件也存在，说明要进行该位点的fisher检验
                if pos in data1_dict:
                    m1, u1 = data1_dict[pos]  # 获取第一个文件的两个数

                    if m1 >= M0 or m2 >= M0:  # 两个methy读段数都要>=2才可能进行下一步
                        if (m1 + u1 >= M2) and (m2 + u2 >= M2):  # 每个文件的两个数之和都要>=4才进行检验
                            cont_table = np.array([[m1, u1], [m2, u2]])  # 建立2*2列联表
                            # 计算两个甲基化率
                            ratio1 = m1 / (m1 + u1) if (m1 + u1) > 0 else 0
                            ratio2 = m2 / (m2 + u2) if (m2 + u2) > 0 else 0
                            meth_diff = abs(ratio1 - ratio2)
                            change = "1" if ratio1 >= ratio2 else "0"  # 根据两个文件的甲基化比率判断
                            # 突变型的甲基化率是否升高了
                            # 调用库函数进行fisher检验，获取到该次检验的pvalue
                            _, pvalue = fisher_exact(cont_table, alternative='two-sided')
                            # 为pvalue保留7位有效数字
                            pvalue = float(f"{pvalue:.7g}")
                            # 为四个文件分别录入需要的数据，其中只有显著的才录入sig_results中
                            all_results.append([pos, pvalue, change, meth_diff])
                            fet_results.append([pos, m1, u1, m2, u2, pvalue, meth_diff])
                            combine_results.append([pos, m1, u1, m2, u2, meth_diff])

                            if pvalue < Th_pValue:
                                sig_results.append([pos, pvalue, change, meth_diff])
        # 至此，这组文件（特定x_y特定甲基化类型特定染色体）所有需要进行的fisher检验已完成
        # 保存结果到磁盘中
        if all_results:
            pd.DataFrame(all_results, columns=["Position", "Pvalue", "Methylation_Change", "MethDiff"]).to_csv(all_output, sep='\t',
                                                                                                   index=False)
            pd.DataFrame(sig_results, columns=["Position", "Pvalue", "Methylation_Change", "MethDiff"]).to_csv(sig_output, sep='\t',
                                                                                                   index=False)
            pd.DataFrame(fet_results,
                         columns=["Position", "Methy1", "Unmethy1", "Methy2", "Unmethy2", "Pvalue", "MethDiff"]).to_csv(
                stats_filename, sep='\t', index=False)
            pd.DataFrame(combine_results, columns=["Position", "Methy1", "Unmethy1", "Methy2", "Unmethy2", "MethDiff"]).to_csv(
                combine_output, sep='\t', index=False)

        print(f"        染色体 {chr_num} 处理完成！共处理 {len(all_results)} 个位点，其中 {len(sig_results)} 个显著")
        return True

    except Exception as e:
        print(f"        处理甲基化类型 {methylation_type}，染色体 {chr_num} 时发生错误: {e}")
        return False

def process_methylation_type_with_collection_pvfilter(file1_path, file2_path, methylation_type, output_dir,
                                             dir1_name, dir2_name, chr_num, meth_diff_threshold=0.0):
    """处理甲基化类型和染色体的Fisher检验 处理1次检验的一个染色体
    参数解释：file_path-当前处理的both文件的相对路径
    methylation_type-当前处理的甲基化类型
    output_dir-输出文件夹./output_x_y/
    dir_name-输入的两个基因型目录路径的最后一个部分
    chr_num-处理的第几个染色体，而不是染色体号

    返回此次的pvalue>0.05的df
        all_results.append([pos, pvalue, change])
    """

    print(f"      处理甲基化类型 {methylation_type}，染色体 {chr_num}...")

    # 计算当前染色体的数据列
    mOrder = 3 * (chr_num - 1)  # 这里的chr_num是第几个染色体，而不是染色体号！！！
    usecols = [mOrder, mOrder + 1, mOrder + 2]

    try:
        # 分块读取第一个文件的所有数据到字典
        data1_dict = {}
        for chunk in pd.read_csv(file1_path, sep=r'\s+', header=None, usecols=usecols, chunksize=100000):
            chunk.columns = ['pos', 'methy', 'unmethy']
            # 用zip每次获取对应列的一个元素，共三个元素放到(pos,methy,unmethy)元组中
            for pos, methy, unmethy in zip(chunk['pos'].values, chunk['methy'].values, chunk['unmethy'].values):
                # 按照如下格式将数据放入对应字典中
                data1_dict[pos] = (methy, unmethy)
        # 至此第一个文件的所有位点数据都录入了data1_dict字典中

        # 分块读取第二个文件并查找共同位点，这是因为如果同时将两个文件的所有数据读到内存里，可能会占用太大的内存
        # 而第一个文件必须全部载入内存，因为我们需要对它进行快速的随机查找以确认某个位点是否在两个文件中都有，都有的话就得检验
        # 这里reader由于有chunksize这个参数，所以read_csv的返回值是一个迭代器，遍历它每次可以每次最多返回100000行数据
        reader2 = pd.read_csv(file2_path, sep=r'\s+', header=None, usecols=usecols, chunksize=100000)

    except Exception as e:
        print(f"        加载数据失败: {e}")
        return False

    # 创建输出文件夹:output_x_y/甲基化类型/
    methylation_dir = os.path.join(output_dir, methylation_type)
    os.makedirs(methylation_dir, exist_ok=True)

    # 创建输出文件路径
    stats_filename = os.path.join(methylation_dir,
                                  f"FET_results_{methylation_type}_{dir1_name}_{dir2_name}_Chr{chr_num}.txt")
    all_output = os.path.join(methylation_dir, f"all_simple_Chr{chr_num}.txt")
    sig_output = os.path.join(methylation_dir, f"sig_simple_Chr{chr_num}.txt")
    combine_output = os.path.join(methylation_dir, f"combineResult_Chr{chr_num}.txt")

    # 设置参数
    M0, M2, Th_pValue = 2, 4, 0.05

    try:
        # 创建四个列表用于存储各类数据以便后续输出
        all_results, sig_results, fet_results, combine_results = [], [], [], []

        # 处理第二个文件的每个数据块
        for chunk2 in reader2:
            chunk2.columns = ['pos', 'methy', 'unmethy']
            # 对于第二个文件的每个数据块，遍历每一行三个值，存储到pos,m2,u2中
            for pos, m2, u2 in zip(chunk2['pos'].values, chunk2['methy'].values, chunk2['unmethy'].values):
                # 若该位点在第一个文件也存在，说明要进行该位点的fisher检验
                if pos in data1_dict:
                    m1, u1 = data1_dict[pos]  # 获取第一个文件的两个数

                    if m1 >= M0 or m2 >= M0:  # 两个methy读段数都要>=2才可能进行下一步
                        if (m1 + u1 >= M2) and (m2 + u2 >= M2):  # 每个文件的两个数之和都要>=4才进行检验
                            cont_table = np.array([[m1, u1], [m2, u2]])  # 建立2*2列联表
                            # 计算两个甲基化率
                            ratio1 = m1 / (m1 + u1) if (m1 + u1) > 0 else 0
                            ratio2 = m2 / (m2 + u2) if (m2 + u2) > 0 else 0
                            meth_diff = abs(ratio1 - ratio2)
                            change = "1" if ratio1 > ratio2 else "0"  # 根据两个文件的甲基化比率判断
                            # 突变型的甲基化率是否升高了
                            # 调用库函数进行fisher检验，获取到该次检验的pvalue
                            _, pvalue = fisher_exact(cont_table, alternative='two-sided')
                            # 为pvalue保留7位有效数字
                            pvalue = float(f"{pvalue:.7g}")
                            # 为四个文件分别录入需要的数据，其中只有显著的才录入sig_results中
                            all_results.append([pos, pvalue, change, meth_diff])
                            fet_results.append([pos, m1, u1, m2, u2, pvalue, meth_diff])
                            combine_results.append([pos, m1, u1, m2, u2, meth_diff])

                            if pvalue < Th_pValue:
                                sig_results.append([pos, pvalue, change, meth_diff])
        # 至此，这组文件（特定x_y特定甲基化类型特定染色体）所有需要进行的fisher检验已完成
        # 保存结果到磁盘中
        if all_results:
            # 转换为DataFrame方便筛选
            all_df = pd.DataFrame(all_results, columns=["Position", "Pvalue", "Methylation_Change", "MethDiff"])
            sig_df = pd.DataFrame(sig_results, columns=["Position", "Pvalue", "Methylation_Change", "MethDiff"])
            fet_df = pd.DataFrame(fet_results,
                                  columns=["Position", "Methy1", "Unmethy1", "Methy2", "Unmethy2", "Pvalue", "MethDiff"])
            combine_df = pd.DataFrame(combine_results, columns=["Position", "Methy1", "Unmethy1", "Methy2", "Unmethy2", "MethDiff"])

            all_df_ndmp = all_df[all_df['Pvalue'] > 0.05]  # 保留pvalue>0.05的信息，后续补充到qvalue列表里

            # 只保留 p≤0.05 的位点
            all_df = all_df[all_df['Pvalue'] <= 0.05]  # 筛选pvalue<=0.05的进行FDR检验
            fet_df = fet_df[fet_df['Pvalue'] <= 0.05]  # 筛选pvalue<=0.05的进行FDR检验
            #fet_df_ndmp = all_df[all_df['Pvalue'] > 0.05]
            # sig_df 本来就是 p<0.05，不需要再筛选

            # 也筛选 combine
            positions_to_keep = set(all_df['Position'].values)
            combine_df = combine_df[combine_df['Position'].isin(positions_to_keep)]  # 筛选pvalue<=0.05的进行FDR检验
            #combine_df_ndmp = combine_df[~combine_df['Position'].isin(positions_to_keep)]
            # 保存筛选后的结果
            all_df.to_csv(all_output, sep='\t', index=False)
            sig_df.to_csv(sig_output, sep='\t', index=False)
            fet_df.to_csv(stats_filename, sep='\t', index=False)
            combine_df.to_csv(combine_output, sep='\t', index=False)

        print(
            f"        染色体 {chr_num} 处理完成！原始检验 {len(all_results)} 个位点，筛选后(p≤0.05) {len(all_df)} 个位点，p>0.05的共{len(all_df_ndmp)}个位点")
        return all_df_ndmp
    #                   all_results.append([pos, pvalue, change])
    #                   fet_results.append([pos, m1, u1, m2, u2, pvalue])
    #               combine_results.append([pos, m1, u1, m2, u2])

    except Exception as e:
        print(f"        处理甲基化类型 {methylation_type}，染色体 {chr_num} 时发生错误: {e}")
        return False

def merge_fet_results_and_fdr(output_dir, replicate_x, replicate_y, mtype3, all_dfs_ndmp_dict, n_chromosomes, is_twostep_context=False):
    """合并output文件夹中的所有FET结果并进行FDR校正
    其中FET文件的格式如：pos, m1, u1, m2, u2, pvalue
    这里输入的output_dir是output_x_y/甲基化类型
    success_dfs_dict[methylation_type][chr_num]访问的是
                    对应次pvalue>0.05的df
                      all_results([pos, pvalue, change])
    n_chromosomes是总染色体数
                      """
    print(f"\n    合并 {output_dir} 中的FET结果并进行FDR校正...")

    if not os.path.exists(output_dir):
        print(f"    错误：目录 {output_dir} 不存在！")
        return False, get_dmp_threshold(mtype3)

    # 搜索所有FET结果文件
    # 这里**表示任意深度的子目录，所以这里glob在递归搜索（recursive=True）时就会在output_dir目录下、以及其
    # 所有子目录下搜寻符合条件的文件，并将符合条件的文件路径（从output_dir开始的相对路径）以列表的形式返回
    file_pattern = os.path.join(output_dir, "**", "FET_results_*_Chr*.txt") # 这里的染色体号实际上是both文件中的第几个三列
    fet_files = glob.glob(file_pattern, recursive=True)

    if not fet_files:
        print(f"    警告：在 {output_dir} 中未找到FET结果文件")
        return False, get_dmp_threshold(mtype3)

    print(f"    找到 {len(fet_files)} 个FET结果文件")

    # 创建列表用于收集所有p值和相关信息
    all_data = []

    for file_path in sorted(fet_files):
        # 提取甲基化类型和染色体信息
        # 这里第一个捕获组用于捕获甲基化类型，.*为零次或多次任意字符，用于匹配replicatex_replicatey，第二个捕获组用于捕获染色体号
        #                                                   # 这里的染色体号实际上是both文件中的第几个三列
        methy_match = re.search(r'/FET_results_([^_]+)_.*_Chr(\d+)\.txt$', file_path.replace('\\', '/'))
        if not methy_match:
            continue

        # 获取到甲基化类型和染色体序号（这里的染色体号是both文件中的第几个三列）
        methylation_type = mtype3
        chr_num = int(methy_match.group(2))

        try:
            df = pd.read_csv(file_path, sep='\t', header=0)
            if 'Position' in df.columns and 'Pvalue' in df.columns:
                # 调整列顺序：Chromosome, Methylation_Type, Position, Pvalue
                df_subset = df[['Position', 'Pvalue', 'MethDiff']].copy()
                df_subset['Chromosome'] = chr_num
                df_subset['Methylation_Type'] = methylation_type

                # 重新排列列顺序
                df_subset = df_subset[['Chromosome', 'Methylation_Type', 'Position', 'Pvalue', 'MethDiff']]
                all_data.append(df_subset) # 将当前FET的所有信息加入all_data中
        except Exception as e:
            print(f"    警告：读取 {file_path} 失败: {e}")
            continue

    if not all_data:
        print(f"    错误：没有成功读取任何数据")
        return False, get_dmp_threshold(mtype3)

    # 合并所有数据（由于上面是直接将每个文件对应的dataframe追加到all_data末尾，所以列表中每个元素都是一个dataframe，因此需要合并）
    combined_df = pd.concat(all_data, ignore_index=True)
    # 将数据根据染色体序号排序（这里的染色体号实际上是both文件中的第几个三列）
    combined_df = combined_df.sort_values(['Methylation_Type', 'Chromosome', 'Position'])
    # 该dataframe格式如下：'Chromosome', 'Methylation_Type', 'Position', 'Pvalue'
    print(f"    合并后总共 {len(combined_df)} 个位点")

    # 计算FDR校正的q值并将其向右增加到combined_df中
    pvalues = combined_df['Pvalue'].values
    qvalues = calculate_qvalues(pvalues, 1.0)
    combined_df['Qvalue'] = qvalues

    fixed_dmp_threshold = get_dmp_threshold(mtype3)
    dmp_threshold_to_use = fixed_dmp_threshold
    auto_q_info = None

    # 自动q-value阈值只用于“两步法”上下文：
    # 即该context已经先按 p<=AUTO_QVALUE_P_CUTOFF 预筛选，再在预筛子集内计算FDR。
    # 非两步法context不启用，避免把全量FDR的q值曲线误用于该规则。
    should_estimate_auto_q = (AUTO_QVALUE_TWOSTEP or AUTO_QVALUE_REPORT_ONLY) and bool(is_twostep_context)
    if should_estimate_auto_q:
        estimated_q, auto_q_info = estimate_auto_qvalue_threshold_twostep(
            pvalues=combined_df['Pvalue'].values,
            qvalues=combined_df['Qvalue'].values,
            pvalue_cutoff=AUTO_QVALUE_P_CUTOFF,
            fallback_q=fixed_dmp_threshold,
            min_candidates=AUTO_QVALUE_MIN_CANDIDATES,
            use_smooth=AUTO_QVALUE_USE_SMOOTH,
            smooth_sigma=AUTO_QVALUE_SMOOTH_SIGMA,
        )

        if AUTO_QVALUE_TWOSTEP and not AUTO_QVALUE_REPORT_ONLY:
            dmp_threshold_to_use = estimated_q
            mode_msg = "用于DMP判定"
        else:
            dmp_threshold_to_use = fixed_dmp_threshold
            mode_msg = "仅报告，不改变DMP判定"

        print(
            f"    自动q-value阈值({mtype3}, 两步法context): "
            f"estimated={estimated_q:.6g}, used={dmp_threshold_to_use:.6g} "
            f"({mode_msg}, status={auto_q_info.get('auto_q_status')}, "
            f"p_at_max={auto_q_info.get('auto_q_pvalue_at_max')}, "
            f"diff_at_max={auto_q_info.get('auto_q_diff_at_max')}, "
            f"n={auto_q_info.get('auto_q_n_candidates')})"
        )
    elif AUTO_QVALUE_TWOSTEP or AUTO_QVALUE_REPORT_ONLY:
        print(f"    {mtype3} 不是两步法context：不进行自动q-value阈值估计，使用固定阈值 {fixed_dmp_threshold}")

    combined_df['Qvalue_Threshold_Used'] = dmp_threshold_to_use
    combined_df['Qvalue_Threshold_Mode'] = (
        "auto_twostep" if (AUTO_QVALUE_TWOSTEP and not AUTO_QVALUE_REPORT_ONLY and bool(is_twostep_context))
        else "fixed"
    )

    # 最终列顺序：Chromosome, Methylation_Type, Position, Pvalue, Qvalue
    combined_df = combined_df[['Chromosome', 'Methylation_Type', 'Position', 'Pvalue', 'MethDiff', 'Qvalue',
                               'Qvalue_Threshold_Used', 'Qvalue_Threshold_Mode']]

    print(f"开始处理 {mtype3} 的 ndmp 数据，")

    # 首先检查 all_dfs_ndmp_dict 中该甲基化类型是否有数据，因为有可能该生物种类的该甲基化类型不需要调用filter版本
    if mtype3 in all_dfs_ndmp_dict and all_dfs_ndmp_dict[mtype3]:
        print(f"该甲基化类型中，ndmp的dfs数量为{len(all_dfs_ndmp_dict[mtype3])}个")
        dfs_ndmp = []
        for chr_num11 in range(1, n_chromosomes + 1):
            # 只处理该甲基化类型中存在的染色体
            if chr_num11 in all_dfs_ndmp_dict[mtype3]:
                df_ndmp = all_dfs_ndmp_dict[mtype3][chr_num11]

                # 确保是有效的 DataFrame
                if isinstance(df_ndmp, pd.DataFrame) and not df_ndmp.empty:
                    df_ndmp = df_ndmp.copy()

                    df_ndmp['Chromosome'] = chr_num11
                    df_ndmp['Methylation_Type'] = mtype3

                    if 'change' in df_ndmp.columns:
                        df_ndmp.drop(columns=['change'], inplace=True)

                    df_ndmp['Qvalue'] = 1
                    df_ndmp['Qvalue_Threshold_Used'] = dmp_threshold_to_use
                    df_ndmp['Qvalue_Threshold_Mode'] = (
                        "auto_twostep" if (AUTO_QVALUE_TWOSTEP and not AUTO_QVALUE_REPORT_ONLY and bool(is_twostep_context))
                        else "fixed"
                    )
                    df_ndmp = df_ndmp[['Chromosome', 'Methylation_Type', 'Position', 'Pvalue', 'MethDiff', 'Qvalue',
                                       'Qvalue_Threshold_Used', 'Qvalue_Threshold_Mode']]

                    dfs_ndmp.append(df_ndmp)
                    print(f"  为 {mtype3}--{chr_num11} 添加所需列，共 {len(df_ndmp)} 个位点")

        if dfs_ndmp:
            total_ndmp = pd.concat(dfs_ndmp, ignore_index=True)
            combined_df = pd.concat([combined_df, total_ndmp], ignore_index=True)
            print(f"  共合并 {len(total_ndmp)} 个 ndmp 位点")
        else:
            print(f"  {mtype3} 类型未找到任何 ndmp 数据")
    else:
        print(f"  {mtype3} 类型在 ndmp 字典中不存在")
    # 至此对于比较多的甲基化类型已加上前面舍弃的pvalue>0.05的位点的信息，
    #   在前面要将其存储在all_dfs_ndmp_dict[甲基化][染色体号]里，此时才能实现该目的

    # 统计显著性结果，计算这 m*n*3次检验中的1次中，pvalues和qvalues显著的分别有多少个
    n_pval_sig = np.sum(pvalues <= 0.05)
    dmp_threshold = dmp_threshold_to_use
    n_qval_sig = np.sum(qvalues <= dmp_threshold)
    # 计算显著的比例
    print(f"    P值显著位点: {n_pval_sig} ({n_pval_sig / len(pvalues) * 100:.1f}%)")
    print(f"    Q值显著位点: {n_qval_sig} ({n_qval_sig / len(qvalues) * 100:.1f}%)")

    # 保存合并的p值列表（用于外部FDR工具）
    pvalue_file = os.path.join(output_dir, f"united_pvalues_wt_replicate{replicate_y}_vs_mut_replicate{replicate_x}.csv")
    with open(pvalue_file, 'w') as f:
        for pvalue in pvalues:
            f.write(f"{pvalue}\n")

    # 保存完整的FDR校正结果（output_dir是output_x_y/甲基化类型）
    #  包含完整的pvalues和qvalues（格式为：Chromosome, Methylation_Type, Position, Pvalue, Qvalue）
    fdr_file = os.path.join(output_dir, f"FDR_corrected_results_wt_replicate{replicate_y}_vs_mut_replicate{replicate_x}.txt")
    combined_df.to_csv(fdr_file, sep='\t', index=False)

    # 保存显著位点（根据甲基化类型使用固定阈值或两步法自动阈值）
    sig_df = combined_df[combined_df['Qvalue'] <= dmp_threshold_to_use]
    if not sig_df.empty: # 若有显著的位点，将显著的那部分数据输出（output_dir是output_x_y/甲基化类型）
        sig_file = os.path.join(output_dir, f"FDR_significant_results_wt_replicate{replicate_y}_vs_mut_replicate{replicate_x}.txt")
        sig_df.to_csv(sig_file, sep='\t', index=False)
        print(f"    显著位点结果保存至: {sig_file}")

    if should_estimate_auto_q:
        save_auto_qvalue_report(
            output_dir=output_dir,
            replicate_x=replicate_x,
            replicate_y=replicate_y,
            methylation_type=mtype3,
            info=auto_q_info,
            threshold_used=dmp_threshold_to_use,
            report_only=AUTO_QVALUE_REPORT_ONLY or not AUTO_QVALUE_TWOSTEP,
        )

    print(f"    P值列表保存至: {pvalue_file}")
    print(f"    FDR结果保存至: {fdr_file}")
    return True, dmp_threshold_to_use


# 在你的代码中，将calculate_qvalues函数替换为：
def     calculate_qvalues(pvalues, pi=1.0):
    """使用Storey方法计算Q值"""

    pvalues = np.array(pvalues, dtype=float)
    if len(pvalues) == 0:
        return np.array([])

    pvalues = np.array(pvalues)
    pvalues_clean = pvalues.copy()
    pvalues_clean[np.isnan(pvalues_clean)] = 1.0
    pvalues_clean = np.clip(pvalues_clean, 0, 1)

    m = len(pvalues_clean)
    sorted_indices = np.argsort(pvalues_clean)
    sorted_pvalues = pvalues_clean[sorted_indices]

    # Storey方法的关键：包含π因子
    q_values_sorted = np.zeros_like(sorted_pvalues)

    # 如果需要自动估计π
    if pi is None:
        # 简化的π估计
        lrange = np.linspace(0.05, 0.95, max(int(m / 100.0), 10))
        pil = np.mean(sorted_pvalues[:, np.newaxis] > lrange, axis=0)
        pilr = pil / (1.0 - lrange)
        pi = 1.0
        if pilr[-1] < 1.0:
            pi = pilr[-1]

    # Storey方法计算
    q_values_sorted = pi * m * sorted_pvalues / np.arange(1, m + 1)
    q_values_sorted[-1] = min(q_values_sorted[-1], 1.0)

    # 单调性调整
    for i in range(m - 2, -1, -1):
        q_values_sorted[i] = min(q_values_sorted[i], q_values_sorted[i + 1])

    # 恢复原顺序
    q_values = np.zeros_like(pvalues_clean)
    q_values[sorted_indices] = q_values_sorted
    q_values[np.isnan(pvalues)] = np.nan

    return q_values

def perform_sliding_window_on_dmp_files(output_dir, replicate_x, replicate_y,):
    """对m*n*3次处理中的一次的结果(N)DMP文件进行滑动窗口分析"""

    print(f"\n    开始对DMP文件进行滑动窗口分析...")

    # 只处理DMP文件
    dmp_patterns = [
        f"DMP_wt_replicate{replicate_y}_mut_replicate{replicate_x}_Chr*.txt",
        f"N-DMP_wt_replicate{replicate_y}_mut_replicate{replicate_x}_Chr*.txt"
    ]

    for pattern in dmp_patterns:
        dmp_files = glob.glob(os.path.join(output_dir, pattern)) #找到符合当前正则表达式格式的所有文件的路径，并放入dmp_files列表中

        for dmp_file in dmp_files: # 遍历列表，每次获取一个DMP或N-DMP文件
            try:

                df = pd.read_csv(dmp_file, sep=r'\s+', header=None, skiprows=1,
                                 names=['position', 'pvalue', 'change']) # 位点号、pvalue、change

                # 检验数据非空
                if df.empty or len(df) == 0:
                    continue

                # 数据类型转换和清理
                df = df.dropna()
                df['position'] = df['position'].astype(int)
                df['change'] = df['change'].astype(int)

                print(f"      处理文件: {os.path.basename(dmp_file)} ({len(df)} 个位点)")

                # 设置输出前缀，其实就是将.txt去掉之后的 (N)DMP_replicate_wt{replicate_y}_mut_replicate{replicate_x}_Chr* 这一部分
                base_name = os.path.basename(dmp_file).replace('.txt', '')

                # 传入DataFrame给sliding_window_analysis
                # 前者格式为window_start,window_end,count_change_1,count_change_0,total_count
                # 后者没有count_change，而有standardized_count，即每个区间显著位点数与最大的那个显著位点数的比值
                sliding_results, std_results = sliding_window_analysis(
                    df, # df是m*n*3次中的一次的某一条染色体的dmp文件
                    window_size=1000000,
                    step_ratio=0.05,
                    save_files=True,
                    output_identifier=base_name,
                    outputdir1=output_dir
                )

                print(f"      完成 {base_name}: {len(sliding_results)} 个窗口")

            except Exception as e:
                print(f"      处理 {dmp_file} 失败: {e}")
                continue

    print(f"    滑动窗口分析完成")

#下面这个版本用于处理pv<0.05筛选后再进行FDR检验的情况，比如植物的CHH和CHG
def perform_sliding_window_on_dmp_files_after_filter(output_dir, replicate_x, replicate_y,all_dfs_ndmp_dict=None, methylation_type=None):
    """对m*n*3次处理中的一次的结果(N)DMP文件进行滑动窗口分析"""

    print(f"\n    开始对DMP文件进行滑动窗口分析...")

    # 只处理DMP文件
    dmp_patterns = [
        f"DMP_wt_replicate{replicate_y}_mut_replicate{replicate_x}_Chr*.txt",
        f"N-DMP_wt_replicate{replicate_y}_mut_replicate{replicate_x}_Chr*.txt"
    ]

    for pattern in dmp_patterns:
        dmp_files = glob.glob(os.path.join(output_dir, pattern)) #找到符合当前正则表达式格式的所有文件的路径，并放入dmp_files列表中
        is_ndmp = pattern.startswith("N-DMP") # 判断当前处理的文件是否是NDMP文件

        for dmp_file in dmp_files: # 遍历列表，每次获取一个DMP或N-DMP文件
            try:

                df = pd.read_csv(dmp_file, sep=' ', header=None, skiprows=1,
                                 names=['position', 'pvalue', 'change']) # 位点号、pvalue、change

                #如果当前文件是ndmp文件，由于不能改all_simple文件，因为那个要用来合并成一列用于FDR检验，所以就在此加上原来pv>0.05的
                                # 那些位点的信息，这样的话，后续生成滑动窗口文件的时候，统计ndmp信息的时候才不会缺漏
                if is_ndmp and all_dfs_ndmp_dict and methylation_type:
                    # 从文件名提取染色体号
                    chr_match = re.search(r'Chr(\d+)\.txt$', dmp_file)
                    if chr_match:
                        chr_num = int(chr_match.group(1))
                        if (methylation_type in all_dfs_ndmp_dict and
                                chr_num in all_dfs_ndmp_dict[methylation_type]):
                            ndmp_df = all_dfs_ndmp_dict[methylation_type][chr_num]
                            if isinstance(ndmp_df, pd.DataFrame) and not ndmp_df.empty:
                                ndmp_df = ndmp_df.rename(columns={
                                    'Position': 'position',
                                    'Pvalue': 'pvalue',
                                    'Methylation_Change': 'change'
                                })
                                # 合并数据
                                df = pd.concat([df, ndmp_df[['position', 'pvalue', 'change']]],
                                               ignore_index=True)
                                df = df.drop_duplicates(subset=['position'])
                                df = df.sort_values('position').reset_index(drop=True)
                                print(f"      为 {os.path.basename(dmp_file)} 合并了之前FDR检验时忽略的 {len(ndmp_df)} 个NDMP位点")
                # 检验数据非空
                if df.empty or len(df) == 0:
                    continue

                # 数据类型转换和清理
                df = df.dropna()
                df['position'] = df['position'].astype(int)
                df['change'] = df['change'].astype(int)

                print(f"      处理文件: {os.path.basename(dmp_file)} ({len(df)} 个位点)")

                # 设置输出前缀，其实就是将.txt去掉之后的 (N)DMP_wt_replicate{replicate_y}_mut_replicate{replicate_x}_Chr* 这一部分
                base_name = os.path.basename(dmp_file).replace('.txt', '')

                # 传入DataFrame给sliding_window_analysis
                # 前者格式为window_start,window_end,count_change_1,count_change_0,total_count
                # 后者没有count_change，而有standardized_count，即每个区间显著位点数与最大的那个显著位点数的比值
                sliding_results, std_results = sliding_window_analysis(
                    df, # df是m*n*3次中的一次的某一条染色体的dmp文件
                    window_size=1000000,
                    step_ratio=0.05,
                    save_files=True,
                    output_identifier=base_name,
                    outputdir1=output_dir
                )

                print(f"      完成 {base_name}: {len(sliding_results)} 个窗口")

            except Exception as e:
                print(f"      处理 {dmp_file} 失败: {e}")
                continue

    print(f"    滑动窗口分析完成")


def generate_dmp_files(dir1,dir2,output_dir, replicate_x, replicate_y, fdr_threshold=0.05, mtype1="CpG",
                      all_dfs_ndmp_dict=None,unfilter_mtypes=["CpG"],n_chromosomes = 5, meth_diff_threshold=0.0,
                      skip_dmr=False, skip_window=False):
    """参数：output_dir是output_x_y/甲基化类型/,
        replicate_x,replicate_y是处理的组号，后面是qvalues的阈值和目前处理的甲基化类型"""

    print(f"\n    生成 {output_dir} 的DMP文件...")

    def safe_float_convert(value):
        """将各种格式的数据转换为浮点数"""
        try:
            # 如果已经是数字类型，均转换为浮点数返回
            if isinstance(value, (int, float)):
                return float(value)
            # 如果是字符串，去除空白后转换为浮点数
            if isinstance(value, str):
                value = value.strip()
                if value:
                    return float(value)
            return None
        except (ValueError, TypeError):
            return None

    # 1. 读取FDR校正结果文件(这里output_dir是output_x_y/甲基化类型/)
    #  该文件格式为 Chromosome, Methylation_Type, Position, Pvalue, Qvalue
    fdr_file = os.path.join(output_dir, f"FDR_corrected_results_wt_replicate{replicate_y}_vs_mut_replicate{replicate_x}.txt")
    if not os.path.exists(fdr_file):
        print(f"    错误：FDR结果文件不存在 {fdr_file}")
        return False

    try:
        fdr_df = pd.read_csv(fdr_file, sep='\t') # 该文件格式为 Chromosome, Methylation_Type, Position, Pvalue, Qvalue
        print(f"    读取FDR结果文件，共 {len(fdr_df)} 个位点")
    except Exception as e:
        print(f"    错误：读取FDR文件失败 {e}")
        return False

    # 2. 读取各染色体的all_simple文件
    methylation_change_data = {}
    total_read_lines = 0

    # (这里output_dir是output_x_y/甲基化类型/)
    mtype_dir = output_dir

    # 将该m*n*3次中的一次处理的所有all_simple文件的文件名获取到列表中
    all_simple_files = [f for f in os.listdir(mtype_dir) if f.startswith('all_simple_Chr') and f.endswith('.txt')]

    # all_simple文件格式为："Position", "Pvalue", "Methylation_Change"
    # 遍历所给一次处理的所有all_simple文件
    for file in all_simple_files:
        chr_match = re.search(r'Chr(\d+)\.txt$', file)
        if not chr_match:
            continue

        chr_num = int(chr_match.group(1))# 这里的chr_num也指的是both文件的第几个三列
        file_path = os.path.join(mtype_dir, file)

        try:
            with open(file_path, 'r') as f:  # 读取的all_simple文件格式为："Position", "Pvalue", "Methylation_Change"
                lines = f.readlines() # readlines用于返回文件的所有行组成的列表，每个元素是一行

            file_valid_lines = 0
            for line_num, line in enumerate(lines, 1): #遍历每一行，其中line_num从1开始枚举
                line = line.strip()
                if not line or line.startswith("Position"): # 跳过空行或首行
                    continue

                parts = line.split('\t')    #根据之前导出的时候设定的分隔符进行分割，将每次的三个数据项放到列表parts里
                if len(parts) >= 3:
                    # 由于读入的是字符串，所以得进行类型转换，将pvalue转换为浮点数，其余两个转换为整型
                    position = int(safe_float_convert(parts[0]))
                    pvalue = safe_float_convert(parts[1])
                    change = int(safe_float_convert(parts[2]))
                    meth_diff = safe_float_convert(parts[3]) if len(parts) >= 4 else 0.0

                    # 检查所有值是否有效
                    if (position is not None and
                            pvalue is not None and
                            change is not None):

                        # 检验change是否在0,1之中
                        if change in [0, 1]:
                            # 存储该映射: (chr, mtype, position) -> change 到字典中
                            methylation_change_data[(chr_num, mtype1, position)] = (change, meth_diff)
                            file_valid_lines += 1
                            total_read_lines += 1

            print(f"    {mtype_dir}/Chr{chr_num}: 读取 {file_valid_lines} 个有效位点")

        except Exception as e:
            print(f"    警告：读取 {file_path} 失败: {e}")
            continue

    print(f"    总共读取甲基化变化方向数据: {total_read_lines} 个位点")

    # 3. 合并数据，利用前面记录在字典中的(chr, mtype, position) -> change,将change这个属性添加到fdr_df的副本中
    combined_data = []
    missing_change = 0
    match_debug = defaultdict(int) #创建一个带默认值的字典，当访问不存在的键的时候，该字典会自动调用int类内置的调用函数
                                                # 为该键创建一个键值对，并将值置为int()的返回值——0

    for _, row in fdr_df.iterrows(): # fdr_df文件格式为 Chromosome, Methylation_Type, Position, Pvalue, Qvalue
        chr_num = row['Chromosome']
        mtype = row['Methylation_Type']
        # 统一转换为浮点数进行匹配
        position = safe_float_convert(row['Position'])
        qvalue = safe_float_convert(row['Qvalue'])

        # 若位点或q值缺失，则跳过当前行
        if position is None or qvalue is None:
            missing_change += 1
            continue

        # 查找对应的甲基化变化方向
        change_key = (chr_num, mtype, position)
        if change_key in methylation_change_data:
            change, meth_diff = methylation_change_data[change_key]  # 该change_data中存储该映射: (chr, mtype, position) -> change
            combined_data.append({
                'chromosome': chr_num,
                'methylation_type': mtype,
                'position': int(position),  # 输出时转为整数
                'qvalue': qvalue,
                'change': change,
                'meth_diff': meth_diff # 记录绝对甲基化差异，不过这个时候combined_data还是个列表而非dataframe
            })
            match_debug[chr_num] += 1 # 这里的chr_num也是both文件中第几个三列，从fisher检验读取Both文件就开始是这样了
        else:
            missing_change += 1

    print(f"    合并后有 {len(combined_data)} 个完整位点")
    print(f"    各染色体匹配情况: {dict(match_debug)}") # 其实就是每个染色体有多少条信息
    if missing_change > 0:
        print(f"    警告：{missing_change} 个位点缺少甲基化变化方向信息")

    # 4. 按染色体分组生成DMP文件
    chr_groups = defaultdict(list) # 这里的默认值字典在初始化时会调用list()函数，用[]空列表作为值
    for item in combined_data: # 每个item是一个字典，五个键值对的值为chr_num,mtype,position,qvalue,change
        chr_groups[item['chromosome']].append(item)
            # chr_groups这个字典将每个染色体的每一条数据（字典形式）作为一个元素记录到chr_groups->chr_num->list1这个列表中

    if not chr_groups:
        print(f"    错误：没有找到任何可生成DMP文件的数据")
        return False

    print(f"    将为以下染色体生成DMP文件: {sorted(chr_groups.keys())}")

    dir1_name = f"replicate{replicate_x}"
    dir2_name = f"replicate{replicate_y}"

    total_dmp = total_ndmp = total_hyper = total_hypo = 0

    # 为每个染色体生成文件
    for chr_num in sorted(chr_groups.keys()): # 这里的keys是染色体序号，序号-1就是一开始的全染色体映射里的索引，因为那个是从零开始算的
        chr_data = chr_groups[chr_num]  # chr_groups这个字典将每个染色体的每一条数据（字典形式）作为一个元素记录下来
                                    #所以这里的chr_data还是一个字典
        chr_data.sort(key=lambda x: x['position']) #lambda表达式匿名函数应用到列表的每个元素上，获取返回值作为排序依据
                                #  这里的匿名函数相当于对于chr_data中的每个元素（是一个字典）应用该函数
                                                #   def get_position(x):
                                            #     return x['position'] 返回值是位点号，所以就是根据位点号排序

        # 生成文件名（这里output是甲基化类型目录）
        dmp_file = os.path.join(output_dir, f"DMP_wt_{dir2_name}_mut_{dir1_name}_Chr{chr_num}.txt")
        ndmp_file = os.path.join(output_dir, f"N-DMP_wt_{dir2_name}_mut_{dir1_name}_Chr{chr_num}.txt")
        hyper_file = os.path.join(output_dir, f"hyper_DMP_wt_{dir2_name}_mut_{dir1_name}_Chr{chr_num}.txt")
        hypo_file = os.path.join(output_dir, f"hypo_DMP_wt_{dir2_name}_mut_{dir1_name}_Chr{chr_num}.txt")

        # 分类数据
        dmp_data = []
        ndmp_data = []
        hyper_data = []
        hypo_data = []

        # 遍历chr_data这个字典
        for item in chr_data: # chr_data是一个字典，五个键值对的值为chr_num,mtype,position,qvalue,change
            position = item['position']
            qvalue = item['qvalue']
            change = item['change']

            if qvalue <= fdr_threshold and item.get('meth_diff', 0.0) >= meth_diff_threshold: #判断是否显著，同时满足甲基化差异阈值
                dmp_data.append((position, qvalue, change))
                if change == 1:
                    hyper_data.append((position, qvalue, change))
                elif change == 0:
                    hypo_data.append((position, qvalue, change))
            else:
                ndmp_data.append((position, qvalue, change))

        if mtype1 not in unfilter_mtypes and all_dfs_ndmp_dict and mtype1 in all_dfs_ndmp_dict:
            if chr_num in all_dfs_ndmp_dict[mtype1]:
                ndmp_df = all_dfs_ndmp_dict[mtype1][chr_num]
                if isinstance(ndmp_df, pd.DataFrame) and not ndmp_df.empty:
                    # 获取已存在的位点集合，避免重复添加
                    existing_positions = set(item[0] for item in ndmp_data)
                    added_count = 0
                    for _, row in ndmp_df.iterrows():
                        pos = int(row['Position'])
                        pval = float(row['Pvalue'])
                        chg = int(row['Methylation_Change'])
                        # 只添加不重复的位点
                        if pos not in existing_positions:
                            ndmp_data.append((pos, pval, chg))
                            added_count += 1
                    if added_count > 0:
                        print(f"      为 Chr{chr_num} 添加了 {added_count} 个pvalue>0.05的NDMP位点")


        # 写入文件
        def write_dmp_file(filename, data):
            with open(filename, 'w') as f:
                f.write("first line\n")
                for pos, qval, chg in data:
                    f.write(f"{pos} {qval} {chg}\n")

        write_dmp_file(dmp_file, dmp_data)
        write_dmp_file(ndmp_file, ndmp_data)
        write_dmp_file(hyper_file, hyper_data)
        write_dmp_file(hypo_file, hypo_data)

        # print(f"    跳过 Chr{chr_num} 的单 pair DMR 候选区域识别（已移除分支 DMR 输出；final DMR 将基于 common DMP 重新计算）")

        print(f"    Chr{chr_num}: DMP={len(dmp_data)}, N-DMP={len(ndmp_data)}, Hyper={len(hyper_data)}, Hypo={len(hypo_data)}")

        # generate_dmp_files(output_dir, replicate_x, replicate_y, fdr_threshold=0.05, mtype1="CpG",
        #                    all_dfs_ndmp_dict=None, unfilter_mtypes=["CpG"]):
        # 统计
        total_dmp += len(dmp_data)
        total_ndmp += len(ndmp_data)
        total_hyper += len(hyper_data)
        total_hypo += len(hypo_data)

    print(f"    DMP文件生成完成！")
    print(f"    总计: DMP={total_dmp}, N-DMP={total_ndmp}, Hyper={total_hyper}, Hypo={total_hypo}")

    bothfile1 = os.path.join(dir1,f"{replicate_x}-bothMeUnme_diffChromo_NOREPEATED_methy_sites_{mtype1}.txt")

    bothfile2 = os.path.join(dir2,f"{replicate_y}-bothMeUnme_diffChromo_NOREPEATED_methy_sites_{mtype1}.txt")

    # print(f"    跳过 {mtype1} 的单 pair DMR reads 汇总/Fisher/FDR（已移除分支 DMR 输出；final DMR 将在 and_output/dmr_analysis_* 中重新计算）")

    if skip_window:
        print(f"    跳过 {mtype1} 的单次 DMP/N-DMP 滑动窗口分析 (--skip-window)")
    else:
        # 根据甲基化类型决定使用哪个滑动窗口函数
        if mtype1 not in unfilter_mtypes:
            print(f"使用新版本的滑动窗口分析")
            # 这里本来也要区分甲基化类型
            perform_sliding_window_on_dmp_files_after_filter(
                output_dir, replicate_x, replicate_y,
                all_dfs_ndmp_dict=all_dfs_ndmp_dict,
                methylation_type=mtype1
            )
            # perform_sliding_window_on_dmp_files(
            #     output_dir, replicate_x, replicate_y
            # )
        else:
            print(f"    使用标准版本的滑动窗口分析")
            perform_sliding_window_on_dmp_files(output_dir, replicate_x, replicate_y)
    return True



def regenerate_pair_dmp_outputs_from_fdr(
        m: int,
        n: int,
        dir1: str,
        dir2: str,
        work_dir: str = ".",
        meth_diff_threshold: float = 0.0,
        contexts=("CpG", "CHH", "CHG"),
):
    """Regenerate pairwise DMP/N-DMP/hyper/hypo files after auto-methdiff.

    Pairwise DMP files are initially produced before the auto-methdiff threshold
    is known.  When auto-methdiff changes the operational MethDiff threshold,
    those files would otherwise remain stale even though auto-vote and final
    common-DMP calling use the new threshold from the complete FDR tables.

    This function rebuilds the user-facing pairwise partition directly from:
      * complete FDR_corrected_results tables;
      * each row's Qvalue_Threshold_Used (or the context fallback);
      * the final MethDiff threshold;
      * raw WT/MUT methylation rates for the hyper/hypo direction.

    It also writes a validation summary.  A mismatch between the FDR universe
    and the raw WT/MUT intersection is treated as an error rather than being
    silently dropped.
    """
    work_path = Path(work_dir)
    meth_diff_threshold = float(meth_diff_threshold)
    summary_records = []

    wt_dir = Path(dir2)
    mut_dir = Path(dir1)
    wt_group_name = wt_dir.name
    mut_group_name = mut_dir.name

    print(
        "\n重新生成 pairwise DMP 输出："
        f"q通过 且 abs(MethDiff)>={meth_diff_threshold:.6g}"
    )

    for mut_idx in range(1, int(m) + 1):
        for wt_idx in range(1, int(n) + 1):
            wt_file = _mianjifa_find_replicate_file(
                wt_dir, wt_group_name, wt_idx
            )
            mut_file = _mianjifa_find_replicate_file(
                mut_dir, mut_group_name, mut_idx
            )

            for ctx in contexts:
                ctx_dir = (
                    work_path /
                    f"output_wt{wt_idx}_mut{mut_idx}" /
                    ctx
                )
                if not ctx_dir.is_dir():
                    continue

                fdr_file = ctx_dir / (
                    "FDR_corrected_results_"
                    f"wt_replicate{wt_idx}_vs_mut_replicate{mut_idx}.txt"
                )
                if not fdr_file.exists():
                    print(f"  ⚠ 跳过，FDR文件不存在：{fdr_file}")
                    continue

                fdr = pd.read_csv(fdr_file, sep="\t")
                required = {
                    "Chromosome", "Methylation_Type", "Position",
                    "Qvalue", "MethDiff",
                }
                missing = required - set(fdr.columns)
                if missing:
                    raise ValueError(
                        f"{fdr_file} 缺少必要列：{sorted(missing)}"
                    )

                frame = pd.DataFrame({
                    "chr": fdr["Chromosome"].map(
                        _mianjifa_normalize_chr
                    ),
                    "position": pd.to_numeric(
                        fdr["Position"], errors="coerce"
                    ),
                    "qvalue": pd.to_numeric(
                        fdr["Qvalue"], errors="coerce"
                    ),
                    "abs_methdiff_fdr": pd.to_numeric(
                        fdr["MethDiff"], errors="coerce"
                    ).abs(),
                })
                frame = frame.dropna(
                    subset=["position", "qvalue", "abs_methdiff_fdr"]
                ).copy()
                frame["position"] = frame["position"].astype(np.int64)

                fallback_q = float(get_dmp_threshold(ctx))
                if "Qvalue_Threshold_Used" in fdr.columns:
                    qthreshold = pd.to_numeric(
                        fdr.loc[frame.index, "Qvalue_Threshold_Used"],
                        errors="coerce",
                    ).fillna(fallback_q)
                    threshold_source = "Qvalue_Threshold_Used"
                else:
                    qthreshold = pd.Series(
                        fallback_q, index=frame.index, dtype=float
                    )
                    threshold_source = "fixed_context_threshold"
                frame["qthreshold"] = qthreshold.to_numpy(dtype=float)

                wt_raw = _mianjifa_read_raw_methylation_file(wt_file, ctx)
                mut_raw = _mianjifa_read_raw_methylation_file(mut_file, ctx)
                wt_raw = wt_raw.rename(
                    columns={"pos": "position", "meth_rate": "wt_rate"}
                )[["chr", "position", "wt_rate"]]
                mut_raw = mut_raw.rename(
                    columns={"pos": "position", "meth_rate": "mut_rate"}
                )[["chr", "position", "mut_rate"]]

                merged = frame.merge(
                    wt_raw,
                    on=["chr", "position"],
                    how="left",
                    validate="one_to_one",
                ).merge(
                    mut_raw,
                    on=["chr", "position"],
                    how="left",
                    validate="one_to_one",
                )

                missing_raw = int(
                    merged[["wt_rate", "mut_rate"]].isna().any(axis=1).sum()
                )
                if missing_raw:
                    examples = merged.loc[
                        merged[["wt_rate", "mut_rate"]].isna().any(axis=1),
                        ["chr", "position"],
                    ].head().to_dict("records")
                    raise RuntimeError(
                        f"{fdr_file}: {missing_raw} 个FDR位点在原始WT/MUT中"
                        f"无法匹配，示例={examples}"
                    )

                merged["signed_methdiff_raw"] = (
                    merged["mut_rate"] - merged["wt_rate"]
                )
                merged["change"] = (
                    merged["signed_methdiff_raw"] > 0
                ).astype(int)
                merged["q_pass"] = (
                    merged["qvalue"] <= merged["qthreshold"]
                )
                merged["methdiff_pass"] = (
                    merged["abs_methdiff_fdr"] >= meth_diff_threshold
                )
                merged["pair_dmp"] = (
                    merged["q_pass"] & merged["methdiff_pass"]
                )

                # Remove the preliminary partition before writing the final one.
                for pattern in (
                    "DMP_wt_replicate*_mut_replicate*_Chr*.txt",
                    "N-DMP_wt_replicate*_mut_replicate*_Chr*.txt",
                    "hyper_DMP_wt_replicate*_mut_replicate*_Chr*.txt",
                    "hypo_DMP_wt_replicate*_mut_replicate*_Chr*.txt",
                ):
                    for old_file in ctx_dir.glob(pattern):
                        old_file.unlink()

                total_dmp = 0
                total_ndmp = 0
                total_hyper = 0
                total_hypo = 0

                def _write_partition(path: Path, part: pd.DataFrame):
                    with path.open("w", encoding="utf-8") as handle:
                        handle.write("first line\n")
                        for row in part.itertuples(index=False):
                            handle.write(
                                f"{int(row.position)} {float(row.qvalue):.12g} "
                                f"{int(row.change)}\n"
                            )

                for chrom, chrom_df in merged.groupby("chr", sort=True):
                    chrom_df = chrom_df.sort_values("position")
                    dmp = chrom_df[chrom_df["pair_dmp"]].copy()
                    ndmp = chrom_df[~chrom_df["pair_dmp"]].copy()
                    hyper = dmp[dmp["change"] == 1].copy()
                    hypo = dmp[dmp["change"] == 0].copy()

                    prefix = (
                        f"wt_replicate{wt_idx}_mut_replicate{mut_idx}_"
                        f"Chr{chrom}.txt"
                    )
                    _write_partition(ctx_dir / f"DMP_{prefix}", dmp)
                    _write_partition(ctx_dir / f"N-DMP_{prefix}", ndmp)
                    _write_partition(ctx_dir / f"hyper_DMP_{prefix}", hyper)
                    _write_partition(ctx_dir / f"hypo_DMP_{prefix}", hypo)

                    total_dmp += len(dmp)
                    total_ndmp += len(ndmp)
                    total_hyper += len(hyper)
                    total_hypo += len(hypo)

                if total_dmp + total_ndmp != len(merged):
                    raise RuntimeError(
                        f"wt{wt_idx}-mut{mut_idx} {ctx}: DMP+N-DMP "
                        f"({total_dmp}+{total_ndmp}) != FDR rows ({len(merged)})"
                    )

                summary_records.append({
                    "context": ctx,
                    "wt_replicate": wt_idx,
                    "mut_replicate": mut_idx,
                    "fdr_rows": int(len(merged)),
                    "q_pass_count": int(merged["q_pass"].sum()),
                    "methdiff_pass_count": int(
                        merged["methdiff_pass"].sum()
                    ),
                    "pair_dmp_count": int(total_dmp),
                    "pair_ndmp_count": int(total_ndmp),
                    "hyper_count": int(total_hyper),
                    "hypo_count": int(total_hypo),
                    "q_threshold_source": threshold_source,
                    "meth_diff_threshold_used": meth_diff_threshold,
                    "wt_file": wt_file.name,
                    "mut_file": mut_file.name,
                    "status": "ok",
                })
                print(
                    f"  ✓ wt{wt_idx}-mut{mut_idx} {ctx}: "
                    f"DMP={total_dmp}, N-DMP={total_ndmp}"
                )

    summary = pd.DataFrame(summary_records)
    summary_dir = work_path / "and_output" / "auto_methdiff_thresholds"
    summary_dir.mkdir(parents=True, exist_ok=True)
    summary_file = summary_dir / "pair_dmp_regeneration_summary.tsv"
    summary.to_csv(summary_file, sep="\t", index=False)
    print(f"pairwise DMP重建汇总表保存至: {summary_file}")
    return summary

# 将此函数集成到 process_replicate_pair 中
def process_replicate_pair(replicate_x, replicate_y, files1, files2, dir1, dir2, dir1_name, dir2_name, unfilter_mtypes, work_dir=".", meth_diff_threshold=0.0, skip_dmr=False, skip_window=False):
    """处理一对组的所有甲基化类型，是处理1*3次检验
    replicate_x,replicate_y为both文件的序号
    files为一个字典链路：序号->甲基化类型->对应文件名
    dir1,dir2为两个基因型数据的目录名
    dir1_name，dir2_name是所给两个dir1,dir2路径的最后一个部分"""

    print(f"\n  处理组对 (wt{replicate_y}, mut{replicate_x})...")

    # 创建该组输出目录
    output_dir = os.path.join(work_dir, f"output_wt{replicate_y}_mut{replicate_x}")
    os.makedirs(output_dir, exist_ok=True)

    methylation_types = ['CpG', 'CHH', 'CHG']

    all_dfs_ndmp_dict = {}

    pair_success_count = 0
    pair_total_tests = 0

    # 循环处理每种甲基化类型
    for methylation_type in methylation_types:
        if methylation_type not in all_dfs_ndmp_dict:
            all_dfs_ndmp_dict[methylation_type] = {}

        success_count = total_tests = 0
        # 这里 in files1[replicate_x] 是判断methylation_type是否在files1[replicate_x] 这个字典的键中存在，存在的话就说明
        # 对应的 i-both...methylation_type.txt文件存在
        if (methylation_type not in files1[replicate_x] or
                methylation_type not in files2[replicate_y]):
            print(f"    跳过甲基化类型 {methylation_type}：因为文件不存在")
            continue
        # 否则获取到当前需要处理的文件的相对路径
        file1_path = os.path.join(dir1, files1[replicate_x][methylation_type])
        file2_path = os.path.join(dir2, files2[replicate_y][methylation_type])

        # 获取两个文件中的染色体数量
        n_chromosomes_1 = get_column_count(file1_path)
        n_chromosomes_2 = get_column_count(file2_path)

        if n_chromosomes_1 is None or n_chromosomes_2 is None:
            print(f"    无法获取 {methylation_type} 的染色体数量")
            continue

        if n_chromosomes_1 != n_chromosomes_2:
            print(f"    {methylation_type} 染色体数量不一致：{n_chromosomes_1} vs {n_chromosomes_2}")
            continue

        # 到这里就说明要处理的两个both文件中都有染色体数据且总列数相同
        # 不过其实可以不用判断的，因为一开始newtoboth已经保证列数肯定相同
        n_chromosomes = n_chromosomes_1  # 获取总染色体数
        print(f"    处理甲基化类型 {methylation_type}，共 {n_chromosomes} 条染色体")

        # 处理当前x_y组的当前甲基化类型的文件的每条染色体，这里的chr_num是第几个染色体，而不是染色体号
        for chr_num in range(1, n_chromosomes + 1):
            if methylation_type in unfilter_mtypes:
                success = process_methylation_type_with_collection(
                    file1_path, file2_path, methylation_type, output_dir,
                    dir1_name, dir2_name, chr_num, meth_diff_threshold=meth_diff_threshold
                )  # 处理1次检验的一个染色体
            else:
                # 这里要区分甲基化类型
                all_df_ndmp = process_methylation_type_with_collection_pvfilter(
                    file1_path, file2_path, methylation_type, output_dir,
                    dir1_name, dir2_name, chr_num, meth_diff_threshold=meth_diff_threshold
                )   #返回此次的一个pvalue>0.05的df
                    #  all_results([pos, pvalue, change])
                if chr_num not in all_dfs_ndmp_dict[methylation_type]:
                    all_dfs_ndmp_dict[methylation_type][chr_num] = all_df_ndmp
                    if (isinstance(all_df_ndmp, pd.DataFrame)):
                        print(f"成功添加{methylation_type}-{chr_num}-特定df到字典中，此次df中有{len(all_df_ndmp)}行数据")
                success = isinstance(all_df_ndmp, pd.DataFrame)
            total_tests += 1
            if success:  # 上一步处理正常进行就会返回True，否则返回False
                success_count += 1

        pair_success_count += success_count
        pair_total_tests += total_tests

        # 合并FET结果并进行FDR校正
        #  其中 FET结果的格式如：pos, m1, u1, m2, u2, pvalue
        if success_count > 0:
            output_dir1 = os.path.join(output_dir, methylation_type)
            # 获取一个output_x_y/甲基化类型/目录下的所有FET文件，将其concat起来，然后FDR检验为qvalues,并导出到磁盘，最终格式为：
            #                                        Chromosome, Methylation_Type, Position, Pvalue, Qvalue
            merge_ok, dmp_threshold = merge_fet_results_and_fdr(
                output_dir1, replicate_x, replicate_y, methylation_type,
                all_dfs_ndmp_dict, n_chromosomes,
                is_twostep_context=(methylation_type not in unfilter_mtypes)
            )
            if not merge_ok:
                print(f"    {methylation_type} 的FDR合并失败，跳过DMP文件生成")
                continue
            # 生成DMP文件：固定阈值或两步法自动阈值由 merge_fet_results_and_fdr 返回
            generate_dmp_files(dir1,dir2,output_dir1, replicate_x, replicate_y, fdr_threshold=dmp_threshold, mtype1=methylation_type,all_dfs_ndmp_dict=all_dfs_ndmp_dict
                               ,unfilter_mtypes=unfilter_mtypes,n_chromosomes=n_chromosomes, meth_diff_threshold=meth_diff_threshold,
                               skip_dmr=skip_dmr, skip_window=skip_window)

    print(f"  组对 (wt{replicate_y}, mut{replicate_x}) 处理完成！成功 {pair_success_count}/{pair_total_tests} 次检验")
    return pair_success_count, pair_total_tests    # 这里是一个染色体就算一次检验


def process_all_combinations(dir1, dir2, m, n, unfilter_mtypes, work_dir=".", meth_diff_threshold=0.0,
                             skip_dmr=False, skip_window=False, threads=1):
    """处理所有组合，进行m*n*3次检验；支持 replicate-pair 级并行。"""

    print(f"扫描文件目录...")
    files1 = scan_sample_files_by_replicates(dir1, m)
    files2 = scan_sample_files_by_replicates(dir2, n)

    print(f"目录1 ({dir1}) 找到 {len(files1)} 组文件")
    print(f"目录2 ({dir2}) 找到 {len(files2)} 组文件")

    missing_replicates1 = [i for i in range(1, m + 1) if i not in files1]
    missing_replicates2 = [i for i in range(1, n + 1) if i not in files2]
    if missing_replicates1:
        print(f"警告：目录1缺少这些组: {missing_replicates1}")
    if missing_replicates2:
        print(f"警告：目录2缺少这些组: {missing_replicates2}")

    available_replicates1 = [i for i in range(1, m + 1) if i in files1]
    available_replicates2 = [i for i in range(1, n + 1) if i in files2]

    total_combinations = len(available_replicates1) * len(available_replicates2)
    print(f"\n开始处理 {total_combinations} 个组合...")

    dir1_name = sanitize_filename(os.path.basename(dir1.rstrip(os.sep)))
    dir2_name = sanitize_filename(os.path.basename(dir2.rstrip(os.sep)))

    total_success = total_tests = 0
    start_time = time.time()
    threads = max(1, int(threads))

    tasks = []
    config = _build_parallel_config()
    for replicate_x in available_replicates1:
        for replicate_y in available_replicates2:
            tasks.append({
                "replicate_x": replicate_x,
                "replicate_y": replicate_y,
                "files1": files1,
                "files2": files2,
                "dir1": dir1,
                "dir2": dir2,
                "dir1_name": dir1_name,
                "dir2_name": dir2_name,
                "unfilter_mtypes": unfilter_mtypes,
                "work_dir": work_dir,
                "meth_diff_threshold": meth_diff_threshold,
                "skip_dmr": skip_dmr,
                "skip_window": skip_window,
                "config": config,
            })

    if threads <= 1 or len(tasks) <= 1:
        for i, task in enumerate(tasks, 1):
            print(f"\n进度: {i}/{total_combinations}")
            success_count, test_count = process_replicate_pair(
                task["replicate_x"], task["replicate_y"], files1, files2,
                dir1, dir2, dir1_name, dir2_name, unfilter_mtypes,
                work_dir=work_dir,
                meth_diff_threshold=meth_diff_threshold,
                skip_dmr=skip_dmr,
                skip_window=skip_window
            )
            total_success += success_count
            total_tests += test_count
    else:
        max_workers = min(threads, len(tasks))
        log_dir = os.path.join(work_dir, "parallel_logs")
        os.makedirs(log_dir, exist_ok=True)
        for task in tasks:
            task["log_file"] = os.path.join(
                log_dir,
                f"pair_wt{task['replicate_y']}_mut{task['replicate_x']}.log"
            )

        print(f"启用并行 replicate-pair 处理: workers={max_workers}, tasks={len(tasks)}")
        with ProcessPoolExecutor(max_workers=max_workers) as executor:
            future_to_task = {executor.submit(_run_replicate_pair_worker, task): task for task in tasks}
            for i, future in enumerate(as_completed(future_to_task), 1):
                task = future_to_task[future]
                label = f"wt{task['replicate_y']}_mut{task['replicate_x']}"
                try:
                    res = future.result()
                    total_success += res["success_count"]
                    total_tests += res["test_count"]
                    print(
                        f"[DONE] {i}/{len(tasks)} {label}: 成功 "
                        f"{res['success_count']}/{res['test_count']} 次检验, "
                        f"{res['elapsed']:.2f}s, log={res['log_file']}"
                    )
                except Exception as e:
                    print(f"[FAILED] {i}/{len(tasks)} {label}: {e}")
                    raise

    end_time = time.time()
    print(f"\n所有处理完成！")
    print(f"总计: {total_success}/{total_tests} 次成功检验")
    print(f"用时: {end_time - start_time:.2f} 秒")
    print(f"单次比较结果保存在 ./output_x_y/甲基化类型/ 目录中")

    # 兼容原代码行为：若没有可计数测试但输出流程未抛异常，不判定为失败。
    return total_success == total_tests

def _dmp_required_count(total_count, auto_vote_threshold=None):
    """Return the integer final-DMP vote requirement used for a context."""
    total_count = int(total_count)
    if total_count <= 0:
        return 0
    if auto_vote_threshold is not None:
        return max(1, min(int(auto_vote_threshold), total_count))
    return max(1, min(int(np.floor(VOTE_THRESHOLD * total_count + 0.5)), total_count))


def _dmp_boundary_abs_methdiff(support_methdiffs, base_required):
    """Boundary abs(MethDiff): the base_required-th largest supporting abs(MethDiff)."""
    try:
        base_required = int(base_required)
    except Exception:
        return np.nan
    if base_required <= 0:
        return np.nan
    values = []
    for v in support_methdiffs or []:
        try:
            v = abs(float(v))
        except Exception:
            continue
        if np.isfinite(v):
            values.append(v)
    if len(values) < base_required:
        return np.nan
    values.sort(reverse=True)
    return float(values[base_required - 1])


def _write_lowdiff_summary_record(summary_file, record):
    """Update lowdiff summary TSV, replacing the current context row if present."""
    os.makedirs(os.path.dirname(summary_file), exist_ok=True)
    new_df = pd.DataFrame([record])
    if os.path.exists(summary_file):
        try:
            old_df = pd.read_csv(summary_file, sep="\t")
            if "context" in old_df.columns:
                old_df = old_df[old_df["context"].astype(str) != str(record.get("context"))]
            out_df = pd.concat([old_df, new_df], ignore_index=True)
        except Exception:
            out_df = new_df
    else:
        out_df = new_df
    out_df.to_csv(summary_file, sep="\t", index=False)


def apply_dmp_lowdiff_strict_vote_to_result_df(
        result_df,
        site_statistics,
        common_sites,
        methytype,
        and_output_dir,
        total_groups,
        auto_vote_threshold=None):
    """Apply/report final-DMP low-difference strict vote for one methylation context.

    The baseline/provisional DMP set has already passed the normal final vote.
    For each provisional DMP, compute boundary_abs_methdiff as the
    base_required-th largest abs(MethDiff) among supporting pairwise tests.
    If boundary_abs_methdiff <= DMP_LOWDIFF_CUTOFF, require a stricter vote:
        low_required = ceil((base_required + total_groups) / 2)

    In report-only mode, diagnostics are written but result_df is not filtered.
    """
    if result_df is None or result_df.empty:
        return result_df

    diag_dir = os.path.join(and_output_dir, "lowdiff_strict_vote")
    os.makedirs(diag_dir, exist_ok=True)

    total_groups = int(total_groups)
    base_required = _dmp_required_count(total_groups, auto_vote_threshold=auto_vote_threshold)
    low_required = int(math.ceil((base_required + total_groups) / 2.0))
    low_required = max(base_required, min(low_required, total_groups))
    cutoff = float(DMP_LOWDIFF_CUTOFF)

    print(
        f"\nDMP lowdiff strict vote ({methytype}): "
        f"enabled={DMP_LOWDIFF_STRICT_VOTE}, "
        f"report_only={DMP_LOWDIFF_STRICT_VOTE_REPORT_ONLY}, "
        f"cutoff={cutoff:.6g}, base_required={base_required}, "
        f"low_required={low_required}, total_groups={total_groups}"
    )
    print(
        "  规则：provisional final DMP 中，若 boundary_abs_methdiff <= cutoff "
        "且 support_count < low_required，则从最终DMP中移除。"
    )

    records = []
    remove_site_ids = set()
    for site_id in common_sites:
        stats = site_statistics.get(site_id, {})
        support_count = int(stats.get("sig_count", 0))
        total_count = int(stats.get("total_count", total_groups))
        support_methdiffs = stats.get("support_methdiffs", [])
        boundary = _dmp_boundary_abs_methdiff(support_methdiffs, base_required)

        clean_support = []
        for v in support_methdiffs:
            try:
                v = abs(float(v))
            except Exception:
                continue
            if np.isfinite(v):
                clean_support.append(v)

        max_support = float(np.max(clean_support)) if clean_support else np.nan
        min_support = float(np.min(clean_support)) if clean_support else np.nan
        mean_support = float(np.mean(clean_support)) if clean_support else np.nan
        is_lowdiff = bool(np.isfinite(boundary) and boundary <= cutoff)
        removed = bool(is_lowdiff and support_count < low_required)
        if removed:
            remove_site_ids.add(site_id)

        chr_num = stats.get("chromosome", np.nan)
        position = stats.get("position", np.nan)
        records.append({
            "Chromosome": chr_num,
            "Methylation_Type": methytype,
            "Position": position,
            "site_id": site_id,
            "support_count": support_count,
            "total_count": total_count,
            "base_required": base_required,
            "low_required": low_required,
            "lowdiff_cutoff": cutoff,
            "boundary_abs_methdiff": boundary,
            "max_support_abs_methdiff": max_support,
            "min_support_abs_methdiff": min_support,
            "mean_support_abs_methdiff": mean_support,
            "is_lowdiff_candidate": is_lowdiff,
            "passed_base_vote": True,
            "passed_lowdiff_strict_vote": not removed,
            "removed_by_lowdiff_strict_vote": removed,
            "report_only": bool(DMP_LOWDIFF_STRICT_VOTE_REPORT_ONLY),
        })

    diag_df = pd.DataFrame(records)
    diag_df = diag_df.sort_values(["Chromosome", "Position"], na_position="last")
    diag_file = os.path.join(diag_dir, f"{methytype}-lowdiff_strict_vote_diagnostics.tsv")
    diag_df.to_csv(diag_file, sep="\t", index=False)

    provisional_file = os.path.join(diag_dir, f"{methytype}-provisional_final_DMPs_before_lowdiff.tsv")
    result_df.to_csv(provisional_file, sep="\t", index=False)

    n_provisional = int(len(result_df))
    n_lowdiff = int(diag_df["is_lowdiff_candidate"].sum()) if not diag_df.empty else 0
    n_removed = int(diag_df["removed_by_lowdiff_strict_vote"].sum()) if not diag_df.empty else 0

    if DMP_LOWDIFF_STRICT_VOTE and not DMP_LOWDIFF_STRICT_VOTE_REPORT_ONLY:
        filtered_df = result_df[~result_df["site_id"].isin(remove_site_ids)].copy()
    else:
        filtered_df = result_df.copy()

    summary_record = {
        "context": methytype,
        "total_groups": total_groups,
        "base_required": base_required,
        "low_required": low_required,
        "lowdiff_cutoff": cutoff,
        "n_provisional_final_dmp": n_provisional,
        "n_lowdiff_candidates": n_lowdiff,
        "n_removed_by_lowdiff_strict_vote": n_removed,
        "n_final_after_lowdiff": int(len(filtered_df)),
        "enabled": bool(DMP_LOWDIFF_STRICT_VOTE),
        "report_only": bool(DMP_LOWDIFF_STRICT_VOTE_REPORT_ONLY),
        "diagnostic_file": diag_file,
        "provisional_file": provisional_file,
    }
    summary_file = os.path.join(diag_dir, "lowdiff_strict_vote_summary.tsv")
    _write_lowdiff_summary_record(summary_file, summary_record)

    print(f"  provisional final DMP: {n_provisional}")
    print(f"  lowdiff candidates: {n_lowdiff}")
    print(f"  removed by lowdiff strict vote: {n_removed}")
    print(f"  final after lowdiff strict vote: {len(filtered_df)}")
    print(f"  lowdiff诊断表保存至: {diag_file}")
    print(f"  provisional DMP备份保存至: {provisional_file}")
    print(f"  lowdiff汇总表保存至: {summary_file}")

    return filtered_df


def bayes_deciding(sig_count, nonsig_count, auto_vote_threshold=None):
    """根据跨replicate组合中的显著支持次数进行最终DMP/DMR投票判定。

    若 auto_vote_threshold 为整数，则直接使用该 required_count；
    否则沿用 --vote-threshold 对应的 round-half-up required_count。
    """
    total_count = sig_count + nonsig_count
    if total_count <= 0:
        return 0

    if auto_vote_threshold is not None:
        required_count = int(auto_vote_threshold)
    else:
        required_count = int(np.floor(VOTE_THRESHOLD * total_count + 0.5))
    final_decision = 1 if sig_count >= required_count else 0
    # print(f"\n判定（投票阈值={VOTE_THRESHOLD * 100:.1f}%）")
    # print(f"  判定结果：{'显著' if final_decision else '不显著'}")
    # print(f"  支持比例：{support_ratio * 100:.1f}%")

    return final_decision

def find_common_significant_sites(output_dirs=None, methytype2='CpG', dir1=None, dir2=None, work_dir=".", meth_diff_threshold=0.0):
    """
    找出在所有组合检验中都显著的位点，并获取对应位点的相关信息
    参数：
        output_dirs: 输出目录列表，如果为None则自动扫描
        methytype2: 甲基化类型
    """

    print("\n寻找所有组合中共同显著的位点...")

    and_output_dir = os.path.join(work_dir, "and_output")
    os.makedirs(and_output_dir, exist_ok=True)

    # 1. 自动扫描所有output_x_y目录
    if output_dirs is None:
        output_dirs = glob.glob(os.path.join(work_dir,f"output_*_*/{methytype2}"))
        # 找出所有output_x_y/methytype2/这种格式的目录，一共m*n个
        output_dirs = [d for d in output_dirs if os.path.isdir(d)]

    if not output_dirs: # 没找到的话
        print("未找到任何输出目录")
        return None

    print(f"找到 {len(output_dirs)} 个输出目录")

    auto_vote_threshold = AUTO_DMP_VOTE_THRESHOLDS.get(methytype2)
    if AUTO_VOTE_THRESHOLD_REPORT_ONLY:
        auto_vote_threshold = None
    if auto_vote_threshold is not None:
        print(f"使用 {methytype2} DMP自动投票阈值: {auto_vote_threshold}")
    else:
        print(f"{methytype2} DMP自动投票阈值不可用或仅报告，使用比例阈值 {VOTE_THRESHOLD}")

    # 2. 一次性读取所有FDR_correct文件到内存，读取所有FDR_corrected文件（包含所有位点）,因为如果只读取显著位点信息的话，
            # 假设某个位点在1_1 sig文件中有出现，但在后面一次检验中没出现，就不知道是因为不显著没出现还是这次输入的原始数据里就没有该位点
    valid_dirs = [] # 收集在后续操作中起到作用的目录的路径
    site_statistics = {}  # 建立如此映射：{site_id: {'sig': 0, 'total': 0}}
    all_dataframes = {} # 最终可以通过all_dataframes[目录]->FDR_correct对应的df
    dir_to_replicate = {} #记录目录对应的replicate编号
    for output_dir in output_dirs: # 遍历所有output_x_y/methytype2/
        # 获取到当前甲基化目录下的FDR_correct文件，
            # 格式为：'Chromosome', 'Methylation_Type', 'Position', 'Pvalue', 'Qvalue'
        # 其中qvalue不同于之前都是小于阈值的，而是所有的
        match = re.search(r'output_wt(\d+)_mut(\d+)', output_dir)
        if not match:
            continue

        replicate_y = int(match.group(1))
        replicate_x = int(match.group(2))


        fdr_all_files = glob.glob(os.path.join(output_dir, "FDR_corrected_results_*.txt"))

        # 检验存在性
        if not fdr_all_files:
            print(f"  警告：{output_dir} 中未找到FDR_corrected文件")
            continue

        fdr_all_file = fdr_all_files[0] # 这是因为glob.glob的返回值是一个列表，所以用[0]获取到实际存在的文件路径

        try:
            df = pd.read_csv(fdr_all_file, sep=r'\s+') # 读取该文件到df中FDR_corrected格式为：
                                        # 'Chromosome', 'Methylation_Type', 'Position', 'Pvalue', 'Qvalue'
            if len(df) > 0:
                # 为每一行创建唯一标识符
                #  具体为：对df中的每一行，将其中的染色体号-甲基化类型-位点号提取出来，放在右侧作为新的一列site_id
                #     这是因为后面需要统计所有位点在多次检验中的显著与否，通过site_id就能仅用一个属性区分不同的位点了，不必去
                #   连续判断多个属性是否相等来筛选或访问（染色体号、甲基化类型、位点号） （注意这里的染色体号是both文件中的第几个三列）
                df['site_id'] = df.apply(
                    lambda row: f"{int(row['Chromosome'])}-{row['Methylation_Type']}-{int(row['Position'])}",
                    axis=1
                )
                all_dataframes[output_dir] = df # all_dataframes[目录]->FDR_correct对应的df
                valid_dirs.append(output_dir)
                dir_to_replicate[output_dir] = (replicate_x, replicate_y)  # 记录编号
                # 获取当前甲基化类型的默认DMP阈值；若FDR文件含有每次比较实际使用的阈值，则优先使用该列
                dmp_threshold = get_dmp_threshold(methytype2)
                for _, row in df.iterrows(): # 遍历每一行，每一行是一次检验中FDR_correct中所有不管显著还是不显著的位点信息
                    site_id = row['site_id']

                    if site_id not in site_statistics:
                        site_statistics[site_id] = {
                            'sig_count': 0,
                            'total_count': 0,
                            'chromosome': int(row['Chromosome']),
                            'methylation_type': row['Methylation_Type'],
                            'position': int(row['Position']),
                            'support_methdiffs': [],
                            'support_qvalues': [],
                        }

                    site_statistics[site_id]['total_count'] += 1

                    threshold_used = row.get('Qvalue_Threshold_Used', dmp_threshold)
                    try:
                        threshold_used = float(threshold_used)
                    except Exception:
                        threshold_used = dmp_threshold

                    try:
                        qvalue_for_support = float(row['Qvalue'])
                    except Exception:
                        qvalue_for_support = np.nan
                    try:
                        methdiff_for_support = abs(float(row.get('MethDiff', 0.0)))
                    except Exception:
                        methdiff_for_support = np.nan

                    support_pass = (
                        np.isfinite(qvalue_for_support)
                        and qvalue_for_support <= threshold_used
                        and np.isfinite(methdiff_for_support)
                        and methdiff_for_support >= meth_diff_threshold
                    )
                    if support_pass:
                        site_statistics[site_id]['sig_count'] += 1
                        site_statistics[site_id]['support_methdiffs'].append(methdiff_for_support)
                        site_statistics[site_id]['support_qvalues'].append(qvalue_for_support)


                print(f"  {output_dir}: {len(df)} 个位点")
            else:
                print(f"  {output_dir}: 无显著位点")
        except Exception as e:
            print(f"  错误：读取 {fdr_all_file} 失败: {e}")
            continue
    valid_dirs.sort(key=lambda d: dir_to_replicate[d])
    print(f"\n统计完成，共 {len(site_statistics)} 个不同位点")

    if not site_statistics:
        print("没有找到任何有效的位点信息结果")
        return None

    # 3. 获取到贝叶斯方法判定为显著的位点，放到common_sites里
    common_sites = []
    for site_id, stats in site_statistics.items():
        if stats['total_count'] != len(valid_dirs):
            stats['total_count'] = len(valid_dirs)
        sig_count = stats['sig_count']  # 获取到当前位点显著检验次数
        nonsig_count = stats['total_count'] - sig_count  # 获取到当前位点非显著检验次数
        is_significant = bayes_deciding(sig_count, nonsig_count, auto_vote_threshold=auto_vote_threshold)
        if is_significant:
            common_sites.append(site_id)
    if not common_sites:
        print("没有在所有组合中都显著的位点")
        return None
    else:
        print(f"\n贝叶斯判定后，共 {len(common_sites)} 个显著位点")

    # 4. 读取每个目录的甲基化变化方向信息
    print("正在读取甲基化变化方向信息...")
    methylation_change_by_dir = {}  # {output_dir: {site_id: change}}

    for output_dir in valid_dirs: # 遍历上述有效的该甲基化类型的目录
        methylation_change_by_dir[output_dir] = {} # 创建当前目录的字典元素，
                                        # 构建 output_dir->site_id->change的链路
        # 查找该目录下所有的all_simple_Chr文件，因为这个文件中有change信息，其格式为pos, pvalue, change，
                                                # 而染色体号从文件名获取
        all_simple_files = glob.glob(os.path.join(output_dir, "all_simple_Chr*.txt"))
                                                        # 获取当前目录下所有的all_simple文件的路径并组成列表

        for file_path in all_simple_files: # 遍历all_simple文件，格式为pos, pvalue, change
            chr_match = re.search(r'Chr(\d+)\.txt$', file_path) #创建正则表达式及捕获组，用于获取染色体号
            if not chr_match:
                continue
            chr_num = int(chr_match.group(1)) # 通过捕获组获取到当前文件的染色体号

            try:
                df = pd.read_csv(
                    file_path,
                    sep='\t',
                    usecols=[0, 2],
                    names=['position', 'change'],
                    dtype={'position': float, 'change': float}, # float啥字符串都可以转换，不会发生错误，是比较稳妥的方案
                    skiprows=1
                )

                # 此时将pos和change都转换为int类型就不会出问题，否则str转int有可能出问题，比如"100.0"转int就会报错
                df['position'] = df['position'].astype(int)
                df['change'] = df['change'].astype(int)

                positions = df['position'].values
                changes = df['change'].values

                for position, change in zip(positions, changes):
                    site_id = f"{chr_num}-{methytype2}-{position}"
                    methylation_change_by_dir[output_dir][site_id] = change

            except Exception as e:
                print(f"  警告：读取 {file_path} 失败: {e}")
                continue

        print(f"  {output_dir}: 读取了 {len(methylation_change_by_dir[output_dir])} 个位点的变化方向")

    # 5. 处理共同位点信息
    print("正在处理共同位点详细信息...")
    common_site_details = []

    # 为每个DataFrame创建索引以加快查询
    indexed_dfs = {} # 建立output_dir->df1（site_id变作索引后）的链路
    for output_dir in valid_dirs:
        df = all_dataframes[output_dir] # 获取每个目录中的FDR_correct对应的df（已经加上了site_id这列）
        indexed_dfs[output_dir] = df.set_index('site_id') # 将site_id设置为索引，并将新的df1返回作为字典的值

    # 逐一处理common位点
    for i, site_id in enumerate(common_sites):
        if i % 10000 == 0:  # 每处理10000个位点打印进度
            print(f"  已处理 {i}/{len(common_sites)} 个位点")

        site_info = {'site_id': site_id}
        chr_num, mtype, pos = site_id.split('-')
        site_info['Chromosome'] = int(chr_num)
        site_info['Methylation_Type'] = mtype
        site_info['Position'] = int(pos)

        # 收集q值的列表，收集不同output_x_y/methytype2/目录下相同染色体相同位点号(符合当前site_id信息的那个位点)在所有检验中的q值
        qvalues = []
        sig_qvalues_for_mean = []
        # 收集不同output_x_y/methytype2/目录下相同染色体相同位点号(符合当前site_id信息的那个位点)在所有检验中的变化方向
        change_values = []
        qvalue_dict = {}
        for output_dir in valid_dirs:
            replicate_x, replicate_y = dir_to_replicate[output_dir]
            col_name = f'qvalue_{os.path.basename(dir2.rstrip("/"))}{replicate_y}_{os.path.basename(dir1.rstrip("/"))}{replicate_x}'   # 新增：列名，前面是野生型的序号，后面是突变型的序号
            indexed_df = indexed_dfs[output_dir] # output_dir->df1（site_id变作索引后）,获取到索引为site_id的当前甲基化的某个FDR_correct文件
                                # 内容格式为：'Chromosome', 'Methylation_Type', 'Position', 'Pvalue', 'Qvalue'
            if site_id in indexed_df.index:
                row_for_site = indexed_df.loc[site_id]
                if isinstance(row_for_site, pd.DataFrame):
                    row_for_site = row_for_site.iloc[0]
                qval = row_for_site['Qvalue']
                qvalues.append(qval)
                qvalue_dict[col_name] = qval

                threshold_used = row_for_site.get('Qvalue_Threshold_Used', get_dmp_threshold(methytype2))
                try:
                    threshold_used = float(threshold_used)
                except Exception:
                    threshold_used = get_dmp_threshold(methytype2)
                if float(qval) <= threshold_used:
                    sig_qvalues_for_mean.append(float(qval))

                # 获取该位点在该目录中的change值，(output_dir->site_id->change的链路)
                if site_id in methylation_change_by_dir[output_dir]:
                    change_values.append(methylation_change_by_dir[output_dir][site_id])
            else:
                qvalue_dict[col_name] = 1.0

        if qvalues:  # 确保有q值数据
            if len(sig_qvalues_for_mean) > 0:
                site_info['Sig_Mean_Qvalue'] = np.mean(sig_qvalues_for_mean)
            else:
                site_info['Sig_Mean_Qvalue'] = 1
            # site_info['Max_Qvalue'] = np.max(qvalues)
            # site_info['Min_Qvalue'] = np.min(qvalues)
            site_info['Num_Comparisons'] = len(qvalues)

            # 投票计算甲基化变化方向
            if change_values:
                # 统计change==1的次数
                num_hyper = sum(change_values)
                total_comparisons = len(change_values)
                hyper_ratio = num_hyper / total_comparisons

                # 多数投票：>= 50%则记为1（高甲基化），否则为0
                site_info['Methylation_Change'] = 1 if hyper_ratio >= 0.5 else 0
                site_info['Hyper_Count'] = num_hyper  # 高甲基化次数
                site_info['Hypo_Count'] = total_comparisons - num_hyper  # 低甲基化次数
                site_info['Hyper_Ratio'] = hyper_ratio  # 高甲基化比例
            else:
                # 如果没有change信息，标记为缺失（不过按理来说一次检验是正好对应一个qvalue和一个change的）
                site_info['Methylation_Change'] = -1  # -1表示无法确定
                site_info['Hyper_Count'] = 0
                site_info['Hypo_Count'] = 0
                site_info['Hyper_Ratio'] = 0

            site_info.update(qvalue_dict) # 添加所有replicate的qvalue列

            common_site_details.append(site_info) # 这个列表中记录了一个个字典，每个字典是一个site_id对应的各个信息，格式：
        # site_id-染色体号-甲基化类型-位点号-总检验次数-change-hypercount-hypocount-hyperratio-Sig_Mean_Qvalue-所有qvalues
            # （此处染色体号是both文件中的第几个3列）

    print(f"  完成处理 {len(common_site_details)} 个位点")

    # 6. 生成结果DataFrame
    result_df = pd.DataFrame(common_site_details)
    result_df = result_df.sort_values(['Methylation_Type', 'Chromosome', 'Position']) # 排序

    # 可选：final DMP lowdiff strict vote 后处理。
    # 这里按甲基化类型分别执行；cutoff可以全局固定为0.3，但support_count、base_required、
    # boundary_abs_methdiff、诊断表都必须按CpG/CHG/CHH分别计算。
    if DMP_LOWDIFF_STRICT_VOTE or DMP_LOWDIFF_STRICT_VOTE_REPORT_ONLY:
        result_df = apply_dmp_lowdiff_strict_vote_to_result_df(
            result_df=result_df,
            site_statistics=site_statistics,
            common_sites=common_sites,
            methytype=methytype2,
            and_output_dir=and_output_dir,
            total_groups=len(valid_dirs),
            auto_vote_threshold=auto_vote_threshold,
        )

    # 调整列顺序，将Methylation_Change放在更靠前的位置
    column_order = [
        'Chromosome', 'Methylation_Type', 'Position',
        'Methylation_Change', 'Hyper_Ratio', 'Hyper_Count', 'Hypo_Count', 'Num_Comparisons',
        'Sig_Mean_Qvalue'
    ]
    # 获取所有replicate列名并排序
    # 新的列名格式：dir2_y_dir1_x，需要获取包含下划线分隔数字的列
    replicate_columns = sorted(
        [col for col in result_df.columns if col.startswith('qvalue_')],
        key=lambda x: tuple(map(int, re.findall(r'\d+', x)))  # 提取所有数字并排序
    )
    column_order = column_order + replicate_columns
    result_df = result_df[column_order]

    # 保存结果
    output_file = os.path.join(and_output_dir, f"{methytype2}-final_significant_sites_DMPs.txt")
    result_df.to_csv(output_file, sep='\t', index=False)
    print(f"\n共同显著位点已保存至: {output_file}")

    # 输出统计信息
    print("\n共同显著位点统计:")
    mtype = methytype2
    mtype_df = result_df
    count = len(mtype_df)
    if count == 0:
        print(f"  {mtype}: 0 个位点")
        return result_df
    hyper_count = len(mtype_df[mtype_df['Methylation_Change'] == 1])
    hypo_count = len(mtype_df[mtype_df['Methylation_Change'] == 0])
    unknown_count = len(mtype_df[mtype_df['Methylation_Change'] == -1])

    print(f"  {mtype}: {count} 个位点")
    print(f"    - 高甲基化(Change=1): {hyper_count} ({hyper_count / count * 100:.1f}%)")
    print(f"    - 低甲基化(Change=0): {hypo_count} ({hypo_count / count * 100:.1f}%)")
    if unknown_count > 0:
        print(f"    - 未知(Change=-1): {unknown_count} ({unknown_count / count * 100:.1f}%)")

    return result_df # 其格式为：'Chromosome', 'Methylation_Type', 'Position',
                            # 'Methylation_Change', 'Hyper_Ratio', 'Hyper_Count', 'Hypo_Count',
                           # 'Sig_Mean_Qvalue', 'Max_Qvalue', 'Min_Qvalue', 'Num_Comparisons'

def sliding_window_analysis(
        input_data,
        window_size = 1000000,
        step_ratio = 0.05,
        save_files = False,
        output_identifier = None,
        outputdir1="./"
):
    """
    对甲基化位点数据进行滑动窗口分析
    参数：
    input_data : 是m*n*3次中一次的某一条染色体的(N)dmp文件转换为的DataFrame，包含列: ['position', 'pvalue', 'change']
    window_size :滑动窗口的大小
    step_ratio : 步长比例（窗口大小的百分比）
    save_files : 是否保存结果文件
    output_identifier : 输出文件前缀，如果save_files=True且未提供则自动生成
    """

    # 1. 数据读取和预处理
    df = input_data.copy()

    # 确保列名正确
    expected_cols = ['position', 'pvalue', 'change']
    if not all(col in df.columns for col in expected_cols):
        raise ValueError(f"DataFrame必须包含列: {expected_cols}")

    if output_identifier is None:
        output_identifier = "sliding_window_analysis"

    # 数据验证和清理
    df = df.dropna()
    if len(df) == 0:
        raise ValueError("没有有效的数据行")

    # 按位点号排序
    df = df.sort_values('position').reset_index(drop=True) #重排dataframe索引，且不将原来的索引变作新的一列

    print(f"数据预处理完成，共 {len(df)} 个位点")

    # 2. 滑动窗口分析
    positions = df['position'].values #获取所有位点号（前面已经排了序了）
    changes = df['change'].values #获取所有甲基化变化方向

    last_pos = positions[-1] #获取到最后一个位点号
    step_size = int(window_size * step_ratio) #计算得到窗口每一步位移的大小

    # 生成所有窗口的起始位置
    window_starts = np.arange(0, last_pos + step_size, step_size) #中间+step_size是为了使得最后一个区间覆盖last_pos这个位点

    print(f"窗口配置: 大小={window_size}, 步长={step_size}, 窗口数={len(window_starts)}")

    #创建列表用于存储后续每个区间的首尾位点号和不同甲基化变化方向的数量以及总显著位点数
    results = []

    for i, start in enumerate(window_starts): # start每次赋值为每个窗口的起始点位点号-1（比如第一次是0）
        if i % 1000 == 0:  # 进度提示
            print(f"处理进度: {i}/{len(window_starts)}")

        end = start + window_size  # 根据窗口宽度计算得到窗口的终止点

        # 使用numpy的searchsorted进行快速查找（相当于二分查找）
        left_idx = np.searchsorted(positions, start, side='right') # 此处意为查找positions数组中第一个大于start的元素的下标
        right_idx = np.searchsorted(positions, end, side='right') # 同理，查找positions数组中第一个大于end的元素的下标

        # 在窗口内的位点
        window_changes = changes[left_idx:right_idx] # 获取到pos在[start,end)宽度为window_size的那些位点的change构成的数组
                                                    # 注意这里的left_idx和right_idx都是下标
                            #left_idx:right_idx这个范围内都是pos位点号元素值>start&&<=end的，
                         # 但是总的列表长度数量基本上不是window_size，会少很多，因为很多位点应该没有数据
        # 统计不同类别的数量
        num_change_1 = np.sum(window_changes == 1)   # 统计hyper的数量
        num_change_0 = np.sum(window_changes == 0)  # 统计hypo的数量，或者用 len(window_changes) - num_change_1

        results.append({
            'window_start': start + 1,  # start+1才是每个窗口真正的位点号的起始点
            'window_end': end,   #start+1到end刚好window_size个位点被统计
            'count_change_1': num_change_1,
            'count_change_0': num_change_0,
            'total_count': num_change_1 + num_change_0 # 当前行对应的区间内的所有位点的甲基化变化方向（突变型相对于野生型）
                                                # 也是当前行对应区间的所有位点的总显著位点数量，因为一次检验会有一个变化方向
                                            #  而因为是从DMP文件中读来的，所有都是显著的
        })

    # 转换为DataFrame
    sliding_results = pd.DataFrame(results)

    # 3. 标准化处理
    max_count = sliding_results['total_count'].max()  # 所有区间中最大的那个总显著位点数
    if max_count == 0:
        max_count = 1  # 避免除零

    standardized_results = sliding_results.copy()
    standardized_results['standardized_count'] = sliding_results['total_count'] / max_count #求得当前区间总显著位点数与所有区间中
                                                                                        # 最大的那个显著位点数的比值
    # 选择需要的列用于标准化输出
    standardized_results = standardized_results[[
        'window_start', 'window_end', 'total_count', 'standardized_count'
    ]]

    print(f"滑动窗口分析完成，共生成 {len(sliding_results)} 个窗口")
    print(f"最大计数: {max_count}")

    # 4. 保存文件
    if save_files:
        # 滑动窗口结果
        sliding_file = f"slidingW_{output_identifier}.txt"
        sliding_file = os.path.join(outputdir1, sliding_file)
        sliding_results[['window_start', 'window_end', 'count_change_1', 'count_change_0']].to_csv(
            sliding_file,
            sep='\t',
            index=False,
            header=False
        )

        # 标准化结果
        std_file = f"noTitle_allDMCs_new_Standardized_slidingW_{output_identifier}.txt"
        std_file = os.path.join(outputdir1, std_file)
        # 格式化输出以匹配原始C++风格的对齐
        standardized_results.to_csv(
            std_file,
            sep='\t',
            index=False,
            header=False,
            float_format='%.6f'  # 控制浮点数精度
        )

        print(f"结果已保存:")
        print(f"  滑动窗口: {sliding_file}")
        print(f"  标准化结果: {std_file}")

    return sliding_results, standardized_results
    # 前者格式为window_start,window_end,count_change_1,count_change_0,total_count
                        #后者没有count_change，而有standardized_count，即每个区间显著位点数与最大的那个显著位点数的比值

def process_common_sites_sliding_window(common_sites_df=None,
                                        window_size=1000000,
                                        step_ratio=0.05,
                                        methytype='CpG',
                                        work_dir="."):
    """
    对共同显著位点进行滑动窗口分析
    参数：
    common_sites_df : 共同显著位点数据，如果为None则自动加载
    window_size : 滑动窗口大小
    step_ratio : 步长比例
    methytype : 甲基化类型
    """

    print(f"\n开始对共同显著位点进行滑动窗口分析...")

    and_output_dir = os.path.join(work_dir, "and_output")
    os.makedirs(and_output_dir, exist_ok=True)

    # 1. 加载共同显著位点数据
    if common_sites_df is None:
        common_file = os.path.join(and_output_dir, f"{methytype}-final_significant_sites_DMPs.txt")
        if not os.path.exists(common_file):
            print(f"错误：共同显著位点文件不存在 {common_file}")
            return None
        common_sites_df = pd.read_csv(common_file, sep='\t')
        print(f"从文件加载共同显著位点: {len(common_sites_df)} 个位点")

    if common_sites_df.empty:
        print("没有共同显著位点数据")
        return None

    # 2. 按染色体分组
    results = {}

    # 按染色体分组处理
    chr_groups = common_sites_df.groupby('Chromosome') #获取一个迭代器，可以迭代获取每个染色体号及其对应的子df

    for chr_num, chr_data in chr_groups: # 迭代获取每个染色体号及其对应的子df
        print(f"\n  处理染色体 {chr_num}: {len(chr_data)} 个位点")

        # 准备滑动窗口分析的数据，保持和all_simple_chr文件相同的格式以确保顺利进行
        window_data = pd.DataFrame({
            'position': chr_data['Position'].astype(int),
            'pvalue': chr_data['Sig_Mean_Qvalue'],
            'change': chr_data['Methylation_Change']
        })

        # 排序
        window_data = window_data.sort_values('position').reset_index(drop=True)

        # 执行滑动窗口分析，调用之前的函数就行
        try:
            sliding_results, std_results = sliding_window_analysis(
                window_data,
                window_size=window_size,
                step_ratio=step_ratio,
                save_files=True,
                output_identifier=f"common_sites_{methytype}_Chr{chr_num}",
                outputdir1=and_output_dir
            )

            # 收集结果文件
            results[chr_num] = {
                'sliding_results': sliding_results,
                'standardized_results': std_results,
                'input_data': window_data
            }

            print(f"    完成染色体 {chr_num}: {len(sliding_results)} 个窗口")

        except Exception as e:
            print(f"    错误：处理染色体 {chr_num} 失败: {e}")
            continue

    print(f"结果保存在以 'common_sites_{methytype}_Chr' 开头的文件中")

    return results

def find_max_total_in_outputs(output_dirs, methylation_type):
    """
    查找所有输出目录中指定甲基化类型的最大total值
    参数：
        output_dirs: 输出目录列表
        methylation_type: 甲基化类型（CpG, CHH, CHG）
    返回：
        max_total: 最大total值
    """
    max_total = 0

    for out_dir in output_dirs:
        mtype_dir = os.path.join(out_dir, methylation_type)
        if not os.path.exists(mtype_dir):
            continue

        # 查找所有DMP标准化文件
        std_files = glob.glob(os.path.join(mtype_dir, "noTitle_allDMCs_new_Standardized_slidingW_DMP_*.txt"))

        for std_file in std_files:
            try:
                df = pd.read_csv(
                    std_file,
                    sep=r'\s+',
                    header=None,
                    names=['start', 'end', 'total', 'normalized']
                )
                if not df.empty:
                    file_max = df['total'].max()
                    if file_max > max_total:
                        max_total = file_max
            except Exception as e:
                print(f"    警告: 读取 {std_file} 时出错: {e}")
                continue

    print(f"  {methylation_type} 类型所有染色体中最大的total值为: {max_total}")
    return max_total


def plot_methylation_sliding_windows(output_dir=None, chr_series=None,work_dir="."):
    """
    对所有滑动窗口结果进行可视化并保存到磁盘
    使用全局max_total进行标准化，使不同染色体具有可比性
    参数：
    output_dir : 指定输出目录，如果为None则自动扫描所有output_x_y目录
    chr_series : 染色体映射Series
    """

    matplotlib.use('Agg')  # 使用非交互式后端

    # 设置中文字体
    #plt.rcParams['font.sans-serif'] = ['SimHei']
    plt.rcParams['axes.unicode_minus'] = False

    print("\n开始生成甲基化滑动窗口可视化图表")

    # 扫描所有输出目录
    if output_dir is None:
        output_dirs = glob.glob(os.path.join(work_dir, "output_*_*"))
        output_dirs = [d for d in output_dirs if os.path.isdir(d)]
    else:
        output_dirs = [output_dir]

    if not output_dirs:
        print("未找到输出目录")
        return

    total_plots = 0
    methylation_types = ['CpG', 'CHH', 'CHG']

    # 为每种甲基化类型分别计算全局max_total
    max_totals = {}
    for mtype in methylation_types:
        max_totals[mtype] = find_max_total_in_outputs(output_dirs, mtype)
        if max_totals[mtype] == 0:
            print(f"  警告: {mtype} 类型未找到有效的total值，将使用1作为默认值")
            max_totals[mtype] = 1

    for out_dir in output_dirs:
        print(f"\n处理目录: {out_dir}")

        for mtype in methylation_types:
            mtype_dir = os.path.join(out_dir, mtype)
            if not os.path.exists(mtype_dir):
                continue

            print(f"  处理甲基化类型: {mtype}（使用全局max_total={max_totals[mtype]}）")

            # 查找当前甲基化目录下所有DMP滑动窗口文件
            dmp_sliding_files = glob.glob(os.path.join(mtype_dir, "slidingW_DMP_*.txt"))

            # 按前缀分组处理，同一前缀的所有染色体绘制在一张大图上
            prefix_groups = {}
            for dmp_sliding_file in dmp_sliding_files:
                basename = os.path.basename(dmp_sliding_file)
                match = re.search(r'slidingW_DMP_(.+)_Chr(\d+)\.txt$', basename)
                if not match:
                    continue

                prefix = match.group(1)  # replicate1_replicate2
                chr_num = match.group(2)

                if prefix not in prefix_groups:
                    prefix_groups[prefix] = []
                prefix_groups[prefix].append(chr_num)

            # 对每个前缀组进行处理
            for prefix, chr_nums in prefix_groups.items():
                # 对染色体编号排序
                chr_nums = sorted(chr_nums, key=lambda x: int(x))

                # 收集所有染色体的数据
                all_chrom_data = []
                chrom_names = []

                for chr_num in chr_nums:
                    # 构建对应的标准化文件路径
                    chr_name = get_chr_name(chr_num, chr_series)

                    # 修正：使用正确的变量名构建文件路径
                    dmp_sliding_file = os.path.join(mtype_dir, f"slidingW_DMP_{prefix}_Chr{chr_num}.txt")
                    dmp_std_file = os.path.join(mtype_dir,
                                                f"noTitle_allDMCs_new_Standardized_slidingW_DMP_{prefix}_Chr{chr_num}.txt")
                    ndmp_std_file = os.path.join(mtype_dir,
                                                 f"noTitle_allDMCs_new_Standardized_slidingW_N-DMP_{prefix}_Chr{chr_num}.txt")

                    if not all(os.path.exists(f) for f in [dmp_sliding_file, dmp_std_file, ndmp_std_file]):
                        print(f"    警告: Chr{chr_num} 的文件不完整，跳过")
                        continue

                    try:
                        # 读取数据
                        sliding_df = pd.read_csv(
                            dmp_sliding_file,
                            sep=r'\s+',
                            header=None,
                            names=['start', 'end', 'hyper', 'hypo']
                        )

                        dmp_std_df = pd.read_csv(
                            dmp_std_file,
                            sep=r'\s+',
                            header=None,
                            names=['start', 'end', 'total', 'normalized']
                        )

                        ndmp_std_df = pd.read_csv(
                            ndmp_std_file,
                            sep=r'\s+',
                            header=None,
                            names=['start', 'end', 'total', 'ndmp_normalized']
                        )

                        if sliding_df.empty or dmp_std_df.empty or ndmp_std_df.empty:
                            print(f"    警告: Chr{chr_num} 的数据为空，跳过")
                            continue

                        # 使用全局max_total重新计算比率
                        max_total = max_totals[mtype]

                        x = (sliding_df['start'] + sliding_df['end']) / 2

                        # 重新计算所有比率（使用全局max_total）
                        y_dmp = (sliding_df['hyper'] + sliding_df['hypo']) / max_total
                        y_hyper = sliding_df['hyper'] / max_total
                        y_hypo = sliding_df['hypo'] / max_total
                        y_ndmp = ndmp_std_df['ndmp_normalized']  # NDMP保持原样

                        # 处理可能的长度不一致
                        max_len = max(len(x), len(y_dmp), len(y_hyper), len(y_hypo), len(y_ndmp))
                        x = x.reindex(range(max_len), fill_value=0)
                        y_dmp = y_dmp.reindex(range(max_len), fill_value=0)
                        y_hyper = y_hyper.reindex(range(max_len), fill_value=0)
                        y_hypo = y_hypo.reindex(range(max_len), fill_value=0)
                        y_ndmp = y_ndmp.reindex(range(max_len), fill_value=0)

                        # 存储数据
                        all_chrom_data.append({
                            'x': x,
                            'y_dmp': y_dmp,
                            'y_hyper': y_hyper,
                            'y_hypo': y_hypo,
                            'y_ndmp': y_ndmp
                        })
                        chrom_names.append(chr_name)

                        print(f"    成功加载染色体 {chr_name} 的数据")

                    except Exception as e:
                        print(f"    处理 Chr{chr_num} 时出错: {e}")
                        continue

                # 如果有数据，绘制大图
                if all_chrom_data:
                    try:
                        # 创建大图，每个染色体一个子图
                        n_chromosomes = len(all_chrom_data)
                        fig, axes = plt.subplots(n_chromosomes, 1, figsize=(15, 3 * n_chromosomes))

                        # 如果只有一个染色体，axes不是数组，需要转换为数组
                        if n_chromosomes == 1:
                            axes = [axes]

                        # 设置总标题
                        fig.suptitle(f'{mtype} Methylation Analysis - {prefix} (Global Normalized)',
                                     fontsize=16, fontfamily='DejaVu Sans')

                        # 绘制每个染色体的子图
                        for idx, (chrom_data, chrom_name) in enumerate(zip(all_chrom_data, chrom_names)):
                            ax = axes[idx]

                            # 绘制所有数据线
                            ax.plot(chrom_data['x'], chrom_data['y_dmp'], label='DMP', color='red', linewidth=1.5)
                            ax.plot(chrom_data['x'], chrom_data['y_hyper'], label='Hyper-ratio', color='green',
                                    linewidth=1)
                            ax.plot(chrom_data['x'], chrom_data['y_hypo'], label='Hypo-ratio', color='blue',
                                    linewidth=1)
                            ax.plot(chrom_data['x'], chrom_data['y_ndmp'], label='NDMP', color='darkgray',
                                    linewidth=1.5)

                            ax.set_ylim(bottom=0)
                            ax.set_ylabel('Ratio', fontsize=10, fontfamily='DejaVu Sans')

                            # 设置子图标题
                            ax.set_title(f'{chrom_name}', fontsize=18, fontfamily='DejaVu Sans', pad=20,y=-0.4)

                            # 添加网格
                            ax.grid(True, alpha=0.3)

                            # 只在第一个子图添加完整图例
                            if idx == 0:
                                ax.legend(loc='upper right', ncol=2, fontsize=8, framealpha=0.7)

                        # 调整布局
                        plt.tight_layout(rect=[0, 0, 1, 0.96])  # 为总标题留出空间

                        # 保存图片
                        plot_filename = os.path.join(mtype_dir,
                                                     f"methylation_plot_{mtype}_{prefix}_all_chromosomes.png")
                        plt.savefig(plot_filename, dpi=300, bbox_inches='tight')
                        plt.close()

                        total_plots += 1
                        print(f"    成功生成大图: {mtype}_{prefix} -> {os.path.basename(plot_filename)}")

                    except Exception as e:
                        print(f"    绘制大图时出错: {e}")
                        continue

    print(f"\n图表生成完成！共生成 {total_plots} 张大图")


def plot_common_sites_sliding_windows(methytype='CpG', chr_series=None, work_dir="."):
    """
    对共同显著位点的滑动窗口结果进行可视化并保存到磁盘
    使用全局max_total进行标准化，将所有染色体绘制在一张大图上
    """
    matplotlib.use('Agg')

    # 设置全局字体为DejaVu Sans
    plt.rcParams['font.family'] = 'DejaVu Sans'
    plt.rcParams['axes.unicode_minus'] = False

    print(f"\n开始生成共同显著位点({methytype})的滑动窗口可视化图表...")

    and_output_dir = os.path.join(work_dir, "and_output")

    # 先找到该甲基化类型的全局max_total
    print(f"  正在查找 {methytype} 的全局max_total...")
    max_total = 0
    std_files = glob.glob(
        os.path.join(and_output_dir, f"noTitle_allDMCs_new_Standardized_slidingW_common_sites_{methytype}_Chr*.txt"))

    for std_file in std_files:
        try:
            df = pd.read_csv(std_file, sep=r'\s+', header=None, names=['start', 'end', 'total', 'normalized'])
            if not df.empty:
                file_max = df['total'].max()
                if file_max > max_total:
                    max_total = file_max
        except Exception as e:
            print(f"    警告: 读取 {std_file} 时出错: {e}")

    if max_total == 0:
        print(f"  警告: 未找到有效的total值，使用默认值1")
        max_total = 1
    else:
        print(f"  {methytype} 的全局max_total为: {max_total}")

    # 读取DMR数据 - 根据甲基化类型选择对应的DMR文件
    print("  正在读取DMR数据...")
    dmr_file = os.path.join(and_output_dir, f"{methytype}-final_significant_regions_DMRs.txt")
    dmr_data = {}
    if os.path.exists(dmr_file):
        try:
            dmr_df = pd.read_csv(dmr_file, sep=r'\s+')

            for _, row in dmr_df.iterrows():
                try:
                    chrom = str(row['Chromosome'])  # 基于列名访问
                    direction = int(row['Direction'])  # 基于列名访问
                    start = int(row['DMR_start'])  # 基于列名访问
                    end = int(row['DMR_end'])  # 基于列名访问

                    # 计算中点
                    mid = (start + end) / 2

                    # 提取染色体数字部分
                    chrom_num = str(chrom).replace('Chr', '').replace('chr', '')

                    if chrom_num not in dmr_data:
                        dmr_data[chrom_num] = []

                    dmr_data[chrom_num].append((mid, direction))
                except (ValueError, IndexError):
                    continue
            print(f"  成功加载 {sum(len(dmrs) for dmrs in dmr_data.values())} 个DMR")
        except Exception as e:
            print(f"  读取DMR文件时出错: {e}")
    else:
        print(f"  警告: DMR文件 {dmr_file} 不存在")

    # 动态获取染色体列表
    sliding_files = glob.glob(os.path.join(and_output_dir, f"slidingW_common_sites_{methytype}_Chr*.txt"))

    # 从文件名中提取染色体编号并排序
    chromosomes = []
    for file in sliding_files:
        match = re.search(r'slidingW_common_sites_.+_Chr(\d+)\.txt$', os.path.basename(file))
        if match:
            chr_num = match.group(1)
            if chr_num not in chromosomes:
                chromosomes.append(chr_num)

    # 按数字顺序排序染色体
    chromosomes.sort(key=int)

    if not chromosomes:
        print(f"未找到共同显著位点({methytype})的滑动窗口文件")
        return

    print(f"找到 {len(chromosomes)} 个染色体的滑动窗口文件")

    all_chrom_data = []
    chrom_nums1 = []

    for chr_num in chromosomes:
        sliding_file = os.path.join(and_output_dir, f"slidingW_common_sites_{methytype}_Chr{chr_num}.txt")
        std_file = os.path.join(and_output_dir,
                                f"noTitle_allDMCs_new_Standardized_slidingW_common_sites_{methytype}_Chr{chr_num}.txt")

        if not all(os.path.exists(f) for f in [sliding_file, std_file]):
            print(f"    警告: Chr{chr_num} 的文件不完整，跳过")
            continue

        try:
            # 读取数据
            sliding_df = pd.read_csv(sliding_file, sep=r'\s+', header=None, names=['start', 'end', 'hyper', 'hypo'])
            std_df = pd.read_csv(std_file, sep=r'\s+', header=None, names=['start', 'end', 'total', 'normalized'])

            if sliding_df.empty or std_df.empty:
                print(f"    警告: Chr{chr_num} 的数据为空，跳过")
                continue

            if len(sliding_df) != len(std_df):
                print(f"    警告: Chr{chr_num} 的数据长度不一致，跳过")
                continue

            # 新增：读取 output_1_1 的 NDMP 数据
            ndmp_file = os.path.join(work_dir, "output_wt1_mut1", methytype,
                                     f"noTitle_allDMCs_new_Standardized_slidingW_N-DMP_wt_replicate1_mut_replicate1_Chr{chr_num}.txt")
            y_ndmp = None  # 初始化为None
            if os.path.exists(ndmp_file):
                try:
                    ndmp_df = pd.read_csv(ndmp_file, sep=r'\s+', header=None,
                                          names=['start', 'end', 'total', 'ndmp_normalized'])
                    if not ndmp_df.empty:
                        y_ndmp = ndmp_df['ndmp_normalized']
                        print(f"    成功读取 output_1_1 的 NDMP 数据: Chr{chr_num}")
                except Exception as e:
                    print(f"    警告: 读取 output_1_1 NDMP 数据失败 (Chr{chr_num}): {e}")
            else:
                print(f"    提示: output_1_1 的 NDMP 文件不存在 (Chr{chr_num})")

            # 使用全局max_total重新计算比率
            x = (sliding_df['start'] + sliding_df['end']) / 2
            y_total = (sliding_df['hyper'] + sliding_df['hypo']) / max_total
            y_hyper = sliding_df['hyper'] / max_total
            y_hypo = sliding_df['hypo'] / max_total

            max_len = max(len(x), len(y_total), len(y_hyper), len(y_hypo))
            if y_ndmp is not None:
                max_len = max(max_len, len(y_ndmp))
            x = x.reindex(range(max_len), fill_value=0)
            y_total = y_total.reindex(range(max_len), fill_value=0)
            y_hyper = y_hyper.reindex(range(max_len), fill_value=0)
            y_hypo = y_hypo.reindex(range(max_len), fill_value=0)
            if y_ndmp is not None:
                y_ndmp = y_ndmp.reindex(range(max_len), fill_value=0)

            # 存储数据
            chr_real_name = get_chr_name(chr_num, chr_series)
            all_chrom_data.append({
                'x': x,
                'y_total': y_total,
                'y_hyper': y_hyper,
                'y_hypo': y_hypo,
                'y_ndmp': y_ndmp
            })
            chrom_nums1.append(chr_num)
            print(f"    成功加载染色体 {chr_real_name} 的数据")

        except Exception as e:
            print(f"    处理 Chr{chr_num} 时出错: {e}")
            continue

    # 如果有数据，绘制大图
    if all_chrom_data:
        try:
            # 创建大图，包含多个子图
            n_chromosomes = len(all_chrom_data)
            fig, axes = plt.subplots(n_chromosomes, 1, figsize=(15, 3 * n_chromosomes))

            # 如果只有一个染色体，axes不是数组，需要转换为数组
            if n_chromosomes == 1:
                axes = [axes]

            # 设置总标题
            fig.suptitle(f'Distribution of Common Significant Sites - {methytype} context',
                         fontsize=16, fontfamily='DejaVu Sans')

            # 绘制每个染色体的子图
            for idx, (chrom_data, chrom_num1) in enumerate(zip(all_chrom_data, chrom_nums1)):
                ax = axes[idx]

                # 绘制DMP数据
                ax.plot(chrom_data['x'], chrom_data['y_total'], label='DMP', color='red', linewidth=2)
                ax.plot(chrom_data['x'], chrom_data['y_hyper'], label='hyper-methylation', color='green', linewidth=1.5)
                ax.plot(chrom_data['x'], chrom_data['y_hypo'], label='hypo-methylation', color='blue', linewidth=1.5)


                if chrom_data['y_ndmp'] is not None:
                    ax.plot(chrom_data['x'], chrom_data['y_ndmp'], label='NDMP',
                            color='darkgray', linewidth=1.5)
                # 统一设置Y轴范围为0到1.2
                ax.set_ylim(0, 1.2)

                # 获取染色体号（从chrom_name中提取）
                chrom_num1 = chrom_num1.replace('chr', '').replace('Chr', '')

                if idx == 0:  # 只在第一个子图时打印一次
                    print(f"  调试: dmr_data keys = {list(dmr_data.keys())}")
                print(f"    {get_chr_name(chrom_num1,chr_series)} -> chrom_num = '{chrom_num1}', 在dmr_data中: {chrom_num1 in dmr_data}")

                # 添加DMR标记
                if chrom_num1 in dmr_data:
                    for mid, direction in dmr_data[chrom_num1]:
                        # 根据direction选择颜色：1=hyper(绿色), 0=hypo(蓝色)
                        color = 'green' if direction == 1 else 'blue'
                        # 在DMR中点位置添加竖线，显示在Y轴1.0到1.2的范围内
                        ax.axvline(x=mid, ymin=0.9, ymax=1, color=color, linewidth=2, alpha=0.7)

                # 设置子图标题和标签
                ax.text(0.5, -0.2, f"{get_chr_name(chrom_num1,chr_series)}",
                        transform=ax.transAxes,
                        fontfamily='DejaVu Sans',
                        ha='center', va='top',
                        fontsize=15)

                # 添加网格
                ax.grid(True, alpha=0.3)

                # 只在第一个子图添加完整图例
                if idx == 0:
                    # 创建自定义图例条目，包括DMR标记
                    from matplotlib.lines import Line2D
                    legend_elements = [
                        Line2D([0], [0], color='red', linewidth=2, label='DMP'),
                        Line2D([0], [0], color='green', linewidth=1.5, label='hyper-DMP'),
                        Line2D([0], [0], color='blue', linewidth=1.5, label='hypo-DMP'),
                        Line2D([0], [0], color='green', linewidth=2, label='hyper-DMR'),
                        Line2D([0], [0], color='blue', linewidth=2, label='hypo-DMR')
                    ]
                    # 如果有NDMP数据，添加到图例
                    if chrom_data['y_ndmp'] is not None:
                        legend_elements.append(
                            Line2D([0], [0], color='darkgray', linewidth=1.5,
                                  label='NDMP')
                        )
                    ax.legend(handles=legend_elements, loc='upper right', bbox_to_anchor=(1.0, 0.95),
                              ncol=5, fontsize=8, framealpha=0.7)

            # 调整布局
            plt.tight_layout(rect=[0, 0, 1, 0.96])  # 为总标题留出空间

            # 保存图片
            plot_filename = os.path.join(and_output_dir, f"common_sites_plot_{methytype}_all_chromosomes.png")
            plt.savefig(plot_filename, dpi=300, bbox_inches='tight')
            plt.close()

            print(f"成功生成大图: {methytype} -> {os.path.basename(plot_filename)}")

        except Exception as e:
            print(f"绘制大图时出错: {e}")
    else:
        print(f"未找到 {methytype} 的有效数据")



def rename_chromosome_files(chr_series, work_dir="."):
    """
    批量重命名所有包含染色体号的文件，将Chr数字替换为真实染色体名称

    参数：
        chr_series: 染色体映射Series (染色体名称 -> 索引)
        work_dir: 工作目录
    """
    print("\n开始批量重命名染色体文件...")

    if chr_series is None or len(chr_series) == 0:
        print("错误：染色体映射为空，跳过重命名")
        return

    # 创建反向映射：数字索引 -> 染色体名称
    # chr_series 的 index 是染色体名称，value 是数字索引
    index_to_chr = {i: chr_series.index[i] for i in range(len(chr_series))}

    print(f"染色体映射: {index_to_chr}")

    # 定义需要搜索的目录模式
    search_dirs = []

    # 添加所有 output_x_y 目录
    output_dirs = glob.glob(os.path.join(work_dir, "output_*_*"))
    search_dirs.extend([d for d in output_dirs if os.path.isdir(d)])

    # 添加 and_output 目录
    and_output_dir = os.path.join(work_dir, "and_output")
    if os.path.exists(and_output_dir):
        search_dirs.append(and_output_dir)

    if not search_dirs:
        print("未找到需要处理的目录")
        return

    print(f"将在 {len(search_dirs)} 个目录中搜索文件")

    # 统计信息
    total_renamed = 0
    failed_renames = 0

    # 定义匹配染色体号的正则表达式模式
    # 匹配 Chr 后跟数字的模式，如 Chr1, Chr12 等
    chr_pattern = re.compile(r'(.*?)Chr(\d+)(.*?)$')

    # 遍历所有目录
    for search_dir in search_dirs:
        print(f"\n处理目录: {search_dir}")

        # 递归遍历目录中的所有文件
        for root, dirs, files in os.walk(search_dir):
            for filename in files:
                # 检查文件名是否包含 Chr数字 模式
                match = chr_pattern.match(filename)

                if match:
                    prefix = match.group(1)  # Chr 之前的部分
                    chr_num = int(match.group(2))  # 染色体数字
                    suffix = match.group(3)  # Chr数字 之后的部分

                    # 根据 chr_series 获取真实染色体名称
                    # 注意：文件名中的 Chr1 对应 index 0
                    chr_index = chr_num - 1

                    if chr_index not in index_to_chr:
                        print(f"  警告: Chr{chr_num} 不在映射表中，跳过文件 {filename}")
                        continue

                    real_chr_name = index_to_chr[chr_index]

                    # 构造新文件名
                    # 如果染色体名称本身包含 'chr' 前缀，直接使用
                    # 否则使用 Chr 前缀
                    if real_chr_name.lower().startswith('chr'):
                        chr_part = real_chr_name
                    else:
                        chr_part = f"Chr{real_chr_name}"

                    new_filename = f"{prefix}{chr_part}{suffix}"

                    # 如果新旧文件名相同，跳过
                    if filename == new_filename:
                        continue

                    # 构造完整路径
                    old_path = os.path.join(root, filename)
                    new_path = os.path.join(root, new_filename)

                    # 检查目标文件是否已存在
                    if os.path.exists(new_path):
                        print(f"  警告: 目标文件已存在，跳过重命名: {filename} -> {new_filename}")
                        failed_renames += 1
                        continue

                    # 执行重命名
                    try:
                        os.rename(old_path, new_path)
                        total_renamed += 1
                        print(f"  {filename} -> {new_filename}")
                    except Exception as e:
                        print(f" 重命名失败: {filename} -> {new_filename}, 错误: {e}")
                        failed_renames += 1

    # 输出统计信息
    print(f"\n重命名完成！")
    print(f"  成功重命名: {total_renamed} 个文件")
    if failed_renames > 0:
        print(f"  失败或跳过: {failed_renames} 个文件")


def convert_output_to_csv(work_dir="."):
    """
    将 and_output 目录下的 final DMP 和 final DMR 文件转换为 CSV 格式（逗号分隔）

    参数：
        work_dir: 工作目录
    """

    and_output_dir = os.path.join(work_dir, "and_output")

    if not os.path.exists(and_output_dir):
        print(f"错误：目录 {and_output_dir} 不存在")
        return 0

    # 定义需要转换的文件模式
    file_patterns = [
        "*-final_significant_sites_DMPs.txt",  # final DMP 文件
        "*-final_significant_regions_DMRs.txt"  # final DMR 文件
    ]

    converted_count = 0
    failed_count = 0

    for pattern in file_patterns:
        # 搜索匹配的文件
        matching_files = glob.glob(os.path.join(and_output_dir, pattern))

        for txt_file in matching_files:
            try:
                # 读取制表符分隔的文件
                df = pd.read_csv(txt_file, sep=r'\s+')

                if df.empty:
                    print(f"  跳过空文件: {os.path.basename(txt_file)}")
                    continue

                # 生成 CSV 文件路径（将 .txt 替换为 .csv）
                csv_file = txt_file.replace('.txt', '.csv')

                # 保存为逗号分隔的 CSV 文件
                df.to_csv(csv_file, sep=',', index=False)

                print(f"  成功转换： {os.path.basename(txt_file)} -> {os.path.basename(csv_file)}")
                converted_count += 1

            except Exception as e:
                print(f"  转换失败: {os.path.basename(txt_file)}, 错误: {e}")
                failed_count += 1

    print(f"\\n转换完成！")
    print(f"  成功转换: {converted_count} 个文件")
    if failed_count > 0:
        print(f"  转换失败: {failed_count} 个文件")

    return converted_count


def convert_chromosome_to_names(chr_series, work_dir="."):
    """
    将 and_output 目录下 final DMP 和 final DMR 文件中的 Chromosome 列
    从数值转换为真实的染色体名称

    参数：
        chr_series: 染色体映射 Series (染色体名称 -> 索引，索引从0开始)
        work_dir: 工作目录

    返回：
        成功转换的文件数量
    """

    if chr_series is None or len(chr_series) == 0:
        print("错误：染色体映射为空，跳过转换")
        return 0

    and_output_dir = os.path.join(work_dir, "and_output")

    if not os.path.exists(and_output_dir):
        print(f"错误：目录 {and_output_dir} 不存在")
        return 0

    # 创建数值到染色体名称的映射
    # 文件中的 Chromosome 列是从1开始的数字，需要转换为 chr_series 中的染色体名称
    # chr_series.index[0] 对应文件中的 1
    index_to_chr = {i + 1: chr_series.index[i] for i in range(len(chr_series))}

    print(f"染色体映射表: {index_to_chr}")

    # 定义需要处理的文件模式
    file_patterns = [
        "*-final_significant_sites_DMPs.txt",  # final DMP 文件
        "*-final_significant_regions_DMRs.txt"  # final DMR 文件
    ]

    converted_count = 0
    failed_count = 0

    for pattern in file_patterns:
        # 搜索匹配的文件
        matching_files = glob.glob(os.path.join(and_output_dir, pattern))

        for file_path in matching_files:
            try:
                # 读取文件
                df = pd.read_csv(file_path, sep=r'\s+')

                if df.empty:
                    print(f"  跳过空文件: {os.path.basename(file_path)}")
                    continue

                # 检查是否有 Chromosome 列
                if 'Chromosome' not in df.columns:
                    print(f"  警告: {os.path.basename(file_path)} 中没有 Chromosome 列，跳过")
                    continue

                # 保存原始的 Chromosome 列用于调试
                original_chrs = df['Chromosome'].unique()

                # 转换 Chromosome 列
                # 先确保是整数类型
                df['Chromosome'] = df['Chromosome'].astype(int)

                # 使用映射转换为染色体名称
                df['Chromosome'] = df['Chromosome'].map(index_to_chr)

                # 检查是否有未成功映射的值
                if df['Chromosome'].isna().any():
                    unmapped_count = df['Chromosome'].isna().sum()
                    print(f"  警告: {os.path.basename(file_path)} 中有 {unmapped_count} 个染色体编号无法映射")
                    # 可选：删除无法映射的行
                    df = df.dropna(subset=['Chromosome'])

                # 保存回原文件（覆盖）
                df.to_csv(file_path, sep='\t', index=False)

                print(f"    成功转换: {os.path.basename(file_path)}")
                print(f"    原始编号: {sorted(original_chrs)}")
                print(f"    转换后: {sorted(df['Chromosome'].unique())}")
                converted_count += 1

            except Exception as e:
                print(f"  ✗ 转换失败: {os.path.basename(file_path)}")
                print(f"    错误信息: {e}")
                import traceback
                traceback.print_exc()
                failed_count += 1

    # 输出统计信息
    print(f"\n染色体编号转换完成！")
    print(f"  成功转换: {converted_count} 个文件")
    if failed_count > 0:
        print(f"  转换失败: {failed_count} 个文件")

    return converted_count


def main():
    global DMR_QVALUE_THRESHOLD, VOTE_THRESHOLD, DMR_ENGINE
    global AUTO_DMP_VOTE_THRESHOLD, AUTO_DMR_VOTE_THRESHOLD, AUTO_VOTE_THRESHOLD_REPORT_ONLY
    global AUTO_DMP_VOTE_THRESHOLDS, AUTO_DMR_VOTE_THRESHOLDS
    global AUTO_QVALUE_TWOSTEP, AUTO_QVALUE_REPORT_ONLY
    global AUTO_QVALUE_P_CUTOFF, AUTO_QVALUE_MIN_CANDIDATES
    global AUTO_QVALUE_USE_SMOOTH, AUTO_QVALUE_SMOOTH_SIGMA
    global DMP_LOWDIFF_STRICT_VOTE, DMP_LOWDIFF_CUTOFF, DMP_LOWDIFF_STRICT_VOTE_REPORT_ONLY

    start_time = time.time()
    parser = argparse.ArgumentParser(
        description="MultiDMPcaller: pairwise Fisher tests, FDR correction, DMP/DMR calling"
    )
    # ===== 兼容旧版位置参数：python script.py n m dir2 dir1 biotype =====
    parser.add_argument("n_pos", type=int, nargs="?", help=argparse.SUPPRESS)
    parser.add_argument("m_pos", type=int, nargs="?", help=argparse.SUPPRESS)
    parser.add_argument("dir2_pos", nargs="?", help=argparse.SUPPRESS)
    parser.add_argument("dir1_pos", nargs="?", help=argparse.SUPPRESS)
    parser.add_argument("biotype_pos", type=int, nargs="?", choices=[0, 1, 2], help=argparse.SUPPRESS)

    # ===== 新版命名参数：python script.py --n 2 --m 2 --dir2 wt --dir1 mut --biotype 0 =====
    parser.add_argument("--n", dest="n_opt", metavar="N", type=int, help="对照组/野生型重复数")
    parser.add_argument("--m", dest="m_opt", metavar="M", type=int, help="实验组/突变型重复数")
    parser.add_argument("--dir-wt", dest="dir2_opt", metavar="DIR_WT", help="对照组/野生型样本目录")
    parser.add_argument("--dir-mut", dest="dir1_opt", metavar="DIR_MUT", help="实验组/突变型样本目录")
    parser.add_argument("--biotype", dest="biotype_opt", metavar="BIOTYPE", type=int, choices=[0, 1, 2], help="0=动物, 1=植物, 2=不过滤")
    parser.add_argument(
        "--meth-diff",
        type=float,
        default=0.0,
        help="DMP最终筛选所需的最小绝对甲基化差异，范围0-1；例如0.25表示25%%。默认0.0，兼容旧行为。"
    )
    parser.add_argument(
        "--auto-meth-diff",
        action="store_true",
        help="启用 mianjifa auto-methdiff：在 pairwise FDR完成后，根据q显著且尚未经过MethDiff筛选的位点之raw MethDiff分布估计全局阈值；该阈值同时用于auto-vote支持构建和final common DMP。默认关闭。"
    )
    parser.add_argument(
        "--auto-meth-diff-report-only",
        action="store_true",
        help="仅输出 mianjifa auto-methdiff 阈值诊断表和分布图，不改变实际 DMP 判定阈值。默认关闭。"
    )
    parser.add_argument(
        "--auto-meth-diff-cut-percent",
        type=float,
        default=0.05,
        help="mianjifa auto-methdiff 从0向两侧切除的直方图面积比例。默认0.05。"
    )
    parser.add_argument(
        "--auto-meth-diff-fallback",
        type=float,
        default=0.3,
        help="mianjifa auto-methdiff 估计失败时回退使用的 methdiff 阈值。默认0.3。"
    )
    parser.add_argument(
        "--auto-meth-diff-aggregate",
        choices=["median", "mean", "max", "min"],
        default="median",
        help="将各比较左右阈值聚合为一个全局 abs(methdiff) 阈值的方法。默认median。"
    )
    parser.add_argument(
        "--q-cpg",
        type=float,
        default=DMP_QVALUE_THRESHOLDS["CpG"],
        help="CpG DMP的q-value阈值，范围0-1。默认使用原代码阈值。"
    )
    parser.add_argument(
        "--q-chg",
        type=float,
        default=DMP_QVALUE_THRESHOLDS["CHG"],
        help="CHG DMP的q-value阈值，范围0-1。默认使用原代码阈值。"
    )
    parser.add_argument(
        "--q-chh",
        type=float,
        default=DMP_QVALUE_THRESHOLDS["CHH"],
        help="CHH DMP的q-value阈值，范围0-1。默认使用原代码阈值。"
    )
    parser.add_argument(
        "--dmr-q",
        type=float,
        default=DMR_QVALUE_THRESHOLD,
        help="DMR的q-value阈值，范围0-1。默认使用原代码阈值。"
    )
    parser.add_argument(
        "--auto-qvalue-twostep",
        action="store_true",
        help="仅对两步法context自动估计DMP q-value阈值：在p值预筛选子集内寻找max(qvalue-pvalue)点，并用该点qvalue作为该pair/context阈值。默认关闭。"
    )
    parser.add_argument(
        "--auto-qvalue-report-only",
        action="store_true",
        help="仅输出两步法自动q-value阈值诊断表，不改变DMP判定阈值。默认关闭。"
    )
    parser.add_argument(
        "--auto-qvalue-p-cutoff",
        type=float,
        default=AUTO_QVALUE_P_CUTOFF,
        help="两步法自动q-value阈值估计使用的p-value候选上限，默认0.05。"
    )
    parser.add_argument(
        "--auto-qvalue-min-candidates",
        type=int,
        default=AUTO_QVALUE_MIN_CANDIDATES,
        help="估计自动q-value阈值所需的最少候选位点数；不足时回退到固定q阈值。默认10。"
    )
    parser.add_argument(
        "--auto-qvalue-use-smooth",
        action="store_true",
        help="使用平滑后的(qvalue-pvalue)曲线寻找最大差值点。默认关闭；建议正式分析使用默认raw diff。"
    )
    parser.add_argument(
        "--auto-qvalue-smooth-sigma",
        type=float,
        default=AUTO_QVALUE_SMOOTH_SIGMA,
        help="--auto-qvalue-use-smooth 启用时的Gaussian sigma。默认4。"
    )
    parser.add_argument(
        "--vote-threshold",
        type=float,
        default=VOTE_THRESHOLD,
        help="最终DMP/DMR跨replicate组合投票阈值，范围(0,1]；默认0.6667，即原代码2/3规则。"
    )
    parser.add_argument(
        "--auto-dmp-vote-threshold",
        action="store_true",
        help="自动估计 final DMP 的整数投票 required_count。默认关闭。"
    )
    parser.add_argument(
        "--auto-dmr-vote-threshold",
        action="store_true",
        help="自动估计 final DMR 的整数投票 required_count。默认关闭。"
    )
    parser.add_argument(
        "--auto-vote-threshold-report-only",
        action="store_true",
        help="仅计算并输出自动投票阈值和分布图，不改变 final DMP/DMR 判定。默认关闭。"
    )
    parser.add_argument(
        "--dmp-lowdiff-strict-vote",
        action="store_true",
        help="启用 final DMP 低差异严格投票后处理：先按正常q-value+投票得到provisional DMP；若boundary abs(MethDiff)<=cutoff，则要求更高投票数。默认关闭，保持旧行为。"
    )
    parser.add_argument(
        "--dmp-lowdiff-cutoff",
        type=float,
        default=DMP_LOWDIFF_CUTOFF,
        help="final DMP lowdiff strict vote 的 boundary abs(MethDiff) cutoff。默认0.3。"
    )
    parser.add_argument(
        "--dmp-lowdiff-strict-vote-report-only",
        action="store_true",
        help="仅输出 final DMP lowdiff strict vote 诊断表和summary，不改变最终DMP文件。默认关闭。"
    )
    parser.add_argument(
        "--skip-dmr",
        action="store_true",
        help="跳过所有 DMR 相关步骤，仅输出 DMP/final DMP 结果。默认关闭，兼容原行为。"
    )
    parser.add_argument(
        "--skip-window",
        action="store_true",
        help="跳过所有滑动窗口和可视化绘图步骤，仅输出表格结果。默认关闭，兼容原行为。"
    )
    parser.add_argument(
        "--dmr-engine",
        choices=["python", "cpp"],
        default="python",
        help="DMR候选区域识别引擎：python=原始Python实现；cpp=使用dmr_step1 + dmr_step2_dynamic。默认python，便于验证兼容性。"
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=1,
        help="并行 worker 进程数。默认1为全流程串行；threads>1 时自动并行所有安全阶段，包括 newtoboth 文件转换、replicate-pair 处理和 common DMR 汇总。建议小数据2-4，大数据根据内存/I/O谨慎设置。"
    )
    args = parser.parse_args()

    def choose_arg(opt_value, pos_value, name):
        """同时兼容新版 --xxx 参数和旧版位置参数。"""
        if opt_value is not None and pos_value is not None and opt_value != pos_value:
            parser.error(f"参数冲突：--{name}={opt_value} 与旧版位置参数 {pos_value} 不一致")
        if opt_value is not None:
            return opt_value
        if pos_value is not None:
            return pos_value
        parser.error(f"缺少必要参数：--{name}")

    args.n = choose_arg(args.n_opt, args.n_pos, "n")
    args.m = choose_arg(args.m_opt, args.m_pos, "m")
    args.dir2 = choose_arg(args.dir2_opt, args.dir2_pos, "dir2")
    args.dir1 = choose_arg(args.dir1_opt, args.dir1_pos, "dir1")
    args.biotype = choose_arg(args.biotype_opt, args.biotype_pos, "biotype")

    m = args.m
    n = args.n
    dir1 = args.dir1
    dir2 = args.dir2
    biotype = args.biotype
    meth_diff_threshold = args.meth_diff

    numeric_thresholds = {
        "--meth-diff": meth_diff_threshold,
        "--q-cpg": args.q_cpg,
        "--q-chg": args.q_chg,
        "--q-chh": args.q_chh,
        "--dmr-q": args.dmr_q,
        "--auto-qvalue-p-cutoff": args.auto_qvalue_p_cutoff,
        "--auto-meth-diff-fallback": args.auto_meth_diff_fallback,
        "--dmp-lowdiff-cutoff": args.dmp_lowdiff_cutoff,
    }
    for name, value in numeric_thresholds.items():
        if not 0 <= value <= 1:
            print(f"错误：{name} 必须在 0 到 1 之间，例如 0.25 表示 25%")
            sys.exit(1)
    if not 0 < args.vote_threshold <= 1:
        print("错误：--vote-threshold 必须在 (0, 1] 之间，例如 0.667 表示2/3多数投票")
        sys.exit(1)
    if not 0 < args.auto_meth_diff_cut_percent < 1:
        print("错误：--auto-meth-diff-cut-percent 必须在 (0, 1) 之间，例如0.05")
        sys.exit(1)
    if args.threads < 1:
        print("错误：--threads 必须是 >= 1 的整数")
        sys.exit(1)
    if args.auto_qvalue_min_candidates < 1:
        print("错误：--auto-qvalue-min-candidates 必须是 >= 1 的整数")
        sys.exit(1)
    if args.auto_qvalue_smooth_sigma <= 0:
        print("错误：--auto-qvalue-smooth-sigma 必须大于0")
        sys.exit(1)

    # 将命令行参数写回全局阈值配置，保持原有 get_dmp_threshold() 调用链不变
    DMP_QVALUE_THRESHOLDS["CpG"] = args.q_cpg
    DMP_QVALUE_THRESHOLDS["CHG"] = args.q_chg
    DMP_QVALUE_THRESHOLDS["CHH"] = args.q_chh
    DMR_QVALUE_THRESHOLD = args.dmr_q
    VOTE_THRESHOLD = args.vote_threshold
    AUTO_DMP_VOTE_THRESHOLD = args.auto_dmp_vote_threshold
    AUTO_DMR_VOTE_THRESHOLD = args.auto_dmr_vote_threshold
    AUTO_VOTE_THRESHOLD_REPORT_ONLY = args.auto_vote_threshold_report_only
    DMR_ENGINE = args.dmr_engine
    AUTO_QVALUE_TWOSTEP = args.auto_qvalue_twostep
    AUTO_QVALUE_REPORT_ONLY = args.auto_qvalue_report_only
    AUTO_QVALUE_P_CUTOFF = args.auto_qvalue_p_cutoff
    AUTO_QVALUE_MIN_CANDIDATES = args.auto_qvalue_min_candidates
    AUTO_QVALUE_USE_SMOOTH = args.auto_qvalue_use_smooth
    AUTO_QVALUE_SMOOTH_SIGMA = args.auto_qvalue_smooth_sigma
    DMP_LOWDIFF_STRICT_VOTE = args.dmp_lowdiff_strict_vote
    DMP_LOWDIFF_CUTOFF = args.dmp_lowdiff_cutoff
    DMP_LOWDIFF_STRICT_VOTE_REPORT_ONLY = args.dmp_lowdiff_strict_vote_report_only

    # 验证目录存在，且确保组数合法
    if not os.path.exists(dir1):
        print(f"错误：目录 '{dir1}' 不存在！")
        sys.exit(1)
    if not os.path.exists(dir2):
        print(f"错误：目录 '{dir2}' 不存在！")
        sys.exit(1)
    if m <= 0 or n <= 0:
        print("错误：组数必须大于0")
        sys.exit(1)

    print(f"\n参数确认:")
    print(f"突变型目录: {dir1} (包含 {m} 组文件)")
    print(f"野生型目录: {dir2} (包含 {n} 组文件)")
    print(f"DMP甲基化差异阈值: {meth_diff_threshold} ({meth_diff_threshold * 100:.1f}%)")
    print(
        "mianjifa auto-methdiff: "
        f"启用判定={args.auto_meth_diff}, "
        f"仅报告={args.auto_meth_diff_report_only}, "
        f"cut_percent={args.auto_meth_diff_cut_percent}, "
        f"aggregate={args.auto_meth_diff_aggregate}, "
        f"fallback={args.auto_meth_diff_fallback}"
    )
    print(
        "DMP q-value阈值: "
        f"CpG={DMP_QVALUE_THRESHOLDS['CpG']}, "
        f"CHG={DMP_QVALUE_THRESHOLDS['CHG']}, "
        f"CHH={DMP_QVALUE_THRESHOLDS['CHH']}"
    )
    print(f"DMR q-value阈值: {DMR_QVALUE_THRESHOLD}")
    print(
        "两步法自动q-value阈值: "
        f"启用判定={AUTO_QVALUE_TWOSTEP}, "
        f"仅报告={AUTO_QVALUE_REPORT_ONLY}, "
        f"p_cutoff={AUTO_QVALUE_P_CUTOFF}, "
        f"min_candidates={AUTO_QVALUE_MIN_CANDIDATES}, "
        f"use_smooth={AUTO_QVALUE_USE_SMOOTH}"
    )
    print(f"最终DMP/DMR投票阈值: {VOTE_THRESHOLD} ({VOTE_THRESHOLD * 100:.1f}%)")
    print(
        "自动投票阈值: "
        f"DMP={AUTO_DMP_VOTE_THRESHOLD}, "
        f"DMR={AUTO_DMR_VOTE_THRESHOLD}, "
        f"report_only={AUTO_VOTE_THRESHOLD_REPORT_ONLY}"
    )
    print(
        "DMP lowdiff strict vote后处理: "
        f"启用判定={DMP_LOWDIFF_STRICT_VOTE}, "
        f"仅报告={DMP_LOWDIFF_STRICT_VOTE_REPORT_ONLY}, "
        f"cutoff={DMP_LOWDIFF_CUTOFF}"
    )
    if (DMP_LOWDIFF_STRICT_VOTE or DMP_LOWDIFF_STRICT_VOTE_REPORT_ONLY) and meth_diff_threshold > 0:
        print(
            "提示：--dmp-lowdiff-strict-vote 是 final DMP 后处理；"
            f"当前 --meth-diff={meth_diff_threshold:.6g} 会先在pair支持层过滤，"
            "lowdiff strict vote 将在该结果基础上继续应用。"
        )
    print(f"是否跳过DMR分析: {args.skip_dmr}")
    print(f"是否跳过滑动窗口/绘图: {args.skip_window}")
    print(f"并行 worker 进程数: {args.threads}")

    print("\n第一阶段：newtoboth进行中")
    # 将bismark新格式数据转换为both格式
    chr_series = newtoboth(m, n, dir1, dir2, threads=args.threads, work_dir=".")
    if biotype == 0:
        unfilter_mtypes = ["CHH", "CHG"]
    elif biotype == 1:
        unfilter_mtypes = ["CpG"]
    elif biotype == 2:
        unfilter_mtypes = ["CHH", "CHG", "CpG"]
    else:
        print("错误：生物类型必须是0、1或2")
        sys.exit(1)
    print(f"不需要p值预过滤的甲基化类型: {unfilter_mtypes}")
    success = process_all_combinations(
        dir1, dir2, m, n, unfilter_mtypes,
        meth_diff_threshold=meth_diff_threshold,
        skip_dmr=args.skip_dmr,
        skip_window=args.skip_window,
        threads=args.threads
    )  # process_all_combinations是进行m*n*3次检验

    if success:  # 全部检验都成功
        print("\n所有检验和FDR校正均成功完成！")

        if AUTO_QVALUE_TWOSTEP or AUTO_QVALUE_REPORT_ONLY:
            plot_all_auto_qvalue_panels(
                work_dir=".",
                m=m,
                n=n,
                unfilter_mtypes=unfilter_mtypes
            )

        if args.auto_meth_diff or args.auto_meth_diff_report_only:
            print(
                "\n提示：mianjifa auto-methdiff 读取完整pairwise FDR表中的q显著位点，"
                "不使用预先经过MethDiff过滤的DMP文件；估计出的阈值将继续用于"
                "auto-vote支持构建和final common DMP。"
            )
            auto_methdiff_threshold, auto_methdiff_summary = estimate_mianjifa_auto_methdiff_threshold(
                mut_dir=dir1,
                wt_dir=dir2,
                work_dir=".",
                cut_fraction=args.auto_meth_diff_cut_percent,
                fallback=args.auto_meth_diff_fallback,
                aggregate=args.auto_meth_diff_aggregate,
                output_dir=os.path.join(".", "and_output", "auto_methdiff_thresholds"),
                report_only=args.auto_meth_diff_report_only
            )
            if args.auto_meth_diff and not args.auto_meth_diff_report_only:
                meth_diff_threshold = auto_methdiff_threshold
                print(
                    "启用 mianjifa auto-methdiff：auto-vote与final common DMP "
                    f"统一使用 methdiff 阈值 = {meth_diff_threshold:.6g}"
                )
                regenerate_pair_dmp_outputs_from_fdr(
                    m=m,
                    n=n,
                    dir1=dir1,
                    dir2=dir2,
                    work_dir=".",
                    meth_diff_threshold=meth_diff_threshold,
                    contexts=("CpG", "CHH", "CHG"),
                )
            else:
                print(
                    f"mianjifa auto-methdiff report-only：估计阈值 = {auto_methdiff_threshold:.6g}；"
                    f"实际仍使用 --meth-diff = {meth_diff_threshold:.6g}"
                )

        methylation_types = ['CpG', 'CHH', 'CHG']
        if AUTO_DMP_VOTE_THRESHOLD or AUTO_VOTE_THRESHOLD_REPORT_ONLY:
            update_auto_dmp_vote_thresholds(
                m=m, n=n, work_dir=".",
                meth_diff_threshold=meth_diff_threshold
            )
        for mtype in methylation_types:
            common_sites_df = find_common_significant_sites(methytype2=mtype, dir1=dir1, dir2=dir2, meth_diff_threshold=meth_diff_threshold)
            if common_sites_df is not None and not common_sites_df.empty:
                print(f"找到 {len(common_sites_df)} 个 {mtype} 类型的共同显著位点")

                # 进行滑动窗口分析
                if args.skip_window:
                    print(f"跳过 {mtype} 共同显著位点滑动窗口分析 (--skip-window)")
                else:
                    print(f"开始对 {mtype} 共同显著位点进行滑动窗口分析...")
                    results = process_common_sites_sliding_window(
                        common_sites_df=common_sites_df,
                        methytype=mtype
                    )
                    print(f"完成 {mtype} 类型的滑动窗口分析")
            else:
                print(f"未找到 {mtype} 类型的共同显著位点")

        if args.skip_dmr:
            print("跳过最终 common DMR 分析流程 (--skip-dmr)")
        else:
            print("开始 DMR 分析流程")
            process_common_sites_dmr_and_summarize(
                dir1=dir1,
                dir2=dir2,
                m=m,
                n=n,
                methylation_types=methylation_types,
                threads=args.threads,
            )

        # 生成所有滑动窗口的可视化图表
        if args.skip_window:
            print("跳过所有滑动窗口可视化绘图 (--skip-window)")
        else:
            plot_methylation_sliding_windows(chr_series=chr_series)
            for mtype111 in ["CpG", "CHH", "CHG"]:
                plot_common_sites_sliding_windows(mtype111, chr_series=chr_series)
        convert_chromosome_to_names(chr_series=chr_series, work_dir=".")
        rename_chromosome_files(chr_series=chr_series, work_dir=".")
        convert_output_to_csv(work_dir=".")
        print(f"- DMP 结果：output_x_y/甲基化类型/")
        print(f"- 共同显著位点：and_output/")
        print(f"- 最终显著 DMR：and_output/*-final_significant_regions_DMRs.txt")
    else:
        print("\n部分检验失败，请检查输出信息。")
    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"总耗时: {elapsed_time:.2f} 秒")


if __name__ == "__main__":
    main()



