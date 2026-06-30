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

# DMP q-value significance threshold configuration
DMP_QVALUE_THRESHOLDS = {
    'CpG': 0.05,
    'CHH': 0.045,
    'CHG': 0.04,
}

# DMR q-value significance threshold
DMR_QVALUE_THRESHOLD = 0.05

# Final DMP/DMR voting threshold across replicate combinations; defaults to 2/3 for backward compatibility
VOTE_THRESHOLD = 2 / 3

# Automatically estimate DMP/DMR voting thresholds. Disabled by default so --vote-threshold is preserved when no new option is provided.
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
    """Return the DMP q-value threshold for the specified methylation type."""
    return DMP_QVALUE_THRESHOLDS.get(methylation_type, 0.05)


def calc_meth_diff(m1, u1, m2, u2):
    """Calculate the absolute methylation-level difference between two samples; the range is 0-1."""
    total1 = m1 + u1
    total2 = m2 + u2
    ratio1 = m1 / total1 if total1 > 0 else 0
    ratio2 = m2 / total2 if total2 > 0 else 0
    return abs(ratio1 - ratio2)


# =========================================================
# Mianjifa auto-methdiff module (embedded, no GUI)
# Adapted from mianjifa_automethdiff.py. It only:
#   1) estimates one global abs(methdiff) threshold from pairwise DMP distributions;
#   2) saves mxn MethDiff distribution plots for each methylation context.
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
        raise ValueError(f"Folder name does not start with output_: {name}")
    body = name[len("output_"):]
    m = re.match(r"^(.+?)(\d+)[_-](.+?)(\d+)$", body)
    if not m:
        raise ValueError(f"Cannot parse output folder name: {name}; expected a format like output_wt1_mut1")
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
    raise FileNotFoundError(f"{parent}  does not contain the original group folder: {dirname}")


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
        raise FileNotFoundError(f"{methylation_type_dir}  contains no DMP* files")
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
            f"{methylation_type_dir}  contains no FDR_corrected_results_*.txt"
        )
    if len(fdr_files) > 1:
        raise RuntimeError(
            f"{methylation_type_dir}  contains multiple FDR files; cannot determine a unique file: "
            f"{[p.name for p in fdr_files]}"
        )

    fdr_file = fdr_files[0]
    try:
        df = pd.read_csv(fdr_file, sep="\t")
        if len(df.columns) < 6:
            df = pd.read_csv(fdr_file, sep=r"\s+", engine="python")
    except Exception as exc:
        raise RuntimeError(f"Failed to read FDR file {fdr_file}: {exc}") from exc

    required = {"Chromosome", "Position", "Qvalue"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(
            f"{fdr_file} is missing required columns {sorted(missing)}; actual columns={list(df.columns)}"
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
        raise FileNotFoundError(f"{group_dir}  not found under  {rep_id}-*.txt")

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
    print(f"Saved auto-methdiff plot: {output_png}")
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
        print(f"Auto methdiff threshold estimation failed: no output_*_* folders found; falling back to {fallback}")
        print(f"Auto methdiff diagnostic table saved to: {out_file}")
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
        print(f"Auto methdiff threshold estimation failed: no methylation-context directories with complete FDR tables were found; falling back to {fallback}")
        print(f"Auto methdiff diagnostic table saved to: {out_file}")
        return float(fallback), summary_df

    print("\nEstimating the auto-methdiff threshold: mianjifa q-significant raw MethDiff area-cut method")
    print(f"  work_dir = {selected_dir}")
    print(f"  mut_dir = {mut_dir}")
    print(f"  wt_dir = {wt_dir}")
    print(f"  cut_fraction = {cut_fraction}")
    print(f"  aggregate = {aggregate}")

    all_summary = []
    for meth_type in meth_types:
        print(f"\nProcessing auto-methdiff methylation context: {meth_type}")
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

    print(f"Auto methdiff diagnostic table saved to: {summary_tsv}")
    print(f"Auto methdiff aggregate-threshold table saved to: {aggregate_tsv}")
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

    print(f"    Auto q-value threshold diagnostic table saved to: {out_file}")

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
        print(f"  [auto-q plot] read failed {filepath}: {e}")
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
        print("  [auto-q plot] no max-diff results to save")
        return None

    out_df = pd.DataFrame(rows)
    out_csv = os.path.join(save_dir, f"max_pvalue_adjustedP_diff_points_{len(rows)}comparisons.csv")
    out_df.to_csv(out_csv, index=False)

    print(f"  [auto-q plot] max-diff table saved to: {out_csv}")
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

    print(f"  [auto-q plot] FDR panel plot saved to : {save_path}")
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

    print(f"  [auto-q plot] diff panel plot saved to : {out_png}")
    return out_png


def plot_auto_qvalue_panels_for_context(work_dir, m, n, methylation_type, unfilter_mtypes):
    """
    Summarize FDR_corrected_results_*.txt files from all pairwise comparisons for one methylation context,
    and generate two plots with the same visual style as auto_qvalue_and_plot.py:
      1. p-value / adjusted p-value curve with the max-difference point;
      2. adjusted p-value minus p-value curve with the max-difference point.
    """
    if not (AUTO_QVALUE_TWOSTEP or AUTO_QVALUE_REPORT_ONLY):
        return

    # Use auto-qvalue only for two-step-FDR contexts; skip non-two-step contexts to avoid misinterpretation.
    if methylation_type in unfilter_mtypes:
        print(f"  [auto-q plot] {methylation_type} is not a two-step-FDR context; skipping auto-qvalue plotting")
        return

    input_files, ordered_keys, display_name = _autoq_collect_context_fdr_files(
        work_dir=work_dir,
        m=m,
        n=n,
        methylation_type=methylation_type
    )

    if not ordered_keys:
        print(f"  [auto-q plot] No {methylation_type} FDR_corrected_results file found; skipping")
        return

    save_dir = os.path.join(work_dir, "and_output", "auto_qvalue_plots", methylation_type)
    os.makedirs(save_dir, exist_ok=True)

    _autoq_setup_sci_style()
    color_pairs = _autoq_build_color_pairs(ordered_keys)
    diff_cache = {}

    for key in ordered_keys:
        rank, p, adj = _autoq_read_pq_file(input_files[key])
        if rank is None:
            print(f"  [auto-q plot] read failed; skipping: {input_files[key]}")
            continue

        rank_sub, p_sub, adj_sub = _autoq_keep_pvalue_below_cutoff(
            p,
            adj,
            cutoff=AUTO_QVALUE_P_CUTOFF
        )
        if rank_sub is None:
            print(f"  [auto-q plot] {display_name[key]} has no sites with p < {AUTO_QVALUE_P_CUTOFF} sites; skipping")
            continue

        info = _autoq_compute_max_diff_point(rank_sub, p_sub, adj_sub)
        if info is None:
            print(f"  [auto-q plot] {display_name[key]} max-diff calculation failed; skipping")
            continue

        diff_cache[key] = info

    if not diff_cache:
        print(f"  [auto-q plot] {methylation_type} has no plottable data")
        return

    print(f"  [auto-q plot] {methylation_type}: detected {len(diff_cache)} pairwise FDR files")
    _autoq_save_max_diff_table(diff_cache, ordered_keys, display_name, save_dir)
    _autoq_plot_fdr_panels(diff_cache, ordered_keys, display_name, color_pairs, save_dir, AUTO_QVALUE_P_CUTOFF)
    _autoq_plot_diff_panels(diff_cache, ordered_keys, display_name, color_pairs, save_dir)


def plot_all_auto_qvalue_panels(work_dir, m, n, unfilter_mtypes):
    if not (AUTO_QVALUE_TWOSTEP or AUTO_QVALUE_REPORT_ONLY):
        return

    print("\nGenerating auto q-value threshold diagnostic plots...")
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

    # For a unimodal distribution, GMM has no reliable valley; use the fallback integer threshold implied by 2/3 or the current vote-threshold.
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
        print(f"  WARNING: failed to save the auto vote-threshold distribution plot {out_file}: {e}")
        try:
            plt.close()
        except Exception:
            pass


def _write_auto_vote_summary(records, output_dir, filename):
    os.makedirs(output_dir, exist_ok=True)
    out_file = os.path.join(output_dir, filename)
    pd.DataFrame(records).to_csv(out_file, sep='\t', index=False)
    print(f"Auto vote-threshold summary table saved to: {out_file}")


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

    group_label = f"{m}x{n}"
    total_groups = int(m) * int(n)
    meth_diff_threshold = float(meth_diff_threshold)
    print(
        f"Analyzing DMP auto vote threshold: {m} rows x {n}  columns of comparison groups"
        f"(with {total_groups} groups)"
    )
    print(
        "  data source = completeFDRtable;pairsupport condition = "
        f"q passes and abs(MethDiff)>={meth_diff_threshold:.6g}"
    )

    group_dirs = []
    for mut_idx in range(1, int(m) + 1):
        for wt_idx in range(1, int(n) + 1):
            dir_path = base_path / f"output_wt{wt_idx}_mut{mut_idx}"
            if not dir_path.is_dir():
                print(f"  [WARN] WARNING: directory does not exist {dir_path}")
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
        print(f"\n> Processing DMP {ctx} context ...")
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
                print(f"  [WARN] {ctx_dir}  contains no complete FDR file")
                continue
            if len(fdr_files) > 1:
                print(
                    f"  [WARN] {ctx_dir}  contains multiple FDR files; using the first one: "
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
                print(f"  [WARN] Failed to read FDR file {fdr_file}: {exc}")
                continue

            required = {
                "Chromosome", "Position", "Pvalue", "Qvalue", "MethDiff"
            }
            missing = required - set(df.columns)
            if missing:
                print(
                    f"  [WARN] {fdr_file} is missing columns {sorted(missing)}; skipping"
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
                f"  No {ctx} sites,"
                "Using --vote-threshold fallback"
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
            print(f"  {ctx} formed no valid group columns; using the proportional threshold as fallback")
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
            print(f"  [WARN] {ctx} failed to save the auto vote-threshold plot {img_file}: {e}")
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
            f"  [OK] {ctx}: recommended threshold = {best_t};"
            f"MethDifffiltering threshold={meth_diff_threshold:.6g};"
            f"plot={img_file}"
        )

    _write_auto_vote_summary(
        records, str(plot_path), "DMP_vote_threshold_summary.tsv"
    )
    pd.DataFrame(pair_filter_records).to_csv(
        plot_path / "DMP_pair_filtering_summary.tsv",
        sep="\t", index=False,
    )
    print(
        "DMP pair-filtering summary table saved to: "
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
    Automatically read significant DMR files from m*n comparison groups under
    and_output/dmr_analysis_wt<wt>_mut<mut>/<context>/dmr_fisher_significant_*.txt.
    DMRs with qvalue <= DMR_QVALUE_THRESHOLD are merged, support counts are calculated,
    and a bimodal GMM is used to automatically select the DMR voting threshold.
    The DMR support-count distribution plot and DMR_vote_threshold_summary.tsv are also saved.

    Special cases:
      1. Unimodal distribution: use the fallback required_count implied by the current --vote-threshold.
      2. Truncated distribution: if the low-support range has no data, decrease the GMM threshold by 1.
    """
    if base_dir is None:
        base_path = Path(__file__).parent
    else:
        base_path = Path(base_dir)

    # Support two cases:
    # 1) base_dir is the run root, DMR files are under base_dir/and_output/dmr_analysis_...
    # 2) base_dir is already and_output, DMR files are under base_dir/dmr_analysis_...
    and_output_path = base_path / 'and_output'
    if not and_output_path.is_dir():
        and_output_path = base_path

    plot_path = Path(output_dir) if output_dir is not None else and_output_path / 'auto_vote_thresholds'
    plot_path.mkdir(parents=True, exist_ok=True)

    group_label = f"{m}x{n}"
    total_groups = int(m) * int(n)
    print(f"Analyzing DMR auto vote threshold: {m} rows x {n} columns of comparison groups (total {total_groups} groups)")

    group_dirs = []
    for mut_idx in range(1, int(m) + 1):
        for wt_idx in range(1, int(n) + 1):
            dir_path = and_output_path / f"dmr_analysis_wt{wt_idx}_mut{mut_idx}"
            if not dir_path.is_dir():
                print(f"  [WARN] WARNING: directory does not exist {dir_path}")
            group_dirs.append(dir_path)

    if chromosomes is None or chromosomes == "all":
        chrom_filter = None
    else:
        chrom_filter = {str(c).lower() for c in chromosomes}

    contexts = ["CpG", "CHH", "CHG"]
    thresholds = {}
    records = []

    for ctx in contexts:
        print(f"\n> Processing DMR {ctx} context ...")
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
                    print(f"  [WARN] WARNING: failed to read DMR file {file}: {e}")
                    continue

        if not all_records:
            print(f"  No significant {ctx} DMRs found; using --vote-threshold as fallback")
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
            print(f"  {ctx} formed no valid group columns; using --vote-threshold as fallback")
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

        # ========== Special handling: unimodal distribution ==========
        if len(nonzero_counts) == 1:
            single_value = int(nonzero_counts.index[0])
            best_t = int(fallback)
            valley_x = np.nan
            means = np.array([float(single_value), float(single_value)])
            status = 'fallback_single_support_count'
            method = 'fallback_vote_threshold'

            print(
                f"  [WARN] Detected a unimodal distribution: all DMR support counts are {single_value},"
                f"using fallback threshold t={best_t}"
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

        # ========== Normal bimodal fitting ==========
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

            # ========== Special handling: truncated distribution ==========
            half_limit = total_groups // 2
            low_range_zero = all(counts.loc[i] == 0 for i in range(1, half_limit + 1))

            if low_range_zero:
                original_t = int(best_t)
                best_t = max(1, int(best_t) - 1)
                status = 'ok_truncated_adjusted'
                method = f'{method}_truncated_adjust_{original_t}_to_{best_t}'
                print(
                    f"  [WARN] Detected a truncated distribution: support counts from 1 to {half_limit} have no DMR,"
                    f"adjust threshold from {original_t} to {best_t}"
                )

            best_t = max(1, min(int(best_t), total_groups))

        # ========== Plotting ==========
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
            print(f"  [WARN] WARNING: {ctx} failed to save the auto vote-threshold plot {img_file}: {e}")
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

        print(f"  [OK] {ctx}: recommended threshold = {best_t}; plot saved to  {img_file}")

    _write_auto_vote_summary(records, str(plot_path), 'DMR_vote_threshold_summary.tsv')
    return thresholds


def update_auto_dmp_vote_thresholds(m, n, work_dir='.', meth_diff_threshold=0.0):
    global AUTO_DMP_VOTE_THRESHOLDS
    if not (AUTO_DMP_VOTE_THRESHOLD or AUTO_VOTE_THRESHOLD_REPORT_ONLY):
        AUTO_DMP_VOTE_THRESHOLDS = {'CpG': None, 'CHH': None, 'CHG': None}
        print("DMP auto vote threshold is disabled; using --vote-threshold")
        return AUTO_DMP_VOTE_THRESHOLDS
    AUTO_DMP_VOTE_THRESHOLDS = compute_dmp_vote_thresholds(
        m=m, n=n, chromosomes='all', base_dir=work_dir,
        output_dir=os.path.join(work_dir, 'and_output', 'auto_vote_thresholds'),
        meth_diff_threshold=meth_diff_threshold
    )
    print(f"DMP auto vote thresholds: {AUTO_DMP_VOTE_THRESHOLDS}")
    return AUTO_DMP_VOTE_THRESHOLDS


def update_auto_dmr_vote_thresholds(m, n, work_dir='.'):
    global AUTO_DMR_VOTE_THRESHOLDS
    if not (AUTO_DMR_VOTE_THRESHOLD or AUTO_VOTE_THRESHOLD_REPORT_ONLY):
        AUTO_DMR_VOTE_THRESHOLDS = {'CpG': None, 'CHH': None, 'CHG': None}
        print("DMR auto vote threshold is disabled; using --vote-threshold")
        return AUTO_DMR_VOTE_THRESHOLDS
    AUTO_DMR_VOTE_THRESHOLDS = compute_dmr_vote_thresholds(
        m=m, n=n, chromosomes='all', base_dir=work_dir,
        output_dir=os.path.join(work_dir, 'and_output', 'auto_vote_thresholds')
    )
    print(f"DMR auto vote threshold: {AUTO_DMR_VOTE_THRESHOLDS}")
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
    """Worker for common DMR reads summation/Fisher/FDR for one pair x context task."""
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
        print(f"Processing file {filepath}")
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
    Full workflow: process common-site DMR candidates, summarize methylation reads, and aggregate results.

    Parallelization strategy:
      1. common significant sites -> common DMR candidates: parallelized by methylation type;
      2. common DMR read summarization/Fisher/FDR across replicate pairs: parallelized in summarize_all_dmr_methylation.
    """
    print("Stage 3: process common significant sites into DMR candidates")

    threads = max(1, int(threads))
    config = _build_parallel_config()

    if threads <= 1 or len(methylation_types) <= 1:
        for mtype in methylation_types:
            dmr_results = process_common_sites_to_dmr(methylation_type=mtype, work_dir=work_dir)
            if dmr_results:
                print(f"\n{mtype} DMR analysis completed; {len(dmr_results)} chromosomes have DMR results")
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

        print(f"Enabled parallel common-DMR candidate processing: workers={max_workers}, tasks={len(tasks)}")
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

    # Step 2: summarize DMR methylation over all output_x_y comparisons and aggregate the results
    summarize_all_dmr_methylation(dir1, dir2, m, n, methylation_types, work_dir=work_dir, threads=threads)

def generate_final_significant_dmr(dmr_data_dict, methylation_type, output_dir,m,n,dir1,dir2):
    """
    Generate the final significant DMR file based on the voting/Bayesian decision rule.

    Parameters:
        dmr_data_dict: mapping from chromosome to DMR interval (start, end) to a data list of per-site statistics.
        methylation_type: methylation context.
        output_dir: output directory.
        threshold: voting/Bayesian decision threshold; default is 2/3.
    """
    print(f"\nProcessing {methylation_type}  final DMR aggregation...")

    if not dmr_data_dict:
        print(f"  {methylation_type} has no DMR data")
        return

    auto_vote_threshold = AUTO_DMR_VOTE_THRESHOLDS.get(methylation_type)
    if AUTO_VOTE_THRESHOLD_REPORT_ONLY:
        auto_vote_threshold = None
    if auto_vote_threshold is not None:
        print(f"  Using {methylation_type} DMR auto vote threshold: {auto_vote_threshold}")
    else:
        print(f"  {methylation_type} DMR auto vote threshold is unavailable or report-only; using proportional threshold {VOTE_THRESHOLD}")

    final_dmr_list = []

    # Iterate over all chromosomes
    for chr_num in sorted(dmr_data_dict.keys()):
        print(f"  Processing chromosome {chr_num}...")
        chr_dmr_dict = dmr_data_dict[chr_num]  # Get the mapping for the current chromosome:DMR interval (start, end) -> data list (per-site statistics)
                                                                    # The data list contains m*n elements, one for each test for this region

        # Iterate over all DMR regions on this chromosome
        for dmr_key, data_list in chr_dmr_dict.items():
            start, end = dmr_key

            # Count significant calls
            sig_count = sum(1 for item in data_list if item.qvalue <= DMR_QVALUE_THRESHOLD)
            total_count = len(data_list)

            # Decision
            is_significant = bayes_deciding(sig_count, total_count - sig_count, auto_vote_threshold=auto_vote_threshold)

            if not is_significant:
                continue

            # Calculate mean values
            avg_exp_m = np.mean([item.exp_methy for item in data_list])
            avg_exp_u = np.mean([item.exp_unmethy for item in data_list])
            avg_wild_m = np.mean([item.wild_methy for item in data_list])
            avg_wild_u = np.mean([item.wild_unmethy for item in data_list])

            filtered_values = [item.qvalue for item in data_list if item.qvalue <= DMR_QVALUE_THRESHOLD]
            avg_qvalue = np.mean(filtered_values) if filtered_values else 1

            # Determine direction by majority voting
            direction_votes = [item.direction for item in data_list]
            direction = 1 if sum(direction_votes) >= len(direction_votes) / 2 else 0

            # Calculate probability
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
            # Add all replicate q-values in order
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
        print(f"  {methylation_type} has no significant DMRs")
        return

    # Convert to DataFrame and sort
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
    # Save results
    output_file = os.path.join(output_dir, f"{methylation_type}-final_significant_regions_DMRs.txt")
    final_df.to_csv(output_file, sep='\t', index=False)

    print(f"  {methylation_type} final significant DMRs: {len(final_df)} ")
    print(f"  saved to: {output_file}")

    # Summary statistics
    hyper_count = len(final_df[final_df['Direction'] == 1])
    hypo_count = len(final_df[final_df['Direction'] == 0])
    print(f"    - hypermethylated: {hyper_count} ({hyper_count / len(final_df) * 100:.1f}%)")
    print(f"    - hypomethylated: {hypo_count} ({hypo_count / len(final_df) * 100:.1f}%)")

    return final_df

def collect_dmr_results(methy_dir, methylation_type, all_dmr_results):
    """
    Collect results from one DMR analysis.

    Parameters:
        methy_dir: methylation-type directory, e.g., ./and_output/CpG/.
        methylation_type: methylation context.
        all_dmr_results: aggregation dictionary initialized with the three methylation-type keys.
    """
    # Find all dmr_fisher_Chr*.txt files
    fisher_files = glob.glob(os.path.join(methy_dir, "dmr_fisher_Chr*.txt"))

    for fisher_file in fisher_files:
        # Extract chromosome number
        match = re.search(r'Chr(\d+)\.txt$', fisher_file)
        if not match:
            continue
        chr_num = int(match.group(1))

        # Read file
        try:
            df = pd.read_csv(fisher_file, sep='\t')

            if df.empty:
                continue

            # Initialize the dictionary for this chromosome
            if chr_num not in all_dmr_results[methylation_type]:
                all_dmr_results[methylation_type][chr_num] = defaultdict(list)

            # Collect data for each DMR region
            for _, row in df.iterrows():
                dmr_key = (int(row['DMR_start']), int(row['DMR_end']))

                # Store (exp_m, exp_u, wild_m, wild_u, qvalue, direction)
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
            print(f"    WARNING: read {fisher_file} failed: {e}")
            continue

def summarize_all_dmr_methylation(dir1, dir2, m, n, methylation_types=['CpG', 'CHH', 'CHG'], work_dir=".", threads=1):
    """
    For all output_x_y directories, summarize methylation reads over common DMR regions.
    In the parallel version, each worker writes only to its own dmr_analysis_wt*_mut*/context directory.
    The main process still aggregates all_dmr_results sequentially to avoid concurrent writes to a shared dictionary.
    """
    print("Stage 4: sum methylation reads for common DMRs across all combinations")

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
                    print(f"  Skipping wt{replicate_y}_mut{replicate_x} {mtype}: both-format file does not exist")
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
            print(f"\nProcessing combination ({dir1}{task['replicate_x']}, {dir2}{task['replicate_y']}) {task['methylation_type']}...")

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
                print(f"  ERROR: failed to process {task['methylation_type']} failed: {e}")
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

        print(f"Enabled parallel common-DMR read-sum/Fisher/FDR processing: workers={max_workers}, tasks={len(tasks)}")
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

    # The main process collects results centrally to avoid concurrent modification of all_dmr_results.
    # Note: parallel task completion order is nondeterministic; collecting directly by completion order
    # would append q-values to data_list in a random order, causing qvalue_* columns in the final DMR file to be misaligned.
    # Therefore restore the original serial order here:replicate_x -> replicate_y -> methylation_type.
    completed = sorted(
        completed,
        key=lambda t: (int(t["replicate_x"]), int(t["replicate_y"]), str(t["methylation_type"]))
    )

    print("Collecting Fisher/FDR results for common DMRs across all combinations...")
    for task in completed:
        collect_dmr_results(task["methy_output_dir"], task["methylation_type"], all_dmr_results)

    if AUTO_DMR_VOTE_THRESHOLD or AUTO_VOTE_THRESHOLD_REPORT_ONLY:
        update_auto_dmr_vote_thresholds(m=m, n=n, work_dir=work_dir)

    print("Stage 5: aggregate DMR results and apply the Bayesian decision rule")

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
    Read final_significant_sites_DMPs.txt and process DMRs by chromosome.

    Parameters:
        methylation_type: methylation context (CpG, CHH, CHG).

    Returns:
        dmr_results: dictionary in the form {chr_num: dmr_list_file_path}.
    """
    print(f"\nStarting {methylation_type}  common significant-site DMR analysis...")

    and_output_dir = os.path.join(work_dir, "and_output")
    os.makedirs(and_output_dir, exist_ok=True)

    # 1. Read the final_significant_sites_DMPs file
    common_file = os.path.join(and_output_dir, f"{methylation_type}-final_significant_sites_DMPs.txt")
    if not os.path.exists(common_file):
        print(f"ERROR: file does not exist {common_file}")
        return None

    try:
        df = pd.read_csv(common_file, sep='\t')
        print(f"Successfully read file; {len(df)}  sites")
    except Exception as e:
        print(f"Failed to read file: {e}")
        return None

    if df.empty:
        print(f"WARNING: {methylation_type} file is empty")
        return None

    # 2. Group by chromosome
    chromosomes = sorted(df['Chromosome'].unique(), key=natural_sort_key)
    print(f"Found {len(chromosomes)}  chromosomes: {chromosomes}")

    dmr_results = {}
    total_chromosomes = len(chromosomes)

    # 3. Process each chromosome
    for chr_num in chromosomes:
        print(f"\n  Processing chromosome {chr_num} ({chr_num}/{total_chromosomes})...")

        # Filter records for the current chromosome
        chr_df = df[df['Chromosome'] == chr_num].copy()

        # Extract the required three columns:Position, Sig_Mean_Qvalue, Methylation_Change
        dmp_data = chr_df[['Position', 'Sig_Mean_Qvalue', 'Methylation_Change']].copy()

        # Ensure correct data types
        dmp_data['Position'] = dmp_data['Position'].astype(int)
        dmp_data['Sig_Mean_Qvalue'] = dmp_data['Sig_Mean_Qvalue'].astype(float)
        dmp_data['Methylation_Change'] = dmp_data['Methylation_Change'].astype(int)

        # Sort by position
        dmp_data = dmp_data.sort_values('Position').reset_index(drop=True)

        print(f"    chromosome {chr_num} with {len(dmp_data)}  sites")

        if len(dmp_data) == 0:
            print(f"    skipping:chromosome {chr_num} has no valid sites")
            continue

        # 4. Create a temporary DMP file (format: pos qvalue change)
        temp_dmp_file = os.path.join(and_output_dir,
                                     f"DMP_common_{methylation_type}_Chr{chr_num}.txt")

        # Write the DMP-format file (first line is "first line", followed by pos qvalue change)
        with open(temp_dmp_file, 'w') as f:
            f.write("first line\n")
            for _, row in dmp_data.iterrows():
                f.write(f"{int(float(row['Position']))} {float(row['Sig_Mean_Qvalue'])} {int(float(row['Methylation_Change']))}\n")

        print(f"    Created DMP file: {os.path.basename(temp_dmp_file)}")

        # 5. Call the DMR analysis function
        # Note:run_dmr_pipeline_on_dmp_file requires the chromoNo parameter
        # Pass the total chromosome count
        try:
            dmr_list_file = run_dmr_pipeline_on_dmp_file_auto(
                dmp_file=temp_dmp_file,
                chromoNo=total_chromosomes
            )

            if dmr_list_file:
                dmr_results[chr_num] = dmr_list_file
                print(f"     chromosome {chr_num} DMR analysis completed")
            else:
                print(f"    chromosome {chr_num} DMR analysis failed (possibly no valid DMRs)")

        except Exception as e:
            print(f"     chromosome {chr_num} processing error: {e}")
            import traceback
            traceback.print_exc()
            continue

    print(f"\n{methylation_type} common significant-site DMR analysis completed.")
    print(f"Successfully processed {len(dmr_results)}/{len(chromosomes)}  chromosomes")

    return dmr_results

def summarize_dmr_methylation(methy_dir, replicate_x, replicate_y, file1_path, file2_path, methylation_type='CpG', custom_dmr_dir=None):
    """
    Sum methylated/unmethylated reads over DMR regions and calculate Fisher p-values,
    FDR q-values, and methylation direction. All output files are saved in methy_dir,
    for example ./output_1_1/CpG/.
    custom_dmr_dir specifies a custom DMR-file directory. If it is None, methy_dir is used.
    This separates DMRs generated inside each output_x_y directory from DMRs generated from common sites.
    """
    # This block handles only one test
    print(f"    Starting DMR methylation-read summation...")

    n_chromosomes = get_column_count(file1_path)
    if n_chromosomes is None:
        print("    Unable to obtain the chromosome count; skipping DMR summation")
        return

    chromosomes = [f'Chr{i}' for i in range(1, n_chromosomes + 1)]

    # Collect DMR data for all chromosomes without p-values
    all_dmr_data = []  # (chrom, start, end, exp_m, exp_u, wild_m, wild_u)

    for idx, chrom in enumerate(chromosomes):
        chrom_num = idx + 1
        if custom_dmr_dir is not None:
            # Read common DMR file
            dmr_file = os.path.join(custom_dmr_dir, f"DMR_list_DMP_common_{methylation_type}_Chr{chrom_num}.txt")
        else:
            # Read the DMR file from the output directory
            dmr_file = os.path.join(methy_dir, f"DMR_list_DMP_wt_replicate{replicate_y}_mut_replicate{replicate_x}_Chr{chrom_num}.txt")
        if not os.path.exists(dmr_file):
            print(f"      skipping {chrom}; DMR file does not exist")
            continue

        # Read DMR regions
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
            print(f"      {chrom} has no valid DMR regions")
            continue

        col_start = idx * 3

        # Read experimental-group data
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
            print(f"      Failed to read experimental-group data ({chrom}): {e}")
            continue

        # Read control-group data
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
            print(f"      Failed to read control-group data ({chrom}): {e}")
            continue

        # Sum reads for each DMR region
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
        print("    No valid DMR data; skipping downstream analysis")
        return

    # === Step 1: write dmr_summary_{chrom}.txt ===
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
            print(f"      {chrom} DMR summary completed -> {summary_file}")

    # === Step 2: calculate p-values for each DMR ===
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

    # === Step 3: global FDR correction across chromosomes ===
    pvalues = np.array([item[7] for item in dmr_with_pvals])
    qvalues = calculate_qvalues(pvalues, pi=1.0)

    # === Step 4: organize results by chromosome and add direction ===
    chrom_data_dict = defaultdict(list)
    for i, (chrom, start, end, exp_m, exp_u, wild_m, wild_u, pval) in enumerate(dmr_with_pvals):
        qval = qvalues[i]

        # Calculate methylation direction
        exp_total = exp_m + exp_u
        wild_total = wild_m + wild_u
        if exp_total > 0 and wild_total > 0:
            exp_rate = exp_m / exp_total
            wild_rate = wild_m / wild_total
            direction = 1 if exp_rate > wild_rate else 0
        else:
            direction = 0  # Set to 0 when the direction cannot be determined

        chrom_data_dict[chrom].append((start, end, exp_m, exp_u, wild_m, wild_u, pval, qval, direction))

    # === Step 5: write full Fisher results and significant subset ===
    for chrom in chromosomes:
        if chrom not in chrom_data_dict:
            continue

        # Full results
        fisher_file = os.path.join(methy_dir, f"dmr_fisher_{chrom}.txt")
        with open(fisher_file, 'w') as f:
            f.write("DMR_start\tDMR_end\texp_methy_sum\texp_unmethy_sum\twild_methy_sum\twild_unmethy_sum\tpvalue\tqvalue\tdirection\n")
            for row in chrom_data_dict[chrom]:
                start, end, exp_m, exp_u, wild_m, wild_u, pval, qval, direction = row
                p_out = pval if not np.isnan(pval) else 'nan'
                q_out = qval if not np.isnan(qval) else 'nan'
                f.write(f"{start}\t{end}\t{exp_m}\t{exp_u}\t{wild_m}\t{wild_u}\t{p_out:.6g}\t{q_out:.6g}\t{direction}\n")
        print(f"      {chrom} Fisher + FDR + direction completed -> {fisher_file}")

        # Significant results using the configurable DMR q-value threshold
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
            print(f"        -> significant DMR (q<={DMR_QVALUE_THRESHOLD}): {sig_file}")
        else:
            print(f"        -> {chrom} has no significant DMRs (q<={DMR_QVALUE_THRESHOLD})")

    print("    DMR methylation analysis completed.")

def run_dmr_pipeline_on_dmp_file(dmp_file: str, chromoNo: int = 10):
    """
    Generate DMR regions and the dmr_list file from a DMP file.
    - methylation_matrix_file is used to dynamically obtain the chromosome count, e.g., 1-bothMeUnme_...txt.
    - All output files are saved in the directory containing dmp_file.
    """
    # Core idea:
    # Use sliding windows to find DMP-dense regions
    # Connect adjacent dense regions by jump merging
    # Finally retain DMRs that are long enough and sufficiently significant
    sWinN = 1000  # Sliding-window size(1000 bp)
    M0 = 4  # Minimum number of DMPs in a window(at least 4)
    M1 = 10  # Minimum number of DMPs in the final DMR(at least 10)
    M2 = 10  # Jump step size, preserving the original segmentation behavior,

    # For safety, ensure chromoNo >= 6 because arrayMethy1_script1[5] is used
    chromoNo = max(chromoNo, 6)

    class PositionNoNode:
        def __init__(self, pos=0, end=0, pV=0.0, ratio=0.0, num=0, num2=0, numCom=0, markR=0, DMR_S=0, DMR_E=0):
            self.pos = pos  # Window start position
            self.end = end  # Window end position
            self.pV = pV
            self.ratio = ratio
            self.num = num  # Number of hypermethylated sites
            self.num2 = num2  # Number of hypomethylated sites
            self.numCom = numCom
            self.markR = markR
            self.DMR_S = DMR_S
            self.DMR_E = DMR_E
            self.posV = []
            self.logPV = []
            self.meUnV = []

    output_dir = os.path.dirname(dmp_file)
    base_name = os.path.basename(dmp_file)

    arrayMethy1 = [[] for _ in range(chromoNo)]  # Create one sublist per chromosome

    # === Read DMP sites ===
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
                    tmpNode.num = change  # The lines above build site information for the current DMP row
                    arrayMethy1[0].append(tmpNode)
                except (ValueError, IndexError):
                    continue
    except Exception as e:
        print(f"    Failed to read DMP file {dmp_file}: {e}")
        return None

    if not arrayMethy1[0]:
        return None

    arrayMethy1[0].sort(key=lambda x: x.pos)
    firstP = arrayMethy1[0][0].pos
    lastP = arrayMethy1[0][-1].pos

    # === Sliding window ===
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

    # === Build the sliding-window list ===
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

    # === Standardized output ===
    if arrayMethy1[2]:
        maxCom = max(node.numCom for node in arrayMethy1[2])
        maxCom = max(maxCom, 1)
        std_file = os.path.join(output_dir, f"noTitle_allDMCs_new_Standardized_slidingW_{base_name}")
        with open(std_file, 'w') as f:
            for node in arrayMethy1[2]:
                if node.end <= lastP:
                    std_val = node.numCom / maxCom
                    f.write(f"{node.pos:<20}{node.end:<20}{node.numCom:<20}{std_val:.6f}\n")

    # === Jump merging using the provided code===
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

        # Extend leftward
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

        # Extend rightward
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
    print(f"Identified {len(arrayMethy1[3])}  DMR regions")

    # Write DMR results and boundary files
    with open(dmr_out_file, 'w') as cout05:
        with open(boundary_file, 'w') as bound_out:
            for node in arrayMethy1[3]:
                cout05.write(f"{node.pos:<20}{node.end:<20}{node.numCom:<20}[{node.DMR_S} {node.DMR_E}]\n")
                bound_out.write(f"{node.DMR_S} {node.DMR_E}\n")

    print(f"Generated boundary file: {boundary_file}")

    # === Merge overlapping boundaries with dynamic chromoL===
    chromoL = lastP + 100000  # Dynamic genome length
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

    print(f"    DMR analysis completed: {base_name} -> {final_file}")
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
        print(f"    C++ DMR input file does not exist: {dmp_abs}")
        return None

    # Match the Python behavior for empty / header-only DMP files.
    try:
        with open(dmp_abs, "r") as f:
            lines = f.readlines()
        valid_data_lines = [line for line in lines[1:] if line.strip()]
        if len(lines) < 2 or not valid_data_lines:
            return None
    except Exception as e:
        print(f"    C++ DMR failed to read DMP file {dmp_abs}: {e}")
        return None

    try:
        step1_bin = resolve_cpp_dmr_binary("dmr_step1")
        step2_bin = resolve_cpp_dmr_binary("dmr_step2_dynamic")
    except Exception as e:
        print(f"    Failed to locate C++ DMR executable: {e}")
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

    print(f"    Processing with the C++ DMR engine: {base_name}")

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

    print(f"    C++ DMR analysis completed: {base_name} -> {expected_list_file}")
    return expected_list_file


def run_dmr_pipeline_on_dmp_file_auto(dmp_file: str, chromoNo: int = 10):
    """
    Dispatch DMP -> DMR candidate/list generation according to --dmr-engine.
    """
    if DMR_ENGINE == "cpp":
        return run_dmr_pipeline_on_dmp_file_cpp(dmp_file, chromoNo=chromoNo)
    return run_dmr_pipeline_on_dmp_file(dmp_file, chromoNo=chromoNo)


def process_chr_in_one_file(df):
    """Add the chr prefix to a single input file when needed and return all chromosome information from that file."""

    # Keep chromosome labels that already start with chr unchanged; otherwise add the chr prefix
    def to_chr(chr_str):
        chr_str = str(chr_str).strip()
        if chr_str.lower().startswith('chr'):
            return f"chr{chr_str[3:]}"
        else:
            return f"chr{chr_str}"

    # Ensure all chromosome labels start with chr
    df['Chromosome_Label'] = df['Chromosome_Label'].apply(to_chr)
    # Return the unique set of chromosome labels
    return df['Chromosome_Label'].unique()  # Return a numpy.array

def natural_sort_key(chr_name):
    """
    Natural chromosome sorting key.
    Numeric chromosomes are sorted by numeric value, and alphabetic chromosomes such as X, Y, and M are sorted afterward alphabetically.
    """
    chr_name = str(chr_name).lower()
    # Remove 'chr' prefix
    if chr_name.startswith('chr'):
        suffix = chr_name[3:]
    else:
        suffix = chr_name

    # Try converting to integer
    try:
        # If numeric, return (0, numeric value, '')
        return (0, int(suffix), '')
    except ValueError:
        # If alphabetic(such as X, Y, or M), return (1, 0, letter)
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
            print(f"WARNING: file does not exist {filepath}")
            continue

        try:
            for chunk in pd.read_csv(
                filepath,
                sep=r"\s+",
                header=None,
                usecols=[0],
                names=["Chromosome_Label"],
                dtype={"Chromosome_Label": "string"},
                chunksize=chunksize,
            ):
                chrs = _normalize_chr_series(chunk["Chromosome_Label"])
                all_chromosomes.update(chrs.dropna().unique().tolist())

            print(f"Chromosome scan completed: {filepath}")

        except Exception as e:
            print(f"Error while reading file {filepath}:  {e}")

    unique_chrs = sorted(all_chromosomes, key=natural_sort_key)
    chr_series = pd.Series(range(len(unique_chrs)), index=unique_chrs)

    print(f"Unified chromosome mapping: {chr_series}")
    return chr_series

def scan_all_files_for_chr_mapping(m, n, dir1, dir2):
    """Scan all m+n files and collect chromosome information to create a unified chromosome mapping."""

    all_chromosomes = set()  # Use set uniqueness to store all chromosome labels

    # Scan the m files in the first directory; file-name format is i-dir1.txt,e.g., 3-msv.txt
    for i in range(1, m + 1):
        filepath = os.path.join(dir1, f"{i}-{os.path.basename(dir1)}.txt")
        if not os.path.exists(filepath):
            print(f"WARNING: file does not exist {filepath}")
            continue

        try:
            df = pd.read_csv(filepath,
                             sep=r'\s+',
                             header=None,
                             names=['Chromosome_Label', 'Site_Position', 'Methylated_Reads', 'Unmethylated_Reads', 'Methylation_Context'],
                             dtype=str)
            chromosomes = process_chr_in_one_file(df)  # Get the unique chromosome labels in this file
            all_chromosomes.update(chromosomes)  # Add all unique chromosome labels from this file to all_chromosomes
            # Note: update accepts any iterable and adds its elements to the set
            print(f"file {filepath} contains chromosomes: {sorted(chromosomes, key=natural_sort_key)}")
        except Exception as e:
            print(f"Error while reading file {filepath}:  {e}")

    # Scan files in the second directory
    for j in range(1, n + 1):
        filepath = os.path.join(dir2, f"{j}-{os.path.basename(dir2)}.txt")
        if not os.path.exists(filepath):
            print(f"WARNING: file does not exist {filepath}")
            continue

        try:
            df = pd.read_csv(filepath,
                             sep=r'\s+',
                             header=None,
                             names=['Chromosome_Label', 'Site_Position', 'Methylated_Reads', 'Unmethylated_Reads', 'Methylation_Context'],
                             dtype=str)
            chromosomes = process_chr_in_one_file(df)  # Same as above
            all_chromosomes.update(chromosomes)  # Same as above
            print(f"file {filepath} contains chromosomes: {sorted(chromosomes, key=natural_sort_key)}")
        except Exception as e:
            print(f"Error while reading file {filepath}:  {e}")

    # Create a unified chromosome mapping,
    # At this point, all_chromosomes contains all unique chromosome labels across files in both genotype directories
    unique_chrs = sorted(all_chromosomes, key=natural_sort_key)
    # Build a sorted Series mapping each unique chromosome label across both genotype directories to a numeric index
    chr_series = pd.Series(range(len(unique_chrs)), index=unique_chrs)

    print(f"Unified chromosome mapping: {chr_series}")
    return chr_series

# def single_newtoboth(filepath1, output_dir, num1, chr_series):
    """
    Parameters:
        filepath1 is a new-format file such as 1-wt.txt.
        output_dir is the directory where both-format output files are written, usually the directory containing filepath1.
        num1 is the index of the new-format file currently being processed, i.e., num1-genotype.txt.
        chr_series is the sorted Series mapping all unique chromosome labels to numeric indices.
    """

    df = pd.read_csv(filepath1,
                     sep=r'\s+',
                     header=None,
                     names=['Chromosome_Label', 'Site_Position', 'Methylated_Reads', 'Unmethylated_Reads', 'Methylation_Context'],
                     dtype=str)

    def to_chr(chr_str):
        chr_str = str(chr_str).strip()
        if chr_str.lower().startswith('chr'):
            return f"chr{chr_str[3:]}"
        else:
            return f"chr{chr_str}"

    # Because chr_series creation did not modify source files, labels read here may still lack the chr prefix and need normalization
    df['Chromosome_Label'] = df['Chromosome_Label'].apply(to_chr)
    df['Chromosome_Label'] = df['Chromosome_Label'].map(chr_series)  # Map each chromosome label to its numeric index, e.g., chr1 -> 0 and chr2 -> 1
    # The following three lines convert position and read-count columns to numeric values
    df['Site_Position'] = pd.to_numeric(df['Site_Position'])
    df['Methylated_Reads'] = pd.to_numeric(df['Methylated_Reads'])
    df['Unmethylated_Reads'] = pd.to_numeric(df['Unmethylated_Reads'])
    # Map chromosome labels to chr_series indices, convert position/read counts to numeric values, and convert methylation context CG to CpG
    df['Methylation_Context'] = df['Methylation_Context'].str.replace('CG', 'CpG')
    # Group the three methylation-context datasets for separate processing
    data_groups = df.groupby('Methylation_Context')
    # Iterate over group keys and their corresponding sub-dataframes
    for methy_type, data_ind in data_groups:
        if data_ind.empty:
            continue
        # Get and sort chromosome indices that actually exist for the current methylation type,
        # These are zero-based converted chromosome indices present for the current methylation type
        actual_chrs = sorted(data_ind['Chromosome_Label'].dropna().unique())
        # Also get the total chromosome count regardless of whether this methylation type has data for every chromosome
        chr_count = len(chr_series)
        chr_data_dict = {}
        # mlen stores the maximum number of rows among chromosomes for the current methylation type
        mlen = 0

        for chr_num in actual_chrs:  # Iterate over each zero-based chromosome index for this methylation type
            # Sort all data for the current chromosome by position in ascending order and store in chr_data
            chr_data = data_ind[data_ind['Chromosome_Label'] == chr_num].sort_values('Site_Position').reset_index(drop=True)
            # Store position, methylated reads, and unmethylated reads in chr_data_dict keyed by chromosome index
            # Each key maps to a numpy array containing these three attributes
            chr_data_dict[chr_num] = chr_data[['Site_Position', 'Methylated_Reads', 'Unmethylated_Reads']].values
            # Update mlen if a chromosome has more rows
            mlen = max(mlen, len(chr_data_dict[chr_num]))

        # Create the output matrix: rows are the maximum chromosome row count and columns are three times the number of unique chromosomes across both directories; initialize with zeros
        output_matrix = np.zeros((mlen, chr_count * 3), dtype=np.int32)

        # Iterate over all rows of the output matrix
        for i in range(mlen):
            # Iterate over zero-based chromosome indices present for the current methylation type
            for chr_num in actual_chrs:  # Use a continuous index
                col_start = chr_num * 3  # Calculate column position from the continuous index
                if i < len(chr_data_dict[chr_num]):
                    output_matrix[i, col_start:col_start + 3] = chr_data_dict[chr_num][i]

        # Output via pandas
        output_df = pd.DataFrame(output_matrix)  # Convert the output matrix to a DataFrame for output
        output_file = f"{num1}-bothMeUnme_diffChromo_NOREPEATED_methy_sites_{methy_type}.txt"
        output_path = os.path.join(output_dir, output_file)
        output_df.to_csv(output_path, sep='\t', header=False, index=False)

def single_newtoboth(filepath1, output_dir, num1, chr_series):
    """
    Accelerated single-file conversion for newtoboth.
    The original output format is preserved:
    each methylation type outputs one bothMeUnme_diffChromo_NOREPEATED file,
    and each chromosome occupies three columns: Position, methylated reads, and unmethylated reads.
    """

    col_names = [
        "Chromosome_Label",
        "Site_Position",
        "Methylated_Reads",
        "Unmethylated_Reads",
        "Methylation_Context",
    ]

    df = pd.read_csv(
        filepath1,
        sep=r"\s+",
        header=None,
        names=col_names,
        dtype={
            "Chromosome_Label": "string",
            "Site_Position": "int64",
            "Methylated_Reads": "int32",
            "Unmethylated_Reads": "int32",
            "Methylation_Context": "string",
        },
    )

    chr_map = chr_series.to_dict()

    df["Chromosome_Label"] = _normalize_chr_series(df["Chromosome_Label"]).map(chr_map)
    df = df.dropna(subset=["Chromosome_Label"]).copy()
    df["Chromosome_Label"] = df["Chromosome_Label"].astype(np.int32)

    # Convert only true CG contexts to CpG to avoid unintended replacement in other strings
    df["Methylation_Context"] = df["Methylation_Context"].replace({
        "CG": "CpG",
        "cg": "CpG",
        "cG": "CpG",
        "Cg": "CpG",
    })

    chr_count = len(chr_series)
    out_dtype = np.int32

    for methy_type, data_ind in df.groupby("Methylation_Context", sort=False):
        if data_ind.empty:
            continue

        # Sort once to avoid sorting each chromosome separately
        data_ind = data_ind.sort_values(
            ["Chromosome_Label", "Site_Position"],
            kind="mergesort",
        )

        chr_blocks = []
        mlen = 0

        for chr_num, chr_df in data_ind.groupby("Chromosome_Label", sort=True):
            arr = chr_df[
                ["Site_Position", "Methylated_Reads", "Unmethylated_Reads"]
            ].to_numpy(dtype=out_dtype, copy=True)

            chr_num = int(chr_num)
            chr_blocks.append((chr_num, arr))
            if arr.shape[0] > mlen:
                mlen = arr.shape[0]

        if mlen == 0:
            continue

        output_matrix = np.zeros((mlen, chr_count * 3), dtype=out_dtype)

        # Key acceleration: assign whole chromosome blocks rather than row by row and chromosome by chromosome
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
    Get chromosome name from chromosome index.

    Parameters:
        chr_num: ordinal chromosome index, starting from 1.
        chr_series: chromosome mapping Series.

    Returns:
        Chromosome name string, e.g., 'chr1' or 'chrX'.
    """
    chr_num = int(chr_num)
    if chr_series is not None:
        try:
            # chr_num is 1-based, so subtract 1
            if 0 <= chr_num - 1 < len(chr_series):
                return chr_series.index[chr_num - 1]
        except Exception as e:
            print(f"WARNING: failed to get chromosome name: {e}")

    # If an error occurs or no mapping exists, return the numeric index
    return f"chr{chr_num}"

def newtoboth(m, n, dir1, dir2, threads=1, work_dir="."):
    """Convert raw methylation files to bothMeUnme matrix files.

    Parallel design:
      1. Keep chromosome mapping construction serial to guarantee one shared chr_series.
      2. If threads > 1, convert the m+n independent raw files concurrently.
         Each worker writes only its own replicate-numbered output files, so filenames do not collide.
    """
    # No need to check directory existence here because main already checked it
    # Get the sorted Series mapping all unique chromosome labels across the two genotype directories to numeric indices
    chr_series = scan_all_files_for_chr_mapping_fast(m, n, dir1, dir2)
    print(f"Mapping: {chr_series}")

    threads = max(1, int(threads))
    tasks = []

    # Loop m+n times to process files from both genotype directories. Keep the original task order: all dir1 files first, then all dir2 files.
    for i in range(1, m + 1):
        filepath = os.path.join(dir1, f"{i}-{os.path.basename(dir1)}.txt")
        if not os.path.exists(filepath):
            print(f"WARNING: file does not exist {filepath}")
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
            print(f"WARNING: file does not exist {filepath}")
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
            print(f"Processing file {task['filepath']}")
            single_newtoboth(task["filepath"], task["output_dir"], task["num"], chr_series)
    else:
        max_workers = min(threads, len(tasks))
        log_dir = os.path.join(work_dir, "parallel_logs", "newtoboth")
        os.makedirs(log_dir, exist_ok=True)
        for task in tasks:
            safe_label = sanitize_filename(task["label"])
            task["log_file"] = os.path.join(log_dir, f"newtoboth_{safe_label}.log")

        print(f"Enabled parallel newtoboth file conversion: workers={max_workers}, tasks={len(tasks)}")
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
    """Sanitize special characters that are not allowed in file names."""
    return re.sub(r'[\\/*?:"<>|]', "", name)  # Replace special characters in the name to make it safe


def get_column_count(file_path):
    """Return the file column count divided by 3."""
    try:
        with open(file_path, 'r') as file:
            first_line = file.readline().strip()
            column_count = len(first_line.split())
            return column_count // 3
    except Exception as e:
        print(f"Error while reading file: {e}")
        return None


def parse_filename(filename):
    """Parse a file name and extract the replicate index and methylation type using a regular expression."""
    # Use capturing groups(the parenthesized parts here)to capture the file index and methylation type
    pattern = r'^(\d+)-bothMeUnme_diffChromo_NOREPEATED_methy_sites_(.+)\.txt$'
    match = re.match(pattern, filename)  # Match the file name against this regular expression
    if match:  # After matching, use match.group(i) to get the i-th captured group
        file_id = int(match.group(1))
        methylation_type = match.group(2)
        return file_id, methylation_type
    return None


def scan_sample_files_by_replicates(sample_dir, max_replicates):
    """Search sample_dir for files of the form 1~max_replicates-both...methylation_type.txt."""
    files_by_replicates = {}  # Create a dictionary to store file names by replicate index
    methylation_types = ['CpG', 'CHH', 'CHG']
    # If the directory does not exist, return an empty dictionary; this can be omitted after merging because main checks this at the beginning
    # if not os.path.exists(sample_dir):
    #     return files_by_replicates

    # Iterate over all entries in this folder, including files and subdirectories
    for filename in os.listdir(sample_dir):
        if filename.endswith('.txt'):  # If a txt file is found (it can only be a both file or a new-format file)
            parsed = parse_filename(filename)  # Parse the file to obtain a possible replicate index and methylation type
            if parsed:  # If the replicate index and methylation type were parsed from this file
                file_id, methylation_type = parsed  # A tuple is indicated by the comma; parentheses are only used to avoid ambiguity
                if file_id <= max_replicates and methylation_type in methylation_types:  # Both replicate index and methylation type are valid
                    if file_id not in files_by_replicates:
                        files_by_replicates[file_id] = {}  # Make files_by_replicates become file_id->{sub-dictionary}
                    files_by_replicates[file_id][methylation_type] = filename  # Make the sub-dictionary format methylation type->filename
    return files_by_replicates  # This returns a mapping chain replicate index->methylation type->file name , so the file name can be retrieved by replicate index and methylation type


def process_methylation_type_with_collection(file1_path, file2_path, methylation_type, output_dir,
                                             dir1_name, dir2_name, chr_num, meth_diff_threshold=0.0):
    """
    Run Fisher tests for one methylation type and one chromosome in a single comparison.

    Parameters:
        file_path: relative path of the current both-format file.
        methylation_type: current methylation context.
        output_dir: output directory, e.g., ./output_x_y/.
        dir_name: final path component of the two genotype input directories.
        chr_num: ordinal chromosome index, not the chromosome label.
    """

    print(f"      Processing methylation context {methylation_type}; chromosome {chr_num}...")

    # Calculate data columns for the current chromosome
    mOrder = 3 * (chr_num - 1)  # Here chr_num means the ordinal chromosome index, not the chromosome label
    usecols = [mOrder, mOrder + 1, mOrder + 2]

    try:
        # Read all data from the first file into a dictionary in chunks
        data1_dict = {}
        for chunk in pd.read_csv(file1_path, sep=r'\s+', header=None, usecols=usecols, chunksize=100000):
            chunk.columns = ['pos', 'methy', 'unmethy']
            # Use zip to retrieve one value from each relevant column and pack the three values into a (pos, methy, unmethy) tuple
            for pos, methy, unmethy in zip(chunk['pos'].values, chunk['methy'].values, chunk['unmethy'].values):
                # Put data into the corresponding dictionary in the following format
                data1_dict[pos] = (methy, unmethy)
        # At this point, all site data from the first file have been loaded into data1_dict

        # Read the second file in chunks and look for common sites because loading both files fully may consume too much memory
        # The first file must be fully loaded because fast random lookup is needed to check whether a site exists in both files and should be tested
        # Because chunksize is set, read_csv returns an iterator; each iteration returns at most 100000 rows
        reader2 = pd.read_csv(file2_path, sep=r'\s+', header=None, usecols=usecols, chunksize=100000)

    except Exception as e:
        print(f"        Failed to load data: {e}")
        return False

    # Create output folder:output_x_y/methylation type/
    methylation_dir = os.path.join(output_dir, methylation_type)
    os.makedirs(methylation_dir, exist_ok=True)

    # Create output file path
    stats_filename = os.path.join(methylation_dir,
                                  f"FET_results_{methylation_type}_{dir1_name}_{dir2_name}_Chr{chr_num}.txt")
    all_output = os.path.join(methylation_dir, f"all_simple_Chr{chr_num}.txt")
    sig_output = os.path.join(methylation_dir, f"sig_simple_Chr{chr_num}.txt")
    combine_output = os.path.join(methylation_dir, f"combineResult_Chr{chr_num}.txt")

    # Set parameters
    M0, M2, Th_pValue = 2, 4, 0.05

    try:
        # Create four lists to store different types of data for later output
        all_results, sig_results, fet_results, combine_results = [], [], [], []

        # Process each data chunk from the second file
        for chunk2 in reader2:
            chunk2.columns = ['pos', 'methy', 'unmethy']
            # For each chunk from the second file, iterate over each row and store the three values in pos, m2, and u2
            for pos, m2, u2 in zip(chunk2['pos'].values, chunk2['methy'].values, chunk2['unmethy'].values):
                # If the site also exists in the first file, perform Fisher test for this site
                if pos in data1_dict:
                    m1, u1 = data1_dict[pos]  # Get the two counts from the first file

                    if m1 >= M0 or m2 >= M0:  # Both methylated read counts must be >= 2 before proceeding
                        if (m1 + u1 >= M2) and (m2 + u2 >= M2):  # The total count in each file must be >= 4 before testing
                            cont_table = np.array([[m1, u1], [m2, u2]])  # Build the 2x2 contingency table
                            # Calculate the two methylation rates
                            ratio1 = m1 / (m1 + u1) if (m1 + u1) > 0 else 0
                            ratio2 = m2 / (m2 + u2) if (m2 + u2) > 0 else 0
                            meth_diff = abs(ratio1 - ratio2)
                            change = "1" if ratio1 >= ratio2 else "0"  # Determine direction from methylation ratios in the two files
                            # whether the mutant methylation rate increased
                            # Call the library Fisher test and obtain the p-value
                            _, pvalue = fisher_exact(cont_table, alternative='two-sided')
                            # Keep 7 significant digits for the p-value
                            pvalue = float(f"{pvalue:.7g}")
                            # Append required data to the four output lists; only significant records are appended to sig_results
                            all_results.append([pos, pvalue, change, meth_diff])
                            fet_results.append([pos, m1, u1, m2, u2, pvalue, meth_diff])
                            combine_results.append([pos, m1, u1, m2, u2, meth_diff])

                            if pvalue < Th_pValue:
                                sig_results.append([pos, pvalue, change, meth_diff])
        # At this point, all required Fisher tests for this file pair, methylation type, and chromosome are complete
        # Save results to disk
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

        print(f"        chromosome {chr_num} processing completed; processed {len(all_results)}  sites,including {len(sig_results)}  significant sites")
        return True

    except Exception as e:
        print(f"        Processing methylation context {methylation_type}; chromosome {chr_num} error occurred: {e}")
        return False

def process_methylation_type_with_collection_pvfilter(file1_path, file2_path, methylation_type, output_dir,
                                             dir1_name, dir2_name, chr_num, meth_diff_threshold=0.0):
    """
    Run Fisher tests for one methylation type and one chromosome in a single comparison.

    Parameters:
        file_path: relative path of the current both-format file.
        methylation_type: current methylation context.
        output_dir: output directory, e.g., ./output_x_y/.
        dir_name: final path component of the two genotype input directories.
        chr_num: ordinal chromosome index, not the chromosome label.

    Returns:
        DataFrame for pvalue > 0.05 records from this test.
        all_results rows have the form [pos, pvalue, change].
    """

    print(f"      Processing methylation context {methylation_type}; chromosome {chr_num}...")

    # Calculate data columns for the current chromosome
    mOrder = 3 * (chr_num - 1)  # Here chr_num means the ordinal chromosome index, not the chromosome label
    usecols = [mOrder, mOrder + 1, mOrder + 2]

    try:
        # Read all data from the first file into a dictionary in chunks
        data1_dict = {}
        for chunk in pd.read_csv(file1_path, sep=r'\s+', header=None, usecols=usecols, chunksize=100000):
            chunk.columns = ['pos', 'methy', 'unmethy']
            # Use zip to retrieve one value from each relevant column and pack the three values into a (pos, methy, unmethy) tuple
            for pos, methy, unmethy in zip(chunk['pos'].values, chunk['methy'].values, chunk['unmethy'].values):
                # Put data into the corresponding dictionary in the following format
                data1_dict[pos] = (methy, unmethy)
        # At this point, all site data from the first file have been loaded into data1_dict

        # Read the second file in chunks and look for common sites because loading both files fully may consume too much memory
        # The first file must be fully loaded because fast random lookup is needed to check whether a site exists in both files and should be tested
        # Because chunksize is set, read_csv returns an iterator; each iteration returns at most 100000 rows
        reader2 = pd.read_csv(file2_path, sep=r'\s+', header=None, usecols=usecols, chunksize=100000)

    except Exception as e:
        print(f"        Failed to load data: {e}")
        return False

    # Create output folder:output_x_y/methylation type/
    methylation_dir = os.path.join(output_dir, methylation_type)
    os.makedirs(methylation_dir, exist_ok=True)

    # Create output file path
    stats_filename = os.path.join(methylation_dir,
                                  f"FET_results_{methylation_type}_{dir1_name}_{dir2_name}_Chr{chr_num}.txt")
    all_output = os.path.join(methylation_dir, f"all_simple_Chr{chr_num}.txt")
    sig_output = os.path.join(methylation_dir, f"sig_simple_Chr{chr_num}.txt")
    combine_output = os.path.join(methylation_dir, f"combineResult_Chr{chr_num}.txt")

    # Set parameters
    M0, M2, Th_pValue = 2, 4, 0.05

    try:
        # Create four lists to store different types of data for later output
        all_results, sig_results, fet_results, combine_results = [], [], [], []

        # Process each data chunk from the second file
        for chunk2 in reader2:
            chunk2.columns = ['pos', 'methy', 'unmethy']
            # For each chunk from the second file, iterate over each row and store the three values in pos, m2, and u2
            for pos, m2, u2 in zip(chunk2['pos'].values, chunk2['methy'].values, chunk2['unmethy'].values):
                # If the site also exists in the first file, perform Fisher test for this site
                if pos in data1_dict:
                    m1, u1 = data1_dict[pos]  # Get the two counts from the first file

                    if m1 >= M0 or m2 >= M0:  # Both methylated read counts must be >= 2 before proceeding
                        if (m1 + u1 >= M2) and (m2 + u2 >= M2):  # The total count in each file must be >= 4 before testing
                            cont_table = np.array([[m1, u1], [m2, u2]])  # Build the 2x2 contingency table
                            # Calculate the two methylation rates
                            ratio1 = m1 / (m1 + u1) if (m1 + u1) > 0 else 0
                            ratio2 = m2 / (m2 + u2) if (m2 + u2) > 0 else 0
                            meth_diff = abs(ratio1 - ratio2)
                            change = "1" if ratio1 > ratio2 else "0"  # Determine direction from methylation ratios in the two files
                            # whether the mutant methylation rate increased
                            # Call the library Fisher test and obtain the p-value
                            _, pvalue = fisher_exact(cont_table, alternative='two-sided')
                            # Keep 7 significant digits for the p-value
                            pvalue = float(f"{pvalue:.7g}")
                            # Append required data to the four output lists; only significant records are appended to sig_results
                            all_results.append([pos, pvalue, change, meth_diff])
                            fet_results.append([pos, m1, u1, m2, u2, pvalue, meth_diff])
                            combine_results.append([pos, m1, u1, m2, u2, meth_diff])

                            if pvalue < Th_pValue:
                                sig_results.append([pos, pvalue, change, meth_diff])
        # At this point, all required Fisher tests for this file pair, methylation type, and chromosome are complete
        # Save results to disk
        if all_results:
            # Convert to DataFrame for filtering
            all_df = pd.DataFrame(all_results, columns=["Position", "Pvalue", "Methylation_Change", "MethDiff"])
            sig_df = pd.DataFrame(sig_results, columns=["Position", "Pvalue", "Methylation_Change", "MethDiff"])
            fet_df = pd.DataFrame(fet_results,
                                  columns=["Position", "Methy1", "Unmethy1", "Methy2", "Unmethy2", "Pvalue", "MethDiff"])
            combine_df = pd.DataFrame(combine_results, columns=["Position", "Methy1", "Unmethy1", "Methy2", "Unmethy2", "MethDiff"])

            all_df_ndmp = all_df[all_df['Pvalue'] > 0.05]  # Keep pvalue > 0.05 records for later addition to the q-value list

            # Keep only sites with p <= 0.05
            all_df = all_df[all_df['Pvalue'] <= 0.05]  # Filter pvalue <= 0.05 records for FDR correction
            fet_df = fet_df[fet_df['Pvalue'] <= 0.05]  # Filter pvalue <= 0.05 records for FDR correction
            #fet_df_ndmp = all_df[all_df['Pvalue'] > 0.05]
            # sig_df already contains p < 0.05 records, so no additional filtering is needed

            # Also filter combine
            positions_to_keep = set(all_df['Position'].values)
            combine_df = combine_df[combine_df['Position'].isin(positions_to_keep)]  # Filter pvalue <= 0.05 records for FDR correction
            #combine_df_ndmp = combine_df[~combine_df['Position'].isin(positions_to_keep)]
            # Save the filtered results
            all_df.to_csv(all_output, sep='\t', index=False)
            sig_df.to_csv(sig_output, sep='\t', index=False)
            fet_df.to_csv(stats_filename, sep='\t', index=False)
            combine_df.to_csv(combine_output, sep='\t', index=False)

        print(
            f"        chromosome {chr_num} processing completed; raw tests = {len(all_results)}  sites,after filtering (p <= 0.05): {len(all_df)}  sites,p > 0.05: {len(all_df_ndmp)} sites")
        return all_df_ndmp
    #                   all_results.append([pos, pvalue, change])
    #                   fet_results.append([pos, m1, u1, m2, u2, pvalue])
    #               combine_results.append([pos, m1, u1, m2, u2])

    except Exception as e:
        print(f"        Processing methylation context {methylation_type}; chromosome {chr_num} error occurred: {e}")
        return False

def merge_fet_results_and_fdr(output_dir, replicate_x, replicate_y, mtype3, all_dfs_ndmp_dict, n_chromosomes, is_twostep_context=False):
    """
    Merge all FET result files in one output directory and perform FDR correction.
    FET files have the form: pos, m1, u1, m2, u2, pvalue.
    The input output_dir is output_x_y/methylation_type.
    success_dfs_dict[methylation_type][chr_num] accesses the corresponding pvalue > 0.05 DataFrame:
        all_results([pos, pvalue, change]).
    n_chromosomes is the total chromosome count.
    """
    print(f"\n    Merging {output_dir}  FET results and performing FDR correction...")

    if not os.path.exists(output_dir):
        print(f"    ERROR: directory {output_dir} does not exist.")
        return False, get_dmp_threshold(mtype3)

    # Search all FET result files
    # Here ** means subdirectories at any depth; with recursive=True, glob searches under output_dir and all of its
    # subdirectories for matching files and returns matching relative paths from output_dir as a list
    file_pattern = os.path.join(output_dir, "**", "FET_results_*_Chr*.txt") # Here the chromosome number actually means the index of the three-column block in the both file
    fet_files = glob.glob(file_pattern, recursive=True)

    if not fet_files:
        print(f"    WARNING: no {output_dir}  FET result files found in ")
        return False, get_dmp_threshold(mtype3)

    print(f"    Found {len(fet_files)}  FET result files")

    # Create a list to collect all p-values and related information
    all_data = []

    for file_path in sorted(fet_files):
        # Extract methylation type and chromosome information
        # The first capturing group captures methylation type, .* matches replicatex_replicatey, and the second captures chromosome number
        # # Here the chromosome number actually means the index of the three-column block in the both file
        methy_match = re.search(r'/FET_results_([^_]+)_.*_Chr(\d+)\.txt$', file_path.replace('\\', '/'))
        if not methy_match:
            continue

        # Get methylation type and chromosome index(chromosome numberboth)
        methylation_type = mtype3
        chr_num = int(methy_match.group(2))

        try:
            df = pd.read_csv(file_path, sep='\t', header=0)
            if 'Position' in df.columns and 'Pvalue' in df.columns:
                # Adjust column order:Chromosome, Methylation_Type, Position, Pvalue
                df_subset = df[['Position', 'Pvalue', 'MethDiff']].copy()
                df_subset['Chromosome'] = chr_num
                df_subset['Methylation_Type'] = methylation_type

                # Reorder columns
                df_subset = df_subset[['Chromosome', 'Methylation_Type', 'Position', 'Pvalue', 'MethDiff']]
                all_data.append(df_subset) # Append all information for the current FET file to all_data
        except Exception as e:
            print(f"    WARNING: read {file_path} failed: {e}")
            continue

    if not all_data:
        print(f"    ERROR: no data were successfully read")
        return False, get_dmp_threshold(mtype3)

    # Merge all data(Since each file-specific DataFrame was appended to all_data, each list element is a DataFrame and must be concatenated)
    combined_df = pd.concat(all_data, ignore_index=True)
    # Sort data by chromosome index(Here the chromosome number actually means the index of the three-column block in the both file)
    combined_df = combined_df.sort_values(['Methylation_Type', 'Chromosome', 'Position'])
    # The DataFrame format is as follows:'Chromosome', 'Methylation_Type', 'Position', 'Pvalue'
    print(f"    After merging, total sites = {len(combined_df)}  sites")

    # Calculate FDR-corrected q-values and append them to combined_df
    pvalues = combined_df['Pvalue'].values
    qvalues = calculate_qvalues(pvalues, 1.0)
    combined_df['Qvalue'] = qvalues

    fixed_dmp_threshold = get_dmp_threshold(mtype3)
    dmp_threshold_to_use = fixed_dmp_threshold
    auto_q_info = None

    # Use automatic q-value thresholds only in two-step contexts:
    # i.e., this context has first been prefiltered by p <= AUTO_QVALUE_P_CUTOFF and then FDR-corrected within the prefiltered subset.
    # Disable this for non-two-step contexts to avoid applying this rule to the full-FDR q-value curve.
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
            mode_msg = "used for DMP calling"
        else:
            dmp_threshold_to_use = fixed_dmp_threshold
            mode_msg = "report-only; DMP calling unchanged"

        print(
            f"    Auto q-value threshold ({mtype3}, two-step-FDR context): "
            f"estimated={estimated_q:.6g}, used={dmp_threshold_to_use:.6g} "
            f"({mode_msg}, status={auto_q_info.get('auto_q_status')}, "
            f"p_at_max={auto_q_info.get('auto_q_pvalue_at_max')}, "
            f"diff_at_max={auto_q_info.get('auto_q_diff_at_max')}, "
            f"n={auto_q_info.get('auto_q_n_candidates')})"
        )
    elif AUTO_QVALUE_TWOSTEP or AUTO_QVALUE_REPORT_ONLY:
        print(f"    {mtype3} is not a two-step-FDR context; skipping auto q-value threshold estimation and using fixed threshold {fixed_dmp_threshold}")

    combined_df['Qvalue_Threshold_Used'] = dmp_threshold_to_use
    combined_df['Qvalue_Threshold_Mode'] = (
        "auto_twostep" if (AUTO_QVALUE_TWOSTEP and not AUTO_QVALUE_REPORT_ONLY and bool(is_twostep_context))
        else "fixed"
    )

    # Final column order:Chromosome, Methylation_Type, Position, Pvalue, Qvalue
    combined_df = combined_df[['Chromosome', 'Methylation_Type', 'Position', 'Pvalue', 'MethDiff', 'Qvalue',
                               'Qvalue_Threshold_Used', 'Qvalue_Threshold_Mode']]

    print(f"Starting {mtype3}  NDMP data,")

    # First check whether all_dfs_ndmp_dict has data for this methylation type because some organisms/contexts may not require the filter version
    if mtype3 in all_dfs_ndmp_dict and all_dfs_ndmp_dict[mtype3]:
        print(f"Number of NDMP dataframes for this methylation context: {len(all_dfs_ndmp_dict[mtype3])}")
        dfs_ndmp = []
        for chr_num11 in range(1, n_chromosomes + 1):
            # Process only chromosomes present for this methylation type
            if chr_num11 in all_dfs_ndmp_dict[mtype3]:
                df_ndmp = all_dfs_ndmp_dict[mtype3][chr_num11]

                # Ensure it is a valid DataFrame
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
                    print(f"  for {mtype3}--{chr_num11} added required columns; {len(df_ndmp)}  sites")

        if dfs_ndmp:
            total_ndmp = pd.concat(dfs_ndmp, ignore_index=True)
            combined_df = pd.concat([combined_df, total_ndmp], ignore_index=True)
            print(f"  Merged {len(total_ndmp)}  NDMP sites")
        else:
            print(f"  {mtype3} context has no NDMP data")
    else:
        print(f"  {mtype3} context is not present in the NDMP dictionary")
    # At this point, for high-count methylation types, information for previously discarded pvalue > 0.05 sites has been added back,
    # They must have been stored in all_dfs_ndmp_dict[methylation_type][chromosome] earlier for this to work

    # Summarize significance counts for one of the m*n*3 tests: numbers significant by p-value and q-value
    n_pval_sig = np.sum(pvalues <= 0.05)
    dmp_threshold = dmp_threshold_to_use
    n_qval_sig = np.sum(qvalues <= dmp_threshold)
    # Calculate significant proportions
    print(f"    P-value significant sites: {n_pval_sig} ({n_pval_sig / len(pvalues) * 100:.1f}%)")
    print(f"    Q-value significant sites: {n_qval_sig} ({n_qval_sig / len(qvalues) * 100:.1f}%)")

    # Save merged p-value list for external FDR tools
    pvalue_file = os.path.join(output_dir, f"united_pvalues_wt_replicate{replicate_y}_vs_mut_replicate{replicate_x}.csv")
    with open(pvalue_file, 'w') as f:
        for pvalue in pvalues:
            f.write(f"{pvalue}\n")

    # Save full FDR-corrected results(output_dir is output_x_y/methylation_type)
    # Includes complete p-values and q-values(Format::Chromosome, Methylation_Type, Position, Pvalue, Qvalue)
    fdr_file = os.path.join(output_dir, f"FDR_corrected_results_wt_replicate{replicate_y}_vs_mut_replicate{replicate_x}.txt")
    combined_df.to_csv(fdr_file, sep='\t', index=False)

    # Save significant sites(Use fixed thresholds or two-step automatic thresholds according to methylation type)
    sig_df = combined_df[combined_df['Qvalue'] <= dmp_threshold_to_use]
    if not sig_df.empty: # If significant sites exist, write the significant subset(output_dir is output_x_y/methylation_type)
        sig_file = os.path.join(output_dir, f"FDR_significant_results_wt_replicate{replicate_y}_vs_mut_replicate{replicate_x}.txt")
        sig_df.to_csv(sig_file, sep='\t', index=False)
        print(f"    Significant-site results saved to: {sig_file}")

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

    print(f"    P-value list saved to: {pvalue_file}")
    print(f"    FDR results saved to: {fdr_file}")
    return True, dmp_threshold_to_use


# In this code, replace calculate_qvalues with:
def     calculate_qvalues(pvalues, pi=1.0):
    """Calculate q-values using the Storey method."""

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

    # Key point of the Storey method: include the pi factor
    q_values_sorted = np.zeros_like(sorted_pvalues)

    # If pi needs to be estimated automatically
    if pi is None:
        # Simplified pi estimation
        lrange = np.linspace(0.05, 0.95, max(int(m / 100.0), 10))
        pil = np.mean(sorted_pvalues[:, np.newaxis] > lrange, axis=0)
        pilr = pil / (1.0 - lrange)
        pi = 1.0
        if pilr[-1] < 1.0:
            pi = pilr[-1]

    # Storey method calculation
    q_values_sorted = pi * m * sorted_pvalues / np.arange(1, m + 1)
    q_values_sorted[-1] = min(q_values_sorted[-1], 1.0)

    # Monotonicity adjustment
    for i in range(m - 2, -1, -1):
        q_values_sorted[i] = min(q_values_sorted[i], q_values_sorted[i + 1])

    # Restore original order
    q_values = np.zeros_like(pvalues_clean)
    q_values[sorted_indices] = q_values_sorted
    q_values[np.isnan(pvalues)] = np.nan

    return q_values

def perform_sliding_window_on_dmp_files(output_dir, replicate_x, replicate_y,):
    """Run sliding-window analysis on one result (N)DMP file from one of the m*n*3 tests."""

    print(f"\n    Starting sliding-window analysis for DMP files...")

    # Process only DMP files
    dmp_patterns = [
        f"DMP_wt_replicate{replicate_y}_mut_replicate{replicate_x}_Chr*.txt",
        f"N-DMP_wt_replicate{replicate_y}_mut_replicate{replicate_x}_Chr*.txt"
    ]

    for pattern in dmp_patterns:
        dmp_files = glob.glob(os.path.join(output_dir, pattern)) # Find all file paths matching the current regex and put them into dmp_files

        for dmp_file in dmp_files: # Iterate over the list and get one DMP or N-DMP file each time
            try:

                df = pd.read_csv(dmp_file, sep=r'\s+', header=None, skiprows=1,
                                 names=['position', 'pvalue', 'change']) # Position, p-value, and change

                # Check that the data are non-empty
                if df.empty or len(df) == 0:
                    continue

                # Convert data types and clean data
                df = df.dropna()
                df['position'] = df['position'].astype(int)
                df['change'] = df['change'].astype(int)

                print(f"      Processing file: {os.path.basename(dmp_file)} ({len(df)}  sites)")

                # Set output prefix, i.e., the file name with .txt removed (N)DMP_replicate_wt{replicate_y}_mut_replicate{replicate_x}_Chr* part
                base_name = os.path.basename(dmp_file).replace('.txt', '')

                # Pass the DataFrame to sliding_window_analysis
                # The former format is window_start, window_end, count_change_1, count_change_0, total_count
                # The latter has no count_change but has standardized_count, the ratio of significant-site count in each interval to the maximum interval count
                sliding_results, std_results = sliding_window_analysis(
                    df, # df is the DMP file for one chromosome from one of the m*n*3 tests
                    window_size=1000000,
                    step_ratio=0.05,
                    save_files=True,
                    output_identifier=base_name,
                    outputdir1=output_dir
                )

                print(f"      Completed {base_name}: {len(sliding_results)} windows")

            except Exception as e:
                print(f"      Processing {dmp_file} failed: {e}")
                continue

    print(f"    Sliding-window analysis completed")

# This version handles cases where p < 0.05 filtering is applied before FDR correction, e.g., plant CHH and CHG
def perform_sliding_window_on_dmp_files_after_filter(output_dir, replicate_x, replicate_y,all_dfs_ndmp_dict=None, methylation_type=None):
    """Run sliding-window analysis on one result (N)DMP file from one of the m*n*3 tests."""

    print(f"\n    Starting sliding-window analysis for DMP files...")

    # Process only DMP files
    dmp_patterns = [
        f"DMP_wt_replicate{replicate_y}_mut_replicate{replicate_x}_Chr*.txt",
        f"N-DMP_wt_replicate{replicate_y}_mut_replicate{replicate_x}_Chr*.txt"
    ]

    for pattern in dmp_patterns:
        dmp_files = glob.glob(os.path.join(output_dir, pattern)) # Find all file paths matching the current regex and put them into dmp_files
        is_ndmp = pattern.startswith("N-DMP") # Check whether the current file is an NDMP file

        for dmp_file in dmp_files: # Iterate over the list and get one DMP or N-DMP file each time
            try:

                df = pd.read_csv(dmp_file, sep=' ', header=None, skiprows=1,
                                 names=['position', 'pvalue', 'change']) # Position, p-value, and change

                # If the current file is an NDMP file, add the previously pvalue > 0.05 records here because all_simple cannot be modified since it is used for FDR merging
                                # so that NDMP information is not missing when sliding-window files are generated later
                if is_ndmp and all_dfs_ndmp_dict and methylation_type:
                    # Extract chromosome number from file name
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
                                # Merge data
                                df = pd.concat([df, ndmp_df[['position', 'pvalue', 'change']]],
                                               ignore_index=True)
                                df = df.drop_duplicates(subset=['position'])
                                df = df.sort_values('position').reset_index(drop=True)
                                print(f"      for {os.path.basename(dmp_file)} merged NDMP sites ignored by the previous FDR test: {len(ndmp_df)}  NDMP sites")
                # Check that the data are non-empty
                if df.empty or len(df) == 0:
                    continue

                # Convert data types and clean data
                df = df.dropna()
                df['position'] = df['position'].astype(int)
                df['change'] = df['change'].astype(int)

                print(f"      Processing file: {os.path.basename(dmp_file)} ({len(df)}  sites)")

                # Set output prefix, i.e., the file name with .txt removed (N)DMP_wt_replicate{replicate_y}_mut_replicate{replicate_x}_Chr* part
                base_name = os.path.basename(dmp_file).replace('.txt', '')

                # Pass the DataFrame to sliding_window_analysis
                # The former format is window_start, window_end, count_change_1, count_change_0, total_count
                # The latter has no count_change but has standardized_count, the ratio of significant-site count in each interval to the maximum interval count
                sliding_results, std_results = sliding_window_analysis(
                    df, # df is the DMP file for one chromosome from one of the m*n*3 tests
                    window_size=1000000,
                    step_ratio=0.05,
                    save_files=True,
                    output_identifier=base_name,
                    outputdir1=output_dir
                )

                print(f"      Completed {base_name}: {len(sliding_results)} windows")

            except Exception as e:
                print(f"      Processing {dmp_file} failed: {e}")
                continue

    print(f"    Sliding-window analysis completed")


def generate_dmp_files(dir1,dir2,output_dir, replicate_x, replicate_y, fdr_threshold=0.05, mtype1="CpG",
                      all_dfs_ndmp_dict=None,unfilter_mtypes=["CpG"],n_chromosomes = 5, meth_diff_threshold=0.0,
                      skip_dmr=False, skip_window=False):
    """
    Parameters:
        output_dir is output_x_y/methylation_type/.
        replicate_x and replicate_y are replicate IDs, followed by the q-value threshold and current methylation type.
    """

    print(f"\n    Generating {output_dir}  DMP files...")

    def safe_float_convert(value):
        """Convert data in multiple formats to a floating-point number."""
        try:
            # If already numeric, convert to float and return
            if isinstance(value, (int, float)):
                return float(value)
            # If string, strip whitespace and convert to float
            if isinstance(value, str):
                value = value.strip()
                if value:
                    return float(value)
            return None
        except (ValueError, TypeError):
            return None

    # 1. Read the FDR-corrected result file(output_dir is output_x_y/methylation_type/)
    # The file format is Chromosome, Methylation_Type, Position, Pvalue, Qvalue
    fdr_file = os.path.join(output_dir, f"FDR_corrected_results_wt_replicate{replicate_y}_vs_mut_replicate{replicate_x}.txt")
    if not os.path.exists(fdr_file):
        print(f"    ERROR: FDR result file does not exist {fdr_file}")
        return False

    try:
        fdr_df = pd.read_csv(fdr_file, sep='\t') # The file format is Chromosome, Methylation_Type, Position, Pvalue, Qvalue
        print(f"    Read FDR result file; {len(fdr_df)}  sites")
    except Exception as e:
        print(f"    ERROR: failed to read FDR file {e}")
        return False

    # 2. Read all_simple files for each chromosome
    methylation_change_data = {}
    total_read_lines = 0

    # (output_dir is output_x_y/methylation_type/)
    mtype_dir = output_dir

    # Get the file names of all all_simple files generated for one of the m*n*3 tests
    all_simple_files = [f for f in os.listdir(mtype_dir) if f.startswith('all_simple_Chr') and f.endswith('.txt')]

    # The all_simple file format is:"Position", "Pvalue", "Methylation_Change"
    # Iterate over all all_simple files for the given test
    for file in all_simple_files:
        chr_match = re.search(r'Chr(\d+)\.txt$', file)
        if not chr_match:
            continue

        chr_num = int(chr_match.group(1))# Here chr_num also means the index of the three-column block in the both file
        file_path = os.path.join(mtype_dir, file)

        try:
            with open(file_path, 'r') as f:  # The read all_simple file format is:"Position", "Pvalue", "Methylation_Change"
                lines = f.readlines() # readlines returns all file lines as a list, with each element being one line

            file_valid_lines = 0
            for line_num, line in enumerate(lines, 1): # Iterate over each line, with line_num enumerated from 1
                line = line.strip()
                if not line or line.startswith("Position"): # Skip empty lines or the first line
                    continue

                parts = line.split('\t')    # Split by the delimiter used during export and put the three fields into parts
                if len(parts) >= 3:
                    # Because the input is string, convert pvalue to float and the other two values to integers
                    position = int(safe_float_convert(parts[0]))
                    pvalue = safe_float_convert(parts[1])
                    change = int(safe_float_convert(parts[2]))
                    meth_diff = safe_float_convert(parts[3]) if len(parts) >= 4 else 0.0

                    # Check whether all values are valid
                    if (position is not None and
                            pvalue is not None and
                            change is not None):

                        # Check whether change is in {0, 1}
                        if change in [0, 1]:
                            # Store this mapping: (chr, mtype, position) -> change in the dictionary
                            methylation_change_data[(chr_num, mtype1, position)] = (change, meth_diff)
                            file_valid_lines += 1
                            total_read_lines += 1

            print(f"    {mtype_dir}/Chr{chr_num}: read {file_valid_lines}  valid sites")

        except Exception as e:
            print(f"    WARNING: read {file_path} failed: {e}")
            continue

    print(f"    Total methylation-change direction records read: {total_read_lines}  sites")

    # 3. Merge data using the previously recorded dictionary mapping(chr, mtype, position) -> change,add the change attribute to a copy of fdr_df
    combined_data = []
    missing_change = 0
    match_debug = defaultdict(int) # Create a defaultdict; when a missing key is accessed, it automatically calls int()
                                                # Create a key-value pair and set its value to int(), i.e., 0

    for _, row in fdr_df.iterrows(): # fdr_df fileFormat: Chromosome, Methylation_Type, Position, Pvalue, Qvalue
        chr_num = row['Chromosome']
        mtype = row['Methylation_Type']
        # Convert uniformly to float for matching
        position = safe_float_convert(row['Position'])
        qvalue = safe_float_convert(row['Qvalue'])

        # Skip the current row if position or q-value is missing
        if position is None or qvalue is None:
            missing_change += 1
            continue

        # Look up the corresponding methylation-change direction
        change_key = (chr_num, mtype, position)
        if change_key in methylation_change_data:
            change, meth_diff = methylation_change_data[change_key]  # change_data stores this mapping: (chr, mtype, position) -> change
            combined_data.append({
                'chromosome': chr_num,
                'methylation_type': mtype,
                'position': int(position),  # Convert to integer for output
                'qvalue': qvalue,
                'change': change,
                'meth_diff': meth_diff # Record absolute methylation difference; at this point combined_data is still a list rather than a DataFrame
            })
            match_debug[chr_num] += 1 # Here chr_num is also the three-column block index in the both file; this convention starts when the both file is read for Fisher tests
        else:
            missing_change += 1

    print(f"    After merging, {len(combined_data)}  complete sites")
    print(f"    Per-chromosome matching summary: {dict(match_debug)}") # This is the number of records for each chromosome
    if missing_change > 0:
        print(f"    WARNING: {missing_change}  sitesmissing methylation-change direction information")

    # 4. Generate DMP files by chromosome group
    chr_groups = defaultdict(list) # This defaultdict calls list() during initialization and uses an empty list as the value
    for item in combined_data: # Each item is a dictionary with values for chr_num, mtype, position, qvalue, and change
        chr_groups[item['chromosome']].append(item)
            # chr_groups stores each record for each chromosome as a dictionary element under chr_groups[chr_num]

    if not chr_groups:
        print(f"    ERROR:no data available for generating DMP files")
        return False

    print(f"    DMP files will be generated for chromosomes: {sorted(chr_groups.keys())}")

    dir1_name = f"replicate{replicate_x}"
    dir2_name = f"replicate{replicate_y}"

    total_dmp = total_ndmp = total_hyper = total_hypo = 0

    # Generate a file for each chromosome
    for chr_num in sorted(chr_groups.keys()): # The keys are chromosome indices; key - 1 corresponds to the zero-based index in the initial chromosome mapping
        chr_data = chr_groups[chr_num]  # The chr_groups dictionary stores each record for each chromosome(dictionary form)as an element
                                    # Therefore chr_data is still a dictionary here
        chr_data.sort(key=lambda x: x['position']) # The lambda anonymous function is applied to each list element and its return value is used as the sort key
                                # This anonymous function is equivalent to applying the function to each dictionary element in chr_data
                                                #   def get_position(x):
                                            # return x['position'] returns the position, so records are sorted by position

        # Generate file name(output is the methylation-type directory here)
        dmp_file = os.path.join(output_dir, f"DMP_wt_{dir2_name}_mut_{dir1_name}_Chr{chr_num}.txt")
        ndmp_file = os.path.join(output_dir, f"N-DMP_wt_{dir2_name}_mut_{dir1_name}_Chr{chr_num}.txt")
        hyper_file = os.path.join(output_dir, f"hyper_DMP_wt_{dir2_name}_mut_{dir1_name}_Chr{chr_num}.txt")
        hypo_file = os.path.join(output_dir, f"hypo_DMP_wt_{dir2_name}_mut_{dir1_name}_Chr{chr_num}.txt")

        # Classify data
        dmp_data = []
        ndmp_data = []
        hyper_data = []
        hypo_data = []

        # Iterate over the chr_data dictionary
        for item in chr_data: # chr_data is a dictionary, with five key-value fields:chr_num,mtype,position,qvalue,change
            position = item['position']
            qvalue = item['qvalue']
            change = item['change']

            if qvalue <= fdr_threshold and item.get('meth_diff', 0.0) >= meth_diff_threshold: # Determine significance while also requiring the methylation-difference threshold
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
                    # Get the set of existing sites to avoid duplicate additions
                    existing_positions = set(item[0] for item in ndmp_data)
                    added_count = 0
                    for _, row in ndmp_df.iterrows():
                        pos = int(row['Position'])
                        pval = float(row['Pvalue'])
                        chg = int(row['Methylation_Change'])
                        # Add only non-duplicate sites
                        if pos not in existing_positions:
                            ndmp_data.append((pos, pval, chg))
                            added_count += 1
                    if added_count > 0:
                        print(f"      for Chr{chr_num} added {added_count}  pvalue > 0.05 NDMP sites")


        # Write file
        def write_dmp_file(filename, data):
            with open(filename, 'w') as f:
                f.write("first line\n")
                for pos, qval, chg in data:
                    f.write(f"{pos} {qval} {chg}\n")

        write_dmp_file(dmp_file, dmp_data)
        write_dmp_file(ndmp_file, ndmp_data)
        write_dmp_file(hyper_file, hyper_data)
        write_dmp_file(hypo_file, hypo_data)

        # print(f" Skip Chr{chr_num} single-pair DMR candidate detection; branch DMR output has been removed and final DMR will be recalculated from common DMPs")

        print(f"    Chr{chr_num}: DMP={len(dmp_data)}, N-DMP={len(ndmp_data)}, Hyper={len(hyper_data)}, Hypo={len(hypo_data)}")

        # generate_dmp_files(output_dir, replicate_x, replicate_y, fdr_threshold=0.05, mtype1="CpG",
        #                    all_dfs_ndmp_dict=None, unfilter_mtypes=["CpG"]):
        # Count
        total_dmp += len(dmp_data)
        total_ndmp += len(ndmp_data)
        total_hyper += len(hyper_data)
        total_hypo += len(hypo_data)

    print(f"    DMP file generation completed.")
    print(f"    Total: DMP={total_dmp}, N-DMP={total_ndmp}, Hyper={total_hyper}, Hypo={total_hypo}")

    bothfile1 = os.path.join(dir1,f"{replicate_x}-bothMeUnme_diffChromo_NOREPEATED_methy_sites_{mtype1}.txt")

    bothfile2 = os.path.join(dir2,f"{replicate_y}-bothMeUnme_diffChromo_NOREPEATED_methy_sites_{mtype1}.txt")

    # print(f" Skip {mtype1} single-pair DMR read summarization/Fisher/FDR; branch DMR output has been removed and final DMR will be recalculated in and_output/dmr_analysis_*")

    if skip_window:
        print(f"    skipping {mtype1} single-run  DMP/N-DMP sliding-window analysis (--skip-window)")
    else:
        # Choose the sliding-window function according to methylation type
        if mtype1 not in unfilter_mtypes:
            print(f"Using the new sliding-window analysis implementation")
            # Methylation type also needs to be distinguished here
            perform_sliding_window_on_dmp_files_after_filter(
                output_dir, replicate_x, replicate_y,
                all_dfs_ndmp_dict=all_dfs_ndmp_dict,
                methylation_type=mtype1
            )
            # perform_sliding_window_on_dmp_files(
            #     output_dir, replicate_x, replicate_y
            # )
        else:
            print(f"    Using the standard sliding-window analysis implementation")
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
        "\nRebuilding Generating pairwise DMP output:"
        f"q passes and abs(MethDiff)>={meth_diff_threshold:.6g}"
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
                    print(f"  [WARN] skipping; FDR file does not exist: {fdr_file}")
                    continue

                fdr = pd.read_csv(fdr_file, sep="\t")
                required = {
                    "Chromosome", "Methylation_Type", "Position",
                    "Qvalue", "MethDiff",
                }
                missing = required - set(fdr.columns)
                if missing:
                    raise ValueError(
                        f"{fdr_file} is missing required columns:{sorted(missing)}"
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
                        f"{fdr_file}: {missing_raw} FDR sites in the original WT/MUT inputs "
                        f"could not be matched,examples={examples}"
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
                    f"  [OK] wt{wt_idx}-mut{mut_idx} {ctx}: "
                    f"DMP={total_dmp}, N-DMP={total_ndmp}"
                )

    summary = pd.DataFrame(summary_records)
    summary_dir = work_path / "and_output" / "auto_methdiff_thresholds"
    summary_dir.mkdir(parents=True, exist_ok=True)
    summary_file = summary_dir / "pair_dmp_regeneration_summary.tsv"
    summary.to_csv(summary_file, sep="\t", index=False)
    print(f"pairwise DMPrebuild summary tablesaved to: {summary_file}")
    return summary

# Integrate this function into process_replicate_pair
def process_replicate_pair(replicate_x, replicate_y, files1, files2, dir1, dir2, dir1_name, dir2_name, unfilter_mtypes, work_dir=".", meth_diff_threshold=0.0, skip_dmr=False, skip_window=False):
    """
    Process all methylation types for one replicate pair, corresponding to 1*3 tests.
    replicate_x and replicate_y are the indices of both-format files.
    files is a mapping chain: replicate index -> methylation type -> corresponding file name.
    dir1 and dir2 are the two genotype data directories.
    dir1_name and dir2_name are the final path components of the given dir1 and dir2 paths.
    """

    print(f"\n  Processing pair (wt{replicate_y}, mut{replicate_x})...")

    # Create the output directory for this comparison
    output_dir = os.path.join(work_dir, f"output_wt{replicate_y}_mut{replicate_x}")
    os.makedirs(output_dir, exist_ok=True)

    methylation_types = ['CpG', 'CHH', 'CHG']

    all_dfs_ndmp_dict = {}

    pair_success_count = 0
    pair_total_tests = 0

    # Loop over each methylation type
    for methylation_type in methylation_types:
        if methylation_type not in all_dfs_ndmp_dict:
            all_dfs_ndmp_dict[methylation_type] = {}

        success_count = total_tests = 0
        # Here "in files1[replicate_x]" checks whether methylation_type exists as a key in files1[replicate_x]; if it exists,
        # the corresponding i-both...methylation_type.txt file exists
        if (methylation_type not in files1[replicate_x] or
                methylation_type not in files2[replicate_y]):
            print(f"    skipping methylation context {methylation_type}: file does not exist")
            continue
        # Otherwise, get the relative paths of the files to be processed
        file1_path = os.path.join(dir1, files1[replicate_x][methylation_type])
        file2_path = os.path.join(dir2, files2[replicate_y][methylation_type])

        # Get the chromosome counts from the two files
        n_chromosomes_1 = get_column_count(file1_path)
        n_chromosomes_2 = get_column_count(file2_path)

        if n_chromosomes_1 is None or n_chromosomes_2 is None:
            print(f"    unable to get {methylation_type}  chromosome count")
            continue

        if n_chromosomes_1 != n_chromosomes_2:
            print(f"    {methylation_type} chromosome counts are inconsistent:{n_chromosomes_1} vs {n_chromosomes_2}")
            continue

        # At this point, both both-files contain chromosome data and have the same total column count
        # This check is not strictly necessary because newtoboth already guarantees matching column counts
        n_chromosomes = n_chromosomes_1  # Get the total chromosome count
        print(f"    Processing methylation context {methylation_type},with {n_chromosomes}  chromosomes")

        # Process each chromosome for the current x_y pair and methylation type; chr_num is the ordinal chromosome index, not the chromosome label
        for chr_num in range(1, n_chromosomes + 1):
            if methylation_type in unfilter_mtypes:
                success = process_methylation_type_with_collection(
                    file1_path, file2_path, methylation_type, output_dir,
                    dir1_name, dir2_name, chr_num, meth_diff_threshold=meth_diff_threshold
                )  # Process one chromosome for one test
            else:
                # Methylation type must be distinguished here
                all_df_ndmp = process_methylation_type_with_collection_pvfilter(
                    file1_path, file2_path, methylation_type, output_dir,
                    dir1_name, dir2_name, chr_num, meth_diff_threshold=meth_diff_threshold
                )   # Return a DataFrame of pvalue > 0.05 records for this test
                    #  all_results([pos, pvalue, change])
                if chr_num not in all_dfs_ndmp_dict[methylation_type]:
                    all_dfs_ndmp_dict[methylation_type][chr_num] = all_df_ndmp
                    if (isinstance(all_df_ndmp, pd.DataFrame)):
                        print(f"successfully added {methylation_type}-{chr_num} specific DataFrame to the dictionary; this DataFrame has {len(all_df_ndmp)} rows")
                success = isinstance(all_df_ndmp, pd.DataFrame)
            total_tests += 1
            if success:  # The previous step returns True if processing succeeded, otherwise False
                success_count += 1

        pair_success_count += success_count
        pair_total_tests += total_tests

        # Merge FET results and perform FDR correction
        # The FET result format is:pos, m1, u1, m2, u2, pvalue
        if success_count > 0:
            output_dir1 = os.path.join(output_dir, methylation_type)
            # Get all FET files under one output_x_y/methylation_type directory, concatenate them, calculate FDR q-values, and export to disk; final format:
            #                                        Chromosome, Methylation_Type, Position, Pvalue, Qvalue
            merge_ok, dmp_threshold = merge_fet_results_and_fdr(
                output_dir1, replicate_x, replicate_y, methylation_type,
                all_dfs_ndmp_dict, n_chromosomes,
                is_twostep_context=(methylation_type not in unfilter_mtypes)
            )
            if not merge_ok:
                print(f"    {methylation_type} ofFDRmerge failed; skippingDMPfile generation")
                continue
            # Generate DMP files: fixed thresholds or two-step automatic thresholds are returned by merge_fet_results_and_fdr
            generate_dmp_files(dir1,dir2,output_dir1, replicate_x, replicate_y, fdr_threshold=dmp_threshold, mtype1=methylation_type,all_dfs_ndmp_dict=all_dfs_ndmp_dict
                               ,unfilter_mtypes=unfilter_mtypes,n_chromosomes=n_chromosomes, meth_diff_threshold=meth_diff_threshold,
                               skip_dmr=skip_dmr, skip_window=skip_window)

    print(f"  pair (wt{replicate_y}, mut{replicate_x}) processing completed: successfully {pair_success_count}/{pair_total_tests} tests")
    return pair_success_count, pair_total_tests    # Here one chromosome is treated as one test


def process_all_combinations(dir1, dir2, m, n, unfilter_mtypes, work_dir=".", meth_diff_threshold=0.0,
                             skip_dmr=False, skip_window=False, threads=1):
    """Process all combinations and perform m*n*3 tests; replicate-pair-level parallelization is supported."""

    print(f"Scanning file directories...")
    files1 = scan_sample_files_by_replicates(dir1, m)
    files2 = scan_sample_files_by_replicates(dir2, n)

    print(f"directory1 ({dir1}) Found {len(files1)} replicate files")
    print(f"directory2 ({dir2}) Found {len(files2)} replicate files")

    missing_replicates1 = [i for i in range(1, m + 1) if i not in files1]
    missing_replicates2 = [i for i in range(1, n + 1) if i not in files2]
    if missing_replicates1:
        print(f"WARNING: directory1is missing these replicates: {missing_replicates1}")
    if missing_replicates2:
        print(f"WARNING: directory2is missing these replicates: {missing_replicates2}")

    available_replicates1 = [i for i in range(1, m + 1) if i in files1]
    available_replicates2 = [i for i in range(1, n + 1) if i in files2]

    total_combinations = len(available_replicates1) * len(available_replicates2)
    print(f"\nStarting {total_combinations}  combinations...")

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
            print(f"\nprogress: {i}/{total_combinations}")
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

        print(f"Enabled parallel replicate-pair Processing: workers={max_workers}, tasks={len(tasks)}")
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
                        f"[DONE] {i}/{len(tasks)} {label}: successfully "
                        f"{res['success_count']}/{res['test_count']} tests, "
                        f"{res['elapsed']:.2f}s, log={res['log_file']}"
                    )
                except Exception as e:
                    print(f"[FAILED] {i}/{len(tasks)} {label}: {e}")
                    raise

    end_time = time.time()
    print(f"\nAll processing completed!")
    print(f"Total: {total_success}/{total_tests} successful tests")
    print(f"Elapsed time: {end_time - start_time:.2f} seconds")
    print(f"Single-comparison results are saved in ./output_x_y/<methylation_context>/ directories")

    # Backward compatibility: if no countable tests exist but output did not raise an exception, do not mark it as failed.
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
        "  Rule:provisional final DMP in,if boundary_abs_methdiff <= cutoff "
        "and support_count < low_required, then remove the site from the final DMP set."
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
    print(f"  lowdiffdiagnostic tablesaved to: {diag_file}")
    print(f"  provisional DMPbackupsaved to: {provisional_file}")
    print(f"  lowdiffsummary tablesaved to: {summary_file}")

    return filtered_df


def bayes_deciding(sig_count, nonsig_count, auto_vote_threshold=None):
    """
    Make the final DMP/DMR voting decision from the number of significant supporting replicate combinations.

    If auto_vote_threshold is an integer, use that required_count directly;
    otherwise use the round-half-up required_count implied by --vote-threshold.
    """
    total_count = sig_count + nonsig_count
    if total_count <= 0:
        return 0

    if auto_vote_threshold is not None:
        required_count = int(auto_vote_threshold)
    else:
        required_count = int(np.floor(VOTE_THRESHOLD * total_count + 0.5))
    final_decision = 1 if sig_count >= required_count else 0
    # print(f"\nDecision(voting threshold={VOTE_THRESHOLD * 100:.1f}%)")
    # print(f" Decisionresult:{'significant' if final_decision else 'not significant'}")
    # print(f" support ratio:{support_ratio * 100:.1f}%")

    return final_decision

def find_common_significant_sites(output_dirs=None, methytype2='CpG', dir1=None, dir2=None, work_dir=".", meth_diff_threshold=0.0):
    """
    Find sites that are significant across the combination tests and collect related information.

    Parameters:
        output_dirs: list of output directories; if None, directories are scanned automatically.
        methytype2: methylation context.
    """

    print("\nSearching for common significant sites across all combinations...")

    and_output_dir = os.path.join(work_dir, "and_output")
    os.makedirs(and_output_dir, exist_ok=True)

    # 1. Automatically scan all output_x_y directories
    if output_dirs is None:
        output_dirs = glob.glob(os.path.join(work_dir,f"output_*_*/{methytype2}"))
        # Find all directories of the form output_x_y/methytype2, for a total of m*n
        output_dirs = [d for d in output_dirs if os.path.isdir(d)]

    if not output_dirs: # If none are found
        print("No output directories were found")
        return None

    print(f"Found {len(output_dirs)} output directories")

    auto_vote_threshold = AUTO_DMP_VOTE_THRESHOLDS.get(methytype2)
    if AUTO_VOTE_THRESHOLD_REPORT_ONLY:
        auto_vote_threshold = None
    if auto_vote_threshold is not None:
        print(f"Using {methytype2} DMP auto vote thresholds: {auto_vote_threshold}")
    else:
        print(f"{methytype2} DMPauto vote threshold is unavailable or report-only,using proportional threshold {VOTE_THRESHOLD}")

    # 2. Read all FDR_corrected files into memory at once, including all sites; if only significant-site files were read,,
            # when a site appears in the 1_1 significant file but not in a later test, it would be unclear whether it was nonsignificant or absent from the raw input for that test
    valid_dirs = [] # Collect paths to directories used in later operations
    site_statistics = {}  # Build this mapping:{site_id: {'sig': 0, 'total': 0}}
    all_dataframes = {} # Finally, the corresponding FDR_corrected DataFrame can be accessed as all_dataframes[directory]
    dir_to_replicate = {} # Record replicate IDs corresponding to each directory
    for output_dir in output_dirs: # Iterate over all output_x_y/methytype2 directories
        # Get the FDR_corrected file under the current methylation directory,
            # Format::'Chromosome', 'Methylation_Type', 'Position', 'Pvalue', 'Qvalue'
        # Here qvalue includes all sites rather than only values below the threshold
        match = re.search(r'output_wt(\d+)_mut(\d+)', output_dir)
        if not match:
            continue

        replicate_y = int(match.group(1))
        replicate_x = int(match.group(2))


        fdr_all_files = glob.glob(os.path.join(output_dir, "FDR_corrected_results_*.txt"))

        # Check existence
        if not fdr_all_files:
            print(f"  WARNING: {output_dir} not found inFDR_correctedfile")
            continue

        fdr_all_file = fdr_all_files[0] # Because glob.glob returns a list, use [0] to get the actual existing file path

        try:
            df = pd.read_csv(fdr_all_file, sep=r'\s+') # Read this file into df; FDR_corrected format is:
                                        # 'Chromosome', 'Methylation_Type', 'Position', 'Pvalue', 'Qvalue'
            if len(df) > 0:
                # Create a unique identifier for each row
                # Specifically, for each row in df, extract chromosome, methylation type, and position and add them as a new site_id column
                # This is needed to count whether each site is significant across tests; site_id distinguishes sites with one attribute, avoiding repeated
                # checks of multiple attributes for filtering or access (chromosome, methylation type, position); note that chromosome number here is the three-column block index in the both file
                df['site_id'] = df.apply(
                    lambda row: f"{int(row['Chromosome'])}-{row['Methylation_Type']}-{int(row['Position'])}",
                    axis=1
                )
                all_dataframes[output_dir] = df # all_dataframes[directory]->FDR_correctcorresponding df
                valid_dirs.append(output_dir)
                dir_to_replicate[output_dir] = (replicate_x, replicate_y)  # Record IDs
                # Get the default DMP threshold for the current methylation type; if the FDR file has the actual threshold used for each comparison, prefer that column
                dmp_threshold = get_dmp_threshold(methytype2)
                for _, row in df.iterrows(): # Iterate over each row; each row is one site from FDR_corrected for one test, regardless of significance
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


                print(f"  {output_dir}: {len(df)}  sites")
            else:
                print(f"  {output_dir}: no significant sites")
        except Exception as e:
            print(f"  ERROR:read {fdr_all_file} failed: {e}")
            continue
    valid_dirs.sort(key=lambda d: dir_to_replicate[d])
    print(f"\nstatistics completed,with {len(site_statistics)} distinct sites")

    if not site_statistics:
        print("no valid site-information results were found")
        return None

    # 3. Get sites called significant by the voting/Bayesian rule and put them into common_sites
    common_sites = []
    for site_id, stats in site_statistics.items():
        if stats['total_count'] != len(valid_dirs):
            stats['total_count'] = len(valid_dirs)
        sig_count = stats['sig_count']  # Get the significant-test count for the current site
        nonsig_count = stats['total_count'] - sig_count  # Get the nonsignificant-test count for the current site
        is_significant = bayes_deciding(sig_count, nonsig_count, auto_vote_threshold=auto_vote_threshold)
        if is_significant:
            common_sites.append(site_id)
    if not common_sites:
        print("no sites were significant in all combinations")
        return None
    else:
        print(f"\nafter Bayesian decision,with {len(common_sites)}  significant sitessites")

    # 4. Read methylation-change direction information from each directory
    print("reading methylation-change direction information...")
    methylation_change_by_dir = {}  # {output_dir: {site_id: change}}

    for output_dir in valid_dirs: # Iterate over valid directories for this methylation type
        methylation_change_by_dir[output_dir] = {} # Create the dictionary entry for the current directory,
                                        # Build the output_dir -> site_id -> change mapping
        # Find all all_simple_Chr files in this directory because they contain change information in the format pos, pvalue, change,
                                                # and chromosome number is obtained from the file name
        all_simple_files = glob.glob(os.path.join(output_dir, "all_simple_Chr*.txt"))
                                                        # Get all all_simple file paths under the current directory as a list

        for file_path in all_simple_files: # Iterate over all_simple files; format is pos, pvalue, change
            chr_match = re.search(r'Chr(\d+)\.txt$', file_path) # Create the regex and capturing group used to get chromosome number
            if not chr_match:
                continue
            chr_num = int(chr_match.group(1)) # Get the chromosome number of the current file from the capturing group

            try:
                df = pd.read_csv(
                    file_path,
                    sep='\t',
                    usecols=[0, 2],
                    names=['position', 'change'],
                    dtype={'position': float, 'change': float}, # Converting via float can handle more string formats robustly
                    skiprows=1
                )

                # Convert both pos and change to int at this point to avoid errors such as int("100.0")
                df['position'] = df['position'].astype(int)
                df['change'] = df['change'].astype(int)

                positions = df['position'].values
                changes = df['change'].values

                for position, change in zip(positions, changes):
                    site_id = f"{chr_num}-{methytype2}-{position}"
                    methylation_change_by_dir[output_dir][site_id] = change

            except Exception as e:
                print(f"  WARNING: read {file_path} failed: {e}")
                continue

        print(f"  {output_dir}: read {len(methylation_change_by_dir[output_dir])}  siteschange-direction records")

    # 5. Process common-site information
    print("processing detailed common-site information...")
    common_site_details = []

    # Create an index for each DataFrame to speed up lookup
    indexed_dfs = {} # Build the output_dir -> df1 mapping after site_id is set as the index
    for output_dir in valid_dirs:
        df = all_dataframes[output_dir] # Get the FDR_corrected DataFrame from each directory, already with the site_id column added
        indexed_dfs[output_dir] = df.set_index('site_id') # Set site_id as the index and return the new df1 as the dictionary value

    # Process common sites one by one
    for i, site_id in enumerate(common_sites):
        if i % 10000 == 0:  # Print progress every 10000 sites
            print(f"  Processed {i}/{len(common_sites)}  sites")

        site_info = {'site_id': site_id}
        chr_num, mtype, pos = site_id.split('-')
        site_info['Chromosome'] = int(chr_num)
        site_info['Methylation_Type'] = mtype
        site_info['Position'] = int(pos)

        # Collect q-values for the same chromosome and position across output_x_y/methytype2 directories for the current site_id
        qvalues = []
        sig_qvalues_for_mean = []
        # Collect change directions for the same chromosome and position across output_x_y/methytype2 directories for the current site_id
        change_values = []
        qvalue_dict = {}
        for output_dir in valid_dirs:
            replicate_x, replicate_y = dir_to_replicate[output_dir]
            col_name = f'qvalue_{os.path.basename(dir2.rstrip("/"))}{replicate_y}_{os.path.basename(dir1.rstrip("/"))}{replicate_x}'   # New column name: control/wild-type index first, mutant index second
            indexed_df = indexed_dfs[output_dir] # output_dir->df1(after site_id is used as the index),get the current methylation FDR_corrected file indexed by site_id
                                # Content format::'Chromosome', 'Methylation_Type', 'Position', 'Pvalue', 'Qvalue'
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

                # Get the change value of this site in this directory, (output_dir->site_id->changemapping)
                if site_id in methylation_change_by_dir[output_dir]:
                    change_values.append(methylation_change_by_dir[output_dir][site_id])
            else:
                qvalue_dict[col_name] = 1.0

        if qvalues:  # Ensure q-value data are available
            if len(sig_qvalues_for_mean) > 0:
                site_info['Sig_Mean_Qvalue'] = np.mean(sig_qvalues_for_mean)
            else:
                site_info['Sig_Mean_Qvalue'] = 1
            # site_info['Max_Qvalue'] = np.max(qvalues)
            # site_info['Min_Qvalue'] = np.min(qvalues)
            site_info['Num_Comparisons'] = len(qvalues)

            # Use voting to calculate methylation-change direction
            if change_values:
                # Count change == 1
                num_hyper = sum(change_values)
                total_comparisons = len(change_values)
                hyper_ratio = num_hyper / total_comparisons

                # Majority voting: record 1 (hypermethylated) if >= 50%, otherwise 0
                site_info['Methylation_Change'] = 1 if hyper_ratio >= 0.5 else 0
                site_info['Hyper_Count'] = num_hyper  # Hypermethylated count
                site_info['Hypo_Count'] = total_comparisons - num_hyper  # Hypomethylated count
                site_info['Hyper_Ratio'] = hyper_ratio  # Hypermethylated ratio
            else:
                # If no change information is available, mark it as missing; in principle, each test should have exactly one q-value and one change value
                site_info['Methylation_Change'] = -1  # -1 means undetermined
                site_info['Hyper_Count'] = 0
                site_info['Hypo_Count'] = 0
                site_info['Hyper_Ratio'] = 0

            site_info.update(qvalue_dict) # Add all replicate q-value columns

            common_site_details.append(site_info) # This list stores dictionaries; each dictionary contains information for one site_id with the following format:
        # site_id-chromosome number-methylation type-position-total test count-change-hypercount-hypocount-hyperratio-Sig_Mean_Qvalue-all q-values
            # Here chromosome number is the index of the three-column block in the both file

    print(f"  finished processing {len(common_site_details)}  sites")

    # 6. Generate result DataFrame
    result_df = pd.DataFrame(common_site_details)
    result_df = result_df.sort_values(['Methylation_Type', 'Chromosome', 'Position']) # Sort

    # Optional final-DMP low-difference strict-vote post-processing.
    # Execute separately by methylation type. The cutoff can be fixed globally at 0.3, but support_count, base_required,
    # boundary_abs_methdiff, and diagnostic tables must be calculated separately for CpG/CHG/CHH.
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

    # Adjust column order and move Methylation_Change earlier
    column_order = [
        'Chromosome', 'Methylation_Type', 'Position',
        'Methylation_Change', 'Hyper_Ratio', 'Hyper_Count', 'Hypo_Count', 'Num_Comparisons',
        'Sig_Mean_Qvalue'
    ]
    # Get and sort all replicate column names
    # New column-name format: dir2_y_dir1_x; need to get columns containing underscore-separated numbers
    replicate_columns = sorted(
        [col for col in result_df.columns if col.startswith('qvalue_')],
        key=lambda x: tuple(map(int, re.findall(r'\d+', x)))  # Extract all numbers and sort
    )
    column_order = column_order + replicate_columns
    result_df = result_df[column_order]

    # Save results
    output_file = os.path.join(and_output_dir, f"{methytype2}-final_significant_sites_DMPs.txt")
    result_df.to_csv(output_file, sep='\t', index=False)
    print(f"\nCommon significant sites have been saved to: {output_file}")

    # Print summary statistics
    print("\nCommon significant-site statistics:")
    mtype = methytype2
    mtype_df = result_df
    count = len(mtype_df)
    if count == 0:
        print(f"  {mtype}: 0  sites")
        return result_df
    hyper_count = len(mtype_df[mtype_df['Methylation_Change'] == 1])
    hypo_count = len(mtype_df[mtype_df['Methylation_Change'] == 0])
    unknown_count = len(mtype_df[mtype_df['Methylation_Change'] == -1])

    print(f"  {mtype}: {count}  sites")
    print(f"    - hypermethylated(Change=1): {hyper_count} ({hyper_count / count * 100:.1f}%)")
    print(f"    - hypomethylated(Change=0): {hypo_count} ({hypo_count / count * 100:.1f}%)")
    if unknown_count > 0:
        print(f"    - unknown(Change=-1): {unknown_count} ({unknown_count / count * 100:.1f}%)")

    return result_df # Its format is:'Chromosome', 'Methylation_Type', 'Position',
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
    Run sliding-window analysis on methylation-site data.

    Parameters:
        input_data: DataFrame converted from one chromosome-level (N)DMP file from one of the m*n*3 tests;
            it contains columns ['position', 'pvalue', 'change'].
        window_size: sliding-window size.
        step_ratio: step-size ratio as a percentage of the window size.
        save_files: whether to save result files.
        output_identifier: output-file prefix; if save_files=True and this is not provided, it is generated automatically.
    """

    # 1. Data loading and preprocessing
    df = input_data.copy()

    # Ensure column names are correct
    expected_cols = ['position', 'pvalue', 'change']
    if not all(col in df.columns for col in expected_cols):
        raise ValueError(f"DataFramemust contain columns: {expected_cols}")

    if output_identifier is None:
        output_identifier = "sliding_window_analysis"

    # Validate and clean data
    df = df.dropna()
    if len(df) == 0:
        raise ValueError("no valid data rows")

    # Sort by position
    df = df.sort_values('position').reset_index(drop=True) # Reset the DataFrame index without keeping the old index as a new column

    print(f"Data preprocessing completed,with {len(df)}  sites")

    # 2. Sliding-window analysis
    positions = df['position'].values # Get all positions; they have already been sorted
    changes = df['change'].values # Get all methylation-change directions

    last_pos = positions[-1] # Get the last position
    step_size = int(window_size * step_ratio) # Calculate the step size for each window shift

    # Generate start positions for all windows
    window_starts = np.arange(0, last_pos + step_size, step_size) # The +step_size ensures the last interval covers last_pos

    print(f"Window configuration: size={window_size}, step size={step_size}, number of windows={len(window_starts)}")

    # Create a list to store start/end positions, counts of each methylation-change direction, and total significant-site counts for each interval
    results = []

    for i, start in enumerate(window_starts): # Each start is set to the window start position minus 1 (e.g., 0 for the first window)
        if i % 1000 == 0:  # Progress message
            print(f"processing progress: {i}/{len(window_starts)}")

        end = start + window_size  # Calculate the window end from the window width

        # Use numpy.searchsorted for fast lookup, equivalent to binary search
        left_idx = np.searchsorted(positions, start, side='right') # Find the index of the first element in positions greater than start
        right_idx = np.searchsorted(positions, end, side='right') # Similarly, find the index of the first element in positions greater than end

        # Sites inside the window
        window_changes = changes[left_idx:right_idx] # Get the change array for positions in the [start, end) window of width window_size
                                                    # Note that left_idx and right_idx are both indices
                            # The left_idx:right_idx range contains positions with values > start and <= end,
                         # but the number of records is usually much smaller than window_size because many positions have no data
        # Count categories
        num_change_1 = np.sum(window_changes == 1)   # Count hypermethylated sites
        num_change_0 = np.sum(window_changes == 0)  # Count hypomethylated sites, or use len(window_changes) - num_change_1

        results.append({
            'window_start': start + 1,  # start + 1 is the true genomic start position of each window
            'window_end': end,   # Positions from start + 1 to end cover exactly window_size bases
            'count_change_1': num_change_1,
            'count_change_0': num_change_0,
            'total_count': num_change_1 + num_change_0 # Methylation-change directions for all sites in the current interval (mutant relative to wild type)
                                                # This is also the total significant-site count for the current interval because each test has one change direction
                                            # Since records are read from a DMP file, all are significant
        })

    # Convert to DataFrame
    sliding_results = pd.DataFrame(results)

    # 3. Standardization
    max_count = sliding_results['total_count'].max()  # Maximum total significant-site count across all intervals
    if max_count == 0:
        max_count = 1  # Avoid division by zero

    standardized_results = sliding_results.copy()
    standardized_results['standardized_count'] = sliding_results['total_count'] / max_count # Calculate the ratio of the current interval total significant-site count to the
                                                                                        # maximum significant-site count across intervals
    # Select required columns for standardized output
    standardized_results = standardized_results[[
        'window_start', 'window_end', 'total_count', 'standardized_count'
    ]]

    print(f"Sliding-window analysis completed,generated {len(sliding_results)} windows")
    print(f"maximum count: {max_count}")

    # 4. Save files
    if save_files:
        # Sliding-window results
        sliding_file = f"slidingW_{output_identifier}.txt"
        sliding_file = os.path.join(outputdir1, sliding_file)
        sliding_results[['window_start', 'window_end', 'count_change_1', 'count_change_0']].to_csv(
            sliding_file,
            sep='\t',
            index=False,
            header=False
        )

        # Standardized results
        std_file = f"noTitle_allDMCs_new_Standardized_slidingW_{output_identifier}.txt"
        std_file = os.path.join(outputdir1, std_file)
        # Format output to match the original C++-style alignment
        standardized_results.to_csv(
            std_file,
            sep='\t',
            index=False,
            header=False,
            float_format='%.6f'  # Control floating-point precision
        )

        print(f"Results saved:")
        print(f"  sliding window: {sliding_file}")
        print(f"  standardized results: {std_file}")

    return sliding_results, standardized_results
    # The former format is window_start, window_end, count_change_1, count_change_0, total_count
                        # The latter has no count_change but has standardized_count, the ratio of significant-site count in each interval to the maximum interval count

def process_common_sites_sliding_window(common_sites_df=None,
                                        window_size=1000000,
                                        step_ratio=0.05,
                                        methytype='CpG',
                                        work_dir="."):
    """
    Run sliding-window analysis on common significant sites.

    Parameters:
        common_sites_df: common significant-site data; if None, it is loaded automatically.
        window_size: sliding-window size.
        step_ratio: step-size ratio.
        methytype: methylation context.
    """

    print(f"\nStarting sliding-window analysis for common significant sites...")

    and_output_dir = os.path.join(work_dir, "and_output")
    os.makedirs(and_output_dir, exist_ok=True)

    # 1. Load common significant-site data
    if common_sites_df is None:
        common_file = os.path.join(and_output_dir, f"{methytype}-final_significant_sites_DMPs.txt")
        if not os.path.exists(common_file):
            print(f"ERROR:common significant-site file does not exist {common_file}")
            return None
        common_sites_df = pd.read_csv(common_file, sep='\t')
        print(f"Loaded common significant sites from file: {len(common_sites_df)}  sites")

    if common_sites_df.empty:
        print("No common significant-site data")
        return None

    # 2. Group by chromosome
    results = {}

    # Process by chromosome group
    chr_groups = common_sites_df.groupby('Chromosome') # Get an iterator over chromosome labels and their corresponding sub-dataframes

    for chr_num, chr_data in chr_groups: # Iterate over chromosome labels and their corresponding sub-dataframes
        print(f"\n  Processing chromosome {chr_num}: {len(chr_data)}  sites")

        # Prepare data for sliding-window analysis in the same format as all_simple_chr files to ensure compatibility
        window_data = pd.DataFrame({
            'position': chr_data['Position'].astype(int),
            'pvalue': chr_data['Sig_Mean_Qvalue'],
            'change': chr_data['Methylation_Change']
        })

        # Sort
        window_data = window_data.sort_values('position').reset_index(drop=True)

        # Run sliding-window analysis by calling the existing function
        try:
            sliding_results, std_results = sliding_window_analysis(
                window_data,
                window_size=window_size,
                step_ratio=step_ratio,
                save_files=True,
                output_identifier=f"common_sites_{methytype}_Chr{chr_num}",
                outputdir1=and_output_dir
            )

            # Collect result files
            results[chr_num] = {
                'sliding_results': sliding_results,
                'standardized_results': std_results,
                'input_data': window_data
            }

            print(f"    completed chromosome {chr_num}: {len(sliding_results)} windows")

        except Exception as e:
            print(f"    ERROR:Processing chromosome {chr_num} failed: {e}")
            continue

    print(f"Results are saved in files starting with 'common_sites_{methytype}_Chr' ")

    return results

def find_max_total_in_outputs(output_dirs, methylation_type):
    """
    Find the maximum total value for a specified methylation type across all output directories.

    Parameters:
        output_dirs: list of output directories.
        methylation_type: methylation context (CpG, CHH, CHG).

    Returns:
        max_total: maximum total value.
    """
    max_total = 0

    for out_dir in output_dirs:
        mtype_dir = os.path.join(out_dir, methylation_type)
        if not os.path.exists(mtype_dir):
            continue

        # Find all standardized DMP files
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
                print(f"    WARNING: read {std_file}:  {e}")
                continue

    print(f"  {methylation_type} largest value across all chromosomes for contexttotalis: {max_total}")
    return max_total


def plot_methylation_sliding_windows(output_dir=None, chr_series=None,work_dir="."):
    """
    Visualize all sliding-window results and save the figures to disk.
    Global max_total is used for normalization so that different chromosomes are comparable.

    Parameters:
        output_dir: specified output directory; if None, all output_x_y directories are scanned automatically.
        chr_series: chromosome mapping Series.
    """

    matplotlib.use('Agg')  # Use a non-interactive backend

    # Set Chinese font
    #plt.rcParams['font.sans-serif'] = ['SimHei']
    plt.rcParams['axes.unicode_minus'] = False

    print("\nStarting methylation sliding-window visualization plot generation")

    # Scan all output directories
    if output_dir is None:
        output_dirs = glob.glob(os.path.join(work_dir, "output_*_*"))
        output_dirs = [d for d in output_dirs if os.path.isdir(d)]
    else:
        output_dirs = [output_dir]

    if not output_dirs:
        print("No output directories found")
        return

    total_plots = 0
    methylation_types = ['CpG', 'CHH', 'CHG']

    # Calculate global max_total separately for each methylation type
    max_totals = {}
    for mtype in methylation_types:
        max_totals[mtype] = find_max_total_in_outputs(output_dirs, mtype)
        if max_totals[mtype] == 0:
            print(f"  WARNING: {mtype} context has no validtotalvalue,will use1as the default value")
            max_totals[mtype] = 1

    for out_dir in output_dirs:
        print(f"\nProcessing directory: {out_dir}")

        for mtype in methylation_types:
            mtype_dir = os.path.join(out_dir, mtype)
            if not os.path.exists(mtype_dir):
                continue

            print(f"  Processing methylation context: {mtype}(using global max_total={max_totals[mtype]})")

            # Find all DMP sliding-window files under the current methylation directory
            dmp_sliding_files = glob.glob(os.path.join(mtype_dir, "slidingW_DMP_*.txt"))

            # Group by prefix and draw all chromosomes with the same prefix in one large figure
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

            # Process each prefix group
            for prefix, chr_nums in prefix_groups.items():
                # Sort chromosome numbers
                chr_nums = sorted(chr_nums, key=lambda x: int(x))

                # Collect data for all chromosomes
                all_chrom_data = []
                chrom_names = []

                for chr_num in chr_nums:
                    # Construct the corresponding standardized-file path
                    chr_name = get_chr_name(chr_num, chr_series)

                    # Fix: use the correct variable name to construct the file path
                    dmp_sliding_file = os.path.join(mtype_dir, f"slidingW_DMP_{prefix}_Chr{chr_num}.txt")
                    dmp_std_file = os.path.join(mtype_dir,
                                                f"noTitle_allDMCs_new_Standardized_slidingW_DMP_{prefix}_Chr{chr_num}.txt")
                    ndmp_std_file = os.path.join(mtype_dir,
                                                 f"noTitle_allDMCs_new_Standardized_slidingW_N-DMP_{prefix}_Chr{chr_num}.txt")

                    if not all(os.path.exists(f) for f in [dmp_sliding_file, dmp_std_file, ndmp_std_file]):
                        print(f"    WARNING: Chr{chr_num} file is incomplete; skipping")
                        continue

                    try:
                        # Read data
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
                            print(f"    WARNING: Chr{chr_num} data are empty; skipping")
                            continue

                        # Recalculate ratios using global max_total
                        max_total = max_totals[mtype]

                        x = (sliding_df['start'] + sliding_df['end']) / 2

                        # Recalculate all ratios using global max_total
                        y_dmp = (sliding_df['hyper'] + sliding_df['hypo']) / max_total
                        y_hyper = sliding_df['hyper'] / max_total
                        y_hypo = sliding_df['hypo'] / max_total
                        y_ndmp = ndmp_std_df['ndmp_normalized']  # Keep NDMP unchanged

                        # Handle possible length mismatch
                        max_len = max(len(x), len(y_dmp), len(y_hyper), len(y_hypo), len(y_ndmp))
                        x = x.reindex(range(max_len), fill_value=0)
                        y_dmp = y_dmp.reindex(range(max_len), fill_value=0)
                        y_hyper = y_hyper.reindex(range(max_len), fill_value=0)
                        y_hypo = y_hypo.reindex(range(max_len), fill_value=0)
                        y_ndmp = y_ndmp.reindex(range(max_len), fill_value=0)

                        # Storedata
                        all_chrom_data.append({
                            'x': x,
                            'y_dmp': y_dmp,
                            'y_hyper': y_hyper,
                            'y_hypo': y_hypo,
                            'y_ndmp': y_ndmp
                        })
                        chrom_names.append(chr_name)

                        print(f"    successfully loadedchromosome {chr_name}  data")

                    except Exception as e:
                        print(f"    Processing Chr{chr_num}:  {e}")
                        continue

                # If data are available, draw the large figure
                if all_chrom_data:
                    try:
                        # Create a large figure with one subplot per chromosome
                        n_chromosomes = len(all_chrom_data)
                        fig, axes = plt.subplots(n_chromosomes, 1, figsize=(15, 3 * n_chromosomes))

                        # If there is only one chromosome, axes is not an array and must be converted to an array
                        if n_chromosomes == 1:
                            axes = [axes]

                        # Set figure title
                        fig.suptitle(f'{mtype} Methylation Analysis - {prefix} (Global Normalized)',
                                     fontsize=16, fontfamily='DejaVu Sans')

                        # Draw the subplot for each chromosome
                        for idx, (chrom_data, chrom_name) in enumerate(zip(all_chrom_data, chrom_names)):
                            ax = axes[idx]

                            # Draw all data lines
                            ax.plot(chrom_data['x'], chrom_data['y_dmp'], label='DMP', color='red', linewidth=1.5)
                            ax.plot(chrom_data['x'], chrom_data['y_hyper'], label='Hyper-ratio', color='green',
                                    linewidth=1)
                            ax.plot(chrom_data['x'], chrom_data['y_hypo'], label='Hypo-ratio', color='blue',
                                    linewidth=1)
                            ax.plot(chrom_data['x'], chrom_data['y_ndmp'], label='NDMP', color='darkgray',
                                    linewidth=1.5)

                            ax.set_ylim(bottom=0)
                            ax.set_ylabel('Ratio', fontsize=10, fontfamily='DejaVu Sans')

                            # Set subplot title
                            ax.set_title(f'{chrom_name}', fontsize=18, fontfamily='DejaVu Sans', pad=20,y=-0.4)

                            # Add grid
                            ax.grid(True, alpha=0.3)

                            # Add the full legend only to the first subplot
                            if idx == 0:
                                ax.legend(loc='upper right', ncol=2, fontsize=8, framealpha=0.7)

                        # Adjust layout
                        plt.tight_layout(rect=[0, 0, 1, 0.96])  # Leave space for the figure title

                        # Save figure
                        plot_filename = os.path.join(mtype_dir,
                                                     f"methylation_plot_{mtype}_{prefix}_all_chromosomes.png")
                        plt.savefig(plot_filename, dpi=300, bbox_inches='tight')
                        plt.close()

                        total_plots += 1
                        print(f"    successfully generated large plot: {mtype}_{prefix} -> {os.path.basename(plot_filename)}")

                    except Exception as e:
                        print(f"    error while drawing large plot: {e}")
                        continue

    print(f"\nPlot generation completed!generated {total_plots} large plots")


def plot_common_sites_sliding_windows(methytype='CpG', chr_series=None, work_dir="."):
    """
    Visualize sliding-window results for common significant sites and save the figures to disk.
    Global max_total is used for normalization, and all chromosomes are drawn in one large figure.
    """
    matplotlib.use('Agg')

    # Set the global font to DejaVu Sans
    plt.rcParams['font.family'] = 'DejaVu Sans'
    plt.rcParams['axes.unicode_minus'] = False

    print(f"\nStarting common significant-site({methytype})sliding-window visualization plot generation...")

    and_output_dir = os.path.join(work_dir, "and_output")

    # First find the global max_total for this methylation type
    print(f"  searching for {methytype}  global max_total...")
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
            print(f"    WARNING: read {std_file}:  {e}")

    if max_total == 0:
        print(f"  WARNING: no validtotalvalue,using default value1")
        max_total = 1
    else:
        print(f"  {methytype}  global max_totalfor: {max_total}")

    # Read DMR data and choose the corresponding DMR file according to methylation type
    print("  reading DMR data...")
    dmr_file = os.path.join(and_output_dir, f"{methytype}-final_significant_regions_DMRs.txt")
    dmr_data = {}
    if os.path.exists(dmr_file):
        try:
            dmr_df = pd.read_csv(dmr_file, sep=r'\s+')

            for _, row in dmr_df.iterrows():
                try:
                    chrom = str(row['Chromosome'])  # Access by column name
                    direction = int(row['Direction'])  # Access by column name
                    start = int(row['DMR_start'])  # Access by column name
                    end = int(row['DMR_end'])  # Access by column name

                    # Calculate midpoint
                    mid = (start + end) / 2

                    # Extract numeric part of chromosome label
                    chrom_num = str(chrom).replace('Chr', '').replace('chr', '')

                    if chrom_num not in dmr_data:
                        dmr_data[chrom_num] = []

                    dmr_data[chrom_num].append((mid, direction))
                except (ValueError, IndexError):
                    continue
            print(f"  successfully loaded {sum(len(dmrs) for dmrs in dmr_data.values())} DMR")
        except Exception as e:
            print(f"  readDMR file error: {e}")
    else:
        print(f"  WARNING: DMR file {dmr_file} does not exist")

    # Get chromosome list dynamically
    sliding_files = glob.glob(os.path.join(and_output_dir, f"slidingW_common_sites_{methytype}_Chr*.txt"))

    # Extract chromosome numbers from file names and sort
    chromosomes = []
    for file in sliding_files:
        match = re.search(r'slidingW_common_sites_.+_Chr(\d+)\.txt$', os.path.basename(file))
        if match:
            chr_num = match.group(1)
            if chr_num not in chromosomes:
                chromosomes.append(chr_num)

    # Sort chromosomes numerically
    chromosomes.sort(key=int)

    if not chromosomes:
        print(f"No common significant sites found for({methytype})ofsliding windowfile")
        return

    print(f"Found {len(chromosomes)}  chromosomesofsliding windowfile")

    all_chrom_data = []
    chrom_nums1 = []

    for chr_num in chromosomes:
        sliding_file = os.path.join(and_output_dir, f"slidingW_common_sites_{methytype}_Chr{chr_num}.txt")
        std_file = os.path.join(and_output_dir,
                                f"noTitle_allDMCs_new_Standardized_slidingW_common_sites_{methytype}_Chr{chr_num}.txt")

        if not all(os.path.exists(f) for f in [sliding_file, std_file]):
            print(f"    WARNING: Chr{chr_num} file is incomplete; skipping")
            continue

        try:
            # Read data
            sliding_df = pd.read_csv(sliding_file, sep=r'\s+', header=None, names=['start', 'end', 'hyper', 'hypo'])
            std_df = pd.read_csv(std_file, sep=r'\s+', header=None, names=['start', 'end', 'total', 'normalized'])

            if sliding_df.empty or std_df.empty:
                print(f"    WARNING: Chr{chr_num} data are empty; skipping")
                continue

            if len(sliding_df) != len(std_df):
                print(f"    WARNING: Chr{chr_num} data lengths are inconsistent; skipping")
                continue

            # New: read NDMP data from output_1_1
            ndmp_file = os.path.join(work_dir, "output_wt1_mut1", methytype,
                                     f"noTitle_allDMCs_new_Standardized_slidingW_N-DMP_wt_replicate1_mut_replicate1_Chr{chr_num}.txt")
            y_ndmp = None  # Initialize as None
            if os.path.exists(ndmp_file):
                try:
                    ndmp_df = pd.read_csv(ndmp_file, sep=r'\s+', header=None,
                                          names=['start', 'end', 'total', 'ndmp_normalized'])
                    if not ndmp_df.empty:
                        y_ndmp = ndmp_df['ndmp_normalized']
                        print(f"    successfully read output_1_1 NDMP data: Chr{chr_num}")
                except Exception as e:
                    print(f"    WARNING: read output_1_1 NDMP data failed (Chr{chr_num}): {e}")
            else:
                print(f"    Note: output_1_1 of NDMP filedoes not exist (Chr{chr_num})")

            # Recalculate ratios using global max_total
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

            # Storedata
            chr_real_name = get_chr_name(chr_num, chr_series)
            all_chrom_data.append({
                'x': x,
                'y_total': y_total,
                'y_hyper': y_hyper,
                'y_hypo': y_hypo,
                'y_ndmp': y_ndmp
            })
            chrom_nums1.append(chr_num)
            print(f"    successfully loadedchromosome {chr_real_name}  data")

        except Exception as e:
            print(f"    Processing Chr{chr_num}:  {e}")
            continue

    # If data are available, draw the large figure
    if all_chrom_data:
        try:
            # Create a large figure with multiple subplots
            n_chromosomes = len(all_chrom_data)
            fig, axes = plt.subplots(n_chromosomes, 1, figsize=(15, 3 * n_chromosomes))

            # If there is only one chromosome, axes is not an array and must be converted to an array
            if n_chromosomes == 1:
                axes = [axes]

            # Set figure title
            fig.suptitle(f'Distribution of Common Significant Sites - {methytype} context',
                         fontsize=16, fontfamily='DejaVu Sans')

            # Draw the subplot for each chromosome
            for idx, (chrom_data, chrom_num1) in enumerate(zip(all_chrom_data, chrom_nums1)):
                ax = axes[idx]

                # Draw DMP data
                ax.plot(chrom_data['x'], chrom_data['y_total'], label='DMP', color='red', linewidth=2)
                ax.plot(chrom_data['x'], chrom_data['y_hyper'], label='hyper-methylation', color='green', linewidth=1.5)
                ax.plot(chrom_data['x'], chrom_data['y_hypo'], label='hypo-methylation', color='blue', linewidth=1.5)


                if chrom_data['y_ndmp'] is not None:
                    ax.plot(chrom_data['x'], chrom_data['y_ndmp'], label='NDMP',
                            color='darkgray', linewidth=1.5)
                # Set a unified y-axis range from 0 to 1.2
                ax.set_ylim(0, 1.2)

                # Get chromosome number extracted from chrom_name
                chrom_num1 = chrom_num1.replace('chr', '').replace('Chr', '')

                if idx == 0:  # Print only once for the first subplot
                    print(f"  debug: dmr_data keys = {list(dmr_data.keys())}")
                print(f"    {get_chr_name(chrom_num1,chr_series)} -> chrom_num = '{chrom_num1}', indmr_datain: {chrom_num1 in dmr_data}")

                # Add DMR markers
                if chrom_num1 in dmr_data:
                    for mid, direction in dmr_data[chrom_num1]:
                        # Choose color by direction: 1 = hyper, 0 = hypo
                        color = 'green' if direction == 1 else 'blue'
                        # Add vertical lines at DMR midpoints, shown within y-axis range 1.0 to 1.2
                        ax.axvline(x=mid, ymin=0.9, ymax=1, color=color, linewidth=2, alpha=0.7)

                # Set subplot titles and labels
                ax.text(0.5, -0.2, f"{get_chr_name(chrom_num1,chr_series)}",
                        transform=ax.transAxes,
                        fontfamily='DejaVu Sans',
                        ha='center', va='top',
                        fontsize=15)

                # Add grid
                ax.grid(True, alpha=0.3)

                # Add the full legend only to the first subplot
                if idx == 0:
                    # Create custom legend entries including DMR markers
                    from matplotlib.lines import Line2D
                    legend_elements = [
                        Line2D([0], [0], color='red', linewidth=2, label='DMP'),
                        Line2D([0], [0], color='green', linewidth=1.5, label='hyper-DMP'),
                        Line2D([0], [0], color='blue', linewidth=1.5, label='hypo-DMP'),
                        Line2D([0], [0], color='green', linewidth=2, label='hyper-DMR'),
                        Line2D([0], [0], color='blue', linewidth=2, label='hypo-DMR')
                    ]
                    # If NDMP data are available, add them to the legend
                    if chrom_data['y_ndmp'] is not None:
                        legend_elements.append(
                            Line2D([0], [0], color='darkgray', linewidth=1.5,
                                  label='NDMP')
                        )
                    ax.legend(handles=legend_elements, loc='upper right', bbox_to_anchor=(1.0, 0.95),
                              ncol=5, fontsize=8, framealpha=0.7)

            # Adjust layout
            plt.tight_layout(rect=[0, 0, 1, 0.96])  # Leave space for the figure title

            # Save figure
            plot_filename = os.path.join(and_output_dir, f"common_sites_plot_{methytype}_all_chromosomes.png")
            plt.savefig(plot_filename, dpi=300, bbox_inches='tight')
            plt.close()

            print(f"successfully generated large plot: {methytype} -> {os.path.basename(plot_filename)}")

        except Exception as e:
            print(f"error while drawing large plot: {e}")
    else:
        print(f"No {methytype} valid data")



def rename_chromosome_files(chr_series, work_dir="."):
    """
    Batch-rename files containing chromosome indices by replacing Chr-number labels with real chromosome names.

    Parameters:
        chr_series: chromosome mapping Series (chromosome name -> index).
        work_dir: working directory.
    """
    print("\nStarting batch renaming of chromosome files...")

    if chr_series is None or len(chr_series) == 0:
        print("ERROR:chromosome mapping is empty; skipping renaming")
        return

    # Create reverse mapping: numeric index -> chromosome name
    # chr_series index is chromosome name and value is numeric index
    index_to_chr = {i: chr_series.index[i] for i in range(len(chr_series))}

    print(f"chromosome mapping: {index_to_chr}")

    # Define directory patterns to search
    search_dirs = []

    # Add all output_x_y directories
    output_dirs = glob.glob(os.path.join(work_dir, "output_*_*"))
    search_dirs.extend([d for d in output_dirs if os.path.isdir(d)])

    # Add and_output directory
    and_output_dir = os.path.join(work_dir, "and_output")
    if os.path.exists(and_output_dir):
        search_dirs.append(and_output_dir)

    if not search_dirs:
        print("No directories requiring processing were found")
        return

    print(f"Will search files in {len(search_dirs)} directories")

    # Summary statistics
    total_renamed = 0
    failed_renames = 0

    # Define regex pattern for matching chromosome numbers
    # Match patterns where Chr is followed by digits, such as Chr1 or Chr12
    chr_pattern = re.compile(r'(.*?)Chr(\d+)(.*?)$')

    # Iterate over all directories
    for search_dir in search_dirs:
        print(f"\nProcessing directory: {search_dir}")

        # Recursively iterate over all files in the directory
        for root, dirs, files in os.walk(search_dir):
            for filename in files:
                # Check whether the file name contains the Chr + digits pattern
                match = chr_pattern.match(filename)

                if match:
                    prefix = match.group(1)  # Part before Chr
                    chr_num = int(match.group(2))  # Chromosome number
                    suffix = match.group(3)  # Part after Chr + digits

                    # Get the real chromosome name from chr_series
                    # Note: Chr1 in file names corresponds to index 0
                    chr_index = chr_num - 1

                    if chr_index not in index_to_chr:
                        print(f"  WARNING: Chr{chr_num} is not in the mapping table; skippingfile {filename}")
                        continue

                    real_chr_name = index_to_chr[chr_index]

                    # Construct new file name
                    # If the chromosome name itself contains the 'chr' prefix, use it directly
                    # Otherwise use the Chr prefix
                    if real_chr_name.lower().startswith('chr'):
                        chr_part = real_chr_name
                    else:
                        chr_part = f"Chr{real_chr_name}"

                    new_filename = f"{prefix}{chr_part}{suffix}"

                    # Skip if the new and old file names are identical
                    if filename == new_filename:
                        continue

                    # Construct full path
                    old_path = os.path.join(root, filename)
                    new_path = os.path.join(root, new_filename)

                    # Check whether the target file already exists
                    if os.path.exists(new_path):
                        print(f"  WARNING: target file already exists; skipping renaming: {filename} -> {new_filename}")
                        failed_renames += 1
                        continue

                    # Perform renaming
                    try:
                        os.rename(old_path, new_path)
                        total_renamed += 1
                        print(f"  {filename} -> {new_filename}")
                    except Exception as e:
                        print(f" rename failed: {filename} -> {new_filename}, ERROR: {e}")
                        failed_renames += 1

    # Print summary statistics
    print(f"\nrenaming completed!")
    print(f"  successfully renamed: {total_renamed} files")
    if failed_renames > 0:
        print(f"  failed or skipped: {failed_renames} files")


def convert_output_to_csv(work_dir="."):
    """
    Convert final DMP and final DMR files under and_output to comma-separated CSV files.

    Parameters:
        work_dir: working directory.
    """

    and_output_dir = os.path.join(work_dir, "and_output")

    if not os.path.exists(and_output_dir):
        print(f"ERROR: directory {and_output_dir} does not exist")
        return 0

    # Define file patterns to convert
    file_patterns = [
        "*-final_significant_sites_DMPs.txt",  # final DMP files
        "*-final_significant_regions_DMRs.txt"  # final DMR files
    ]

    converted_count = 0
    failed_count = 0

    for pattern in file_patterns:
        # Search for matching files
        matching_files = glob.glob(os.path.join(and_output_dir, pattern))

        for txt_file in matching_files:
            try:
                # Read tab-delimited file
                df = pd.read_csv(txt_file, sep=r'\s+')

                if df.empty:
                    print(f"  Skipping empty file: {os.path.basename(txt_file)}")
                    continue

                # Generate CSV file path by replacing .txt with .csv
                csv_file = txt_file.replace('.txt', '.csv')

                # Save as comma-separated CSV file
                df.to_csv(csv_file, sep=',', index=False)

                print(f"  conversion successful: {os.path.basename(txt_file)} -> {os.path.basename(csv_file)}")
                converted_count += 1

            except Exception as e:
                print(f"  conversion failed: {os.path.basename(txt_file)}, ERROR: {e}")
                failed_count += 1

    print(f"\\nconversion completed!")
    print(f"  conversion successful: {converted_count} files")
    if failed_count > 0:
        print(f"  conversion failed: {failed_count} files")

    return converted_count


def convert_chromosome_to_names(chr_series, work_dir="."):
    """
    Convert the Chromosome column in final DMP and final DMR files under and_output
    from numeric indices to real chromosome names.

    Parameters:
        chr_series: chromosome mapping Series (chromosome name -> index, zero-based index).
        work_dir: working directory.

    Returns:
        Number of files successfully converted.
    """

    if chr_series is None or len(chr_series) == 0:
        print("ERROR:chromosome mapping is empty; skippingconversion")
        return 0

    and_output_dir = os.path.join(work_dir, "and_output")

    if not os.path.exists(and_output_dir):
        print(f"ERROR: directory {and_output_dir} does not exist")
        return 0

    # Create numeric-to-chromosome-name mapping
    # The Chromosome column in files is 1-based numeric and must be converted to chromosome names in chr_series
    # chr_series.index[0] corresponds to 1 in files
    index_to_chr = {i + 1: chr_series.index[i] for i in range(len(chr_series))}

    print(f"chromosome mapping table: {index_to_chr}")

    # Define file patterns to process
    file_patterns = [
        "*-final_significant_sites_DMPs.txt",  # final DMP files
        "*-final_significant_regions_DMRs.txt"  # final DMR files
    ]

    converted_count = 0
    failed_count = 0

    for pattern in file_patterns:
        # Search for matching files
        matching_files = glob.glob(os.path.join(and_output_dir, pattern))

        for file_path in matching_files:
            try:
                # Read file
                df = pd.read_csv(file_path, sep=r'\s+')

                if df.empty:
                    print(f"  Skipping empty file: {os.path.basename(file_path)}")
                    continue

                # Check whether the Chromosome column exists
                if 'Chromosome' not in df.columns:
                    print(f"  WARNING: {os.path.basename(file_path)} does not contain Chromosome column; skipping")
                    continue

                # Save the original Chromosome column for debugging
                original_chrs = df['Chromosome'].unique()

                # Convert the Chromosome column
                # Ensure integer type first
                df['Chromosome'] = df['Chromosome'].astype(int)

                # Convert to chromosome names using the mapping
                df['Chromosome'] = df['Chromosome'].map(index_to_chr)

                # Check whether any values were not mapped successfully
                if df['Chromosome'].isna().any():
                    unmapped_count = df['Chromosome'].isna().sum()
                    print(f"  WARNING: {os.path.basename(file_path)} contains {unmapped_count}  chromosome IDs could not be mapped")
                    # Optional: remove rows that could not be mapped
                    df = df.dropna(subset=['Chromosome'])

                # Save back to the original file, overwriting it
                df.to_csv(file_path, sep='\t', index=False)

                print(f"    conversion successful: {os.path.basename(file_path)}")
                print(f"    original IDs: {sorted(original_chrs)}")
                print(f"    after conversion: {sorted(df['Chromosome'].unique())}")
                converted_count += 1

            except Exception as e:
                print(f"  [FAIL] conversion failed: {os.path.basename(file_path)}")
                print(f"    error message: {e}")
                import traceback
                traceback.print_exc()
                failed_count += 1

    # Print summary statistics
    print(f"\nchromosome ID conversion completed!")
    print(f"  conversion successful: {converted_count} files")
    if failed_count > 0:
        print(f"  conversion failed: {failed_count} files")

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
    # ===== Backward-compatible legacy positional arguments:python script.py n m dir2 dir1 biotype =====
    parser.add_argument("n_pos", type=int, nargs="?", help=argparse.SUPPRESS)
    parser.add_argument("m_pos", type=int, nargs="?", help=argparse.SUPPRESS)
    parser.add_argument("dir2_pos", nargs="?", help=argparse.SUPPRESS)
    parser.add_argument("dir1_pos", nargs="?", help=argparse.SUPPRESS)
    parser.add_argument("biotype_pos", type=int, nargs="?", choices=[0, 1, 2], help=argparse.SUPPRESS)

    # ===== New named arguments:python script.py --n 2 --m 2 --dir2 wt --dir1 mut --biotype 0 =====
    parser.add_argument("--wt-reps", dest="n_opt", metavar="WT_REPS", type=int, help="control group / wild-type replicate count")
    parser.add_argument("--mut-reps", dest="m_opt", metavar="MUT_REPS", type=int, help="experimental group / mutant replicate count")
    parser.add_argument("--dir-wt", dest="dir2_opt", metavar="DIR_WT", help="control group / wild-type sample directory")
    parser.add_argument("--dir-mut", dest="dir1_opt", metavar="DIR_MUT", help="experimental group / mutant sample directory")
    parser.add_argument("--biotype", dest="biotype_opt", metavar="BIOTYPE", type=int, choices=[0, 1, 2], help="0=animal, 1=plant, 2=no filtering")
    parser.add_argument(
        "--meth-diff",
        type=float,
        default=0.0,
        help="Minimum absolute methylation difference required for final DMP filtering, in the range 0-1; for example, 0.25 means 25%%. Default: 0.0 for backward compatibility."
    )
    parser.add_argument(
        "--auto-meth-diff",
        action="store_true",
        help="Enable mianjifa auto-methdiff. After pairwise FDR correction, estimate one global threshold from the raw MethDiff distribution of q-significant sites before MethDiff filtering. The estimated threshold is also used for auto-vote support construction and final common-DMP calling. Disabled by default."
    )
    parser.add_argument(
        "--auto-meth-diff-report-only",
        action="store_true",
        help="Only output mianjifa auto-methdiff diagnostic tables and distribution plots; do not change the actual DMP calling threshold. Disabled by default."
    )
    parser.add_argument(
        "--auto-meth-diff-cut-percent",
        type=float,
        default=0.05,
        help="Histogram area fraction cut outward from zero for mianjifa auto-methdiff. Default: 0.05."
    )
    parser.add_argument(
        "--auto-meth-diff-fallback",
        type=float,
        default=0.3,
        help="Fallback MethDiff threshold used when mianjifa auto-methdiff estimation fails. Default: 0.3."
    )
    parser.add_argument(
        "--auto-meth-diff-aggregate",
        choices=["median", "mean", "max", "min"],
        default="median",
        help="Method for aggregating left/right thresholds across comparisons into one global abs(MethDiff) threshold. Default: median."
    )
    parser.add_argument(
        "--q-cpg",
        type=float,
        default=DMP_QVALUE_THRESHOLDS["CpG"],
        help="CpG DMP q-value threshold, in the range 0-1. Default: use the original code threshold."
    )
    parser.add_argument(
        "--q-chg",
        type=float,
        default=DMP_QVALUE_THRESHOLDS["CHG"],
        help="CHG DMP q-value threshold, in the range 0-1. Default: use the original code threshold."
    )
    parser.add_argument(
        "--q-chh",
        type=float,
        default=DMP_QVALUE_THRESHOLDS["CHH"],
        help="CHH DMP q-value threshold, in the range 0-1. Default: use the original code threshold."
    )
    parser.add_argument(
        "--dmr-q",
        type=float,
        default=DMR_QVALUE_THRESHOLD,
        help="DMR q-value threshold, in the range 0-1. Default: use the original code threshold."
    )
    parser.add_argument(
        "--auto-qvalue-twostep",
        action="store_true",
        help="Estimate the DMP q-value threshold only for two-step-FDR contexts: within the p-value-prefiltered subset, find the maximum qvalue-pvalue point and use its q-value as the pair/context threshold. Disabled by default."
    )
    parser.add_argument(
        "--auto-qvalue-report-only",
        action="store_true",
        help="Only output the two-step auto q-value threshold diagnostic table; do not change the DMP calling threshold. Disabled by default."
    )
    parser.add_argument(
        "--auto-qvalue-p-cutoff",
        type=float,
        default=AUTO_QVALUE_P_CUTOFF,
        help="P-value candidate upper limit used for two-step auto q-value threshold estimation. Default: 0.05."
    )
    parser.add_argument(
        "--auto-qvalue-min-candidates",
        type=int,
        default=AUTO_QVALUE_MIN_CANDIDATES,
        help="Minimum number of candidate sites required for auto q-value threshold estimation; fall back to the fixed q-threshold if the number is insufficient. Default: 10."
    )
    parser.add_argument(
        "--auto-qvalue-use-smooth",
        action="store_true",
        help="Use the smoothed (qvalue-pvalue) curve to find the maximum-difference point. Disabled by default; for formal analysis, keeping the default raw-difference setting is recommended."
    )
    parser.add_argument(
        "--auto-qvalue-smooth-sigma",
        type=float,
        default=AUTO_QVALUE_SMOOTH_SIGMA,
        help="Gaussian sigma used when --auto-qvalue-use-smooth is enabled. Default: 4."
    )
    parser.add_argument(
        "--vote-threshold",
        type=float,
        default=VOTE_THRESHOLD,
        help="Final DMP/DMR voting threshold across replicate combinations, in the range (0,1]. Default: 0.6667, i.e. the original 2/3 rule."
    )
    parser.add_argument(
        "--auto-dmp-vote-threshold",
        action="store_true",
        help="Automatically estimate the integer required_count for final DMP voting. Disabled by default."
    )
    parser.add_argument(
        "--auto-dmr-vote-threshold",
        action="store_true",
        help="Automatically estimate the integer required_count for final DMR voting. Disabled by default."
    )
    parser.add_argument(
        "--auto-vote-threshold-report-only",
        action="store_true",
        help="Only calculate and output auto vote thresholds and distribution plots; do not change final DMP/DMR calling. Disabled by default."
    )
    parser.add_argument(
        "--dmp-lowdiff-strict-vote",
        action="store_true",
        help="Enable final-DMP low-difference strict-vote post-processing: first obtain provisional DMPs using the normal q-value + voting rule; if boundary abs(MethDiff) <= cutoff, require a higher vote count. Disabled by default to preserve legacy behavior."
    )
    parser.add_argument(
        "--dmp-lowdiff-cutoff",
        type=float,
        default=DMP_LOWDIFF_CUTOFF,
        help="Boundary abs(MethDiff) cutoff for final-DMP lowdiff strict voting. Default: 0.3."
    )
    parser.add_argument(
        "--dmp-lowdiff-strict-vote-report-only",
        action="store_true",
        help="Only output the final-DMP lowdiff strict-vote diagnostic table and summary; do not change the final DMP file. Disabled by default."
    )
    parser.add_argument(
        "--skip-dmr",
        action="store_true",
        help="Skip all DMR-related steps and only output DMP/final-DMP results. Disabled by default for backward compatibility."
    )
    parser.add_argument(
        "--skip-window",
        action="store_true",
        help="Skip all sliding-window and visualization plotting steps; only output table results. Disabled by default for backward compatibility."
    )
    parser.add_argument(
        "--dmr-engine",
        choices=["python", "cpp"],
        default="python",
        help="DMR candidate-region detection engine: python=original Python implementation; cpp=use dmr_step1 + dmr_step2_dynamic. Default: python, to facilitate compatibility validation."
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=1,
        help="Number of parallel worker processes. Default: 1, which means the full workflow is serial. When threads > 1, all safe stages are parallelized automatically, including newtoboth file conversion, replicate-pair processing, and common DMR aggregation. Recommended: 2-4 for small datasets; adjust carefully according to memory and I/O for large datasets."
    )
    args = parser.parse_args()

    def choose_arg(opt_value, pos_value, name):
        """Support both new --xxx arguments and legacy positional arguments."""
        if opt_value is not None and pos_value is not None and opt_value != pos_value:
            parser.error(f"parameter conflict:--{name}={opt_value} with the legacy positional argument {pos_value} inconsistent")
        if opt_value is not None:
            return opt_value
        if pos_value is not None:
            return pos_value
        parser.error(f"missing required parameter:--{name}")

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
            print(f"ERROR: {name} must be in 0 to 1 range, for example 0.25 means 25%")
            sys.exit(1)
    if not 0 < args.vote_threshold <= 1:
        print("ERROR: --vote-threshold must be in (0, 1] range, for example 0.667 means 2/3majority voting")
        sys.exit(1)
    if not 0 < args.auto_meth_diff_cut_percent < 1:
        print("ERROR: --auto-meth-diff-cut-percent must be in (0, 1) range, for example 0.05")
        sys.exit(1)
    if args.threads < 1:
        print("ERROR: --threads must be >= 1 integer")
        sys.exit(1)
    if args.auto_qvalue_min_candidates < 1:
        print("ERROR: --auto-qvalue-min-candidates must be >= 1 integer")
        sys.exit(1)
    if args.auto_qvalue_smooth_sigma <= 0:
        print("ERROR: --auto-qvalue-smooth-sigma must be greater than 0")
        sys.exit(1)

    # Write command-line parameters back to global threshold configuration to preserve the existing get_dmp_threshold() call chain
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

    # Validate directory existence and ensure replicate counts are valid
    if not os.path.exists(dir1):
        print(f"ERROR: directory '{dir1}' does not exist.")
        sys.exit(1)
    if not os.path.exists(dir2):
        print(f"ERROR: directory '{dir2}' does not exist.")
        sys.exit(1)
    if m <= 0 or n <= 0:
        print("ERROR:the number of groups must be greater than 0")
        sys.exit(1)

    print(f"\nParameter confirmation:")
    print(f"mutant directory: {dir1} (containing {m} replicate files)")
    print(f"wild-type directory: {dir2} (containing {n} replicate files)")
    print(f"DMP methylation-difference threshold: {meth_diff_threshold} ({meth_diff_threshold * 100:.1f}%)")
    print(
        "mianjifa auto-methdiff: "
        f"enabled for calling={args.auto_meth_diff}, "
        f"report-only={args.auto_meth_diff_report_only}, "
        f"cut_percent={args.auto_meth_diff_cut_percent}, "
        f"aggregate={args.auto_meth_diff_aggregate}, "
        f"fallback={args.auto_meth_diff_fallback}"
    )
    print(
        "DMP q-value threshold: "
        f"CpG={DMP_QVALUE_THRESHOLDS['CpG']}, "
        f"CHG={DMP_QVALUE_THRESHOLDS['CHG']}, "
        f"CHH={DMP_QVALUE_THRESHOLDS['CHH']}"
    )
    print(f"DMR q-value threshold: {DMR_QVALUE_THRESHOLD}")
    print(
        "two-step auto q-value threshold: "
        f"enabled for calling={AUTO_QVALUE_TWOSTEP}, "
        f"report-only={AUTO_QVALUE_REPORT_ONLY}, "
        f"p_cutoff={AUTO_QVALUE_P_CUTOFF}, "
        f"min_candidates={AUTO_QVALUE_MIN_CANDIDATES}, "
        f"use_smooth={AUTO_QVALUE_USE_SMOOTH}"
    )
    print(f"final DMP/DMR vote threshold: {VOTE_THRESHOLD} ({VOTE_THRESHOLD * 100:.1f}%)")
    print(
        "auto vote threshold: "
        f"DMP={AUTO_DMP_VOTE_THRESHOLD}, "
        f"DMR={AUTO_DMR_VOTE_THRESHOLD}, "
        f"report_only={AUTO_VOTE_THRESHOLD_REPORT_ONLY}"
    )
    print(
        "DMP lowdiff strict-vote post-processing: "
        f"enabled for calling={DMP_LOWDIFF_STRICT_VOTE}, "
        f"report-only={DMP_LOWDIFF_STRICT_VOTE_REPORT_ONLY}, "
        f"cutoff={DMP_LOWDIFF_CUTOFF}"
    )
    if (DMP_LOWDIFF_STRICT_VOTE or DMP_LOWDIFF_STRICT_VOTE_REPORT_ONLY) and meth_diff_threshold > 0:
        print(
            "Note: --dmp-lowdiff-strict-vote is final DMP post-processing;"
            f"current --meth-diff={meth_diff_threshold:.6g} will first filter at the pair-support layer,"
            "lowdiff strict vote will then be applied on top of those results."
        )
    print(f"skip DMR analysis: {args.skip_dmr}")
    print(f"skip sliding-window/plotting: {args.skip_window}")
    print(f"parallel worker process count: {args.threads}")

    print("\nStage 1: running newtoboth")
    # Convert Bismark new-format data to both format
    chr_series = newtoboth(m, n, dir1, dir2, threads=args.threads, work_dir=".")
    if biotype == 0:
        unfilter_mtypes = ["CHH", "CHG"]
    elif biotype == 1:
        unfilter_mtypes = ["CpG"]
    elif biotype == 2:
        unfilter_mtypes = ["CHH", "CHG", "CpG"]
    else:
        print("ERROR: biotype must be 0, 1, or 2")
        sys.exit(1)
    print(f"methylation contexts not requiring p-value prefiltering: {unfilter_mtypes}")
    success = process_all_combinations(
        dir1, dir2, m, n, unfilter_mtypes,
        meth_diff_threshold=meth_diff_threshold,
        skip_dmr=args.skip_dmr,
        skip_window=args.skip_window,
        threads=args.threads
    )  # process_all_combinations performs m*n*3 tests

    if success:  # All tests succeeded
        print("\nAll tests and FDR corrections completed successfully!")

        if AUTO_QVALUE_TWOSTEP or AUTO_QVALUE_REPORT_ONLY:
            plot_all_auto_qvalue_panels(
                work_dir=".",
                m=m,
                n=n,
                unfilter_mtypes=unfilter_mtypes
            )

        if args.auto_meth_diff or args.auto_meth_diff_report_only:
            print(
                "\nNote:mianjifa auto-methdiff reads q-significant sites from the complete pairwise FDR tables,"
                "does not use previouslyMethDiff-filtered DMP files; the estimated threshold will continue to be used for "
                "auto-vote support construction and final common-DMP calling."
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
                    "Enable mianjifa auto-methdiff: auto-vote and final common DMP "
                    f"uniformly using methdiff threshold = {meth_diff_threshold:.6g}"
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
                    f"mianjifa auto-methdiff report-only:estimated threshold = {auto_methdiff_threshold:.6g};"
                    f"actually still using --meth-diff = {meth_diff_threshold:.6g}"
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
                print(f"Found {len(common_sites_df)}  {mtype} common significant sites")

                # Perform sliding-window analysis
                if args.skip_window:
                    print(f"skipping {mtype} common significant site sliding-window analysis (--skip-window)")
                else:
                    print(f"Starting {mtype} common significant sites for sliding-window analysis...")
                    results = process_common_sites_sliding_window(
                        common_sites_df=common_sites_df,
                        methytype=mtype
                    )
                    print(f"Completed {mtype} context sliding-window analysis")
            else:
                print(f"No {mtype} common significant sites")

        if args.skip_dmr:
            print("Skipping final common DMR analysis workflow (--skip-dmr)")
        else:
            print("Starting the DMR analysis workflow")
            process_common_sites_dmr_and_summarize(
                dir1=dir1,
                dir2=dir2,
                m=m,
                n=n,
                methylation_types=methylation_types,
                threads=args.threads,
            )

        # Generate visualizations for all sliding-window results
        if args.skip_window:
            print("skip all sliding-window visualization plotting (--skip-window)")
        else:
            plot_methylation_sliding_windows(chr_series=chr_series)
            for mtype111 in ["CpG", "CHH", "CHG"]:
                plot_common_sites_sliding_windows(mtype111, chr_series=chr_series)
        convert_chromosome_to_names(chr_series=chr_series, work_dir=".")
        rename_chromosome_files(chr_series=chr_series, work_dir=".")
        convert_output_to_csv(work_dir=".")
        print(f"- DMP results: output_x_y/<methylation_context>/")
        print(f"- common significant sites: and_output/")
        print(f"- final significant DMR: and_output/*-final_significant_regions_DMRs.txt")
    else:
        print("\nSome tests failed,please check the output messages.")
    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"total elapsed time: {elapsed_time:.2f} seconds")


if __name__ == "__main__":
    main()



