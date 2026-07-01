# MultiDMPcaller

MultiDMPcaller is a one-stop DNA methylation analysis tool for detecting and visualizing **differentially methylated positions (DMPs)** and **differentially methylated regions (DMRs)** between two biological groups, such as a control/wild-type group and an experimental/mutant group. It supports the three major cytosine methylation contexts, **CpG**, **CHG**, and **CHH**, and performs replicate-aware analysis through all pairwise comparisons between the two groups.

---

## 1. Installation

### 1.1 Recommended environment

Python 3.10 or later is recommended.

```bash
conda create -n multidmpcaller python=3.10 -y
conda activate multidmpcaller
pip install -r requirements.txt
```

Alternatively, with `venv`:

```bash
python -m venv multidmpcaller_env
source multidmpcaller_env/bin/activate
pip install -r requirements.txt
```

### 1.2 Optional C++ DMR engine

If you want to use the accelerated C++ DMR engine, make sure the two executable files are placed in the same directory as the main Python script:

```text
dmr_step1
dmr_step2_dynamic
```

Then make them executable:

```bash
chmod +x dmr_step1 dmr_step2_dynamic
```

The Python DMR engine remains available and is the default for compatibility.

---

## Documentation

- [Input format guide](docs/Input_format.md)
- [Parameter guide](docs/Parameter_guide.md)
- [Troubleshooting guide](docs/Troubleshooting.md)
- [Frequently asked questions (FAQ)](docs/FAQ.md)

---

## 2. Input data format

### 2.1 Directory structure

Input files should be placed in two directories: one for the control/wild-type group and one for the experimental/mutant group.

Example for a 2 vs 2 comparison:

```text
MultiDMPcaller/
├── MultiDMPcaller.py
├── wt/
│   ├── 1-wt.txt
│   └── 2-wt.txt
└── mut/
    ├── 1-mut.txt
    └── 2-mut.txt
```

The directory name and file suffix must be consistent. For example, if the wild-type directory is named `wt`, the files should be named `1-wt.txt`, `2-wt.txt`, etc. If the mutant directory is named `mut`, the files should be named `1-mut.txt`, `2-mut.txt`, etc.

### 2.2 File naming rules

For each group directory, sample files must follow this pattern:

```text
{replicate_index}-{directory_name}.txt
```

Examples:

```text
1-wt.txt
2-wt.txt
1-mut.txt
2-mut.txt
```

Replicate indices should start from 1 and be consecutive.

### 2.3 File content format

Each input file should be a plain text file without a header. Each line represents one methylation site and contains five columns separated by spaces or tabs:

```text
Chromosome    Position    Methylated_reads    Unmethylated_reads    Context
```

Example:

```text
chr1    1005    23    5     CpG
chr1    1030    18    9     CHH
chr2    5002    7     32    CHG
```

Notes:

- `Context` should be one of `CpG`, `CHG`, or `CHH`. The program also normalizes common CpG/CG-style labels where applicable.
- Chromosome labels such as `1`, `chr1`, and `Chr1` are automatically normalized internally.
- You do not need to split files by methylation context. MultiDMPcaller automatically processes CpG, CHG, and CHH contexts separately.

---

## 3. Basic command-line usage

The recommended interface uses named arguments:

```bash
python MultiDMPcaller.py \
  --wt-reps 2 \
  --mut-reps 2 \
  --dir-wt wt \
  --dir-mut mut \
  --biotype 1
```

Required arguments:

| Argument | Meaning |
| :--- | :--- |
| `--wt-reps` | Number of control/wild-type replicates. |
| `--mut-reps` | Number of experimental/mutant replicates. |
| `--dir-wt` | Directory containing the control/wild-type sample files. |
| `--dir-mut` | Directory containing the experimental/mutant sample files. |
| `--biotype` | Organism/data mode. `0`: animal; `1`: plant; `2`: no p-value prefiltering for all contexts. |

For clarity, the current README documents the readable argument names `--dir-wt` and `--dir-mut`. Legacy names such as `--dir1` and `--dir2` are not recommended for public use.

---

## 4. Example runs

### 4.1 Quick test

For a small 1 vs 1 test dataset:

```bash
python MultiDMPcaller.py \
  --wt-reps 1 \
  --mut-reps 1 \
  --dir-wt test_data/wt \
  --dir-mut test_data/mut \
  --biotype 1
```

### 4.2 Standard plant WGBS analysis

For a 2 vs 2 plant dataset:

```bash
python MultiDMPcaller.py \
  --wt-reps 2 \
  --mut-reps 2 \
  --dir-wt wt \
  --dir-mut mut \
  --biotype 1 \
  --threads 4 \
  --dmr-engine cpp \
  --auto-qvalue-twostep \
  --auto-dmp-vote-threshold \
  --auto-dmr-vote-threshold \
  --dmp-lowdiff-strict-vote
```

### 4.3 Standard animal WGBS analysis

For a 2 vs 2 animal dataset:

```bash
python MultiDMPcaller.py \
  --wt-reps 2 \
  --mut-reps 2 \
  --dir-wt liver \
  --dir-mut brain \
  --biotype 0 \
  --threads 4 \
  --dmr-engine cpp \
  --auto-qvalue-twostep \
  --auto-dmp-vote-threshold \
  --auto-dmr-vote-threshold \
  --dmp-lowdiff-strict-vote
```

### 4.4 DMP-only analysis

To skip DMR calling and only generate DMP outputs:

```bash
python MultiDMPcaller.py \
  --wt-reps 2 \
  --mut-reps 2 \
  --dir-wt wt \
  --dir-mut mut \
  --biotype 1 \
  --skip-dmr
```

### 4.5 Table-only analysis without visualization

To skip sliding-window visualization and only generate table outputs:

```bash
python MultiDMPcaller.py \
  --wt-reps 2 \
  --mut-reps 2 \
  --dir-wt wt \
  --dir-mut mut \
  --biotype 1 \
  --skip-window
```

---

## 5. Main parameters

### 5.1 Significance thresholds

| Argument | Default | Meaning |
| :--- | :--- | :--- |
| `--q-cpg` | `0.05` | DMP q-value threshold for CpG sites. |
| `--q-chg` | `0.04` | DMP q-value threshold for CHG sites. |
| `--q-chh` | `0.045` | DMP q-value threshold for CHH sites. |
| `--dmr-q` | `0.05` | DMR q-value threshold. |
| `--meth-diff` | `0.0` | Optional hard absolute methylation-difference filter for final DMP calling. For example, `0.2` means 20 percentage points. |

### 5.2 Adaptive q-value and voting modules

| Argument | Meaning |
| :--- | :--- |
| `--auto-qvalue-twostep` | Enable adaptive q-value threshold estimation for two-step FDR contexts. |
| `--auto-dmp-vote-threshold` | Automatically estimate the final DMP voting requirement across replicate comparisons. |
| `--auto-dmr-vote-threshold` | Automatically estimate the final DMR voting requirement across replicate comparisons. |

### 5.3 Low-difference strict voting

| Argument | Meaning |
| :--- | :--- |
| `--dmp-lowdiff-strict-vote` | Enable stricter final-DMP voting for low-difference candidates. |
| `--dmp-lowdiff-cutoff` | Boundary absolute MethDiff cutoff used to define low-difference final DMP candidates. Default: `0.3`. |

### 5.4 Runtime and workflow control

| Argument | Default | Meaning |
| :--- | :--- | :--- |
| `--threads` | `1` | Number of parallel worker processes. Recommended values are usually 2–4 for ordinary workstations. |
| `--dmr-engine` | `python` | DMR candidate-region engine. Use `cpp` for the accelerated C++ engine if `dmr_step1` and `dmr_step2_dynamic` are available. |
| `--skip-dmr` | disabled | Skip all DMR-related steps. |
| `--skip-window` | disabled | Skip sliding-window and visualization steps. |

---

## 6. Output files

MultiDMPcaller writes intermediate and final results to the working directory. The main final outputs are saved in the `and_output/` directory.

For each methylation context (`CpG`, `CHG`, and `CHH`), the key outputs include:

```text
and_output/{context}-final_significant_sites_DMPs.txt
and_output/{context}-final_significant_sites_DMPs.csv
and_output/{context}-final_significant_regions_DMRs.txt
and_output/{context}-final_significant_regions_DMRs.csv
```

Examples:

```text
and_output/CpG-final_significant_sites_DMPs.txt
and_output/CHG-final_significant_sites_DMPs.txt
and_output/CHH-final_significant_sites_DMPs.txt
and_output/CpG-final_significant_regions_DMRs.txt
```

Visualization files are also generated unless `--skip-window` is used. Typical visualization outputs include context-specific chromosomal DMP/DMR distribution plots.

---

## 7. Output interpretation

### 7.1 Final DMP table

A final DMP table contains significant methylation sites supported across replicate comparisons. Typical columns include:

| Column | Meaning |
| :--- | :--- |
| `Chromosome` | Chromosome identifier. |
| `Methylation_Type` | Methylation context: CpG, CHG, or CHH. |
| `Position` | Genomic coordinate of the DMP. |
| `Methylation_Change` | Direction of methylation change. `1` usually denotes hyper-methylation in the experimental/mutant group; `0` denotes hypo-methylation. |
| `Hyper_Count` / `Hypo_Count` | Number of pairwise comparisons supporting hyper- or hypo-methylation. |
| `Num_Comparisons` | Total number of pairwise comparisons. |
| `Sig_Mean_Qvalue` | Mean q-value across significant pairwise comparisons. |

### 7.2 Final DMR table

A final DMR table contains significant methylated regions inferred from DMP-enriched candidate regions and replicate voting. Typical columns include:

| Column | Meaning |
| :--- | :--- |
| `Chromosome` | Chromosome identifier. |
| `Methylation_Type` | Methylation context. |
| `DMR_start` / `DMR_end` | Genomic boundaries of the DMR. |
| `Length` | Region length. |
| `Direction` | Direction of methylation change. |
| `Sig_count` / `Total_count` | Number of supporting comparisons and total comparisons. |
| `Sig_probability` | Support proportion across comparisons. |
| `Sig_Avg_qvalue` | Average q-value across supporting comparisons. |

---

## 8. Web server

A public web server is available at:

```text
https://ciebioinfo.nwafu.edu.cn/
```

The web interface supports file upload, parameter selection, job submission, job-ID based status query, result visualization, and ZIP result download. For large jobs, keep the returned Job ID so that results can be retrieved later.

---

## 9. Notes and recommendations

1. Make sure the number of replicate files matches `--wt-reps` and `--mut-reps`.
2. Make sure the directory names match the file suffixes. For example, files inside `wt/` should be named `1-wt.txt`, `2-wt.txt`, etc.
3. For plant WGBS data, use `--biotype 1`; for animal WGBS data, use `--biotype 0`.
4. For large datasets, use `--threads 2` to `--threads 4` first. Increasing threads may increase memory and I/O pressure.
5. Use `--dmr-engine cpp` only when the C++ executables are available and executable.
6. For formal analyses, it is recommended to keep all log files and the full `and_output/` directory for reproducibility.

---

## 10. Citation

If you use MultiDMPcaller, please cite:

```text
Yuan Q#., Zhao H#., Zhang Z#., Yue C., Zhang B., Zou Q. and Yu J*. 2026. MultiDMPcaller: A one-stop software for detection and visualization of differentially methylated positions and regions. Bioinformatics (under review).
```

