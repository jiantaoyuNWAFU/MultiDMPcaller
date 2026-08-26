# MultiDMPcaller

MultiDMPcaller is a one-stop DNA methylation analysis tool for detecting and visualizing **differentially methylated positions (DMPs)** and **differentially methylated regions (DMRs)** between two biological groups, such as a control/wild-type group and an experimental/mutant group. It supports the three major cytosine methylation contexts, **CpG**, **CHG**, and **CHH**, and performs replicate-aware analysis through all pairwise comparisons between the two groups.

---

## 1. Installation

### 1.1 Recommended environment

Python 3.10 is recommended and was used for validation.

Create and activate a Conda environment:

```bash
conda create -n multidmpcaller python=3.10 -y
conda activate multidmpcaller
pip install -r requirements.txt
```

Alternatively, create a virtual environment with `venv`:

```bash
python -m venv multidmpcaller_env
```

Activate the environment on Linux or macOS:

```bash
source multidmpcaller_env/bin/activate
```

Activate it in Windows PowerShell:

```powershell
.\multidmpcaller_env\Scripts\Activate.ps1
```

Activate it in Windows Command Prompt:

```bat
multidmpcaller_env\Scripts\activate.bat
```

Then install the dependencies:

```bash
pip install -r requirements.txt
```

### 1.2 Optional C++ DMR engine

MultiDMPcaller provides both Python and C++ implementations for DMR candidate-region detection. The Python engine is the default and does not require external executables.

The repository includes the following precompiled executables and C++ source files:

```text
dmr_step1
dmr_step2_dynamic
dmr_step1.cpp
dmr_step2_dynamic.cpp
```

On Linux or macOS, the supplied executables can be made executable with:

```bash
chmod +x dmr_step1 dmr_step2_dynamic
```

The C++ programs can also be compiled from source with a compatible C++ compiler. For example, on Linux or macOS:

```bash
g++ -O3 -std=c++17 -o dmr_step1 dmr_step1.cpp
g++ -O3 -std=c++17 -o dmr_step2_dynamic dmr_step2_dynamic.cpp
chmod +x dmr_step1 dmr_step2_dynamic
```

On Windows with MinGW-w64 or another compatible `g++` installation:

```powershell
g++ -O3 -std=c++17 -o dmr_step1.exe dmr_step1.cpp
g++ -O3 -std=c++17 -o dmr_step2_dynamic.exe dmr_step2_dynamic.cpp
```

On Windows, make sure the directory containing the compiled `.exe` files is included in `PATH` before running MultiDMPcaller.

Enable the accelerated engine with:

```bash
--dmr-engine cpp
```

The Python engine remains available through:

```bash
--dmr-engine python
```

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

Example for a 2 × 2 comparison:

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

The directory name and file suffix must be consistent. For example, if the wild-type directory is named `wt`, the files should be named `1-wt.txt`, `2-wt.txt`, and so on. If the mutant directory is named `mut`, the files should be named `1-mut.txt`, `2-mut.txt`, and so on.

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

Each input file should be a plain-text file without a header. Each line represents one methylation site and contains five columns separated by spaces or tabs:

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

- `Context` should be one of `CpG`, `CHG`, or `CHH`. Common `CG`-style labels are normalized to `CpG` where applicable.
- Chromosome labels such as `1`, `chr1`, and `Chr1` are normalized internally.
- Input files do not need to be split by methylation context. MultiDMPcaller processes CpG, CHG, and CHH contexts separately.

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
| `--biotype` | Organism/data mode. `0`: animal; `1`: plant; `2`: no p-value prefiltering for any context. |

Legacy positional invocation is retained internally for backward compatibility. New analyses should use the named arguments documented above.

---

## 4. Example runs

### 4.1 Quick test

For a small 1 × 1 test dataset:

```bash
python MultiDMPcaller.py \
  --wt-reps 1 \
  --mut-reps 1 \
  --dir-wt test_data/wt \
  --dir-mut test_data/mut \
  --biotype 1
```

### 4.2 Example adaptive plant WGBS workflow

The following example explicitly enables the optional adaptive q-value, automatic voting, and low-difference strict-voting modules. These modules are disabled by default.

```bash
python MultiDMPcaller.py \
  --wt-reps 2 \
  --mut-reps 2 \
  --dir-wt wt \
  --dir-mut mut \
  --biotype 1 \
  --q-cpg 0.05 \
  --q-chg 0.04 \
  --q-chh 0.045 \
  --dmr-q 0.05 \
  --methy-diff-dmp 0.0 \
  --methy-diff-dmr 0.0 \
  --vote-threshold 2/3 \
  --processes 4 \
  --dmr-engine cpp \
  --auto-qvalue-twostep \
  --auto-dmp-vote-threshold \
  --auto-dmr-vote-threshold \
  --dmp-lowdiff-strict-vote \
  --dmp-lowdiff-cutoff 0.3
```

### 4.3 Example adaptive animal WGBS workflow

```bash
python MultiDMPcaller.py \
  --wt-reps 2 \
  --mut-reps 2 \
  --dir-wt liver \
  --dir-mut brain \
  --biotype 0 \
  --q-cpg 0.05 \
  --q-chg 0.04 \
  --q-chh 0.045 \
  --dmr-q 0.05 \
  --methy-diff-dmp 0.0 \
  --methy-diff-dmr 0.0 \
  --vote-threshold 2/3 \
  --processes 4 \
  --dmr-engine cpp \
  --auto-qvalue-twostep \
  --auto-dmp-vote-threshold \
  --auto-dmr-vote-threshold \
  --dmp-lowdiff-strict-vote \
  --dmp-lowdiff-cutoff 0.3
```

### 4.4 DMP-only analysis

To skip DMR calling and generate only DMP outputs:

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

To skip sliding-window visualization and generate table outputs only:

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

### 5.1 Significance and methylation-difference thresholds

| Argument | Default | Meaning |
| :--- | :--- | :--- |
| `--q-cpg` | `0.05` | Fixed DMP q-value threshold for CpG sites. |
| `--q-chg` | `0.04` | Fixed DMP q-value threshold for CHG sites. |
| `--q-chh` | `0.045` | Fixed DMP q-value threshold for CHH sites. |
| `--dmr-q` | `0.05` | DMR q-value threshold for each pairwise regional comparison. |
| `--methy-diff-dmp` | `0.0` | Minimum absolute site-level methylation difference required for a pairwise comparison to contribute DMP support. The value must be between `0` and `1`; for example, `0.2` means 20 percentage points. |
| `--methy-diff-dmr` | `0.0` | Minimum absolute regional methylation difference required for each pairwise DMR support. The regional difference is calculated from aggregated methylated and unmethylated read counts within the candidate region. |

The legacy `--meth-diff` and `--dmr-meth-diff` names are retained internally as hidden compatibility aliases. New commands should use `--methy-diff-dmp` and `--methy-diff-dmr`.

### 5.2 Adaptive q-value and voting modules

| Argument | Default | Meaning |
| :--- | :--- | :--- |
| `--vote-threshold` | `0.6667` | Fixed support proportion used for final DMP and DMR voting. Fractions such as `2/3` are accepted. This value is also used as the fallback when automatic voting cannot produce a valid threshold. |
| `--auto-qvalue-twostep` | disabled | Estimate the DMP q-value threshold for contexts using two-step FDR correction. |
| `--auto-dmp-vote-threshold` | disabled | Automatically estimate the integer support count required for final DMP calling. |
| `--auto-dmr-vote-threshold` | disabled | Automatically estimate the integer support count required for final DMR calling. |

If an automatic DMP or DMR voting model raises an exception during fitting or threshold calculation, MultiDMPcaller falls back to the value specified by `--vote-threshold`. If the parameter is not explicitly supplied, the default two-thirds rule is used.

### 5.3 Low-difference strict voting

| Argument | Default | Meaning |
| :--- | :--- | :--- |
| `--dmp-lowdiff-strict-vote` | disabled | Apply an additional strict-voting rule to provisional final DMPs with relatively small boundary methylation differences. |
| `--dmp-lowdiff-cutoff` | `0.3` | Boundary absolute MethDiff cutoff used to identify provisional DMPs that require stricter voting. |

The hard DMP filter specified by `--methy-diff-dmp` is applied first at the pairwise-support level. Low-difference strict voting is then applied to the remaining provisional DMPs after the ordinary q-value and voting procedure.

### 5.4 Runtime and workflow control

| Argument | Default | Meaning |
| :--- | :--- | :--- |
| `--processes` | `1` | Maximum number of parallel worker processes used by safe parallel stages, including raw-file conversion, replicate-pair processing, and common DMR aggregation. The actual number of preprocessing workers cannot exceed the number of replicate files. Start with `2`–`4` for ordinary workstations and increase it only when sufficient memory and disk I/O bandwidth are available. |
| `--dmr-engine` | `python` | DMR candidate-region engine. Use `cpp` after compatible `dmr_step1` and `dmr_step2_dynamic` executables have been prepared. |
| `--skip-dmr` | disabled | Skip all DMR-related steps. |
| `--skip-window` | disabled | Skip sliding-window and visualization steps. |

### 5.5 Low-memory and cross-platform preprocessing

Before statistical testing, MultiDMPcaller converts each raw five-column methylation file into context-specific matrix files.

The current preprocessing implementation:

- reads large input files in chunks;
- performs stable sorting within each chunk;
- writes temporary sorted runs by methylation context and chromosome;
- uses a deterministic pure-Python k-way merge;
- writes final matrices in row blocks;
- does not depend on the GNU/Linux `sort` command;
- uses `--processes` directly for process-based parallel conversion of replicate files.

The implementation is designed to reduce peak memory usage and avoid operating-system-specific sorting commands.

Two optional environment variables control only the I/O block sizes:

| Environment variable | Default | Meaning |
| :--- | :--- | :--- |
| `MULTIDMPCALLER_NEWTOBOTH_CHUNKSIZE` | `1000000` | Number of input rows read per chunk. |
| `MULTIDMPCALLER_NEWTOBOTH_BLOCK_ROWS` | `100000` | Number of matrix rows written per output block. |

Most users do not need to set these variables. Lower values may reduce peak memory use but can increase temporary-file and disk-I/O overhead. Input directories must be writable and should have sufficient free disk space for generated matrices and temporary spill files.

---

## 6. Output files

MultiDMPcaller writes intermediate and final results to the working directory. The main final outputs are saved in the `and_output/` directory.

When results are available for a methylation context, the key final outputs include:

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

Visualization files are generated unless `--skip-window` is used. Typical visualization outputs include context-specific chromosomal DMP/DMR distribution plots.

---

## 7. Output interpretation

### 7.1 Final DMP table

A final DMP table contains significant methylation sites supported across replicate comparisons. Typical columns include:

| Column | Meaning |
| :--- | :--- |
| `Chromosome` | Chromosome identifier. |
| `Methylation_Type` | Methylation context: CpG, CHG, or CHH. |
| `Position` | Genomic coordinate of the DMP. |
| `Methylation_Change` | Direction of methylation change. `1` indicates higher methylation in the experimental/mutant group, whereas `0` indicates lower methylation in the experimental/mutant group. |
| `Hyper_Ratio` | Proportion of supporting pairwise comparisons indicating hypermethylation in the experimental/mutant group. |
| `Hyper_Count` / `Hypo_Count` | Number of pairwise comparisons supporting hyper- or hypomethylation. |
| `Num_Comparisons` | Total number of pairwise comparisons. |
| `Sig_Mean_Qvalue` | Mean q-value across significant pairwise comparisons. |
| `qvalue_*` | Pairwise q-values for the corresponding wild-type and mutant replicate comparison. |

### 7.2 Final DMR table

A final DMR table contains significant methylated regions inferred from DMP-enriched candidate regions and replicate voting. Typical columns include:

| Column | Meaning |
| :--- | :--- |
| `Chromosome` | Chromosome identifier. |
| `Methylation_Type` | Methylation context. |
| `DMR_start` / `DMR_end` | Genomic boundaries of the DMR. |
| `Length` | Region length. |
| `Direction` | Regional methylation direction. `1` indicates hypermethylation and `0` indicates hypomethylation in the experimental/mutant group. |
| `Sig_count` / `Total_count` | Number of supporting comparisons and total comparisons. |
| `Sig_probability` | Support proportion across comparisons. |
| `Sig_Avg_qvalue` | Average q-value across supporting comparisons. |

---

## 8. Web server

A public web server is available at:

[https://ciebioinfo.nwafu.edu.cn/](https://ciebioinfo.nwafu.edu.cn/)

The web interface supports file upload, parameter selection, job submission, job-ID-based status queries, result visualization, and ZIP result download. Keep the returned Job ID so that results can be retrieved later.

Current service limits:

| Item | Current setting |
| :--- | :--- |
| Maximum worker processes per job | `4` |
| Maximum concurrent jobs | `4` |
| Upload limit | No fixed software-level cap |
| Result retention period | `72 h` |

“No fixed software-level cap” does not imply unlimited upload size. Practical limits may still depend on browser behavior, network conditions, reverse-proxy configuration, available storage, and server resources.

---

## 9. Code and data availability

The source code, documentation, Python implementation, C++ DMR source files, and release history of MultiDMPcaller are available in this GitHub repository.

The human simulation benchmark datasets used in the study are available from Zenodo:

[Human simulation benchmark datasets for MultiDMPcaller](https://zenodo.org/records/21121883)

Publicly available experimental datasets used in the study are identified by their original database accession numbers in the manuscript and associated benchmark documentation.

---

## 10. Notes and recommendations

1. Make sure the number of replicate files matches `--wt-reps` and `--mut-reps`.
2. Make sure directory names match the file suffixes. For example, files inside `wt/` should be named `1-wt.txt`, `2-wt.txt`, and so on.
3. Use `--biotype 1` for plant WGBS data and `--biotype 0` for animal WGBS data.
4. For large datasets, start with `--processes 2` to `--processes 4`. Increasing the number of worker processes can increase memory use and disk-I/O pressure.
5. Use `--dmr-engine cpp` only after compatible C++ executables have been prepared and made discoverable by the program.
6. Keep log files, command lines, software versions, and the complete `and_output/` directory for reproducibility.
7. Ensure sufficient free disk space is available before processing whole-genome bisulfite sequencing datasets.

---

## 11. Citation

If you use MultiDMPcaller, please cite:

```text
Yuan Q#., Zhao H#., Zhang Z#., Yue C., Zhang B., Xue S., Zou Q. and Yu J*. 2026. MultiDMPcaller: A one-stop software for detection and visualization of differentially methylated positions and regions. Bioinformatics (accepted).
```
