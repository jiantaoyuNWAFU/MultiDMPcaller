# Parameter Guide

This page explains the main command-line parameters of **MultiDMPcaller** and provides practical recommendations for common workflows.

## 1. Required parameters

### `--wt-reps`

Number of control/wild-type biological replicates.

```bash
--wt-reps 2
```

### `--mut-reps`

Number of experimental/mutant biological replicates.

```bash
--mut-reps 2
```

### `--dir-wt`

Directory containing the control/wild-type files.

```bash
--dir-wt wt
```

### `--dir-mut`

Directory containing the experimental/mutant files.

```bash
--dir-mut mut
```

### `--biotype`

Organism/data mode:

```text
0 = animal
1 = plant
2 = no p-value prefiltering for any context
```

| Data type | Recommended setting |
| :--- | :--- |
| Animal WGBS data | `--biotype 0` |
| Plant WGBS data | `--biotype 1` |
| Generic mode without p-value prefiltering | `--biotype 2` |

## 2. Basic examples

### 2.1 Plant WGBS

```bash
python MultiDMPcaller.py \
  --wt-reps 2 \
  --mut-reps 2 \
  --dir-wt wt \
  --dir-mut mut \
  --biotype 1
```

### 2.2 Animal WGBS

```bash
python MultiDMPcaller.py \
  --wt-reps 2 \
  --mut-reps 2 \
  --dir-wt liver \
  --dir-mut brain \
  --biotype 0
```

### 2.3 DMP-only analysis

```bash
python MultiDMPcaller.py \
  --wt-reps 2 \
  --mut-reps 2 \
  --dir-wt wt \
  --dir-mut mut \
  --biotype 1 \
  --skip-dmr
```

### 2.4 Table-only analysis

```bash
python MultiDMPcaller.py \
  --wt-reps 2 \
  --mut-reps 2 \
  --dir-wt wt \
  --dir-mut mut \
  --biotype 1 \
  --skip-window
```

## 3. Significance thresholds

### `--q-cpg`

Fixed DMP q-value threshold for CpG sites.

```bash
--q-cpg 0.05
```

### `--q-chg`

Fixed DMP q-value threshold for CHG sites.

```bash
--q-chg 0.04
```

### `--q-chh`

Fixed DMP q-value threshold for CHH sites.

```bash
--q-chh 0.045
```

### `--dmr-q`

Q-value threshold for each pairwise DMR comparison.

```bash
--dmr-q 0.05
```

Values must be between `0` and `1`. More stringent values generally reduce the number of retained candidates.

## 4. Methylation-difference thresholds

### `--methy-diff-dmp`

Minimum absolute site-level methylation difference required for a pairwise comparison to contribute DMP support.

Default:

```bash
--methy-diff-dmp 0.0
```

Example:

```bash
--methy-diff-dmp 0.2
```

A value of `0.2` means 20 percentage points.

This filter is applied at the pairwise-support layer before final DMP voting.

### `--methy-diff-dmr`

Minimum absolute regional methylation difference required for each pairwise DMR support.

Default:

```bash
--methy-diff-dmr 0.0
```

Example:

```bash
--methy-diff-dmr 0.2
```

The regional difference is calculated from aggregated methylated and unmethylated read counts within the candidate region.

The legacy names `--meth-diff` and `--dmr-meth-diff` are retained only as hidden compatibility aliases. New commands should use `--methy-diff-dmp` and `--methy-diff-dmr`.

## 5. Voting threshold

### `--vote-threshold`

Fixed support proportion used for final DMP and DMR voting.

Default:

```bash
--vote-threshold 2/3
```

Decimal values such as `0.6666666667` are also accepted.

For `m` control replicates and `n` experimental replicates, MultiDMPcaller performs `m × n` pairwise comparisons and converts the support proportion into an integer required count.

This threshold is used when automatic voting is disabled. It also serves as the fallback if automatic DMP or DMR voting raises an exception or cannot produce a valid threshold.

## 6. Adaptive q-value module

### `--auto-qvalue-twostep`

Estimate a data-adaptive DMP q-value threshold for contexts using two-step FDR correction.

```bash
--auto-qvalue-twostep
```


### Additional controls

| Argument | Default | Meaning |
| :--- | :--- | :--- |
| `--auto-qvalue-p-cutoff` | `0.05` | Candidate p-value upper limit used during estimation. |
| `--auto-qvalue-min-candidates` | `10` | Minimum candidate count required before automatic estimation is attempted. |
| `--auto-qvalue-use-smooth` | disabled | Use a smoothed q-value-minus-p-value curve. |
| `--auto-qvalue-smooth-sigma` | `4` | Gaussian smoothing sigma when smoothing is enabled. |

For formal analysis, the default raw-difference setting is recommended unless smoothing has been explicitly justified.

## 7. Automatic voting modules

### `--auto-dmp-vote-threshold`

Automatically estimate the integer support count required for final DMP calling.

### `--auto-dmr-vote-threshold`

Automatically estimate the integer support count required for final DMR calling.


If GMM fitting or threshold calculation raises an exception, the program falls back to the support requirement implied by `--vote-threshold`. A convergence warning alone does not necessarily trigger fallback.

## 8. Low-difference strict voting

### `--dmp-lowdiff-strict-vote`

Apply an additional strict-voting rule to provisional final DMPs with relatively small boundary methylation differences.

### `--dmp-lowdiff-cutoff`

Boundary absolute MethDiff cutoff used to identify provisional DMPs requiring stricter voting.

Default:

```bash
--dmp-lowdiff-cutoff 0.3
```


The `--methy-diff-dmp` hard filter is applied first at the pairwise-support layer. Low-difference strict voting is applied later to provisional final DMPs.

## 9. Runtime and workflow control

### `--processes`

Maximum number of parallel worker processes.

Default:

```bash
--processes 1
```

Safe parallel stages include raw-file conversion, replicate-pair processing, and common DMR aggregation. The actual number of active worker processes may be lower when fewer independent tasks are available.

Start with:

```bash
--processes 2
```

or:

```bash
--processes 4
```

Increasing the number of worker processes can increase memory use and disk-I/O pressure.

### `--dmr-engine`

DMR candidate-region engine.

Default:

```bash
--dmr-engine python
```

Accelerated C++ mode:

```bash
--dmr-engine cpp
```

The repository includes:

```text
dmr_step1
dmr_step2_dynamic
dmr_step1.cpp
dmr_step2_dynamic.cpp
```

On Linux or macOS:

```bash
g++ -O3 -std=c++17 -o dmr_step1 dmr_step1.cpp
g++ -O3 -std=c++17 -o dmr_step2_dynamic dmr_step2_dynamic.cpp
chmod +x dmr_step1 dmr_step2_dynamic
```

On Windows with MinGW-w64 or a compatible `g++`:

```powershell
g++ -O3 -std=c++17 -o dmr_step1.exe dmr_step1.cpp
g++ -O3 -std=c++17 -o dmr_step2_dynamic.exe dmr_step2_dynamic.cpp
```

On Windows, add the directory containing the compiled `.exe` files to `PATH`.

### `--skip-dmr`

Skip all DMR-related steps and generate DMP outputs only.

### `--skip-window`

Skip sliding-window and visualization steps and generate table outputs only.

## 10. Low-memory preprocessing

The raw-file conversion stage:

- reads input files in chunks;
- writes output matrices in row blocks;
- uses `--processes` for process-based parallel conversion of replicate files.

Optional environment variables:

| Variable | Default | Meaning |
| :--- | :--- | :--- |
| `MULTIDMPCALLER_NEWTOBOTH_CHUNKSIZE` | `1000000` | Number of input rows read per chunk. |
| `MULTIDMPCALLER_NEWTOBOTH_BLOCK_ROWS` | `100000` | Number of matrix rows written per block. |

Most users should keep the defaults. Lower values can reduce peak memory but increase temporary-file and disk-I/O overhead.

Input directories must be writable and have enough free disk space for generated matrices and temporary spill files.

## 11. Suggested parameter sets

### 11.1 Quick sanity test

```bash
python MultiDMPcaller.py \
  --wt-reps 1 \
  --mut-reps 1 \
  --dir-wt test_data/wt \
  --dir-mut test_data/mut \
  --biotype 1 \
  --skip-dmr \
  --skip-window
```

### 11.2 Adaptive plant workflow

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
  --auto-dmr-vote-threshold
```

### 11.3 Adaptive animal workflow

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
  --auto-dmr-vote-threshold
```

### 11.4 Conservative high-confidence workflow

```bash
python MultiDMPcaller.py \
  --wt-reps 2 \
  --mut-reps 2 \
  --dir-wt wt \
  --dir-mut mut \
  --biotype 1 \
  --methy-diff-dmp 0.2 \
  --methy-diff-dmr 0.2 \
  --vote-threshold 0.75 \
  --processes 4
```

### 11.5 Exploratory workflow

```bash
python MultiDMPcaller.py \
  --wt-reps 2 \
  --mut-reps 2 \
  --dir-wt wt \
  --dir-mut mut \
  --biotype 1 \
  --methy-diff-dmp 0 \
  --methy-diff-dmr 0 \
  --vote-threshold 0.5 \
  --processes 4
```

## 12. Reproducibility checklist

Record:

- MultiDMPcaller version or GitHub commit hash;
- full command line;
- `--biotype`;
- all q-value thresholds;
- `--methy-diff-dmp` and `--methy-diff-dmr`;
- `--vote-threshold`;
- whether adaptive q-value estimation was enabled;
- whether automatic DMP or DMR voting was enabled;
- whether low-difference strict voting was enabled;
- whether the Python or C++ DMR engine was used;
- number of worker processes;
- operating system, architecture, and Python version;
- checksums of the input files when possible.
