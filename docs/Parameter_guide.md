# Parameter Guide

This page explains the main command-line parameters of **MultiDMPcaller** and provides practical recommendations for common use cases.

## 1. Required parameters

### `--wt-reps`

Number of control/wild-type replicates.

Example:

```bash
--wt-reps 2
```

This means the control/wild-type directory should contain:

```text
1-wt.txt
2-wt.txt
```

if the directory is named `wt`.

### `--mut-reps`

Number of experimental/mutant replicates.

Example:

```bash
--mut-reps 2
```

This means the experimental/mutant directory should contain:

```text
1-mut.txt
2-mut.txt
```

if the directory is named `mut`.

### `--dir-wt`

Directory containing the control/wild-type replicate files.

Example:

```bash
--dir-wt wt
```

### `--dir-mut`

Directory containing the experimental/mutant replicate files.

Example:

```bash
--dir-mut mut
```

### `--biotype`

Organism/data mode.

```text
0 = animal
1 = plant
2 = no p-value prefiltering for all contexts
```

Recommended use:

| Data type | Recommended setting |
|---|---|
| Animal WGBS data | `--biotype 0` |
| Plant WGBS data | `--biotype 1` |
| Generic/no p-value prefiltering mode | `--biotype 2` |

## 2. Basic examples

### 2.1 Plant WGBS example

```bash
python MultiDMPcaller.py \
  --wt-reps 2 \
  --mut-reps 2 \
  --dir-wt wt \
  --dir-mut mut \
  --biotype 1
```

### 2.2 Animal WGBS example

```bash
python MultiDMPcaller.py \
  --wt-reps 2 \
  --mut-reps 2 \
  --dir-wt liver \
  --dir-mut brain \
  --biotype 0
```

### 2.3 DMP-only example

```bash
python MultiDMPcaller.py \
  --wt-reps 2 \
  --mut-reps 2 \
  --dir-wt wt \
  --dir-mut mut \
  --biotype 1 \
  --skip-dmr
```

### 2.4 Table-only example without visualization

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

DMP q-value threshold for CpG sites.

Default:

```bash
--q-cpg 0.05
```

### `--q-chg`

DMP q-value threshold for CHG sites.

Default:

```bash
--q-chg 0.04
```

### `--q-chh`

DMP q-value threshold for CHH sites.

Default:

```bash
--q-chh 0.045
```

### `--dmr-q`

DMR q-value threshold.

Default:

```bash
--dmr-q 0.05
```

### Practical notes

More stringent q-value thresholds reduce the number of detected DMPs/DMRs and may increase precision, but they can also reduce recall.

More permissive q-value thresholds may detect more candidates, but they can also increase false positives.

## 4. Methylation-difference threshold

### `--meth-diff`

Optional hard absolute methylation-difference filter for final DMP calling.

Default:

```bash
--meth-diff 0.0
```

Example:

```bash
--meth-diff 0.2
```

This means that a candidate DMP must have an absolute methylation-level difference of at least 0.2, corresponding to 20 percentage points.

### Practical notes

- Use a lower value for exploratory analysis or weak-effect datasets.
- Use a higher value for high-confidence candidate screening.
- Common exploratory values include `0`, `0.1`, and `0.2`.
- More stringent values such as `0.3` may be useful when focusing on stronger methylation differences.

## 5. Voting threshold

### `--vote-threshold`

Final DMP/DMR voting threshold across replicate pairwise comparisons.

For `m` control/wild-type replicates and `n` experimental/mutant replicates, MultiDMPcaller performs `m × n` pairwise comparisons.

The final result is selected based on the support proportion across these pairwise comparisons.

Default:

```bash
--vote-threshold 0.6666666666666666
```

This corresponds approximately to a two-thirds majority rule.

### Example

For a 2 vs 2 design, there are 4 pairwise comparisons.

A two-thirds threshold means that a DMP/DMR must be supported by enough pairwise comparisons to pass the required support count derived from the threshold.

### Practical notes

- Lower thresholds, such as `0.5`, are more permissive and may increase recall.
- Higher thresholds, such as `0.75` or `1.0`, are more stringent and may increase confidence but reduce recall.
- The default two-thirds threshold is intended as a balanced setting.

## 6. Adaptive q-value module

### `--auto-qvalue-twostep`

Enable adaptive q-value threshold estimation for two-step FDR contexts.

Example:

```bash
--auto-qvalue-twostep
```

### `--auto-qvalue-report-only`

Generate diagnostic tables for adaptive q-value estimation without changing final calling thresholds.

Example:

```bash
--auto-qvalue-report-only
```

### Recommended use

If you are testing the effect of adaptive q-value estimation for the first time, use report-only mode first:

```bash
python MultiDMPcaller.py \
  --wt-reps 2 \
  --mut-reps 2 \
  --dir-wt wt \
  --dir-mut mut \
  --biotype 1 \
  --auto-qvalue-twostep \
  --auto-qvalue-report-only
```

Then inspect the diagnostic output before deciding whether to apply adaptive thresholds in formal analysis.

## 7. Automatic voting-threshold modules

### `--auto-dmp-vote-threshold`

Automatically estimate the final DMP voting requirement across replicate comparisons.

### `--auto-dmr-vote-threshold`

Automatically estimate the final DMR voting requirement across replicate comparisons.

### `--auto-vote-threshold-report-only`

Report auto-estimated voting thresholds without applying them to final DMP/DMR calling.

### Recommended use

For diagnostic runs:

```bash
--auto-dmp-vote-threshold \
--auto-dmr-vote-threshold \
--auto-vote-threshold-report-only
```

This allows users to inspect the estimated voting thresholds before changing the final calling behavior.

## 8. Low-difference strict voting

### `--dmp-lowdiff-strict-vote`

Enable stricter final-DMP voting for low-difference candidates.

### `--dmp-lowdiff-cutoff`

Boundary absolute MethDiff cutoff used to define low-difference final DMP candidates.

Default:

```bash
--dmp-lowdiff-cutoff 0.3
```

### `--dmp-lowdiff-strict-vote-report-only`

Generate diagnostic output for low-difference strict voting without changing the final DMP file.

### Recommended use

For first-time testing, use:

```bash
--dmp-lowdiff-strict-vote-report-only
```

This provides diagnostic output while preserving the original final DMP file.

## 9. Runtime and workflow control

### `--threads`

Number of parallel worker processes.

Default:

```bash
--threads 1
```

Recommended starting values:

```bash
--threads 2
--threads 4
```

Using too many threads may increase memory and I/O pressure, especially on shared servers.

### `--dmr-engine`

DMR candidate-region engine.

Default:

```bash
--dmr-engine python
```

To use the accelerated C++ DMR engine:

```bash
--dmr-engine cpp
```

Before using the C++ engine, make sure the following files are available and executable:

```text
dmr_step1
dmr_step2_dynamic
```

Run:

```bash
chmod +x dmr_step1 dmr_step2_dynamic
```

### `--skip-dmr`

Skip all DMR-related steps.

This is useful for DMP-only analysis or quick testing.

Example:

```bash
--skip-dmr
```

### `--skip-window`

Skip sliding-window and visualization steps.

This is useful when only tabular output is needed.

Example:

```bash
--skip-window
```

## 10. Suggested parameter sets

### 10.1 Quick sanity test

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

### 10.2 Standard plant analysis

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
  --dmp-lowdiff-strict-vote-report-only
```

### 10.3 Standard animal analysis

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
  --dmp-lowdiff-strict-vote-report-only
```

### 10.4 Conservative high-confidence analysis

```bash
python MultiDMPcaller.py \
  --wt-reps 2 \
  --mut-reps 2 \
  --dir-wt wt \
  --dir-mut mut \
  --biotype 1 \
  --meth-diff 0.2 \
  --vote-threshold 0.75 \
  --threads 4
```

### 10.5 Exploratory analysis

```bash
python MultiDMPcaller.py \
  --wt-reps 2 \
  --mut-reps 2 \
  --dir-wt wt \
  --dir-mut mut \
  --biotype 1 \
  --meth-diff 0 \
  --vote-threshold 0.5 \
  --threads 4
```

## 11. Recommended reporting for reproducibility

When reporting results generated by MultiDMPcaller, we recommend recording:

- MultiDMPcaller version or GitHub commit hash
- Full command line
- `--biotype`
- q-value thresholds
- `--meth-diff`
- `--vote-threshold`
- Whether adaptive q-value estimation was enabled
- Whether auto-vote estimation was enabled
- Whether low-difference strict voting was enabled
- Whether C++ or Python DMR engine was used
- Number of threads
- Operating system and Python version
