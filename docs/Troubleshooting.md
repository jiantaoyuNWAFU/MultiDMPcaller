# Troubleshooting Guide

This guide summarizes common problems that may occur when running **MultiDMPcaller**, together with possible causes and suggested solutions.

## 1. Installation and environment issues

### 1.1 `ModuleNotFoundError: No module named 'pandas'`, `scipy`, `sklearn`, or other packages

#### Possible cause

The Python environment has not been activated, or required packages have not been installed.

#### Solution

Create and activate the recommended conda environment:

```bash
conda create -n multidmpcaller python=3.10 -y
conda activate multidmpcaller
pip install -r requirements.txt
```

Alternatively, use `venv`:

```bash
python -m venv multidmpcaller_env
source multidmpcaller_env/bin/activate
pip install -r requirements.txt
```

Then check:

```bash
python --version
python -c "import pandas, scipy, sklearn, matplotlib; print('OK')"
```

### 1.2 The program uses the wrong Python environment

#### Possible cause

Multiple Python environments are installed, and the command is executed under a different environment.

#### Solution

Check the active Python interpreter:

```bash
which python
python --version
which pip
```

If using conda, activate the environment before running:

```bash
conda activate multidmpcaller
```

### 1.3 Matplotlib or display-related errors on a server

#### Possible cause

The server does not have a graphical display environment.

#### Solution

If you do not need visualization, use:

```bash
--skip-window
```

Example:

```bash
python MultiDMPcaller.py \
  --wt-reps 2 \
  --mut-reps 2 \
  --dir-wt wt \
  --dir-mut mut \
  --biotype 1 \
  --skip-window
```

If the error persists, check whether your Python environment has a working non-interactive matplotlib backend.

## 2. C++ DMR engine issues

### 2.1 `dmr_step1` or `dmr_step2_dynamic` not found

#### Possible cause

The accelerated C++ DMR engine files are missing or not placed in the same directory as the main Python script.

#### Solution

Make sure these two files are located together with `MultiDMPcaller.py`:

```text
dmr_step1
dmr_step2_dynamic
```

Then run with:

```bash
--dmr-engine cpp
```

If the files are not available, use the default Python DMR engine by omitting `--dmr-engine cpp`.

### 2.2 `Permission denied` when using the C++ DMR engine

#### Possible cause

The C++ executable files do not have executable permission.

#### Solution

On Linux or macOS, run:

```bash
chmod +x dmr_step1 dmr_step2_dynamic
```

Then rerun the command.

### 2.3 The C++ DMR engine fails, but DMP outputs are generated

#### Possible cause

The DMP calling step and DMR calling step are separate. DMP outputs may be generated successfully even if the DMR engine fails.

#### Solution

First verify DMP-only analysis:

```bash
python MultiDMPcaller.py \
  --wt-reps 2 \
  --mut-reps 2 \
  --dir-wt wt \
  --dir-mut mut \
  --biotype 1 \
  --skip-dmr
```

If DMP-only analysis works, the problem is likely specific to DMR calling or the C++ DMR engine. Try the Python DMR engine by omitting:

```bash
--dmr-engine cpp
```

## 3. Input file and directory problems

### 3.1 Input files cannot be found

#### Possible cause

The directory name and file suffix are inconsistent.

MultiDMPcaller expects files to follow this pattern:

```text
{replicate_index}-{directory_name}.txt
```

For example, if the directory is named `wt`, expected files are:

```text
wt/1-wt.txt
wt/2-wt.txt
```

If the directory is named `mut`, expected files are:

```text
mut/1-mut.txt
mut/2-mut.txt
```

#### Solution

Check the directory structure:

```bash
tree wt mut
```

If necessary, rename files or create symbolic links:

```bash
ln -s original_control_rep1.txt wt/1-wt.txt
ln -s original_control_rep2.txt wt/2-wt.txt
ln -s original_treatment_rep1.txt mut/1-mut.txt
ln -s original_treatment_rep2.txt mut/2-mut.txt
```

### 3.2 The number of input files does not match `--wt-reps` or `--mut-reps`

#### Possible cause

The number of replicate files in the input directory does not match the command-line arguments.

#### Solution

Check the files:

```bash
ls -lh wt
ls -lh mut
```

For example, if you specify:

```bash
--wt-reps 2 --mut-reps 3
```

you should have:

```text
wt/1-wt.txt
wt/2-wt.txt
mut/1-mut.txt
mut/2-mut.txt
mut/3-mut.txt
```

### 3.3 Replicate indices are not consecutive

#### Possible cause

Replicate indices should start from 1 and be consecutive. For example, `1-wt.txt` and `3-wt.txt` without `2-wt.txt` may cause problems.

#### Solution

Rename or link the files so that indices are consecutive:

```text
1-wt.txt
2-wt.txt
3-wt.txt
```

### 3.4 The input file has the wrong number of columns

#### Possible cause

Each input row should contain five fields:

```text
Chromosome    Position    Methylated_reads    Unmethylated_reads    Context
```

#### Solution

Check the first few lines:

```bash
head -5 wt/1-wt.txt
head -5 mut/1-mut.txt
```

Check the number of columns:

```bash
awk '{print NF}' wt/1-wt.txt | sort | uniq -c
awk '{print NF}' mut/1-mut.txt | sort | uniq -c
```

Valid data rows should have five fields.

### 3.5 The input file contains a header line

#### Possible cause

MultiDMPcaller expects input files without a header.

#### Solution

Remove the header:

```bash
tail -n +2 input_with_header.txt > input_no_header.txt
```

Then use the no-header file as input.

### 3.6 Methylated or unmethylated read counts are not numeric

#### Possible cause

Columns 3 and 4 should contain numeric read counts.

#### Solution

Check the file:

```bash
awk 'NF>=4 && ($3 !~ /^[0-9.]+$/ || $4 !~ /^[0-9.]+$/) {print NR, $0; exit}' wt/1-wt.txt
```

If problematic rows exist, fix or remove them before running MultiDMPcaller.

### 3.7 Invalid methylation context labels

#### Possible cause

The context column contains labels other than `CpG`, `CHG`, or `CHH`.

#### Solution

Check the context labels:

```bash
awk '{print $5}' wt/1-wt.txt | sort | uniq -c
awk '{print $5}' mut/1-mut.txt | sort | uniq -c
```

Recommended labels:

```text
CpG
CHG
CHH
```

If your upstream tool uses `CG`, convert it to `CpG` if needed:

```bash
awk 'BEGIN{OFS="\t"} {$5=($5=="CG" ? "CpG" : $5); print}' input.txt > converted.txt
```

### 3.8 Chromosome names differ between replicates

#### Possible cause

Some files use `chr1`, while others use `1`, `Chr1`, or other naming conventions.

MultiDMPcaller normalizes common chromosome labels internally, but consistent naming is still recommended.

#### Solution

Inspect chromosome names:

```bash
awk '{print $1}' wt/1-wt.txt | sort | uniq | head
awk '{print $1}' mut/1-mut.txt | sort | uniq | head
```

If needed, standardize chromosome labels before analysis.

## 4. Parameter-related problems

### 4.1 No or very few DMPs are detected

#### Possible causes

- The selected organism/data mode is not appropriate.
- The q-value thresholds are too strict.
- The methylation-difference threshold is too strict.
- The vote threshold is too strict.
- Input coverage is low.
- The biological difference between groups is weak.
- Chromosome names or context labels are inconsistent.

#### Suggested checks

Start with a DMP-only exploratory run:

```bash
python MultiDMPcaller.py \
  --wt-reps 2 \
  --mut-reps 2 \
  --dir-wt wt \
  --dir-mut mut \
  --biotype 1 \
  --meth-diff 0 \
  --q-cpg 0.05 \
  --q-chg 0.05 \
  --q-chh 0.05 \
  --vote-threshold 0.5 \
  --skip-dmr \
  --skip-window
```

If this produces DMPs, gradually increase stringency.

### 4.2 Too many DMPs are detected

#### Possible causes

- `--meth-diff` is too low.
- `--vote-threshold` is too permissive.
- There may be replicate-level quality differences.
- The dataset may contain strong global methylation differences.

#### Suggested solution

Use stricter thresholds:

```bash
--meth-diff 0.2
--vote-threshold 0.6666666666666666
```

You may also inspect replicate-level methylation distributions before final interpretation.

### 4.3 Animal data were accidentally run with plant mode

#### Possible cause

Animal WGBS data should generally use:

```bash
--biotype 0
```

If animal data are run with plant mode, the context-specific two-step strategy may not match the intended data mode.

#### Solution

Rerun using:

```bash
--biotype 0
```

### 4.4 Plant data were accidentally run with animal mode

#### Possible cause

Plant WGBS data should generally use:

```bash
--biotype 1
```

This is especially important for CHG and CHH methylation contexts.

#### Solution

Rerun using:

```bash
--biotype 1
```

### 4.5 Threshold-estimation options changed the final output unexpectedly

#### Possible cause

Options such as `--auto-qvalue-twostep`, `--auto-dmp-vote-threshold`, and `--auto-dmr-vote-threshold` may change thresholds used in final calling.

#### Solution

For formal analysis, record whether these options were enabled. If you need fully fixed thresholds, omit automatic threshold-estimation options and explicitly specify:

```bash
--q-cpg 0.05
--q-chg 0.04
--q-chh 0.045
--dmr-q 0.05
--vote-threshold 0.6666666666666666
```

## 5. Runtime and memory problems

### 5.1 The program is slow

#### Possible causes

- Large WGBS datasets contain tens of millions of methylation sites.
- DMR calling can be time-consuming.
- Visualization can add extra runtime.
- Too many threads may increase I/O pressure on some systems.

#### Suggested solutions

Use a moderate number of threads first:

```bash
--threads 4
```

Use the C++ DMR engine if available:

```bash
--dmr-engine cpp
```

Skip DMR calling during DMP-only tests:

```bash
--skip-dmr
```

Skip visualization if only tables are needed:

```bash
--skip-window
```

### 5.2 The job is killed by the system

#### Possible cause

The job may exceed available memory or runtime limits.

#### Suggested solutions

Use fewer threads:

```bash
--threads 2
```

Skip optional steps:

```bash
--skip-dmr
--skip-window
```

Run the job on a machine with more memory if possible.

### 5.3 The terminal disconnects during a long run

#### Possible cause

Long jobs may continue running or be terminated depending on the shell/session environment.

#### Suggested solution

Use `nohup` for long-running analyses:

```bash
nohup python MultiDMPcaller.py \
  --wt-reps 2 \
  --mut-reps 2 \
  --dir-wt wt \
  --dir-mut mut \
  --biotype 1 \
  --threads 4 \
  --dmr-engine cpp \
  > run.log 2>&1 &
```

Monitor the log:

```bash
tail -f run.log
```

## 6. Output-related problems

### 6.1 I cannot find the final DMP files

#### Solution

Check the `and_output/` directory:

```bash
ls -lh and_output/*final_significant_sites_DMPs*
```

Expected files include:

```text
and_output/CpG-final_significant_sites_DMPs.txt
and_output/CHG-final_significant_sites_DMPs.txt
and_output/CHH-final_significant_sites_DMPs.txt
```

### 6.2 I cannot find the final DMR files

#### Possible cause

DMR calling may have been skipped, failed, or produced no significant DMRs.

#### Solution

Check:

```bash
ls -lh and_output/*final_significant_regions_DMRs*
```

If you used `--skip-dmr`, DMR files will not be generated.

If DMR calling failed, inspect the log file.

### 6.3 Visualization outputs are missing

#### Possible cause

Visualization was skipped with:

```bash
--skip-window
```

#### Solution

Rerun without `--skip-window` if visualization files are needed.

### 6.4 Pairwise outputs exist, but final DMP/DMR files are empty

#### Possible cause

Some sites or regions may be significant in individual pairwise comparisons but fail to meet the final voting/support threshold.

#### Solution

Try a more permissive exploratory vote threshold:

```bash
--vote-threshold 0.5
```

Then compare with stricter settings.

### 6.5 Results are different from other tools

#### Explanation

Different methylation analysis tools use different statistical tests, replicate models, filtering criteria, multiple-testing correction procedures, and DMR definitions. Exact agreement with methylKit, DSS, DMRcaller, methylSig, DMRcate, or other tools should not be expected.

MultiDMPcaller uses all pairwise comparisons between groups and applies a final voting strategy across comparisons.

## 7. Quick self-check before running

Before running MultiDMPcaller, check the following items.

### 7.1 Check directory structure

```bash
tree wt mut
```

Expected example:

```text
wt/
├── 1-wt.txt
└── 2-wt.txt
mut/
├── 1-mut.txt
└── 2-mut.txt
```

### 7.2 Check first lines of input files

```bash
head -5 wt/1-wt.txt
head -5 mut/1-mut.txt
```

Expected format:

```text
chr1    1005    23    5     CpG
chr1    1030    18    9     CHH
chr2    5002    7     32    CHG
```

### 7.3 Check column counts

```bash
awk '{print NF}' wt/1-wt.txt | sort | uniq -c
awk '{print NF}' mut/1-mut.txt | sort | uniq -c
```

Most valid rows should have five columns.

### 7.4 Check methylation contexts

```bash
awk '{print $5}' wt/1-wt.txt | sort | uniq -c
awk '{print $5}' mut/1-mut.txt | sort | uniq -c
```

Expected labels:

```text
CpG
CHG
CHH
```

### 7.5 Run a small DMP-only test

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

## 8. Information to include when reporting an issue

When opening a GitHub issue or contacting the developers, please include:

- MultiDMPcaller version or GitHub commit hash
- Full command line
- Operating system
- Python version
- Whether conda or venv was used
- Full error log
- Input directory structure
- First 5 lines of representative input files
- Output of `awk '{print NF}' file | sort | uniq -c`
- Output of `awk '{print $5}' file | sort | uniq -c`
- Whether the problem still occurs with `--skip-dmr`
- Whether the problem still occurs with `--skip-window`
