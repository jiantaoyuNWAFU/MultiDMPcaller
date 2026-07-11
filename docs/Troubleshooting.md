# Troubleshooting Guide

This guide summarizes common problems that may occur when running **MultiDMPcaller**, together with likely causes and suggested solutions.

## 1. Installation and environment issues

### 1.1 `ModuleNotFoundError` for `pandas`, `scipy`, `sklearn`, or another package

#### Possible cause

The intended Python environment has not been activated, or the required packages have not been installed.

#### Solution

Create and activate the recommended Conda environment:

```bash
conda create -n multidmpcaller python=3.10 -y
conda activate multidmpcaller
pip install -r requirements.txt
```

Alternatively, create a virtual environment with `venv`:

```bash
python -m venv multidmpcaller_env
```

Activate it on Linux or macOS:

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

Check the installation:

```bash
python --version
python -c "import pandas, scipy, sklearn, matplotlib; print('OK')"
```

Python 3.10 is recommended and was used for validation.

### 1.2 The program uses the wrong Python environment

#### Possible cause

Multiple Python installations or environments are available, and the command is being executed with a different interpreter.

#### Solution

On Linux or macOS, check:

```bash
which python
which pip
python --version
```

In Windows PowerShell, check:

```powershell
Get-Command python
Get-Command pip
python --version
```

In Windows Command Prompt, check:

```bat
where python
where pip
python --version
```

Activate the intended environment before running MultiDMPcaller.

### 1.3 Matplotlib or display-related errors on a server

#### Possible cause

The server does not have a graphical display environment.

#### Solution

If visualization is not required, use:

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

If the error persists, verify that the Python environment has a working non-interactive Matplotlib backend.

## 2. C++ DMR engine issues

### 2.1 `dmr_step1` or `dmr_step2_dynamic` cannot be found

#### Possible cause

The C++ executables are missing, are not executable, or are not discoverable by the program.

The repository includes both source files:

```text
dmr_step1.cpp
dmr_step2_dynamic.cpp
```

and precompiled Linux executables:

```text
dmr_step1
dmr_step2_dynamic
```

#### Solution on Linux or macOS

Place the executables beside `MultiDMPcaller.py` and make them executable:

```bash
chmod +x dmr_step1 dmr_step2_dynamic
```

They can also be compiled from source:

```bash
g++ -O3 -std=c++17 -o dmr_step1 dmr_step1.cpp
g++ -O3 -std=c++17 -o dmr_step2_dynamic dmr_step2_dynamic.cpp
chmod +x dmr_step1 dmr_step2_dynamic
```

#### Solution on Windows

Compile the source files with MinGW-w64 or another compatible `g++` installation:

```powershell
g++ -O3 -std=c++17 -o dmr_step1.exe dmr_step1.cpp
g++ -O3 -std=c++17 -o dmr_step2_dynamic.exe dmr_step2_dynamic.cpp
```

Add the directory containing the compiled `.exe` files to `PATH` before running MultiDMPcaller.

Enable the C++ engine with:

```bash
--dmr-engine cpp
```

The Python DMR engine remains the default:

```bash
--dmr-engine python
```

### 2.2 `Permission denied` when using the C++ DMR engine

#### Possible cause

On Linux or macOS, the executable permission is missing.

#### Solution

```bash
chmod +x dmr_step1 dmr_step2_dynamic
```

Then rerun the command.

### 2.3 The C++ DMR engine fails, but DMP outputs are generated

#### Explanation

DMP calling and DMR calling are separate stages. DMP outputs can be generated successfully even if the DMR engine fails.

#### Suggested checks

First run a DMP-only analysis:

```bash
python MultiDMPcaller.py \
  --wt-reps 2 \
  --mut-reps 2 \
  --dir-wt wt \
  --dir-mut mut \
  --biotype 1 \
  --skip-dmr
```

If the DMP-only run succeeds, inspect the C++ DMR logs and verify that both executables match the current operating system and architecture. You can also compare behavior with:

```bash
--dmr-engine python
```

## 3. Input file and directory problems

### 3.1 Input files cannot be found

#### Possible cause

The directory name and file suffix are inconsistent.

MultiDMPcaller expects:

```text
{replicate_index}-{directory_name}.txt
```

For example:

```text
wt/1-wt.txt
wt/2-wt.txt
mut/1-mut.txt
mut/2-mut.txt
```

#### Solution

Rename the files so that their suffix matches the directory name. On Linux or macOS, symbolic links can also be used:

```bash
ln -s original_control_rep1.txt wt/1-wt.txt
ln -s original_control_rep2.txt wt/2-wt.txt
ln -s original_treatment_rep1.txt mut/1-mut.txt
ln -s original_treatment_rep2.txt mut/2-mut.txt
```

### 3.2 The number of input files does not match `--wt-reps` or `--mut-reps`

If the command contains:

```bash
--wt-reps 2 --mut-reps 3
```

then the input directories must contain:

```text
wt/1-wt.txt
wt/2-wt.txt
mut/1-mut.txt
mut/2-mut.txt
mut/3-mut.txt
```

Replicate indices must start from `1` and be consecutive.

### 3.3 The input file has the wrong number of columns

Each data row must contain five whitespace- or tab-separated fields:

```text
Chromosome    Position    Methylated_reads    Unmethylated_reads    Context
```

On Linux or macOS, check with:

```bash
awk '{print NF}' wt/1-wt.txt | sort | uniq -c
```

Valid rows should have five fields.

### 3.4 The input file contains a header line

MultiDMPcaller expects input files without a header.

On Linux or macOS, remove the first line with:

```bash
tail -n +2 input_with_header.txt > input_no_header.txt
```

### 3.5 Read counts are not valid non-negative integers

Columns 3 and 4 must contain non-negative integer read counts.

Invalid examples include `NA`, negative values, or decimal values.

On Linux or macOS, check with:

```bash
awk 'NF>=4 && ($3 !~ /^[0-9]+$/ || $4 !~ /^[0-9]+$/) {print NR, $0; exit}' wt/1-wt.txt
```

### 3.6 Invalid methylation-context labels

Recommended labels are:

```text
CpG
CHG
CHH
```

Common `CG`-style labels are normalized to `CpG`, but consistent labels are recommended.

### 3.7 Chromosome names differ between replicates

Chromosome labels such as `1`, `chr1`, and `Chr1` are normalized internally. Consistent naming is still recommended for easier interpretation and debugging.

### 3.8 The input directory is not writable or the disk is full

The preprocessing stage writes context-specific matrices and temporary spill files into the input directories. The directories therefore must be writable and must have sufficient free disk space.

On Linux or macOS, check:

```bash
df -h .
ls -ld wt mut
```

On Windows PowerShell, check the target drive with:

```powershell
Get-PSDrive
```

## 4. Parameter-related problems

### 4.1 No or very few DMPs are detected

Possible causes include overly strict q-value, methylation-difference, or voting thresholds; low coverage; weak biological effects; or an inappropriate `--biotype` setting.

Try an exploratory DMP-only run:

```bash
python MultiDMPcaller.py \
  --wt-reps 2 \
  --mut-reps 2 \
  --dir-wt wt \
  --dir-mut mut \
  --biotype 1 \
  --methy-diff-dmp 0 \
  --q-cpg 0.05 \
  --q-chg 0.05 \
  --q-chh 0.05 \
  --vote-threshold 0.5 \
  --skip-dmr \
  --skip-window
```

If this produces DMPs, increase stringency gradually.

### 4.2 Too many DMPs are detected

Use a more stringent site-level difference or voting threshold, for example:

```bash
--methy-diff-dmp 0.2
--vote-threshold 2/3
```

Also inspect replicate-level methylation distributions and sample quality.

### 4.3 Too many or too few DMR supports are detected

The DMR-specific methylation-difference filter is controlled by:

```bash
--methy-diff-dmr
```

For example:

```bash
--methy-diff-dmr 0.2
```

requires an absolute regional methylation difference of at least 20 percentage points for a pairwise DMR support. The regional difference is calculated from aggregated methylated and unmethylated read counts within the candidate region.

### 4.4 Animal data were run with plant mode, or plant data were run with animal mode

Use:

```text
--biotype 0    animal
--biotype 1    plant
--biotype 2    no p-value prefiltering for any context
```

### 4.5 Automatic threshold estimation changed the final output

The following options can change final calling thresholds:

```text
--auto-qvalue-twostep
--auto-dmp-vote-threshold
--auto-dmr-vote-threshold
```

For fully fixed thresholds, omit the automatic options and explicitly set:

```bash
--q-cpg 0.05
--q-chg 0.04
--q-chh 0.045
--dmr-q 0.05
--vote-threshold 2/3
```


### 4.6 Automatic voting fails during GMM fitting or threshold calculation

If the automatic DMP or DMR voting model raises an exception, MultiDMPcaller falls back to the support requirement implied by `--vote-threshold`. If this parameter is not supplied, the default two-thirds rule is used.

A convergence warning alone does not necessarily trigger this fallback; the fallback applies when fitting or threshold calculation raises an actual exception.

## 5. Runtime, memory, and preprocessing problems

### 5.1 The program is slow

Large WGBS files can contain tens of millions of sites. Start with:

```bash
--threads 2
```

or:

```bash
--threads 4
```

The `--threads` value controls safe parallel stages, including raw-file conversion, replicate-pair processing, and common DMR aggregation. More workers can increase disk-I/O pressure and may not always improve runtime.

Other options include:

```bash
--dmr-engine cpp
--skip-dmr
--skip-window
```

### 5.2 The job is killed by the system

Possible causes include insufficient memory, insufficient disk space, excessive parallelism, or external runtime limits.

Try:

```bash
--threads 1
```

or:

```bash
--threads 2
```

The low-memory preprocessing stage can also be tuned with:

```text
MULTIDMPCALLER_NEWTOBOTH_CHUNKSIZE
MULTIDMPCALLER_NEWTOBOTH_BLOCK_ROWS
```

Defaults are `1000000` input rows per chunk and `100000` output rows per block. Lower values can reduce peak memory but increase temporary-file and I/O overhead.

### 5.3 Temporary `newtoboth_*` directories remain after interruption

Normal completion and ordinary exceptions trigger automatic cleanup. A forced termination such as `kill -9`, a power failure, or an external scheduler kill may prevent cleanup.

After confirming that no MultiDMPcaller process is running, residual directories can be located on Linux or macOS with:

```bash
find wt mut -maxdepth 1 -type d -name 'newtoboth_*'
```

Remove only confirmed stale directories.

### 5.4 The terminal disconnects during a long run

On Linux or macOS, use `nohup`:

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

On Windows, use a persistent PowerShell session, Windows Terminal, Task Scheduler, or another process-management method appropriate for the local environment.

## 6. Output-related problems

### 6.1 Final DMP files cannot be found

Check:

```bash
ls -lh and_output/*final_significant_sites_DMPs*
```

Expected names include:

```text
and_output/CpG-final_significant_sites_DMPs.txt
and_output/CHG-final_significant_sites_DMPs.txt
and_output/CHH-final_significant_sites_DMPs.txt
```

A context-specific file may be absent if that context was not present or no output was produced for it.

### 6.2 Final DMR files cannot be found

Check:

```bash
ls -lh and_output/*final_significant_regions_DMRs*
```

DMR files are not generated when `--skip-dmr` is used. If DMR calling failed, inspect the run log and C++ engine logs if applicable.

### 6.3 Visualization outputs are missing

Visualization is not generated when:

```bash
--skip-window
```

is used.

### 6.4 Pairwise outputs exist, but final DMP or DMR files are empty

A site or region can be significant in one pairwise comparison but fail the final support requirement. For exploratory analysis, compare with a more permissive setting such as:

```bash
--vote-threshold 0.5
```

### 6.5 Results differ from other tools

Different methylation-analysis tools use different statistical tests, replicate models, filtering rules, multiple-testing procedures, and DMR definitions. Exact agreement with methylKit, DSS, DMRcaller, methylSig, DMRcate, or other tools is not expected.

MultiDMPcaller performs all pairwise comparisons between the two groups and then applies final support voting.

## 7. Quick self-check before running

Confirm that:

- the replicate counts match `--wt-reps` and `--mut-reps`;
- file names match their directory names;
- each row has five columns and no header;
- read counts are non-negative integers;
- context labels are valid;
- input directories are writable;
- sufficient disk space is available;
- the selected `--biotype` is appropriate.

A small DMP-only test can be run with:

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

Please include:

- MultiDMPcaller version or GitHub commit hash;
- full command line;
- operating system and architecture;
- Python version;
- whether Conda or `venv` was used;
- full error log;
- input directory structure;
- first five lines of representative input files;
- replicate counts and methylation contexts;
- available memory and free disk space;
- whether the problem persists with `--threads 1`;
- whether it persists with `--skip-dmr` or `--skip-window`;
- whether the Python or C++ DMR engine was used.
