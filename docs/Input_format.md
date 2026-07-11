# Input Format Guide

This page describes the required input format for **MultiDMPcaller**.

## 1. Overview

MultiDMPcaller requires one input file per biological replicate. Files from the control/wild-type group and experimental/mutant group must be placed in two separate directories.

Each input file must be a plain-text file without a header. Each row represents one methylation site.

## 2. Directory structure

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

The directory name and file suffix must match.

If the control directory is named `wt`, the expected files are:

```text
1-wt.txt
2-wt.txt
3-wt.txt
```

If the experimental directory is named `mut`, the expected files are:

```text
1-mut.txt
2-mut.txt
3-mut.txt
```

## 3. File naming rules

Each replicate file must follow:

```text
{replicate_index}-{directory_name}.txt
```

Examples:

```text
wt/1-wt.txt
wt/2-wt.txt
mut/1-mut.txt
mut/2-mut.txt
```

Replicate indices must start from `1` and be consecutive.

Incorrect examples include:

```text
wt/wt1.txt
wt/sample1.txt
wt/1-control.txt
mut/treated_rep1.txt
```

If the original files use different names, rename them before running MultiDMPcaller. On Linux or macOS, symbolic links can also be used:

```bash
mkdir -p wt mut
ln -s original_control_rep1.txt wt/1-wt.txt
ln -s original_control_rep2.txt wt/2-wt.txt
ln -s original_treatment_rep1.txt mut/1-mut.txt
ln -s original_treatment_rep2.txt mut/2-mut.txt
```

## 4. File content format

Each input file must contain five whitespace- or tab-separated columns:

```text
Chromosome    Position    Methylated_reads    Unmethylated_reads    Context
```

Example:

```text
chr1    1005    23    5     CpG
chr1    1030    18    9     CHH
chr2    5002    7     32    CHG
```

## 5. Column descriptions

| Column | Name | Requirement |
| :--- | :--- | :--- |
| 1 | `Chromosome` | Chromosome or scaffold identifier, for example `chr1`, `1`, or `Chr1`. |
| 2 | `Position` | Non-negative integer genomic coordinate. |
| 3 | `Methylated_reads` | Non-negative integer number of methylated reads. |
| 4 | `Unmethylated_reads` | Non-negative integer number of unmethylated reads. |
| 5 | `Context` | Methylation context. Recommended values: `CpG`, `CHG`, and `CHH`. |

Missing values such as `NA`, negative counts, and decimal read counts are not valid.

## 6. Header line

Input files must not contain a header line.

Correct:

```text
chr1    1005    23    5     CpG
chr1    1030    18    9     CHH
```

Incorrect:

```text
Chromosome    Position    Methylated_reads    Unmethylated_reads    Context
chr1          1005        23                  5                    CpG
```

On Linux or macOS, a header can be removed with:

```bash
tail -n +2 input_with_header.txt > input_no_header.txt
```

## 7. Methylation-context labels

Recommended labels are:

```text
CpG
CHG
CHH
```

Common `CG`-style labels are normalized to `CpG` where applicable. For reproducibility, use one consistent label style in all replicates.

## 8. Chromosome naming

Labels such as `1`, `chr1`, and `Chr1` are normalized internally. A consistent naming style across all replicates is nevertheless recommended.

Scaffold identifiers are also accepted, provided the same identifiers are used consistently across samples.

## 9. Read-count requirements

Columns 3 and 4 must contain non-negative integers.

Valid examples:

```text
chr1    1005    23    5     CpG
chr1    1030    0     17    CHH
```

Invalid examples:

```text
chr1    1005    NA    5     CpG
chr1    1030    -1    17    CHH
chr2    5002    7.5   32    CHG
```

On Linux or macOS, check with:

```bash
awk 'NF>=4 && ($3 !~ /^[0-9]+$/ || $4 !~ /^[0-9]+$/) {print NR, $0; exit}' wt/1-wt.txt
```

## 10. Should input files be sorted?

Sorting by chromosome and position is recommended for consistent file organization and easier manual inspection. When available, use the sorting function provided by your upstream methylation-processing workflow.

## 11. Duplicate sites

For best results, each chromosome-position-context combination should appear once per replicate. MultiDMPcaller does not use duplicate rows as a substitute for read-count aggregation. If duplicate records exist, resolve them with the upstream methylation-calling or preprocessing workflow before analysis.

## 12. Do I need separate files for CpG, CHG, and CHH?

No. One replicate file can contain all three contexts. MultiDMPcaller separates CpG, CHG, and CHH internally.

## 13. Writable directories and disk space

The preprocessing stage writes context-specific matrix files and temporary spill files into the input directories. Therefore:

- both input directories must be writable;
- sufficient free disk space must be available;
- the same directory should not be used simultaneously by two independent runs.

Normal completion and ordinary exceptions trigger automatic cleanup of temporary directories. A forced termination can leave stale `newtoboth_*` directories.

## 14. Minimum self-check before running

Verify that:

1. replicate indices start from `1` and are consecutive;
2. each file name matches its directory name;
3. each row has exactly five columns;
4. there is no header line;
5. positions and read counts are non-negative integers;
6. context labels are valid;
7. chromosome names are consistent enough for interpretation;
8. the input directories are writable;
9. sufficient free disk space is available.

On Linux or macOS, useful checks include:

```bash
head -5 wt/1-wt.txt
head -5 mut/1-mut.txt
awk '{print NF}' wt/1-wt.txt | sort | uniq -c
awk '{print $5}' wt/1-wt.txt | sort | uniq -c
```

## 15. Example commands

For a 2 × 2 plant dataset:

```bash
python MultiDMPcaller.py \
  --wt-reps 2 \
  --mut-reps 2 \
  --dir-wt wt \
  --dir-mut mut \
  --biotype 1
```

For a 2 × 2 animal dataset:

```bash
python MultiDMPcaller.py \
  --wt-reps 2 \
  --mut-reps 2 \
  --dir-wt wt \
  --dir-mut mut \
  --biotype 0
```
