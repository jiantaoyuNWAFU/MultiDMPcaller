# Input Format Guide

This page describes the required input file format for **MultiDMPcaller**.

## 1. Overview

MultiDMPcaller requires one input file per biological replicate. Files from the control/wild-type group and experimental/mutant group should be placed in two separate directories.

Each input file should be a plain text file without a header. Each row represents one methylation site.

## 2. Directory structure

Input files should be placed in two directories:

- One directory for the control/wild-type group
- One directory for the experimental/mutant group

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

The directory name and file suffix must be consistent.

For example, if the control directory is named `wt`, the files inside it should be named:

```text
1-wt.txt
2-wt.txt
3-wt.txt
```

If the experimental directory is named `mut`, the files inside it should be named:

```text
1-mut.txt
2-mut.txt
3-mut.txt
```

## 3. File naming rules

Each replicate file should follow this pattern:

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

Correct examples:

```text
wt/1-wt.txt
wt/2-wt.txt
mut/1-mut.txt
mut/2-mut.txt
```

Incorrect examples:

```text
wt/wt1.txt
wt/sample1.txt
wt/1-control.txt
mut/treated_rep1.txt
```

If your original files have different names, you can either rename them or create symbolic links.

Example:

```bash
mkdir -p wt mut
ln -s original_control_rep1.txt wt/1-wt.txt
ln -s original_control_rep2.txt wt/2-wt.txt
ln -s original_treatment_rep1.txt mut/1-mut.txt
ln -s original_treatment_rep2.txt mut/2-mut.txt
```

## 4. File content format

Each input file should contain five whitespace- or tab-separated columns:

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

| Column | Name | Description |
|---|---|---|
| 1 | `Chromosome` | Chromosome or scaffold identifier. Examples: `chr1`, `1`, `Chr1`. |
| 2 | `Position` | Genomic coordinate of the methylation site. |
| 3 | `Methylated_reads` | Number of reads supporting methylation at this site. |
| 4 | `Unmethylated_reads` | Number of reads supporting unmethylation at this site. |
| 5 | `Context` | Methylation context. Recommended values: `CpG`, `CHG`, `CHH`. |

## 6. Header line

Input files should not contain a header line.

Correct:

```text
chr1    1005    23    5     CpG
chr1    1030    18    9     CHH
chr2    5002    7     32    CHG
```

Incorrect:

```text
Chromosome    Position    Methylated_reads    Unmethylated_reads    Context
chr1          1005        23                  5                    CpG
```

If your file contains a header, remove it:

```bash
tail -n +2 input_with_header.txt > input_no_header.txt
```

## 7. Methylation context labels

Recommended context labels:

```text
CpG
CHG
CHH
```

The program also normalizes common CpG/CG-style labels where applicable. Nevertheless, for clarity, we recommend using `CpG` consistently instead of mixing `CG` and `CpG`.

Check context labels:

```bash
awk '{print $5}' wt/1-wt.txt | sort | uniq -c
```

If your file uses `CG` and you want to convert it to `CpG`:

```bash
awk 'BEGIN{OFS="\t"} {$5=($5=="CG" ? "CpG" : $5); print}' input.txt > converted.txt
```

## 8. Chromosome naming

Chromosome labels such as `1`, `chr1`, and `Chr1` are normalized internally.

However, we recommend using a consistent chromosome naming style across all input files.

Check chromosome names:

```bash
awk '{print $1}' wt/1-wt.txt | sort | uniq | head
awk '{print $1}' mut/1-mut.txt | sort | uniq | head
```

## 9. Read count requirements

Columns 3 and 4 should be numeric read counts.

Valid examples:

```text
chr1    1005    23    5     CpG
chr1    1030    0     17    CHH
```

Invalid examples:

```text
chr1    1005    NA    5     CpG
chr1    1030    18    NA    CHH
```

Check for non-numeric read counts:

```bash
awk 'NF>=4 && ($3 !~ /^[0-9.]+$/ || $4 !~ /^[0-9.]+$/) {print NR, $0; exit}' wt/1-wt.txt
```

## 10. Should input files be sorted?

Sorting by chromosome and position is recommended, especially for large datasets and downstream interpretation.

Example sorting command:

```bash
sort -k1,1 -k2,2n input.txt > input.sorted.txt
```

## 11. Do I need separate files for CpG, CHG, and CHH?

No. You can provide one file per replicate containing all contexts. MultiDMPcaller automatically processes CpG, CHG, and CHH contexts separately.

## 12. Minimum self-check before running

Before running a full analysis, check the input files:

```bash
head -5 wt/1-wt.txt
head -5 mut/1-mut.txt

awk '{print NF}' wt/1-wt.txt | sort | uniq -c
awk '{print NF}' mut/1-mut.txt | sort | uniq -c

awk '{print $5}' wt/1-wt.txt | sort | uniq -c
awk '{print $5}' mut/1-mut.txt | sort | uniq -c
```

A valid row should look like:

```text
chr1    1005    23    5     CpG
```

## 13. Example command after preparing input files

For a 2 vs 2 plant dataset:

```bash
python MultiDMPcaller.py \
  --wt-reps 2 \
  --mut-reps 2 \
  --dir-wt wt \
  --dir-mut mut \
  --biotype 1
```

For a 2 vs 2 animal dataset:

```bash
python MultiDMPcaller.py \
  --wt-reps 2 \
  --mut-reps 2 \
  --dir-wt wt \
  --dir-mut mut \
  --biotype 0
```
