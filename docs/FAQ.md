# Frequently Asked Questions (FAQ)

This page answers common questions about installing, running, and interpreting results from **MultiDMPcaller**.

## General

### Q1. What is MultiDMPcaller used for?

MultiDMPcaller is a downstream DNA methylation analysis tool for detecting and visualizing differentially methylated positions (DMPs) and differentially methylated regions (DMRs) between two biological groups.

It supports CpG, CHG, and CHH contexts and performs all pairwise comparisons between the two groups before applying final support voting.

### Q2. Can MultiDMPcaller detect both DMPs and DMRs?

Yes. Main final outputs are written to `and_output/`:

```text
and_output/{context}-final_significant_sites_DMPs.txt
and_output/{context}-final_significant_sites_DMPs.csv
and_output/{context}-final_significant_regions_DMRs.txt
and_output/{context}-final_significant_regions_DMRs.csv
```

### Q3. Does MultiDMPcaller support plant and animal methylomes?

Yes:

```text
--biotype 0    animal
--biotype 1    plant
--biotype 2    no p-value prefiltering for any context
```

## Input format

### Q4. What is the required input format?

Each replicate file must be a headerless plain-text file with five whitespace- or tab-separated columns:

```text
Chromosome    Position    Methylated_reads    Unmethylated_reads    Context
```

### Q5. What values are allowed in the read-count columns?

Columns 3 and 4 must be non-negative integers. Missing values, negative values, and decimal counts are not valid.

### Q6. What methylation-context labels are accepted?

Recommended labels are:

```text
CpG
CHG
CHH
```

Common `CG`-style labels are normalized to `CpG` where applicable.

### Q7. Do I need separate files for CpG, CHG, and CHH?

No. One replicate file can contain all three contexts.

### Q8. Should input files be sorted?

Sorting by chromosome and genomic position is recommended for consistent file organization and easier manual inspection.

### Q9. Can I use Bismark output directly?

Usually not without conversion. Convert the upstream output to the required five-column format first.

## Directory and file naming

### Q10. What directory structure is required?

Example for a 2 × 2 analysis:

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

### Q11. What are the file naming rules?

Use:

```text
{replicate_index}-{directory_name}.txt
```

Replicate indices must start from `1` and be consecutive.

### Q12. Can I use different input-directory names?

Yes, but each file suffix must match its directory name. For example:

```text
control/1-control.txt
control/2-control.txt
treated/1-treated.txt
treated/2-treated.txt
```

## Command-line usage

### Q13. What is the recommended basic command?

Plant example:

```bash
python MultiDMPcaller.py \
  --wt-reps 2 \
  --mut-reps 2 \
  --dir-wt wt \
  --dir-mut mut \
  --biotype 1
```

Animal example:

```bash
python MultiDMPcaller.py \
  --wt-reps 2 \
  --mut-reps 2 \
  --dir-wt liver \
  --dir-mut brain \
  --biotype 0
```

### Q14. What do `--wt-reps` and `--mut-reps` mean?

They specify the numbers of control/wild-type and experimental/mutant biological replicates.

### Q15. What is the direction of methylation change?

Output direction is interpreted as experimental/mutant relative to control/wild-type. In final DMP and DMR tables, `1` indicates higher methylation in the experimental/mutant group and `0` indicates lower methylation.

## Parameter selection

### Q16. What do `--q-cpg`, `--q-chg`, and `--q-chh` mean?

They are fixed DMP q-value thresholds. Defaults are:

```text
--q-cpg 0.05
--q-chg 0.04
--q-chh 0.045
```

### Q17. What does `--dmr-q` mean?

It specifies the q-value threshold for each pairwise DMR comparison. The default is:

```text
--dmr-q 0.05
```

### Q18. What does `--methy-diff-dmp` mean?

It specifies the minimum absolute site-level methylation difference required for a pairwise comparison to contribute DMP support.

```bash
--methy-diff-dmp 0.2
```

means 20 percentage points. The default is `0.0`.

### Q19. What does `--methy-diff-dmr` mean?

It specifies the minimum absolute regional methylation difference required for each pairwise DMR support.

```bash
--methy-diff-dmr 0.2
```

means 20 percentage points. The regional difference is calculated from aggregated methylated and unmethylated read counts within the candidate region. The default is `0.0`.

### Q20. What does `--vote-threshold` mean?

It defines the support proportion required across all `m × n` pairwise comparisons.

```bash
--vote-threshold 2/3
```

corresponds to the default two-thirds rule. Decimal values are also accepted.

### Q21. What does `--auto-qvalue-twostep` do?

It estimates a data-adaptive DMP q-value threshold for contexts using two-step FDR correction.


### Q22. What do the automatic voting options do?

```text
--auto-dmp-vote-threshold
--auto-dmr-vote-threshold
```

automatically estimate the integer support counts required for final DMP and DMR calling.


### Q23. What happens if GMM-based automatic voting fails?

If fitting or threshold calculation raises an exception, MultiDMPcaller falls back to the support requirement implied by `--vote-threshold`. If the parameter is not supplied, the default two-thirds rule is used.

A convergence warning alone does not necessarily trigger fallback.

### Q24. What is low-difference strict voting?

It is an optional post-voting rule for provisional final DMPs with relatively small boundary methylation differences.

Relevant options are:

```text
--dmp-lowdiff-strict-vote
--dmp-lowdiff-cutoff
```

The default cutoff is `0.3`.

### Q25. What is the difference between `--methy-diff-dmp` and `--dmp-lowdiff-cutoff`?

`--methy-diff-dmp` is a hard filter applied at the pairwise-support layer. `--dmp-lowdiff-cutoff` is used later to identify provisional final DMPs that require stricter voting.

## DMR engine and performance

### Q26. When should I use `--dmr-engine cpp`?

Use it when compatible C++ executables are available. The repository includes:

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

On Windows, compile the source files and add the directory containing `dmr_step1.exe` and `dmr_step2_dynamic.exe` to `PATH`.

The default engine is:

```bash
--dmr-engine python
```

### Q27. What does `--threads` control?

It controls safe parallel stages, including raw-file conversion, replicate-pair processing, and common DMR aggregation. The preprocessing worker count cannot exceed the number of replicate files.

Start with `2` to `4` workers and increase only when memory and disk-I/O bandwidth are sufficient.

### Q28. How does the low-memory preprocessing work?

It reads raw files in chunks and writes final matrices in row blocks. The preprocessing stage can also use `--threads` for replicate-level parallel conversion.

Optional tuning variables are:

```text
MULTIDMPCALLER_NEWTOBOTH_CHUNKSIZE
MULTIDMPCALLER_NEWTOBOTH_BLOCK_ROWS
```

Most users should keep the defaults.

### Q29. Why must the input directories be writable?

The preprocessing stage writes context-specific matrices and temporary spill files into the input directories. It also needs sufficient free disk space.

### Q30. When should I use `--skip-dmr` or `--skip-window`?

Use `--skip-dmr` for DMP-only analysis or quick testing. Use `--skip-window` when only table outputs are required.

## Output interpretation

### Q31. Which files contain final DMP results?

```text
and_output/{context}-final_significant_sites_DMPs.txt
and_output/{context}-final_significant_sites_DMPs.csv
```

### Q32. Which files contain final DMR results?

```text
and_output/{context}-final_significant_regions_DMRs.txt
and_output/{context}-final_significant_regions_DMRs.csv
```

### Q33. Why are pairwise results and final results different?

Pairwise results come from individual replicate comparisons. Final DMPs and DMRs must satisfy support voting across all pairwise comparisons.

### Q34. Why are my results different from methylKit, DSS, DMRcaller, or another tool?

Different tools use different statistical models, filtering rules, multiple-testing corrections, replicate strategies, and DMR definitions. Exact one-to-one agreement is not expected.

## Web server and data availability

### Q35. Can I use the web server instead of local installation?

Yes. The public server is available at:

[https://ciebioinfo.nwafu.edu.cn/](https://ciebioinfo.nwafu.edu.cn/)

It supports upload, parameter selection, job submission, job-ID-based status queries, visualization, and ZIP result download.

### Q36. What are the current web-server limits?

| Item | Current setting |
| :--- | :--- |
| Maximum threads per job | `4` |
| Maximum concurrent jobs | `4` |
| Upload limit | No fixed software-level cap |
| Result retention period | `72 h` |

“No fixed software-level cap” does not imply unlimited upload size. Practical limits can still depend on browser behavior, network conditions, reverse-proxy configuration, available storage, and server resources.

Keep the returned Job ID until the result package has been downloaded.

### Q37. Where are the human simulation benchmark datasets available?

They are available from Zenodo:

[Human simulation benchmark datasets for MultiDMPcaller](https://zenodo.org/records/21121883)

## Reproducibility and issue reporting

### Q38. What should I keep for reproducibility?

Keep:

- the exact command line;
- the software version or GitHub commit hash;
- input checksums;
- the complete `and_output/` directory;
- full logs;
- all threshold settings;
- operating-system and Python-version information.

### Q39. What should I include when reporting a problem?

Include:

- MultiDMPcaller version or commit hash;
- full command line;
- operating system and architecture;
- Python version;
- full error log;
- input directory structure;
- representative input lines;
- replicate counts and contexts;
- available memory and free disk space;
- whether the problem persists with `--threads 1`;
- whether it persists with `--skip-dmr` or `--skip-window`;
- whether the Python or C++ DMR engine was used.
