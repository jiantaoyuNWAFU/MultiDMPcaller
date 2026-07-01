# Frequently Asked Questions (FAQ)

This page answers common questions about installing, running, and interpreting results from **MultiDMPcaller**.

## General

### Q1. What is MultiDMPcaller used for?

MultiDMPcaller is a downstream DNA methylation analysis tool for detecting and visualizing differentially methylated positions (DMPs) and differentially methylated regions (DMRs) between two biological groups, such as a control/wild-type group and an experimental/mutant group.

It supports the three major cytosine methylation contexts:

- CpG
- CHG
- CHH

It also performs replicate-aware analysis by conducting all pairwise comparisons between the two groups and then applying a voting strategy to identify final DMPs and DMRs.

### Q2. Can MultiDMPcaller detect both DMPs and DMRs?

Yes. MultiDMPcaller reports both site-level DMPs and region-level DMRs.

The main final output files are saved in the `and_output/` directory:

```text
and_output/{context}-final_significant_sites_DMPs.txt
and_output/{context}-final_significant_sites_DMPs.csv
and_output/{context}-final_significant_regions_DMRs.txt
and_output/{context}-final_significant_regions_DMRs.csv
```

For example:

```text
and_output/CpG-final_significant_sites_DMPs.txt
and_output/CHG-final_significant_sites_DMPs.txt
and_output/CHH-final_significant_sites_DMPs.txt
and_output/CpG-final_significant_regions_DMRs.txt
```

### Q3. Does MultiDMPcaller support plant and animal methylomes?

Yes. MultiDMPcaller supports both plant and animal methylome data.

Use the `--biotype` parameter to specify the organism/data mode:

```text
0 = animal
1 = plant
2 = no p-value prefiltering for all contexts
```

In general:

- Use `--biotype 0` for animal WGBS data.
- Use `--biotype 1` for plant WGBS data.
- Use `--biotype 2` when no p-value prefiltering is desired for all contexts.

## Input format

### Q4. What is the required input format?

Each input file should be a plain text file without a header. Each row should contain five whitespace- or tab-separated columns:

```text
Chromosome    Position    Methylated_reads    Unmethylated_reads    Context
```

Example:

```text
chr1    1005    23    5     CpG
chr1    1030    18    9     CHH
chr2    5002    7     32    CHG
```

### Q5. Can I include a header line in the input file?

No. The current input format expects no header line. If your file has a header, remove it before running MultiDMPcaller.

For example:

```bash
tail -n +2 original_file.txt > no_header_file.txt
```

### Q6. What methylation context labels are accepted?

The recommended context labels are:

```text
CpG
CHG
CHH
```

The program also normalizes common CpG/CG-style labels where applicable. However, for clarity and reproducibility, we recommend using `CpG`, `CHG`, and `CHH` consistently in all input files.

### Q7. Do I need to split input files by methylation context?

No. You can provide one input file per replicate containing CpG, CHG, and CHH sites together. MultiDMPcaller automatically processes the three methylation contexts separately.

### Q8. Should chromosome names be written as `1` or `chr1`?

Both styles can be used. Chromosome labels such as `1`, `chr1`, and `Chr1` are normalized internally.

However, it is still recommended to use a consistent chromosome naming style across all replicates to reduce confusion during downstream interpretation.

### Q9. Can I use Bismark output directly?

Usually not directly. MultiDMPcaller requires five columns:

```text
Chromosome    Position    Methylated_reads    Unmethylated_reads    Context
```

If your upstream tool outputs a different format, convert it to the required five-column format before running MultiDMPcaller.

## Directory and file naming

### Q10. What directory structure is required?

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

### Q11. What are the file naming rules?

For each group directory, replicate files should follow this pattern:

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

### Q12. Can I use different names for the input directories?

Yes, but the file suffix must match the directory name.

For example, if your wild-type directory is named `control`, the files should be named:

```text
control/1-control.txt
control/2-control.txt
```

If your mutant directory is named `treated`, the files should be named:

```text
treated/1-treated.txt
treated/2-treated.txt
```

Then run:

```bash
python MultiDMPcaller.py \
  --wt-reps 2 \
  --mut-reps 2 \
  --dir-wt control \
  --dir-mut treated \
  --biotype 1
```

## Command-line usage

### Q13. What is the recommended basic command?

For a 2 vs 2 plant WGBS dataset:

```bash
python MultiDMPcaller.py \
  --wt-reps 2 \
  --mut-reps 2 \
  --dir-wt wt \
  --dir-mut mut \
  --biotype 1
```

For a 2 vs 2 animal WGBS dataset:

```bash
python MultiDMPcaller.py \
  --wt-reps 2 \
  --mut-reps 2 \
  --dir-wt liver \
  --dir-mut brain \
  --biotype 0
```

### Q14. What do `--wt-reps` and `--mut-reps` mean?

`--wt-reps` specifies the number of control/wild-type replicates.

`--mut-reps` specifies the number of experimental/mutant replicates.

For example, if you have:

```text
wt/1-wt.txt
wt/2-wt.txt
mut/1-mut.txt
mut/2-mut.txt
mut/3-mut.txt
```

then use:

```bash
--wt-reps 2 --mut-reps 3
```

### Q15. What is the difference between `--dir-wt` and `--dir-mut`?

`--dir-wt` points to the directory containing control/wild-type samples.

`--dir-mut` points to the directory containing experimental/mutant samples.

The direction of methylation change in the output is interpreted relative to the experimental/mutant group versus the control/wild-type group.

## Parameter selection

### Q16. What do `--q-cpg`, `--q-chg`, and `--q-chh` mean?

These parameters define DMP q-value thresholds for different methylation contexts:

```text
--q-cpg    DMP q-value threshold for CpG sites
--q-chg    DMP q-value threshold for CHG sites
--q-chh    DMP q-value threshold for CHH sites
```

Default values:

```text
--q-cpg 0.05
--q-chg 0.04
--q-chh 0.045
```

### Q17. What does `--dmr-q` mean?

`--dmr-q` specifies the q-value threshold for DMR calling.

Default:

```text
--dmr-q 0.05
```

### Q18. What does `--meth-diff` mean?

`--meth-diff` specifies an optional hard absolute methylation-difference filter for final DMP calling.

For example:

```bash
--meth-diff 0.2
```

means that a site must have at least a 0.2 absolute methylation-level difference, corresponding to 20 percentage points, to be retained as a final DMP.

### Q19. What does `--vote-threshold` mean?

MultiDMPcaller performs all pairwise comparisons between control and experimental replicates. For `m` control replicates and `n` experimental replicates, there are `m × n` pairwise comparisons.

The final DMPs and DMRs are selected using a support/voting rule across these pairwise comparisons.

For example:

```bash
--vote-threshold 0.6666666666666666
```

corresponds to an approximately two-thirds majority rule.

### Q20. What does `--auto-qvalue-twostep` do?

`--auto-qvalue-twostep` enables adaptive q-value threshold estimation for two-step FDR contexts.

This option is intended to make q-value threshold selection more data-adaptive in applicable contexts.

### Q21. What do `--auto-dmp-vote-threshold` and `--auto-dmr-vote-threshold` do?

These options automatically estimate final voting requirements across replicate comparisons:

```text
--auto-dmp-vote-threshold    Automatically estimate the final DMP voting requirement.
--auto-dmr-vote-threshold    Automatically estimate the final DMR voting requirement.
```

### Q22. What is low-difference strict voting?

Low-difference strict voting is an optional post-filter for final DMP candidates with relatively small methylation differences.

Relevant parameters:

```text
--dmp-lowdiff-strict-vote
--dmp-lowdiff-cutoff
```

The cutoff defines the boundary absolute MethDiff used to identify low-difference candidates. For example:

```bash
--dmp-lowdiff-cutoff 0.3
```

### Q23. When should I use `--dmr-engine cpp`?

Use `--dmr-engine cpp` when the accelerated C++ DMR engine is available.

Before using it, make sure the two executable files are in the same directory as the main Python script:

```text
dmr_step1
dmr_step2_dynamic
```

On Linux or macOS, make them executable:

```bash
chmod +x dmr_step1 dmr_step2_dynamic
```

If the C++ engine is unavailable or you want maximum compatibility, use the default Python DMR engine.

### Q24. When should I use `--skip-dmr`?

Use `--skip-dmr` when you only need DMP outputs or when you want to perform a quick DMP-only test.

Example:

```bash
python MultiDMPcaller.py \
  --wt-reps 2 \
  --mut-reps 2 \
  --dir-wt wt \
  --dir-mut mut \
  --biotype 1 \
  --skip-dmr
```

### Q25. When should I use `--skip-window`?

Use `--skip-window` when you want to skip sliding-window visualization and only generate table outputs.

This can reduce runtime for large datasets.

## Output interpretation

### Q26. Which files contain the final DMP results?

Final DMP results are saved in:

```text
and_output/{context}-final_significant_sites_DMPs.txt
and_output/{context}-final_significant_sites_DMPs.csv
```

For example:

```text
and_output/CpG-final_significant_sites_DMPs.txt
and_output/CHG-final_significant_sites_DMPs.txt
and_output/CHH-final_significant_sites_DMPs.txt
```

### Q27. Which files contain the final DMR results?

Final DMR results are saved in:

```text
and_output/{context}-final_significant_regions_DMRs.txt
and_output/{context}-final_significant_regions_DMRs.csv
```

For example:

```text
and_output/CpG-final_significant_regions_DMRs.txt
```

### Q28. What does `Methylation_Change` mean?

In the final DMP table, `Methylation_Change` indicates the direction of methylation change.

In typical output interpretation:

- `1` denotes hyper-methylation in the experimental/mutant group.
- `0` denotes hypo-methylation in the experimental/mutant group.

### Q29. Why are my results different from methylKit, DSS, DMRcaller, or other tools?

Different tools use different statistical models, filtering rules, multiple-testing correction strategies, replicate-handling strategies, and DMR definitions. Therefore, exact one-to-one agreement is not expected.

MultiDMPcaller uses all pairwise comparisons between the two groups and then applies a final voting strategy, which may retain loci that are consistently supported across replicate comparisons.

### Q30. Why are pairwise results and final results different?

Pairwise results are generated from individual replicate comparisons. Final DMPs/DMRs are selected based on support across all pairwise comparisons. Therefore, a site or region significant in one pairwise comparison may not appear in the final output if it does not meet the required voting/support threshold.

## Web server

### Q31. Can I use the web server instead of local installation?

Yes. A public web server is available at:

```text
https://ciebioinfo.nwafu.edu.cn/
```

The web interface supports file upload, parameter selection, job submission, job-ID based status query, result visualization, and ZIP result download.

### Q32. What should I keep after submitting a web-server job?

Keep the returned Job ID. It is required for checking job status and retrieving results later.

## Reproducibility and issue reporting

### Q33. What information should I keep for reproducibility?

For formal analyses, we recommend keeping:

- The exact command line
- The software version or GitHub commit hash
- The input files or their checksums
- The full `and_output/` directory
- The full log file
- Parameter settings
- Runtime environment information

### Q34. What information should I include when reporting a problem?

Please include:

- MultiDMPcaller version or GitHub commit hash
- Full command line
- Operating system
- Python version
- Full error log
- Input directory structure
- First 5 lines of representative input files
- Whether the problem still occurs with `--skip-dmr` or `--skip-window`
