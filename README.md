## Introduction  
   MultiDMPcaller is a comprehensive tool designed for analyzing both **Differentially Methylated Positions (DMPs)** and **Differentially Methylated Regions (DMRs)** between two groups of samples (e.g., mutant vs. wild-type). It systematically supports the analysis of three cytosine methylation contexts (CpG, CHH, and CHG) and outputs significantly differential sites and regions along with detailed statistical metrics.

## Environment Setup  
```bash
pip install -r requirements.txt
```  


## Running Instructions  

### Data File Format  
1. **File Naming Rules**  
   - Sample files must be placed in two directories (e.g., `msv` and `wt`; directory names can be customized).
   - File name format: `{sample index}-{directory name}.txt`
     Examples: 1-msv.txt (first sample in the mutant group), 2-wt.txt (second sample in the wild-type group).

   **File Content Format**  
   No header. Each line represents a methylation site, separated by spaces or tabs,with exactly 5 columns:  
   ```text
   Chromosome   Position   Methylated_reads   Unmethylated_reads   Context
Note: The program is highly robust with chromosome naming. Both numeric formats (e.g., 1, 2) and prefix formats (e.g., chr1, Chr2) are fully supported and will be automatically normalized.
   ```  
   Example content:  
   ```
   chr1  1005  23  5   CpG
   chr1  1030  18  9   CHH
   chr2  5002  7   32  CHG
   ```  
### Run via Command-Line Arguments  
```bash
python MultiDMPcaller.py <n> <m> <dir_wt> <dir_mut> <biotype>
```  

#### Parameter Description (in order):
1. <n>: Number of samples in the wild-type (WT) group.
2. <m>: Number of samples in the mutant group.  
3. <dir_wt>: Directory name of the wild-type group (e.g., wt).  
4. <dir_mut>: Directory name of the mutant group (e.g., msv).  
5. <biotype>: Biological type filtering logic (0 = animal / 1 = plant / 2 = no filtering).  


### Example Run (with Test Data)  
Using the test data included in the repository ( `wt/1-wt.txt`,`msv/1-msv.txt`), the running command is:  
```bash
python MultiDMPcaller.py 1 1 wt msv 1
```  
(Explanation: 1 wild-type sample, 1 mutant sample, wild-type directory is wt, mutant directory is msv, biological type is plant) .

## Output Results  
Analysis results are saved in the `and_output` folder within the working directory. For user convenience, the core results are provided in both **.txt** (tab-separated) and **.csv** (comma-separated) formats, which can be easily opened in Excel or R.
The core files include:
- `{methylation type}-final_significant_sites_DMPs.csv / .txt`: Significantly differential methylation positions (DMPs) across all replicate combinations, including exact chromosomal positions, methylation change direction, and statistical significance (Q-values).
- `{methylation type}-final_significant_regions_DMRs.csv / .txt`: Significantly differential DMR regions classified by methylation context, including information such as chromosomal boundaries, length, average methylation levels, and significance probability.


## Notes  
1. Ensure that the directory names are consistent with the prefixes in the file names (e.g., all files in the `msv` directory must be named `x-msv.txt`).  
2. The number of samples must strictly match the number of files in the directory (e.g., if 3 samples are input, there must be exactly 3 files in the directory).  
3. If there is no valid data for a specific methylation type, the program will automatically skip it and prompt a message without affecting the analysis of other types.
