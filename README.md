## Introduction  
   MultiDMPcaller is a comprehensive tool designed for analyzing both **Differentially Methylated Positions (DMPs)** and **Differentially Methylated Regions (DMRs)** between two groups of samples (e.g., mutant vs. wild-type). It systematically supports the analysis of three cytosine methylation contexts (CpG, CHH, and CHG) and outputs significantly differential sites and regions along with detailed statistical metrics.

## Environment Setup  
```bash
pip install -r requirements.txt
```  


## Running Instructions  

### 1. Data File Format  
1. **File Naming Rules**  
   - Sample files must be placed in two directories (e.g., `msv` and `wt`; directory names can be customized).
   - File name format: `{sample index}-{directory name}.txt`
     Examples: 1-msv.txt (first sample in the mutant group), 2-wt.txt (second sample in the wild-type group).

   **File Content Format**  
   No header. Each line represents a methylation site, separated by spaces or tabs,with exactly 5 columns:  
   ```text
   Chromosome   Position   Methylated_reads   Unmethylated_reads   Context
   ```
   > **Note:** The program is highly robust with chromosome naming. Both numeric formats (e.g., `1`, `2`) and prefix formats (e.g., `chr1`, `Chr2`) are fully supported and will be automatically normalized.
   **Example Content:**
   ```text
   chr1  1005  23  5   CpG
   chr1  1030  18  9   CHH
   chr2  5002  7   32  CHG
   ```

### 2. Run via Command-Line Arguments  
```bash
python MultiDMPcaller.py <n> <m> <dir_wt> <dir_mut> <biotype>
```  

#### Parameter Description (in order):
| Parameter | Description |
| :--- | :--- |
| `<n>` | Number of replicates in the **Wild-type (WT)** group |
| `<m>` | Number of replicates in the **Mutant** group |
| `<dir_wt>` | Directory name of the WT group (e.g., `wt`) |
| `<dir_mut>` | Directory name of the Mutant group (e.g., `msv`) |
| `<biotype>` | Biological filter: `0`=Animal, `1`=Plant, `2`=No filter |



### 3. Example Run (with Test Data)  
To analyze **1 Wild-type replicate** (in `wt/`) and **1 Mutant replicate** (in `msv/`) for a **plant** genome:
```bash
python MultiDMPcaller.py 1 1 wt msv 1
```  
**Explanation:**
- `n=1`: 1 WT replicate.
- `m=1`: 1 Mutant replicate.
- `dir_wt=wt`: WT files are in `wt/`.
- `dir_mut=msv`: Mutant files are in `msv/`.
- `biotype=1`: Plant mode.


## Output Results  
Analysis results are saved in the `and_output` folder within the working directory. 
The core results are provided in both **.txt** (tab-separated) and **.csv** (comma-separated) formats.

**Core Files:**
*   **`{context}-final_significant_sites_DMPs.csv / .txt`**
    *   Contains significantly differential methylation positions (DMPs) across replicate combinations.
    *   Includes: Chromosomal positions, methylation change direction, and statistical significance (Q-values).

*   **`{context}-final_significant_regions_DMRs.csv / .txt`**
    *   Contains significantly differential DMR regions classified by methylation context.
    *   Includes: Chromosomal boundaries, length, average methylation levels, and significance probability.


## Notes  
1. Ensure that the directory names are consistent with the prefixes in the file names (e.g., all files in the `msv` directory must be named `x-msv.txt`).  
2. The number of samples must strictly match the number of files in the directory (e.g., if 3 samples are input, there must be exactly 3 files in the directory).  
3. If there is no valid data for a specific methylation type, the program will automatically skip it and prompt a message without affecting the analysis of other types.

## Cite MultiDMPcaller
If you use MultiDMPcaller, please cite
```text
Yuan Q., Zhao H., Zhang Z., Yue C., Zhang B., Zou Q. and Yu J. 2026. MultiDMPcaller: A one-stop software for detection and visualization of differentially methylated positions and regions. Bioinformatics (under review).
   ```
