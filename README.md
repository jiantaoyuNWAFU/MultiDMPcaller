## Introduction  
MultiDMPcaller is a tool for analyzing differentially methylated regions (DMRs) between two groups of samples (e.g., mutant vs wild-type). It supports the analysis of three methylation types (CpG, CHH, and CHG) and outputs significantly differential DMR regions along with statistical information.  


## Environment Setup  
```bash
pip install -r requirements.txt
```  


## Running Instructions  

### Data File Format  
1. **File Naming Rules**  
   - Sample files must be placed in two directories (e.g., `msv` and `wt`; directory names can be customized).  
   - File name format: `{sample index}-{directory name}.txt`  
     Examples: `1-msv.txt` (first sample in the first group), `2-wt.txt` (second sample in the second group).  

   **File Content Format**  
   No header. Each line represents a methylation site, separated by spaces or tabs, with 5 columns:  
   ```
   Chromosome number  Position number  Methylated read count  Unmethylated read count  Methylation type (CpG/CHH/CHG)
   ```  
   Example content:  
   ```
   1  1005  23  5  CpG
   1  1030  18  9  CHH
   2  5002  7  32  CHG
   ```  
### Run via Command-Line Arguments  
```bash
python MultiDMPcaller.py <m> <n> <dir1> <dir2> <biotype>
```  

#### Parameter Description (in order):
1. `<m>`: Number of samples in the first group's directory  
2. `<n>`: Number of samples in the second group's directory  
3. `<dir1>`: Directory name of the first group of samples (e.g., `wt` in the example)  
4. `<dir2>`: Directory name of the second group of samples (e.g., `msv` in the example)  
5. `<biotype>`: Biological type (`0`=animal / `1`=plant / `2`=no filtering)  


### Example Run (with Test Data)  
Using the test data included in the repository ( `wt/1-wt.txt`,`msv/1-msv.txt`), the running command is:  
```bash
python MultiDMPcaller.py 1 1 wt msv 1
```  
(Explanation: Number of samples in the second group = 1, number of samples in the first group = 1, first group directory = wt, second group directory = msv, biological type = plant)


## Output Results  
Analysis results are saved in the `and_output` folder in the working directory. The core result files include:  
- `{methylation type}-final_significant_regions_DMRs.txt`: Significantly differential DMR regions classified by methylation type, including information such as chromosomal location, length, methylation level, and significance probability.  


## Notes  
1. Ensure that the directory names are consistent with the prefixes in the file names (e.g., all files in the `msv` directory must be named `x-msv.txt`).  
2. The number of samples must strictly match the number of files in the directory (e.g., if 3 samples are input, there must be exactly 3 files in the directory).  
3. If there is no valid data for a specific methylation type, the program will automatically skip it and prompt a message without affecting the analysis of other types.
