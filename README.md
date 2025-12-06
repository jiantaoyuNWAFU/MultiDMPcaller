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

2. Run the core code:  
```bash
python MultiDMPcaller.py
```  

3. Follow the prompts to input the data directories, number of samples, and biological type (0=animal / 1=plant / 2=no filtering).  


## Example Data  
`msv/1-msv.txt` and `wt/1-wt.txt` are test data files, containing basic methylation site formats, which can be directly used for test runs.  


## Output Results  
Analysis results are saved in the `and_output` folder in the working directory. The core result files include:  
- `{methylation type}-final_significant_regions_DMRs.txt`: Significantly differential DMR regions classified by methylation type, including information such as chromosomal location, length, methylation level, and significance probability.  


## Notes  
1. Ensure that the directory names are consistent with the prefixes in the file names (e.g., all files in the `msv` directory must be named `x-msv.txt`).  
2. The number of samples must strictly match the number of files in the directory (e.g., if 3 samples are input, there must be exactly 3 files in the directory).  
3. If there is no valid data for a specific methylation type, the program will automatically skip it and prompt a message without affecting the analysis of other types.
