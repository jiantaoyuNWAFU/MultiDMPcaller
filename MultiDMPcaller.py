import sys
import os
import re
import pandas as pd
import numpy as np
from scipy.stats import fisher_exact
import glob
import time
from collections import defaultdict
from pathlib import Path
from typing import Union, Tuple, Optional, List
import matplotlib
import matplotlib.pyplot as plt
from scipy.special import betainc
import bisect
from collections import namedtuple

# Configuration of significance thresholds for DMP Q-values
DMP_QVALUE_THRESHOLDS = {
    'CpG': 0.05,
    'CHH': 0.045,
    'CHG': 0.04,
}

# Significance threshold for DMR Q-values
DMR_QVALUE_THRESHOLD = 0.05

def get_dmp_threshold(methylation_type):
    """Get the DMP Q-value significance threshold for the specified methylation type"""
    return DMP_QVALUE_THRESHOLDS.get(methylation_type, 0.05)

DmrRecord = namedtuple('DmrRecord', ['exp_methy', 'exp_unmethy', 'wild_methy', 'wild_unmethy', 'qvalue', 'direction'])

def process_common_sites_dmr_and_summarize(dir1, dir2, m, n, methylation_types=['CpG', 'CHH', 'CHG'],work_dir="."):
    """
    Complete workflow: Process DMRs of common significant sites + cyclic summation analysis + summarization

    parameter：
        dir1: Directory of the first sample group (mutant)
        dir2: Directory of the second sample group (wild-type)
        m: Number of mutant samples
        n: Number of wild-type samples
        methylation_types: List of methylation types
    """
    # Step 1: Process common significant sites to generate DMRs
    print("Phase 3: Process DMRs of Common Significant Sites")

    for mtype in methylation_types:
        dmr_results = process_common_sites_to_dmr(methylation_type=mtype,work_dir=work_dir)
        if dmr_results:
            print(f"\n{mtype} type DMR analysis completed. DMR results obtained for {len(dmr_results)} chromosomes")

    # Step 2: Perform summation analysis and summarization for DMRs in all output_x_y directories
    summarize_all_dmr_methylation(dir1, dir2, m, n, methylation_types,work_dir=work_dir)

def generate_final_significant_dmr(dmr_data_dict, methylation_type, output_dir,m,n,dir1,dir2):
    """
    Generate the final significant DMR file (based on Bayesian decision)

    parameter：
        dmr_data_dict: Mapping dictionary: Chromosome -> DMR region (start, end) -> Data list (statistical data of each site in multiple tests)
        methylation_type: methylation type
        output_dir: Output directory
        threshold: Bayesian decision threshold (default: 2/3)
    """
    print(f"\nProcessing final DMR summarization for {methylation_type} ...")

    if not dmr_data_dict:
        print(f"  No DMR data for {methylation_type}")
        return

    final_dmr_list = []

    # Iterate over all chromosomes
    for chr_num in sorted(dmr_data_dict.keys()):
        print(f"  Processing chromosome {chr_num}...")
        chr_dmr_dict = dmr_data_dict[chr_num]  # Get the mapping in the current chromosome: DMR region (start, end) -> Data list (statistical data of each site in multiple tests)
                                                                    # The data list contains m*n elements, representing the statistical information of the region in each test

        # Iterate through all DMRs on this chromosome
        for dmr_key, data_list in chr_dmr_dict.items():
            start, end = dmr_key

            # Count significant occurrences
            sig_count = sum(1 for item in data_list if item.qvalue <= 0.05)
            total_count = len(data_list)

            # Determine significance
            is_significant = bayes_deciding(sig_count, total_count - sig_count)

            if not is_significant:
                continue

            # Calculate average values
            avg_exp_m = np.mean([item.exp_methy for item in data_list])
            avg_exp_u = np.mean([item.exp_unmethy for item in data_list])
            avg_wild_m = np.mean([item.wild_methy for item in data_list])
            avg_wild_u = np.mean([item.wild_unmethy for item in data_list])

            filtered_values = [item.qvalue for item in data_list if item.qvalue <= 0.05]
            avg_qvalue = np.mean(filtered_values) if filtered_values else 1

            # Determine direction by voting (majority vote)
            direction_votes = [item.direction for item in data_list]
            direction = 1 if sum(direction_votes) >= len(direction_votes) / 2 else 0

            # Calculate probability
            prob = sig_count / total_count

            dmr_record = {
                'Chromosome': chr_num,
                'Methylation_Type': methylation_type,
                'DMR_start': start,
                'DMR_end': end,
                'Length': end - start + 1,
                'Avg_exp_methy': avg_exp_m,
                'Avg_exp_unmethy': avg_exp_u,
                'Avg_wild_methy': avg_wild_m,
                'Avg_wild_unmethy': avg_wild_u,
                'Sig_Avg_qvalue': avg_qvalue,
                'Direction': direction,
                'Sig_count': sig_count,
                'Total_count': total_count,
                'Sig_probability': prob
            }
            # Add q-values for all replicates (in sequential order)
            idx = 0
            for x in range(1, m + 1):
                for y in range(1, n + 1):
                    col_name = f'qvalue_{os.path.basename(dir2.rstrip("/"))}{y}_{os.path.basename(dir1.rstrip("/"))}{x}'
                    if idx < len(data_list):
                        dmr_record[col_name] = data_list[idx].qvalue
                    else:
                        dmr_record[col_name] = 1.0
                    idx += 1
            final_dmr_list.append(dmr_record)

    if not final_dmr_list:
        print(f"  {methylation_type} No significant DMRs found")
        return

    # Convert to DataFrame and sort
    final_df = pd.DataFrame(final_dmr_list)
    final_df = final_df.sort_values(['Chromosome', 'DMR_start'])
    column_order = [
        'Chromosome', 'Methylation_Type', 'DMR_start', 'DMR_end', 'Length',
        'Direction', 'Sig_count', 'Total_count', 'Sig_probability',
        'Avg_exp_methy', 'Avg_exp_unmethy', 'Avg_wild_methy', 'Avg_wild_unmethy',
        'Sig_Avg_qvalue'
    ]
    replicate_columns = sorted(
        [col for col in final_df.columns if col.startswith('qvalue_')],
        key=lambda x: tuple(map(int, re.findall(r'\d+', x)))
    )
    column_order = column_order + replicate_columns
    final_df = final_df[column_order]
    # Save results
    output_file = os.path.join(output_dir, f"{methylation_type}-final_significant_regions_DMRs.txt")
    final_df.to_csv(output_file, sep='\t', index=False)

    print(f"  Final significant DMRs for {methylation_type} : {len(final_df)} ")
    print(f"  Saved to: {output_file}")

    # Summary statistics
    hyper_count = len(final_df[final_df['Direction'] == 1])
    hypo_count = len(final_df[final_df['Direction'] == 0])
    print(f"    - Hypermethylated: {hyper_count} ({hyper_count / len(final_df) * 100:.1f}%)")
    print(f"    - Hypomethylated : {hypo_count} ({hypo_count / len(final_df) * 100:.1f}%)")

    return final_df

def collect_dmr_results(methy_dir, methylation_type, all_dmr_results):
    """
    Collect results from a single DMR analysis

    Parameters：
        methy_dir: Directory for methylation type (e.g., ./and_output/CpG/)
        methylation_type: Methylation type
        all_dmr_results: Summary dictionary; initially contains only three methylation type keys
    """
    # Find all dmr_fisher_Chr*.txt files
    fisher_files = glob.glob(os.path.join(methy_dir, "dmr_fisher_Chr*.txt"))

    for fisher_file in fisher_files:
        # Extract chromosome number
        match = re.search(r'Chr(\d+)\.txt$', fisher_file)
        if not match:
            continue
        chr_num = int(match.group(1))

        # Read file
        try:
            df = pd.read_csv(fisher_file, sep='\t')

            if df.empty:
                continue

            # Initialize the dictionary for this chromosome
            if chr_num not in all_dmr_results[methylation_type]:
                all_dmr_results[methylation_type][chr_num] = defaultdict(list)

            # Collect data for each DMR
            for _, row in df.iterrows():
                dmr_key = (int(row['DMR_start']), int(row['DMR_end']))

                # Store (exp_m, exp_u, wild_m, wild_u, qvalue, direction)
                dmr_data = DmrRecord(
                    exp_methy=int(row['exp_methy_sum']),
                    exp_unmethy=int(row['exp_unmethy_sum']),
                    wild_methy=int(row['wild_methy_sum']),
                    wild_unmethy=int(row['wild_unmethy_sum']),
                    qvalue=float(row['qvalue']) if not pd.isna(row['qvalue']) else 1.0,
                    direction=int(row['direction'])
                )

                all_dmr_results[methylation_type][chr_num][dmr_key].append(dmr_data)

        except Exception as e:
            print(f"    Warning: Failed to read {fisher_file} : {e}")
            continue

def summarize_all_dmr_methylation(dir1, dir2, m, n, methylation_types=['CpG', 'CHH', 'CHG'],work_dir="."):
    """
    Perform summation analysis of methylation reads using common DMRs across all output_x_y directories
    """
    print("Phase 4: Summation analysis of methylation reads for common DMRs across all combinations")

    and_output_dir = os.path.join(work_dir, "and_output")
    os.makedirs(and_output_dir, exist_ok=True)

    # Store methylation data for all DMRs
    all_dmr_results = {mtype: {} for mtype in methylation_types}

    # 1. Loop through all output_x_y combinations
    for replicate_x in range(1, m + 1):
        for replicate_y in range(1, n + 1):
            print(f"\nProcessing combination ({dir1}{replicate_x}, {dir2}{replicate_y})...")

            # Get paths for both files
            file1_path = os.path.join(dir1, f"{replicate_x}-bothMeUnme_diffChromo_NOREPEATED_methy_sites_{{}}.txt")
            file2_path = os.path.join(dir2, f"{replicate_y}-bothMeUnme_diffChromo_NOREPEATED_methy_sites_{{}}.txt")

            # 2. Loop through the three methylation types
            for mtype in methylation_types:
                # Check if both files exist
                f1 = file1_path.format(mtype)
                f2 = file2_path.format(mtype)
                if not os.path.exists(f1) or not os.path.exists(f2):
                    print(f"  Skipping {mtype}:both files do not exist")
                    continue

                print(f"  Processing {mtype} type...")

                # Create output directory (saved under and_output)
                methy_output_dir = os.path.join(and_output_dir, f"dmr_analysis_wt{replicate_y}_mut{replicate_x}", mtype)
                os.makedirs(methy_output_dir, exist_ok=True)

                # 3. Call the modified summarize_dmr_methylation function, passing custom_dmr_dir
                try:
                    summarize_dmr_methylation(
                        methy_dir=methy_output_dir,  # Results are saved to a new directory under and_output
                        replicate_x=replicate_x,
                        replicate_y=replicate_y,
                        file1_path=f1,
                        file2_path=f2,
                        methylation_type=mtype,
                        custom_dmr_dir=and_output_dir  #  Read common DMRs from and_output
                    )

                    # 4. Read the generated dmr_fisher files and collect results
                    collect_dmr_results(methy_output_dir, mtype, all_dmr_results)

                except Exception as e:
                    print(f"  Error:Failed to process {mtype} : {e}")
                    import traceback
                    traceback.print_exc()
                    continue

    # 5. Summarize all results to generate final significant DMRs
    print("Phase 5: Summarize DMR results and perform decision")

    for mtype in methylation_types:
        generate_final_significant_dmr(
            all_dmr_results[mtype],
            methylation_type=mtype,
            output_dir=and_output_dir,
            m=m,
            n=n,
            dir1=dir1,
            dir2=dir2
        )

def process_common_sites_to_dmr(methylation_type='CpG',work_dir="."):
    """
    Read the final_significant_sites_DMPs.txt file and process DMRs grouped by chromosome

    parameter：

    methylation_type: Methylation type (CpG, CHH, CHG)
    Returns：
        dmr_results: Dictionary {chr_num: dmr_list_file_path}
    """
    print(f"\nStarting common significant site DMR analysis for {methylation_type} ...")

    and_output_dir = os.path.join(work_dir, "and_output")
    os.makedirs(and_output_dir, exist_ok=True)

    # 1.Read the final_significant_sites_DMPs file
    common_file = os.path.join(and_output_dir, f"{methylation_type}-final_significant_sites_DMPs.txt")
    if not os.path.exists(common_file):
        print(f"Error: File does not exist {common_file}")
        return None

    try:
        df = pd.read_csv(common_file, sep='\t')
        print(f"Successfully read file, total {len(df)} sites")
    except Exception as e:
        print(f"Failed to read file: {e}")
        return None

    if df.empty:
        print(f"Warning:{methylation_type} file is empty")
        return None

    # 2. Group by chromosome
    chromosomes = sorted(df['Chromosome'].unique(), key=natural_sort_key)
    print(f"Found {len(chromosomes)} chromosomes: {chromosomes}")

    dmr_results = {}
    total_chromosomes = len(chromosomes)

    # 3. Process each chromosome
    for chr_num in chromosomes:
        print(f"\n  Processing chromosome {chr_num} ({chr_num}/{total_chromosomes})...")

        # Filter data for the current chromosome
        chr_df = df[df['Chromosome'] == chr_num].copy()

        # Extract the required three columns:Position, Sig_Mean_Qvalue, Methylation_Change
        dmp_data = chr_df[['Position', 'Sig_Mean_Qvalue', 'Methylation_Change']].copy()

        # Ensure correct data types
        dmp_data['Position'] = dmp_data['Position'].astype(int)
        dmp_data['Sig_Mean_Qvalue'] = dmp_data['Sig_Mean_Qvalue'].astype(float)
        dmp_data['Methylation_Change'] = dmp_data['Methylation_Change'].astype(int)

        # Sort by position
        dmp_data = dmp_data.sort_values('Position').reset_index(drop=True)

        print(f"    Chromosome {chr_num} has a total of {len(dmp_data)} sites")

        if len(dmp_data) == 0:
            print(f"    Skipping: Chromosome {chr_num} has no valid sites")
            continue

        # 4. Create temporary DMP file (format: pos qvalue change)
        temp_dmp_file = os.path.join(and_output_dir,
                                     f"DMP_common_{methylation_type}_Chr{chr_num}.txt")

        # Write to DMP format file (first line is "first line", followed by pos qvalue change)
        with open(temp_dmp_file, 'w') as f:
            f.write("first line\n")
            for _, row in dmp_data.iterrows():
                f.write(f"{int(float(row['Position']))} {float(row['Sig_Mean_Qvalue'])} {int(float(row['Methylation_Change']))}\n")

        print(f"    Created DMP file: {os.path.basename(temp_dmp_file)}")

        # 5. Call DMR analysis function
        # Note: run_dmr_pipeline_on_dmp_file requires the chromoNo parameter
        # pass the total number of chromosomes
        try:
            dmr_list_file = run_dmr_pipeline_on_dmp_file(
                dmp_file=temp_dmp_file,
                chromoNo=total_chromosomes
            )

            if dmr_list_file:
                dmr_results[chr_num] = dmr_list_file
                print(f"     Chromosome {chr_num} DMR analysis complete")
            else:
                print(f"    Chromosome {chr_num} DMR analysis failed (possibly no valid DMRs)")

        except Exception as e:
            print(f"    Error processing chromosome {chr_num}: {e}")
            import traceback
            traceback.print_exc()
            continue

    print(f"\nCommon significant site DMR analysis for {methylation_type} complete!")
    print(f"Successfully processed {len(dmr_results)}/{len(chromosomes)} chromosomes")

    return dmr_results

def summarize_dmr_methylation(methy_dir, replicate_x, replicate_y, file1_path, file2_path, methylation_type='CpG', custom_dmr_dir=None):
    """
    Sum methylation reads for DMR regions, and calculate Fisher's p-value, FDR q-value, and methylation direction.
    All output files will be saved in the methy_dir directory (e.g., ./output_1_1/CpG/).
    NOTE!!! custom_dmr_dir: Directory containing custom DMR files (if None, methy_dir is used). This distinguishes DMRs generated by each output_x_y from those generated by common_sites.
    """
    # Processing a single test here
    print(f"    Starting DMR methylation read summation analysis...")

    n_chromosomes = get_column_count(file1_path)
    if n_chromosomes is None:
        print("    Unable to get the number of chromosomes, skipping DMR summation")
        return

    chromosomes = [f'Chr{i}' for i in range(1, n_chromosomes + 1)]

    # Collect DMR data for all chromosomes (excluding p-values)
    all_dmr_data = []  # (chrom, start, end, exp_m, exp_u, wild_m, wild_u)

    for idx, chrom in enumerate(chromosomes):
        chrom_num = idx + 1
        if custom_dmr_dir is not None:
            # Read common DMR file
            dmr_file = os.path.join(custom_dmr_dir, f"DMR_list_DMP_common_{methylation_type}_Chr{chrom_num}.txt")
        else:
            # Read DMR file from output directory
            dmr_file = os.path.join(methy_dir, f"DMR_list_DMP_wt_replicate{replicate_y}_mut_replicate{replicate_x}_Chr{chrom_num}.txt")
        if not os.path.exists(dmr_file):
            print(f"      Skipping {chrom}, DMR file does not exist")
            continue

        # Read DMR regions
        dmr_list = []
        with open(dmr_file, 'r') as f:
            for line in f:
                parts = line.strip().split()
                if len(parts) >= 2:
                    try:
                        start = int(parts[0])
                        end = int(parts[1])
                        dmr_list.append((start, end))
                    except ValueError:
                        continue

        if not dmr_list:
            print(f"      {chrom} has no valid DMR regions")
            continue

        col_start = idx * 3

        # Read experimental group data
        exp_sites, exp_methy_dict, exp_unmethy_dict = [], {}, {}
        try:
            with open(file1_path, 'r') as f:
                for line in f:
                    parts = line.strip().split()
                    if len(parts) < col_start + 3:
                        continue
                    try:
                        site = int(float(parts[col_start]))
                        methy = int(float(parts[col_start + 1]))
                        unmethy = int(float(parts[col_start + 2]))
                        if site > 0:
                            exp_sites.append(site)
                            exp_methy_dict[site] = methy
                            exp_unmethy_dict[site] = unmethy
                    except:
                        continue
            exp_sites.sort()
        except Exception as e:
            print(f"      Failed to read experimental group ({chrom}): {e}")
            continue

        # Read control group data
        wild_sites, wild_methy_dict, wild_unmethy_dict = [], {}, {}
        try:
            with open(file2_path, 'r') as f:
                for line in f:
                    parts = line.strip().split()
                    if len(parts) < col_start + 3:
                        continue
                    try:
                        site = int(float(parts[col_start]))
                        methy = int(float(parts[col_start + 1]))
                        unmethy = int(float(parts[col_start + 2]))
                        if site > 0:
                            wild_sites.append(site)
                            wild_methy_dict[site] = methy
                            wild_unmethy_dict[site] = unmethy
                    except:
                        continue
            wild_sites.sort()
        except Exception as e:
            print(f"      Failed to read control group ({chrom}): {e}")
            continue

        # Summarize reads for each DMR region
        for start, end in dmr_list:
            l1 = bisect.bisect_left(exp_sites, start)
            r1 = bisect.bisect_right(exp_sites, end)
            exp_m_sum = sum(exp_methy_dict[site] for site in exp_sites[l1:r1])
            exp_u_sum = sum(exp_unmethy_dict[site] for site in exp_sites[l1:r1])

            l2 = bisect.bisect_left(wild_sites, start)
            r2 = bisect.bisect_right(wild_sites, end)
            wild_m_sum = sum(wild_methy_dict[site] for site in wild_sites[l2:r2])
            wild_u_sum = sum(wild_unmethy_dict[site] for site in wild_sites[l2:r2])

            all_dmr_data.append((chrom, start, end, exp_m_sum, exp_u_sum, wild_m_sum, wild_u_sum))

    if not all_dmr_data:
        print("    No valid DMR data found, skipping subsequent analysis")
        return

    # === Step 1: Output dmr_summary_{chrom}.txt ===
    chrom_summary_dict = defaultdict(list)
    for item in all_dmr_data:
        chrom, start, end, exp_m, exp_u, wild_m, wild_u = item
        chrom_summary_dict[chrom].append((start, end, exp_m, exp_u, wild_m, wild_u))

    for chrom in chromosomes:
        if chrom in chrom_summary_dict:
            summary_file = os.path.join(methy_dir, f"dmr_summary_{chrom}.txt")
            with open(summary_file, 'w') as f:
                f.write("DMR_start\tDMR_end\texp_methy_sum\texp_unmethy_sum\twild_methy_sum\twild_unmethy_sum\n")
                for row in chrom_summary_dict[chrom]:
                    f.write(f"{row[0]}\t{row[1]}\t{row[2]}\t{row[3]}\t{row[4]}\t{row[5]}\n")
            print(f"      {chrom} DMR summary complete → {summary_file}")

    # === Step 2: Calculate p-values for each DMR ===
    dmr_with_pvals = []
    for (chrom, start, end, exp_m, exp_u, wild_m, wild_u) in all_dmr_data:
        if (exp_m + exp_u < 4) or (wild_m + wild_u < 4):
            pval = np.nan
        else:
            try:
                _, pval = fisher_exact([[exp_m, exp_u], [wild_m, wild_u]], alternative='two-sided')
                pval = float(pval)
            except:
                pval = np.nan
        dmr_with_pvals.append((chrom, start, end, exp_m, exp_u, wild_m, wild_u, pval))

    # === Step 3: Global FDR correction (cross-chromosome) ===
    pvalues = np.array([item[7] for item in dmr_with_pvals])
    qvalues = calculate_qvalues(pvalues, pi=1.0)

    # === Step 4: Organize results by chromosome and determine direction ===
    chrom_data_dict = defaultdict(list)
    for i, (chrom, start, end, exp_m, exp_u, wild_m, wild_u, pval) in enumerate(dmr_with_pvals):
        qval = qvalues[i]

        # Calculate methylation direction
        exp_total = exp_m + exp_u
        wild_total = wild_m + wild_u
        if exp_total > 0 and wild_total > 0:
            exp_rate = exp_m / exp_total
            wild_rate = wild_m / wild_total
            direction = 1 if exp_rate > wild_rate else 0
        else:
            direction = 0  # Set to 0 if undeterminable

        chrom_data_dict[chrom].append((start, end, exp_m, exp_u, wild_m, wild_u, pval, qval, direction))

    # === Step 5: Output complete Fisher results and significant subsets ===
    for chrom in chromosomes:
        if chrom not in chrom_data_dict:
            continue

        # Complete results
        fisher_file = os.path.join(methy_dir, f"dmr_fisher_{chrom}.txt")
        with open(fisher_file, 'w') as f:
            f.write("DMR_start\tDMR_end\texp_methy_sum\texp_unmethy_sum\twild_methy_sum\twild_unmethy_sum\tpvalue\tqvalue\tdirection\n")
            for row in chrom_data_dict[chrom]:
                start, end, exp_m, exp_u, wild_m, wild_u, pval, qval, direction = row
                p_out = pval if not np.isnan(pval) else 'nan'
                q_out = qval if not np.isnan(qval) else 'nan'
                f.write(f"{start}\t{end}\t{exp_m}\t{exp_u}\t{wild_m}\t{wild_u}\t{p_out:.6g}\t{q_out:.6g}\t{direction}\n")
        print(f"      {chrom} Fisher + FDR + direction complete → {fisher_file}")

        # Significant results (q < 0.05)
        sig_rows = [
            row for row in chrom_data_dict[chrom]
            if not np.isnan(row[7]) and row[7] <= 0.05
        ]
        if sig_rows:
            sig_file = os.path.join(methy_dir, f"dmr_fisher_significant_{chrom}.txt")
            with open(sig_file, 'w') as f:
                f.write("DMR_start\tDMR_end\texp_methy_sum\texp_unmethy_sum\twild_methy_sum\twild_unmethy_sum\tpvalue\tqvalue\tdirection\n")
                for row in sig_rows:
                    start, end, exp_m, exp_u, wild_m, wild_u, pval, qval, direction = row
                    f.write(f"{start}\t{end}\t{exp_m}\t{exp_u}\t{wild_m}\t{wild_u}\t{pval:.6g}\t{qval:.6g}\t{direction}\n")
            print(f"        → Significant DMRs (q<0.05): {sig_file}")
        else:
            print(f"        → {chrom} has no significant DMRs (q<0.05)")

    print("    DMR methylation analysis complete.")

def run_dmr_pipeline_on_dmp_file(dmp_file: str, chromoNo: int = 10):
    """
    Generate DMR regions from the DMP file and output the dmr_list file
    - methylation_matrix_file: Used to dynamically obtain the number of chromosomes (e.g., 1-bothMeUnme_...txt)
    - All output files are saved in the directory of the dmp_file.
    """
    # The core concept is as follows:
    # Use a sliding window to identify DMP-dense regions
    # Connect adjacent dense regions via "jump merging"
    # Finally, filter out sufficiently long and significant DMRs
    sWinN = 1000  # Sliding window size（1000 bp）
    M0 = 4  # Minimum number of DMPs within a window (at least 4)
    M1 = 10  # Minimum number of DMPs in the final DMR (at least 10)
    M2 = 10  # Jump step size, retaining original segmentation behavior

    # For safety, ensure chromoNo >= 6 (since arrayMethy1_script1[5] is used)
    chromoNo = max(chromoNo, 6)

    class PositionNoNode:
        def __init__(self, pos=0, end=0, pV=0.0, ratio=0.0, num=0, num2=0, numCom=0, markR=0, DMR_S=0, DMR_E=0):
            self.pos = pos  # Window start position
            self.end = end  # Window end position
            self.pV = pV
            self.ratio = ratio
            self.num = num  # Number of hypermethylated sites
            self.num2 = num2  # Number of hypomethylated sites
            self.numCom = numCom
            self.markR = markR
            self.DMR_S = DMR_S
            self.DMR_E = DMR_E
            self.posV = []
            self.logPV = []
            self.meUnV = []

    output_dir = os.path.dirname(dmp_file)
    base_name = os.path.basename(dmp_file)

    arrayMethy1 = [[] for _ in range(chromoNo)]  # Create sublists matching the number of chromosomes

    # === Read DMP sites ===
    try:
        with open(dmp_file, 'r') as fin1:
            lines = fin1.readlines()
            if not lines or len(lines) < 2:
                return None
            for line in lines[1:]:
                parts = line.strip().split()
                if len(parts) < 3:
                    continue
                try:
                    pos = int(parts[0])
                    pValue = float(parts[1])
                    change = int(parts[2])
                    tmpNode = PositionNoNode()
                    tmpNode.pos = pos
                    tmpNode.pV = pValue
                    tmpNode.num = change  # The above lines establish the site information corresponding to the current DMP row
                    arrayMethy1[0].append(tmpNode)
                except (ValueError, IndexError):
                    continue
    except Exception as e:
        print(f"    Failed to read DMP file {dmp_file}: {e}")
        return None

    if not arrayMethy1[0]:
        return None

    arrayMethy1[0].sort(key=lambda x: x.pos)
    firstP = arrayMethy1[0][0].pos
    lastP = arrayMethy1[0][-1].pos

    # === Sliding window ===
    ratioI = 0.1
    cirN0 = (firstP - 1) // sWinN
    cir2 = 0
    while True:
        startP = cirN0 * sWinN + int(sWinN * ratioI * cir2) + 1
        endP = cirN0 * sWinN + int(sWinN * (ratioI * cir2 + 1))
        if startP > lastP:
            break

        numS1 = numS2 = 0
        for node in arrayMethy1[0]:
            if node.pos > endP:
                break
            if startP <= node.pos <= endP:
                if node.num == 1:
                    numS1 += 1
                else:
                    numS2 += 1

        tmpNode = PositionNoNode()
        tmpNode.pos = startP
        tmpNode.end = endP
        tmpNode.num = numS1
        tmpNode.num2 = numS2
        arrayMethy1[1].append(tmpNode)
        cir2 += 1

    # === Construct sliding window list ===
    arrayMethy1[2] = []
    for node in arrayMethy1[1]:
        new_node = PositionNoNode()
        new_node.pos = node.pos
        new_node.end = node.end
        new_node.num = node.num
        new_node.num2 = node.num2
        new_node.numCom = node.num + node.num2
        new_node.markR = 0
        arrayMethy1[2].append(new_node)

    # === Standardize output ===
    if arrayMethy1[2]:
        maxCom = max(node.numCom for node in arrayMethy1[2])
        maxCom = max(maxCom, 1)
        std_file = os.path.join(output_dir, f"noTitle_allDMCs_new_Standardized_slidingW_{base_name}")
        with open(std_file, 'w') as f:
            for node in arrayMethy1[2]:
                if node.end <= lastP:
                    std_val = node.numCom / maxCom
                    f.write(f"{node.pos:<20}{node.end:<20}{node.numCom:<20}{std_val:.6f}\n")

    # === Jump merging ===
    dmr_out_file = os.path.join(output_dir, f"DMR_{base_name}")
    boundary_file = os.path.join(output_dir, f"DMR_boundaries_{base_name}")

    arrayMethy1[3] = []
    for i, it02 in enumerate(arrayMethy1[2]):
        it03 = it02
        it04 = it02
        num1 = it02.pos
        num2 = it02.end
        num3 = it02.numCom

        num01 = num1
        num02 = num2
        num03 = num3
        num04 = 0

        flag1 = False
        flag2 = False
        flag3 = False
        flag4 = False
        flag5 = False
        flag6 = False
        flag7 = False

        # Expand to the left
        while num03 >= M0 and num01 >= firstP:
            if not flag3:
                flag1 = True
                flag3 = True
                num07 = num02
                num09 = num01
                num10 = num02

            num04 += num03
            num06 = num01
            it03.markR = 1

            idx = arrayMethy1[2].index(it03)
            if idx > 0:
                for k in range(M2):
                    idx -= 1
                    if idx <= 0:
                        flag5 = True
                        break
                    it03 = arrayMethy1[2][idx]
                if flag5:
                    break
                else:
                    num01 = it03.pos
                    num03 = it03.numCom
            else:
                break

        # Expand to the right
        if flag1:
            idx = arrayMethy1[2].index(it04)
            for k in range(M2):
                idx += 1
                if idx >= len(arrayMethy1[2]) - 1:
                    flag6 = True
                    break
                it04 = arrayMethy1[2][idx]

            if not flag6:
                while idx < len(arrayMethy1[2]):
                    it04 = arrayMethy1[2][idx]
                    num02 = it04.end
                    num03 = it04.numCom
                    # M0:4
                    if num03 >= M0 and num02 <= lastP:
                        if not flag4:
                            flag2 = True
                            flag4 = True

                        num07 = num02
                        num04 += num03
                        it04.markR = 1

                        for k in range(M2):
                            idx += 1
                            if idx >= len(arrayMethy1[2]) - 1:
                                flag7 = True
                                break
                            it04 = arrayMethy1[2][idx]

                        if flag7:
                            break
                    else:
                        break

        if flag1 or flag2:
            tmpNode = PositionNoNode()
            tmpNode.pos = num06
            tmpNode.end = num07
            tmpNode.numCom = num04
            tmpNode.DMR_S = num09
            tmpNode.DMR_E = num10
            arrayMethy1[3].append(tmpNode)

    print("=" * 60)
    print(f"Identified {len(arrayMethy1[3])} DMR regions")

    # Output DMR results and boundary files
    with open(dmr_out_file, 'w') as cout05:
        with open(boundary_file, 'w') as bound_out:
            for node in arrayMethy1[3]:
                cout05.write(f"{node.pos:<20}{node.end:<20}{node.numCom:<20}[{node.DMR_S} {node.DMR_E}]\n")
                bound_out.write(f"{node.DMR_S} {node.DMR_E}\n")

    print(f"Generated boundary file: {boundary_file}")

    # === Merge overlapping boundaries (dynamic chromoL) ===
    chromoL = lastP + 100000  # Dynamic genome length
    chromoArray = [0] * chromoL
    try:
        with open(boundary_file, 'r') as f_in:
            for line in f_in:
                parts = line.strip().split()
                if len(parts) >= 2:
                    try:
                        s, e = int(parts[0]), int(parts[1])
                        for idx in range(s, min(e + 1, chromoL)):
                            chromoArray[idx] = 1
                    except:
                        continue
    except:
        return boundary_file

    arrayMethy1_script1 = [[] for _ in range(chromoNo)]
    i = 0
    while i < chromoL:
        start = None
        for j in range(i, chromoL):
            if chromoArray[j] == 1:
                start = j
                break
        if start is None:
            break
        end = start
        while end + 1 < chromoL and chromoArray[end + 1] == 1:
            end += 1
        node = PositionNoNode()
        node.DMR_S = start
        node.DMR_E = end
        arrayMethy1_script1[0].append(node)
        i = end + 1

    no_overlap_file = os.path.join(output_dir, f"boundaries_noOverlapping_{base_name}")
    with open(no_overlap_file, 'w') as f_no:
        for node in arrayMethy1_script1[0]:
            f_no.write(f"{node.DMR_S}    {node.DMR_E}\n")

    arrayMethy1_script1[5] = arrayMethy1[0]
    final_dmr_list = []
    for region in arrayMethy1_script1[0]:
        s, e = region.DMR_S, region.DMR_E
        hyper = hypo = 0
        for site in arrayMethy1_script1[5]:
            if site.pos < s:
                continue
            if site.pos > e:
                break
            if site.num == 1:
                hyper += 1
            else:
                hypo += 1

        total = hyper + hypo
        length = e - s + 1
        if total >= M1 and length >= sWinN:
            node = PositionNoNode()
            node.pos = s
            node.end = e
            node.num = hyper
            node.num2 = hypo
            node.numCom = total
            final_dmr_list.append(node)

    final_file = os.path.join(output_dir, f"DMR_list_{base_name}")
    final_dmr_list.sort(key=lambda x: x.pos)
    with open(final_file, 'w') as f_final:
        for node in final_dmr_list:
            f_final.write(f"{node.pos:<20}{node.end:<20}{node.num:<20}{node.num2:<20}{node.numCom:<20}\n")

    print(f"    DMR analysis complete: {base_name} → {final_file}")
    return final_file


def process_chr_in_one_file(df):
    """Modify the prefix of a single input file to 'chr' and return all chromosome information for that file"""

    # This function ensures that chromosome names starting with 'chr' remain unchanged, while those without the prefix will have 'chr' added as the new chromosome name
    def to_chr(chr_str):
        chr_str = str(chr_str).strip()
        if chr_str.lower().startswith('chr'):
            return f"chr{chr_str[3:]}"
        else:
            return f"chr{chr_str}"

    # Change all chromosome names to start with 'chr'
    df['Chromosome'] = df['Chromosome'].apply(to_chr)
    # Return a unique set of all chromosome numbers
    return df['Chromosome'].unique()  # Return a numpy.array

def natural_sort_key(chr_name):
    """
    Custom chromosome sorting function
    Numeric chromosomes are sorted numerically; letter chromosomes (X, Y, M, etc.) are placed at the end and sorted alphabetically
    """
    chr_name = str(chr_name).lower()
    # Remove 'chr' prefix
    if chr_name.startswith('chr'):
        suffix = chr_name[3:]
    else:
        suffix = chr_name

    # Attempt to convert to integer
    try:
        # If numeric, return (0, numeric value, '')
        return (0, int(suffix), '')
    except ValueError:
        # If a letter (e.g., X, Y, M), return (1, 0, letter)
        return (1, 0, suffix)

def scan_all_files_for_chr_mapping(m, n, dir1, dir2):
    """Scan all m+n files, using the previous function to collect all chromosome information and generate a unified mapping"""

    all_chromosomes = set()  # Use the uniqueness property of sets to store all chromosome numbers

    # Scan m files in the first directory, named in the format i-dir1.txt (e.g., 3-msv.txt)
    for i in range(1, m + 1):
        filepath = os.path.join(dir1, f"{i}-{os.path.basename(dir1)}.txt")
        if not os.path.exists(filepath):
            print(f"Warning: File does not exist {filepath}")
            continue

        try:
            df = pd.read_csv(filepath,
                             sep='\s+',
                             header=None,
                             names=['Chromosome', 'Position', 'Methylated_reads', 'Unmethylated_reads', 'Methylation_type'],
                             dtype=str)
            chromosomes = process_chr_in_one_file(df)  # Get the unique set of all chromosome numbers in this file
            all_chromosomes.update(chromosomes)  # Add all unique chromosome numbers from this file to the all_chromosomes set
            # Note: update accepts any iterable and adds its elements to the set
            print(f"File {filepath} contains chromosomes: {sorted(chromosomes, key=natural_sort_key)}")
        except Exception as e:
            print(f"Error reading file {filepath} : {e}")

    # Scan files in the second directory
    for j in range(1, n + 1):
        filepath = os.path.join(dir2, f"{j}-{os.path.basename(dir2)}.txt")
        if not os.path.exists(filepath):
            print(f"Warning: File does not exist {filepath}")
            continue

        try:
            df = pd.read_csv(filepath,
                             sep='\s+',
                             header=None,
                             names=['Chromosome', 'Position', 'Methylated_reads', 'Unmethylated_reads', 'Methylation_type'],
                             dtype=str)
            chromosomes = process_chr_in_one_file(df)  # 同上
            all_chromosomes.update(chromosomes)  # 同上
            print(f"File {filepath} contains chromosomes: {sorted(chromosomes, key=natural_sort_key)}")
        except Exception as e:
            print(f"Error reading file {filepath} : {e}")

    # Create a unified chromosome mapping,
    # Note: The all_chromosomes set now contains unique chromosome numbers from all files in both genotype folders
    unique_chrs = sorted(all_chromosomes, key=natural_sort_key)
    # Build a Series (or array) mapping unique, sorted chromosome numbers from both genotype directories to numeric values
    chr_series = pd.Series(range(len(unique_chrs)), index=unique_chrs)

    print(f"Unified chromosome mapping: {chr_series}")
    return chr_series

def single_newtoboth(filepath1, output_dir, num1, chr_series):
    '''Parameters here: filepath1 (e.g., 1-wt.txt) is the new format file being processed,
    output_dir is the directory where the resulting 'both' file is outputted, usually the same directory as filepath1,
    num1 indicates which new format file is currently being processed (i.e., currently processing num1-genotype.txt),
    chr_series is a sorted Series mapping all unique chromosome numbers to numeric values'''

    df = pd.read_csv(filepath1,
                     sep='\s+',
                     header=None,
                     names=['Chromosome', 'Position', 'Methylated_reads', 'Unmethylated_reads', 'Methylation_type'],
                     dtype=str)

    def to_chr(chr_str):
        chr_str = str(chr_str).strip()
        if chr_str.lower().startswith('chr'):
            return f"chr{chr_str[3:]}"
        else:
            return f"chr{chr_str}"

    # Since the source file was not modified when building chr_series earlier, the imported chromosome numbers might still lack the 'chr' prefix and need to be processed
    df['Chromosome'] = df['Chromosome'].apply(to_chr)
    df['Chromosome'] = df['Chromosome'].map(chr_series)  # Modify each chromosome number to its corresponding numeric value, e.g., chr1->0, chr2->1
    # The following three lines convert the position and the two read counts to numeric values
    df['Position'] = pd.to_numeric(df['Position'])
    df['Methylated_reads'] = pd.to_numeric(df['Methylated_reads'])
    df['Unmethylated_reads'] = pd.to_numeric(df['Unmethylated_reads'])
    # Replace chromosome numbers in the df with corresponding indices from chr_series, convert positions and read counts to numeric types, and change methylation type 'CG' to 'CpG'
    df['Methylation_type'] = df['Methylation_type'].str.replace('CG', 'CpG')
    # Group the three types of data for subsequent separate operations
    data_groups = df.groupby('Methylation_type')
    # Iterate based on the grouping keys and their corresponding sub-dataframes
    for methy_type, data_ind in data_groups:
        if data_ind.empty:
            continue
        # Get and sort the actually existing chromosome numbers for the current methylation type,
        # These are the chromosome numbers present in the current methylation type. NOTE: !!! These are the converted zero-based indices !!!
        actual_chrs = sorted(data_ind['Chromosome'].dropna().unique())
        # Simultaneously get the total number of chromosomes, regardless of whether the methylation type has corresponding chromosome information
        chr_count = len(chr_series)
        chr_data_dict = {}
        # mlen is used to record the maximum number of chromosome data entries among all chromosomes for the current methylation type
        mlen = 0

        for chr_num in actual_chrs:  # Iterate through each chromosome number (zero-based value) for this methylation type
            # Sort all data for the current chromosome number in ascending order by position, and store it in chr_data
            chr_data = data_ind[data_ind['Chromosome'] == chr_num].sort_values('Position').reset_index(drop=True)
            # Store the three parameters of interest (Position, Methylated_reads, Unmethylated_reads) in the chr_data_dict with the chromosome number as the key
            # The value corresponding to each key is a numpy array composed of these three attributes
            chr_data_dict[chr_num] = chr_data[['Position', 'Methylated_reads', 'Unmethylated_reads']].values
            # Update mlen with any potentially larger chromosome data entry count
            mlen = max(mlen, len(chr_data_dict[chr_num]))

        # Create an output matrix, with rows equal to the maximum chromosome data entry count, and columns equal to the total number of unique chromosomes in both directories * 3, initially filled with 0
        output_matrix = np.zeros((mlen, chr_count * 3), dtype=np.int32)

        # Iterate through all rows in the output matrix
        for i in range(mlen):
            # Iterate through the chromosome numbers present in the current methylation type (zero-based numerical values)
            for chr_num in actual_chrs:  # Use continuous index
                col_start = chr_num * 3  # Calculate column position based on continuous index
                if i < len(chr_data_dict[chr_num]):
                    output_matrix[i, col_start:col_start + 3] = chr_data_dict[chr_num][i]

        # Output(export using pandas)
        output_df = pd.DataFrame(output_matrix)  # Convert the output matrix to a DataFrame for easier output
        output_file = f"{num1}-bothMeUnme_diffChromo_NOREPEATED_methy_sites_{methy_type}.txt"
        output_path = os.path.join(output_dir, output_file)
        output_df.to_csv(output_path, sep='\t', header=False, index=False)

def get_chr_name(chr_num, chr_series):
    """
    Get chromosome name based on chromosome index
    Parameters：
        chr_num: The nth chromosome (1-based)
        chr_series: Chromosome mapping Series
    Returns：
        Chromosome name string (e.g., 'chr1', 'chrX')
    """
    chr_num = int(chr_num)
    if chr_series is not None:
        try:
            # chr_num is 1-based, so subtract 1
            if 0 <= chr_num - 1 < len(chr_series):
                return chr_series.index[chr_num - 1]
        except Exception as e:
            print(f"Warning: Failed to get chromosome name: {e}")

    # Return numeric ID if an error occurs or no mapping exists
    return f"chr{chr_num}"

def newtoboth(m, n, dir1, dir2):
    # No need to check if directories exist, as the main function logic has already verified this
    # Obtain a sorted Series mapping unique chromosome numbers to numeric values across both genotype directories
    chr_series = scan_all_files_for_chr_mapping(m, n, dir1, dir2)
    print(f"Mapping relationship: {chr_series}")
    # Loop m+n times to process files for both genotypes respectively
    for i in range(1, m + 1):
        filepath = os.path.join(dir1, f"{i}-{os.path.basename(dir1)}.txt")
        if not os.path.exists(filepath):
            print(f"Warning: File does not exist {filepath}")
            continue
        print(f"Processing file {filepath}")
        # Convert i-dir1.txt in the dir1 folder to the 'both' format
        single_newtoboth(filepath, dir1, i, chr_series)
    for j in range(1, n + 1):
        filepath = os.path.join(dir2, f"{j}-{os.path.basename(dir2)}.txt")
        if not os.path.exists(filepath):
            print(f"Warning: File does not exist {filepath}")
            continue
        print(f"Processing file{filepath}")
        # Convert j-dir2.txt in the dir2 folder to the 'both' format
        single_newtoboth(filepath, dir2, j, chr_series)
    return chr_series

def sanitize_filename(name):
    """Sanitize filenames by removing special characters not allowed in file names"""
    return re.sub(r'[\\/*?:"<>|]', "", name)  # Replace special characters in the name with an empty string to sanitize it


def get_column_count(file_path):
    """Get the number of columns in the file and return the column count divided by 3"""
    try:
        with open(file_path, 'r') as file:
            first_line = file.readline().strip()
            column_count = len(first_line.split())
            return column_count // 3
    except Exception as e:
        print(f"Error occurred while reading file: {e}")
        return None


def parse_filename(filename):
    """Parse a single filename, using regular expressions to extract the file ID and methylation type"""
    # Use capture groups (the parts in parentheses) to capture the corresponding file ID and methylation type
    pattern = r'^(\d+)-bothMeUnme_diffChromo_NOREPEATED_methy_sites_(.+)\.txt$'
    match = re.match(pattern, filename)  # Match the specific filename against this regular expression
    if match:  # Upon matching, use match.group(i) to get the content of the i-th capture group
        file_id = int(match.group(1))
        methylation_type = match.group(2)
        return file_id, methylation_type
    return None


def scan_sample_files_by_replicates(sample_dir, max_replicates):
    """Search for files with IDs 1~max_replicates ending with '-both...methylation_type.txt' in the sample_dir directory"""
    files_by_replicates = {}  # Create a dictionary to store corresponding filenames based on the prefix ID
    methylation_types = ['CpG', 'CHH', 'CHG']
    # Return an empty dictionary if the directory does not exist; this step can be omitted after merging since it is checked initially in the main function
    # if not os.path.exists(sample_dir):
    #     return files_by_replicates

    # Iterate through all contents (including files and subdirectories) in the folder
    for filename in os.listdir(sample_dir):
        if filename.endswith('.txt'):  # If a txt file is found (it can only be a 'both' file or a new format file)
            parsed = parse_filename(filename)  # Parse the file to obtain the potentially existing ID and methylation type
            if parsed:  # If the ID and methylation type are successfully parsed for the file
                file_id, methylation_type = parsed  # (Tuple; the comma indicates a tuple, while parentheses merely avoid ambiguity)
                if file_id <= max_replicates and methylation_type in methylation_types:  # Both ID and methylation type are valid
                    if file_id not in files_by_replicates:
                        files_by_replicates[file_id] = {}  # Make the files_by_replicates dictionary map file_id -> {sub-dictionary}
                    files_by_replicates[file_id][methylation_type] = filename  # Format the sub-dictionary as methylation_type -> filename
    return files_by_replicates  # Returns a mapping of ID -> methylation_type -> filename, allowing filename retrieval via ID and methylation type


def process_methylation_type_with_collection(file1_path, file2_path, methylation_type, output_dir,
                                             dir1_name, dir2_name, chr_num):
    """Process Fisher's exact test for a methylation type and chromosome (processes one chromosome for one test)
    Parameter explanation: file_path - Relative path of the current 'both' file being processed
    methylation_type - The methylation type currently being processed
    output_dir - Output folder ./output_x_y/
    dir_name - The final part of the two input genotype directory paths
    chr_num - The nth chromosome being processed, not the chromosome name/number"""

    print(f"      Processing methylation type {methylation_type}, chromosome {chr_num}...")

    # Calculate data columns for the current chromosome
    mOrder = 3 * (chr_num - 1)  # Here chr_num refers to the nth chromosome, not the chromosome name/number!!!
    usecols = [mOrder, mOrder + 1, mOrder + 2]

    try:
        # Read all data from the first file into a dictionary in chunks
        data1_dict = {}
        for chunk in pd.read_csv(file1_path, sep='\s+', header=None, usecols=usecols, chunksize=100000):
            chunk.columns = ['pos', 'methy', 'unmethy']
            # Use zip to get one element from each column at a time, placing the three elements into a (pos, methy, unmethy) tuple
            for pos, methy, unmethy in zip(chunk['pos'].values, chunk['methy'].values, chunk['unmethy'].values):
                # Place data into the corresponding dictionary in the following format
                data1_dict[pos] = (methy, unmethy)
        # At this point, all site data from the first file has been loaded into the data1_dict dictionary

        # Read the second file in chunks to find common sites; loading all data from both files into memory simultaneously might consume too much memory
        # The first file must be fully loaded into memory because we need fast random access to confirm if a site exists in both files; if it does, it must be tested
        # Because 'reader' uses the 'chunksize' parameter, the return value of read_csv is an iterator, yielding up to 100,000 rows of data per iteration
        reader2 = pd.read_csv(file2_path, sep='\s+', header=None, usecols=usecols, chunksize=100000)

    except Exception as e:
        print(f"        Failed to load data: {e}")
        return False

    # Create output folder: output_x_y/methylation_type/
    methylation_dir = os.path.join(output_dir, methylation_type)
    os.makedirs(methylation_dir, exist_ok=True)

    # Create output file paths
    stats_filename = os.path.join(methylation_dir,
                                  f"FET_results_{methylation_type}_{dir1_name}_{dir2_name}_Chr{chr_num}.txt")
    all_output = os.path.join(methylation_dir, f"all_simple_Chr{chr_num}.txt")
    sig_output = os.path.join(methylation_dir, f"sig_simple_Chr{chr_num}.txt")
    combine_output = os.path.join(methylation_dir, f"combineResult_Chr{chr_num}.txt")

    # Set parameters
    M0, M2, Th_pValue = 2, 4, 0.05

    try:
        # Create four lists to store various data for subsequent output
        all_results, sig_results, fet_results, combine_results = [], [], [], []

        # Process each data chunk from the second file
        for chunk2 in reader2:
            chunk2.columns = ['pos', 'methy', 'unmethy']
            # For each data chunk of the second file, iterate through the three values in each row and store them in pos, m2, u2
            for pos, m2, u2 in zip(chunk2['pos'].values, chunk2['methy'].values, chunk2['unmethy'].values):
                # If the site also exists in the first file, perform Fisher's exact test for this site
                if pos in data1_dict:
                    m1, u1 = data1_dict[pos]  # Get the two read counts from the first file

                    if m1 >= M0 or m2 >= M0:  # Both methylated read counts must be >= 2 to proceed
                        if (m1 + u1 >= M2) and (m2 + u2 >= M2):  # The sum of the two read counts in each file must be >= 4 to perform the test
                            cont_table = np.array([[m1, u1], [m2, u2]])  # Construct a 2x2 contingency table
                            # Calculate the two methylation ratios
                            ratio1 = m1 / (m1 + u1) if (m1 + u1) > 0 else 0
                            ratio2 = m2 / (m2 + u2) if (m2 + u2) > 0 else 0
                            change = "1" if ratio1 >= ratio2 else "0"  # Determine based on the methylation ratios of the two files
                            # Whether the methylation ratio of the mutant increased
                            # Call library function to perform Fisher's exact test and obtain the p-value
                            _, pvalue = fisher_exact(cont_table, alternative='two-sided')
                            # Retain 7 significant digits for the p-value
                            pvalue = float(f"{pvalue:.7g}")
                            # Enter the required data into four lists respectively; only significant ones are entered into sig_results
                            all_results.append([pos, pvalue, change])
                            fet_results.append([pos, m1, u1, m2, u2, pvalue])
                            combine_results.append([pos, m1, u1, m2, u2])

                            if pvalue < Th_pValue:
                                sig_results.append([pos, pvalue, change])
        # At this point, all required Fisher's exact tests for this set of files (specific x_y, methylation type, and chromosome) are complete
        # Save results to disk
        if all_results:
            pd.DataFrame(all_results, columns=["Position", "Pvalue", "Methylation_Change"]).to_csv(all_output, sep='\t',
                                                                                                   index=False)
            pd.DataFrame(sig_results, columns=["Position", "Pvalue", "Methylation_Change"]).to_csv(sig_output, sep='\t',
                                                                                                   index=False)
            pd.DataFrame(fet_results,
                         columns=["Position", "Methy1", "Unmethy1", "Methy2", "Unmethy2", "Pvalue"]).to_csv(
                stats_filename, sep='\t', index=False)
            pd.DataFrame(combine_results, columns=["Position", "Methy1", "Unmethy1", "Methy2", "Unmethy2"]).to_csv(
                combine_output, sep='\t', index=False)

        print(f"        Chromosome {chr_num} processing complete! Processed a total of {len(all_results)} sites, of which {len(sig_results)} are significant")
        return True

    except Exception as e:
        print(f"        Error occurred while processing methylation type {methylation_type}, chromosome {chr_num} : {e}")
        return False

def process_methylation_type_with_collection_pvfilter(file1_path, file2_path, methylation_type, output_dir,
                                             dir1_name, dir2_name, chr_num):
    """Process Fisher's exact test for a methylation type and chromosome (processes one chromosome for one test)
    Parameter explanation: file_path - Relative path of the current 'both' file being processed
    methylation_type - The methylation type currently being processed
    output_dir - Output folder ./output_x_y/
    dir_name - The final part of the two input genotype directory paths
    chr_num - The nth chromosome being processed, not the chromosome name/number

    Return the dataframe containing p-values > 0.05 for this iteration
        all_results.append([pos, pvalue, change])
    """

    print(f"      Processing methylation type {methylation_type}, chromosome {chr_num}...")

    # Calculate data columns for the current chromosome
    mOrder = 3 * (chr_num - 1)  # Here chr_num refers to the nth chromosome, not the chromosome name/number!!!
    usecols = [mOrder, mOrder + 1, mOrder + 2]

    try:
        # Read all data from the first file into a dictionary in chunks
        data1_dict = {}
        for chunk in pd.read_csv(file1_path, sep='\s+', header=None, usecols=usecols, chunksize=100000):
            chunk.columns = ['pos', 'methy', 'unmethy']
            # Use zip to get one element from each column at a time, placing the three elements into a (pos, methy, unmethy) tuple
            for pos, methy, unmethy in zip(chunk['pos'].values, chunk['methy'].values, chunk['unmethy'].values):
                # Place data into the corresponding dictionary in the following format
                data1_dict[pos] = (methy, unmethy)
        # At this point, all site data from the first file has been loaded into the data1_dict dictionary

        # Read the second file in chunks to find common sites; loading all data from both files into memory simultaneously might consume too much memory
        # The first file must be fully loaded into memory because we need fast random access to confirm if a site exists in both files; if it does, it must be tested
        # Because 'reader' uses the 'chunksize' parameter, the return value of read_csv is an iterator, yielding up to 100,000 rows of data per iteration
        reader2 = pd.read_csv(file2_path, sep='\s+', header=None, usecols=usecols, chunksize=100000)

    except Exception as e:
        print(f"        Failed to load data: {e}")
        return False

    # Create output folder: output_x_y/methylation_type/
    methylation_dir = os.path.join(output_dir, methylation_type)
    os.makedirs(methylation_dir, exist_ok=True)

    # Create output file paths
    stats_filename = os.path.join(methylation_dir,
                                  f"FET_results_{methylation_type}_{dir1_name}_{dir2_name}_Chr{chr_num}.txt")
    all_output = os.path.join(methylation_dir, f"all_simple_Chr{chr_num}.txt")
    sig_output = os.path.join(methylation_dir, f"sig_simple_Chr{chr_num}.txt")
    combine_output = os.path.join(methylation_dir, f"combineResult_Chr{chr_num}.txt")

    # Set parameters
    M0, M2, Th_pValue = 2, 4, 0.05

    try:
        # Create four lists to store various data for subsequent output
        all_results, sig_results, fet_results, combine_results = [], [], [], []

        # Process each data chunk from the second file
        for chunk2 in reader2:
            chunk2.columns = ['pos', 'methy', 'unmethy']
            # For each data chunk of the second file, iterate through the three values in each row and store them in pos, m2, u2
            for pos, m2, u2 in zip(chunk2['pos'].values, chunk2['methy'].values, chunk2['unmethy'].values):
                # If the site also exists in the first file, perform Fisher's exact test for this site
                if pos in data1_dict:
                    m1, u1 = data1_dict[pos]  # Get the two read counts from the first file

                    if m1 >= M0 or m2 >= M0:  # Both methylated read counts must be >= 2 to proceed
                        if (m1 + u1 >= M2) and (m2 + u2 >= M2):  # The sum of the two read counts in each file must be >= 4 to perform the test
                            cont_table = np.array([[m1, u1], [m2, u2]])  # Construct a 2x2 contingency table
                            # Calculate the two methylation ratios
                            ratio1 = m1 / (m1 + u1) if (m1 + u1) > 0 else 0
                            ratio2 = m2 / (m2 + u2) if (m2 + u2) > 0 else 0
                            change = "1" if ratio1 > ratio2 else "0"  # Determine based on the methylation ratios of the two files
                            # Whether the methylation ratio of the mutant increased
                            # Call library function to perform Fisher's exact test and obtain the p-value
                            _, pvalue = fisher_exact(cont_table, alternative='two-sided')
                            # Retain 7 significant digits for the p-value
                            pvalue = float(f"{pvalue:.7g}")
                            # Enter the required data into four lists respectively; only significant ones are entered into sig_results
                            all_results.append([pos, pvalue, change])
                            fet_results.append([pos, m1, u1, m2, u2, pvalue])
                            combine_results.append([pos, m1, u1, m2, u2])

                            if pvalue < Th_pValue:
                                sig_results.append([pos, pvalue, change])
        # At this point, all required Fisher's exact tests for this set of files (specific x_y, methylation type, and chromosome) are complete
        # Save results to disk
        if all_results:
            # Convert to DataFrame for easier filtering
            all_df = pd.DataFrame(all_results, columns=["Position", "Pvalue", "Methylation_Change"])
            sig_df = pd.DataFrame(sig_results, columns=["Position", "Pvalue", "Methylation_Change"])
            fet_df = pd.DataFrame(fet_results,
                                  columns=["Position", "Methy1", "Unmethy1", "Methy2", "Unmethy2", "Pvalue"])
            combine_df = pd.DataFrame(combine_results, columns=["Position", "Methy1", "Unmethy1", "Methy2", "Unmethy2"])

            all_df_ndmp = all_df[all_df['Pvalue'] > 0.05]  # Retain information for p-value > 0.05, to be added to the q-value list later

            # Only keep sites with p <= 0.05
            all_df = all_df[all_df['Pvalue'] <= 0.05]  # Filter for p-value <= 0.05 to perform FDR correction
            fet_df = fet_df[fet_df['Pvalue'] <= 0.05]
            #fet_df_ndmp = all_df[all_df['Pvalue'] > 0.05]
            # sig_df is already p < 0.05, no further filtering needed

            # Also filter combine
            positions_to_keep = set(all_df['Position'].values)
            combine_df = combine_df[combine_df['Position'].isin(positions_to_keep)]  # Filter for p-value <= 0.05 to perform FDR correction
            #combine_df_ndmp = combine_df[~combine_df['Position'].isin(positions_to_keep)]
            # Save the filtered results
            all_df.to_csv(all_output, sep='\t', index=False)
            sig_df.to_csv(sig_output, sep='\t', index=False)
            fet_df.to_csv(stats_filename, sep='\t', index=False)
            combine_df.to_csv(combine_output, sep='\t', index=False)

        print(
            f"        Chromosome {chr_num} processing complete! Original tests: {len(all_results)} sites; after filtering(p≤0.05) {len(all_df)} sites;p>0.05:{len(all_df_ndmp)}sites")
        return all_df_ndmp
    #                   all_results.append([pos, pvalue, change])
    #                   fet_results.append([pos, m1, u1, m2, u2, pvalue])
    #               combine_results.append([pos, m1, u1, m2, u2])

    except Exception as e:
        print(f"        Error occurred while processing methylation type {methylation_type}, chromosome {chr_num} : {e}")
        return False

def merge_fet_results_and_fdr(output_dir, replicate_x, replicate_y, mtype3, all_dfs_ndmp_dict,n_chromosomes):
    """Merge all FET results in the output directory and perform FDR correction
    The format of the FET file is: pos, m1, u1, m2, u2, pvalue
    The output_dir input here is output_x_y/methylation_type
    success_dfs_dict[methylation_type][chr_num] accesses
                    the dataframe with p-value > 0.05 for the corresponding iteration
                      all_results([pos, pvalue, change])
    n_chromosomes is the total number of chromosomes
                      """
    print(f"\n    Merging FET results in {output_dir} and performing FDR correction...")

    if not os.path.exists(output_dir):
        print(f"    Error: Directory {output_dir} does not exist!")
        return False

    # Search for all FET result files
    # Here ** represents subdirectories of any depth, so when glob searches recursively (recursive=True), it will look in the output_dir and all its
    # subdirectories for files meeting the conditions, and return their paths (relative to output_dir) as a list
    file_pattern = os.path.join(output_dir, "**", "FET_results_*_Chr*.txt") # The chromosome number here is actually the nth set of three columns in the 'both' file
    fet_files = glob.glob(file_pattern, recursive=True)

    if not fet_files:
        print(f"    Warning: No FET result files found in {output_dir} ")
        return False

    print(f"   Found {len(fet_files)}FET result files")

    # Create a list to collect all p-values and related information
    all_data = []

    for file_path in sorted(fet_files):
        # Extract methylation type and chromosome information
        # Here, the first capture group captures the methylation type, . matches zero or more arbitrary characters for replicatex_replicatey, and the second capture group captures the chromosome number
        #                                                   # The chromosome number here is actually the nth set of three columns in the 'both' file
        methy_match = re.search(r'/FET_results_([^_]+)_.*_Chr(\d+)\.txt$', file_path.replace('\\', '/'))
        if not methy_match:
            continue

        # Obtained the methylation type and chromosome index (the chromosome index here is the nth set of three columns in the 'both' file)
        methylation_type = mtype3
        chr_num = int(methy_match.group(2))

        try:
            df = pd.read_csv(file_path, sep='\t', header=0)
            if 'Position' in df.columns and 'Pvalue' in df.columns:
                # Adjust column order: Chromosome, Methylation_Type, Position, Pvalue
                df_subset = df[['Position', 'Pvalue']].copy()
                df_subset['Chromosome'] = chr_num
                df_subset['Methylation_Type'] = methylation_type

                # Rearrange column order
                df_subset = df_subset[['Chromosome', 'Methylation_Type', 'Position', 'Pvalue']]
                all_data.append(df_subset) # Add all current FET information to all_data
        except Exception as e:
            print(f"    Warning: Failed to read {file_path} : {e}")
            continue

    if not all_data:
        print(f"    Error: Failed to read any data successfully")
        return False

    # Merge all data (since the dataframes corresponding to each file were appended directly to all_data, each element in the list is a dataframe, so concatenation is required)
    combined_df = pd.concat(all_data, ignore_index=True)
    # Sort the data by chromosome index (the chromosome index here is actually the nth set of three columns in the 'both' file)
    combined_df = combined_df.sort_values(['Methylation_Type', 'Chromosome', 'Position'])
    # The format of this dataframe is: 'Chromosome', 'Methylation_Type', 'Position', 'Pvalue'
    print(f"    Merged a total of {len(combined_df)} sites")

    # Calculate FDR-corrected q-values and append them as a new column to combined_df
    pvalues = combined_df['Pvalue'].values
    qvalues = calculate_qvalues(pvalues, 1.0)
    combined_df['Qvalue'] = qvalues

    # Final column order: Chromosome, Methylation_Type, Position, Pvalue, Qvalue
    combined_df = combined_df[['Chromosome', 'Methylation_Type', 'Position', 'Pvalue', 'Qvalue']]

    print(f"Starting to process ndmp data for {mtype3} ")

    # First check if there is data for this methylation type in all_dfs_ndmp_dict, as it's possible this methylation type for this species doesn't require calling the filter version
    if mtype3 in all_dfs_ndmp_dict and all_dfs_ndmp_dict[mtype3]:
        print(f"In this methylation type, the number of ndmp dataframes is{len(all_dfs_ndmp_dict[mtype3])}")
        dfs_ndmp = []
        for chr_num11 in range(1, n_chromosomes + 1):
            # Only process chromosomes present in this methylation type
            if chr_num11 in all_dfs_ndmp_dict[mtype3]:
                df_ndmp = all_dfs_ndmp_dict[mtype3][chr_num11]

                # Ensure it is a valid DataFrame
                if isinstance(df_ndmp, pd.DataFrame) and not df_ndmp.empty:
                    df_ndmp = df_ndmp.copy()

                    df_ndmp['Chromosome'] = chr_num11
                    df_ndmp['Methylation_Type'] = mtype3

                    if 'change' in df_ndmp.columns:
                        df_ndmp.drop(columns=['change'], inplace=True)

                    df_ndmp['Qvalue'] = 1
                    df_ndmp = df_ndmp[['Chromosome', 'Methylation_Type', 'Position', 'Pvalue', 'Qvalue']]

                    dfs_ndmp.append(df_ndmp)
                    print(f" Added required columns for {mtype3}--{chr_num11}, total {len(df_ndmp)} sites")

        if dfs_ndmp:
            total_ndmp = pd.concat(dfs_ndmp, ignore_index=True)
            combined_df = pd.concat([combined_df, total_ndmp], ignore_index=True)
            print(f"  Merged a total of {len(total_ndmp)} ndmp sites")
        else:
            print(f"  No ndmp data found for {mtype3} type")
    else:
        print(f" {mtype3} type does not exist in the ndmp dictionary")
    # At this point, for methylation types with larger data volumes, the information for sites with p-value > 0.05 previously discarded has been added,
    #   which had to be stored in all_dfs_ndmp_dict[methylation][chromosome_number] earlier to achieve this

    # Tally significant results, calculating how many significant p-values and q-values there are in this 1 out of mn3 tests
    n_pval_sig = np.sum(pvalues <= 0.05)
    dmp_threshold = get_dmp_threshold(mtype3)
    n_qval_sig = np.sum(qvalues <= dmp_threshold)
    # Calculate the proportion of significant sites
    print(f"   P-value significant sites: {n_pval_sig} ({n_pval_sig / len(pvalues) * 100:.1f}%)")
    print(f"   Q-value significant sites: {n_qval_sig} ({n_qval_sig / len(qvalues) * 100:.1f}%)")

    # Save the merged p-value list (for external FDR tools)
    pvalue_file = os.path.join(output_dir, f"united_pvalues_wt_replicate{replicate_y}_vs_mut_replicate{replicate_x}.csv")
    with open(pvalue_file, 'w') as f:
        for pvalue in pvalues:
            f.write(f"{pvalue}\n")

    # Save complete FDR correction results (output_dir is output_x_y/methylation_type)
    #  Includes complete p-values and q-values (format: Chromosome, Methylation_Type, Position, Pvalue, Qvalue)
    fdr_file = os.path.join(output_dir, f"FDR_corrected_results_wt_replicate{replicate_y}_vs_mut_replicate{replicate_x}.txt")
    combined_df.to_csv(fdr_file, sep='\t', index=False)

    # Save significant sites (q-value < 0.05)
    # Save significant sites (using different thresholds based on methylation type)
    dmp_threshold = get_dmp_threshold(mtype3)
    sig_df = combined_df[combined_df['Qvalue'] <= dmp_threshold]
    if not sig_df.empty: # If there are significant sites, output the significant portion of the data (output_dir is output_x_y/methylation_type)
        sig_file = os.path.join(output_dir, f"FDR_significant_results_wt_replicate{replicate_y}_vs_mut_replicate{replicate_x}.txt")
        sig_df.to_csv(sig_file, sep='\t', index=False)
        print(f"  Significant site results saved to: {sig_file}")

    print(f" P-value list saved to: {pvalue_file}")
    print(f"  FDR results saved to: {fdr_file}")
    return True


def calculate_qvalues(pvalues, pi=1.0):
    """Calculate Q-values using the Storey method"""

    pvalues = np.array(pvalues, dtype=float)
    if len(pvalues) == 0:
        return np.array([])

    pvalues = np.array(pvalues)
    pvalues_clean = pvalues.copy()
    pvalues_clean[np.isnan(pvalues_clean)] = 1.0
    pvalues_clean = np.clip(pvalues_clean, 0, 1)

    m = len(pvalues_clean)
    sorted_indices = np.argsort(pvalues_clean)
    sorted_pvalues = pvalues_clean[sorted_indices]

    # Key to the Storey method: incorporating the π factor
    q_values_sorted = np.zeros_like(sorted_pvalues)

    # If automatic estimation of π is needed
    if pi is None:
        # Simplified estimation of π
        lrange = np.linspace(0.05, 0.95, max(int(m / 100.0), 10))
        pil = np.mean(sorted_pvalues[:, np.newaxis] > lrange, axis=0)
        pilr = pil / (1.0 - lrange)
        pi = 1.0
        if pilr[-1] < 1.0:
            pi = pilr[-1]

    # Storey method calculation
    q_values_sorted = pi * m * sorted_pvalues / np.arange(1, m + 1)
    q_values_sorted[-1] = min(q_values_sorted[-1], 1.0)

    # Monotonicity adjustment
    for i in range(m - 2, -1, -1):
        q_values_sorted[i] = min(q_values_sorted[i], q_values_sorted[i + 1])

    # Restore original order
    q_values = np.zeros_like(pvalues_clean)
    q_values[sorted_indices] = q_values_sorted
    q_values[np.isnan(pvalues)] = np.nan

    return q_values

def perform_sliding_window_on_dmp_files(output_dir, replicate_x, replicate_y,):
    """Perform sliding window analysis on the (N)DMP file resulting from 1 of the mn3 processing iterations"""

    print(f"\n    Starting sliding window analysis on DMP files...")

    # Process only DMP files
    dmp_patterns = [
        f"DMP_wt_replicate{replicate_y}_mut_replicate{replicate_x}_Chr*.txt",
        f"N-DMP_wt_replicate{replicate_y}_mut_replicate{replicate_x}_Chr*.txt"
    ]

    for pattern in dmp_patterns:
        dmp_files = glob.glob(os.path.join(output_dir, pattern)) #Find the paths of all files matching the current regular expression format and put them into the dmp_files list

        for dmp_file in dmp_files: # Iterate through the list, obtaining one DMP or N-DMP file at a time
            try:

                df = pd.read_csv(dmp_file, sep='\s+', header=None, skiprows=1,
                                 names=['position', 'pvalue', 'change'])

                # Check if data is not empty
                if df.empty or len(df) == 0:
                    continue

                # Data type conversion and cleanup
                df = df.dropna()
                df['position'] = df['position'].astype(int)
                df['change'] = df['change'].astype(int)

                print(f"      Processing file: {os.path.basename(dmp_file)} ({len(df)} sites)")

                # Set the output prefix, which essentially removes '.txt' to retain the (N)DMP_replicate_wt{replicate_y}_mut_replicate{replicate_x}_Chr* portion
                base_name = os.path.basename(dmp_file).replace('.txt', '')

                # Pass the DataFrame to sliding_window_analysis
                # The former format is window_start, window_end, count_change_1, count_change_0, total_count
                # The latter lacks count_change but includes standardized_count, which is the ratio of significant sites in each interval to the maximum number of significant sites across all intervals
                sliding_results, std_results = sliding_window_analysis(
                    df, # df represents the DMP file for a specific chromosome from 1 of the mn3 processing iterations
                    window_size=1000000,
                    step_ratio=0.05,
                    save_files=True,
                    output_identifier=base_name,
                    outputdir1=output_dir
                )

                print(f"      Completed {base_name}: {len(sliding_results)} windows")

            except Exception as e:
                print(f"      Failed to process {dmp_file} : {e}")
                continue

    print(f"   Sliding window analysis complete")

#The following version handles cases where FDR correction is performed after filtering for p-value < 0.05, such as for plant CHH and CHG contexts
def perform_sliding_window_on_dmp_files_after_filter(output_dir, replicate_x, replicate_y,all_dfs_ndmp_dict=None, methylation_type=None):
    """Perform sliding window analysis on the (N)DMP file resulting from 1 of the m*n*3 processing iterations"""

    print(f"\n    Starting sliding window analysis on DMP files...")

    # Process only DMP files
    dmp_patterns = [
        f"DMP_wt_replicate{replicate_y}_mut_replicate{replicate_x}_Chr*.txt",
        f"N-DMP_wt_replicate{replicate_y}_mut_replicate{replicate_x}_Chr*.txt"
    ]

    for pattern in dmp_patterns:
        dmp_files = glob.glob(os.path.join(output_dir, pattern)) #Find the paths of all files matching the current regular expression format and put them into the dmp_files list
        is_ndmp = pattern.startswith("N-DMP") # Determine if the file currently being processed is an NDMP file

        for dmp_file in dmp_files: # Iterate through the list, obtaining one DMP or N-DMP file at a time
            try:

                df = pd.read_csv(dmp_file, sep=' ', header=None, skiprows=1,
                                 names=['position', 'pvalue', 'change'])

                #If the current file is an ndmp file, we cannot modify the all_simple file (as it will be used to merge into a single column for FDR correction), so we append the original information of sites with p-value > 0.05 here
                                # This ensures no omissions when calculating ndmp information during the subsequent generation of sliding window files
                if is_ndmp and all_dfs_ndmp_dict and methylation_type:
                    # Extract chromosome number from the filename
                    chr_match = re.search(r'Chr(\d+)\.txt$', dmp_file)
                    if chr_match:
                        chr_num = int(chr_match.group(1))
                        if (methylation_type in all_dfs_ndmp_dict and
                                chr_num in all_dfs_ndmp_dict[methylation_type]):
                            ndmp_df = all_dfs_ndmp_dict[methylation_type][chr_num]
                            if isinstance(ndmp_df, pd.DataFrame) and not ndmp_df.empty:
                                ndmp_df = ndmp_df.rename(columns={
                                    'Position': 'position',
                                    'Pvalue': 'pvalue',
                                    'Methylation_Change': 'change'
                                })
                                # Merge data
                                df = pd.concat([df, ndmp_df[['position', 'pvalue', 'change']]],
                                               ignore_index=True)
                                df = df.drop_duplicates(subset=['position'])
                                df = df.sort_values('position').reset_index(drop=True)
                                print(f" Merged {len(ndmp_df)} NDMP sites previously ignored during FDR correction for {os.path.basename(dmp_file)}")
                # Check if data is not empty
                if df.empty or len(df) == 0:
                    continue

                # Data type conversion and cleanup
                df = df.dropna()
                df['position'] = df['position'].astype(int)
                df['change'] = df['change'].astype(int)

                print(f"     process file: {os.path.basename(dmp_file)} ({len(df)} sites)")

                # Set the output prefix, which essentially removes '.txt' to retain the (N)DMP_wt_replicate{replicate_y}_mut_replicate{replicate_x}_Chr* portion
                base_name = os.path.basename(dmp_file).replace('.txt', '')

                # Pass the DataFrame to sliding_window_analysis
                # The former format is window_start, window_end, count_change_1, count_change_0, total_count
                # The latter lacks count_change but includes standardized_count, which is the ratio of significant sites in each interval to the maximum number of significant sites across all intervals
                sliding_results, std_results = sliding_window_analysis(
                    df, # df is the DMP file for a specific chromosome from 1 of the mn3 processing iterations
                    window_size=1000000,
                    step_ratio=0.05,
                    save_files=True,
                    output_identifier=base_name,
                    outputdir1=output_dir
                )

                print(f"      Completed {base_name}: {len(sliding_results)} windows")

            except Exception as e:
                print(f"      Failed to process {dmp_file} : {e}")
                continue

    print(f"    Sliding window analysis complete ")


def generate_dmp_files(dir1,dir2,output_dir, replicate_x, replicate_y, fdr_threshold=0.05, mtype1="CpG",
                      all_dfs_ndmp_dict=None,unfilter_mtypes=["CpG"],n_chromosomes = 5):
    """Parameters: output_dir is output_x_y/methylation_type/,
       replicate_x, replicate_y are the processed group numbers, followed by the q-value threshold and the currently processed methylation type"""

    print(f"\n    Generating DMP files for {output_dir} ...")

    def safe_float_convert(value):
        """Convert data of various formats into floating-point numbers"""
        try:
            # If already a numeric type, convert to float and return
            if isinstance(value, (int, float)):
                return float(value)
            # If a string, strip whitespace and convert to float
            if isinstance(value, str):
                value = value.strip()
                if value:
                    return float(value)
            return None
        except (ValueError, TypeError):
            return None

    # 1. Read the FDR correction results file (here output_dir is output_x_y/methylation_type/)
    #  The format of this file is Chromosome, Methylation_Type, Position, Pvalue, Qvalue
    fdr_file = os.path.join(output_dir, f"FDR_corrected_results_wt_replicate{replicate_y}_vs_mut_replicate{replicate_x}.txt")
    if not os.path.exists(fdr_file):
        print(f"   Error: FDR results file does not exist: {fdr_file}")
        return False

    try:
        fdr_df = pd.read_csv(fdr_file, sep='\t') # The format of this file is Chromosome, Methylation_Type, Position, Pvalue, Qvalue
        print(f"    Read FDR results file, total {len(fdr_df)} sites")
    except Exception as e:
        print(f"   Error: Failed to read FDR file {e}")
        return False

    # 2. Read the all_simple files for each chromosome
    methylation_change_data = {}
    total_read_lines = 0

    # (here output_dir is output_x_y/methylation_type/)
    mtype_dir = output_dir

    # Retrieve the filenames of all all_simple files processed in 1 of the m*n*3 processing iterations into a list
    all_simple_files = [f for f in os.listdir(mtype_dir) if f.startswith('all_simple_Chr') and f.endswith('.txt')]

    # The format of the all_simple file is: "Position", "Pvalue", "Methylation_Change"
    # Iterate through all all_simple files from the given processing iteration
    for file in all_simple_files:
        chr_match = re.search(r'Chr(\d+)\.txt$', file)
        if not chr_match:
            continue

        chr_num = int(chr_match.group(1))# chr_num here also refers to the nth set of three columns in the 'both' file
        file_path = os.path.join(mtype_dir, file)

        try:
            with open(file_path, 'r') as f:  # The read all_simple file format is: "Position", "Pvalue", "Methylation_Change"
                lines = f.readlines() # readlines is used to return a list composed of all lines in the file, where each element is a line

            file_valid_lines = 0
            for line_num, line in enumerate(lines, 1): # Iterate through each line, with line_num enumerating starting from 1
                line = line.strip()
                if not line or line.startswith("Position"): # Skip empty lines or the header row
                    continue

                parts = line.split('\t')    #Split based on the delimiter set during the previous export, placing the three data items into the 'parts' list each time
                if len(parts) >= 3:
                    # Since the input consists of strings, type conversion is needed: convert pvalue to float, and the other two to integers
                    position = int(safe_float_convert(parts[0]))
                    pvalue = safe_float_convert(parts[1])
                    change = int(safe_float_convert(parts[2]))

                    # Check if all values are valid
                    if (position is not None and
                            pvalue is not None and
                            change is not None):

                        # Verify if 'change' is either 0 or 1
                        if change in [0, 1]:
                            # Store this mapping: (chr, mtype, position) -> change in the dictionary
                            methylation_change_data[(chr_num, mtype1, position)] = change
                            file_valid_lines += 1
                            total_read_lines += 1

            print(f"    {mtype_dir}/Chr{chr_num}: Read {file_valid_lines} valid sites")

        except Exception as e:
            print(f"    Warning: Failed to read {file_path} : {e}")
            continue

    print(f"   Total read methylation change direction data: {total_read_lines} sites")

    # 3. Merge data: use the previously recorded (chr, mtype, position) -> change dictionary mapping to add the 'change' attribute to a copy of fdr_df
    combined_data = []
    missing_change = 0
    match_debug = defaultdict(int) #Create a dictionary with a default value; when accessing a non-existent key, the dictionary will automatically call the built-in function of the int class
                                                # to create a key-value pair for that key, setting the value to the return value of int() — which is 0

    for _, row in fdr_df.iterrows(): # The fdr_df file format is Chromosome, Methylation_Type, Position, Pvalue, Qvalue
        chr_num = row['Chromosome']
        mtype = row['Methylation_Type']
        # Uniformly convert to floating-point numbers for matching
        position = safe_float_convert(row['Position'])
        qvalue = safe_float_convert(row['Qvalue'])

        # If position or q-value is missing, skip the current row
        if position is None or qvalue is None:
            missing_change += 1
            continue

        # Look up the corresponding methylation change direction
        change_key = (chr_num, mtype, position)
        if change_key in methylation_change_data:
            change = methylation_change_data[change_key]  # This change_data stores the mapping: (chr, mtype, position) -> change
            combined_data.append({
                'chromosome': chr_num,
                'methylation_type': mtype,
                'position': int(position),
                'qvalue': qvalue,
                'change': change # This is equivalent to adding the 'change' attribute to a copy of fdr_df, but at this point, combined_data is still a list rather than a dataframe
            })
            match_debug[chr_num] += 1 # The chr_num here also represents the nth set of three columns in the 'both' file; this has been the case since reading the 'Both' file for Fisher's exact test
        else:
            missing_change += 1

    print(f"    Merged into {len(combined_data)} complete sites")
    print(f"    Chromosome matching status: {dict(match_debug)}") # This actually indicates how many pieces of information each chromosome has
    if missing_change > 0:
        print(f"    Warning:{missing_change} sites are missing methylation change direction information")

    # 4. Group by chromosome to generate DMP files
    chr_groups = defaultdict(list) # This defaultdict calls the list() function upon initialization, using an empty list [] as the default value
    for item in combined_data: # Each item is a dictionary where the values of the five key-value pairs are chr_num, mtype, position, qvalue, change
        chr_groups[item['chromosome']].append(item)
            # The chr_groups dictionary records each piece of data (in dictionary form) for each chromosome as an element in the list at chr_groups -> chr_num

    if not chr_groups:
        print(f"    Error: No data found to generate DMP files")
        return False

    print(f"    Will generate DMP files for the following chromosomes: {sorted(chr_groups.keys())}")

    dir1_name = f"replicate{replicate_x}"
    dir2_name = f"replicate{replicate_y}"

    total_dmp = total_ndmp = total_hyper = total_hypo = 0

    # Generate files for each chromosome
    for chr_num in sorted(chr_groups.keys()): # The keys here are chromosome indices; index - 1 is the index in the initial total chromosome mapping, since that one is zero-based
        chr_data = chr_groups[chr_num]  # The chr_groups dictionary records each piece of data (in dictionary form) for each chromosome as an element
                                    #So chr_data here is still a dictionary
        chr_data.sort(key=lambda x: x['position']) #A lambda anonymous function is applied to each element of the list, using the return value as the basis for sorting
                                #  The anonymous function here is equivalent to applying this function to each element (which is a dictionary) in chr_data
                                                #   def get_position(x):
                                            #     return x['position'] The return value is the position, so it sorts based on the position

        # Generate filenames (here 'output' refers to the methylation type directory)
        dmp_file = os.path.join(output_dir, f"DMP_wt_{dir2_name}_mut_{dir1_name}_Chr{chr_num}.txt")
        ndmp_file = os.path.join(output_dir, f"N-DMP_wt_{dir2_name}_mut_{dir1_name}_Chr{chr_num}.txt")
        hyper_file = os.path.join(output_dir, f"hyper_DMP_wt_{dir2_name}_mut_{dir1_name}_Chr{chr_num}.txt")
        hypo_file = os.path.join(output_dir, f"hypo_DMP_wt_{dir2_name}_mut_{dir1_name}_Chr{chr_num}.txt")

        # Categorize data
        dmp_data = []
        ndmp_data = []
        hyper_data = []
        hypo_data = []

        # Iterate through chr_data
        for item in chr_data: # chr_data contains items that are dictionaries, where the values of the five key-value pairs are chr_num, mtype, position, qvalue, change
            position = item['position']
            qvalue = item['qvalue']
            change = item['change']

            if qvalue <= fdr_threshold: #Determine if significant
                dmp_data.append((position, qvalue, change))
                if change == 1:
                    hyper_data.append((position, qvalue, change))
                elif change == 0:
                    hypo_data.append((position, qvalue, change))
            else:
                ndmp_data.append((position, qvalue, change))

        if mtype1 not in unfilter_mtypes and all_dfs_ndmp_dict and mtype1 in all_dfs_ndmp_dict:
            if chr_num in all_dfs_ndmp_dict[mtype1]:
                ndmp_df = all_dfs_ndmp_dict[mtype1][chr_num]
                if isinstance(ndmp_df, pd.DataFrame) and not ndmp_df.empty:
                    # Get the set of existing sites to avoid duplicate additions
                    existing_positions = set(item[0] for item in ndmp_data)
                    added_count = 0
                    for _, row in ndmp_df.iterrows():
                        pos = int(row['Position'])
                        pval = float(row['Pvalue'])
                        chg = int(row['Methylation_Change'])
                        # Add only unique sites
                        if pos not in existing_positions:
                            ndmp_data.append((pos, pval, chg))
                            added_count += 1
                    if added_count > 0:
                        print(f"    Added {added_count} NDMP sites with p-value > 0.05 for Chr{chr_num}")


        # Write to files
        def write_dmp_file(filename, data):
            with open(filename, 'w') as f:
                f.write("first line\n")
                for pos, qval, chg in data:
                    f.write(f"{pos} {qval} {chg}\n")

        write_dmp_file(dmp_file, dmp_data)
        write_dmp_file(ndmp_file, ndmp_data)
        write_dmp_file(hyper_file, hyper_data)
        write_dmp_file(hypo_file, hypo_data)

        run_dmr_pipeline_on_dmp_file(dmp_file,chromoNo=n_chromosomes)

        print(f"    Chr{chr_num}: DMP={len(dmp_data)}, N-DMP={len(ndmp_data)}, Hyper={len(hyper_data)}, Hypo={len(hypo_data)}")

        # generate_dmp_files(output_dir, replicate_x, replicate_y, fdr_threshold=0.05, mtype1="CpG",
        #                    all_dfs_ndmp_dict=None, unfilter_mtypes=["CpG"]):
        # Statistics
        total_dmp += len(dmp_data)
        total_ndmp += len(ndmp_data)
        total_hyper += len(hyper_data)
        total_hypo += len(hypo_data)

    print(f"    DMP file generation complete!")
    print(f"    Total: DMP={total_dmp}, N-DMP={total_ndmp}, Hyper={total_hyper}, Hypo={total_hypo}")

    bothfile1 = os.path.join(dir1,f"{replicate_x}-bothMeUnme_diffChromo_NOREPEATED_methy_sites_{mtype1}.txt")

    bothfile2 = os.path.join(dir2,f"{replicate_y}-bothMeUnme_diffChromo_NOREPEATED_methy_sites_{mtype1}.txt")

    summarize_dmr_methylation(output_dir, replicate_x, replicate_y, bothfile1, bothfile2, mtype1,custom_dmr_dir=None)

    # Decide which sliding window function to use based on the methylation type
    if mtype1 not in unfilter_mtypes:
        # Methylation types should also be differentiated here
        perform_sliding_window_on_dmp_files_after_filter(
            output_dir, replicate_x, replicate_y,
            all_dfs_ndmp_dict=all_dfs_ndmp_dict,
            methylation_type=mtype1
        )
        # perform_sliding_window_on_dmp_files(
        #     output_dir, replicate_x, replicate_y
        # )
    else:
        #print(f"  Using the standard version of sliding window analysis")
        perform_sliding_window_on_dmp_files(output_dir, replicate_x, replicate_y)
    return True


# Integrate this function into process_replicate_pair
def process_replicate_pair(replicate_x, replicate_y, files1, files2, dir1, dir2, dir1_name, dir2_name,unfilter_mtypes,work_dir="."):
    """Process all methylation types for a pair of replicates; represents processing 1 set of 3 tests
    replicate_x, replicate_y are the indices of the 'both' files
    'files' is a dictionary mapping: index -> methylation type -> corresponding filename
    dir1, dir2 are the directory names for the two genotype datasets
    dir1_name, dir2_name are the final parts of the provided dir1, dir2 paths"""

    print(f"\n  Processing replicate pair (wt{replicate_y}, mut{replicate_x})...")

    # Create the output directory for this pair
    output_dir = os.path.join(work_dir, f"output_wt{replicate_y}_mut{replicate_x}")
    os.makedirs(output_dir, exist_ok=True)

    methylation_types = ['CpG', 'CHH', 'CHG']

    all_dfs_ndmp_dict = {}

    # Loop through to process each methylation type
    for methylation_type in methylation_types:
        if methylation_type not in all_dfs_ndmp_dict:
            all_dfs_ndmp_dict[methylation_type] = {}

        success_count = total_tests = 0
        # Here, 'in files1[replicate_x]' checks if methylation_type exists as a key in the files1[replicate_x] dictionary; if it does, it indicates
        # the corresponding i-both...methylation_type.txt file exists
        if (methylation_type not in files1[replicate_x] or
                methylation_type not in files2[replicate_y]):
            print(f"   Skipping methylation type {methylation_type}: file does not exist")
            continue
        # Otherwise, get the relative path of the file currently needing to be processed
        file1_path = os.path.join(dir1, files1[replicate_x][methylation_type])
        file2_path = os.path.join(dir2, files2[replicate_y][methylation_type])

        # Get the number of chromosomes in the two files
        n_chromosomes_1 = get_column_count(file1_path)
        n_chromosomes_2 = get_column_count(file2_path)

        if n_chromosomes_1 is None or n_chromosomes_2 is None:
            print(f" Unable to get the number of chromosomes for {methylation_type}")
            continue

        if n_chromosomes_1 != n_chromosomes_2:
            print(f"   Inconsistent chromosome counts for {methylation_type}: {n_chromosomes_1} vs {n_chromosomes_2}")
            continue

        # At this point, both 'both' files to be processed contain chromosome data and have the same total number of columns
        # However, this check is not strictly necessary, as newtoboth guarantees identical column counts from the start
        n_chromosomes = n_chromosomes_1  # Get the total number of chromosomes
        print(f"    Processing methylation type {methylation_type}, total {n_chromosomes} chromosomes")

        # Process each chromosome for the current methylation type in the current x_y group; chr_num here is the nth chromosome, not the chromosome name
        for chr_num in range(1, n_chromosomes + 1):
            if methylation_type in unfilter_mtypes:
                success = process_methylation_type_with_collection(
                    file1_path, file2_path, methylation_type, output_dir,
                    dir1_name, dir2_name, chr_num
                )  # Process one chromosome for one test
            else:
                # Differentiate methylation types here
                all_df_ndmp = process_methylation_type_with_collection_pvfilter(
                    file1_path, file2_path, methylation_type, output_dir,
                    dir1_name, dir2_name, chr_num
                )   #Return the dataframe with p-value > 0.05 for this iteration
                    #  all_results([pos, pvalue, change])
                if chr_num not in all_dfs_ndmp_dict[methylation_type]:
                    all_dfs_ndmp_dict[methylation_type][chr_num] = all_df_ndmp
                    if (isinstance(all_df_ndmp, pd.DataFrame)):
                        print(f"Successfully added the specific dataframe for {methylation_type}-{chr_num} to the dictionary; this dataframe contains {len(all_df_ndmp)} rows")
                success = isinstance(all_df_ndmp, pd.DataFrame)
            total_tests += 1
            if success:  # Returns True if the previous step proceeded normally, False otherwise
                success_count += 1

        # Merge FET results and perform FDR correction
        #  The format of the FET results is: pos, m1, u1, m2, u2, pvalue
        if success_count > 0:
            output_dir1 = os.path.join(output_dir, methylation_type)
            # Retrieve all FET files in an output_x_y/methylation_type/ directory, concatenate them, perform FDR correction to get q-values, and export to disk. The final format is:
            #                                        Chromosome, Methylation_Type, Position, Pvalue, Qvalue
            merge_fet_results_and_fdr(output_dir1, replicate_x, replicate_y, methylation_type,all_dfs_ndmp_dict,n_chromosomes)
            # Generate DMP files
            # Get the corresponding DMP threshold based on the methylation type
            dmp_threshold = get_dmp_threshold(methylation_type)
            generate_dmp_files(dir1,dir2,output_dir1, replicate_x, replicate_y, fdr_threshold=dmp_threshold, mtype1=methylation_type,all_dfs_ndmp_dict=all_dfs_ndmp_dict
                               ,unfilter_mtypes=unfilter_mtypes,n_chromosomes=n_chromosomes)

    print(f"  Replicate pair (wt{replicate_y}, mut{replicate_x}) processing complete!{success_count}/{total_tests} tests successful")
    return success_count, total_tests ,    # Here, one chromosome counts as one test


def process_all_combinations(dir1, dir2, m, n,unfilter_mtypes,work_dir="."):
    """Process all combinations, performing mn3 tests"""

    print(f"Scanning file directories...")
    # 'files' is a mapping: index -> methylation type -> filename, allowing retrieval of the corresponding filename in the directory via index and methylation type
    files1 = scan_sample_files_by_replicates(dir1, m)
    files2 = scan_sample_files_by_replicates(dir2, n)

    print(f"Directory 1 ({dir1}): found {len(files1)} sets of files")
    print(f"Directory 2 ({dir2}): found {len(files2)} sets of files")

    # Check for missing 'both' files; if any, record their indices
    missing_replicates1 = [i for i in range(1, m + 1) if i not in files1]
    missing_replicates2 = [i for i in range(1, n + 1) if i not in files2]
    # Output the indices of missing 'both' files
    if missing_replicates1:
        print(f"Warning: Directory 1 is missing these sets: {missing_replicates1}")
    if missing_replicates2:
        print(f"Warning: Directory 2 is missing these sets: {missing_replicates2}")

    # Record the indices of the acquired 'both' files
    available_replicates1 = [i for i in range(1, m + 1) if i in files1]
    available_replicates2 = [i for i in range(1, n + 1) if i in files2]

    # Calculate the total number of combination groups needed
    total_combinations = len(available_replicates1) * len(available_replicates2)
    print(f"\nStarting to process {total_combinations} combinations...")

    # rstrip(char) removes trailing specified characters from the right; here the specific character is os.sep, i.e., the current system path separator: \
    # After ensuring no trailing path separator, use the basename function to get the last part of the given path (if \ is not removed, an empty string is obtained)
    # That is, obtaining the corresponding directory name—the genotype
    dir1_name = sanitize_filename(os.path.basename(dir1.rstrip(os.sep)))
    dir2_name = sanitize_filename(os.path.basename(dir2.rstrip(os.sep)))

    total_success = total_tests = 0
    start_time = time.time()
    # enumerate starts counting from 0, where i is the 0-based index. replicate_x is each element of available_replicates1, here representing each actually existing index sorted in ascending order
    for i, replicate_x in enumerate(available_replicates1):
        # j is the 0-based index. replicate_y is each element of available_replicates2, here representing each actually existing index sorted in ascending order
        for j, replicate_y in enumerate(available_replicates2):
            print(f"\nnProgress: {i * len(available_replicates2) + j + 1}/{total_combinations}")
            # Process the mut_replicate_X-wt_replicate_y group
            success_count, test_count = process_replicate_pair(
                replicate_x, replicate_y, files1, files2, dir1, dir2, dir1_name, dir2_name,unfilter_mtypes,work_dir=work_dir
            )  # process_replicate_pair handles 13 tests
            total_success += success_count # Here, one chromosome still counts as one test
            total_tests += test_count  # Here, one chromosome still counts as one test

    end_time = time.time()
    print(f"\nAll processing complete!")
    print(f"Total: {total_success}/{total_tests} successful tests")
    print(f"Time elapsed: {end_time - start_time:.2f} seconds")

    # Output instructions
    print(f"Results for single comparisons are saved in the ./output_x_y/methylation_type/ directory")

    return total_success == total_tests #Return True if all succeed, otherwise return False

def bayes_deciding(sig_count, nonsig_count):


    prob_gt_half = sig_count/(sig_count+nonsig_count)
    final_decision = 1 if prob_gt_half >= 2/3 else 0
    # print(f"\nDecision (Threshold={recommended_threshold * 100:.0f}%)")
    # print(f"  Decision result: {'Significant' if final_decision else 'Not significant'}")
    # print(f"  Confidence: {prob_gt_half * 100:.1f}%")

    return final_decision

def find_common_significant_sites(output_dirs=None, methytype2='CpG', dir1=None, dir2=None,work_dir="."):
    """
    Identify sites that are significant across all combination tests and retrieve relevant information for these sites
    Parameters：
        output_dirs: List of output directories; if None, scan automatically
        methytype2: Methylation type
    """

    print("\nSearching for common significant sites across all combinations...")

    and_output_dir = os.path.join(work_dir, "and_output")
    os.makedirs(and_output_dir, exist_ok=True)

    # 1. Automatically scan all output_x_y directories
    if output_dirs is None:
        output_dirs = glob.glob(os.path.join(work_dir,f"output_*_*/{methytype2}"))
        # Find all directories formatted as output_x_y/methytype2/; there are mn in total
        output_dirs = [d for d in output_dirs if os.path.isdir(d)]

    if not output_dirs: # If not found
        print("No output directories found")
        return None

    print(f"Found {len(output_dirs)} output directories")

    # 2. Read all FDR_corrected files into memory at once, including all sites. If only significant sites were read,
            # and a site appeared in the 1_1 sig file but not in a subsequent test, it would be unclear whether it was absent due to non-significance or because it wasn't in the original input data for that iteration
    valid_dirs = [] # Collect paths of directories that play a role in subsequent operations
    site_statistics = {}  # Build the mapping: {site_id: {'sig': 0, 'total': 0}}
    all_dataframes = {} # Ultimately access the df corresponding to FDR_correct via all_dataframes[directory]
    dir_to_replicate = {} #Record the replicate index corresponding to the directory
    for output_dir in output_dirs: # Iterate through all output_x_y/methytype2/
        # Get the FDR_correct file in the current methylation directory,
            # format is: 'Chromosome', 'Methylation_Type', 'Position', 'Pvalue', 'Qvalue'
        # Here, unlike before where q-values were all below the threshold, it includes all q-values
        match = re.search(r'output_wt(\d+)_mut(\d+)', output_dir)
        if not match:
            continue

        replicate_y = int(match.group(1))
        replicate_x = int(match.group(2))


        fdr_all_files = glob.glob(os.path.join(output_dir, "FDR_corrected_results_*.txt"))

        # Check for existence
        if not fdr_all_files:
            print(f" Warning: FDR_corrected file not found in {output_dir}")
            continue

        fdr_all_file = fdr_all_files[0] # Since glob.glob returns a list, use [0] to get the actually existing file path

        try:
            df = pd.read_csv(fdr_all_file, sep='\s+') # Read this file into df. The FDR_corrected format is:
                                        # 'Chromosome', 'Methylation_Type', 'Position', 'Pvalue', 'Qvalue'
            if len(df) > 0:
                # Create a unique identifier for each row
                #  Specifically: for each row in df, extract the chromosome number, methylation type, and position, and append it to the right as a new column 'site_id'
                #    This is because we need to tally the significance of all sites across multiple tests later; site_id allows distinguishing different sites using just one attribute, eliminating the need to
                #     successively check multiple attributes for equality to filter or access them (chromosome number, methylation type, position). (Note: the chromosome number here is the nth set of three columns in the 'both' file)
                df['site_id'] = df.apply(
                    lambda row: f"{int(row['Chromosome'])}-{row['Methylation_Type']}-{int(row['Position'])}",
                    axis=1
                )
                all_dataframes[output_dir] = df # all_dataframes[directory] -> corresponding df for FDR_correct
                valid_dirs.append(output_dir)
                dir_to_replicate[output_dir] = (replicate_x, replicate_y)  # Record index
                # Get the DMP threshold for the current methylation type
                dmp_threshold = get_dmp_threshold(methytype2)
                for _, row in df.iterrows(): # Iterate through each row; each row contains the information of a site, whether significant or not, from FDR_correct in a single test
                    site_id = row['site_id']

                    if site_id not in site_statistics:
                        site_statistics[site_id] = {
                            'sig_count': 0,
                            'total_count': 0,
                            'chromosome': int(row['Chromosome']),
                            'methylation_type': row['Methylation_Type'],
                            'position': int(row['Position'])
                        }

                    site_statistics[site_id]['total_count'] += 1

                    if row['Qvalue'] <= dmp_threshold:
                        site_statistics[site_id]['sig_count'] += 1


                print(f"  {output_dir}: {len(df)} sites")
            else:
                print(f"  {output_dir}: No significant sites")
        except Exception as e:
            print(f"  Error: Failed to read {fdr_all_file} : {e}")
            continue
    valid_dirs.sort(key=lambda d: dir_to_replicate[d])
    print(f"\nStatistics complete, total {len(site_statistics)} distinct sites")

    if not site_statistics:
        print("No valid site information results found")
        return None

    # 3. Obtain sites determined as significant by the Bayesian method and add them to common_sites
    common_sites = []
    for site_id, stats in site_statistics.items():
        if stats['total_count'] != len(valid_dirs):
            stats['total_count'] = len(valid_dirs)
        sig_count = stats['sig_count']  # Obtain the number of times the current site was tested as significant
        nonsig_count = stats['total_count'] - sig_count  # Obtain the number of times the current site was tested as non-significant
        is_significant = bayes_deciding(sig_count, nonsig_count)
        if is_significant:
            common_sites.append(site_id)
    if not common_sites:
        print("No sites were significant across all combinations")
        return None
    else:
        print(f"\nAfterdecision, a total of {len(common_sites)} significant sites")

    # 4. Read the methylation change direction information for each directory
    print("Reading methylation change direction information...")
    methylation_change_by_dir = {}  # {output_dir: {site_id: change}}

    for output_dir in valid_dirs: # Iterate through the above valid directories for this methylation type
        methylation_change_by_dir[output_dir] = {} # Create dictionary element for the current directory,
                                        # Construct the mapping: output_dir -> site_id -> change
        # Search for all all_simple_Chr files in this directory, as they contain 'change' information in the format pos, pvalue, change,
                                                # and the chromosome number is obtained from the filename
        all_simple_files = glob.glob(os.path.join(output_dir, "all_simple_Chr*.txt"))
                                                        # Get the paths of all all_simple files in the current directory and form a list

        for file_path in all_simple_files: # Iterate through all_simple files, format is pos, pvalue, change
            chr_match = re.search(r'Chr(\d+)\.txt$', file_path) # Create regular expression and capture group to extract the chromosome number
            if not chr_match:
                continue
            chr_num = int(chr_match.group(1)) # Obtain the chromosome number of the current file via the capture group

            try:
                df = pd.read_csv(
                    file_path,
                    sep='\t',
                    usecols=[0, 2],
                    names=['position', 'change'],
                    dtype={'position': float, 'change': float}, # float can convert any numeric string without errors, making it a more robust approach
                    skiprows=1
                )

                # Converting pos and change to int at this point avoids issues; otherwise, str to int conversion might fail (e.g., converting "100.0" to int throws an error)
                df['position'] = df['position'].astype(int)
                df['change'] = df['change'].astype(int)

                positions = df['position'].values
                changes = df['change'].values

                for position, change in zip(positions, changes):
                    site_id = f"{chr_num}-{methytype2}-{position}"
                    methylation_change_by_dir[output_dir][site_id] = change

            except Exception as e:
                print(f"  Warning: Failed to read {file_path} : {e}")
                continue

        print(f"  {output_dir}: Read change directions for {len(methylation_change_by_dir[output_dir])} sites")

    # 5. Process common site information
    print("Processing common site detailed information...")
    common_site_details = []

    # Create an index for each DataFrame to speed up queries
    indexed_dfs = {} # Establish the mapping: output_dir -> df1 (with site_id as the index)
    for output_dir in valid_dirs:
        df = all_dataframes[output_dir] # Obtain the df corresponding to FDR_correct in each directory (which already includes the site_id column)
        indexed_dfs[output_dir] = df.set_index('site_id') # Set site_id as the index and return the new df1 as the dictionary value

    # Process common sites one by one
    for i, site_id in enumerate(common_sites):
        if i % 10000 == 0:  # Print progress every 10,000 sites processed
            print(f"  Processed {i}/{len(common_sites)} sites")

        site_info = {'site_id': site_id}
        chr_num, mtype, pos = site_id.split('-')
        site_info['Chromosome'] = int(chr_num)
        site_info['Methylation_Type'] = mtype
        site_info['Position'] = int(pos)

        # List to collect q-values: collects q-values across all tests for the same chromosome and position (matching the current site_id info) in different output_x_y/methytype2/ directories
        qvalues = []
        # Collect change directions across all tests for the same chromosome and position (matching the current site_id info) in different output_x_y/methytype2/ directories
        change_values = []
        qvalue_dict = {}
        for output_dir in valid_dirs:
            replicate_x, replicate_y = dir_to_replicate[output_dir]
            col_name = f'qvalue_{os.path.basename(dir2.rstrip("/"))}{replicate_y}_{os.path.basename(dir1.rstrip("/"))}{replicate_x}'   # Added: column name, prefix is the wild-type index, suffix is the mutant-type index
            indexed_df = indexed_dfs[output_dir] # output_dir -> df1 (with site_id as index), retrieve a specific FDR_correct file for the current methylation type indexed by site_id
                                # Content format is: 'Chromosome', 'Methylation_Type', 'Position', 'Pvalue', 'Qvalue'
            if site_id in indexed_df.index:
                qval = indexed_df.loc[site_id, 'Qvalue']
                qvalues.append(qval)
                qvalue_dict[col_name] = qval
                # Get the change value for this site in this directory (output_dir -> site_id -> change mapping)
                if site_id in methylation_change_by_dir[output_dir]:
                    change_values.append(methylation_change_by_dir[output_dir][site_id])
            else:
                qvalue_dict[col_name] = 1.0

        if qvalues:  # Ensure q-value data exists
            dmp_threshold = get_dmp_threshold(methytype2)
            qvalues_sig = [q for q in qvalues if q <= dmp_threshold]
            if len(qvalues_sig) > 0:
                site_info['Sig_Mean_Qvalue'] = np.mean(qvalues_sig)
            else:
                site_info['Sig_Mean_Qvalue'] = 1
            # site_info['Max_Qvalue'] = np.max(qvalues)
            # site_info['Min_Qvalue'] = np.min(qvalues)
            site_info['Num_Comparisons'] = len(qvalues)

            # Calculate methylation change direction via voting
            if change_values:
                # Count occurrences where change == 1
                num_hyper = sum(change_values)
                total_comparisons = len(change_values)
                hyper_ratio = num_hyper / total_comparisons

                # Majority vote: if >= 50%, record as 1 (hypermethylated), otherwise 0
                site_info['Methylation_Change'] = 1 if hyper_ratio >= 0.5 else 0
                site_info['Hyper_Count'] = num_hyper  # Hypermethylation count
                site_info['Hypo_Count'] = total_comparisons - num_hyper  # Hypomethylation count
                site_info['Hyper_Ratio'] = hyper_ratio  # Hypermethylation ratio
            else:
                # If no change information, mark as missing (though theoretically, one test exactly corresponds to one qvalue and one change)
                site_info['Methylation_Change'] = -1  # -1 indicates undetermined
                site_info['Hyper_Count'] = 0
                site_info['Hypo_Count'] = 0
                site_info['Hyper_Ratio'] = 0

            site_info.update(qvalue_dict) # Add q-value columns for all replicates

            common_site_details.append(site_info) # This list records dictionaries, each containing various information corresponding to a site_id. Format:
        # site_id-Chromosome-Methylation_Type-Position-Total_Tests-change-hypercount-hypocount-hyperratio-Sig_Mean_Qvalue-All_qvalues
            # (The chromosome number here is the nth set of 3 columns in the 'both' file)

    print(f" Finished processing {len(common_site_details)} sites")

    # 6. Generate result DataFrame
    result_df = pd.DataFrame(common_site_details)
    result_df = result_df.sort_values(['Methylation_Type', 'Chromosome', 'Position']) # Sort

    # Adjust column order, placing Methylation_Change further up front
    column_order = [
        'Chromosome', 'Methylation_Type', 'Position',
        'Methylation_Change', 'Hyper_Ratio', 'Hyper_Count', 'Hypo_Count', 'Num_Comparisons',
        'Sig_Mean_Qvalue'
    ]
    # Obtain all replicate column names and sort
    # New column name format: dir2_y_dir1_x; need to retrieve columns containing numbers separated by underscores
    replicate_columns = sorted(
        [col for col in result_df.columns if col.startswith('qvalue_')],
        key=lambda x: tuple(map(int, re.findall(r'\d+', x)))  # Extract all numbers and sort
    )
    column_order = column_order + replicate_columns
    result_df = result_df[column_order]

    # Save results
    output_file = os.path.join(and_output_dir, f"{methytype2}-final_significant_sites_DMPs.txt")
    result_df.to_csv(output_file, sep='\t', index=False)
    print(f"\nCommon significant sites saved to: {output_file}")

    # Output statistical information
    print("\nCommon significant site statistics:")
    mtype = methytype2
    mtype_df = result_df
    count = len(mtype_df)
    hyper_count = len(mtype_df[mtype_df['Methylation_Change'] == 1])
    hypo_count = len(mtype_df[mtype_df['Methylation_Change'] == 0])
    unknown_count = len(mtype_df[mtype_df['Methylation_Change'] == -1])

    print(f"  {mtype}: {count} sites")
    print(f"    - Hypermethylated(Change=1): {hyper_count} ({hyper_count / count * 100:.1f}%)")
    print(f"    - Hypomethylated (Change=0): {hypo_count} ({hypo_count / count * 100:.1f}%)")
    if unknown_count > 0:
        print(f"    - Unknown(Change=-1): {unknown_count} ({unknown_count / count * 100:.1f}%)")

    return result_df # The format is:'Chromosome', 'Methylation_Type', 'Position',
                            # 'Methylation_Change', 'Hyper_Ratio', 'Hyper_Count', 'Hypo_Count',
                           # 'Sig_Mean_Qvalue', 'Max_Qvalue', 'Min_Qvalue', 'Num_Comparisons'

def sliding_window_analysis(
        input_data,
        window_size = 1000000,
        step_ratio = 0.05,
        save_files = False,
        output_identifier = None,
        outputdir1="./"
):
    """
    Perform sliding window analysis on methylation site data
    Parameters：
    input_data: DataFrame converted from the (N)DMP file for a specific chromosome from 1 of the mn3 processing iterations, containing columns: ['position', 'pvalue', 'change']
    window_size: Size of the sliding window
    step_ratio: Step size ratio (percentage of the window size)
    save_files : save_files: Whether to save the result files
    output_identifier: Output file prefix; auto-generated if save_files=True and not provided
    """

    # 1. Data reading and preprocessing
    df = input_data.copy()

    # Ensure correct column names
    expected_cols = ['position', 'pvalue', 'change']
    if not all(col in df.columns for col in expected_cols):
        raise ValueError(f"DataFrame必须包含列: {expected_cols}")

    if output_identifier is None:
        output_identifier = "sliding_window_analysis"

    # Data validation and cleaning
    df = df.dropna()
    if len(df) == 0:
        raise ValueError("No valid data rows")

    # Sort by position
    df = df.sort_values('position').reset_index(drop=True) #Reset dataframe index without adding the original index as a new column

    print(f"Data preprocessing complete, total {len(df)} sites")

    # 2. Sliding window analysis
    positions = df['position'].values #Get all positions (already sorted previously)
    changes = df['change'].values #Get all methylation change directions

    last_pos = positions[-1] #Get the last position
    step_size = int(window_size * step_ratio) #Calculate the step size for the window movement

    # Generate the start positions for all windows
    window_starts = np.arange(0, last_pos + step_size, step_size) #Adding step_size here ensures the last interval covers the last_pos site

    print(f"Window configuration: size={window_size}, step={step_size}, number of windows={len(window_starts)}")

    #Create a list to store the start and end positions, counts of different methylation change directions, and total significant sites for each subsequent interval
    results = []

    for i, start in enumerate(window_starts): # 'start' is assigned the start position of each window minus 1 (e.g., 0 for the first iteration)
        if i % 1000 == 0:  # Progress indicator
            print(f"Processing progress: {i}/{len(window_starts)}")

        end = start + window_size  # Calculate the end point of the window based on the window width

        # Use numpy's searchsorted for fast lookup (equivalent to binary search)
        left_idx = np.searchsorted(positions, start, side='right') #This means finding the index of the first element in the 'positions' array that is greater than 'start'
        right_idx = np.searchsorted(positions, end, side='right') #Similarly, find the index of the first element in the 'positions' array that is greater than 'end'

        # Sites within the window
        window_changes = changes[left_idx:right_idx] # Obtain an array composed of 'change' values for sites where 'pos' is within [start, end) with a width of window_size
                                                    # Note that left_idx and right_idx here are both indices
                            #The range left_idx:right_idx contains elements where the pos value is > start and <= end,
                         # however, the total list length is generally not window_size; it will be much smaller because many positions likely lack data
        # Count the quantities of different categories
        num_change_1 = np.sum(window_changes == 1)   # Count hypermethylated sites
        num_change_0 = np.sum(window_changes == 0)  # Count hypomethylated sites, or use len(window_changes) - num_change_1

        results.append({
            'window_start': start + 1,  # start + 1 is the true starting position of each window
            'window_end': end,   #Exactly window_size positions are counted from start + 1 to end
            'count_change_1': num_change_1,
            'count_change_0': num_change_0,
            'total_count': num_change_1 + num_change_0 # The methylation change directions of all sites in the interval corresponding to the current row (mutant vs. wild type)
                                                # It is also the total number of significant sites in the interval corresponding to the current row, since each test has one change direction
                                            #  And because it is read from the DMP file, all are significant
        })

    # Convert to DataFrame
    sliding_results = pd.DataFrame(results)

    # 3. Standardization processing
    max_count = sliding_results['total_count'].max()  # The maximum total number of significant sites among all intervals
    if max_count == 0:
        max_count = 1  # Avoid division by zero

    standardized_results = sliding_results.copy()
    standardized_results['standardized_count'] = sliding_results['total_count'] / max_count #Calculate the ratio of the total significant sites in the current interval to
                                                                                        # the maximum significant sites among all intervals
    # Select the required columns for standardized output
    standardized_results = standardized_results[[
        'window_start', 'window_end', 'total_count', 'standardized_count'
    ]]

    print(f"Sliding window analysis complete, generated a total of {len(sliding_results)} windows")
    print(f"Maximum count: {max_count}")

    # 4. Save files
    if save_files:
        # Sliding window results
        sliding_file = f"slidingW_{output_identifier}.txt"
        sliding_file = os.path.join(outputdir1, sliding_file)
        sliding_results[['window_start', 'window_end', 'count_change_1', 'count_change_0']].to_csv(
            sliding_file,
            sep='\t',
            index=False,
            header=False
        )

        # Standardized results
        std_file = f"noTitle_allDMCs_new_Standardized_slidingW_{output_identifier}.txt"
        std_file = os.path.join(outputdir1, std_file)
        standardized_results.to_csv(
            std_file,
            sep='\t',
            index=False,
            header=False,
            float_format='%.6f'  # Control floating-point precision
        )

        print(f"Results saved:")
        print(f"  Sliding window: {sliding_file}")
        print(f"  Standardized results: {std_file}")

    return sliding_results, standardized_results
    # The former format is window_start, window_end, count_change_1, count_change_0, total_count
                        #The latter lacks count_change but includes standardized_count, which is the ratio of significant sites in each interval to the maximum number of significant sites across all intervals

def process_common_sites_sliding_window(common_sites_df=None,
                                        window_size=1000000,
                                        step_ratio=0.05,
                                        methytype='CpG',
                                        work_dir="."):
    """
    Perform sliding window analysis on common significant sites
    Parameter：
    common_sites_df : Common significant site data; automatically loaded if None
    window_size : Sliding window size
    step_ratio : Step ratio
    methytype : Methylation type
    """

    print(f"\nStarting sliding window analysis on common significant sites...")

    and_output_dir = os.path.join(work_dir, "and_output")
    os.makedirs(and_output_dir, exist_ok=True)

    # 1. Load common significant site data
    if common_sites_df is None:
        common_file = os.path.join(and_output_dir, f"{methytype}-final_significant_sites_DMPs.txt")
        if not os.path.exists(common_file):
            print(f"Error: Common significant sites file does not exist {common_file}")
            return None
        common_sites_df = pd.read_csv(common_file, sep='\t')
        print(f"Loaded common significant sites from file: {len(common_sites_df)} sites")

    if common_sites_df.empty:
        print("No common significant site data")
        return None

    # 2. Group by chromosome
    results = {}

    # Process grouped by chromosome
    chr_groups = common_sites_df.groupby('Chromosome') #Get an iterator to iteratively obtain each chromosome number and its corresponding sub-dataframe

    for chr_num, chr_data in chr_groups: # Iteratively obtain each chromosome number and its corresponding sub-dataframe
        print(f"\n  Processing chromosome {chr_num}: {len(chr_data)} sites")

        # Prepare data for sliding window analysis, maintaining the same format as the all_simple_chr file to ensure smooth processing
        window_data = pd.DataFrame({
            'position': chr_data['Position'].astype(int),
            'pvalue': chr_data['Sig_Mean_Qvalue'],
            'change': chr_data['Methylation_Change']
        })

        # Sort
        window_data = window_data.sort_values('position').reset_index(drop=True)

        # Execute sliding window analysis by calling the previous function
        try:
            sliding_results, std_results = sliding_window_analysis(
                window_data,
                window_size=window_size,
                step_ratio=step_ratio,
                save_files=True,
                output_identifier=f"common_sites_{methytype}_Chr{chr_num}",
                outputdir1=and_output_dir
            )

            # Collect result files
            results[chr_num] = {
                'sliding_results': sliding_results,
                'standardized_results': std_results,
                'input_data': window_data
            }

            print(f"    Completed chromosome {chr_num}: {len(sliding_results)} windows")

        except Exception as e:
            print(f"    Error: Failed to process chromosome {chr_num} : {e}")
            continue

    print(f"Results are saved in files starting with'common_sites_{methytype}_Chr'")

    return results

def find_max_total_in_outputs(output_dirs, methylation_type):
    """
    Find the maximum total value of the specified methylation type across all output directories
    Parameters：
        output_dirs: List of output directories
        methylation_type: methylation_type: Methylation type (CpG, CHH, CHG)
    Returns：
        max_total: Maximum total value
    """
    max_total = 0

    for out_dir in output_dirs:
        mtype_dir = os.path.join(out_dir, methylation_type)
        if not os.path.exists(mtype_dir):
            continue

        # Find all DMP standardized files
        std_files = glob.glob(os.path.join(mtype_dir, "noTitle_allDMCs_new_Standardized_slidingW_DMP_*.txt"))

        for std_file in std_files:
            try:
                df = pd.read_csv(
                    std_file,
                    sep='\s+',
                    header=None,
                    names=['start', 'end', 'total', 'normalized']
                )
                if not df.empty:
                    file_max = df['total'].max()
                    if file_max > max_total:
                        max_total = file_max
            except Exception as e:
                print(f"    Warning: Error reading {std_file} : {e}")
                continue

    print(f" The maximum total value across all chromosomes for the {methylation_type} type is: {max_total}")
    return max_total


def plot_methylation_sliding_windows(output_dir=None, chr_series=None,work_dir="."):
    """
    Visualize all sliding window results and save them to disk
    Standardize using the global max_total to make different chromosomes visually comparable
    Parameter：
    output_dir : Specify the output directory; if None, automatically scan all output_x_y directories
    chr_series: Chromosome mapping Series
    """

    matplotlib.use('Agg')  # Use non-interactive backend
    #plt.rcParams['font.sans-serif'] = ['SimHei']
    plt.rcParams['axes.unicode_minus'] = False

    print("\nStarting to generate methylation sliding window visualization charts")

    # Scan all output directories
    if output_dir is None:
        output_dirs = glob.glob(os.path.join(work_dir, "output_*_*"))
        output_dirs = [d for d in output_dirs if os.path.isdir(d)]
    else:
        output_dirs = [output_dir]

    if not output_dirs:
        print("No output directories found")
        return

    total_plots = 0
    methylation_types = ['CpG', 'CHH', 'CHG']

    # Calculate the global max_total separately for each methylation type
    max_totals = {}
    for mtype in methylation_types:
        max_totals[mtype] = find_max_total_in_outputs(output_dirs, mtype)
        if max_totals[mtype] == 0:
            print(f"  Warning: No valid total value found for {mtype} type, using 1 as the default value")
            max_totals[mtype] = 1

    for out_dir in output_dirs:
        print(f"\nProcessing directory: {out_dir}")

        for mtype in methylation_types:
            mtype_dir = os.path.join(out_dir, mtype)
            if not os.path.exists(mtype_dir):
                continue

            print(f"  Processing methylation type: {mtype} (using global max_total={max_totals[mtype]})")

            # Find all DMP sliding window files under the current methylation directory
            dmp_sliding_files = glob.glob(os.path.join(mtype_dir, "slidingW_DMP_*.txt"))

            # Process grouped by prefix; all chromosomes with the same prefix are plotted on a single large figure
            prefix_groups = {}
            for dmp_sliding_file in dmp_sliding_files:
                basename = os.path.basename(dmp_sliding_file)
                match = re.search(r'slidingW_DMP_(.+)_Chr(\d+)\.txt$', basename)
                if not match:
                    continue

                prefix = match.group(1)  # replicate1_replicate2
                chr_num = match.group(2)

                if prefix not in prefix_groups:
                    prefix_groups[prefix] = []
                prefix_groups[prefix].append(chr_num)

            # Process each prefix group
            for prefix, chr_nums in prefix_groups.items():
                # Sort chromosome numbers
                chr_nums = sorted(chr_nums, key=lambda x: int(x))

                # Collect data for all chromosomes
                all_chrom_data = []
                chrom_names = []

                for chr_num in chr_nums:
                    # Construct the corresponding standardized file paths
                    chr_name = get_chr_name(chr_num, chr_series)

                    # Correction: Use correct variable names to construct file paths
                    dmp_sliding_file = os.path.join(mtype_dir, f"slidingW_DMP_{prefix}_Chr{chr_num}.txt")
                    dmp_std_file = os.path.join(mtype_dir,
                                                f"noTitle_allDMCs_new_Standardized_slidingW_DMP_{prefix}_Chr{chr_num}.txt")
                    ndmp_std_file = os.path.join(mtype_dir,
                                                 f"noTitle_allDMCs_new_Standardized_slidingW_N-DMP_{prefix}_Chr{chr_num}.txt")

                    if not all(os.path.exists(f) for f in [dmp_sliding_file, dmp_std_file, ndmp_std_file]):
                        print(f"    Warning: Files for Chr{chr_num} are incomplete, skipping")
                        continue

                    try:
                        # Read data
                        sliding_df = pd.read_csv(
                            dmp_sliding_file,
                            sep='\s+',
                            header=None,
                            names=['start', 'end', 'hyper', 'hypo']
                        )

                        dmp_std_df = pd.read_csv(
                            dmp_std_file,
                            sep='\s+',
                            header=None,
                            names=['start', 'end', 'total', 'normalized']
                        )

                        ndmp_std_df = pd.read_csv(
                            ndmp_std_file,
                            sep='\s+',
                            header=None,
                            names=['start', 'end', 'total', 'ndmp_normalized']
                        )

                        if sliding_df.empty or dmp_std_df.empty or ndmp_std_df.empty:
                            print(f"   Warning: Data for Chr{chr_num} is empty, skipping")
                            continue

                        # Recalculate ratios using the global max_total
                        max_total = max_totals[mtype]

                        x = (sliding_df['start'] + sliding_df['end']) / 2

                        # Recalculate all ratios (using the global max_total)
                        y_dmp = (sliding_df['hyper'] + sliding_df['hypo']) / max_total
                        y_hyper = sliding_df['hyper'] / max_total
                        y_hypo = sliding_df['hypo'] / max_total
                        y_ndmp = ndmp_std_df['ndmp_normalized']  # NDMP remains unchanged

                        # Handle possible length inconsistencies
                        max_len = max(len(x), len(y_dmp), len(y_hyper), len(y_hypo), len(y_ndmp))
                        x = x.reindex(range(max_len), fill_value=0)
                        y_dmp = y_dmp.reindex(range(max_len), fill_value=0)
                        y_hyper = y_hyper.reindex(range(max_len), fill_value=0)
                        y_hypo = y_hypo.reindex(range(max_len), fill_value=0)
                        y_ndmp = y_ndmp.reindex(range(max_len), fill_value=0)

                        # Store data
                        all_chrom_data.append({
                            'x': x,
                            'y_dmp': y_dmp,
                            'y_hyper': y_hyper,
                            'y_hypo': y_hypo,
                            'y_ndmp': y_ndmp
                        })
                        chrom_names.append(chr_name)

                        print(f"   Successfully loaded data for chromosome {chr_name}")

                    except Exception as e:
                        print(f"    Error processing Chr{chr_num}: {e}")
                        continue

                # If there is data, plot the large figure
                if all_chrom_data:
                    try:
                        # Create a large figure, one subplot per chromosome
                        n_chromosomes = len(all_chrom_data)
                        fig, axes = plt.subplots(n_chromosomes, 1, figsize=(15, 3 * n_chromosomes))

                        # If there is only one chromosome, axes is not an array and needs to be converted into one
                        if n_chromosomes == 1:
                            axes = [axes]

                        # Set the main title
                        fig.suptitle(f'{mtype} Methylation Analysis - {prefix} (Global Normalized)',
                                     fontsize=16, fontfamily='DejaVu Sans')

                        # Plot the subplot for each chromosome
                        for idx, (chrom_data, chrom_name) in enumerate(zip(all_chrom_data, chrom_names)):
                            ax = axes[idx]

                            # Plot all data lines
                            ax.plot(chrom_data['x'], chrom_data['y_dmp'], label='DMP', color='red', linewidth=1.5)
                            ax.plot(chrom_data['x'], chrom_data['y_hyper'], label='Hyper-ratio', color='green',
                                    linewidth=1)
                            ax.plot(chrom_data['x'], chrom_data['y_hypo'], label='Hypo-ratio', color='blue',
                                    linewidth=1)
                            ax.plot(chrom_data['x'], chrom_data['y_ndmp'], label='NDMP', color='darkgray',
                                    linewidth=1.5)

                            ax.set_ylim(bottom=0)
                            ax.set_ylabel('Ratio', fontsize=10, fontfamily='DejaVu Sans')

                            # Set the subplot title
                            ax.set_title(f'{chrom_name}', fontsize=18, fontfamily='DejaVu Sans', pad=20,y=-0.4)

                            # Add grid
                            ax.grid(True, alpha=0.3)

                            # Add the complete legend only to the first subplot
                            if idx == 0:
                                ax.legend(loc='upper right', ncol=2, fontsize=8, framealpha=0.7)

                        # Adjust layout
                        plt.tight_layout(rect=[0, 0, 1, 0.96])  # Leave space for the main title

                        # Save the figure
                        plot_filename = os.path.join(mtype_dir,
                                                     f"methylation_plot_{mtype}_{prefix}_all_chromosomes.png")
                        plt.savefig(plot_filename, dpi=300, bbox_inches='tight')
                        plt.close()

                        total_plots += 1
                        print(f"    Successfully generated large figure: {mtype}_{prefix} -> {os.path.basename(plot_filename)}")

                    except Exception as e:
                        print(f"    Error plotting large figure: {e}")
                        continue

    print(f"\nChart generation complete! Generated a total of {total_plots} large figures")


def plot_common_sites_sliding_windows(methytype='CpG', chr_series=None, work_dir="."):
    """
    Visualize the sliding window results of common significant sites and save them to disk
    Standardize using the global max_total and plot all chromosomes on a single large figure
    """
    matplotlib.use('Agg')

    plt.rcParams['font.family'] = 'DejaVu Sans'
    plt.rcParams['axes.unicode_minus'] = False

    print(f"\nStarting to generate sliding window visualization charts for common significant sites ({methytype})...")

    and_output_dir = os.path.join(work_dir, "and_output")

    # First find the global max_total for this methylation type
    print(f"  Looking for the global max_total of {methytype}...")
    max_total = 0
    std_files = glob.glob(
        os.path.join(and_output_dir, f"noTitle_allDMCs_new_Standardized_slidingW_common_sites_{methytype}_Chr*.txt"))

    for std_file in std_files:
        try:
            df = pd.read_csv(std_file, sep='\s+', header=None, names=['start', 'end', 'total', 'normalized'])
            if not df.empty:
                file_max = df['total'].max()
                if file_max > max_total:
                    max_total = file_max
        except Exception as e:
            print(f"    Warning: Error reading {std_file} : {e}")

    if max_total == 0:
        print(f" Warning: No valid total value found, using default value 1")
        max_total = 1
    else:
        print(f" The global max_total for {methytype} is: {max_total}")

    # Read DMR data - select the corresponding DMR file based on the methylation type
    print("  Reading DMR data...")
    dmr_file = os.path.join(and_output_dir, f"{methytype}-final_significant_regions_DMRs.txt")
    dmr_data = {}
    if os.path.exists(dmr_file):
        try:
            dmr_df = pd.read_csv(dmr_file, sep='\s+')

            for _, row in dmr_df.iterrows():
                try:
                    chrom = str(row['Chromosome'])  # Access based on column name
                    direction = int(row['Direction'])  # Access based on column name
                    start = int(row['DMR_start'])  # Access based on column name
                    end = int(row['DMR_end'])  # Access based on column name

                    # Calculate the midpoint
                    mid = (start + end) / 2

                    # Extract the numeric part of the chromosome
                    chrom_num = str(chrom).replace('Chr', '').replace('chr', '')

                    if chrom_num not in dmr_data:
                        dmr_data[chrom_num] = []

                    dmr_data[chrom_num].append((mid, direction))
                except (ValueError, IndexError):
                    continue
            print(f"  Successfully loaded {sum(len(dmrs) for dmrs in dmr_data.values())} DMRs")
        except Exception as e:
            print(f"  Error reading DMR file: {e}")
    else:
        print(f"  Warning: DMR file {dmr_file} does not exist")

    # Dynamically obtain the chromosome list
    sliding_files = glob.glob(os.path.join(and_output_dir, f"slidingW_common_sites_{methytype}_Chr*.txt"))

    # Extract chromosome numbers from filenames and sort
    chromosomes = []
    for file in sliding_files:
        match = re.search(r'slidingW_common_sites_.+_Chr(\d+)\.txt$', os.path.basename(file))
        if match:
            chr_num = match.group(1)
            if chr_num not in chromosomes:
                chromosomes.append(chr_num)

    # Sort chromosomes in numerical order
    chromosomes.sort(key=int)

    if not chromosomes:
        print(f"No sliding window files found for common significant sites({methytype})")
        return

    print(f"Found sliding window files for {len(chromosomes)} chromosomes")

    all_chrom_data = []
    chrom_nums1 = []

    for chr_num in chromosomes:
        sliding_file = os.path.join(and_output_dir, f"slidingW_common_sites_{methytype}_Chr{chr_num}.txt")
        std_file = os.path.join(and_output_dir,
                                f"noTitle_allDMCs_new_Standardized_slidingW_common_sites_{methytype}_Chr{chr_num}.txt")

        if not all(os.path.exists(f) for f in [sliding_file, std_file]):
            print(f"    Warning:Files for Chr{chr_num} are incomplete, skipping")
            continue

        try:
            # Read data
            sliding_df = pd.read_csv(sliding_file, sep='\s+', header=None, names=['start', 'end', 'hyper', 'hypo'])
            std_df = pd.read_csv(std_file, sep='\s+', header=None, names=['start', 'end', 'total', 'normalized'])

            if sliding_df.empty or std_df.empty:
                print(f"    Warning: Data for Chr{chr_num} is empty, skipping")
                continue

            if len(sliding_df) != len(std_df):
                print(f"    Warning: Data lengths for Chr{chr_num} are inconsistent, skipping")
                continue

            # Added: Read NDMP data from output_wt1_mut1 (Note: Adjusted to match your directory naming convention 'output_wt1_mut1' for clarity)
            ndmp_file = os.path.join(work_dir, "output_wt1_mut1", methytype,
                                     f"noTitle_allDMCs_new_Standardized_slidingW_N-DMP_wt_replicate1_mut_replicate1_Chr{chr_num}.txt")
            y_ndmp = None  # Initialize as None
            if os.path.exists(ndmp_file):
                try:
                    ndmp_df = pd.read_csv(ndmp_file, sep='\s+', header=None,
                                          names=['start', 'end', 'total', 'ndmp_normalized'])
                    if not ndmp_df.empty:
                        y_ndmp = ndmp_df['ndmp_normalized']
                        print(f"    Successfully read NDMP data from output_wt1_mut1: Chr{chr_num}")
                except Exception as e:
                    print(f"    Warning: Failed to read output_wt1_mut1 NDMP data (Chr{chr_num}): {e}")
            else:
                print(f"   Note: output_wt1_mut1 NDMP file does not exist (Chr{chr_num})")

            # Recalculate ratios using the global max_total
            x = (sliding_df['start'] + sliding_df['end']) / 2
            y_total = (sliding_df['hyper'] + sliding_df['hypo']) / max_total
            y_hyper = sliding_df['hyper'] / max_total
            y_hypo = sliding_df['hypo'] / max_total

            max_len = max(len(x), len(y_total), len(y_hyper), len(y_hypo))
            if y_ndmp is not None:
                max_len = max(max_len, len(y_ndmp))
            x = x.reindex(range(max_len), fill_value=0)
            y_total = y_total.reindex(range(max_len), fill_value=0)
            y_hyper = y_hyper.reindex(range(max_len), fill_value=0)
            y_hypo = y_hypo.reindex(range(max_len), fill_value=0)
            if y_ndmp is not None:
                y_ndmp = y_ndmp.reindex(range(max_len), fill_value=0)

            # Store data
            chr_real_name = get_chr_name(chr_num, chr_series)
            all_chrom_data.append({
                'x': x,
                'y_total': y_total,
                'y_hyper': y_hyper,
                'y_hypo': y_hypo,
                'y_ndmp': y_ndmp
            })
            chrom_nums1.append(chr_num)
            print(f"    Successfully loaded data for chromosome {chr_real_name} ")

        except Exception as e:
            print(f"    Error processing Chr{chr_num} : {e}")
            continue

    # If there is data, plot the large figure
    if all_chrom_data:
        try:
            # Create a large figure containing multiple subplots
            n_chromosomes = len(all_chrom_data)
            fig, axes = plt.subplots(n_chromosomes, 1, figsize=(15, 3 * n_chromosomes))

            # If there is only one chromosome, axes is not an array and needs to be converted to one
            if n_chromosomes == 1:
                axes = [axes]

            # Set the main title
            fig.suptitle(f'Distribution of Common Significant Sites - {methytype} context',
                         fontsize=16, fontfamily='DejaVu Sans')

            # Plot the subplot for each chromosome
            for idx, (chrom_data, chrom_num1) in enumerate(zip(all_chrom_data, chrom_nums1)):
                ax = axes[idx]

                # Plot DMP data
                ax.plot(chrom_data['x'], chrom_data['y_total'], label='DMP', color='red', linewidth=2)
                ax.plot(chrom_data['x'], chrom_data['y_hyper'], label='hyper-methylation', color='green', linewidth=1.5)
                ax.plot(chrom_data['x'], chrom_data['y_hypo'], label='hypo-methylation', color='blue', linewidth=1.5)


                if chrom_data['y_ndmp'] is not None:
                    ax.plot(chrom_data['x'], chrom_data['y_ndmp'], label='NDMP',
                            color='darkgray', linewidth=1.5)
                # Uniformly set the Y-axis range from 0 to 1.2
                ax.set_ylim(0, 1.2)

                # Get the chromosome number (extracted from chrom_name)
                chrom_num1 = chrom_num1.replace('chr', '').replace('Chr', '')

                if idx == 0:  # Print only once during the first subplot iteration
                    print(f"  Debug: dmr_data keys = {list(dmr_data.keys())}")
                print(f"    {get_chr_name(chrom_num1,chr_series)} -> chrom_num = '{chrom_num1}', in dmr_data: {chrom_num1 in dmr_data}")

                # Add DMR markers
                if chrom_num1 in dmr_data:
                    for mid, direction in dmr_data[chrom_num1]:
                        # Select color based on direction: 1=hyper (green), 0=hypo (blue)
                        color = 'green' if direction == 1 else 'blue'
                        # Add a vertical line at the midpoint of the DMR, displayed in the Y-axis range of 1.0 to 1.2
                        ax.axvline(x=mid, ymin=0.9, ymax=1, color=color, linewidth=2, alpha=0.7)

                # Set subplot title and labels
                ax.text(0.5, -0.2, f"{get_chr_name(chrom_num1,chr_series)}",
                        transform=ax.transAxes,
                        fontfamily='DejaVu Sans',
                        ha='center', va='top',
                        fontsize=15)

                # Add grid
                ax.grid(True, alpha=0.3)

                # Add the complete legend only to the first subplot
                if idx == 0:
                    # Create custom legend entries, including DMR markers
                    from matplotlib.lines import Line2D
                    legend_elements = [
                        Line2D([0], [0], color='red', linewidth=2, label='DMP'),
                        Line2D([0], [0], color='green', linewidth=1.5, label='hyper-DMP'),
                        Line2D([0], [0], color='blue', linewidth=1.5, label='hypo-DMP'),
                        Line2D([0], [0], color='green', linewidth=2, label='hyper-DMR'),
                        Line2D([0], [0], color='blue', linewidth=2, label='hypo-DMR')
                    ]
                    # If NDMP data exists, add it to the legend
                    if chrom_data['y_ndmp'] is not None:
                        legend_elements.append(
                            Line2D([0], [0], color='darkgray', linewidth=1.5,
                                  label='NDMP')
                        )
                    ax.legend(handles=legend_elements, loc='upper right', bbox_to_anchor=(1.0, 0.95),
                              ncol=5, fontsize=8, framealpha=0.7)

            # Adjust layout
            plt.tight_layout(rect=[0, 0, 1, 0.96])  # Leave space for the main title

            # Save the figure
            plot_filename = os.path.join(and_output_dir, f"common_sites_plot_{methytype}_all_chromosomes.png")
            plt.savefig(plot_filename, dpi=300, bbox_inches='tight')
            plt.close()

            print(f"Successfully generated large figure: {methytype} -> {os.path.basename(plot_filename)}")

        except Exception as e:
            print(f"绘制大图时出错: {e}")
    else:
        print(f"No valid data found for {methytype} ")



def rename_chromosome_files(chr_series, work_dir="."):
    """
    Batch rename all files containing chromosome numbers, replacing 'Chr' numbers with actual chromosome names

    Parameters：
        chr_series: Chromosome mapping Series (Chromosome name -> index)
        work_dir: Working directory
    """
    print("\nStarting batch renaming of chromosome files...")

    if chr_series is None or len(chr_series) == 0:
        print("Error: Chromosome mapping is empty, skipping renaming")
        return

    # Create reverse mapping: numeric index -> chromosome name
    # The index of chr_series is the chromosome name, and the value is the numeric index
    index_to_chr = {i: chr_series.index[i] for i in range(len(chr_series))}

    print(f"Chromosome mapping: {index_to_chr}")

    # Define the directory patterns to search
    search_dirs = []

    # Add all output_x_y directories
    output_dirs = glob.glob(os.path.join(work_dir, "output_*_*"))
    search_dirs.extend([d for d in output_dirs if os.path.isdir(d)])

    # Add the and_output directory
    and_output_dir = os.path.join(work_dir, "and_output")
    if os.path.exists(and_output_dir):
        search_dirs.append(and_output_dir)

    if not search_dirs:
        print("No directories found to process")
        return

    print(f"Will search for files in {len(search_dirs)} directories")

    # Statistics
    total_renamed = 0
    failed_renames = 0

    # Define the regular expression pattern to match chromosome numbers
    # Match the pattern of 'Chr' followed by digits, such as Chr1, Chr12, etc.
    chr_pattern = re.compile(r'(.*?)Chr(\d+)(.*?)$')

    # Iterate through all directories
    for search_dir in search_dirs:
        print(f"\nProcessing directory: {search_dir}")

        # Recursively traverse all files in the directory
        for root, dirs, files in os.walk(search_dir):
            for filename in files:
                # Check if the filename contains the 'Chr' number pattern
                match = chr_pattern.match(filename)

                if match:
                    prefix = match.group(1)  # The part before 'Chr'
                    chr_num = int(match.group(2))  # The chromosome number
                    suffix = match.group(3)  # The part after the 'Chr' number

                    # Obtain the actual chromosome name based on chr_series
                    # Note: Chr1 in the filename corresponds to index 0
                    chr_index = chr_num - 1

                    if chr_index not in index_to_chr:
                        print(f"  Warning: Chr{chr_num} is not in the mapping table, skipping file {filename}")
                        continue

                    real_chr_name = index_to_chr[chr_index]

                    # Construct the new filename
                    # If the chromosome name itself contains the 'chr' prefix, use it directly
                    # Otherwise, use the 'Chr' prefix
                    if real_chr_name.lower().startswith('chr'):
                        chr_part = real_chr_name
                    else:
                        chr_part = f"Chr{real_chr_name}"

                    new_filename = f"{prefix}{chr_part}{suffix}"

                    # If the old and new filenames are the same, skip
                    if filename == new_filename:
                        continue

                    # Construct full path
                    old_path = os.path.join(root, filename)
                    new_path = os.path.join(root, new_filename)

                    # Check if target file already exists
                    if os.path.exists(new_path):
                        print(f" Warning: Target file already exists, skipping renaming: {filename} -> {new_filename}")
                        failed_renames += 1
                        continue

                    # Execute rename
                    try:
                        os.rename(old_path, new_path)
                        total_renamed += 1
                        print(f"  {filename} -> {new_filename}")
                    except Exception as e:
                        print(f" Rename failed: {filename} -> {new_filename}, Error: {e}")
                        failed_renames += 1

    # Output statistics
    print(f"\nRenaming complete!")
    print(f"  Successfully renamed: {total_renamed} files")
    if failed_renames > 0:
        print(f" Failed or skipped: {failed_renames} files")


def convert_output_to_csv(work_dir="."):
    """
    Convert the final DMP and final DMR files in the and_output directory to CSV format (comma-separated)

    Parameters：
        work_dir: Working directory
    """

    and_output_dir = os.path.join(work_dir, "and_output")

    if not os.path.exists(and_output_dir):
        print(f"Error: Directory {and_output_dir} does not exist")
        return 0

    # Define file patterns to convert
    file_patterns = [
        "*-final_significant_sites_DMPs.txt",  # final DMP files
        "*-final_significant_regions_DMRs.txt"  # final DMR files
    ]

    converted_count = 0
    failed_count = 0

    for pattern in file_patterns:
        # Search for matching files
        matching_files = glob.glob(os.path.join(and_output_dir, pattern))

        for txt_file in matching_files:
            try:
                # Read tab-separated files
                df = pd.read_csv(txt_file, sep=r'\s+')

                if df.empty:
                    print(f"  Skipping empty file: {os.path.basename(txt_file)}")
                    continue

                # Generate CSV file path (replace .txt with .csv)
                csv_file = txt_file.replace('.txt', '.csv')

                # Save as a comma-separated CSV file
                df.to_csv(csv_file, sep=',', index=False)

                print(f"  Successfully converted: {os.path.basename(txt_file)} -> {os.path.basename(csv_file)}")
                converted_count += 1

            except Exception as e:
                print(f" Conversion failed: {os.path.basename(txt_file)}, Error: {e}")
                failed_count += 1

    print(f"\\nConversion complete!")
    print(f"  Successfully converted: {converted_count} files")
    if failed_count > 0:
        print(f"  Conversion failed: {failed_count} files")

    return converted_count


def convert_chromosome_to_names(chr_series, work_dir="."):
    """
    Convert the Chromosome column in the final DMP and final DMR files within the and_output directory
    from numeric values to actual chromosome names

    Parameter：
        chr_series: Chromosome mapping Series (Chromosome name -> index, 0-based index)

    Returns：
        Number of successfully converted files
    """

    if chr_series is None or len(chr_series) == 0:
        print("Error: Chromosome mapping is empty, skipping conversion")
        return 0

    and_output_dir = os.path.join(work_dir, "and_output")

    if not os.path.exists(and_output_dir):
        print(f"Error:Directory {and_output_dir} does not exist")
        return 0

    # Create mapping from numeric values to chromosome names
    # The Chromosome column in the file contains 1-based numbers, which need conversion to chromosome names in chr_series
    # chr_series.index[0] corresponds to 1 in the file
    index_to_chr = {i + 1: chr_series.index[i] for i in range(len(chr_series))}

    print(f"Chromosome mapping table: {index_to_chr}")

    # Define file patterns to process
    file_patterns = [
        "*-final_significant_sites_DMPs.txt",  # final DMP files
        "*-final_significant_regions_DMRs.txt"  # final DMR files
    ]

    converted_count = 0
    failed_count = 0

    for pattern in file_patterns:
        # Search for matching files
        matching_files = glob.glob(os.path.join(and_output_dir, pattern))

        for file_path in matching_files:
            try:
                # Read file
                df = pd.read_csv(file_path, sep=r'\s+')

                if df.empty:
                    print(f" Skipping empty file: {os.path.basename(file_path)}")
                    continue

                # Check if Chromosome column exists
                if 'Chromosome' not in df.columns:
                    print(f"  Warning: No Chromosome column in {os.path.basename(file_path)} , skipping")
                    continue

                # Save the original Chromosome column for debugging
                original_chrs = df['Chromosome'].unique()

                # Convert Chromosome column
                # Ensure integer type first
                df['Chromosome'] = df['Chromosome'].astype(int)

                # Use mapping to convert to chromosome names
                df['Chromosome'] = df['Chromosome'].map(index_to_chr)

                # Check for values that failed to map
                if df['Chromosome'].isna().any():
                    unmapped_count = df['Chromosome'].isna().sum()
                    print(f"  Warning:{unmapped_count} chromosome numbers failed to map in {os.path.basename(file_path)}")
                    # Optional: Drop rows that failed to map
                    df = df.dropna(subset=['Chromosome'])

                # Save back to the original file (overwrite)
                df.to_csv(file_path, sep='\t', index=False)

                print(f"    Successfully converted: {os.path.basename(file_path)}")
                print(f"    Original numbers: {sorted(original_chrs)}")
                print(f"    After conversion: {sorted(df['Chromosome'].unique())}")
                converted_count += 1

            except Exception as e:
                print(f"    Conversion failed: {os.path.basename(file_path)}")
                print(f"    Error information: {e}")
                import traceback
                traceback.print_exc()
                failed_count += 1

    # Output statistics
    print(f"\nChromosome number conversion complete!")
    print(f"  Successfully converted: {converted_count} files")
    if failed_count > 0:
        print(f" Conversion failed: {failed_count} files")

    return converted_count


def main():
    start_time = time.time()
    print("Program Description:")
    print("Perform m*n*3 Fisher's exact tests")
    print("Each group includes 3 methylation types: CpG, CHH, CHG")
    print()

    try:
        # Get command line arguments
        # m = int(sys.argv[2])
        # n = int(sys.argv[1])
        # dir1 = sys.argv[4]
        # dir2 = sys.argv[3]
        # biotype = int(sys.argv[5])
        dir2 = input("Please enter the wild-type sample directory name:").strip()
        dir1 = input("Please enter the mutant-type sample directory name:").strip()
        n = int(input("Please enter the number of wild-type groups (n):"))
        m = int(input("Please enter the number of mutant-type groups (m):"))
        biotype = int(input("Please enter the biological type of the given genotypes (0-Animal, 1-Plant, 2-No filter):"))

    # Since directory names can be entered freely without triggering errors, a ValueError here definitely means the group numbers m or n were entered incorrectly
    except ValueError:
        print("Error: Please enter valid numbers")
        sys.exit(1)

    # Verify directories exist and ensure group numbers are valid
    if not os.path.exists(dir1):
        print(f"Error:Directory '{dir1}' does not exist!")
        sys.exit(1)
    if not os.path.exists(dir2):
        print(f"Error: Directory '{dir2}' does not exist!")
        sys.exit(1)
    if m <= 0 or n <= 0:
        print("Error: Group numbers must be greater than 0")
        sys.exit(1)

    print(f"\nParameter Confirmation:")
    print(f"First directory: {dir1} (contains {m} groups of files)")
    print(f"Second directory: {dir2} (contains {n} groups of files)")
    print(f"Expected to perform: {m * n * 3} Fisher's exact tests")

    print("\nPhase 1: newtoboth in progress")
    # Convert Bismark new format data to both format
    chr_series = newtoboth(m, n, dir1, dir2)
    if biotype == 0:
        unfilter_mtypes = ["CHH","CHG"]
    elif biotype == 1:
        unfilter_mtypes = ["CpG"]
    elif biotype ==2:
        unfilter_mtypes = ["CHH","CHG","CpG"]
    else:
        print("Error: Biological type must be 0, 1, or 2")
        sys.exit(1)
    print(f"Methylation types not requiring p-value pre-filtering: {unfilter_mtypes} {unfilter_mtypes}")
    success = process_all_combinations(dir1, dir2, m, n,unfilter_mtypes)  # process_all_combinations performs m*n*3 tests

    if success: # All tests are successful
        print("\nAll tests and FDR corrections successfully completed!")
        methylation_types = ['CpG', 'CHH', 'CHG']
        for mtype in methylation_types:
            common_sites_df = find_common_significant_sites(methytype2=mtype, dir1=dir1, dir2=dir2)
                             # Its format is: 'Chromosome', 'Methylation_Type', 'Position',
                                    # 'Methylation_Change', 'Hyper_Ratio', 'Hyper_Count', 'Hypo_Count',
                                # 'Sig_Mean_Qvalue', 'Max_Qvalue', 'Min_Qvalue', 'Num_Comparisons'
            if common_sites_df is not None and not common_sites_df.empty:
                print(f"Found {len(common_sites_df)} common significant sites for {mtype} type")

                # Perform sliding window analysis
                print(f"Starting sliding window analysis for {mtype} common significant sites...")

                results = process_common_sites_sliding_window(
                    common_sites_df=common_sites_df,
                    methytype=mtype
                )

                print(f"Completed sliding window analysis for {mtype} type")
            else:
                print(f"No common significant sites found for {mtype} type")



        print("Starting DMR analysis pipeline")
        process_common_sites_dmr_and_summarize(
            dir1=dir1,
            dir2=dir2,
            m=m,
            n=n,
            methylation_types=methylation_types,
        )
        # Generate visualization charts for all sliding windows
        plot_methylation_sliding_windows(chr_series=chr_series)
        for mtype111 in ["CpG","CHH","CHG"]:
            plot_common_sites_sliding_windows(mtype111, chr_series=chr_series)
        convert_chromosome_to_names(chr_series=chr_series, work_dir=".")
        rename_chromosome_files(chr_series=chr_series, work_dir=".")
        convert_output_to_csv(work_dir=".")
        print(f"- DMP results: output_x_y/methylation_type/")
        print(f"- Common significant sites: and_output/")
        print(f"- Final significant DMRs: and_output/-final_significant_regions_DMRs.txt")
    else:
        print("\nSome tests failed. Please check the output information.")
    end_time = time.time()  # Record end time
    elapsed_time = end_time - start_time  # Calculate elapsed time
    print(f"Total time elapsed: {elapsed_time:.2f} seconds")


if __name__ == "__main__":
    main()



