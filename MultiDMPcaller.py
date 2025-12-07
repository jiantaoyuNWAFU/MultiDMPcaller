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

DmrRecord = namedtuple('DmrRecord', ['exp_methy', 'exp_unmethy', 'wild_methy', 'wild_unmethy', 'qvalue', 'direction'])

def process_common_sites_dmr_and_summarize(dir1, dir2, m, n, methylation_types=['CpG', 'CHH', 'CHG'],work_dir="."):
    """
    完整流程：处理 common sites DMR + 循环求和分析 + 汇总  

    参数：
        dir1: 第一个样本目录
        dir2: 第二个样本目录
        m: 突变型样本数
        n: 野生型样本数
        methylation_types: 甲基化类型列表
    """
    # 第一步：处理 common sites 生成 DMR
    print("第三阶段：处理共同显著位点 DMR")

    for mtype in methylation_types:
        dmr_results = process_common_sites_to_dmr(methylation_type=mtype,work_dir=work_dir)
        if dmr_results:
            print(f"\n{mtype} 类型 DMR 分析完成，共 {len(dmr_results)} 个染色体有DMR结果")

    # 第二步：对所有 output_x_y 的 DMR 进行求和分析并汇总
    summarize_all_dmr_methylation(dir1, dir2, m, n, methylation_types,work_dir=work_dir)

def generate_final_significant_dmr(dmr_data_dict, methylation_type, output_dir,m,n,dir1,dir2):
    """
    生成最终的显著 DMR 文件（基于贝叶斯判定）

    参数：
        dmr_data_dict: 染色体->DMR区域(start,end)->数据列表(一系列每个位点的统计数据) 的映射字典
        methylation_type: 甲基化类型
        output_dir: 输出目录
        threshold: 贝叶斯判定阈值（默认 2/3）
    """
    print(f"\n处理 {methylation_type} 的最终 DMR 汇总...")

    if not dmr_data_dict:
        print(f"  {methylation_type} 无 DMR 数据")
        return

    final_dmr_list = []

    # 遍历所有染色体
    for chr_num in sorted(dmr_data_dict.keys()):
        print(f"  处理染色体 {chr_num}...")
        chr_dmr_dict = dmr_data_dict[chr_num]  # 获取当前染色体中的该映射：DMR区域(start,end)->数据列表(一系列每个位点的统计数据)
                                                                    # 此处数据列表中有m*n个元素，代表每次检验中该区域的统计信息

        # 遍历该染色体的所有 DMR 区域
        for dmr_key, data_list in chr_dmr_dict.items():
            start, end = dmr_key

            # 统计显著次数
            sig_count = sum(1 for item in data_list if item.qvalue < 0.05)
            total_count = len(data_list)

            # 判定
            is_significant = bayes_deciding(sig_count, total_count - sig_count)

            if not is_significant:
                continue

            # 计算平均值
            avg_exp_m = np.mean([item.exp_methy for item in data_list])
            avg_exp_u = np.mean([item.exp_unmethy for item in data_list])
            avg_wild_m = np.mean([item.wild_methy for item in data_list])
            avg_wild_u = np.mean([item.wild_unmethy for item in data_list])

            filtered_values = [item.qvalue for item in data_list if item.qvalue < 0.05]
            avg_qvalue = np.mean(filtered_values) if filtered_values else 1

            # 投票决定 direction（多数投票）
            direction_votes = [item.direction for item in data_list]
            direction = 1 if sum(direction_votes) >= len(direction_votes) / 2 else 0

            # 计算概率
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
            # 添加所有replicate的qvalue（按顺序）
            idx = 0
            for x in range(1, m + 1):
                for y in range(1, n + 1):
                    col_name = f'qvalue_{dir2}{y}_{dir1}{x}'
                    if idx < len(data_list):
                        dmr_record[col_name] = data_list[idx].qvalue
                    else:
                        dmr_record[col_name] = 1.0
                    idx += 1
            final_dmr_list.append(dmr_record)

    if not final_dmr_list:
        print(f"  {methylation_type} 无显著 DMR")
        return

    # 转换为 DataFrame 并排序
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
    # 保存结果
    output_file = os.path.join(output_dir, f"{methylation_type}-final_significant_regions_DMRs.txt")
    final_df.to_csv(output_file, sep='\t', index=False)

    print(f"  {methylation_type} 最终显著 DMR: {len(final_df)} 个")
    print(f"  保存至: {output_file}")

    # 统计信息
    hyper_count = len(final_df[final_df['Direction'] == 1])
    hypo_count = len(final_df[final_df['Direction'] == 0])
    print(f"    - 高甲基化: {hyper_count} ({hyper_count / len(final_df) * 100:.1f}%)")
    print(f"    - 低甲基化: {hypo_count} ({hypo_count / len(final_df) * 100:.1f}%)")

    return final_df

def collect_dmr_results(methy_dir, methylation_type, all_dmr_results):
    """
    收集单次 DMR 分析的结果

    参数：
        methy_dir: 甲基化类型目录（如 ./and_output/CpG/）
        methylation_type: 甲基化类型
        all_dmr_results: 汇总字典，输入的时候只有三个甲基化类型键，后面内容是空的
    """
    # 查找所有 dmr_fisher_Chr*.txt 文件
    fisher_files = glob.glob(os.path.join(methy_dir, "dmr_fisher_Chr*.txt"))

    for fisher_file in fisher_files:
        # 提取染色体号
        match = re.search(r'Chr(\d+)\.txt$', fisher_file)
        if not match:
            continue
        chr_num = int(match.group(1))

        # 读取文件
        try:
            df = pd.read_csv(fisher_file, sep='\t')

            if df.empty:
                continue

            # 初始化该染色体的字典
            if chr_num not in all_dmr_results[methylation_type]:
                all_dmr_results[methylation_type][chr_num] = defaultdict(list)

            # 收集每个 DMR 区域的数据
            for _, row in df.iterrows():
                dmr_key = (int(row['DMR_start']), int(row['DMR_end']))

                # 存储 (exp_m, exp_u, wild_m, wild_u, qvalue, direction)
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
            print(f"    警告：读取 {fisher_file} 失败: {e}")
            continue

def summarize_all_dmr_methylation(dir1, dir2, m, n, methylation_types=['CpG', 'CHH', 'CHG'],work_dir="."):
    """
    对所有 output_x_y 目录，使用 common DMR 区域进行甲基化读段求和分析
    """
    print("第四阶段：Common DMR 在各组合中的甲基化读段求和分析")

    and_output_dir = os.path.join(work_dir, "and_output")
    os.makedirs(and_output_dir, exist_ok=True)

    # 存储所有 DMR 的甲基化数据
    all_dmr_results = {mtype: {} for mtype in methylation_types}

    # 1. 循环所有 output_x_y 组合
    for replicate_x in range(1, m + 1):
        for replicate_y in range(1, n + 1):
            print(f"\n处理组合 ({dir1}{replicate_x}, {dir2}{replicate_y})...")

            # 获取 both 文件路径
            file1_path = os.path.join(dir1, f"{replicate_x}-bothMeUnme_diffChromo_NOREPEATED_methy_sites_{{}}.txt")
            file2_path = os.path.join(dir2, f"{replicate_y}-bothMeUnme_diffChromo_NOREPEATED_methy_sites_{{}}.txt")

            # 2. 循环三种甲基化类型
            for mtype in methylation_types:
                # 检查 both 文件是否存在
                f1 = file1_path.format(mtype)
                f2 = file2_path.format(mtype)
                if not os.path.exists(f1) or not os.path.exists(f2):
                    print(f"  跳过 {mtype}：both 文件不存在")
                    continue

                print(f"  处理 {mtype} 类型...")

                # 创建输出目录（保存到 and_output 下）
                methy_output_dir = os.path.join(and_output_dir, f"dmr_analysis_wt{replicate_y}_mut{replicate_x}", mtype)
                os.makedirs(methy_output_dir, exist_ok=True)

                # 3. 调用修改后的 summarize_dmr_methylation，传入 custom_dmr_dir
                try:
                    summarize_dmr_methylation(
                        methy_dir=methy_output_dir,  # 结果保存到 and_output 下的新目录
                        replicate_x=replicate_x,
                        replicate_y=replicate_y,
                        file1_path=f1,
                        file2_path=f2,
                        methylation_type=mtype,
                        custom_dmr_dir=and_output_dir  #  读取 and_output 的 common DMR
                    )

                    # 4. 读取生成的 dmr_fisher 文件，收集结果
                    collect_dmr_results(methy_output_dir, mtype, all_dmr_results)

                except Exception as e:
                    print(f"  错误：处理 {mtype} 失败: {e}")
                    import traceback
                    traceback.print_exc()
                    continue

    # 5. 汇总所有结果，生成最终的显著 DMR
    print("第五阶段：汇总 DMR 结果并进行贝叶斯判定")

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
    读取 final_significant_sites_DMPs.txt 文件，按染色体分组处理DMR

    参数：

    methylation_type: 甲基化类型（CpG, CHH, CHG）
    返回：
        dmr_results: 字典 {chr_num: dmr_list_file_path}
    """
    print(f"\n开始处理 {methylation_type} 的共同显著位点 DMR 分析...")

    and_output_dir = os.path.join(work_dir, "and_output")
    os.makedirs(and_output_dir, exist_ok=True)

    # 1. 读取 final_significant_sites_DMPs 文件
    common_file = os.path.join(and_output_dir, f"{methylation_type}-final_significant_sites_DMPs.txt")
    if not os.path.exists(common_file):
        print(f"错误：文件不存在 {common_file}")
        return None

    try:
        df = pd.read_csv(common_file, sep='\t')
        print(f"成功读取文件，共 {len(df)} 个位点")
    except Exception as e:
        print(f"读取文件失败: {e}")
        return None

    if df.empty:
        print(f"警告：{methylation_type} 文件为空")
        return None

    # 2. 按染色体分组
    chromosomes = sorted(df['Chromosome'].unique(), key=natural_sort_key)
    print(f"找到 {len(chromosomes)} 个染色体: {chromosomes}")

    dmr_results = {}
    total_chromosomes = len(chromosomes)

    # 3. 对每个染色体进行处理
    for chr_num in chromosomes:
        print(f"\n  处理染色体 {chr_num} ({chr_num}/{total_chromosomes})...")

        # 筛选当前染色体的数据
        chr_df = df[df['Chromosome'] == chr_num].copy()

        # 提取需要的三列：Position, Sig_Mean_Qvalue, Methylation_Change
        dmp_data = chr_df[['Position', 'Sig_Mean_Qvalue', 'Methylation_Change']].copy()

        # 确保数据类型正确
        dmp_data['Position'] = dmp_data['Position'].astype(int)
        dmp_data['Sig_Mean_Qvalue'] = dmp_data['Sig_Mean_Qvalue'].astype(float)
        dmp_data['Methylation_Change'] = dmp_data['Methylation_Change'].astype(int)

        # 按位点号排序
        dmp_data = dmp_data.sort_values('Position').reset_index(drop=True)

        print(f"    染色体 {chr_num} 共 {len(dmp_data)} 个位点")

        if len(dmp_data) == 0:
            print(f"    跳过：染色体 {chr_num} 无有效位点")
            continue

        # 4. 创建临时 DMP 文件（格式：pos qvalue change）
        temp_dmp_file = os.path.join(and_output_dir,
                                     f"DMP_common_{methylation_type}_Chr{chr_num}.txt")

        # 写入 DMP 格式文件（第一行是 "first line"，后续是 pos qvalue change）
        with open(temp_dmp_file, 'w') as f:
            f.write("first line\n")
            for _, row in dmp_data.iterrows():
                f.write(f"{int(float(row['Position']))} {float(row['Sig_Mean_Qvalue'])} {int(float(row['Methylation_Change']))}\n")

        print(f"    创建 DMP 文件: {os.path.basename(temp_dmp_file)}")

        # 5. 调用 DMR 分析函数
        # 注意：run_dmr_pipeline_on_dmp_file 需要 chromoNo 参数
        # 我们传入总染色体数
        try:
            dmr_list_file = run_dmr_pipeline_on_dmp_file(
                dmp_file=temp_dmp_file,
                chromoNo=total_chromosomes
            )

            if dmr_list_file:
                dmr_results[chr_num] = dmr_list_file
                print(f"     染色体 {chr_num} DMR 分析完成")
            else:
                print(f"    染色体 {chr_num} DMR 分析失败（可能无有效DMR）")

        except Exception as e:
            print(f"     染色体 {chr_num} 处理出错: {e}")
            import traceback
            traceback.print_exc()
            continue

    print(f"\n{methylation_type} 共同显著位点 DMR 分析完成！")
    print(f"成功处理 {len(dmr_results)}/{len(chromosomes)} 个染色体")

    return dmr_results

def summarize_dmr_methylation(methy_dir, replicate_x, replicate_y, file1_path, file2_path, methylation_type='CpG', custom_dmr_dir=None):
    """
    对 DMR 区域进行甲基化读段求和，并计算 Fisher p 值、FDR q 值和甲基化方向。
    所有输出文件将保存在 methy_dir 目录下（例如：./output_1_1/CpG/）。
    注意！！！custom_dmr_dir: 自定义 DMR 文件所在目录（如果为 None，则使用 methy_dir），这样能区分每个output_x_y自己生成的dmr和common_sites生成的dmr
    """
    # 这里是只处理一次检验的
    print(f"    开始 DMR 甲基化读段求和分析...")

    n_chromosomes = get_column_count(file1_path)
    if n_chromosomes is None:
        print("    无法获取染色体数量，跳过 DMR 求和")
        return

    chromosomes = [f'Chr{i}' for i in range(1, n_chromosomes + 1)]

    # 收集所有染色体的 DMR 数据（不含 p 值）
    all_dmr_data = []  # (chrom, start, end, exp_m, exp_u, wild_m, wild_u)

    for idx, chrom in enumerate(chromosomes):
        chrom_num = idx + 1
        if custom_dmr_dir is not None:
            # 读取 common DMR 文件
            dmr_file = os.path.join(custom_dmr_dir, f"DMR_list_DMP_common_{methylation_type}_Chr{chrom_num}.txt")
        else:
            # 读取 output目录 的 DMR 文件
            dmr_file = os.path.join(methy_dir, f"DMR_list_DMP_wt_replicate{replicate_y}_mut_replicate{replicate_x}_Chr{chrom_num}.txt")
        if not os.path.exists(dmr_file):
            print(f"      跳过 {chrom}，DMR 文件不存在")
            continue

        # 读取 DMR 区域
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
            print(f"      {chrom} 无有效 DMR 区域")
            continue

        col_start = idx * 3

        # 读取实验组数据
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
            print(f"      读取实验组失败 ({chrom}): {e}")
            continue

        # 读取对照组数据
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
            print(f"      读取对照组失败 ({chrom}): {e}")
            continue

        # 汇总每个 DMR 区域的 reads
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
        print("    无任何有效 DMR 数据，跳过后续分析")
        return

    # === 第一步：输出 dmr_summary_{chrom}.txt ===
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
            print(f"      {chrom} DMR 汇总完成 → {summary_file}")

    # === 第二步：为每个 DMR 计算 p 值 ===
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

    # === 第三步：全局 FDR 校正（跨染色体） ===
    pvalues = np.array([item[7] for item in dmr_with_pvals])
    qvalues = calculate_qvalues(pvalues, pi=1.0)

    # === 第四步：按染色体组织结果，添加 direction ===
    chrom_data_dict = defaultdict(list)
    for i, (chrom, start, end, exp_m, exp_u, wild_m, wild_u, pval) in enumerate(dmr_with_pvals):
        qval = qvalues[i]

        # 计算甲基化方向
        exp_total = exp_m + exp_u
        wild_total = wild_m + wild_u
        if exp_total > 0 and wild_total > 0:
            exp_rate = exp_m / exp_total
            wild_rate = wild_m / wild_total
            direction = 1 if exp_rate > wild_rate else 0
        else:
            direction = 0  # 无法判断时设为 0

        chrom_data_dict[chrom].append((start, end, exp_m, exp_u, wild_m, wild_u, pval, qval, direction))

    # === 第五步：输出完整 Fisher 结果 + 显著子集 ===
    for chrom in chromosomes:
        if chrom not in chrom_data_dict:
            continue

        # 完整结果
        fisher_file = os.path.join(methy_dir, f"dmr_fisher_{chrom}.txt")
        with open(fisher_file, 'w') as f:
            f.write("DMR_start\tDMR_end\texp_methy_sum\texp_unmethy_sum\twild_methy_sum\twild_unmethy_sum\tpvalue\tqvalue\tdirection\n")
            for row in chrom_data_dict[chrom]:
                start, end, exp_m, exp_u, wild_m, wild_u, pval, qval, direction = row
                p_out = pval if not np.isnan(pval) else 'nan'
                q_out = qval if not np.isnan(qval) else 'nan'
                f.write(f"{start}\t{end}\t{exp_m}\t{exp_u}\t{wild_m}\t{wild_u}\t{p_out:.6g}\t{q_out:.6g}\t{direction}\n")
        print(f"      {chrom} Fisher + FDR + direction 完成 → {fisher_file}")

        # 显著结果 (q < 0.05)
        sig_rows = [
            row for row in chrom_data_dict[chrom]
            if not np.isnan(row[7]) and row[7] < 0.05
        ]
        if sig_rows:
            sig_file = os.path.join(methy_dir, f"dmr_fisher_significant_{chrom}.txt")
            with open(sig_file, 'w') as f:
                f.write("DMR_start\tDMR_end\texp_methy_sum\texp_unmethy_sum\twild_methy_sum\twild_unmethy_sum\tpvalue\tqvalue\tdirection\n")
                for row in sig_rows:
                    start, end, exp_m, exp_u, wild_m, wild_u, pval, qval, direction = row
                    f.write(f"{start}\t{end}\t{exp_m}\t{exp_u}\t{wild_m}\t{wild_u}\t{pval:.6g}\t{qval:.6g}\t{direction}\n")
            print(f"        → 显著 DMR (q<0.05): {sig_file}")
        else:
            print(f"        → {chrom} 无显著 DMR (q<0.05)")

    print("    DMR 甲基化分析完成。")

def run_dmr_pipeline_on_dmp_file(dmp_file: str, chromoNo: int = 10):
    """
    从 DMP 文件生成 DMR 区域 , 生成dmr_list文件
    - methylation_matrix_file: 用于动态获取染色体数量（如 1-bothMeUnme_...txt）
    - 所有输出文件保存在 dmp_file 所在目录。
    """
    # 总的核心思想如下：
    # 用滑动窗口找到 DMP 密集的区域
    # 通过"跳跃合并"连接相邻的密集区域
    # 最终筛选出足够长、足够显著的 DMR
    sWinN = 1000  # 滑动窗口大小（1000 bp）
    M0 = 4  # 窗口内最少 DMP 数量（至少 4 个）
    M1 = 10  # 最终 DMR 内最少 DMP 数量（至少 10 个）
    M2 = 10  # 跳跃步长，保留原始分割行为，

    # 为安全起见，确保 chromoNo >= 6（因为用到 arrayMethy1_script1[5]）
    chromoNo = max(chromoNo, 6)

    class PositionNoNode:
        def __init__(self, pos=0, end=0, pV=0.0, ratio=0.0, num=0, num2=0, numCom=0, markR=0, DMR_S=0, DMR_E=0):
            self.pos = pos  # 窗口起始位置
            self.end = end  # 窗口结束位置
            self.pV = pV
            self.ratio = ratio
            self.num = num  # 高甲基化数量
            self.num2 = num2  # 低甲基化数量
            self.numCom = numCom
            self.markR = markR
            self.DMR_S = DMR_S
            self.DMR_E = DMR_E
            self.posV = []
            self.logPV = []
            self.meUnV = []

    output_dir = os.path.dirname(dmp_file)
    base_name = os.path.basename(dmp_file)

    arrayMethy1 = [[] for _ in range(chromoNo)]  # 创建染色体个数个子列表

    # === 读取 DMP 位点 ===
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
                    tmpNode.num = change  # 以上几行建立了当前dmp行所对应的位点信息
                    arrayMethy1[0].append(tmpNode)
                except (ValueError, IndexError):
                    continue
    except Exception as e:
        print(f"    DMP 文件读取失败 {dmp_file}: {e}")
        return None

    if not arrayMethy1[0]:
        return None

    arrayMethy1[0].sort(key=lambda x: x.pos)
    firstP = arrayMethy1[0][0].pos
    lastP = arrayMethy1[0][-1].pos

    # === 滑动窗口 ===
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

    # === 构建滑动窗口列表 ===
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

    # === 标准化输出 ===
    if arrayMethy1[2]:
        maxCom = max(node.numCom for node in arrayMethy1[2])
        maxCom = max(maxCom, 1)
        std_file = os.path.join(output_dir, f"noTitle_allDMCs_new_Standardized_slidingW_{base_name}")
        with open(std_file, 'w') as f:
            for node in arrayMethy1[2]:
                if node.end <= lastP:
                    std_val = node.numCom / maxCom
                    f.write(f"{node.pos:<20}{node.end:<20}{node.numCom:<20}{std_val:.6f}\n")

    # === 跳跃合并（使用你提供的代码）===
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

        # 向左扩展
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

        # 向右扩展
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
    print(f"识别到 {len(arrayMethy1[3])} 个DMR区域")

    # 输出DMR结果和边界文件
    with open(dmr_out_file, 'w') as cout05:
        with open(boundary_file, 'w') as bound_out:
            for node in arrayMethy1[3]:
                cout05.write(f"{node.pos:<20}{node.end:<20}{node.numCom:<20}[{node.DMR_S} {node.DMR_E}]\n")
                bound_out.write(f"{node.DMR_S} {node.DMR_E}\n")

    print(f"已生成边界文件: {boundary_file}")

    # === 合并重叠边界（动态 chromoL）===
    chromoL = lastP + 100000  # 动态基因组长度
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

    print(f"    DMR 分析完成: {base_name} → {final_file}")
    return final_file


def process_chr_in_one_file(df):
    """将输入的单个文件的前缀改为chr并将该文件所有的染色体信息返回"""

    # 该函数使得以chr开头的染色体号仍以chr开头，若不以chr开头的会加上前缀chr作为新的染色体号
    def to_chr(chr_str):
        chr_str = str(chr_str).strip()
        if chr_str.lower().startswith('chr'):
            return f"chr{chr_str[3:]}"
        else:
            return f"chr{chr_str}"

    # 将所有染色体号改为chr开头
    df['染色体号'] = df['染色体号'].apply(to_chr)
    # 将所有染色体号的不重复集合返回
    return df['染色体号'].unique()  # 返回numpy.array

def natural_sort_key(chr_name):
    """
    自定义染色体排序函数
    数字染色体按数值排序，字母染色体(X,Y,M等)排在最后按字母排序
    """
    chr_name = str(chr_name).lower()
    # 去除 'chr' 前缀
    if chr_name.startswith('chr'):
        suffix = chr_name[3:]
    else:
        suffix = chr_name

    # 尝试转换为整数
    try:
        # 如果是数字，返回 (0, 数字值, '')
        return (0, int(suffix), '')
    except ValueError:
        # 如果是字母(如X,Y,M)，返回 (1, 0, 字母)
        return (1, 0, suffix)

def scan_all_files_for_chr_mapping(m, n, dir1, dir2):
    """扫描m+n个、即所有文件，利用上一个函数收集所有染色体信息，以生成统一映射"""

    all_chromosomes = set()  # 利用集合不重复的特性存储所有染色体号

    # 扫描第一个目录的m个文件，名称格式为 i-dir1.txt,如3-msv.txt
    for i in range(1, m + 1):
        filepath = os.path.join(dir1, f"{i}-{os.path.basename(dir1)}.txt")
        if not os.path.exists(filepath):
            print(f"警告: 文件不存在 {filepath}")
            continue

        try:
            df = pd.read_csv(filepath,
                             sep='\s+',
                             header=None,
                             names=['染色体号', '位点号', '甲基化读段数', '非甲基化读段数', '甲基化类型'],
                             dtype=str)
            chromosomes = process_chr_in_one_file(df)  # 获取到该文件中的所有染色体号的不重复集合
            all_chromosomes.update(chromosomes)  # 将该文件的所有不重复的染色体号放到all_chromosomes这个集合里
            # 注：update接受任何可迭代对象并将其中元素添加到集合中
            print(f"文件 {filepath} 包含染色体: {sorted(chromosomes, key=natural_sort_key)}")
        except Exception as e:
            print(f"读取文件 {filepath} 时出错: {e}")

    # 扫描第二个目录的文件
    for j in range(1, n + 1):
        filepath = os.path.join(dir2, f"{j}-{os.path.basename(dir2)}.txt")
        if not os.path.exists(filepath):
            print(f"警告: 文件不存在 {filepath}")
            continue

        try:
            df = pd.read_csv(filepath,
                             sep='\s+',
                             header=None,
                             names=['染色体号', '位点号', '甲基化读段数', '非甲基化读段数', '甲基化类型'],
                             dtype=str)
            chromosomes = process_chr_in_one_file(df)  # 同上
            all_chromosomes.update(chromosomes)  # 同上
            print(f"文件 {filepath} 包含染色体: {sorted(chromosomes, key=natural_sort_key)}")
        except Exception as e:
            print(f"读取文件 {filepath} 时出错: {e}")

    # 创建统一的染色体映射，
    # 注意：此时all_chromosomes这个集合中包含两个基因型文件夹中所有文件的不重复染色体号
    unique_chrs = sorted(all_chromosomes, key=natural_sort_key)
    # 建立一个两个基因型目录中所有 不重复的 且排序了的 染色体号->数值的Series，或者说数组
    chr_series = pd.Series(range(len(unique_chrs)), index=unique_chrs)

    print(f"统一染色体映射: {chr_series}")
    return chr_series

def single_newtoboth(filepath1, output_dir, num1, chr_series):
    '''此处参数：filepath1即如1-wt.txt，是处理的新格式文件，
    output_dir是输出结果both文件到的目录，一般是filepath1所在的目录，
    num1是当前处理第几个新格式文件，即当前正处理num1-基因型.txt，
    chr_series是所有 不重复的 且排序了的 染色体号->数值的Series'''

    df = pd.read_csv(filepath1,
                     sep='\s+',
                     header=None,
                     names=['染色体号', '位点号', '甲基化读段数', '非甲基化读段数', '甲基化类型'],
                     dtype=str)

    def to_chr(chr_str):
        chr_str = str(chr_str).strip()
        if chr_str.lower().startswith('chr'):
            return f"chr{chr_str[3:]}"
        else:
            return f"chr{chr_str}"

    # 之前构建chr_series时没有修改源文件，所以此时读进来仍可能有不符合chr开头这个格式的染色体号，所以得处理
    df['染色体号'] = df['染色体号'].apply(to_chr)
    df['染色体号'] = df['染色体号'].map(chr_series)  # 将每个染色体号修改为其对应的数值，如chr1->0,chr2->1
    # 以下三行，将位点号、两个读段数转换为数值
    df['位点号'] = pd.to_numeric(df['位点号'])
    df['甲基化读段数'] = pd.to_numeric(df['甲基化读段数'])
    df['非甲基化读段数'] = pd.to_numeric(df['非甲基化读段数'])
    # 将df中的染色体号换为对应的chr_series中的索引，位点号和两个读段数换成数值类型，此外把甲基化类型CG转换为CpG
    df['甲基化类型'] = df['甲基化类型'].str.replace('CG', 'CpG')
    # 对三类数据进行分组，后续分别操作
    data_groups = df.groupby('甲基化类型')
    # 根据分组的键、对应的子数据表进行遍历
    for methy_type, data_ind in data_groups:
        if data_ind.empty:
            continue
        # 获取当前甲基化类型实际存在的染色体号并排序，
        # 是当前甲基化类型有的染色体号，注意：！！！是转换后从零开始的下标！！！
        actual_chrs = sorted(data_ind['染色体号'].dropna().unique())
        # 同时获取总染色体数，无论该甲基化类型是否具有对应的染色体信息
        chr_count = len(chr_series)
        chr_data_dict = {}
        # mlen用于记录当前甲基化类型中所有 染色体数据条数 中的最大值
        mlen = 0

        for chr_num in actual_chrs:  # 遍历该甲基化类型的每个染色体号（从零开始的数值）
            # 将当前染色体号的所有数据按照位点号进行升序排序，存储在chr_data中
            chr_data = data_ind[data_ind['染色体号'] == chr_num].sort_values('位点号').reset_index(drop=True)
            # 将三个关心的数据 位点号-甲基化读段数-非甲基化读段数 以染色体号为键存储在chr_data_dict字典中
            # 每个键对应的值是该三个属性组成的numpy数组
            chr_data_dict[chr_num] = chr_data[['位点号', '甲基化读段数', '非甲基化读段数']].values
            # 更新可能存在的更大的染色体数据条数到mlen中
            mlen = max(mlen, len(chr_data_dict[chr_num]))

        # 创建输出矩阵，行数为最大的染色体数据条数，列数为两个目录中所有存在的不重复染色体数*3，先用0填充
        output_matrix = np.zeros((mlen, chr_count * 3), dtype=np.int32)

        # 遍历输出矩阵中的所有行
        for i in range(mlen):
            # 遍历当前甲基化类型有的染色体号（从零开始编号的数值）
            for chr_num in actual_chrs:  # 使用连续索引
                col_start = chr_num * 3  # 基于连续索引计算列位置
                if i < len(chr_data_dict[chr_num]):
                    output_matrix[i, col_start:col_start + 3] = chr_data_dict[chr_num][i]

        # 输出（用pandas导出）
        output_df = pd.DataFrame(output_matrix)  # 将输出矩阵转换为DataFrame方便输出
        output_file = f"{num1}-bothMeUnme_diffChromo_NOREPEATED_methy_sites_{methy_type}.txt"
        output_path = os.path.join(output_dir, output_file)
        output_df.to_csv(output_path, sep='\t', header=False, index=False)

def get_chr_name(chr_num, chr_series):
    """
    根据染色体序号获取染色体名称
    参数：
        chr_num: 第几个染色体（从1开始）
        chr_series: 染色体映射 Series
    返回：
        染色体名称字符串（如 'chr1', 'chrX'）
    """
    chr_num = int(chr_num)
    if chr_series is not None:
        try:
            # chr_num 是从1开始的，所以要减1
            if 0 <= chr_num - 1 < len(chr_series):
                return chr_series.index[chr_num - 1]
        except Exception as e:
            print(f"警告：获取染色体名称失败: {e}")

    # 如果出错或没有映射，返回数字编号
    return f"chr{chr_num}"

def newtoboth(m, n, dir1, dir2):
    # 没必要检查目录是否存在了，因为main函数主逻辑检查过了
    # 获取两个基因型目录中所有 不重复的 且排序了的 染色体号->数值的Series
    chr_series = scan_all_files_for_chr_mapping(m, n, dir1, dir2)
    print(f"映射关系: {chr_series}")
    # 循环m+n次，分别处理两个基因型的文件
    for i in range(1, m + 1):
        filepath = os.path.join(dir1, f"{i}-{os.path.basename(dir1)}.txt")
        if not os.path.exists(filepath):
            print(f"警告: 文件不存在 {filepath}")
            continue
        print(f"处理文件{filepath}")
        # 将dir1文件夹中的i-dir1.txt转换为both格式
        single_newtoboth(filepath, dir1, i, chr_series)
    for j in range(1, n + 1):
        filepath = os.path.join(dir2, f"{j}-{os.path.basename(dir2)}.txt")
        if not os.path.exists(filepath):
            print(f"警告: 文件不存在 {filepath}")
            continue
        print(f"处理文件{filepath}")
        # 将dir2文件夹中的j-dir2.txt转换为both格式
        single_newtoboth(filepath, dir2, j, chr_series)
    return chr_series

def sanitize_filename(name):
    """清理文件名中的不允许在文件名中出现的特殊字符"""
    return re.sub(r'[\\/*?:"<>|]', "", name)  # 将name中的特殊字符替换为空，以使其无害化


def get_column_count(file_path):
    """获取文件的列数并返回列数除以3的结果"""
    try:
        with open(file_path, 'r') as file:
            first_line = file.readline().strip()
            column_count = len(first_line.split())
            return column_count // 3
    except Exception as e:
        print(f"读取文件时发生错误: {e}")
        return None


def parse_filename(filename):
    """解析单个文件名，使用正则表达式提取文件编号和甲基化类型"""
    # 用捕获组（这里有括号的部分）捕获对应的文件序号和甲基化类型
    pattern = r'^(\d+)-bothMeUnme_diffChromo_NOREPEATED_methy_sites_(.+)\.txt$'
    match = re.match(pattern, filename)  # 将特定的文件名和此正则表达式进行匹配
    if match:  # 匹配到后，用match.group(i)获取对应第i个捕获组的内容
        file_id = int(match.group(1))
        methylation_type = match.group(2)
        return file_id, methylation_type
    return None


def scan_sample_files_by_replicates(sample_dir, max_replicates):
    """在sample_dir目录中，搜寻出序号 1~max_replicates-both...甲基化类型.txt 文件"""
    files_by_replicates = {}  # 创建字典用于根据前端的序号存储对应的文件名
    methylation_types = ['CpG', 'CHH', 'CHG']
    # 若不存在对应目录，则返回空字典，这步合并后可省略，因为main函数中一开始就判断了
    # if not os.path.exists(sample_dir):
    #     return files_by_replicates

    # 遍历该文件夹中的所有内容（包括文件和子目录）
    for filename in os.listdir(sample_dir):
        if filename.endswith('.txt'):  # 若找到txt文件（只有可能是both文件或新格式文件）
            parsed = parse_filename(filename)  # 解析该文件获取可能存在的 序号,甲基化类型
            if parsed:  # 如果对于该文件解析到了 序号，甲基化类型
                file_id, methylation_type = parsed  # （元组，逗号是元组的标志，括号只是为了为避免歧义的产物）
                if file_id <= max_replicates and methylation_type in methylation_types:  # 序号和甲基化类型都合法
                    if file_id not in files_by_replicates:
                        files_by_replicates[file_id] = {}  # 使得files_by_replicates这个字典成为 file_id->{子字典}
                    files_by_replicates[file_id][methylation_type] = filename  # 使得子字典中格式为 甲基化类型->filename
    return files_by_replicates  # 就是返回一个链路 序号->甲基化类型->文件名 ，可以通过序号和甲基化类型获取到对应文件名


def process_methylation_type_with_collection(file1_path, file2_path, methylation_type, output_dir,
                                             dir1_name, dir2_name, chr_num):
    """处理甲基化类型和染色体的Fisher检验 处理1次检验的一个染色体
    参数解释：file_path-当前处理的both文件的相对路径
    methylation_type-当前处理的甲基化类型
    output_dir-输出文件夹./output_x_y/
    dir_name-输入的两个基因型目录路径的最后一个部分
    chr_num-处理的第几个染色体，而不是染色体号"""

    print(f"      处理甲基化类型 {methylation_type}，染色体 {chr_num}...")

    # 计算当前染色体的数据列
    mOrder = 3 * (chr_num - 1)  # 这里的chr_num是第几个染色体，而不是染色体号！！！
    usecols = [mOrder, mOrder + 1, mOrder + 2]

    try:
        # 分块读取第一个文件的所有数据到字典
        data1_dict = {}
        for chunk in pd.read_csv(file1_path, sep='\s+', header=None, usecols=usecols, chunksize=100000):
            chunk.columns = ['pos', 'methy', 'unmethy']
            # 用zip每次获取对应列的一个元素，共三个元素放到(pos,methy,unmethy)元组中
            for pos, methy, unmethy in zip(chunk['pos'].values, chunk['methy'].values, chunk['unmethy'].values):
                # 按照如下格式将数据放入对应字典中
                data1_dict[pos] = (methy, unmethy)
        # 至此第一个文件的所有位点数据都录入了data1_dict字典中

        # 分块读取第二个文件并查找共同位点，这是因为如果同时将两个文件的所有数据读到内存里，可能会占用太大的内存
        # 而第一个文件必须全部载入内存，因为我们需要对它进行快速的随机查找以确认某个位点是否在两个文件中都有，都有的话就得检验
        # 这里reader由于有chunksize这个参数，所以read_csv的返回值是一个迭代器，遍历它每次可以每次最多返回100000行数据
        reader2 = pd.read_csv(file2_path, sep='\s+', header=None, usecols=usecols, chunksize=100000)

    except Exception as e:
        print(f"        加载数据失败: {e}")
        return False

    # 创建输出文件夹:output_x_y/甲基化类型/
    methylation_dir = os.path.join(output_dir, methylation_type)
    os.makedirs(methylation_dir, exist_ok=True)

    # 创建输出文件路径
    stats_filename = os.path.join(methylation_dir,
                                  f"FET_results_{methylation_type}_{dir1_name}_{dir2_name}_Chr{chr_num}.txt")
    all_output = os.path.join(methylation_dir, f"all_simple_Chr{chr_num}.txt")
    sig_output = os.path.join(methylation_dir, f"sig_simple_Chr{chr_num}.txt")
    combine_output = os.path.join(methylation_dir, f"combineResult_Chr{chr_num}.txt")

    # 设置参数
    M0, M2, Th_pValue = 2, 4, 0.05

    try:
        # 创建四个列表用于存储各类数据以便后续输出
        all_results, sig_results, fet_results, combine_results = [], [], [], []

        # 处理第二个文件的每个数据块
        for chunk2 in reader2:
            chunk2.columns = ['pos', 'methy', 'unmethy']
            # 对于第二个文件的每个数据块，遍历每一行三个值，存储到pos,m2,u2中
            for pos, m2, u2 in zip(chunk2['pos'].values, chunk2['methy'].values, chunk2['unmethy'].values):
                # 若该位点在第一个文件也存在，说明要进行该位点的fisher检验
                if pos in data1_dict:
                    m1, u1 = data1_dict[pos]  # 获取第一个文件的两个数

                    if m1 >= M0 or m2 >= M0:  # 两个methy读段数都要>=2才可能进行下一步
                        if (m1 + u1 >= M2) and (m2 + u2 >= M2):  # 每个文件的两个数之和都要>=4才进行检验
                            cont_table = np.array([[m1, u1], [m2, u2]])  # 建立2*2列联表
                            # 计算两个甲基化率
                            ratio1 = m1 / (m1 + u1) if (m1 + u1) > 0 else 0
                            ratio2 = m2 / (m2 + u2) if (m2 + u2) > 0 else 0
                            change = "1" if ratio1 >= ratio2 else "0"  # 根据两个文件的甲基化比率判断
                            # 突变型的甲基化率是否升高了
                            # 调用库函数进行fisher检验，获取到该次检验的pvalue
                            _, pvalue = fisher_exact(cont_table, alternative='two-sided')
                            # 为pvalue保留7位有效数字
                            pvalue = float(f"{pvalue:.7g}")
                            # 为四个文件分别录入需要的数据，其中只有显著的才录入sig_results中
                            all_results.append([pos, pvalue, change])
                            fet_results.append([pos, m1, u1, m2, u2, pvalue])
                            combine_results.append([pos, m1, u1, m2, u2])

                            if pvalue < Th_pValue:
                                sig_results.append([pos, pvalue, change])
        # 至此，这组文件（特定x_y特定甲基化类型特定染色体）所有需要进行的fisher检验已完成
        # 保存结果到磁盘中
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

        print(f"        染色体 {chr_num} 处理完成！共处理 {len(all_results)} 个位点，其中 {len(sig_results)} 个显著")
        return True

    except Exception as e:
        print(f"        处理甲基化类型 {methylation_type}，染色体 {chr_num} 时发生错误: {e}")
        return False

def process_methylation_type_with_collection_pvfilter(file1_path, file2_path, methylation_type, output_dir,
                                             dir1_name, dir2_name, chr_num):
    """处理甲基化类型和染色体的Fisher检验 处理1次检验的一个染色体
    参数解释：file_path-当前处理的both文件的相对路径
    methylation_type-当前处理的甲基化类型
    output_dir-输出文件夹./output_x_y/
    dir_name-输入的两个基因型目录路径的最后一个部分
    chr_num-处理的第几个染色体，而不是染色体号

    返回此次的pvalue>0.05的df
        all_results.append([pos, pvalue, change])
    """

    print(f"      处理甲基化类型 {methylation_type}，染色体 {chr_num}...")

    # 计算当前染色体的数据列
    mOrder = 3 * (chr_num - 1)  # 这里的chr_num是第几个染色体，而不是染色体号！！！
    usecols = [mOrder, mOrder + 1, mOrder + 2]

    try:
        # 分块读取第一个文件的所有数据到字典
        data1_dict = {}
        for chunk in pd.read_csv(file1_path, sep='\s+', header=None, usecols=usecols, chunksize=100000):
            chunk.columns = ['pos', 'methy', 'unmethy']
            # 用zip每次获取对应列的一个元素，共三个元素放到(pos,methy,unmethy)元组中
            for pos, methy, unmethy in zip(chunk['pos'].values, chunk['methy'].values, chunk['unmethy'].values):
                # 按照如下格式将数据放入对应字典中
                data1_dict[pos] = (methy, unmethy)
        # 至此第一个文件的所有位点数据都录入了data1_dict字典中

        # 分块读取第二个文件并查找共同位点，这是因为如果同时将两个文件的所有数据读到内存里，可能会占用太大的内存
        # 而第一个文件必须全部载入内存，因为我们需要对它进行快速的随机查找以确认某个位点是否在两个文件中都有，都有的话就得检验
        # 这里reader由于有chunksize这个参数，所以read_csv的返回值是一个迭代器，遍历它每次可以每次最多返回100000行数据
        reader2 = pd.read_csv(file2_path, sep='\s+', header=None, usecols=usecols, chunksize=100000)

    except Exception as e:
        print(f"        加载数据失败: {e}")
        return False

    # 创建输出文件夹:output_x_y/甲基化类型/
    methylation_dir = os.path.join(output_dir, methylation_type)
    os.makedirs(methylation_dir, exist_ok=True)

    # 创建输出文件路径
    stats_filename = os.path.join(methylation_dir,
                                  f"FET_results_{methylation_type}_{dir1_name}_{dir2_name}_Chr{chr_num}.txt")
    all_output = os.path.join(methylation_dir, f"all_simple_Chr{chr_num}.txt")
    sig_output = os.path.join(methylation_dir, f"sig_simple_Chr{chr_num}.txt")
    combine_output = os.path.join(methylation_dir, f"combineResult_Chr{chr_num}.txt")

    # 设置参数
    M0, M2, Th_pValue = 2, 4, 0.05

    try:
        # 创建四个列表用于存储各类数据以便后续输出
        all_results, sig_results, fet_results, combine_results = [], [], [], []

        # 处理第二个文件的每个数据块
        for chunk2 in reader2:
            chunk2.columns = ['pos', 'methy', 'unmethy']
            # 对于第二个文件的每个数据块，遍历每一行三个值，存储到pos,m2,u2中
            for pos, m2, u2 in zip(chunk2['pos'].values, chunk2['methy'].values, chunk2['unmethy'].values):
                # 若该位点在第一个文件也存在，说明要进行该位点的fisher检验
                if pos in data1_dict:
                    m1, u1 = data1_dict[pos]  # 获取第一个文件的两个数

                    if m1 >= M0 or m2 >= M0:  # 两个methy读段数都要>=2才可能进行下一步
                        if (m1 + u1 >= M2) and (m2 + u2 >= M2):  # 每个文件的两个数之和都要>=4才进行检验
                            cont_table = np.array([[m1, u1], [m2, u2]])  # 建立2*2列联表
                            # 计算两个甲基化率
                            ratio1 = m1 / (m1 + u1) if (m1 + u1) > 0 else 0
                            ratio2 = m2 / (m2 + u2) if (m2 + u2) > 0 else 0
                            change = "1" if ratio1 > ratio2 else "0"  # 根据两个文件的甲基化比率判断
                            # 突变型的甲基化率是否升高了
                            # 调用库函数进行fisher检验，获取到该次检验的pvalue
                            _, pvalue = fisher_exact(cont_table, alternative='two-sided')
                            # 为pvalue保留7位有效数字
                            pvalue = float(f"{pvalue:.7g}")
                            # 为四个文件分别录入需要的数据，其中只有显著的才录入sig_results中
                            all_results.append([pos, pvalue, change])
                            fet_results.append([pos, m1, u1, m2, u2, pvalue])
                            combine_results.append([pos, m1, u1, m2, u2])

                            if pvalue < Th_pValue:
                                sig_results.append([pos, pvalue, change])
        # 至此，这组文件（特定x_y特定甲基化类型特定染色体）所有需要进行的fisher检验已完成
        # 保存结果到磁盘中
        if all_results:
            # 转换为DataFrame方便筛选
            all_df = pd.DataFrame(all_results, columns=["Position", "Pvalue", "Methylation_Change"])
            sig_df = pd.DataFrame(sig_results, columns=["Position", "Pvalue", "Methylation_Change"])
            fet_df = pd.DataFrame(fet_results,
                                  columns=["Position", "Methy1", "Unmethy1", "Methy2", "Unmethy2", "Pvalue"])
            combine_df = pd.DataFrame(combine_results, columns=["Position", "Methy1", "Unmethy1", "Methy2", "Unmethy2"])

            all_df_ndmp = all_df[all_df['Pvalue'] > 0.05]  # 保留pvalue>0.05的信息，后续

            # 只保留 p≤0.05 的位点
            all_df = all_df[all_df['Pvalue'] <= 0.05]  # 筛选pvalue<=0.05的进行FDR检验
            fet_df = fet_df[fet_df['Pvalue'] <= 0.05]  # 筛选pvalue<=0.05的进行FDR检验
            #fet_df_ndmp = all_df[all_df['Pvalue'] > 0.05]
            # sig_df 本来就是 p<0.05，不需要再筛选

            # 也筛选 combine
            positions_to_keep = set(all_df['Position'].values)
            combine_df = combine_df[combine_df['Position'].isin(positions_to_keep)]  # 筛选pvalue<=0.05的进行FDR检验
            #combine_df_ndmp = combine_df[~combine_df['Position'].isin(positions_to_keep)]
            # 保存筛选后的结果
            all_df.to_csv(all_output, sep='\t', index=False)
            sig_df.to_csv(sig_output, sep='\t', index=False)
            fet_df.to_csv(stats_filename, sep='\t', index=False)
            combine_df.to_csv(combine_output, sep='\t', index=False)

        print(
            f"        染色体 {chr_num} 处理完成！原始检验 {len(all_results)} 个位点，筛选后(p≤0.05) {len(all_df)} 个位点，p>0.05的共{len(all_df_ndmp)}个位点")
        return all_df_ndmp
    #                   all_results.append([pos, pvalue, change])
    #                   fet_results.append([pos, m1, u1, m2, u2, pvalue])
    #               combine_results.append([pos, m1, u1, m2, u2])

    except Exception as e:
        print(f"        处理甲基化类型 {methylation_type}，染色体 {chr_num} 时发生错误: {e}")
        return False

def merge_fet_results_and_fdr(output_dir, replicate_x, replicate_y, mtype3, all_dfs_ndmp_dict,n_chromosomes):
    """合并output文件夹中的所有FET结果并进行FDR校正
    其中FET文件的格式如：pos, m1, u1, m2, u2, pvalue
    这里输入的output_dir是output_x_y/甲基化类型
    success_dfs_dict[methylation_type][chr_num]访问的是
                    对应次pvalue>0.05的df
                      all_results([pos, pvalue, change])
    n_chromosomes是总染色体数
                      """
    print(f"\n    合并 {output_dir} 中的FET结果并进行FDR校正...")

    if not os.path.exists(output_dir):
        print(f"    错误：目录 {output_dir} 不存在！")
        return False

    # 搜索所有FET结果文件
    # 这里**表示任意深度的子目录，所以这里glob在递归搜索（recursive=True）时就会在output_dir目录下、以及其
    # 所有子目录下搜寻符合条件的文件，并将符合条件的文件路径（从output_dir开始的相对路径）以列表的形式返回
    file_pattern = os.path.join(output_dir, "**", "FET_results_*_Chr*.txt") # 这里的染色体号实际上是both文件中的第几个三列
    fet_files = glob.glob(file_pattern, recursive=True)

    if not fet_files:
        print(f"    警告：在 {output_dir} 中未找到FET结果文件")
        return False

    print(f"    找到 {len(fet_files)} 个FET结果文件")

    # 创建列表用于收集所有p值和相关信息
    all_data = []

    for file_path in sorted(fet_files):
        # 提取甲基化类型和染色体信息
        # 这里第一个捕获组用于捕获甲基化类型，.*为零次或多次任意字符，用于匹配replicatex_replicatey，第二个捕获组用于捕获染色体号
        #                                                   # 这里的染色体号实际上是both文件中的第几个三列
        methy_match = re.search(r'/FET_results_([^_]+)_.*_Chr(\d+)\.txt$', file_path.replace('\\', '/'))
        if not methy_match:
            continue

        # 获取到甲基化类型和染色体序号（这里的染色体号是both文件中的第几个三列）
        methylation_type = mtype3
        chr_num = int(methy_match.group(2))

        try:
            df = pd.read_csv(file_path, sep='\t', header=0)
            if 'Position' in df.columns and 'Pvalue' in df.columns:
                # 调整列顺序：Chromosome, Methylation_Type, Position, Pvalue
                df_subset = df[['Position', 'Pvalue']].copy()
                df_subset['Chromosome'] = chr_num
                df_subset['Methylation_Type'] = methylation_type

                # 重新排列列顺序
                df_subset = df_subset[['Chromosome', 'Methylation_Type', 'Position', 'Pvalue']]
                all_data.append(df_subset) # 将当前FET的所有信息加入all_data中
        except Exception as e:
            print(f"    警告：读取 {file_path} 失败: {e}")
            continue

    if not all_data:
        print(f"    错误：没有成功读取任何数据")
        return False

    # 合并所有数据（由于上面是直接将每个文件对应的dataframe追加到all_data末尾，所以列表中每个元素都是一个dataframe，因此需要合并）
    combined_df = pd.concat(all_data, ignore_index=True)
    # 将数据根据染色体序号排序（这里的染色体号实际上是both文件中的第几个三列）
    combined_df = combined_df.sort_values(['Methylation_Type', 'Chromosome', 'Position'])
    # 该dataframe格式如下：'Chromosome', 'Methylation_Type', 'Position', 'Pvalue'
    print(f"    合并后总共 {len(combined_df)} 个位点")

    # 计算FDR校正的q值并将其向右增加到combined_df中
    pvalues = combined_df['Pvalue'].values
    qvalues = calculate_qvalues(pvalues, 1.0)
    combined_df['Qvalue'] = qvalues

    # 最终列顺序：Chromosome, Methylation_Type, Position, Pvalue, Qvalue
    combined_df = combined_df[['Chromosome', 'Methylation_Type', 'Position', 'Pvalue', 'Qvalue']]

    print(f"开始处理 {mtype3} 的 ndmp 数据，")

    # 首先检查 all_dfs_ndmp_dict 中该甲基化类型是否有数据，因为有可能该生物种类的该甲基化类型不需要调用filter版本
    if mtype3 in all_dfs_ndmp_dict and all_dfs_ndmp_dict[mtype3]:
        print(f"该甲基化类型中，ndmp的dfs数量为{len(all_dfs_ndmp_dict[mtype3])}个")
        dfs_ndmp = []
        for chr_num11 in range(1, n_chromosomes + 1):
            # 只处理该甲基化类型中存在的染色体
            if chr_num11 in all_dfs_ndmp_dict[mtype3]:
                df_ndmp = all_dfs_ndmp_dict[mtype3][chr_num11]

                # 确保是有效的 DataFrame
                if isinstance(df_ndmp, pd.DataFrame) and not df_ndmp.empty:
                    df_ndmp = df_ndmp.copy()

                    df_ndmp['Chromosome'] = chr_num11
                    df_ndmp['Methylation_Type'] = mtype3

                    if 'change' in df_ndmp.columns:
                        df_ndmp.drop(columns=['change'], inplace=True)

                    df_ndmp['Qvalue'] = 1
                    df_ndmp = df_ndmp[['Chromosome', 'Methylation_Type', 'Position', 'Pvalue', 'Qvalue']]

                    dfs_ndmp.append(df_ndmp)
                    print(f"  为 {mtype3}--{chr_num11} 添加所需列，共 {len(df_ndmp)} 个位点")

        if dfs_ndmp:
            total_ndmp = pd.concat(dfs_ndmp, ignore_index=True)
            combined_df = pd.concat([combined_df, total_ndmp], ignore_index=True)
            print(f"  共合并 {len(total_ndmp)} 个 ndmp 位点")
        else:
            print(f"  {mtype3} 类型未找到任何 ndmp 数据")
    else:
        print(f"  {mtype3} 类型在 ndmp 字典中不存在")
    # 至此对于比较多的甲基化类型已加上前面舍弃的pvalue>0.05的位点的信息，
    #   在前面要将其存储在all_dfs_ndmp_dict[甲基化][染色体号]里，此时才能实现该目的

    # 统计显著性结果，计算这 m*n*3次检验中的1次中，pvalues和qvalues显著的分别有多少个
    n_pval_sig = np.sum(pvalues < 0.05)
    n_qval_sig = np.sum(qvalues < 0.05)
    # 计算显著的比例
    print(f"    P值显著位点: {n_pval_sig} ({n_pval_sig / len(pvalues) * 100:.1f}%)")
    print(f"    Q值显著位点: {n_qval_sig} ({n_qval_sig / len(qvalues) * 100:.1f}%)")

    # 保存合并的p值列表（用于外部FDR工具）
    pvalue_file = os.path.join(output_dir, f"united_pvalues_wt_replicate{replicate_y}_vs_mut_replicate{replicate_x}.csv")
    with open(pvalue_file, 'w') as f:
        for pvalue in pvalues:
            f.write(f"{pvalue}\n")

    # 保存完整的FDR校正结果（output_dir是output_x_y/甲基化类型）
    #  包含完整的pvalues和qvalues（格式为：Chromosome, Methylation_Type, Position, Pvalue, Qvalue）
    fdr_file = os.path.join(output_dir, f"FDR_corrected_results_wt_replicate{replicate_y}_vs_mut_replicate{replicate_x}.txt")
    combined_df.to_csv(fdr_file, sep='\t', index=False)

    # 保存显著位点（q值<0.05）
    sig_df = combined_df[combined_df['Qvalue'] < 0.05]
    if not sig_df.empty: # 若有显著的位点，将显著的那部分数据输出（output_dir是output_x_y/甲基化类型）
        sig_file = os.path.join(output_dir, f"FDR_significant_results_wt_replicate{replicate_y}_vs_mut_replicate{replicate_x}.txt")
        sig_df.to_csv(sig_file, sep='\t', index=False)
        print(f"    显著位点结果保存至: {sig_file}")

    print(f"    P值列表保存至: {pvalue_file}")
    print(f"    FDR结果保存至: {fdr_file}")
    return True


# 在你的代码中，将calculate_qvalues函数替换为：
def calculate_qvalues(pvalues, pi=1.0):
    """使用Storey方法计算Q值"""

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

    # Storey方法的关键：包含π因子
    q_values_sorted = np.zeros_like(sorted_pvalues)

    # 如果需要自动估计π
    if pi is None:
        # 简化的π估计
        lrange = np.linspace(0.05, 0.95, max(int(m / 100.0), 10))
        pil = np.mean(sorted_pvalues[:, np.newaxis] > lrange, axis=0)
        pilr = pil / (1.0 - lrange)
        pi = 1.0
        if pilr[-1] < 1.0:
            pi = pilr[-1]

    # Storey方法计算
    q_values_sorted = pi * m * sorted_pvalues / np.arange(1, m + 1)
    q_values_sorted[-1] = min(q_values_sorted[-1], 1.0)

    # 单调性调整
    for i in range(m - 2, -1, -1):
        q_values_sorted[i] = min(q_values_sorted[i], q_values_sorted[i + 1])

    # 恢复原顺序
    q_values = np.zeros_like(pvalues_clean)
    q_values[sorted_indices] = q_values_sorted
    q_values[np.isnan(pvalues)] = np.nan

    return q_values

def perform_sliding_window_on_dmp_files(output_dir, replicate_x, replicate_y,):
    """对m*n*3次处理中的一次的结果(N)DMP文件进行滑动窗口分析"""

    print(f"\n    开始对DMP文件进行滑动窗口分析...")

    # 只处理DMP文件
    dmp_patterns = [
        f"DMP_wt_replicate{replicate_y}_mut_replicate{replicate_x}_Chr*.txt",
        f"N-DMP_wt_replicate{replicate_y}_mut_replicate{replicate_x}_Chr*.txt"
    ]

    for pattern in dmp_patterns:
        dmp_files = glob.glob(os.path.join(output_dir, pattern)) #找到符合当前正则表达式格式的所有文件的路径，并放入dmp_files列表中

        for dmp_file in dmp_files: # 遍历列表，每次获取一个DMP或N-DMP文件
            try:

                df = pd.read_csv(dmp_file, sep='\s+', header=None, skiprows=1,
                                 names=['position', 'pvalue', 'change']) # 位点号、pvalue、change

                # 检验数据非空
                if df.empty or len(df) == 0:
                    continue

                # 数据类型转换和清理
                df = df.dropna()
                df['position'] = df['position'].astype(int)
                df['change'] = df['change'].astype(int)

                print(f"      处理文件: {os.path.basename(dmp_file)} ({len(df)} 个位点)")

                # 设置输出前缀，其实就是将.txt去掉之后的 (N)DMP_replicate_wt{replicate_y}_mut_replicate{replicate_x}_Chr* 这一部分
                base_name = os.path.basename(dmp_file).replace('.txt', '')

                # 传入DataFrame给sliding_window_analysis
                # 前者格式为window_start,window_end,count_change_1,count_change_0,total_count
                # 后者没有count_change，而有standardized_count，即每个区间显著位点数与最大的那个显著位点数的比值
                sliding_results, std_results = sliding_window_analysis(
                    df, # df是m*n*3次中的一次的某一条染色体的dmp文件
                    window_size=1000000,
                    step_ratio=0.05,
                    save_files=True,
                    output_identifier=base_name,
                    outputdir1=output_dir
                )

                print(f"      完成 {base_name}: {len(sliding_results)} 个窗口")

            except Exception as e:
                print(f"      处理 {dmp_file} 失败: {e}")
                continue

    print(f"    滑动窗口分析完成")

#下面这个版本用于处理pv<0.05筛选后再进行FDR检验的情况，比如植物的CHH和CHG
def perform_sliding_window_on_dmp_files_after_filter(output_dir, replicate_x, replicate_y,all_dfs_ndmp_dict=None, methylation_type=None):
    """对m*n*3次处理中的一次的结果(N)DMP文件进行滑动窗口分析"""

    print(f"\n    开始对DMP文件进行滑动窗口分析...")

    # 只处理DMP文件
    dmp_patterns = [
        f"DMP_wt_replicate{replicate_y}_mut_replicate{replicate_x}_Chr*.txt",
        f"N-DMP_wt_replicate{replicate_y}_mut_replicate{replicate_x}_Chr*.txt"
    ]

    for pattern in dmp_patterns:
        dmp_files = glob.glob(os.path.join(output_dir, pattern)) #找到符合当前正则表达式格式的所有文件的路径，并放入dmp_files列表中
        is_ndmp = pattern.startswith("N-DMP") # 判断当前处理的文件是否是NDMP文件

        for dmp_file in dmp_files: # 遍历列表，每次获取一个DMP或N-DMP文件
            try:

                df = pd.read_csv(dmp_file, sep=' ', header=None, skiprows=1,
                                 names=['position', 'pvalue', 'change']) # 位点号、pvalue、change

                #如果当前文件是ndmp文件，由于不能改all_simple文件，因为那个要用来合并成一列用于FDR检验，所以就在此加上原来pv>0.05的
                                # 那些位点的信息，这样的话，后续生成滑动窗口文件的时候，统计ndmp信息的时候才不会缺漏
                if is_ndmp and all_dfs_ndmp_dict and methylation_type:
                    # 从文件名提取染色体号
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
                                # 合并数据
                                df = pd.concat([df, ndmp_df[['position', 'pvalue', 'change']]],
                                               ignore_index=True)
                                df = df.drop_duplicates(subset=['position'])
                                df = df.sort_values('position').reset_index(drop=True)
                                print(f"      为 {os.path.basename(dmp_file)} 合并了之前FDR检验时忽略的 {len(ndmp_df)} 个NDMP位点")
                # 检验数据非空
                if df.empty or len(df) == 0:
                    continue

                # 数据类型转换和清理
                df = df.dropna()
                df['position'] = df['position'].astype(int)
                df['change'] = df['change'].astype(int)

                print(f"      处理文件: {os.path.basename(dmp_file)} ({len(df)} 个位点)")

                # 设置输出前缀，其实就是将.txt去掉之后的 (N)DMP_wt_replicate{replicate_y}_mut_replicate{replicate_x}_Chr* 这一部分
                base_name = os.path.basename(dmp_file).replace('.txt', '')

                # 传入DataFrame给sliding_window_analysis
                # 前者格式为window_start,window_end,count_change_1,count_change_0,total_count
                # 后者没有count_change，而有standardized_count，即每个区间显著位点数与最大的那个显著位点数的比值
                sliding_results, std_results = sliding_window_analysis(
                    df, # df是m*n*3次中的一次的某一条染色体的dmp文件
                    window_size=1000000,
                    step_ratio=0.05,
                    save_files=True,
                    output_identifier=base_name,
                    outputdir1=output_dir
                )

                print(f"      完成 {base_name}: {len(sliding_results)} 个窗口")

            except Exception as e:
                print(f"      处理 {dmp_file} 失败: {e}")
                continue

    print(f"    滑动窗口分析完成")


def generate_dmp_files(dir1,dir2,output_dir, replicate_x, replicate_y, fdr_threshold=0.05, mtype1="CpG",
                      all_dfs_ndmp_dict=None,unfilter_mtypes=["CpG"],n_chromosomes = 5):
    """参数：output_dir是output_x_y/甲基化类型/,
        replicate_x,replicate_y是处理的组号，后面是qvalues的阈值和目前处理的甲基化类型"""

    print(f"\n    生成 {output_dir} 的DMP文件...")

    def safe_float_convert(value):
        """将各种格式的数据转换为浮点数"""
        try:
            # 如果已经是数字类型，均转换为浮点数返回
            if isinstance(value, (int, float)):
                return float(value)
            # 如果是字符串，去除空白后转换为浮点数
            if isinstance(value, str):
                value = value.strip()
                if value:
                    return float(value)
            return None
        except (ValueError, TypeError):
            return None

    # 1. 读取FDR校正结果文件(这里output_dir是output_x_y/甲基化类型/)
    #  该文件格式为 Chromosome, Methylation_Type, Position, Pvalue, Qvalue
    fdr_file = os.path.join(output_dir, f"FDR_corrected_results_wt_replicate{replicate_y}_vs_mut_replicate{replicate_x}.txt")
    if not os.path.exists(fdr_file):
        print(f"    错误：FDR结果文件不存在 {fdr_file}")
        return False

    try:
        fdr_df = pd.read_csv(fdr_file, sep='\t') # 该文件格式为 Chromosome, Methylation_Type, Position, Pvalue, Qvalue
        print(f"    读取FDR结果文件，共 {len(fdr_df)} 个位点")
    except Exception as e:
        print(f"    错误：读取FDR文件失败 {e}")
        return False

    # 2. 读取各染色体的all_simple文件
    methylation_change_data = {}
    total_read_lines = 0

    # (这里output_dir是output_x_y/甲基化类型/)
    mtype_dir = output_dir

    # 将该m*n*3次中的一次处理的所有all_simple文件的文件名获取到列表中
    all_simple_files = [f for f in os.listdir(mtype_dir) if f.startswith('all_simple_Chr') and f.endswith('.txt')]

    # all_simple文件格式为："Position", "Pvalue", "Methylation_Change"
    # 遍历所给一次处理的所有all_simple文件
    for file in all_simple_files:
        chr_match = re.search(r'Chr(\d+)\.txt$', file)
        if not chr_match:
            continue

        chr_num = int(chr_match.group(1))# 这里的chr_num也指的是both文件的第几个三列
        file_path = os.path.join(mtype_dir, file)

        try:
            with open(file_path, 'r') as f:  # 读取的all_simple文件格式为："Position", "Pvalue", "Methylation_Change"
                lines = f.readlines() # readlines用于返回文件的所有行组成的列表，每个元素是一行

            file_valid_lines = 0
            for line_num, line in enumerate(lines, 1): #遍历每一行，其中line_num从1开始枚举
                line = line.strip()
                if not line or line.startswith("Position"): # 跳过空行或首行
                    continue

                parts = line.split('\t')    #根据之前导出的时候设定的分隔符进行分割，将每次的三个数据项放到列表parts里
                if len(parts) >= 3:
                    # 由于读入的是字符串，所以得进行类型转换，将pvalue转换为浮点数，其余两个转换为整型
                    position = int(safe_float_convert(parts[0]))
                    pvalue = safe_float_convert(parts[1])
                    change = int(safe_float_convert(parts[2]))

                    # 检查所有值是否有效
                    if (position is not None and
                            pvalue is not None and
                            change is not None):

                        # 检验change是否在0,1之中
                        if change in [0, 1]:
                            # 存储该映射: (chr, mtype, position) -> change 到字典中
                            methylation_change_data[(chr_num, mtype1, position)] = change
                            file_valid_lines += 1
                            total_read_lines += 1

            print(f"    {mtype_dir}/Chr{chr_num}: 读取 {file_valid_lines} 个有效位点")

        except Exception as e:
            print(f"    警告：读取 {file_path} 失败: {e}")
            continue

    print(f"    总共读取甲基化变化方向数据: {total_read_lines} 个位点")

    # 3. 合并数据，利用前面记录在字典中的(chr, mtype, position) -> change,将change这个属性添加到fdr_df的副本中
    combined_data = []
    missing_change = 0
    match_debug = defaultdict(int) #创建一个带默认值的字典，当访问不存在的键的时候，该字典会自动调用int类内置的调用函数
                                                # 为该键创建一个键值对，并将值置为int()的返回值——0

    for _, row in fdr_df.iterrows(): # fdr_df文件格式为 Chromosome, Methylation_Type, Position, Pvalue, Qvalue
        chr_num = row['Chromosome']
        mtype = row['Methylation_Type']
        # 统一转换为浮点数进行匹配
        position = safe_float_convert(row['Position'])
        qvalue = safe_float_convert(row['Qvalue'])

        # 若位点或q值缺失，则跳过当前行
        if position is None or qvalue is None:
            missing_change += 1
            continue

        # 查找对应的甲基化变化方向
        change_key = (chr_num, mtype, position)
        if change_key in methylation_change_data:
            change = methylation_change_data[change_key]  # 该change_data中存储该映射: (chr, mtype, position) -> change
            combined_data.append({
                'chromosome': chr_num,
                'methylation_type': mtype,
                'position': int(position),  # 输出时转为整数
                'qvalue': qvalue,
                'change': change # 相当于将change这个属性添加到fdr_df的副本中，不过这个时候combined_data还是个列表而非dataframe
            })
            match_debug[chr_num] += 1 # 这里的chr_num也是both文件中第几个三列，从fisher检验读取Both文件就开始是这样了
        else:
            missing_change += 1

    print(f"    合并后有 {len(combined_data)} 个完整位点")
    print(f"    各染色体匹配情况: {dict(match_debug)}") # 其实就是每个染色体有多少条信息
    if missing_change > 0:
        print(f"    警告：{missing_change} 个位点缺少甲基化变化方向信息")

    # 4. 按染色体分组生成DMP文件
    chr_groups = defaultdict(list) # 这里的默认值字典在初始化时会调用list()函数，用[]空列表作为值
    for item in combined_data: # 每个item是一个字典，五个键值对的值为chr_num,mtype,position,qvalue,change
        chr_groups[item['chromosome']].append(item)
            # chr_groups这个字典将每个染色体的每一条数据（字典形式）作为一个元素记录到chr_groups->chr_num->list1这个列表中

    if not chr_groups:
        print(f"    错误：没有找到任何可生成DMP文件的数据")
        return False

    print(f"    将为以下染色体生成DMP文件: {sorted(chr_groups.keys())}")

    dir1_name = f"replicate{replicate_x}"
    dir2_name = f"replicate{replicate_y}"

    total_dmp = total_ndmp = total_hyper = total_hypo = 0

    # 为每个染色体生成文件
    for chr_num in sorted(chr_groups.keys()): # 这里的keys是染色体序号，序号-1就是一开始的全染色体映射里的索引，因为那个是从零开始算的
        chr_data = chr_groups[chr_num]  # chr_groups这个字典将每个染色体的每一条数据（字典形式）作为一个元素记录下来
                                    #所以这里的chr_data还是一个字典
        chr_data.sort(key=lambda x: x['position']) #lambda表达式匿名函数应用到列表的每个元素上，获取返回值作为排序依据
                                #  这里的匿名函数相当于对于chr_data中的每个元素（是一个字典）应用该函数
                                                #   def get_position(x):
                                            #     return x['position'] 返回值是位点号，所以就是根据位点号排序

        # 生成文件名（这里output是甲基化类型目录）
        dmp_file = os.path.join(output_dir, f"DMP_wt_{dir2_name}_mut_{dir1_name}_Chr{chr_num}.txt")
        ndmp_file = os.path.join(output_dir, f"N-DMP_wt_{dir2_name}_mut_{dir1_name}_Chr{chr_num}.txt")
        hyper_file = os.path.join(output_dir, f"hyper_DMP_wt_{dir2_name}_mut_{dir1_name}_Chr{chr_num}.txt")
        hypo_file = os.path.join(output_dir, f"hypo_DMP_wt_{dir2_name}_mut_{dir1_name}_Chr{chr_num}.txt")

        # 分类数据
        dmp_data = []
        ndmp_data = []
        hyper_data = []
        hypo_data = []

        # 遍历chr_data这个字典
        for item in chr_data: # chr_data是一个字典，五个键值对的值为chr_num,mtype,position,qvalue,change
            position = item['position']
            qvalue = item['qvalue']
            change = item['change']

            if qvalue < fdr_threshold: #判断是否显著
                dmp_data.append((position, qvalue, change))
                if change == 1:
                    hyper_data.append((position, qvalue, change))
                elif change == 0:
                    hypo_data.append((position, qvalue, change))
            else:
                ndmp_data.append((position, qvalue, change))

        # 写入文件
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
        # 统计
        total_dmp += len(dmp_data)
        total_ndmp += len(ndmp_data)
        total_hyper += len(hyper_data)
        total_hypo += len(hypo_data)

    print(f"    DMP文件生成完成！")
    print(f"    总计: DMP={total_dmp}, N-DMP={total_ndmp}, Hyper={total_hyper}, Hypo={total_hypo}")

    bothfile1 = os.path.join(dir1,f"{replicate_x}-bothMeUnme_diffChromo_NOREPEATED_methy_sites_{mtype1}.txt")

    bothfile2 = os.path.join(dir2,f"{replicate_y}-bothMeUnme_diffChromo_NOREPEATED_methy_sites_{mtype1}.txt")

    summarize_dmr_methylation(output_dir, replicate_x, replicate_y, bothfile1, bothfile2, mtype1,custom_dmr_dir=None)

    # 根据甲基化类型决定使用哪个滑动窗口函数
    if mtype1 not in unfilter_mtypes:
        print(f"使用新版本的滑动窗口分析")
        # 这里本来也要区分甲基化类型
        perform_sliding_window_on_dmp_files_after_filter(
            output_dir, replicate_x, replicate_y,
            all_dfs_ndmp_dict=all_dfs_ndmp_dict,
            methylation_type=mtype1
        )
        # perform_sliding_window_on_dmp_files(
        #     output_dir, replicate_x, replicate_y
        # )
    else:
        print(f"    使用标准版本的滑动窗口分析")
        perform_sliding_window_on_dmp_files(output_dir, replicate_x, replicate_y)
    return True


# 将此函数集成到 process_replicate_pair 中
def process_replicate_pair(replicate_x, replicate_y, files1, files2, dir1, dir2, dir1_name, dir2_name,unfilter_mtypes,work_dir="."):
    """处理一对组的所有甲基化类型，是处理1*3次检验
    replicate_x,replicate_y为both文件的序号
    files为一个字典链路：序号->甲基化类型->对应文件名
    dir1,dir2为两个基因型数据的目录名
    dir1_name，dir2_name是所给两个dir1,dir2路径的最后一个部分"""

    print(f"\n  处理组对 (wt{replicate_y}, mut{replicate_x})...")

    # 创建该组输出目录
    output_dir = os.path.join(work_dir, f"output_wt{replicate_y}_mut{replicate_x}")
    os.makedirs(output_dir, exist_ok=True)

    methylation_types = ['CpG', 'CHH', 'CHG']

    all_dfs_ndmp_dict = {}

    # 循环处理每种甲基化类型
    for methylation_type in methylation_types:
        if methylation_type not in all_dfs_ndmp_dict:
            all_dfs_ndmp_dict[methylation_type] = {}

        success_count = total_tests = 0
        # 这里 in files1[replicate_x] 是判断methylation_type是否在files1[replicate_x] 这个字典的键中存在，存在的话就说明
        # 对应的 i-both...methylation_type.txt文件存在
        if (methylation_type not in files1[replicate_x] or
                methylation_type not in files2[replicate_y]):
            print(f"    跳过甲基化类型 {methylation_type}：因为文件不存在")
            continue
        # 否则获取到当前需要处理的文件的相对路径
        file1_path = os.path.join(dir1, files1[replicate_x][methylation_type])
        file2_path = os.path.join(dir2, files2[replicate_y][methylation_type])

        # 获取两个文件中的染色体数量
        n_chromosomes_1 = get_column_count(file1_path)
        n_chromosomes_2 = get_column_count(file2_path)

        if n_chromosomes_1 is None or n_chromosomes_2 is None:
            print(f"    无法获取 {methylation_type} 的染色体数量")
            continue

        if n_chromosomes_1 != n_chromosomes_2:
            print(f"    {methylation_type} 染色体数量不一致：{n_chromosomes_1} vs {n_chromosomes_2}")
            continue

        # 到这里就说明要处理的两个both文件中都有染色体数据且总列数相同
        # 不过其实可以不用判断的，因为一开始newtoboth已经保证列数肯定相同
        n_chromosomes = n_chromosomes_1  # 获取总染色体数
        print(f"    处理甲基化类型 {methylation_type}，共 {n_chromosomes} 条染色体")

        # 处理当前x_y组的当前甲基化类型的文件的每条染色体，这里的chr_num是第几个染色体，而不是染色体号
        for chr_num in range(1, n_chromosomes + 1):
            if methylation_type in unfilter_mtypes:
                success = process_methylation_type_with_collection(
                    file1_path, file2_path, methylation_type, output_dir,
                    dir1_name, dir2_name, chr_num
                )  # 处理1次检验的一个染色体
            else:
                # 这里要区分甲基化类型
                all_df_ndmp = process_methylation_type_with_collection_pvfilter(
                    file1_path, file2_path, methylation_type, output_dir,
                    dir1_name, dir2_name, chr_num
                )   #返回此次的一个pvalue>0.05的df
                    #  all_results([pos, pvalue, change])
                if chr_num not in all_dfs_ndmp_dict[methylation_type]:
                    all_dfs_ndmp_dict[methylation_type][chr_num] = all_df_ndmp
                    if (isinstance(all_df_ndmp, pd.DataFrame)):
                        print(f"成功添加{methylation_type}-{chr_num}-特定df到字典中，此次df中有{len(all_df_ndmp)}行数据")
                success = isinstance(all_df_ndmp, pd.DataFrame)
            total_tests += 1
            if success:  # 上一步处理正常进行就会返回True，否则返回False
                success_count += 1

        # 合并FET结果并进行FDR校正
        #  其中 FET结果的格式如：pos, m1, u1, m2, u2, pvalue
        if success_count > 0:
            output_dir1 = os.path.join(output_dir, methylation_type)
            # 获取一个output_x_y/甲基化类型/目录下的所有FET文件，将其concat起来，然后FDR检验为qvalues,并导出到磁盘，最终格式为：
            #                                        Chromosome, Methylation_Type, Position, Pvalue, Qvalue
            merge_fet_results_and_fdr(output_dir1, replicate_x, replicate_y, methylation_type,all_dfs_ndmp_dict,n_chromosomes)
            # 生成DMP文件
            generate_dmp_files(dir1,dir2,output_dir1, replicate_x, replicate_y, mtype1=methylation_type,all_dfs_ndmp_dict=all_dfs_ndmp_dict
                               ,unfilter_mtypes=unfilter_mtypes,n_chromosomes=n_chromosomes)

    print(f"  组对 (wt{replicate_y}, mut{replicate_x}) 处理完成！成功 {success_count}/{total_tests} 次检验")
    return success_count, total_tests ,    # 这里是一个染色体就算一次检验


def process_all_combinations(dir1, dir2, m, n,unfilter_mtypes,work_dir="."):
    """处理所有组合，进行m*n*3次检验"""

    print(f"扫描文件目录...")
    # files就是一个链路 序号->甲基化类型->文件名 ，可以通过序号和甲基化类型获取到对应文件夹中的对应文件名
    files1 = scan_sample_files_by_replicates(dir1, m)
    files2 = scan_sample_files_by_replicates(dir2, n)

    print(f"目录1 ({dir1}) 找到 {len(files1)} 组文件")
    print(f"目录2 ({dir2}) 找到 {len(files2)} 组文件")

    # 看看是否有缺漏的both文件，若有，记录其序号
    missing_replicates1 = [i for i in range(1, m + 1) if i not in files1]
    missing_replicates2 = [i for i in range(1, n + 1) if i not in files2]
    # 输出缺漏的both文件的序号
    if missing_replicates1:
        print(f"警告：目录1缺少这些组: {missing_replicates1}")
    if missing_replicates2:
        print(f"警告：目录2缺少这些组: {missing_replicates2}")

    # 记录获取到的both文件的序号
    available_replicates1 = [i for i in range(1, m + 1) if i in files1]
    available_replicates2 = [i for i in range(1, n + 1) if i in files2]

    # 计算总共需要的比对组数
    total_combinations = len(available_replicates1) * len(available_replicates2)
    print(f"\n开始处理 {total_combinations} 个组合...")

    # rstrip(符号)用于去除右侧末尾的连续特定符号，这里特定符号为os.sep即当前系统路径分隔符：\
    # 保证右侧没有路径分隔符之后，用basename函数获取所给路径的最后一个部分（若没有去除\，则获取到的是空字符串）
    # 即获取到对应目录的名称——基因型
    dir1_name = sanitize_filename(os.path.basename(dir1.rstrip(os.sep)))
    dir2_name = sanitize_filename(os.path.basename(dir2.rstrip(os.sep)))

    total_success = total_tests = 0
    start_time = time.time()
    # enumerate从0开始枚举，i对应0开始的序号，replicate_x为available_replicates1的每一个元素，在此为每个从小到大排序的实际存在的序号
    for i, replicate_x in enumerate(available_replicates1):
        # j对应0开始的序号，replicate_y为available_replicates2的每一个元素，在此为每个从小到大排序的实际存在的序号
        for j, replicate_y in enumerate(available_replicates2):
            print(f"\n进度: {i * len(available_replicates2) + j + 1}/{total_combinations}")
            # 进行第 mut_replicate_X-wt_replicate_y 组处理
            success_count, test_count = process_replicate_pair(
                replicate_x, replicate_y, files1, files2, dir1, dir2, dir1_name, dir2_name,unfilter_mtypes,work_dir=work_dir
            )  # process_replicate_pair是处理1*3次检验
            total_success += success_count # 这里仍是一个染色体算一次检验
            total_tests += test_count  # 这里仍是一个染色体算一次检验

    end_time = time.time()
    print(f"\n所有处理完成！")
    print(f"总计: {total_success}/{total_tests} 次成功检验")
    print(f"用时: {end_time - start_time:.2f} 秒")

    # 输出说明
    print(f"单次比较结果保存在 ./output_x_y/甲基化类型/ 目录中")

    return total_success == total_tests #全部成功就return True，否则return False

def bayes_deciding(sig_count, nonsig_count):


    prob_gt_half = sig_count/(sig_count+nonsig_count)
    final_decision = 1 if prob_gt_half > 2/3 else 0
    # print(f"\n判定（阈值={recommended_threshold * 100:.0f}%）")
    # print(f"  判定结果：{'显著' if final_decision else '不显著'}")
    # print(f"  置信度：{prob_gt_half * 100:.1f}%")

    return final_decision

def find_common_significant_sites(output_dirs=None, methytype2='CpG', dir1=None, dir2=None,work_dir="."):
    """
    找出在所有组合检验中都显著的位点，并获取对应位点的相关信息
    参数：
        output_dirs: 输出目录列表，如果为None则自动扫描
        methytype2: 甲基化类型
    """

    print("\n寻找所有组合中共同显著的位点...")

    and_output_dir = os.path.join(work_dir, "and_output")
    os.makedirs(and_output_dir, exist_ok=True)

    # 1. 自动扫描所有output_x_y目录
    if output_dirs is None:
        output_dirs = glob.glob(os.path.join(work_dir,f"output_*_*/{methytype2}"))
        # 找出所有output_x_y/methytype2/这种格式的目录，一共m*n个
        output_dirs = [d for d in output_dirs if os.path.isdir(d)]

    if not output_dirs: # 没找到的话
        print("未找到任何输出目录")
        return None

    print(f"找到 {len(output_dirs)} 个输出目录")

    # 2. 一次性读取所有FDR_correct文件到内存，读取所有FDR_corrected文件（包含所有位点）,因为如果只读取显著位点信息的话，
            # 假设某个位点在1_1 sig文件中有出现，但在后面一次检验中没出现，就不知道是因为不显著没出现还是这次输入的原始数据里就没有该位点
    valid_dirs = [] # 收集在后续操作中起到作用的目录的路径
    site_statistics = {}  # 建立如此映射：{site_id: {'sig': 0, 'total': 0}}
    all_dataframes = {} # 最终可以通过all_dataframes[目录]->FDR_correct对应的df
    dir_to_replicate = {} #记录目录对应的replicate编号
    for output_dir in output_dirs: # 遍历所有output_x_y/methytype2/
        # 获取到当前甲基化目录下的FDR_correct文件，
            # 格式为：'Chromosome', 'Methylation_Type', 'Position', 'Pvalue', 'Qvalue'
        # 其中qvalue不同于之前都是小于阈值的，而是所有的
        match = re.search(r'output_wt(\d+)_mut(\d+)', output_dir)
        if not match:
            continue

        replicate_y = int(match.group(1))
        replicate_x = int(match.group(2))


        fdr_all_files = glob.glob(os.path.join(output_dir, "FDR_corrected_results_*.txt"))

        # 检验存在性
        if not fdr_all_files:
            print(f"  警告：{output_dir} 中未找到FDR_corrected文件")
            continue

        fdr_all_file = fdr_all_files[0] # 这是因为glob.glob的返回值是一个列表，所以用[0]获取到实际存在的文件路径

        try:
            df = pd.read_csv(fdr_all_file, sep='\s+') # 读取该文件到df中FDR_corrected格式为：
                                        # 'Chromosome', 'Methylation_Type', 'Position', 'Pvalue', 'Qvalue'
            if len(df) > 0:
                # 为每一行创建唯一标识符
                #  具体为：对df中的每一行，将其中的染色体号-甲基化类型-位点号提取出来，放在右侧作为新的一列site_id
                #     这是因为后面需要统计所有位点在多次检验中的显著与否，通过site_id就能仅用一个属性区分不同的位点了，不必去
                #   连续判断多个属性是否相等来筛选或访问（染色体号、甲基化类型、位点号） （注意这里的染色体号是both文件中的第几个三列）
                df['site_id'] = df.apply(
                    lambda row: f"{int(row['Chromosome'])}-{row['Methylation_Type']}-{int(row['Position'])}",
                    axis=1
                )
                all_dataframes[output_dir] = df # all_dataframes[目录]->FDR_correct对应的df
                valid_dirs.append(output_dir)
                dir_to_replicate[output_dir] = (replicate_x, replicate_y)  # 记录编号
                for _, row in df.iterrows(): # 遍历每一行，每一行是一次检验中FDR_correct中所有不管显著还是不显著的位点信息
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

                    if row['Qvalue'] < 0.05:
                        site_statistics[site_id]['sig_count'] += 1


                print(f"  {output_dir}: {len(df)} 个位点")
            else:
                print(f"  {output_dir}: 无显著位点")
        except Exception as e:
            print(f"  错误：读取 {fdr_all_file} 失败: {e}")
            continue
    valid_dirs.sort(key=lambda d: dir_to_replicate[d])
    print(f"\n统计完成，共 {len(site_statistics)} 个不同位点")

    if not site_statistics:
        print("没有找到任何有效的位点信息结果")
        return None

    # 3. 获取到贝叶斯方法判定为显著的位点，放到common_sites里
    common_sites = []
    for site_id, stats in site_statistics.items():
        if stats['total_count'] != len(valid_dirs):
            stats['total_count'] = len(valid_dirs)
        sig_count = stats['sig_count']  # 获取到当前位点显著检验次数
        nonsig_count = stats['total_count'] - sig_count  # 获取到当前位点非显著检验次数
        is_significant = bayes_deciding(sig_count, nonsig_count)
        if is_significant:
            common_sites.append(site_id)
    if not common_sites:
        print("没有在所有组合中都显著的位点")
        return None
    else:
        print(f"\n贝叶斯判定后，共 {len(common_sites)} 个显著位点")

    # 4. 读取每个目录的甲基化变化方向信息
    print("正在读取甲基化变化方向信息...")
    methylation_change_by_dir = {}  # {output_dir: {site_id: change}}

    for output_dir in valid_dirs: # 遍历上述有效的该甲基化类型的目录
        methylation_change_by_dir[output_dir] = {} # 创建当前目录的字典元素，
                                        # 构建 output_dir->site_id->change的链路
        # 查找该目录下所有的all_simple_Chr文件，因为这个文件中有change信息，其格式为pos, pvalue, change，
                                                # 而染色体号从文件名获取
        all_simple_files = glob.glob(os.path.join(output_dir, "all_simple_Chr*.txt"))
                                                        # 获取当前目录下所有的all_simple文件的路径并组成列表

        for file_path in all_simple_files: # 遍历all_simple文件，格式为pos, pvalue, change
            chr_match = re.search(r'Chr(\d+)\.txt$', file_path) #创建正则表达式及捕获组，用于获取染色体号
            if not chr_match:
                continue
            chr_num = int(chr_match.group(1)) # 通过捕获组获取到当前文件的染色体号

            try:
                df = pd.read_csv(
                    file_path,
                    sep='\t',
                    usecols=[0, 2],
                    names=['position', 'change'],
                    dtype={'position': float, 'change': float}, # float啥字符串都可以转换，不会发生错误，是比较稳妥的方案
                    skiprows=1
                )

                # 此时将pos和change都转换为int类型就不会出问题，否则str转int有可能出问题，比如"100.0"转int就会报错
                df['position'] = df['position'].astype(int)
                df['change'] = df['change'].astype(int)

                positions = df['position'].values
                changes = df['change'].values

                for position, change in zip(positions, changes):
                    site_id = f"{chr_num}-{methytype2}-{position}"
                    methylation_change_by_dir[output_dir][site_id] = change

            except Exception as e:
                print(f"  警告：读取 {file_path} 失败: {e}")
                continue

        print(f"  {output_dir}: 读取了 {len(methylation_change_by_dir[output_dir])} 个位点的变化方向")

    # 5. 处理共同位点信息
    print("正在处理共同位点详细信息...")
    common_site_details = []

    # 为每个DataFrame创建索引以加快查询
    indexed_dfs = {} # 建立output_dir->df1（site_id变作索引后）的链路
    for output_dir in valid_dirs:
        df = all_dataframes[output_dir] # 获取每个目录中的FDR_correct对应的df（已经加上了site_id这列）
        indexed_dfs[output_dir] = df.set_index('site_id') # 将site_id设置为索引，并将新的df1返回作为字典的值

    # 逐一处理common位点
    for i, site_id in enumerate(common_sites):
        if i % 10000 == 0:  # 每处理10000个位点打印进度
            print(f"  已处理 {i}/{len(common_sites)} 个位点")

        site_info = {'site_id': site_id}
        chr_num, mtype, pos = site_id.split('-')
        site_info['Chromosome'] = int(chr_num)
        site_info['Methylation_Type'] = mtype
        site_info['Position'] = int(pos)

        # 收集q值的列表，收集不同output_x_y/methytype2/目录下相同染色体相同位点号(符合当前site_id信息的那个位点)在所有检验中的q值
        qvalues = []
        # 收集不同output_x_y/methytype2/目录下相同染色体相同位点号(符合当前site_id信息的那个位点)在所有检验中的变化方向
        change_values = []
        qvalue_dict = {}
        for output_dir in valid_dirs:
            replicate_x, replicate_y = dir_to_replicate[output_dir]
            col_name = f'qvalue_{dir2}{replicate_y}_{dir1}{replicate_x}'   # 新增：列名，前面是野生型的序号，后面是突变型的序号
            indexed_df = indexed_dfs[output_dir] # output_dir->df1（site_id变作索引后）,获取到索引为site_id的当前甲基化的某个FDR_correct文件
                                # 内容格式为：'Chromosome', 'Methylation_Type', 'Position', 'Pvalue', 'Qvalue'
            if site_id in indexed_df.index:
                qval = indexed_df.loc[site_id, 'Qvalue']
                qvalues.append(qval)
                qvalue_dict[col_name] = qval
                # 获取该位点在该目录中的change值，(output_dir->site_id->change的链路)
                if site_id in methylation_change_by_dir[output_dir]:
                    change_values.append(methylation_change_by_dir[output_dir][site_id])
            else:
                qvalue_dict[col_name] = 1.0

        if qvalues:  # 确保有q值数据
            qvalues_sig = [q for q in qvalues if q < 0.05]
            if len(qvalues_sig) > 0:
                site_info['Sig_Mean_Qvalue'] = np.mean(qvalues_sig)
            else:
                site_info['Sig_Mean_Qvalue'] = 1
            # site_info['Max_Qvalue'] = np.max(qvalues)
            # site_info['Min_Qvalue'] = np.min(qvalues)
            site_info['Num_Comparisons'] = len(qvalues)

            # 投票计算甲基化变化方向
            if change_values:
                # 统计change==1的次数
                num_hyper = sum(change_values)
                total_comparisons = len(change_values)
                hyper_ratio = num_hyper / total_comparisons

                # 多数投票：>= 50%则记为1（高甲基化），否则为0
                site_info['Methylation_Change'] = 1 if hyper_ratio >= 0.5 else 0
                site_info['Hyper_Count'] = num_hyper  # 高甲基化次数
                site_info['Hypo_Count'] = total_comparisons - num_hyper  # 低甲基化次数
                site_info['Hyper_Ratio'] = hyper_ratio  # 高甲基化比例
            else:
                # 如果没有change信息，标记为缺失（不过按理来说一次检验是正好对应一个qvalue和一个change的）
                site_info['Methylation_Change'] = -1  # -1表示无法确定
                site_info['Hyper_Count'] = 0
                site_info['Hypo_Count'] = 0
                site_info['Hyper_Ratio'] = 0

            site_info.update(qvalue_dict) # 添加所有replicate的qvalue列

            common_site_details.append(site_info) # 这个列表中记录了一个个字典，每个字典是一个site_id对应的各个信息，格式：
        # site_id-染色体号-甲基化类型-位点号-总检验次数-change-hypercount-hypocount-hyperratio-Sig_Mean_Qvalue-所有qvalues
            # （此处染色体号是both文件中的第几个3列）

    print(f"  完成处理 {len(common_site_details)} 个位点")

    # 6. 生成结果DataFrame
    result_df = pd.DataFrame(common_site_details)
    result_df = result_df.sort_values(['Methylation_Type', 'Chromosome', 'Position']) # 排序

    # 调整列顺序，将Methylation_Change放在更靠前的位置
    column_order = [
        'Chromosome', 'Methylation_Type', 'Position',
        'Methylation_Change', 'Hyper_Ratio', 'Hyper_Count', 'Hypo_Count', 'Num_Comparisons',
        'Sig_Mean_Qvalue'
    ]
    # 获取所有replicate列名并排序
    # 新的列名格式：dir2_y_dir1_x，需要获取包含下划线分隔数字的列
    replicate_columns = sorted(
        [col for col in result_df.columns if col.startswith('qvalue_')],
        key=lambda x: tuple(map(int, re.findall(r'\d+', x)))  # 提取所有数字并排序
    )
    column_order = column_order + replicate_columns
    result_df = result_df[column_order]

    # 保存结果
    output_file = os.path.join(and_output_dir, f"{methytype2}-final_significant_sites_DMPs.txt")
    result_df.to_csv(output_file, sep='\t', index=False)
    print(f"\n共同显著位点已保存至: {output_file}")

    # 输出统计信息
    print("\n共同显著位点统计:")
    mtype = methytype2
    mtype_df = result_df
    count = len(mtype_df)
    hyper_count = len(mtype_df[mtype_df['Methylation_Change'] == 1])
    hypo_count = len(mtype_df[mtype_df['Methylation_Change'] == 0])
    unknown_count = len(mtype_df[mtype_df['Methylation_Change'] == -1])

    print(f"  {mtype}: {count} 个位点")
    print(f"    - 高甲基化(Change=1): {hyper_count} ({hyper_count / count * 100:.1f}%)")
    print(f"    - 低甲基化(Change=0): {hypo_count} ({hypo_count / count * 100:.1f}%)")
    if unknown_count > 0:
        print(f"    - 未知(Change=-1): {unknown_count} ({unknown_count / count * 100:.1f}%)")

    return result_df # 其格式为：'Chromosome', 'Methylation_Type', 'Position',
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
    对甲基化位点数据进行滑动窗口分析
    参数：
    input_data : 是m*n*3次中一次的某一条染色体的(N)dmp文件转换为的DataFrame，包含列: ['position', 'pvalue', 'change']
    window_size :滑动窗口的大小
    step_ratio : 步长比例（窗口大小的百分比）
    save_files : 是否保存结果文件
    output_identifier : 输出文件前缀，如果save_files=True且未提供则自动生成
    """

    # 1. 数据读取和预处理
    df = input_data.copy()

    # 确保列名正确
    expected_cols = ['position', 'pvalue', 'change']
    if not all(col in df.columns for col in expected_cols):
        raise ValueError(f"DataFrame必须包含列: {expected_cols}")

    if output_identifier is None:
        output_identifier = "sliding_window_analysis"

    # 数据验证和清理
    df = df.dropna()
    if len(df) == 0:
        raise ValueError("没有有效的数据行")

    # 按位点号排序
    df = df.sort_values('position').reset_index(drop=True) #重排dataframe索引，且不将原来的索引变作新的一列

    print(f"数据预处理完成，共 {len(df)} 个位点")

    # 2. 滑动窗口分析
    positions = df['position'].values #获取所有位点号（前面已经排了序了）
    changes = df['change'].values #获取所有甲基化变化方向

    last_pos = positions[-1] #获取到最后一个位点号
    step_size = int(window_size * step_ratio) #计算得到窗口每一步位移的大小

    # 生成所有窗口的起始位置
    window_starts = np.arange(0, last_pos + step_size, step_size) #中间+step_size是为了使得最后一个区间覆盖last_pos这个位点

    print(f"窗口配置: 大小={window_size}, 步长={step_size}, 窗口数={len(window_starts)}")

    #创建列表用于存储后续每个区间的首尾位点号和不同甲基化变化方向的数量以及总显著位点数
    results = []

    for i, start in enumerate(window_starts): # start每次赋值为每个窗口的起始点位点号-1（比如第一次是0）
        if i % 1000 == 0:  # 进度提示
            print(f"处理进度: {i}/{len(window_starts)}")

        end = start + window_size  # 根据窗口宽度计算得到窗口的终止点

        # 使用numpy的searchsorted进行快速查找（相当于二分查找）
        left_idx = np.searchsorted(positions, start, side='right') # 此处意为查找positions数组中第一个大于start的元素的下标
        right_idx = np.searchsorted(positions, end, side='right') # 同理，查找positions数组中第一个大于end的元素的下标

        # 在窗口内的位点
        window_changes = changes[left_idx:right_idx] # 获取到pos在[start,end)宽度为window_size的那些位点的change构成的数组
                                                    # 注意这里的left_idx和right_idx都是下标
                            #left_idx:right_idx这个范围内都是pos位点号元素值>start&&<=end的，
                         # 但是总的列表长度数量基本上不是window_size，会少很多，因为很多位点应该没有数据
        # 统计不同类别的数量
        num_change_1 = np.sum(window_changes == 1)   # 统计hyper的数量
        num_change_0 = np.sum(window_changes == 0)  # 统计hypo的数量，或者用 len(window_changes) - num_change_1

        results.append({
            'window_start': start + 1,  # start+1才是每个窗口真正的位点号的起始点
            'window_end': end,   #start+1到end刚好window_size个位点被统计
            'count_change_1': num_change_1,
            'count_change_0': num_change_0,
            'total_count': num_change_1 + num_change_0 # 当前行对应的区间内的所有位点的甲基化变化方向（突变型相对于野生型）
                                                # 也是当前行对应区间的所有位点的总显著位点数量，因为一次检验会有一个变化方向
                                            #  而因为是从DMP文件中读来的，所有都是显著的
        })

    # 转换为DataFrame
    sliding_results = pd.DataFrame(results)

    # 3. 标准化处理
    max_count = sliding_results['total_count'].max()  # 所有区间中最大的那个总显著位点数
    if max_count == 0:
        max_count = 1  # 避免除零

    standardized_results = sliding_results.copy()
    standardized_results['standardized_count'] = sliding_results['total_count'] / max_count #求得当前区间总显著位点数与所有区间中
                                                                                        # 最大的那个显著位点数的比值
    # 选择需要的列用于标准化输出
    standardized_results = standardized_results[[
        'window_start', 'window_end', 'total_count', 'standardized_count'
    ]]

    print(f"滑动窗口分析完成，共生成 {len(sliding_results)} 个窗口")
    print(f"最大计数: {max_count}")

    # 4. 保存文件
    if save_files:
        # 滑动窗口结果
        sliding_file = f"slidingW_{output_identifier}.txt"
        sliding_file = os.path.join(outputdir1, sliding_file)
        sliding_results[['window_start', 'window_end', 'count_change_1', 'count_change_0']].to_csv(
            sliding_file,
            sep='\t',
            index=False,
            header=False
        )

        # 标准化结果
        std_file = f"noTitle_allDMCs_new_Standardized_slidingW_{output_identifier}.txt"
        std_file = os.path.join(outputdir1, std_file)
        # 格式化输出以匹配原始C++风格的对齐
        standardized_results.to_csv(
            std_file,
            sep='\t',
            index=False,
            header=False,
            float_format='%.6f'  # 控制浮点数精度
        )

        print(f"结果已保存:")
        print(f"  滑动窗口: {sliding_file}")
        print(f"  标准化结果: {std_file}")

    return sliding_results, standardized_results
    # 前者格式为window_start,window_end,count_change_1,count_change_0,total_count
                        #后者没有count_change，而有standardized_count，即每个区间显著位点数与最大的那个显著位点数的比值

def process_common_sites_sliding_window(common_sites_df=None,
                                        window_size=1000000,
                                        step_ratio=0.05,
                                        methytype='CpG',
                                        work_dir="."):
    """
    对共同显著位点进行滑动窗口分析
    参数：
    common_sites_df : 共同显著位点数据，如果为None则自动加载
    window_size : 滑动窗口大小
    step_ratio : 步长比例
    methytype : 甲基化类型
    """

    print(f"\n开始对共同显著位点进行滑动窗口分析...")

    and_output_dir = os.path.join(work_dir, "and_output")
    os.makedirs(and_output_dir, exist_ok=True)

    # 1. 加载共同显著位点数据
    if common_sites_df is None:
        common_file = os.path.join(and_output_dir, f"{methytype}-final_significant_sites_DMPs.txt")
        if not os.path.exists(common_file):
            print(f"错误：共同显著位点文件不存在 {common_file}")
            return None
        common_sites_df = pd.read_csv(common_file, sep='\t')
        print(f"从文件加载共同显著位点: {len(common_sites_df)} 个位点")

    if common_sites_df.empty:
        print("没有共同显著位点数据")
        return None

    # 2. 按染色体分组
    results = {}

    # 按染色体分组处理
    chr_groups = common_sites_df.groupby('Chromosome') #获取一个迭代器，可以迭代获取每个染色体号及其对应的子df

    for chr_num, chr_data in chr_groups: # 迭代获取每个染色体号及其对应的子df
        print(f"\n  处理染色体 {chr_num}: {len(chr_data)} 个位点")

        # 准备滑动窗口分析的数据，保持和all_simple_chr文件相同的格式以确保顺利进行
        window_data = pd.DataFrame({
            'position': chr_data['Position'].astype(int),
            'pvalue': chr_data['Sig_Mean_Qvalue'],
            'change': chr_data['Methylation_Change']
        })

        # 排序
        window_data = window_data.sort_values('position').reset_index(drop=True)

        # 执行滑动窗口分析，调用之前的函数就行
        try:
            sliding_results, std_results = sliding_window_analysis(
                window_data,
                window_size=window_size,
                step_ratio=step_ratio,
                save_files=True,
                output_identifier=f"common_sites_{methytype}_Chr{chr_num}",
                outputdir1=and_output_dir
            )

            # 收集结果文件
            results[chr_num] = {
                'sliding_results': sliding_results,
                'standardized_results': std_results,
                'input_data': window_data
            }

            print(f"    完成染色体 {chr_num}: {len(sliding_results)} 个窗口")

        except Exception as e:
            print(f"    错误：处理染色体 {chr_num} 失败: {e}")
            continue

    print(f"结果保存在以 'common_sites_{methytype}_Chr' 开头的文件中")

    return results

def find_max_total_in_outputs(output_dirs, methylation_type):
    """
    查找所有输出目录中指定甲基化类型的最大total值
    参数：
        output_dirs: 输出目录列表
        methylation_type: 甲基化类型（CpG, CHH, CHG）
    返回：
        max_total: 最大total值
    """
    max_total = 0

    for out_dir in output_dirs:
        mtype_dir = os.path.join(out_dir, methylation_type)
        if not os.path.exists(mtype_dir):
            continue

        # 查找所有DMP标准化文件
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
                print(f"    警告: 读取 {std_file} 时出错: {e}")
                continue

    print(f"  {methylation_type} 类型所有染色体中最大的total值为: {max_total}")
    return max_total


def plot_methylation_sliding_windows(output_dir=None, chr_series=None,work_dir="."):
    """
    对所有滑动窗口结果进行可视化并保存到磁盘
    使用全局max_total进行标准化，使不同染色体具有可比性
    参数：
    output_dir : 指定输出目录，如果为None则自动扫描所有output_x_y目录
    chr_series : 染色体映射Series
    """

    matplotlib.use('Agg')  # 使用非交互式后端

    # 设置中文字体
    plt.rcParams['font.sans-serif'] = ['SimHei']
    plt.rcParams['axes.unicode_minus'] = False

    print("\n开始生成甲基化滑动窗口可视化图表")

    # 扫描所有输出目录
    if output_dir is None:
        output_dirs = glob.glob(os.path.join(work_dir, "output_*_*"))
        output_dirs = [d for d in output_dirs if os.path.isdir(d)]
    else:
        output_dirs = [output_dir]

    if not output_dirs:
        print("未找到输出目录")
        return

    total_plots = 0
    methylation_types = ['CpG', 'CHH', 'CHG']

    # 为每种甲基化类型分别计算全局max_total
    max_totals = {}
    for mtype in methylation_types:
        max_totals[mtype] = find_max_total_in_outputs(output_dirs, mtype)
        if max_totals[mtype] == 0:
            print(f"  警告: {mtype} 类型未找到有效的total值，将使用1作为默认值")
            max_totals[mtype] = 1

    for out_dir in output_dirs:
        print(f"\n处理目录: {out_dir}")

        for mtype in methylation_types:
            mtype_dir = os.path.join(out_dir, mtype)
            if not os.path.exists(mtype_dir):
                continue

            print(f"  处理甲基化类型: {mtype}（使用全局max_total={max_totals[mtype]}）")

            # 查找当前甲基化目录下所有DMP滑动窗口文件
            dmp_sliding_files = glob.glob(os.path.join(mtype_dir, "slidingW_DMP_*.txt"))

            # 按前缀分组处理，同一前缀的所有染色体绘制在一张大图上
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

            # 对每个前缀组进行处理
            for prefix, chr_nums in prefix_groups.items():
                # 对染色体编号排序
                chr_nums = sorted(chr_nums, key=lambda x: int(x))

                # 收集所有染色体的数据
                all_chrom_data = []
                chrom_names = []

                for chr_num in chr_nums:
                    # 构建对应的标准化文件路径
                    chr_name = get_chr_name(chr_num, chr_series)

                    # 修正：使用正确的变量名构建文件路径
                    dmp_sliding_file = os.path.join(mtype_dir, f"slidingW_DMP_{prefix}_Chr{chr_num}.txt")
                    dmp_std_file = os.path.join(mtype_dir,
                                                f"noTitle_allDMCs_new_Standardized_slidingW_DMP_{prefix}_Chr{chr_num}.txt")
                    ndmp_std_file = os.path.join(mtype_dir,
                                                 f"noTitle_allDMCs_new_Standardized_slidingW_N-DMP_{prefix}_Chr{chr_num}.txt")

                    if not all(os.path.exists(f) for f in [dmp_sliding_file, dmp_std_file, ndmp_std_file]):
                        print(f"    警告: Chr{chr_num} 的文件不完整，跳过")
                        continue

                    try:
                        # 读取数据
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
                            print(f"    警告: Chr{chr_num} 的数据为空，跳过")
                            continue

                        # 使用全局max_total重新计算比率
                        max_total = max_totals[mtype]

                        x = (sliding_df['start'] + sliding_df['end']) / 2

                        # 重新计算所有比率（使用全局max_total）
                        y_dmp = (sliding_df['hyper'] + sliding_df['hypo']) / max_total
                        y_hyper = sliding_df['hyper'] / max_total
                        y_hypo = sliding_df['hypo'] / max_total
                        y_ndmp = ndmp_std_df['ndmp_normalized']  # NDMP保持原样

                        # 处理可能的长度不一致
                        max_len = max(len(x), len(y_dmp), len(y_hyper), len(y_hypo), len(y_ndmp))
                        x = x.reindex(range(max_len), fill_value=0)
                        y_dmp = y_dmp.reindex(range(max_len), fill_value=0)
                        y_hyper = y_hyper.reindex(range(max_len), fill_value=0)
                        y_hypo = y_hypo.reindex(range(max_len), fill_value=0)
                        y_ndmp = y_ndmp.reindex(range(max_len), fill_value=0)

                        # 存储数据
                        all_chrom_data.append({
                            'x': x,
                            'y_dmp': y_dmp,
                            'y_hyper': y_hyper,
                            'y_hypo': y_hypo,
                            'y_ndmp': y_ndmp
                        })
                        chrom_names.append(chr_name)

                        print(f"    成功加载染色体 {chr_name} 的数据")

                    except Exception as e:
                        print(f"    处理 Chr{chr_num} 时出错: {e}")
                        continue

                # 如果有数据，绘制大图
                if all_chrom_data:
                    try:
                        # 创建大图，每个染色体一个子图
                        n_chromosomes = len(all_chrom_data)
                        fig, axes = plt.subplots(n_chromosomes, 1, figsize=(15, 3 * n_chromosomes))

                        # 如果只有一个染色体，axes不是数组，需要转换为数组
                        if n_chromosomes == 1:
                            axes = [axes]

                        # 设置总标题
                        fig.suptitle(f'{mtype} Methylation Analysis - {prefix} (Global Normalized)',
                                     fontsize=16, fontfamily='Times New Roman')

                        # 绘制每个染色体的子图
                        for idx, (chrom_data, chrom_name) in enumerate(zip(all_chrom_data, chrom_names)):
                            ax = axes[idx]

                            # 绘制所有数据线
                            ax.plot(chrom_data['x'], chrom_data['y_dmp'], label='DMP', color='red', linewidth=1.5)
                            ax.plot(chrom_data['x'], chrom_data['y_hyper'], label='Hyper-ratio', color='green',
                                    linewidth=1)
                            ax.plot(chrom_data['x'], chrom_data['y_hypo'], label='Hypo-ratio', color='blue',
                                    linewidth=1)
                            ax.plot(chrom_data['x'], chrom_data['y_ndmp'], label='NDMP', color='darkgray',
                                    linewidth=1.5)

                            ax.set_ylim(bottom=0)
                            ax.set_ylabel('Ratio', fontsize=10, fontfamily='Times New Roman')

                            # 设置子图标题
                            ax.set_title(f'{chrom_name}', fontsize=18, fontfamily='Times New Roman', pad=20,y=-0.4)

                            # 添加网格
                            ax.grid(True, alpha=0.3)

                            # 只在第一个子图添加完整图例
                            if idx == 0:
                                ax.legend(loc='upper right', ncol=2, fontsize=8, framealpha=0.7)

                        # 调整布局
                        plt.tight_layout(rect=[0, 0, 1, 0.96])  # 为总标题留出空间

                        # 保存图片
                        plot_filename = os.path.join(mtype_dir,
                                                     f"methylation_plot_{mtype}_{prefix}_all_chromosomes.png")
                        plt.savefig(plot_filename, dpi=300, bbox_inches='tight')
                        plt.close()

                        total_plots += 1
                        print(f"    成功生成大图: {mtype}_{prefix} -> {os.path.basename(plot_filename)}")

                    except Exception as e:
                        print(f"    绘制大图时出错: {e}")
                        continue

    print(f"\n图表生成完成！共生成 {total_plots} 张大图")


def plot_common_sites_sliding_windows(methytype='CpG', chr_series=None, work_dir="."):
    """
    对共同显著位点的滑动窗口结果进行可视化并保存到磁盘
    使用全局max_total进行标准化，将所有染色体绘制在一张大图上
    """
    matplotlib.use('Agg')

    # 设置全局字体为Times New Roman
    plt.rcParams['font.family'] = 'Times New Roman'
    plt.rcParams['axes.unicode_minus'] = False

    print(f"\n开始生成共同显著位点({methytype})的滑动窗口可视化图表...")

    and_output_dir = os.path.join(work_dir, "and_output")

    # 先找到该甲基化类型的全局max_total
    print(f"  正在查找 {methytype} 的全局max_total...")
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
            print(f"    警告: 读取 {std_file} 时出错: {e}")

    if max_total == 0:
        print(f"  警告: 未找到有效的total值，使用默认值1")
        max_total = 1
    else:
        print(f"  {methytype} 的全局max_total为: {max_total}")

    # 读取DMR数据 - 根据甲基化类型选择对应的DMR文件
    print("  正在读取DMR数据...")
    dmr_file = os.path.join(and_output_dir, f"{methytype}-final_significant_regions_DMRs.txt")
    dmr_data = {}
    if os.path.exists(dmr_file):
        try:
            dmr_df = pd.read_csv(dmr_file, sep='\s+')

            for _, row in dmr_df.iterrows():
                try:
                    chrom = str(row['Chromosome'])  # 基于列名访问
                    direction = int(row['Direction'])  # 基于列名访问
                    start = int(row['DMR_start'])  # 基于列名访问
                    end = int(row['DMR_end'])  # 基于列名访问

                    # 计算中点
                    mid = (start + end) / 2

                    # 提取染色体数字部分
                    chrom_num = str(chrom).replace('Chr', '').replace('chr', '')

                    if chrom_num not in dmr_data:
                        dmr_data[chrom_num] = []

                    dmr_data[chrom_num].append((mid, direction))
                except (ValueError, IndexError):
                    continue
            print(f"  成功加载 {sum(len(dmrs) for dmrs in dmr_data.values())} 个DMR")
        except Exception as e:
            print(f"  读取DMR文件时出错: {e}")
    else:
        print(f"  警告: DMR文件 {dmr_file} 不存在")

    # 动态获取染色体列表
    sliding_files = glob.glob(os.path.join(and_output_dir, f"slidingW_common_sites_{methytype}_Chr*.txt"))

    # 从文件名中提取染色体编号并排序
    chromosomes = []
    for file in sliding_files:
        match = re.search(r'slidingW_common_sites_.+_Chr(\d+)\.txt$', os.path.basename(file))
        if match:
            chr_num = match.group(1)
            if chr_num not in chromosomes:
                chromosomes.append(chr_num)

    # 按数字顺序排序染色体
    chromosomes.sort(key=int)

    if not chromosomes:
        print(f"未找到共同显著位点({methytype})的滑动窗口文件")
        return

    print(f"找到 {len(chromosomes)} 个染色体的滑动窗口文件")

    all_chrom_data = []
    chrom_names = []

    for chr_num in chromosomes:
        sliding_file = os.path.join(and_output_dir, f"slidingW_common_sites_{methytype}_Chr{chr_num}.txt")
        std_file = os.path.join(and_output_dir,
                                f"noTitle_allDMCs_new_Standardized_slidingW_common_sites_{methytype}_Chr{chr_num}.txt")

        if not all(os.path.exists(f) for f in [sliding_file, std_file]):
            print(f"    警告: Chr{chr_num} 的文件不完整，跳过")
            continue

        try:
            # 读取数据
            sliding_df = pd.read_csv(sliding_file, sep='\s+', header=None, names=['start', 'end', 'hyper', 'hypo'])
            std_df = pd.read_csv(std_file, sep='\s+', header=None, names=['start', 'end', 'total', 'normalized'])

            if sliding_df.empty or std_df.empty:
                print(f"    警告: Chr{chr_num} 的数据为空，跳过")
                continue

            if len(sliding_df) != len(std_df):
                print(f"    警告: Chr{chr_num} 的数据长度不一致，跳过")
                continue

            # 新增：读取 output_1_1 的 NDMP 数据
            ndmp_file = os.path.join(work_dir, "output_wt1_mut1", methytype,
                                     f"noTitle_allDMCs_new_Standardized_slidingW_N-DMP_wt_replicate1_mut_replicate1_Chr{chr_num}.txt")
            y_ndmp = None  # 初始化为None
            if os.path.exists(ndmp_file):
                try:
                    ndmp_df = pd.read_csv(ndmp_file, sep='\s+', header=None,
                                          names=['start', 'end', 'total', 'ndmp_normalized'])
                    if not ndmp_df.empty:
                        y_ndmp = ndmp_df['ndmp_normalized']
                        print(f"    成功读取 output_1_1 的 NDMP 数据: Chr{chr_num}")
                except Exception as e:
                    print(f"    警告: 读取 output_1_1 NDMP 数据失败 (Chr{chr_num}): {e}")
            else:
                print(f"    提示: output_1_1 的 NDMP 文件不存在 (Chr{chr_num})")

            # 使用全局max_total重新计算比率
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

            # 存储数据
            chr_name = get_chr_name(chr_num, chr_series)
            all_chrom_data.append({
                'x': x,
                'y_total': y_total,
                'y_hyper': y_hyper,
                'y_hypo': y_hypo,
                'y_ndmp': y_ndmp
            })
            chrom_names.append(chr_name)
            print(f"    成功加载染色体 {chr_name} 的数据")

        except Exception as e:
            print(f"    处理 Chr{chr_num} 时出错: {e}")
            continue

    # 如果有数据，绘制大图
    if all_chrom_data:
        try:
            # 创建大图，包含多个子图
            n_chromosomes = len(all_chrom_data)
            fig, axes = plt.subplots(n_chromosomes, 1, figsize=(15, 3 * n_chromosomes))

            # 如果只有一个染色体，axes不是数组，需要转换为数组
            if n_chromosomes == 1:
                axes = [axes]

            # 设置总标题
            fig.suptitle(f'Distribution of Common Significant Sites - {methytype} context',
                         fontsize=16, fontfamily='Times New Roman')

            # 绘制每个染色体的子图
            for idx, (chrom_data, chrom_name) in enumerate(zip(all_chrom_data, chrom_names)):
                ax = axes[idx]

                # 绘制DMP数据
                ax.plot(chrom_data['x'], chrom_data['y_total'], label='DMP', color='red', linewidth=2)
                ax.plot(chrom_data['x'], chrom_data['y_hyper'], label='hyper-methylation', color='green', linewidth=1.5)
                ax.plot(chrom_data['x'], chrom_data['y_hypo'], label='hypo-methylation', color='blue', linewidth=1.5)


                if chrom_data['y_ndmp'] is not None:
                    ax.plot(chrom_data['x'], chrom_data['y_ndmp'], label='NDMP',
                            color='darkgray', linewidth=1.5)
                # 统一设置Y轴范围为0到1.2
                ax.set_ylim(0, 1.2)

                # 获取染色体编号（从chrom_name中提取）
                chrom_num = chrom_name.replace('chr', '')

                # 添加DMR标记
                if chrom_num in dmr_data:
                    for mid, direction in dmr_data[chrom_num]:
                        # 根据direction选择颜色：1=hyper(绿色), 0=hypo(蓝色)
                        color = 'green' if direction == 1 else 'blue'
                        # 在DMR中点位置添加竖线，显示在Y轴1.0到1.2的范围内
                        ax.axvline(x=mid, ymin=0.9, ymax=1, color=color, linewidth=2, alpha=0.7)

                # 设置子图标题和标签
                ax.text(0.5, -0.2, f"{chrom_name}",
                        transform=ax.transAxes,
                        fontfamily='Times New Roman',
                        ha='center', va='top',
                        fontsize=15)

                # 添加网格
                ax.grid(True, alpha=0.3)

                # 只在第一个子图添加完整图例
                if idx == 0:
                    # 创建自定义图例条目，包括DMR标记
                    from matplotlib.lines import Line2D
                    legend_elements = [
                        Line2D([0], [0], color='red', linewidth=2, label='DMP'),
                        Line2D([0], [0], color='green', linewidth=1.5, label='hyper-DMP'),
                        Line2D([0], [0], color='blue', linewidth=1.5, label='hypo-DMP'),
                        Line2D([0], [0], color='green', linewidth=2, label='hyper-DMR'),
                        Line2D([0], [0], color='blue', linewidth=2, label='hypo-DMR')
                    ]
                    # 如果有NDMP数据，添加到图例
                    if chrom_data['y_ndmp'] is not None:
                        legend_elements.append(
                            Line2D([0], [0], color='darkgray', linewidth=1.5,
                                  label='NDMP')
                        )
                    ax.legend(handles=legend_elements, loc='upper right', bbox_to_anchor=(1.0, 0.95),
                              ncol=5, fontsize=8, framealpha=0.7)

            # 调整布局
            plt.tight_layout(rect=[0, 0, 1, 0.96])  # 为总标题留出空间

            # 保存图片
            plot_filename = os.path.join(and_output_dir, f"common_sites_plot_{methytype}_all_chromosomes.png")
            plt.savefig(plot_filename, dpi=300, bbox_inches='tight')
            plt.close()

            print(f"成功生成大图: {methytype} -> {os.path.basename(plot_filename)}")

        except Exception as e:
            print(f"绘制大图时出错: {e}")
    else:
        print(f"未找到 {methytype} 的有效数据")



def rename_chromosome_files(chr_series, work_dir="."):
    """
    批量重命名所有包含染色体号的文件，将Chr数字替换为真实染色体名称

    参数：
        chr_series: 染色体映射Series (染色体名称 -> 索引)
        work_dir: 工作目录
    """
    print("\n开始批量重命名染色体文件...")

    if chr_series is None or len(chr_series) == 0:
        print("错误：染色体映射为空，跳过重命名")
        return

    # 创建反向映射：数字索引 -> 染色体名称
    # chr_series 的 index 是染色体名称，value 是数字索引
    index_to_chr = {i: chr_series.index[i] for i in range(len(chr_series))}

    print(f"染色体映射: {index_to_chr}")

    # 定义需要搜索的目录模式
    search_dirs = []

    # 添加所有 output_x_y 目录
    output_dirs = glob.glob(os.path.join(work_dir, "output_*_*"))
    search_dirs.extend([d for d in output_dirs if os.path.isdir(d)])

    # 添加 and_output 目录
    and_output_dir = os.path.join(work_dir, "and_output")
    if os.path.exists(and_output_dir):
        search_dirs.append(and_output_dir)

    if not search_dirs:
        print("未找到需要处理的目录")
        return

    print(f"将在 {len(search_dirs)} 个目录中搜索文件")

    # 统计信息
    total_renamed = 0
    failed_renames = 0

    # 定义匹配染色体号的正则表达式模式
    # 匹配 Chr 后跟数字的模式，如 Chr1, Chr12 等
    chr_pattern = re.compile(r'(.*?)Chr(\d+)(.*?)$')

    # 遍历所有目录
    for search_dir in search_dirs:
        print(f"\n处理目录: {search_dir}")

        # 递归遍历目录中的所有文件
        for root, dirs, files in os.walk(search_dir):
            for filename in files:
                # 检查文件名是否包含 Chr数字 模式
                match = chr_pattern.match(filename)

                if match:
                    prefix = match.group(1)  # Chr 之前的部分
                    chr_num = int(match.group(2))  # 染色体数字
                    suffix = match.group(3)  # Chr数字 之后的部分

                    # 根据 chr_series 获取真实染色体名称
                    # 注意：文件名中的 Chr1 对应 index 0
                    chr_index = chr_num - 1

                    if chr_index not in index_to_chr:
                        print(f"  警告: Chr{chr_num} 不在映射表中，跳过文件 {filename}")
                        continue

                    real_chr_name = index_to_chr[chr_index]

                    # 构造新文件名
                    # 如果染色体名称本身包含 'chr' 前缀，直接使用
                    # 否则使用 Chr 前缀
                    if real_chr_name.lower().startswith('chr'):
                        chr_part = real_chr_name
                    else:
                        chr_part = f"Chr{real_chr_name}"

                    new_filename = f"{prefix}{chr_part}{suffix}"

                    # 如果新旧文件名相同，跳过
                    if filename == new_filename:
                        continue

                    # 构造完整路径
                    old_path = os.path.join(root, filename)
                    new_path = os.path.join(root, new_filename)

                    # 检查目标文件是否已存在
                    if os.path.exists(new_path):
                        print(f"  警告: 目标文件已存在，跳过重命名: {filename} -> {new_filename}")
                        failed_renames += 1
                        continue

                    # 执行重命名
                    try:
                        os.rename(old_path, new_path)
                        total_renamed += 1
                        print(f"  {filename} -> {new_filename}")
                    except Exception as e:
                        print(f" 重命名失败: {filename} -> {new_filename}, 错误: {e}")
                        failed_renames += 1

    # 输出统计信息
    print(f"\n重命名完成！")
    print(f"  成功重命名: {total_renamed} 个文件")
    if failed_renames > 0:
        print(f"  失败或跳过: {failed_renames} 个文件")


def convert_output_to_csv(work_dir="."):
    """
    将 and_output 目录下的 final DMP 和 final DMR 文件转换为 CSV 格式（逗号分隔）

    参数：
        work_dir: 工作目录
    """

    and_output_dir = os.path.join(work_dir, "and_output")

    if not os.path.exists(and_output_dir):
        print(f"错误：目录 {and_output_dir} 不存在")
        return 0

    # 定义需要转换的文件模式
    file_patterns = [
        "*-final_significant_sites_DMPs.txt",  # final DMP 文件
        "*-final_significant_regions_DMRs.txt"  # final DMR 文件
    ]

    converted_count = 0
    failed_count = 0

    for pattern in file_patterns:
        # 搜索匹配的文件
        matching_files = glob.glob(os.path.join(and_output_dir, pattern))

        for txt_file in matching_files:
            try:
                # 读取制表符分隔的文件
                df = pd.read_csv(txt_file, sep=r'\s+')

                if df.empty:
                    print(f"  跳过空文件: {os.path.basename(txt_file)}")
                    continue

                # 生成 CSV 文件路径（将 .txt 替换为 .csv）
                csv_file = txt_file.replace('.txt', '.csv')

                # 保存为逗号分隔的 CSV 文件
                df.to_csv(csv_file, sep=',', index=False)

                print(f"  成功转换： {os.path.basename(txt_file)} -> {os.path.basename(csv_file)}")
                converted_count += 1

            except Exception as e:
                print(f"  转换失败: {os.path.basename(txt_file)}, 错误: {e}")
                failed_count += 1

    print(f"\\n转换完成！")
    print(f"  成功转换: {converted_count} 个文件")
    if failed_count > 0:
        print(f"  转换失败: {failed_count} 个文件")

    return converted_count


def convert_chromosome_to_names(chr_series, work_dir="."):
    """
    将 and_output 目录下 final DMP 和 final DMR 文件中的 Chromosome 列
    从数值转换为真实的染色体名称

    参数：
        chr_series: 染色体映射 Series (染色体名称 -> 索引，索引从0开始)
        work_dir: 工作目录

    返回：
        成功转换的文件数量
    """

    if chr_series is None or len(chr_series) == 0:
        print("错误：染色体映射为空，跳过转换")
        return 0

    and_output_dir = os.path.join(work_dir, "and_output")

    if not os.path.exists(and_output_dir):
        print(f"错误：目录 {and_output_dir} 不存在")
        return 0

    # 创建数值到染色体名称的映射
    # 文件中的 Chromosome 列是从1开始的数字，需要转换为 chr_series 中的染色体名称
    # chr_series.index[0] 对应文件中的 1
    index_to_chr = {i + 1: chr_series.index[i] for i in range(len(chr_series))}

    print(f"染色体映射表: {index_to_chr}")

    # 定义需要处理的文件模式
    file_patterns = [
        "*-final_significant_sites_DMPs.txt",  # final DMP 文件
        "*-final_significant_regions_DMRs.txt"  # final DMR 文件
    ]

    converted_count = 0
    failed_count = 0

    for pattern in file_patterns:
        # 搜索匹配的文件
        matching_files = glob.glob(os.path.join(and_output_dir, pattern))

        for file_path in matching_files:
            try:
                # 读取文件
                df = pd.read_csv(file_path, sep=r'\s+')

                if df.empty:
                    print(f"  跳过空文件: {os.path.basename(file_path)}")
                    continue

                # 检查是否有 Chromosome 列
                if 'Chromosome' not in df.columns:
                    print(f"  警告: {os.path.basename(file_path)} 中没有 Chromosome 列，跳过")
                    continue

                # 保存原始的 Chromosome 列用于调试
                original_chrs = df['Chromosome'].unique()

                # 转换 Chromosome 列
                # 先确保是整数类型
                df['Chromosome'] = df['Chromosome'].astype(int)

                # 使用映射转换为染色体名称
                df['Chromosome'] = df['Chromosome'].map(index_to_chr)

                # 检查是否有未成功映射的值
                if df['Chromosome'].isna().any():
                    unmapped_count = df['Chromosome'].isna().sum()
                    print(f"  警告: {os.path.basename(file_path)} 中有 {unmapped_count} 个染色体编号无法映射")
                    # 可选：删除无法映射的行
                    df = df.dropna(subset=['Chromosome'])

                # 保存回原文件（覆盖）
                df.to_csv(file_path, sep='\t', index=False)

                print(f"    成功转换: {os.path.basename(file_path)}")
                print(f"    原始编号: {sorted(original_chrs)}")
                print(f"    转换后: {sorted(df['Chromosome'].unique())}")
                converted_count += 1

            except Exception as e:
                print(f"  ✗ 转换失败: {os.path.basename(file_path)}")
                print(f"    错误信息: {e}")
                import traceback
                traceback.print_exc()
                failed_count += 1

    # 输出统计信息
    print(f"\n染色体编号转换完成！")
    print(f"  成功转换: {converted_count} 个文件")
    if failed_count > 0:
        print(f"  转换失败: {failed_count} 个文件")

    return converted_count


def main():
    start_time = time.time()
    print("程序说明：")
    print("进行m×n×3次Fisher检验")
    print("每组包含3种甲基化类型：CpG, CHH, CHG")
    print()

    try:
        # 获取命令行参数
        m = int(sys.argv[2])
        n = int(sys.argv[1])
        dir1 = sys.argv[4]
        dir2 = sys.argv[3]
        biotype = int(sys.argv[5])
        #dir1 = input("请输入第一个样本目录名：").strip()
        #dir2 = input("请输入第二个样本目录名：").strip()
        #m = int(input("请输入第一个目录的组数 (m)："))
        #n = int(input("请输入第二个目录的组数 (n)："))
        #biotype = int(input("请输入所给基因型源自的生物类型（0-动物，1-植物，2-不过滤）"))

    # 因为目录名随便输都行，不会出现错误的情况，所以如果错了就肯定是组数m,n输入错了
    except ValueError:
        print("错误：请输入有效的数字")
        sys.exit(1)

    # 验证目录存在，且确保组数合法
    if not os.path.exists(dir1):
        print(f"错误：目录 '{dir1}' 不存在！")
        sys.exit(1)
    if not os.path.exists(dir2):
        print(f"错误：目录 '{dir2}' 不存在！")
        sys.exit(1)
    if m <= 0 or n <= 0:
        print("错误：组数必须大于0")
        sys.exit(1)

    print(f"\n参数确认:")
    print(f"第一个目录: {dir1} (包含 {m} 组文件)")
    print(f"第二个目录: {dir2} (包含 {n} 组文件)")
    print(f"预计进行: {m * n * 3} 次Fisher检验")

    print("\n第一阶段：newtoboth进行中")
    # 将bismark新格式数据转换为both格式
    chr_series = newtoboth(m, n, dir1, dir2)
    if biotype == 0:
        unfilter_mtypes = ["CHH","CHG"]
    elif biotype == 1:
        unfilter_mtypes = ["CpG"]
    elif biotype ==2:
        unfilter_mtypes = ["CHH","CHG","CpG"]
    else:
        print("错误：生物类型必须是0、1或2")
        sys.exit(1)
    print(f"不需要p值预过滤的甲基化类型: {unfilter_mtypes}")
    success = process_all_combinations(dir1, dir2, m, n,unfilter_mtypes)  # process_all_combinations是进行m*n*3次检验

    if success: # 全部检验都成功
        print("\n所有检验和FDR校正均成功完成！")
        methylation_types = ['CpG', 'CHH', 'CHG']
        for mtype in methylation_types:
            common_sites_df = find_common_significant_sites(methytype2=mtype, dir1=dir1, dir2=dir2)
                             # 其格式为：'Chromosome', 'Methylation_Type', 'Position',
                                    # 'Methylation_Change', 'Hyper_Ratio', 'Hyper_Count', 'Hypo_Count',
                                # 'Sig_Mean_Qvalue', 'Max_Qvalue', 'Min_Qvalue', 'Num_Comparisons'
            if common_sites_df is not None and not common_sites_df.empty:
                print(f"找到 {len(common_sites_df)} 个 {mtype} 类型的共同显著位点")

                # 进行滑动窗口分析
                print(f"开始对 {mtype} 共同显著位点进行滑动窗口分析...")

                results = process_common_sites_sliding_window(
                    common_sites_df=common_sites_df,
                    methytype=mtype
                )

                print(f"完成 {mtype} 类型的滑动窗口分析")
            else:
                print(f"未找到 {mtype} 类型的共同显著位点")



        print("开始 DMR 分析流程")
        process_common_sites_dmr_and_summarize(
            dir1=dir1,
            dir2=dir2,
            m=m,
            n=n,
            methylation_types=methylation_types,
        )
        # 生成所有滑动窗口的可视化图表
        plot_methylation_sliding_windows(chr_series=chr_series)
        for mtype111 in ["CpG","CHH","CHG"]:
            plot_common_sites_sliding_windows(mtype111, chr_series=chr_series)
        convert_chromosome_to_names(chr_series=chr_series, work_dir=".")
        rename_chromosome_files(chr_series=chr_series, work_dir=".")
        convert_output_to_csv(work_dir=".")
        print(f"- DMP 结果：output_x_y/甲基化类型/")
        print(f"- 共同显著位点：and_output/")
        print(f"- 最终显著 DMR：and_output/*-final_significant_regions_DMRs.txt")
    else:
        print("\n部分检验失败，请检查输出信息。")
    end_time = time.time()  # 记录结束时间
    elapsed_time = end_time - start_time  # 计算耗时
    print(f"总耗时: {elapsed_time:.2f} 秒")


if __name__ == "__main__":
    main()



