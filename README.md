**pDMPcaller**是一款轻量级DNA甲基化差异位点（DMP）和差异区域（DMR）识别工具，专为动植物甲基化组数据分析设计。支持CpG/CHH/CHG三种甲基化类型检测，基于Fisher精确检验结合贝叶斯判定模型，实现高可信度的差异位点筛选，并提供可视化结果输出。


## 核心功能
- 🧬 **多类型甲基化分析**：支持CpG、CHH、CHG位点的差异检测
- 📊 **统计检验**：整合Fisher精确检验与FDR校正，降低假阳性
- 🎯 **DMR识别**：基于滑动窗口算法识别差异甲基化区域
- 📈 **可视化输出**：自动生成甲基化分布图谱、染色体定位图
- 📝 **标准化结果**：输出兼容下游分析的文本格式（BED/CSV）


## 安装指南

### 环境依赖
- Python 3.8及以上版本
- 依赖库：`pandas`、`numpy`、`scipy`、`matplotlib`


### 快速安装
1. 克隆仓库到本地：
   ```bash
   git clone https://github.com/jiantaoyuNWAFU/pDMPcaller.git
   cd pDMPcaller
   ```

2. 安装依赖库：
   ```bash
   pip install -r requirements.txt
   ```


## 快速示例（5分钟上手）
我们提供测试数据（位于`tests/test_data/`），可直接运行验证功能：
```bash
# 1. 进入主程序目录
cd bin

# 2. 运行工具并输入测试参数
python run_pDMPcaller.py
```
运行后按提示输入：
- 第一组样本目录：`../tests/test_data/control`
- 第二组样本目录：`../tests/test_data/treatment`
- 第一组样本数：`1`
- 第二组样本数：`1`
- 生物类型：`1`（植物）

运行完成后，结果会输出到`and_output/`目录。


## 使用教程

### 输入数据格式
需提供**两组样本的甲基化数据目录**（如对照组vs处理组），每组包含若干样本文件（制表符分隔，无表头），列定义如下：

| 列名           | 说明                          | 示例       |
|----------------|-------------------------------|------------|
| 染色体号       | 基因组染色体ID                | chr1       |
| 位点位置       | 甲基化位点的基因组坐标        | 10260657   |
| 甲基化读段数   | 该位点被甲基化的测序读段数    | 5          |
| 非甲基化读段数 | 该位点未被甲基化的测序读段数  | 3          |
| 甲基化类型     | 位点类型（CpG/CHH/CHG）       | CpG        |

示例样本文件（`sample1.txt`）：
```
chr1	10260657	0	3	CpG
chr1	11796376	0	18	CpG
chr2	18978618	0	3	CpG
```


### 运行命令
```bash
# 进入主程序目录
cd bin

# 启动工具
python run_pDMPcaller.py
```


## 输出结果说明

### 核心结果目录
运行完成后生成2个核心目录：

1. **中间结果目录（`output_wtX_mutY/`）**
   - `{type}_DMPs_fisher.txt`：Fisher检验原始结果（含P值）
   - `{type}_DMPs_fdr.txt`：FDR校正后的差异位点
   - `{type}_sliding_window_results.txt`：滑动窗口分析数据

2. **最终结果目录（`and_output/`）**
   - `{type}-final_significant_sites_DMPs.txt`：最终显著DMP位点（核心结果）
   - `{type}-final_significant_regions_DMRs.txt`：显著DMR区域
   - `common_sites_plot_{type}_all_chromosomes.png`：甲基化位点染色体分布可视化图


### 结果字段说明（以`CpG-final_significant_sites_DMPs.txt`为例）
| 列序 | 字段名               | 说明                                                                 |
|------|----------------------|----------------------------------------------------------------------|
| 1    | 染色体号             | 甲基化位点所在染色体                                                 |
| 2    | 位点位置             | 基因组坐标                                                           |
| 3    | 甲基化读段数         | 第一组样本在该位点的甲基化读段总数                                   |
| 4    | 非甲基化读段数       | 第一组样本在该位点的非甲基化读段总数                                 |
| 5    | 甲基化类型           | CpG/CHH/CHG                                                          |
| 6    | Fisher检验P值        | 两组样本甲基化差异的显著性P值                                        |
| 7    | FDR校正后P值         | 多重检验校正后的P值（通常以<0.05为显著阈值）                         |
| 8    | 是否显著             | True=显著差异位点，False=非显著                                     |


## 可视化结果示例
![CpG位点染色体分布](docs/example_output/common_sites_plot_CpG_all_chromosomes.png)
*图：CpG差异甲基化位点在各染色体上的分布*


## 常见问题
1. **运行报错“文件找不到”**
   检查输入的样本目录路径是否正确（建议使用相对路径，如`../tests/test_data/control`）。

2. **结果文件为空**
   可能是样本数据量过小，或两组样本甲基化差异不显著；可尝试调整`src/pDMPcaller.py`中的`fdr_threshold`参数（默认0.05）。

3. **可视化图不显示**
   确保`matplotlib`已正确安装，或在可视化模块代码中添加`plt.show()`语句。


## 许可证
本项目基于MIT许可证开源，详见[LICENSE](LICENSE)文件。


## 联系方式
- 作者：你的姓名/ID
- 邮箱：你的邮箱地址
- 问题反馈：[GitHub Issues](https://github.com/你的GitHub用户名/pDMPcaller/issues)


## 更新日志
### v1.0（2025-12-04）
- 首次发布，支持CpG/CHH/CHG位点的DMP/DMR检测
- 实现基础可视化与标准化结果输出
```
