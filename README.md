# MultiDMPcaller
DNA甲基化差异位点（DMP）/差异区域（DMR）识别工具，支持CpG/CHH/CHG类型分析。

## 环境配置
```bash
pip install -r requirements.txt
```

## 运行方式
1. 将数据放在任意目录，确保格式为：染色体号\t位点位置\t甲基化读段数\t非甲基化读段数\t甲基化类型
2. 运行核心代码：
```bash
cd src
python MultiDMPcaller.py
```
3. 按提示输入数据目录、样本数、生物类型（0=动物/1=植物/2=不过滤）即可。

## 示例数据
`example_data/1-msv-CpG-chr1.txt`,`example_data/1-wt-CpG-chr1.txt` 为测试数据，包含基础甲基化位点格式，可直接用于测试运行。

## 输出说明
运行后生成 `and_output/` 目录，包含显著DMP/DMR位点文件及可视化图表。
