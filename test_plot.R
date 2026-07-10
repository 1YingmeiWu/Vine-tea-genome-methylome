import pandas as pd
import matplotlib.pyplot as plt
import re
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests

def read_bedgraph(file_path):
    """读取bedgraph文件，跳过第一行，返回第四列的数值（过滤非数值数据）"""
    try:
        df = pd.read_csv(
            file_path,
            sep=r'\s+',
            header=None,
            usecols=[3],
            dtype={3: float},
            skiprows=1
        )
        values = df[3].dropna().values
        return values.astype(float)
    except Exception as e:
        print(f"读取文件 {file_path} 出错: {e}")
        return None

def main():
    # -------------------------- 1. 基础配置与数据读取 --------------------------
    input_files = [
        "L.CHG.bedgraph", "L.CHH.bedgraph", "L.CpG.bedgraph",
        "R.CHG.bedgraph", "R.CHH.bedgraph", "R.CpG.bedgraph",
        "S.CHG.bedgraph", "S.CpG.bedgraph", "S.CHH.bedgraph"
    ]
    pal = ["#E64B35", "#4DBBD5", "#00A087"]
    color_map = {"R": pal[0], "S": pal[1], "L": pal[2]}
    group_order = ["R", "S", "L"]  # 样本组
    type_order = ["CpG", "CHG", "CHH"]  # 甲基化类型

    # 读取并整合数据
    data = []
    for file in input_files:
        match = re.match(r'^([RSL])\.([CHGHCpG]+)\.bedgraph$', file)
        if match:
            group = match.group(1)
            methyl_type = match.group(2)
            values = read_bedgraph(file)
            if values is not None and len(values) > 0:
                for val in values:
                    data.append({
                        'group': group,
                        'type': methyl_type,
                        'exp': val
                    })
    df = pd.DataFrame(data)
    if df.empty:
        print("警告：没有有效数据可处理，请检查输入文件！")
        return


    # -------------------------- 2. 输出数据分布统计文件（修复quantile参数） --------------------------
    # 关键修复：删除quantile()中的skipna=True（低版本pandas不支持该参数，默认自动跳过空值）
    group_stats = df.groupby(['group', 'type']).agg(
        number=('exp', lambda x: sum(~pd.isna(x))),  # 有效数据量
        mean=('exp', lambda x: round(x.mean(), 2)),  # 简化：mean()默认skipna=True，低版本兼容
        Q25=('exp', lambda x: round(x.quantile(q=0.25, interpolation='linear'), 2)),  # 删除skipna=True
        Q50=('exp', lambda x: round(x.quantile(q=0.5, interpolation='linear'), 2)),   # 删除skipna=True
        Q75=('exp', lambda x: round(x.quantile(q=0.75, interpolation='linear'), 2))    # 删除skipna=True
    ).reset_index()

    # 保存数据分布文件
    stats_filename = "methylation_multi_group_stat.xls"
    group_stats.to_csv(stats_filename, sep='\t', index=False, quoting=0)
    print(f"数据分布统计文件已保存：{stats_filename}")


    # -------------------------- 3. 输出组间显著性分析文件 --------------------------
    ref_group = group_order[0]  
    sig_results = []

    for methyl_type in type_order:
        single_type_data = df[df['type'] == methyl_type].copy()
        if len(single_type_data) == 0:
            print(f"警告：类型 {methyl_type} 无有效数据，跳过显著性检验")
            continue

        # 提取参考组数据（简化：dropna()显式去空，兼容低版本）
        ref_data = single_type_data[single_type_data['group'] == ref_group]['exp'].dropna()
        if len(ref_data) == 0:
            print(f"警告：参考组 {ref_group}（类型 {methyl_type}）无有效数据，跳过该类型检验")
            continue

        # 循环检验非参考组
        for test_group in group_order:
            if test_group == ref_group:
                continue

            test_data = single_type_data[single_type_data['group'] == test_group]['exp'].dropna()
            if len(test_data) == 0:
                print(f"警告：检验组 {test_group}（类型 {methyl_type}）无有效数据，跳过该组检验")
                continue

            # Wilcoxon检验
            stat, p_value = mannwhitneyu(
                x=test_data,
                y=ref_data,
                alternative='two-sided',
                use_continuity=True
            )

            sig_results.append({
                'type': methyl_type,
                'ref_group': ref_group,
                'test_group': test_group,
                'statistic': round(stat, 4),
                'p_value': round(p_value, 6)
            })

    # 整理显著性结果
    sig_df = pd.DataFrame(sig_results)
    if not sig_df.empty:
        # 多重检验校正
        sig_df['p_adjusted'] = multipletests(
            pvals=sig_df['p_value'],
            alpha=0.05,
            method='bonferroni'
        )[1]

        # 显著性标记
        def get_sig_mark(p):
            if p <= 0.001:
                return '***'
            elif p <= 0.01:
                return '**'
            elif p <= 0.05:
                return '*'
            else:
                return 'ns'
        sig_df['significance'] = sig_df['p_adjusted'].apply(get_sig_mark)

        # 保存显著性文件
        sig_filename = "methylation_intergroup_signif.xls"
        sig_df.to_csv(sig_filename, sep='\t', index=False, quoting=0)
        print(f"组间显著性分析文件已保存：{sig_filename}")
    else:
        print("警告：无有效显著性检验结果，未生成显著性分析文件")


    # -------------------------- 4. 箱线图绘制 --------------------------
    plt.figure(figsize=(8, 6))
    positions = []
    labels = []
    current_pos = 1
    box_width = 0.15
    intra_group_spacing = 0.25
    inter_group_spacing = 0.5

    for methyl_type in type_order:
        type_data = df[df['type'] == methyl_type]
        type_positions = []
        for group in group_order:
            group_data = type_data[type_data['group'] == group]['exp']
            if not group_data.empty:
                plt.boxplot(
                    group_data,
                    positions=[current_pos],
                    widths=box_width,
                    patch_artist=True,
                    boxprops=dict(facecolor=color_map[group], edgecolor='black'),
                    medianprops=dict(color='black'),
                    showfliers=True
                )
                type_positions.append(current_pos)
                current_pos += intra_group_spacing
        positions.append(type_positions)
        labels.append(methyl_type)
        current_pos += inter_group_spacing

    # 计算每个甲基化类型的中心位置作为x轴标签
    x_ticks = [sum(pos)/len(pos) for pos in positions]
    plt.xticks(x_ticks, labels, fontsize=12)
    plt.xlabel('Methylation Type', fontsize=14)
    plt.ylabel('Methylation Level', fontsize=14)
    plt.title('', fontsize=16)

    # 图例放在右下角，显示样本组
    legend_elements = [plt.Rectangle((0,0),1,1, facecolor=color_map[group], edgecolor='black') 
                      for group in group_order]
    plt.legend(legend_elements, group_order, title='Sample', loc='lower right')

    # 保存图片
    plt.tight_layout()
    plt.savefig('methylation_boxplot.pdf', dpi=300, bbox_inches='tight')
    plt.savefig('methylation_boxplot.png', dpi=300, bbox_inches='tight')
    plt.show()
    print("箱线图已保存：methylation_boxplot.pdf / methylation_boxplot.png")

if __name__ == "__main__":
    main()
