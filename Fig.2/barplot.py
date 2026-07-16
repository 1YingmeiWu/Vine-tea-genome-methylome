import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from collections import defaultdict
import numpy as np
from multiprocessing import Pool, cpu_count
import functools

def process_column(args):
    """
    处理单个列的统计（用于多线程）
    """
    col_name, col_data = args
    # 过滤掉0值
    filtered_data = col_data[col_data != 0]
    
    # 统计大于等于70的个数（hyper）
    hyper_count = (filtered_data >= 70).sum()
    # 统计小于70的个数（hypo）
    hypo_count = (filtered_data < 70).sum()
    
    return col_name, {'hyper': hyper_count, 'hypo': hypo_count}

def process_single_file(input_file, output_file, plot_base, title, chunksize=100000, n_threads=None):
    """
    处理单个文件，统计第4列以后每列大于等于70和小于70的个数，生成叠加柱状图
    
    参数:
    input_file -- 输入文件路径
    output_file -- 统计结果输出文件路径
    plot_base -- 图表文件基础名称（不含扩展名）
    title -- 图表标题
    chunksize -- 分块读取大小
    n_threads -- 线程数量
    """
    
    if n_threads is None:
        n_threads = cpu_count()
    
    print(f"处理文件: {input_file}")
    print(f"使用 {n_threads} 线程处理")
    
    try:
        # 分块读取大文件
        if input_file.endswith(('.xls', '.xlsx')):
            # Excel文件一次性读取（通常不会太大）
            if input_file.endswith('.xlsx'):
                df = pd.read_excel(input_file, engine='openpyxl')
            else:
                df = pd.read_excel(input_file, engine='xlrd')
        else:
            # 文本文件分块读取
            print(f"使用分块读取模式，块大小: {chunksize}行")
            
            # 先读取列名
            try:
                header_df = pd.read_csv(input_file, sep='\t', nrows=0, engine='python')
            except:
                header_df = pd.read_csv(input_file, sep=',', nrows=0, engine='python')
            
            value_columns = header_df.columns[3:]
            print(f"分析的列: {list(value_columns)}")
            
            # 初始化统计结果
            column_stats = {col: {'hyper': 0, 'hypo': 0} for col in value_columns}
            
            # 分块处理
            chunk_count = 0
            if input_file.endswith('.tsv') or input_file.endswith('.txt'):
                chunks = pd.read_csv(input_file, sep='\t', chunksize=chunksize, engine='python')
            else:
                chunks = pd.read_csv(input_file, sep=',', chunksize=chunksize, engine='python')
            
            for chunk in chunks:
                chunk_count += 1
                if chunk_count % 10 == 0:
                    print(f"已处理 {chunk_count} 个数据块...")
                
                # 多线程处理当前块的所有列
                with Pool(n_threads) as pool:
                    results = pool.map(process_column, [(col, chunk[col]) for col in value_columns])
                
                # 累加统计结果
                for col_name, stats in results:
                    column_stats[col_name]['hyper'] += stats['hyper']
                    column_stats[col_name]['hypo'] += stats['hypo']
            
            print(f"总共处理 {chunk_count} 个数据块")
            
            # 创建结果DataFrame
            result_df = pd.DataFrame(column_stats).T
            # 确保列顺序：hypo在前，hyper在后（这样hypo在底部）
            result_df = result_df[['hypo', 'hyper']]
        
        # 如果不是分块处理（Excel文件或小文件）
        if 'result_df' not in locals():
            # 检查文件是否有足够的列
            if df.shape[1] < 5:
                print(f"错误: 文件列数不足，至少需要5列")
                return
            
            # 获取第4列及以后的列名（从索引3开始）
            value_columns = df.columns[3:]
            print(f"分析的列: {list(value_columns)}")
            
            # 多线程处理所有列
            print("使用多线程处理列统计...")
            with Pool(n_threads) as pool:
                results = pool.map(process_column, [(col, df[col]) for col in value_columns])
            
            # 创建结果字典
            stats_data = {col_name: stats for col_name, stats in results}
            
            # 创建结果DataFrame
            result_df = pd.DataFrame(stats_data).T
            # 确保列顺序：hypo在前，hyper在后（这样hypo在底部）
            result_df = result_df[['hypo', 'hyper']]
        
        # 保存统计结果
        result_df.to_csv(output_file, sep='\t')
        print(f"统计结果已保存到: {output_file}")
        print("\n统计结果:")
        print(result_df)
        
        # 绘制叠加柱状图
        if not result_df.empty:
            plot_stacked_bar(result_df, plot_base, title)
        else:
            print("没有有效数据可绘制图表")
            
    except Exception as e:
        print(f"处理文件 {input_file} 时出错: {str(e)}")
        import traceback
        traceback.print_exc()

def plot_stacked_bar(result_df, plot_base, title):
    """
    绘制叠加柱状图
    """
    # 设置颜色：蓝色（hypo）在底部，红色（hyper）在上部
    colors = ['#33D6E2', '#F97373']  # 蓝色(hypo), 红色(hyper)
    
    # 绘图 - 不使用网格风格，设置为默认白色背景
    plt.style.use('default')
    fig, ax = plt.subplots(figsize=(max(6, len(result_df) * 0.8), 7))
    
    # 绘制叠加柱状图 - 调整bar_width使柱子变细
    bar_width = 0.6
    
    # 由于列顺序已经是hypo在前，hyper在后，所以hypo会显示在底部
    result_df.plot(kind='bar', stacked=True, ax=ax, 
                  color=colors, width=bar_width)
    
    # 去除所有网格线
    ax.grid(False)
    
    # 添加数量标签
    for i, (idx, row) in enumerate(result_df.iterrows()):
        current_height = 0
        for j, (col_name, value) in enumerate(row.items()):
            if value > 0:
                fontsize = 8 if value > 100000 else 9  # 大数值使用较小字体
                ax.text(i, current_height + value/2, f'{int(value):,}',
                        ha='center', va='center', fontsize=fontsize, fontweight='bold')
                current_height += value
    # 设置标题和标签
    ax.set_title(title, fontsize=15)
    ax.set_xlabel('', fontsize=12)
    ax.set_ylabel('Number of site', fontsize=12)
    
    # 旋转x轴标签
    plt.xticks(rotation=0, ha='center', fontsize=10)
    
    # 调整图例 - 确保图例顺序与柱状图一致
    ax.legend(['Hypo (<0.7)', 'Hyper (≥0.7)'], title='Methylation Type', 
              bbox_to_anchor=(1.05, 1), loc='upper left')
    
    # 调整布局
    plt.tight_layout()
    
    # 保存为PNG和PDF两种格式
    png_file = f"{plot_base}.png"
    pdf_file = f"{plot_base}.pdf"
    plt.savefig(png_file, dpi=300, bbox_inches='tight')
    plt.savefig(pdf_file, dpi=300, bbox_inches='tight')
    print(f"PNG图表已保存到: {png_file}")
    print(f"PDF图表已保存到: {pdf_file}")
    plt.close()

if __name__ == "__main__":
    import sys
    if len(sys.argv) < 5:
        print("用法: python methylation_stats.py <输入文件路径> <统计结果文件路径> <图表基础名称> <图表标题> [块大小] [线程数]")
        print("示例: python methylation_stats.py ./data.txt ./result.txt ./methylation_plot 'Methylation Level Distribution'")
        print("示例: python methylation_stats.py ./big_data.txt ./result.txt ./plot 'Distribution' 100000 8")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    plot_base = sys.argv[3]
    title = sys.argv[4]
    
    # 可选参数
    chunksize = int(sys.argv[5]) if len(sys.argv) > 5 else 100000
    n_threads = int(sys.argv[6]) if len(sys.argv) > 6 else cpu_count()
    
    if not os.path.isfile(input_file):
        print(f"错误: 输入文件 '{input_file}' 不存在")
        sys.exit(1)
    
    process_single_file(input_file, output_file, plot_base, title, chunksize, n_threads)
