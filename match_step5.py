#!/usr/bin/env python3
import sys
import os
import re
import csv

def parse_orthogroups(orthogroup_file):
    """
    解析Orthogroups.tsv文件，提取三种关系：
    1. 多对一：多个橘子基因对应一个拟南芥基因
    2. 一对多：一个橘子基因对应多个拟南芥基因
    3. 一对一：一个橘子基因对应一个拟南芥基因

    返回：橘子基因到拟南芥基因的映射字典
    """
    print("正在解析Orthogroups.tsv文件...")

    # 存储映射关系：橘子基因 -> [拟南芥基因1, 拟南芥基因2, ...]
    citrus_to_arabidopsis = {}

    with open(orthogroup_file, 'r') as f:
        # 读取表头
        header = f.readline().strip()
        columns = header.split('\t')

        # 找到橘子基因列和拟南芥基因列的索引
        citrus_col_idx = None
        arabidopsis_col_idx = None

        for i, col_name in enumerate(columns):
            # 根据您的描述，橘子基因列包含"Cau."
            if 'Cau.' in col_name or 'Citrus' in col_name or 'Orange' in col_name:
                citrus_col_idx = i
            # 拟南芥基因列包含"gene_"或"AT"
            elif 'gene_' in col_name or 'AT' in col_name or 'Arabidopsis' in col_name:
                arabidopsis_col_idx = i

        # 如果没找到，尝试按位置猜测
        if citrus_col_idx is None:
            citrus_col_idx = 1  # 默认第二列是橘子基因
        if arabidopsis_col_idx is None:
            arabidopsis_col_idx = 2  # 默认第三列是拟南芥基因

        print(f"橘子基因列索引: {citrus_col_idx}, 拟南芥基因列索引: {arabidopsis_col_idx}")

        for line_num, line in enumerate(f, start=2):
            line = line.strip()
            if not line:
                continue

            parts = line.split('\t')
            if len(parts) <= max(citrus_col_idx, arabidopsis_col_idx):
                continue

            og_id = parts[0]
            citrus_genes_str = parts[citrus_col_idx].strip()
            arabidopsis_genes_str = parts[arabidopsis_col_idx].strip()

            # 提取橘子基因列表（以Cau.开头）
            citrus_genes = []
            if citrus_genes_str:
                # 处理逗号分隔的列表
                for gene in citrus_genes_str.split(','):
                    gene = gene.strip()
                    if gene.startswith('Cau.'):
                        citrus_genes.append(gene)

            # 提取拟南芥基因列表（以gene_开头或AT开头）
            arabidopsis_genes = []
            if arabidopsis_genes_str:
                # 处理逗号分隔的列表
                for gene in arabidopsis_genes_str.split(','):
                    gene = gene.strip()
                    # 提取gene_AT3G56230或AT3G56230格式
                    if gene.startswith('gene_'):
                        # 提取AT3G56230部分
                        match = re.search(r'gene_(AT[A-Za-z0-9]+)', gene)
                        if match:
                            arabidopsis_genes.append(match.group(1))
                    elif gene.startswith('AT'):
                        arabidopsis_genes.append(gene)

            # 存储映射关系
            for citrus_gene in citrus_genes:
                if citrus_gene not in citrus_to_arabidopsis:
                    citrus_to_arabidopsis[citrus_gene] = []

                for arabidopsis_gene in arabidopsis_genes:
                    if arabidopsis_gene not in citrus_to_arabidopsis[citrus_gene]:
                        citrus_to_arabidopsis[citrus_gene].append(arabidopsis_gene)

    print(f"解析完成，共找到 {len(citrus_to_arabidopsis)} 个橘子基因的映射关系")
    return citrus_to_arabidopsis

def generate_intermediate_file(citrus_to_arabidopsis, output_file):
    """生成中间映射文件"""
    print("生成中间映射文件...")

    with open(output_file, 'w') as f:
        f.write("Citrus_Gene\tArabidopsis_Genes\n")

        for citrus_gene, arabidopsis_genes in sorted(citrus_to_arabidopsis.items()):
            # 用分号分隔多个拟南芥基因
            arabidopsis_str = ";".join(arabidopsis_genes)
            f.write(f"{citrus_gene}\t{arabidopsis_str}\n")

    print(f"中间文件已保存: {output_file}")

def detect_delimiter(file_path):
    """检测文件的分隔符（逗号或制表符）"""
    with open(file_path, 'r') as f:
        first_line = f.readline().strip()
        if ',' in first_line and '\t' not in first_line:
            return ','
        else:
            return '\t'

def process_high_variable_genes(hv_file, citrus_to_arabidopsis, output_file):
    """处理高变基因文件，添加拟南芥基因列，支持CSV和TSV格式"""
    print("处理高变基因文件...")

    # 检测文件分隔符
    delimiter = detect_delimiter(hv_file)
    print(f"检测到文件分隔符: {'逗号' if delimiter == ',' else '制表符'}")

    with open(hv_file, 'r') as f_in, open(output_file, 'w') as f_out:
        # 使用csv模块读取文件，处理可能的引号
        reader = csv.reader(f_in, delimiter=delimiter)
        writer = csv.writer(f_out, delimiter='\t')

        # 处理表头
        header = next(reader)
        header.append("Arabidopsis_Genes")
        writer.writerow(header)

        # 处理数据行
        matched_count = 0
        unmatched_count = 0

        for row in reader:
            if not row:
                continue

            citrus_gene = row[0].strip()
            arabidopsis_genes = citrus_to_arabidopsis.get(citrus_gene, [])

            # 用分号分隔多个拟南芥基因，如果没有匹配到则显示斜杠
            if arabidopsis_genes:
                arabidopsis_str = ";".join(arabidopsis_genes)
                matched_count += 1
            else:
                arabidopsis_str = "/"
                unmatched_count += 1

            # 添加新列
            row.append(arabidopsis_str)

            # 写入新行
            writer.writerow(row)

    print(f"结果文件已保存: {output_file}")
    print(f"匹配成功: {matched_count} 个基因")
    print(f"未匹配: {unmatched_count} 个基因")

def main():
    # 检查输入参数
    if len(sys.argv) != 5:
        print("用法: python map_genes_complex.py <orthogroups.tsv> <high_variable_genes.csv/tsv> <intermediate_file.txt> <output_file.txt>")
        print("示例: python map_genes_complex.py Orthogroups.tsv high_var_genes.csv intermediate.txt result.txt")
        sys.exit(1)

    orthogroup_file = sys.argv[1]
    hv_file = sys.argv[2]
    intermediate_file = sys.argv[3]
    output_file = sys.argv[4]

    # 检查文件是否存在
    for file_path in [orthogroup_file, hv_file]:
        if not os.path.exists(file_path):
            print(f"错误: 文件不存在 - {file_path}")
            sys.exit(1)

    # 步骤1: 解析Orthogroups.tsv文件，提取三种关系
    citrus_to_arabidopsis = parse_orthogroups(orthogroup_file)

    # 步骤2: 生成中间文件
    generate_intermediate_file(citrus_to_arabidopsis, intermediate_file)

    # 步骤3: 处理高变基因文件
    process_high_variable_genes(hv_file, citrus_to_arabidopsis, output_file)

    print("完成!")

if __name__ == "__main__":
    main()
