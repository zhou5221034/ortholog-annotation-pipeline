#!/usr/bin/env python3
import csv
import sys

def match_genes_and_transcripts(gene_file, transcript_file, output_file):
    """
    匹配基因文件和比对文件，将转录本ID添加到基因文件后面
    未匹配到的基因在最后一列显示斜杠 /
    """
    
    # 第一步：读取比对文件，建立基因名到转录本ID的映射
    gene_to_transcript = {}
    
    try:
        with open(transcript_file, 'r', encoding='utf-8') as f:
            for line in f:
                line = line.strip()
                if not line:  # 跳过空行
                    continue
                
                # 使用制表符分割（比对文件看起来是制表符分隔）
                parts = line.split('\t')
                if len(parts) >= 2:
                    transcript_id = parts[0].strip()
                    gene_name = parts[1].strip()
                    
                    # 将基因名作为键，转录本ID作为值存入字典
                    gene_to_transcript[gene_name] = transcript_id
    
    except FileNotFoundError:
        print(f"错误：找不到比对文件 {transcript_file}")
        return
    except Exception as e:
        print(f"读取比对文件时出错：{e}")
        return
    
    print(f"成功读取比对文件，共找到 {len(gene_to_transcript)} 个基因的转录本信息")
    
    # 第二步：读取基因CSV文件并添加转录本信息
    matched_count = 0
    unmatched_count = 0
    
    try:
        with open(gene_file, 'r', encoding='utf-8') as f_in, \
             open(output_file, 'w', newline='', encoding='utf-8') as f_out:
            
            # 读取CSV文件
            reader = csv.reader(f_in)
            writer = csv.writer(f_out)
            
            # 读取表头
            header = next(reader)
            
            # 在表头末尾添加新列
            new_header = header + ['transcript_id']
            writer.writerow(new_header)
            
            # 处理每一行数据
            for row in reader:
                if not row:  # 跳过空行
                    continue
                
                gene_name = row[0]  # 第一列是基因名
                
                # 查找对应的转录本ID，未找到则显示斜杠 /
                transcript_id = gene_to_transcript.get(gene_name, '/')
                
                if transcript_id != '/':
                    matched_count += 1
                else:
                    unmatched_count += 1
                
                # 将转录本ID添加到行末尾
                new_row = row + [transcript_id]
                writer.writerow(new_row)
    
    except FileNotFoundError:
        print(f"错误：找不到基因文件 {gene_file}")
        return
    except Exception as e:
        print(f"处理文件时出错：{e}")
        return
    
    # 第三步：输出统计信息
    total = matched_count + unmatched_count
    print(f"\n匹配完成！")
    print(f"输出文件：{output_file}")
    print(f"总基因数：{total}")
    print(f"成功匹配：{matched_count}")
    print(f"未匹配：{unmatched_count}")
    if total > 0:
        print(f"匹配率：{matched_count/total*100:.2f}%")
    else:
        print("匹配率：0.00%")

def main():
    """主函数，通过命令行参数接收文件"""
    # 检查命令行参数
    if len(sys.argv) != 4:
        print("用法: python match_genes.py <基因CSV文件> <比对TXT文件> <输出文件>")
        print("示例: python match_genes.py Total_markergenes.csv result.txt gene_with_transcript.csv")
        sys.exit(1)
    
    gene_csv_file = sys.argv[1]
    transcript_txt_file = sys.argv[2]
    output_csv_file = sys.argv[3]
    
    print("开始基因匹配程序...")
    print(f"基因文件：{gene_csv_file}")
    print(f"比对文件：{transcript_txt_file}")
    print(f"输出文件：{output_csv_file}")
    print("-" * 50)
    
    # 执行匹配
    match_genes_and_transcripts(gene_csv_file, transcript_txt_file, output_csv_file)
    
    # 显示输出文件的前几行作为示例
    try:
        print("\n输出文件前5行示例：")
        with open(output_csv_file, 'r', encoding='utf-8') as f:
            for i, line in enumerate(f):
                if i < 6:  # 显示表头+前5行数据
                    print(line.strip())
                else:
                    break
    except:
        pass

if __name__ == "__main__":
    main()
