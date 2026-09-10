import sys
import os
import re

def main():
    # 检查输入参数（新增eggNOG文件作为第5个参数）
    if len(sys.argv) != 6:
        print("用法: python add_ko_go_pathway_v3.py <result_file.txt> <ko_pathway_file.txt> <swissprot_file.txt> <go_file.txt> <eggnog_file.txt>")
        print("示例: python add_ko_go_pathway_v3.py result_02.txt KEGG.universal.tsv Swissprot.tsv go_results.txt CitrusAurantium.emapper.annotations")
        sys.exit(1)
    
    result_file = sys.argv[1]
    ko_file = sys.argv[2]
    swissprot_file = sys.argv[3]
    go_file = sys.argv[4]
    eggnog_file = sys.argv[5]  # 新增eggNOG文件路径
    
    # 检查所有文件是否存在
    for file_path in [result_file, ko_file, swissprot_file, go_file, eggnog_file]:
        if not os.path.exists(file_path):
            print(f"错误: 文件不存在 - {file_path}")
            sys.exit(1)
    
    print("正在读取 KO/Pathway 信息...")
    gene_to_ko_info = {}
    with open(ko_file, 'r') as f:
        f.readline()  # 跳过表头
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split('\t')
            if len(parts) < 4:
                continue
            citrus_gene = parts[0].strip()
            ko_id = parts[1].strip()
            pathway_id = parts[2].strip()
            description = parts[3].strip()
            if citrus_gene not in gene_to_ko_info:
                gene_to_ko_info[citrus_gene] = []
            gene_to_ko_info[citrus_gene].append((ko_id, pathway_id, description))
    print(f"已读取 {len(gene_to_ko_info)} 个基因的 KO/Pathway 信息")
    
    print("正在读取 SwissProt 信息...")
    gene_to_swissprot = {}
    with open(swissprot_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split('\t')
            if len(parts) < 2:
                continue
            citrus_gene = parts[0].strip()
            swissprot_match = parts[1].strip()
            swissprot_accession = ""
            if swissprot_match.startswith('sp|'):
                match = re.search(r'sp\|([A-Za-z0-9]+)\|', swissprot_match)
                swissprot_accession = match.group(1) if match else swissprot_match
            elif swissprot_match.startswith('tr|'):
                match = re.search(r'tr\|([A-Za-z0-9]+)\|', swissprot_match)
                swissprot_accession = match.group(1) if match else swissprot_match
            else:
                swissprot_accession = swissprot_match
            gene_to_swissprot[citrus_gene] = {
                'full': swissprot_match,
                'accession': swissprot_accession
            }
    print(f"已读取 {len(gene_to_swissprot)} 个基因的 SwissProt 信息")
    
    print("正在读取 GO 信息...")
    gene_to_go_info = {}
    with open(go_file, 'r') as f:
        f.readline()  # 跳过表头（Query\tGO\tDescription\tLevel）
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split('\t')
            if len(parts) < 4:
                continue
            citrus_gene = parts[0].strip()
            go_id = parts[1].strip()
            go_desc = parts[2].strip()
            go_level = parts[3].strip()
            if citrus_gene not in gene_to_go_info:
                gene_to_go_info[citrus_gene] = []
            gene_to_go_info[citrus_gene].append((go_id, go_desc, go_level))
    print(f"已读取 {len(gene_to_go_info)} 个基因的 GO 信息")
    
    print("正在读取 eggNOG-mapper 注释信息...")
    gene_to_eggnog = {}
    with open(eggnog_file, 'r') as f:
        # 跳过注释行，找到真正的表头
        header = None
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('##'):
                continue  # 跳过注释行
            if line.startswith('#query'):
                header = line[1:].split('\t')  # 去掉开头的 #
                break
        
        if header is None:
            print("错误: 无法找到 eggNOG 文件的表头")
            sys.exit(1)
        
        # 获取各列索引
        col_indices = {col: idx for idx, col in enumerate(header)}
        required_cols = ['query', 'seed_ortholog', 'evalue', 'score', 'eggNOG_OGs', 
                        'max_annot_lvl', 'COG_category', 'Description', 'Preferred_name',
                        'GOs', 'EC', 'KEGG_ko', 'KEGG_Pathway', 'KEGG_Module', 
                        'KEGG_Reaction', 'KEGG_rclass', 'BRITE', 'KEGG_TC', 'CAZy', 
                        'BiGG_Reaction', 'PFAMs']
        
        for col in required_cols:
            if col not in col_indices:
                print(f"警告: eggNOG 文件缺少列 '{col}'")
        
        # 读取数据行
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split('\t')
            if len(parts) <= col_indices.get('query', 0):
                continue
            
            citrus_gene = parts[col_indices['query']].strip()
            gene_to_eggnog[citrus_gene] = {}
            
            # 提取所有需要的列
            for col in required_cols:
                if col in col_indices and col_indices[col] < len(parts):
                    gene_to_eggnog[citrus_gene][col] = parts[col_indices[col]].strip()
                else:
                    gene_to_eggnog[citrus_gene][col] = ''
    
    print(f"已读取 {len(gene_to_eggnog)} 个基因的 eggNOG 注释信息")
    
    # 输出文件
    output_file = "MBOS_eggnog.txt"
    with open(result_file, 'r') as f_in, open(output_file, 'w') as f_out:
        header = f_in.readline().strip()
        
        # 新增eggNOG相关列
        eggnog_cols = [
            'seed_ortholog', 'evalue', 'score', 'eggNOG_OGs', 'max_annot_lvl',
            'COG_category', 'Description', 'Preferred_name', 'eggNOG_GOs', 'EC',
            'KEGG_ko', 'KEGG_Pathway', 'KEGG_Module', 'KEGG_Reaction', 
            'KEGG_rclass', 'BRITE', 'KEGG_TC', 'CAZy', 'BiGG_Reaction', 'PFAMs'
        ]
        
        # 写入新的表头
        f_out.write(header + "\tKO_ID\tPathway_ID\tPathway_Description\tSwissProt_Accession\tSwissProt_Full\tGO_ID\tGO_Description\tGO_Level\t" + "\t".join(eggnog_cols) + "\n")
        
        for line in f_in:
            line = line.strip()
            if not line:
                continue
            parts = line.split('\t')
            if len(parts) < 1:
                continue
            citrus_gene = parts[0].strip()
            
            # 获取KO信息
            ko_info_list = gene_to_ko_info.get(citrus_gene, [])
            ko_str = pathway_str = desc_str = ""
            if ko_info_list:
                ko_ids, pathway_ids, descs = zip(*ko_info_list)
                ko_str = ";".join(ko_ids)
                pathway_str = ";".join(pathway_ids)
                desc_str = ";".join(descs)
            
            # 获取SwissProt信息
            swissprot_info = gene_to_swissprot.get(citrus_gene, {})
            swissprot_accession = swissprot_info.get('accession', '')
            swissprot_full = swissprot_info.get('full', '')
            
            # 获取GO信息
            go_info_list = gene_to_go_info.get(citrus_gene, [])
            go_id_str = go_desc_str = go_level_str = ""
            if go_info_list:
                go_ids, go_descs, go_levels = zip(*go_info_list)
                go_id_str = ";".join(go_ids)
                go_desc_str = ";".join(go_descs)
                go_level_str = ";".join(go_levels)
            
            # 获取eggNOG信息
            eggnog_data = gene_to_eggnog.get(citrus_gene, {})
            eggnog_values = [eggnog_data.get(col, '') for col in eggnog_cols]
            
            # 特殊处理：将eggNOG的GOs列重命名为eggNOG_GOs
            if 'GOs' in eggnog_data:
                eggnog_values[eggnog_cols.index('eggNOG_GOs')] = eggnog_data['GOs']
            
            # 写入所有信息
            f_out.write(f"{line}\t{ko_str}\t{pathway_str}\t{desc_str}\t{swissprot_accession}\t{swissprot_full}\t{go_id_str}\t{go_desc_str}\t{go_level_str}\t" + "\t".join(eggnog_values) + "\n")
    
    print(f"最终文件已保存: {output_file}")
    print("完成!")

if __name__ == "__main__":
    main()
