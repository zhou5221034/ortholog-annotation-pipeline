# ortholog-annotation-pipeline

跨物种同源基因比对与功能注释流程。面向**非模式物种单细胞分析**的下游需求：把物种自身的基因（如柑橘 `Cau.`）映射到**模式物种**（如拟南芥 `AT`）的同源基因上，并汇总 GO / KEGG / eggNOG / SwissProt 多来源注释，最终产出一张可用于富集分析和结果解读的完整注释表。

## 流程总览

```
                     ┌──────────────────────── 蛋白序列准备 ────────────────────────┐
 genome.fa + gff ──► │ LongestGff.pl  →  mapid（基因→最长转录本）                    │
                     │ gffread -x     →  out.cds                                   │
                     │ cds2aa.pl      →  out.pep（CDS 翻译为蛋白）                   │
                     │ GetSeq.pl      →  out.longest.pep（只保留最长转录本蛋白）        │
                     └────────────────────────────┬───────────────────────────────┘
                                                  ▼
                            ┌──────────── 同源基因鉴定（RBH） ────────────┐
                            │ Best2Best.pl → runBest.sh → result         │
                            │ （双向最佳比对，L1 层级结果）                  │
                            └────────────────────┬───────────────────────┘
                                                 ▼
        ┌──────────────────── 功能注释 ────────────────────┐
        │ eggNOG-mapper 注释 → pre_cluster.py              │
        │   ├─ GO.universal.tsv                           │
        │   └─ KEGG.universal.tsv                         │
        │ KeggKingdom.pl / KeggFilt.pl（按物种界过滤 KEGG）   │
        └────────────────────┬────────────────────────────┘
                             ▼
        ┌──────────── 基因匹配与注释合并（step4 → step5 → step6）────────────┐
        │ step4: 标记基因 ←→ 转录本 ID                                      │
        │ step5: 基因 ←→ 模式物种（拟南芥）同源基因（OrthoFinder Orthogroups）  │
        │ step6: 合并 KO/Pathway + SwissProt + GO + eggNOG → MBOS_eggnog.txt │
        └──────────────────────────────────────────────────────────────────┘
```

## 目录结构

```
.
├── Getpep.sh                  # 蛋白序列提取主脚本（串联下面 4 个步骤）
├── LongestGff.pl              # GFF → mapid（每个基因的最长转录本及其长度）
├── cds2aa.pl                  # CDS 翻译为蛋白序列（含基因模型质量检查）
├── GetSeq.pl                  # 按 mapid 从蛋白 fasta 中提取最长转录本序列
├── Best2Best.pl               # 生成 Best2Best(RBH) 运行所需的 config 与 runBest.sh
├── pre_cluster.py             # eggNOG/自建注释 → GO.universal.tsv + KEGG.universal.tsv
├── match_step4.py             # 标记基因表 + 转录本比对结果 → 追加 transcript_id 列
├── match_step5.py             # Orthogroups.tsv → 追加模式物种同源基因列
├── match_step6.py             # 合并 KO/Pathway + SwissProt + GO + eggNOG → 最终注释表
├── demo/                      # 示例数据（柑橘标记基因的完整中间产物与最终结果）
├── eggNOG-mapper/             # eggNOG 注释相关的界别 KO 列表与过滤脚本
│   ├── Filt.sh                # KEGG 过滤入口（串联下面两个脚本）
│   ├── KeggKingdom.pl         # 按物种界（Animals/Plants/...）生成 KO 列表
│   ├── KeggFilt.pl            # 用界别 KO 列表过滤 KEGG.universal.tsv
│   ├── pre_cluster.py         # 与根目录同名脚本一致（流程内独立运行版本）
│   ├── {Animals,Archaea,Bacteria,Fungi,Plants,Protists}.ko.txt
│   └── organism               # KEGG 物种清单（物种缩写 → 界）
└── 跨物种同源基因比对工具方法.docx   # 方法说明文档
```

## 依赖环境

| 依赖 | 用途 | 安装 |
| --- | --- | --- |
| Perl 5 | `*.pl` 脚本 | 系统自带 |
| Python 3 | `*.py` 脚本（仅用标准库） | 系统自带 |
| [gffread](https://github.com/gpertea/gffread) | 从基因组提取 CDS 序列 | `conda install -c bioconda gffread` |
| [Best2Best](https://github.com/zhangjiahao93/Best2Best) | 双向最佳比对（RBH）鉴定同源基因 | 见项目主页 |
| [eggNOG-mapper](https://github.com/eggnogdb/eggnog-mapper) | 功能注释（GO / KEGG / COG 等） | `conda install -c bioconda eggnog-mapper` |
| [DIAMOND](https://github.com/bbuchfink/diamond) | 比对 SwissProt 数据库 | `conda install -c bioconda diamond` |
| [OrthoFinder](https://github.com/davidemms/OrthoFinder) | 生成 `Orthogroups.tsv`（step5 输入） | `conda install -c bioconda orthofinder` |

## 使用说明

### 步骤 1：提取最长转录本蛋白序列

```bash
sh Getpep.sh <genome.fa> <annotation.gff>
```

依次执行 `LongestGff.pl` → `gffread` → `cds2aa.pl` → `GetSeq.pl`，产出：

| 文件 | 说明 |
| --- | --- |
| `mapid` | 基因 ID → 最长转录本 ID → CDS 长度 |
| `out.cds` | 全部转录本的 CDS 序列 |
| `out.pep` | 全部转录本翻译的蛋白序列 |
| `out.longest.pep` | **仅最长转录本**的蛋白序列（后续同源比对的输入） |

### 步骤 2：同源基因鉴定（RBH）

```bash
perl Best2Best.pl <reference.fa> <query.fa> <output_dir> [coverage] [identity]
```

| 参数 | 说明 | 默认值 |
| --- | --- | --- |
| `reference.fa` | 参考物种（如拟南芥）蛋白序列 | 必填 |
| `query.fa` | 目标物种（如柑橘）蛋白序列 | 必填 |
| `output_dir` | 输出目录 | 必填 |
| `coverage` | 比对覆盖率阈值，(0,1] | `0.3` |
| `identity` | 一致性阈值，(0,100] | `30` |

脚本会重命名序列（追加 `_ref` / `_tar` 标记）、生成 `config` 与 `runBest.sh`，执行 `runBest.sh` 后得到双向最佳比对结果 `result`。

### 步骤 3：功能注释

```bash
# eggNOG-mapper 注释结果 → 标准化的 GO / KEGG 表
python pre_cluster.py -t eggnog5 -e <species.emapper.annotations> -s <kegg_species_abbr> -o <outdir>
# 或使用自建注释
python pre_cluster.py -t inhouse -g <go_file> -k <kegg_file> -s <kegg_species_abbr> -o <outdir>

# 按物种界过滤 KEGG 结果
sh eggNOG-mapper/Filt.sh
perl eggNOG-mapper/KeggFilt.pl KEGG.universal.tsv Plants
```

产出 `GO.universal.tsv`（`Query / GO / Description / Level`）与 `KEGG.universal.tsv`（`Query / KO / Pathway / Description`）。

### 步骤 4～6：基因匹配与注释合并

```bash
# step4：标记基因表 ←→ 转录本比对结果，追加 transcript_id
python match_step4.py <marker_genes.csv> <best2best_result.txt> <out_step4.csv>

# step5：追加模式物种同源基因列（基于 OrthoFinder 的 Orthogroups.tsv）
python match_step5.py <Orthogroups.tsv> <step4_out.csv> <intermediate.txt> <out_step5.txt>

# step6：合并全部注释，输出最终表
python match_step6.py <out_step5.txt> KEGG.universal.tsv Swissprot.tsv GO.universal.tsv <species.emapper.annotations>
```

step6 的最终输出为 `MBOS_eggnog.txt`，在标记基因表基础上依次追加：

`transcript_id`、`Arabidopsis_Genes`、`KO_ID`、`Pathway_ID`、`Pathway_Description`、`SwissProt_Accession`、`SwissProt_Full`、`GO_ID`、`GO_Description`、`GO_Level`，以及 eggNOG 的全部注释列（`seed_ortholog`、`evalue`、`COG_category`、`Preferred_name`、`PFAMs` 等 20 列）。

未匹配到的项统一以 `/` 填充，并在终端输出匹配率统计。

## demo 数据

`demo/` 提供了一份完整的柑橘（`Cau.`）标记基因注释示例，可直接串起 step4 → step6：

| 文件 | 行数 | 说明 |
| --- | --- | --- |
| `Total_markergenes_02.csv` | 1556 | 输入：单细胞标记基因表（`gene, p_val, avg_log2FC, pct.1, pct.2, p_val_adj, cluster`） |
| `step3_result.csv` | 12874 | Best2Best RBH 结果（`transcript:AT... ↔ Cau...` + 覆盖率/一致性/L1 层级） |
| `step4_result.csv` | 1556 | step4 输出：追加 `transcript_id` 列 |
| `step5_result.csv` | 1556 | step5 输出：追加 `Arabidopsis_Genes` 列 |
| `Swissprot.tsv` | 32282 | DIAMOND 比对 SwissProt 结果 |
| `MBOS_eggnog.txt` | 1556 | **最终输出**：全部注释合并后的完整注释表 |

## 注意事项

1. **脚本内含硬编码的绝对路径**，移植到自己的环境时需要修改：
   - `Best2Best.pl` → `/gpfs/users/liuzm/software/Best2Best/`（Best2Best 套件位置）
   - `pre_cluster.py` → `/gpfs/users/liuzm/software/SingleCell/data`（GO/KEGG 参考数据库目录，对应 `eggNOG-mapper/data/`）
   - `eggNOG-mapper/Filt.sh`、`KeggFilt.pl` → `/gpfs/users/liuzm/software/eggNOG-mapper/`
2. `match_step5.py` 中的物种判定（`Cau.` / `gene_AT` / `AT`）与 `match_step6.py` 的输出文件名（`MBOS_eggnog.txt`）按当前项目写死，换物种时需相应调整。
3. `cds2aa.pl` 为第三方脚本（作者 Fan Wei，2006，BGI），此处仅作流程依赖收录，版权归原作者所有。
4. `pre_cluster.py` 中 `-s/--species` 需传入 KEGG 物种缩写（见 `eggNOG-mapper/organism`），用于把通路 ID 从 `map` 前缀替换为物种前缀。

## 本仓库未包含的内容

以下内容因体积或版权原因**未纳入版本控制**（`.gitignore` 已排除），需要自行准备：

| 未包含 | 体积 | 说明 |
| --- | --- | --- |
| `eggNOG-mapper/data/` | 102 MB | GO / KEGG 参考数据库（`go_base`、`k_id2map`、`kegg_pathway`、`all_kegg_organs`、`kegg_data/`），第三方数据 |
| `eggnog-mapper-master/` | 48 GB | eggNOG-mapper 第三方项目本体及其数据库、预编译二进制 |
| `uniprot_sprot.dmnd` | 280 MB | DIAMOND 格式的 SwissProt 数据库，可从 UniProt 重新构建 |
| `out.cds` / `out.pep` / `out.longest.pep` | 74 MB | 步骤 1 生成的中间文件，可由 `Getpep.sh` 重新生成 |
| `mapid` | 1.4 MB | `LongestGff.pl` 生成的基因-转录本对应表 |
