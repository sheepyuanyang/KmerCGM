# KmerCGM：基于 K-mer 的候选基因挖掘

> KmerCGM（K-mer-based Candidate Gene Mining）是一个面向短读长测序数据的候选基因挖掘流程。它通过“先做 K-mer 频率减法、后做从头组装”的策略，绕开单一线性参考基因组带来的比对偏差（reference bias），直接从常规 Illumina RNA-seq 数据中恢复高度分化、存在结构变异或参考基因组中缺失的候选基因序列，无需三代长读长测序、目标基因家族富集或完整参考基因组。

**适用场景**：目标基因已通过遗传定位 / BSA 锁定到一定染色体区间，但常规“比对—找差异”流程因参考偏差、大片段插入缺失（Indel）、存在/缺失变异（PAV）或序列高度分化而无法找回候选序列时，尤其适合使用本流程。

---

## 目录

- [1. 方法原理](#1-方法原理)
- [2. 工作流程](#2-工作流程)
- [3. 仓库结构](#3-仓库结构)
- [4. 依赖软件](#4-依赖软件)
- [5. 快速开始](#5-快速开始)
- [6. 示例：YrNAM 复现](#6-示例yrnam-复现)
- [7. 输出文件说明](#7-输出文件说明)
- [8. 主要参数](#8-主要参数)
- [9. 注意事项](#9-注意事项)
- [10. 引用与联系](#10-引用与联系)

## 1. 方法原理

传统图位克隆依赖短读长序列向单一线性参考基因组的高质量比对。当目标基因与参考基因组差异过大（如 NLR 类抗病基因的快速进化、大片段插入缺失、染色体重排），短读长会比对失败（unmapped）或错误比对到旁系同源位点，导致候选基因丢失。

KmerCGM 的核心思路：

1. 先**放宽条件回收**与目标区间、同源区间、未锚定 scaffold（ChrUn）以及完全未比对相关的 reads，避免候选信息过早丢失；
2. 将 reads 切成固定长度 K-mer，直接对“携带目标特征的样本（Sample-T）”和“对照样本（Sample-W）”做**精确的频率减法**，得到特征特异 K-mer；
3. 在组装之前完成去污染和泛基因组交叉过滤，最后才用**从头组装**重建候选转录本。

这种“先减法、后组装”的架构把比对误差尽可能排除在组装之前，同时完整保留了非参考等位基因信息。

## 2. 工作流程

KmerCGM 共 6 步（见论文 Figure 1）：

| 步骤 | 内容 | 主要工具 |
| --- | --- | --- |
| Step 1 | 目标关联 reads 回收：保留目标区间内、同源区间内、比对到 ChrUn、未比对 reads 及存在错配的非目标 reads，丢弃“完美比对到非目标区域”的 reads | samtools |
| Step 2 | 特征特异 K-mer 提取：Sample-T 与 Sample-W 共有 K-mer 被减去，保留 Sample-T 特异且深度达到阈值的 K-mer（默认 51-mer） | Jellyfish |
| Step 3 | 非宿主污染去除：与多物种植物转录组/CDS 库比对，剔除不属于植物的 K-mer | BLAST+ |
| Step 4 | 小麦泛基因组交叉过滤：将候选 K-mer 逐一比对到多个小麦参考基因组，按“未锚定噪音规则”和“同源得分规则”剔除离靶噪音 | BLAST+ |
| Step 5 | 目标 reads 回收与从头组装：提取含 Step 4 保留 K-mer 的 reads，用 Trinity 组装成非冗余 contig，再二次泛基因组验证 | Trinity、BLAST+ |
| Step 6 | 候选筛选与人工核查：将 Sample-T / Sample-W 的 reads 重新比对到 contig，进行 SNP/Indel 检测，结合表达模式与基因家族注释人工筛选 | GATK4 HaplotypeCaller、IGV、NLR-Annotator 等 |

## 3. 仓库结构

```text
Github/
├── CGM.pl                                   # 主流程脚本（Step 3 去污染—Step 5 组装；Step 1/2 示例命令见代码注释）
├── scripts/                                 # 【待补充】Step 1、Step 2 及 Step 6 的完整可执行脚本（当前目录为空）
├── example_YrNAM/                           # YrNAM 案例完整示例（输入配置文件）
│   ├── File.list                            # 样本清单（BAM/FASTQ 路径 + 样本类型）
│   ├── target_region.bed                    # 目标区间及同源区间
│   ├── blastdb.list                         # 小麦泛基因组 blast 数据库 + 目标/同源区间坐标
│   ├── plant.trans.list                     # 植物转录组/CDS 去污染数据库
│   └── unmap_chr.bed                        # 未锚定 scaffold（ChrUnknown）区间
├── Aegilops_bicornis_TB01.mRNA.pos.csv      # 【待补充】Ae. bicornis TB01 mRNA/CDS 坐标表，用途见下文
└── ReadMe_CN.md                             # 本说明文件
```

## 4. 依赖软件

核心流程（CGM.pl）依赖：

- Perl 5（模块：`List::Util`）
- samtools（读取 BAM、区间提取）
- Jellyfish 2.2.6（K-mer 计数，K=51）
- BLAST+ 2.8.0（`blastn`）
- Trinity 2.6.6（从头组装）
- gzip / `gunzip`（FASTQ 解压）

上游预处理与 Step 6 分析（不在 CGM.pl 内）：

- HISAT2 2.1.0（RNA-seq 比对）
- GATK4 4.2.4.0（`MarkDuplicates` 去重复、`HaplotypeCaller` 变异检测）
- NLR-Annotator（NLR 基因家族注释）
- IGV（候选 contig 人工核查）

【待补充：经过完整测试的操作系统 / Linux 发行版、容器镜像或 conda 环境文件；建议的运行内存与核数（示例 YrNAM 中 Trinity 使用 300 G 内存）。】

## 5. 快速开始

### 5.1 准备 4 个输入文件

CGM.pl 要求以下 4 个文件**与脚本放在同一目录**，且均为 Tab 分隔的纯文本：

**① File.list —— 样本清单**

```text
#Type  Run      BAM_file                  FQ1                       FQ2
W      SRRxxx   /path/to/W.sort.bam       /path/to/W_1.fastq.gz     /path/to/W_2.fastq.gz
M1     SRRyyy   /path/to/M1.sort.bam      /path/to/M1_1.fastq.gz    /path/to/M1_2.fastq.gz
```

- 第 1 列 Type：`W` = 野生型 / 对照（**只能有一个**）；`M1`、`M2`… = 携带目标基因 / 目标特征的样本（可多个）。
- 第 2 列 Run：样本编号 / 测序 run 名；同一 Type 允许对应多个 Run（同一批样本多次测序）。
- 第 3 列 BAM：比对后**已排序去重**的 BAM（建议 HISAT2 比对 + GATK4 MarkDuplicates 生成）。
- 第 4/5 列：双端 FASTQ（gzip 压缩）。

**② target_region.bed —— 目标区间及同源区间**

```text
Chr1B  1509827  5592386
Chr1A  1337814  4834683
Chr1D  174680   2945825
```

行数决定运行模式：**1 行 = 二倍体模式；>1 行 = 多倍体模式**（六倍体小麦示例为 3 行：目标染色体 1B + 两个同源染色体 1A/1D）。

**③ blastdb.list —— 小麦泛基因组 blast 数据库及目标/同源区间**

```text
#Genome  Blastdb  TargetChromosome  Start  End  HomeologousChromosome1  Start1  End1  HomeologousChromosome2  Start2  End2
CSv2.1   /path/to/iwgsc_refseqv2.1_assembly  Chr1B  1509827  5592386  Chr1A  1337814  4834683  Chr1D  174680  2945825
```

- 每行一个参考基因组（泛基因组库）；Step 4 将候选 K-mer 逐一比对到这些库进行交叉过滤。
- 多倍体模式需要 11 列；二倍体模式只使用前 5 列（Genome、Blastdb、TargetChromosome、Start、End）。
- 第一行 `#` 注释必须保留。

**④ plant.trans.list —— 植物转录组/CDS 去污染数据库**

```text
CS         /path/to/CS_v2.1_cds
Emmer      /path/to/Emmer.cds
ae_bicornis /path/to/ae.bicornis.TB01.annotation_gene.fa
...
```

每行两列（名称 + blast 数据库路径）。Step 3 按文件中的顺序逐库比对，剔除与任何植物转录本都不匹配的 K-mer（示例中含小麦、野生二粒小麦、节节麦、黑麦、Ae. bicornis、玉米、拟南芥、番茄、草莓、水稻、大麦等 11 个库）。

**（可选）unmap_chr.bed —— 未锚定 scaffold 区间**

存在该文件时，Step 1 会额外保留比对到这些区间（如 `ChrUnknown`）的 reads。

### 5.2 运行

```bash
cd /path/to/repo            # 确保 4 个输入文件与 CGM.pl 在同一目录
perl CGM.pl
```

所有输出写入 `KmerCGM_result/` 目录。

> **注意**：当前版本 CGM.pl 中 Step 1（reads 提取）与 Step 2（Jellyfish 计数）以注释形式保留示例命令，实际运行时需要先用对应脚本生成中间文件（`KmerCGM_result/<Type>.step1.1.fq` / `KmerCGM_result/<Type>.step1.2.fq` 和 `KmerCGM_result/Kmer.step3.<Type>.target.fa`）。【待补充：scripts/ 目录中 Step 1、Step 2、Step 3 提取及 Step 6 的完整脚本；补充后请更新本节为完整一键流程】

## 6. 示例：YrNAM 复现

**案例背景**：小麦条锈病抗性基因 YrNAM，定位在 1B 染色体 1.5–5.6 Mb（标记 Xpsp3000–Xsdauw79 之间）。该基因序列与参考基因组差异极大，原始研究中需要三代全长转录组才能克隆；KmerCGM 仅用 Illumina RNA-seq 即可找回。

**样本构成**（`example_YrNAM/File.list`）：

- `W`（对照）：SRR21456679；
- 6 个 EMS 突变体：M168（SRR21456680）、M160（SRR21456681）、M110（SRR21456682）、M38（SRR21456684）、M19（SRR21456685）、M16（SRR21456686）。

**运行方式**：将 `example_YrNAM/` 下 5 个配置文件复制到 CGM.pl 所在目录（或直接在该目录运行），先按本机实际路径修改 BAM/FASTQ 与 blast 数据库路径（当前为作者服务器绝对路径），然后执行：

```bash
perl CGM.pl
```

**预期结果**：【待补充：各突变体 Step 3/4/5 的中间数量统计、最终候选 contig 分组结果文件、以及“6 个突变体收敛到 YrNAM”的示例输出，便于用户核对复现是否成功。】

## 7. 输出文件说明

所有输出在 `KmerCGM_result/` 目录下：

| 文件 | 说明 |
| --- | --- |
| `<Type>.step1.1.fq` / `<Type>.step1.2.fq` | Step 1 保留的 R1/R2 reads |
| `<Type>.step2.jf` / `<Type>.step2.txt` | Jellyfish 计数结果（K-mer + 频数） |
| `Kmer.step3.<Type>.target.fa` / `.fix.fa` | 特征特异 K-mer 及去污染后的 K-mer |
| `Kmer.<Genome>.step4.<Type>.blast.txt` | 各泛基因组 blast 结果（outfmt 6） |
| `Kmer.step4.<Type>.target.fix.select.fa` | Step 4 筛选后保留的 K-mer |
| `<Type>.step5.1.fq` / `<Type>.step5.2.fq` | 含目标 K-mer 的 reads（供组装） |
| `trinity_OUT_<Type>/trinity.fasta` | Trinity 组装结果（候选转录本 contig） |
| `<pre>.<db>.fa` / `<pre>.<db>.blast.txt` | 去污染中间文件 |

## 8. 主要参数

| 参数 | 默认值 | 说明 |
| --- | --- | --- |
| 读长 `$readslen` | 150 | CGM.pl 顶部硬编码，需按数据实际读长修改 |
| K-mer 长度 | 51 | `jellyfish count -m 51` |
| Jellyfish | `-s 2G -C -t 2` | canonical 计数，内存 2 G，线程 2 |
| Step 2 最低深度 | ≥3 | 论文示例值；脚本注释中曾用 ≥5，请以 scripts/ 实际脚本为准 |
| Step 4 blastn | `-outfmt 6 -max_target_seqs 6 -num_threads 4` | 泛基因组过滤 |
| Step 3 去污染 blastn | `-word_size 21 -dust no -evalue 1e-5 -max_target_seqs 1 -num_threads 10` | 植物库去污染 |
| Trinity | `--max_memory 300G --CPU 4` | 示例内存，按机器调整 |

## 9. 注意事项

- **最适合低背景差异设计**：EMS 突变体与野生型背景几乎一致，目标变异在 K-mer 减法后极易突出；多突变体取交集可将候选收敛到极少数基因（如 YrNAM 案例收敛到 2 个候选组）。
- **高背景群体需结合生物学过滤**：F2/RIL/NIL 或亲本直接比较会产生大量背景多态。建议结合基因家族注释（如抗病 NLR 基因）、表达模式或更精细的遗传定位进一步筛选（Pm5e、Pm1a 案例即如此）。
- **计算资源**：Trinity 组装内存需求大（示例 300 G），建议在高性能计算节点运行。【待补充：各步骤实测运行时间与资源占用】
- **BAM 与 FASTQ 均需提供**：BAM 用于判断 reads 是否落在目标/同源/未比对区间，FASTQ 用于提取实际序列；两者样本一一对应。
- 本仓库路径（BAM、FASTQ、blast 数据库）均为作者服务器绝对路径，复现时需全部替换为本地路径。

## 10. 引用与联系

**论文**：

> Yuan Yang, Na Li, Xiaowu Liu, Haiyan Jia\*, Zhengqiang Ma\*. KmerCGM: K-mer-Based Candidate Gene Mining to Reduce Reference Bias in Short-Read Sequencing Data. Plant Methods.
>
> 【待补充：DOI、卷/期/页码及正式发表信息】

**仓库与版本**：

- GitHub 仓库：【待补充：仓库 URL】
- Zenodo / 存档 DOI：【待补充】
- License：【待补充】

**BibTeX 模板**（补全后即可使用）：

```bibtex
@article{kmercgm,
  title   = {KmerCGM: K-mer-Based Candidate Gene Mining to Reduce Reference Bias in Short-Read Sequencing Data},
  author  = {Yang, Yuan and Li, Na and Liu, Xiaowu and Jia, Haiyan and Ma, Zhengqiang},
  journal = {Plant Methods},
  year    = {2026},
  doi     = {},   % 待补充
  url     = {}    % 待补充
}
```

**致谢**：本研究受南京农业大学生物信息中心与高性能计算公共平台支持。

**联系方式**：zqm2@njau.edu.cn（马正强）；hyjia@njau.edu.cn（贾海燕）。
