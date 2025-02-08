import os
import random
import re
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import Counter, defaultdict
from scipy.stats import fisher_exact
from tqdm import tqdm
from statsmodels.stats.multitest import multipletests
import numpy as np

# IUPAC nucleotide ambiguity codes (用于将明确的碱基集合映射回 IUPAC 碱基)
IUPAC_MAP = {
    frozenset(['A']): 'A',
    frozenset(['C']): 'C',
    frozenset(['G']): 'G',
    frozenset(['T']): 'T',
    frozenset(['A','G']): 'R',
    frozenset(['C','T']): 'Y',
    frozenset(['G','C']): 'S',
    frozenset(['A','T']): 'W',
    frozenset(['G','T']): 'K',
    frozenset(['A','C']): 'M',
    frozenset(['C','G','T']): 'B',
    frozenset(['A','G','T']): 'D',
    frozenset(['A','C','T']): 'H',
    frozenset(['A','C','G']): 'V',
    frozenset(['A','C','G','T']): 'N'
}

# 反向：简并碱基 -> 具体碱基集合
IUPAC_CODES = {
    'A': 'A', 'C': 'C', 'G': 'G', 'T': 'T',
    'R': 'AG', 'Y': 'CT', 'S': 'GC', 'W': 'AT',
    'K': 'GT', 'M': 'AC', 'B': 'CGT', 'D': 'AGT',
    'H': 'ACT', 'V': 'ACG', 'N': 'ACGT'
}

def parse_bed_file(bed_file_path, flank_size=7):
    locations = []
    with open(bed_file_path) as bed_file:
        for line in bed_file:
            if line.startswith('#') or line.strip() == '':
                continue
            parts = line.strip().split()
            chrom, start, end = parts[0], int(parts[1]), int(parts[2])
            strand = parts[3] if len(parts) > 3 else '+'
            start = max(0, start - flank_size)
            end += flank_size
            locations.append({
                'chrom': chrom,
                'start': start,
                'end': end,
                'strand': strand
            })
    return locations

def extract_sequences(genomic_locations, genome_fasta):
    sequences = []
    for loc in genomic_locations:
        chrom = loc['chrom']
        start = loc['start']
        end = loc['end']
        strand = loc['strand']
        if chrom not in genome_fasta:
            continue
        sequence = genome_fasta[chrom].seq[start:end]
        if strand == '-':
            sequence = sequence.reverse_complement()
        sequences.append(SeqRecord(sequence, id=f"{chrom}:{start}-{end}", description=""))
    return sequences

def generate_random_background(genome_fasta, num_sequences, seq_length):
    """
       Use only the longest chromosome as the source for background sequence extraction.
    """
    background = []
    valid_chroms = [chrom for chrom in genome_fasta if len(genome_fasta[chrom]) >= seq_length]
    if not valid_chroms:
        return background

    max_length = max(len(genome_fasta[chrom]) for chrom in valid_chroms)
    longest_chroms = [chrom for chrom in valid_chroms if len(genome_fasta[chrom]) == max_length]

    for _ in range(num_sequences):
        selected_chrom = random.choice(longest_chroms)
        chrom_seq = genome_fasta[selected_chrom].seq
        max_start = len(chrom_seq) - seq_length
        start = random.randint(0, max_start)
        end = start + seq_length
        background.append(
            SeqRecord(
                chrom_seq[start:end],
                id=f"{selected_chrom}:{start}-{end}",
                description=""
            )
        )
    return background

def motif_to_regex(motif):
    regex = []
    for base in motif.upper():
        bases = IUPAC_CODES.get(base, base)
        if len(bases) == 1:
            regex.append(bases)
        else:
            regex.append(f'[{bases}]')
    return ''.join(regex)

def scan_sequences_for_motifs(sequences, motifs):
    motif_regex = {m: re.compile(motif_to_regex(m)) for m in motifs}
    matches = []
    for seq_rec in sequences:
        seq_str = str(seq_rec.seq).upper()
        for motif, pattern in motif_regex.items():
            for match in pattern.finditer(seq_str):
                matches.append({
                    'motif': motif,
                    'position': match.start(),
                    'sequence': seq_str,
                    'location': seq_rec.id
                })
    return matches

def kmer_counting(sequences, k):
    """高效k-mer计数"""
    counts = defaultdict(int)
    for seq in sequences:
        seq_str = str(seq.seq).upper()
        for i in range(len(seq_str) - k + 1):
            kmer = seq_str[i:i + k]
            # 跳过含N的序列，以免统计不稳定
            if 'N' in kmer:
                continue
            counts[kmer] += 1
    return counts

def hierarchical_motif_discovery(reg_seqs, bg_seqs, min_len, max_len, top_n):
    """分层候选发现：从短到长逐层筛选并取top富集的k-mer"""
    candidates = []
    for k in range(min_len, max_len + 1):
        reg_counts = kmer_counting(reg_seqs, k)
        bg_counts = kmer_counting(bg_seqs, k)

        enrich_scores = {}
        # +1e-6 避免除零
        bg_denominator = (len(bg_seqs) + 1e-6)
        reg_denominator = (len(reg_seqs) + 1e-6)

        for kmer in reg_counts:
            reg_freq = reg_counts[kmer] / reg_denominator
            bg_freq = bg_counts.get(kmer, 0) / bg_denominator
            enrich = reg_freq / (bg_freq + 1e-9)  # 避免除0
            enrich_scores[kmer] = enrich

        # 排序取top_n
        sorted_kmers = sorted(enrich_scores.items(), key=lambda x: -x[1])
        top_kmers_this_length = [kmer for (kmer, _) in sorted_kmers[:top_n]]
        candidates += top_kmers_this_length

    # 去重
    return list(set(candidates))

def count_motif_occurrences(seq_rec, motif):
    """使用正则表达式高效计数"""
    seq = str(seq_rec.seq).upper()
    pattern = re.compile(motif_to_regex(motif))
    return len(pattern.findall(seq))

def precise_enrichment_analysis(reg_seqs, bg_seqs, candidates, alpha):
    """对候选k-mer做Fisher精确检验，返回显著者"""
    results = []
    for motif in tqdm(candidates, desc="Analyzing motifs"):
        k = len(motif)
        reg_obs = sum(count_motif_occurrences(seq, motif) for seq in reg_seqs)
        bg_obs = sum(count_motif_occurrences(seq, motif) for seq in bg_seqs)

        reg_total = sum(len(seq.seq) - k + 1 for seq in reg_seqs)
        bg_total = sum(len(seq.seq) - k + 1 for seq in bg_seqs)

        # Fisher检验
        table = [[reg_obs, bg_obs], [reg_total - reg_obs, bg_total - bg_obs]]
        _, p = fisher_exact(table, alternative='greater')
        results.append((motif, p, reg_obs, bg_obs, reg_total, bg_total))

    # 多重检验校正
    pvals = [x[1] for x in results]
    _, adj_pvals, _, _ = multipletests(pvals, alpha=alpha, method='fdr_bh')

    # 返回通过阈值的motif
    final = []
    for (motif, raw_p, reg_obs, bg_obs, reg_total, bg_total), adj_p in zip(results, adj_pvals):
        if adj_p <= alpha:
            final.append((motif, raw_p, adj_p, reg_obs, bg_obs, reg_total, bg_total))
    return final

def print_significant_results(results):
    """格式化输出结果"""
    print(f"\n{'Motif':<12}{'P-value':<12}{'Adj.P':<12}{'Reg/Obs':<12}{'Bg/Obs':<12}{'Fold':<8}")
    for res in sorted(results, key=lambda x: x[1]):
        motif, p, adj_p, reg_obs, bg_obs, reg_total, bg_total = res
        fold = ((reg_obs+1e-9)/(reg_total+1e-9)) / ((bg_obs+1e-9)/(bg_total+1e-9))  # 避免inf
        print(f"{motif:<12}{p:.2e}  {adj_p:.2e}  {reg_obs}/{reg_total:<10}  "
              f"{bg_obs}/{bg_total:<10}  {fold:.1f}x")

# =============== 新增：k-mer聚类与IUPAC共识 ===============

def hamming_distance(s1, s2):
    """仅适用于长度相同的k-mer"""
    if len(s1) != len(s2):
        return float('inf')
    return sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2))

def cluster_kmers(kmers, max_distance=1):
    """
    将同长度的k-mer做简单单链式聚类:
    若某k-mer与已有某簇中任何成员的海明距离≤max_distance，就并入该簇；否则另起新簇。
    返回一个列表，其中每个元素是k-mer的列表(一簇)。
    """
    clusters = []
    for kmer in kmers:
        added = False
        for clust in clusters:
            # 只要与此cluster中任意一个k-mer距离足够近，就并入
            if any(hamming_distance(kmer, member) <= max_distance for member in clust):
                clust.append(kmer)
                added = True
                break
        if not added:
            clusters.append([kmer])
    return clusters

def bases_to_iupac(bases_set):
    """将给定位置的碱基集合(例如{'A','T'})转换为对应的IUPAC码"""
    return IUPAC_MAP[frozenset(bases_set)]

def get_iupac_consensus_for_cluster(cluster):
    """
    对一簇长度相同的k-mer逐位统计，得到IUPAC共识序列。
    cluster是一个列表，如 ['GATC','GATT','GACC', ...]
    """
    if not cluster:
        return ""

    length = len(cluster[0])
    consensus = []
    for pos in range(length):
        # 收集该pos上所有k-mer的碱基
        bases_here = {k[pos] for k in cluster}
        # 转成IUPAC
        iupac_base = bases_to_iupac(bases_here)
        consensus.append(iupac_base)
    return "".join(consensus)

def cluster_significant_kmers(significant_motifs, max_distance=1):
    """
    将显著k-mer按长度分组, 各组内做聚类, 最后输出每簇的IUPAC共识
    :param significant_motifs: precise_enrichment_analysis 返回的结果列表
                               格式: [(motif, p, adj_p, reg_obs, bg_obs, reg_total, bg_total), ...]
    :param max_distance: 最大海明距离阈值
    :return: dict, key为k-mer长度, value为 [ (cluster成员列表, IUPAC共识), ... ]
    """
    # 仅取出 motif 字段
    motifs = [x[0] for x in significant_motifs]
    # 按k-mer长度分组
    length_groups = defaultdict(list)
    for m in motifs:
        length_groups[len(m)].append(m)

    consensus_dict = {}
    for k_len, kmers in length_groups.items():
        # 对同长度kmers做聚类
        c = cluster_kmers(kmers, max_distance=max_distance)
        # 生成共识
        clusters_with_consensus = []
        for clust in c:
            consensus_seq = get_iupac_consensus_for_cluster(clust)
            clusters_with_consensus.append((clust, consensus_seq))
        consensus_dict[k_len] = clusters_with_consensus
    return consensus_dict

# =============== 主流程 ===============
def motif_enrichment_pipeline(bed_file, fasta_file,
                              min_len=4, max_len=8,
                              top_per_length=10,
                              fdr_threshold=0.01,
                              num_bg=5000,
                              max_cluster_distance=1):
    """
    一个示例的主流程：
    1. 解析 bed+fasta，提取目标区序列 & 背景序列
    2. 层级扫描获取候选 k-mer
    3. Fisher 精确检验得到显著k-mer
    4. 对显著k-mer做聚类并生成IUPAC共识
    """
    # 1) 加载数据
    genome = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))
    reg_seqs = extract_sequences(parse_bed_file(bed_file), genome)
    if not reg_seqs:
        print("No valid regions found!")
        return
    bg_seqs = generate_random_background(genome, num_bg, len(reg_seqs[0].seq))

    # 2) 分层 k-mer 候选发现
    candidates = hierarchical_motif_discovery(reg_seqs, bg_seqs,
                                              min_len=min_len,
                                              max_len=max_len,
                                              top_n=top_per_length)

    # 3) 对候选做精确富集分析
    significant_motifs = precise_enrichment_analysis(reg_seqs, bg_seqs,
                                                     candidates,
                                                     alpha=fdr_threshold)

    # 打印显著结果
    print_significant_results(significant_motifs)

    # 4) 针对显著 k-mer 做聚类并生成IUPAC共识
    consensus_result = cluster_significant_kmers(significant_motifs,
                                                 max_distance=max_cluster_distance)

    # 输出聚类与共识结果
    print("\n=== Clustered IUPAC consensus (by k-mer length) ===")
    for k_len, cluster_info in consensus_result.items():
        print(f"\n--- k-mer length = {k_len} ---")
        for i, (cluster_list, consensus_seq) in enumerate(cluster_info, 1):
            print(f"Cluster {i}:")
            print(f"  Members: {cluster_list}")
            print(f"  IUPAC consensus: {consensus_seq}")

# =============== 直接运行的示例 ===============
if __name__ == "__main__":
    # 你可以根据实际文件名修改调用
    motif_enrichment_pipeline(
        bed_file="methyl_sites.bed",
        fasta_file="ecoli.fasta",
        min_len=4,
        max_len=8,
        top_per_length=10,
        fdr_threshold=0.01,
        num_bg=5000,
        max_cluster_distance=1  # 海明距离阈值，可根据需要调节
    )