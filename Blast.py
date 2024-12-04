from collections import defaultdict


def create_kmer_index(sequence, k):
    """ 创建k-mer索引，返回一个字典，其中键是k-mer，值是k-mer出现的所有位置 """
    kmer_index = defaultdict(list)
    for i in range(len(sequence) - k + 1):
        kmer = sequence[i:i + k]
        kmer_index[kmer].append(i)
    return kmer_index


def find_potential_matches(query, db_sequence, k):
    """ 在数据库序列中查找潜在的匹配序列 """
    kmer_index = create_kmer_index(db_sequence, k)
    potential_matches = []

    for i in range(len(query) - k + 1):
        query_kmer = query[i:i + k]
        if query_kmer in kmer_index:
            potential_matches.extend(kmer_index[query_kmer])

    return sorted(set(potential_matches))


def local_alignment(seq1, seq2, match=1, mismatch=-1, gap=-2):
    """ 简单的局部比对（Smith-Waterman算法） """
    len1, len2 = len(seq1), len(seq2)
    # 创建动态规划矩阵
    dp = [[0] * (len2 + 1) for _ in range(len1 + 1)]

    # 填充动态规划表
    for i in range(1, len1 + 1):
        for j in range(1, len2 + 1):
            match_score = match if seq1[i - 1] == seq2[j - 1] else mismatch
            dp[i][j] = max(0, dp[i - 1][j - 1] + match_score, dp[i - 1][j] + gap, dp[i][j - 1] + gap)

    # 找到最大值及其位置
    max_score = 0
    end_pos = (0, 0)
    for i in range(len1 + 1):
        for j in range(len2 + 1):
            if dp[i][j] > max_score:
                max_score = dp[i][j]
                end_pos = (i, j)

    # 回溯获得最佳比对
    aligned_seq1, aligned_seq2 = [], []
    i, j = end_pos
    while i > 0 and j > 0 and dp[i][j] > 0:
        if dp[i][j] == dp[i - 1][j - 1] + (match if seq1[i - 1] == seq2[j - 1] else mismatch):
            aligned_seq1.append(seq1[i - 1])
            aligned_seq2.append(seq2[j - 1])
            i -= 1
            j -= 1
        elif dp[i][j] == dp[i - 1][j] + gap:
            aligned_seq1.append(seq1[i - 1])
            aligned_seq2.append("-")
            i -= 1
        else:
            aligned_seq1.append("-")
            aligned_seq2.append(seq2[j - 1])
            j -= 1

    # 反转结果
    return ''.join(reversed(aligned_seq1)), ''.join(reversed(aligned_seq2)), max_score


def heuristic_blast(query, db_sequence, k=6):
    """ 启发式BLAST算法实现（种子匹配 + 局部比对） """
    potential_matches = find_potential_matches(query, db_sequence, k)
    alignments = []

    for start_pos in potential_matches:
        # 提取数据库序列的子串
        db_subseq = db_sequence[start_pos:start_pos + len(query)]
        # 局部比对
        aligned_query, aligned_db, score = local_alignment(query, db_subseq)
        if score > 0:
            alignments.append((aligned_query, aligned_db, score, start_pos))

    # 按照得分排序
    alignments.sort(key=lambda x: x[2], reverse=True)

    return alignments


# 示例数据
query_sequence = "GATTACA"
db_sequence = "CGATGATTACAGGATTAG"

# 启发式BLAST比对
alignments = heuristic_blast(query_sequence, db_sequence, k=4)

# 打印比对结果
for aligned_query, aligned_db, score, pos in alignments:
    print(f"Match at position {pos}:")
    print(f"Query:  {aligned_query}")
    print(f"DB:     {aligned_db}")
    print(f"Score:  {score}")
    print()
