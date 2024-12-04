import re
from Bio.Seq import Seq
from Bio.SeqIO import parse

#  1.参考基因组的转换


def preprocess_reference(reference_path):
    """
    预处理参考基因组，生成两种版本：
    - 所有 C 转为 T 的参考序列
    - 所有 G 转为 A 的反向互补参考序列
    """
    with open(reference_path, 'r') as file:
        ref_seq = "".join(line.strip() for line in file if not line.startswith('>'))

    # 生成正向参考序列的 C->T 转换
    forward_ref = ref_seq.replace('C', 'T')

    # 生成反向互补序列的 G->A 转换
    reverse_ref = str(Seq(ref_seq).reverse_complement()).replace('G', 'A')

    return forward_ref, reverse_ref
# 2. 预处理参考序列


def preprocess_reads(reads_path):
    """
    读取输入的测序数据，生成正向和反向互补序列。
    """
    reads = []
    with open(reads_path, 'r') as file:
        for line in file:
            if not line.startswith('>'):
                seq = line.strip()
                reads.append((seq, str(Seq(seq).reverse_complement())))
    return reads

# 3.模糊对比


def align_read_to_reference(read, reference, allow_mismatch=True):
    """
    将一个 read 与参考序列进行比对，允许 C-T/G-A 模糊匹配。
    返回匹配位置和得分。
    """
    best_score = -1
    best_position = -1

    for i in range(len(reference) - len(read) + 1):
        # 提取比对窗口
        window = reference[i:i+len(read)]

        # 比对得分计算
        score = 0
        for r_base, w_base in zip(read, window):
            if r_base == w_base:
                score += 1
            elif allow_mismatch and ((r_base == 'C' and w_base == 'T') or (r_base == 'G' and w_base == 'A')):
                score += 1

        # 更新最佳匹配
        if score > best_score:
            best_score = score
            best_position = i

    return best_position, best_score

# 4.全局对比
def align_reads_to_reference(reads, forward_ref, reverse_ref):
    """
    将所有 reads 与参考序列进行比对。
    """
    results = []

    for read, rev_read in reads:
        # 与正向参考序列比对
        fwd_pos, fwd_score = align_read_to_reference(read, forward_ref)
        rev_pos, rev_score = align_read_to_reference(rev_read, reverse_ref)

        # 记录最佳比对结果
        if fwd_score >= rev_score:
            results.append((read, "forward", fwd_pos, fwd_score))
        else:
            results.append((read, "reverse", rev_pos, rev_score))

    return results

# 5.主程序及计算


def main(reference_path, reads_path, output_path):
    # 1. 预处理参考基因组
    forward_ref, reverse_ref = preprocess_reference(reference_path)

    # 2. 预处理测序数据
    reads = preprocess_reads(reads_path)

    # 3. 比对 reads 到参考序列
    results = align_reads_to_reference(reads, forward_ref, reverse_ref)

    # 4. 保存结果
    with open(output_path, 'w') as output_file:
        for read, strand, position, score in results:
            output_file.write(f"{read}\t{strand}\t{position}\t{score}\n")

    print("比对完成，结果已保存到:", output_path)


# 执行主程序
if __name__ == "__main__":
    reference_path = "reference.fa"  # 输入参考基因组文件
    reads_path = "reads.fa"          # 输入测序 reads 文件
    output_path = "alignment_results.txt"  # 输出结果文件

    main(reference_path, reads_path, output_path)



