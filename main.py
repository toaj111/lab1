NEGATIVE_INFINITY = -10**9

def dna_complement(seq):
    mapping = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    result = ""
    for ch in seq:
        result += mapping[ch]
    return result[::-1]

def kmp_search(seq, pattern):
    positions = []
    n, m = len(seq), len(pattern)
    if m == 0:
        return positions
    lps = [0] * m
    j = 0
    for i in range(1, m):
        while j > 0 and pattern[i] != pattern[j]:
            j = lps[j-1]
        if pattern[i] == pattern[j]:
            j += 1
        lps[i] = j
    j = 0
    for i in range(n):
        while j > 0 and seq[i] != pattern[j]:
            j = lps[j-1]
        if text[i] == pattern[j]:
            j += 1
        if j == m:
            positions.append(i - m + 1)
            j = lps[j-1]
    return positions

def compute_hash_matrix(seq, base=25, mod=10000000007):
    L = len(seq)
    hash_mat = [[0 for _ in range(L)] for _ in range(L)]
    power = [1] * L
    for i in range(1, L):
        power[i] = (power[i-1] * base) % mod
    hash_mat[0][0] = ord(seq[0]) - ord('A') + 1
    for i in range(1, L):
        hash_mat[0][i] = (hash_mat[0][i-1] * base + (ord(seq[i]) - ord('A') + 1)) % mod
    for sub_len in range(1, L+1):
        for i in range(1, L - sub_len + 1):
            prev_val = hash_mat[i-1][i+sub_len-2]
            adjust = (power[sub_len-1] * (ord(seq[i-1]) - ord('A') + 1)) % mod
            diff = (prev_val - adjust + mod) % mod
            hash_mat[i][i+sub_len-1] = (diff * base + (ord(seq[i+sub_len-1]) - ord('A') + 1)) % mod
    return hash_mat

class AlignmentCell:
    def __init__(self, score, ref_index, streak):
        self.score = score
        self.ref_index = ref_index
        self.streak = streak

class RepeatBlock:
    def __init__(self, block, ref_start, block_len, occurrences, inverted, q_start, q_end):
        self.block = block
        self.ref_start = ref_start
        self.block_len = block_len
        self.occurrences = occurrences
        self.inverted = inverted
        self.q_start = q_start
        self.q_end = q_end

class RepeatAnalyzer:
    def __init__(self, ref_seq, query_seq, unit_size=1):
        self.ref_seq = ref_seq
        self.query_seq = query_seq
        self.unit = unit_size
        self.repeat_blocks = []  
        self.dp = []           
        self.best_dp = []      
        self.bt_matrix = []   
        self.bt_path = []      

    def init_analysis(self):
        self.ref_len = len(self.ref_seq) - self.unit + 1
        self.query_len = len(self.query_seq) - self.unit + 1
        self.ref_hash = compute_hash_matrix(self.ref_seq)
        self.query_hash = compute_hash_matrix(self.query_seq)
        self.ref_rev_hash = compute_hash_matrix(dna_complement(self.ref_seq))
        
        self.dp = [[AlignmentCell(NEGATIVE_INFINITY, -1, 0) for _ in range(self.ref_len)]
                   for _ in range(self.query_len)]
        self.dp[0][0] = AlignmentCell(0, 0, 1)
        self.best_dp = [AlignmentCell(NEGATIVE_INFINITY, -1, -1) for _ in range(self.query_len)]
        self.best_dp[0] = AlignmentCell(0, 0, -1)
        self.bt_matrix = [[NEGATIVE_INFINITY for _ in range(self.ref_len)] for _ in range(self.query_len)]
        self.bt_path = [AlignmentCell(0, -1, 0) for _ in range(self.query_len)]

    def exact_match(self, qi, rj):
        if qi + self.unit > len(self.query_seq) or rj + self.unit > len(self.ref_seq):
            return False
        return self.query_seq[qi:qi+self.unit] == self.ref_seq[rj:rj+self.unit]

    def inverse_match(self, qi, rj):
        idx1 = self.ref_len - rj - 1
        idx2 = self.ref_len - rj + self.unit - 2
        if idx1 < 0 or idx1 >= len(self.ref_rev_hash):
            return False
        if idx2 < 0 or idx2 >= len(self.ref_rev_hash[idx1]):
            return False
        if qi + self.unit > len(self.query_seq):
            return False
        return self.query_hash[qi][qi + self.unit - 1] == self.ref_rev_hash[idx1][idx2]

    def run_alignment(self):
        score_main = 50
        score_new = -10
        score_main_inv = 49
        score_new_inv = -10
        score_cont = 1

        for qi in range(1, self.query_len):
            for rj in range(self.ref_len):
                if self.exact_match(qi, rj):
                    new_val = self.best_dp[qi-1].score + score_main + score_new
                    self.dp[qi][rj] = AlignmentCell(new_val, self.best_dp[qi-1].ref_index, 1)
                    if rj - 1 >= 0 and self.exact_match(qi-1, rj-1):
                        candidate = self.dp[qi-1][rj-1].score + score_main + self.dp[qi-1][rj-1].streak * score_cont
                        if self.dp[qi][rj].score < candidate:
                            self.dp[qi][rj] = AlignmentCell(candidate, rj-1, self.dp[qi-1][rj-1].streak + 1)
                    if self.best_dp[qi].score < self.dp[qi][rj].score:
                        self.best_dp[qi] = AlignmentCell(self.dp[qi][rj].score, rj, self.dp[qi][rj].ref_index)
                elif self.inverse_match(qi, rj):
                    new_val = self.best_dp[qi-1].score + score_main_inv + score_new_inv
                    self.dp[qi][rj] = AlignmentCell(new_val, self.best_dp[qi-1].ref_index, 1)
                    if rj + 1 < self.ref_len and self.inverse_match(qi-1, rj+1):
                        candidate = self.dp[qi-1][rj+1].score + score_main_inv + self.dp[qi-1][rj+1].streak * score_cont
                        if self.dp[qi][rj].score < candidate:
                            self.dp[qi][rj] = AlignmentCell(candidate, rj+1, self.dp[qi-1][rj+1].streak + 1)
                    if self.best_dp[qi].score <= self.dp[qi][rj].score:
                        self.best_dp[qi] = AlignmentCell(self.dp[qi][rj].score, rj, self.dp[qi][rj].ref_index)

    def trace_back(self):
        qi = self.query_len - 1
        rj = self.ref_len - 1
        while qi >= 0 and rj >= 0:
            self.bt_matrix[qi][rj] = self.dp[qi][rj].score
            self.bt_path[qi] = AlignmentCell(self.dp[qi][rj].score, rj, self.dp[qi][rj].ref_index)
            rj = self.dp[qi][rj].ref_index
            qi -= 1

        q_begin = 0
        while q_begin < self.query_len:
            q_end = q_begin + 1
            while q_end < self.query_len:
                temp = self.bt_path[q_end].ref_index if self.bt_path[q_end].ref_index != -1 else 0
                if abs(self.dp[q_end][temp].streak) == 1:
                    break
                q_end += 1
            ref_position = self.bt_path[q_begin].ref_index
            if ref_position < 0:
                ref_position = 0
            if self.inverse_match(q_begin, ref_position):
                block_obj = RepeatBlock(
                    self.query_seq[q_begin : q_end],
                    ref_position,
                    q_end - q_begin,
                    1,
                    True,
                    q_begin,
                    q_end - 1
                )
            else:
                block_obj = RepeatBlock(
                    self.query_seq[q_begin : q_end],
                    ref_position + (q_end - q_begin),
                    q_end - q_begin,
                    1,
                    False,
                    q_begin,
                    q_end - 1
                )
            self.repeat_blocks.append(block_obj)
            q_begin = q_end

    def merge_blocks(self):
        marker = [[0 for _ in range(self.ref_len)] for _ in range(self.query_len)]
        for r in range(self.ref_len - 1, -1, -1):
            q = 0
            while q < self.query_len and self.bt_matrix[q][r] == NEGATIVE_INFINITY:
                q += 1
            if q < self.query_len:
                marker[q][r] = 1

        for block in self.repeat_blocks:
            for j in range(block.q_end, block.q_start - 1, -1):
                temp = self.best_dp[j].ref_index if self.best_dp[j].ref_index != -1 else 0
                if temp < self.ref_len and marker[j][temp] == 1:
                    block.block_len -= 1

        idx = len(self.repeat_blocks) - 1
        while idx >= 0:
            jdx = 0
            merged = False
            while jdx < idx:
                if (abs(self.repeat_blocks[idx].ref_start - self.repeat_blocks[jdx].ref_start) < 10 and
                    abs(self.repeat_blocks[idx].block_len - self.repeat_blocks[jdx].block_len) < 10 and
                    self.repeat_blocks[idx].inverted == self.repeat_blocks[jdx].inverted):
                    self.repeat_blocks[jdx].occurrences += 1
                    del self.repeat_blocks[idx]
                    merged = True
                    break
                jdx += 1
            if not merged:
                idx -= 1
        self.repeat_blocks = [blk for blk in self.repeat_blocks if blk.block_len >= 10]

    def output_results(self):
        header = "{:>5} | {:>8} | {:>4} | {:>5} | {:>5}".format("Idx", "Ref_Pos", "Len", "Count", "Inv")
        print(header)
        for i, blk in enumerate(self.repeat_blocks):
            inv_str = "True" if blk.inverted else "False"
            print("{:5d} | {:8d} | {:4d} | {:5d} | {:>5}".format(i+1, blk.ref_start, blk.block_len, blk.occurrences, inv_str))

    def analyze(self):
        self.init_analysis()
        self.run_alignment()
        self.trace_back()
        self.merge_blocks()
        self.output_results()


def main():
    try:
        with open("ref.txt", "r") as ref_file:
            ref_data = ref_file.readline().strip()
    except Exception as err:
        print("无法读取 ref.txt")
        return

    try:
        with open("query.txt", "r") as qry_file:
            qry_data = qry_file.readline().strip()
    except Exception as err:
        print("无法读取 query.txt")
        return

    analyzer = RepeatAnalyzer(ref_data, qry_data, unit_size=1)
    analyzer.analyze()

if __name__ == '__main__':
    main()

