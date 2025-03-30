from collections import defaultdict, deque
from typing import Dict, List, Tuple

class RepeatUnit:
    def __init__(self, start_idx: int, unit_size: int, repeats: int, is_rev: bool):
        self.start_idx = start_idx    # 在参考序列中的起始索引
        self.unit_size = unit_size    # 重复单元的长度（不含重复）
        self.repeats = repeats        # 重复次数
        self.is_rev = is_rev          # 是否为反向互补序列

    def __repr__(self):
        return f"RepeatUnit(start_idx={self.start_idx}, unit_size={self.unit_size}, repeats={self.repeats}, is_rev={self.is_rev})"

class GraphStructure:
    def __init__(self):
        self.adjacency: Dict[int, Dict[int, RepeatUnit]] = defaultdict(dict)

    def insert_edge(self, src: int, dst: int, unit: RepeatUnit):
        self.adjacency[src][dst] = unit

def reverse_complement(seq: str) -> str:
    comp = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join(comp[base] for base in reversed(seq))

def discover_edges(ref_seq: str, query_seq: str) -> List[Tuple[int, RepeatUnit]]:
    found_edges = []
    substring_map = defaultdict(list)
    ref_length = len(ref_seq)
    ref_rc = reverse_complement(ref_seq)

    # 构建子串映射表
    for sub_len in range(1, ref_length + 1):
        for pos in range(ref_length - sub_len + 1):
            sub_str = ref_seq[pos:pos+sub_len]
            rev_sub = ref_rc[ref_length - pos - sub_len : ref_length - pos]
            substring_map[sub_str].append((pos, sub_len, False))
            substring_map[rev_sub].append((pos, sub_len, True))

    # 查询序列处理
    qry_length = len(query_seq)
    for sub_len in range(1, qry_length + 1):
        for start in range(qry_length - sub_len + 1):
            segment = query_seq[start:start+sub_len]
            if segment not in substring_map:
                continue
            for ref_info in substring_map[segment]:
                ref_start, length_val, is_rev = ref_info
                count = 1
                found_edges.append((start, RepeatUnit(ref_start, length_val, count, is_rev)))
                nxt = start + sub_len
                while nxt <= qry_length - sub_len:
                    if query_seq[nxt:nxt+sub_len] == segment:
                        count += 1
                        found_edges.append((start, RepeatUnit(ref_start, length_val, count, is_rev)))
                    else:
                        break
                    nxt += sub_len
    return found_edges

def build_graph_structure(graph: GraphStructure, ref_seq: str, query_seq: str):
    edge_list = discover_edges(ref_seq, query_seq)
    for q_start, unit in edge_list:
        end_point = q_start + unit.unit_size * unit.repeats
        graph.insert_edge(q_start, end_point, unit)

def breadth_first_search(graph: GraphStructure, query_length: int) -> List[RepeatUnit]:
    path_units = []
    predecessor = [(-1, 0)] * (query_length + 5)
    state = [0] * (query_length + 5)
    queue = deque([(0, 0)])
    state[0] = 1

    while queue:
        curr_q, curr_ref = queue.popleft()
        if curr_q == query_length:
            break
        for nxt, unit in graph.adjacency.get(curr_q, {}).items():
            valid_flag = 0
            if curr_ref == unit.start_idx:
                valid_flag = 1
            elif curr_ref == unit.start_idx + unit.unit_size:
                valid_flag = 2
            if valid_flag and state[nxt] == 0:
                state[nxt] = 1
                if valid_flag == 1:
                    predecessor[nxt] = (curr_q, unit.repeats - 1)
                    queue.append((nxt, curr_ref + unit.unit_size))
                else:
                    predecessor[nxt] = (curr_q, unit.repeats)
                    queue.append((nxt, curr_ref))
        state[curr_q] = 2

    # 回溯路径
    cur_idx = query_length
    prev_idx = predecessor[cur_idx][0]
    while prev_idx != -1:
        edge_unit = graph.adjacency[prev_idx][cur_idx]
        edge_unit.repeats = predecessor[cur_idx][1]
        if edge_unit.repeats > 0:
            path_units.append(edge_unit)
        cur_idx = prev_idx
        prev_idx = predecessor[cur_idx][0]
    return path_units

def main():
    query_string = "CTGCAACGTTCGTGGTTCATGTTTGAGCGATAGGCCGAAACTAACCGTGCATGCAACGTTAGTGGATCATTGTGGAACTATAGACTCAAACTAAGCGAGCTTGCAACGTTAGTGGACCCTTTTTGAGCTATAGACGAAAACGGACCGAGGCTGCAAGGTTAGTGGATCATTTTTCAGTTTTAGACACAAACAAACCGAGCCATCAACGTTAGTCGATCATTTTTGTGCTATTGACCATATCTCAGCGAGCCTGCAACGTGAGTGGATCATTCTTGAGCTCTGGACCAAATCTAACCGTGCCAGCAACGCTAGTGGATAATTTTGTTGCTATAGACCAACACTAATCGAGACTGCCTCGTTAGTGCATCATTTTTGCGCCATAGACCATAGCTAAGCGAGCCTGCCTCGTTAGTGCATCATTTTTGCGCCATAGACCATAGCTAAGCGAGCCTGCCTCGTTAGTGCATCATTTTTGCGCCATAGACCATAGCTAAGCGAGCCTGCCTCGTTAGTGCATCATTTTTGCGCCATAGACCATAGCTAAGCGAGCCTGCCTCGTTAGTGCATCATTTTTGCGCCATAGACCATAGCTAAGCGAGCTAGACCAACACTAATCGAGACTGCCTCGTTAGTGCATCATTTTTGCGCCATAGACCATAGCTAAGCGAGCTAGACCAACACTAATCGAGACTGCCTCGTTAGTGCATCATTTTTGCGCCATAGACCATAGCTAAGCGAGCTAGACCAACACTAATCGAGACTGCCTCGTTAGTGCATCATTTTTGCGCCATAGACCATAGCTAAGCGAGCGCTCGCTTAGCTATGGTCTATGGCGCAAAAATGATGCACTAACGAGGCAGTCTCGATTAGTGTTGGTCTATAGCAACAAAATTATCCACTAGCGTTGCTGGCTCGCTTAGCTATGGTCTATGGCGCAAAAATGATGCACTAACGAGGCAGTCTCGATTAGTGTTGGTCTATAGCAACAAAATTATCCACTAGCGTTGCTGCTTACCATCGGACCTCCACGAATCTGAAAAGTTTTAATTTCCGAGCGATACTTACGACCGGACCTCCACGAATCAGAAAGGGTTCACTATCCGCTCGATACATACGATCGGACCTCCACGACTCTGTAAGGTTTCAAAATCCGCACGATAGTTACGACCGTACCTCTACGAATCTATAAGGTTTCAATTTCCGCTGGATCCTTACGATCGGACCTCCTCGAATCTGCAAGGTTTCAATATCCGCTCAATGGTTACGGACGGACCTCCACGCATCTTAAAGGTTAAAATAGGCGCTCGGTACTTACGATCGGACCTCTCCGAATCTCAAAGGTTTCAATATCCGCTTGATACTTACGATCGCAACACCACGGATCTGAAAGGTTTCAATATCCACTCTATA"
    ref_string = "CTGCAACGTTCGTGGTTCATGTTTGAGCGATAGGCCGAAACTAACCGTGCATGCAACGTTAGTGGATCATTGTGGAACTATAGACTCAAACTAAGCGAGCTTGCAACGTTAGTGGACCCTTTTTGAGCTATAGACGAAAACGGACCGAGGCTGCAAGGTTAGTGGATCATTTTTCAGTTTTAGACACAAACAAACCGAGCCATCAACGTTAGTCGATCATTTTTGTGCTATTGACCATATCTCAGCGAGCCTGCAACGTGAGTGGATCATTCTTGAGCTCTGGACCAAATCTAACCGTGCCAGCAACGCTAGTGGATAATTTTGTTGCTATAGACCAACACTAATCGAGACTGCCTCGTTAGTGCATCATTTTTGCGCCATAGACCATAGCTAAGCGAGCCTTACCATCGGACCTCCACGAATCTGAAAAGTTTTAATTTCCGAGCGATACTTACGACCGGACCTCCACGAATCAGAAAGGGTTCACTATCCGCTCGATACATACGATCGGACCTCCACGACTCTGTAAGGTTTCAAAATCCGCACGATAGTTACGACCGTACCTCTACGAATCTATAAGGTTTCAATTTCCGCTGGATCCTTACGATCGGACCTCCTCGAATCTGCAAGGTTTCAATATCCGCTCAATGGTTACGGACGGACCTCCACGCATCTTAAAGGTTAAAATAGGCGCTCGGTACTTACGATCGGACCTCTCCGAATCTCAAAGGTTTCAATATCCGCTTGATACTTACGATCGCAACACCACGGATCTGAAAGGTTTCAATATCCACTCTATA"

    graph_obj = GraphStructure()
    build_graph_structure(graph_obj, ref_string, query_string)
    path = breadth_first_search(graph_obj, len(query_string))
    
    print("检测结果：")
    print("长度\t次数\t参考位置\t方向\t序列")
    for unit in reversed(path):
        seq = ref_string[unit.start_idx:unit.start_idx+unit.unit_size]
        direction = "反向" if unit.is_rev else "正向"
        print(f"{unit.unit_size}\t{unit.repeats}\t{unit.start_idx + unit.unit_size}\t{direction}\t{seq}")

if __name__ == '__main__':
    main()