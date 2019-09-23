from Bio import SeqIO
from copy import deepcopy


def k_mer_pairs(read, k):
    size = len(read)
    if size < k:
        return None, None
    for i in range(size - k):
        yield read[i:(i + k)], read[(i + 1):(i + k + 1)]


class DBGraph(object):
    def __init__(self, file_path: str, fst_format: str, k_mer_len: int = 55):
        self.file = SeqIO.parse(file_path, fst_format)
        self.edges = {}
        self.adj = {}
        self.k_mer_len = k_mer_len
        self.read_len = 0

        read = ''
        for read in self.file:
            read = str(read.seq)
            for k_mer_1, k_mer_2 in k_mer_pairs(read, self.k_mer_len):
                self.add_edge(k_mer_1, k_mer_2)
        if read:
            self.read_len = len(read)

        self.inverse_adj = {v: set() for v in self.adj.keys()}
        not_visited = set(self.adj.keys())
        while not_visited:
            start = not_visited.pop()
            stack = [start]
            while stack:
                vert = stack.pop()
                not_visited -= {vert}
                neighbors = self.adj[vert]
                for neighbor in neighbors:
                    self.inverse_adj[neighbor].add(vert)
                    if neighbor in not_visited:
                        stack.append(neighbor)

    def get_mean_cover(self) -> float:
        sum_covers, amount = 0, 0
        for edge_id in self.edges.keys():
            coverage = self.edges[edge_id][1]
            sum_covers += len(self.edges[edge_id][0]) * coverage
            amount += 1
        return sum_covers / amount

    def update_coverage(self):
        for e_id in self.edges:
            coverage = self.edges[e_id][1]
            edge_length = len(self.edges[e_id][0])
            self.edges[e_id][1] = coverage / edge_length

    def add_edge(self, k_mer_1: str, k_mer_2: str) -> None:
        if k_mer_1 in self.adj:
            if k_mer_2 in self.adj[k_mer_1]:
                edge_id = self.adj[k_mer_1][k_mer_2]
                self.edges[edge_id][1] += 1
            else:
                edge_id = len(self.edges)
                self.adj[k_mer_1][k_mer_2] = edge_id
                self.edges[edge_id] = [k_mer_1 + k_mer_2[-1], 1]
                if k_mer_2 not in self.adj:
                    self.adj[k_mer_2] = dict()
        else:
            edge_id = len(self.edges)
            self.adj[k_mer_1] = {k_mer_2: edge_id}
            self.edges[edge_id] = [k_mer_1 + k_mer_2[-1], 1]
            if k_mer_2 not in self.adj:
                self.adj[k_mer_2] = dict()

    def delete_zero_deg(self) -> None:
        zero_deg = set()
        for vert in self.adj.keys():
            in_degree = len(self.adj[vert].keys())
            out_degree = len(self.inverse_adj[vert])
            if (in_degree == 0) and (out_degree == 0):
                zero_deg.add(vert)
        for vert in zero_deg:
            self.adj.pop(vert)
            self.inverse_adj.pop(vert)

    def compression(self) -> None:
        k_mer_len = self.k_mer_len
        target_vertices = set()
        for vert in self.adj.keys():
            in_degree = len(self.adj[vert])
            out_degree = len(self.inverse_adj[vert])
            if (in_degree == 1) and (out_degree == 1):
                target_vertices.add(vert)
        while target_vertices:
            vert = target_vertices.pop()
            prev_vert = list(self.inverse_adj[vert])[0]
            in_edge = self.adj[prev_vert][vert]
            neighbor = list(self.adj[vert].keys())[0]
            out_edge = self.adj[vert][neighbor]
            k_mer_1_cov = self.edges[in_edge][1]
            k_mer_2_cov = self.edges[out_edge][1]
            mean_covering = (k_mer_1_cov + k_mer_2_cov)
            self.edges[in_edge][1] = mean_covering
            self.edges[in_edge][0] += self.edges[out_edge][0][k_mer_len:]
            self.adj[prev_vert][neighbor] = self.adj[prev_vert].pop(vert)
            self.adj.pop(vert)
            self.inverse_adj[neighbor] -= {vert}
            self.inverse_adj[neighbor].add(prev_vert)
            self.inverse_adj.pop(vert)

    def simplify_tail(self):
        lower_bound = .5
        while True:
            deletions = 0
            mean_cover = self.get_mean_cover()
            adj_copy = deepcopy(self.adj)
            for from_vert in adj_copy.keys():
                for to_vert in adj_copy[from_vert].keys():
                    in_degree = len(self.adj[to_vert])
                    out_degree = len(self.inverse_adj[from_vert])
                    if (in_degree == 0) or (out_degree == 0):
                        edge_id = self.adj[from_vert][to_vert]
                        coverage: int = self.edges[edge_id][1]
                        current_cover = len(self.edges[edge_id][0]) * coverage
                        if current_cover < lower_bound * mean_cover:
                            self.adj[from_vert].pop(to_vert)
                            self.inverse_adj[to_vert] -= {from_vert}
                            self.edges.pop(edge_id)
                            deletions += 1

            if deletions == 0:
                break

            self.delete_zero_deg()
            self.compression()

    def simplify_bad_cover(self):
        lower_bound = 0
        while True:
            deletions = 0
            mean_cover = self.get_mean_cover()
            adj_copy = deepcopy(self.adj)
            for from_vert in adj_copy.keys():
                for to_vert in adj_copy[from_vert].keys():
                    edge_id = self.adj[from_vert][to_vert]
                    coverage: int = self.edges[edge_id][1]
                    current_cover = len(self.edges[edge_id][0]) * coverage
                    if current_cover < lower_bound * mean_cover:
                        self.adj[from_vert].pop(to_vert)
                        self.inverse_adj[to_vert] -= {from_vert}
                        self.edges.pop(edge_id)
                        deletions += 1

            if deletions == 0:
                break

            self.delete_zero_deg()
            self.compression()

    def get_true_length(self) -> dict:
        result = dict()
        for from_vert in self.adj.keys():
            for to_vert in self.adj[from_vert].keys():
                in_degree = len(self.adj[to_vert])
                out_degree = len(self.inverse_adj[from_vert])
                e_id = self.adj[from_vert][to_vert]
                if (in_degree == 0) or (out_degree == 0):
                    # true_len = len(self.edges[e_id][0]) - self.k_mer_len
                    true_len = len(self.edges[e_id][0])
                    # for uniform coverage
                    normalize_len = 1 / (1 - 1 / 2 * self.read_len / true_len)
                    result[e_id] = true_len / normalize_len
                else:
                    result[e_id] = len(self.edges[e_id][0])
        return result
