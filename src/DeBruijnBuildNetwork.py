from Bio import SeqIO
import networkx as nx
from copy import deepcopy
from collections import Counter
import Util


# noinspection PyCallingNonCallable
class DBGraph(nx.MultiDiGraph):
    def __init__(
            self,
            file_path: str,
            read_format: str,
            k_mer_len: int = 55,
            **attr):
        super().__init__(**attr)
        self.file_path = file_path
        self.format = read_format
        self.k_mer_len = k_mer_len
        self.read_len = 0
        self.coded_nodes = {}

    def build(self):
        file_with_reads = SeqIO.parse(self.file_path, self.format)
        read = ''
        for read in file_with_reads:
            read = str(read.seq)
            if len(read) <= self.k_mer_len:
                continue
            if set(read.lower()) != set('atgc'):
                continue
            for k_mer_1, k_mer_2 in Util.k_mer_pairs(read, self.k_mer_len):
                contig = k_mer_1 + k_mer_2[-1]
                if (k_mer_1, k_mer_2) in self.edges:
                    self.edges[(k_mer_1, k_mer_2, 0)]['coverage'] += 1
                else:
                    self.add_edge(k_mer_1, k_mer_2)
                    self.edges[(k_mer_1, k_mer_2, 0)]['coverage'] = 1
                    self.edges[(k_mer_1, k_mer_2, 0)]['contig'] = contig
        for _, info in self.edges.items():
            info['coverage'] *= len(info['contig'])
        if read:
            self.read_len = len(read)

        self.coded_nodes = {v: idx for idx, v in enumerate(self.nodes)}

    def get_edge_substring(self, edge: tuple) -> str:
        u, v, _ = edge
        in_degree = self.in_degree(u)
        out_degree = self.out_degree(v)
        k = self.k_mer_len
        if in_degree == 0 and out_degree == 0:
            # |0--->0|
            return self.edges[edge]['contig']
        elif out_degree == 0:
            # ...0--->0|
            return self.edges[edge]['contig'][k // 2:]
        elif in_degree == 0:
            # |0--->0...
            return self.edges[edge]['contig'][:-k // 2]
        else:
            # |...0--->0...|
            return self.edges[edge]['contig'][k // 2: -k // 2]

    def get_mean_cover(self) -> float:
        sum_covers, amount = 0, 0
        for e, info in self.edges:
            coverage = info['coverage']
            sum_covers += len(info['contig']) * coverage
            amount += 1
        return sum_covers / amount

    def compression(self) -> None:
        k_mer_len = self.k_mer_len
        while True:
            is_simplified = False
            all_verts = list(self.adj.keys())
            for vertex in all_verts:
                out_degree = self.in_degree(vertex)
                in_degree = self.out_degree(vertex)
                if (in_degree == 1) and (out_degree == 1):
                    is_simplified = True
                    prev_vert = list(self.in_edges(vertex))[0][0]
                    next_vert = list(self.out_edges(vertex))[0][1]
                    cov_1 = self.edges[(prev_vert, vertex, 0)]['coverage']
                    cov_2 = self.edges[(vertex, next_vert, 0)]['coverage']
                    contig_1 = self.edges[(prev_vert, vertex, 0)]['contig']
                    contig_2 = self.edges[(vertex, next_vert, 0)]['contig']
                    edge_idx = self.add_edge(prev_vert, next_vert)
                    new_edge = self.edges[(prev_vert, next_vert, edge_idx)]
                    new_edge['coverage'] = cov_1 + cov_2
                    new_edge['contig'] = contig_1 + contig_2[k_mer_len:]
                    self.remove_node(vertex)
            if not is_simplified:
                break

    def delete_zero_deg(self) -> None:
        vertices = list(self.nodes)
        for vertex in vertices:
            if self.degree(vertex) == 0:
                self.remove_node(vertex)

    def simplify_tails(self) -> None:
        lbnd = .1
        while True:
            mean_cover = self.get_mean_cover()
            target_tails = []
            for e, info in self.edges.items():
                is_start = (self.in_degree(e[0]) + self.out_degree(e[0])) == 1
                is_end = (self.out_degree(e[1]) + self.in_degree(e[1])) == 1
                if is_start or is_end and info['coverage'] < lbnd * mean_cover:
                    target_tails.append(e)
            if not target_tails:
                break
            for u, v, key in target_tails:
                self.remove_edge(u, v, key)
            self.compression()

        self.delete_zero_deg()

    def simplify_bad_cover(self) -> None:
        lbnd = 0.1
        while True:
            mean_cover = self.get_mean_cover()
            target_edges = []
            for e, info in self.edges.items():
                if info['coverage'] < lbnd * mean_cover:
                    target_edges.append(e)
            if not target_edges:
                break
            for u, v, key in target_edges:
                self.remove_edge(u, v, key)
            self.compression()

        self.delete_zero_deg()

    def __get_length_coefficients(self) -> dict:
        result = dict()
        for edge, info in self.edges.items():
            from_vertex, to_vertex, _ = edge
            out_degree = self.out_degree(to_vertex)
            in_degree = self.in_degree(from_vertex)
            result[edge] = 1
            if (in_degree == 0) or (out_degree == 0):
                # for uniform coverage
                length = len(info['contig'])
                result[edge] -= 1 / 2 * (self.read_len - 1) / length
        return result

    def __get_expected_bridges_cover(self) -> float:
        db_cpy = deepcopy(self)
        all_src = []
        all_dst = []
        for u in self.nodes:
            in_deg = self.in_degree(u)
            out_deg = self.out_degree(u)
            if in_deg == 0:
                all_src.append(u)
            if out_deg == 0:
                all_dst.append(u)
        for u in all_src:
            db_cpy.add_edge('src', u)
        for v in all_dst:
            db_cpy.add_edge(v, 'dst')

        repeats_of_edges = Counter([e[:2] for e in db_cpy.edges])
        # multi edge is not bridge
        db_undirected = nx.Graph()
        for e in repeats_of_edges:
            db_undirected.add_edge(*e)
        bridges = list(nx.bridges(db_undirected))
        sum_bridge_cover, bridges_amount, amount = 0, 0, 0
        for e in bridges:
            if repeats_of_edges[e] > 1 or repeats_of_edges[e[::-1]] > 1:
                continue

            if 'dst' in e or 'src' in e:
                continue
            # bridge is not multi-edge, so key=0
            e = (*e, 0)
            if e not in self.edges:
                e = list(e)
                e[0], e[1] = e[1], e[0]
            e = tuple(e)
            sum_bridge_cover += self.edges[e]['coverage']
            bridges_amount += 1

        return sum_bridge_cover / bridges_amount if bridges_amount else 1

    def normalize_coverage(self) -> None:
        norm_length = self.__get_length_coefficients()

        for e, info in self.edges.items():
            info['coverage'] /= norm_length[e]

    def get_heaviest_path(self) -> list:
        path = []
        all_src = []
        for u in self.nodes:
            in_deg = self.in_degree(u)
            if in_deg == 0:
                all_src.append(u)

        all_start_edges = list()
        for src in all_src:
            for u in self.adj[src]:
                for idx, info in self.adj[src][u].items():
                    edge, cov = (src, u, idx), info['coverage']
                    all_start_edges.append((cov, edge))

        _, current_edge = max(all_start_edges)
        visited = set()
        while current_edge:
            path.append(current_edge)
            visited.add(current_edge)
            _, current_node, idx = current_edge
            possible_next_edges = list()
            for next_node in self.adj[current_node]:
                for idx, info in self.adj[current_node][next_node].items():
                    edge = (current_node, next_node, idx)
                    cov = info['coverage']
                    possible_next_edges.append((cov, edge))
            if possible_next_edges:
                _, current_edge = max(possible_next_edges)
            if current_edge in visited:
                break

        return path

    def print_encoded_adjacency_list(self):
        cod = self.coded_nodes
        for v, neighbours in self.adj.items():
            print(cod[v], ': ', end='')
            all_neighbours = []
            for u, items in neighbours.items():
                all_neighbours.extend([cod[u]] * len(items))
            print(
                '{' +
                ', '.join(map(str, [u for u in all_neighbours])) +
                '}'
            )
