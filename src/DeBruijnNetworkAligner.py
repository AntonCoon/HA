from src import DeBruijnBuildNetwork
from src import Util
from networkx import MultiDiGraph
from Bio import SeqIO
from bisect import bisect_left
from bisect import bisect_right
from collections import defaultdict
from tqdm import tqdm
from collections import Counter


class AlignedDBNetwork(MultiDiGraph):
    def __init__(self):
        super().__init__()

    def get_edge_substring(self, edge: tuple) -> str:
        return self.edges[edge]['contig']

    def get_edge_coverage(self, edge: tuple) -> str:
        return self.edges[edge]['coverage']


class NetworkAligner(object):
    def __init__(self, db_graph: DeBruijnBuildNetwork.DBGraph):
        self.db_graph = db_graph
        self.edge_alignment = dict()
        self.read_alignment = dict()
        self.aligned_db_graph = None
        self.reference = str()
        self.encoded_nodes = None
        self.__buckets = list()
        self.__positions = list()
        self.__edges_in_bucket = None

    @staticmethod
    def align(ref: str, contig: str) -> tuple:
        arg_max = 0
        max_value = 0
        for start in range(len(ref) - len(contig) + 1):
            score = sum([ri == ci for ri, ci in zip(ref[start:], contig)])
            max_value, arg_max = max([(max_value, arg_max), (score, start)])
        return arg_max, arg_max + len(contig)

    def align_db_graph(self) -> None:
        reference_path = self.db_graph.get_heaviest_path()
        reference = Util.get_haplotype_by_path(self.db_graph, reference_path)
        self.reference = reference
        title = 'edge alignment'
        for edge, info in tqdm(self.db_graph.edges.items(), desc=title):
            # align whole contig from edge but set position according to rules
            # in get_edge_substring function
            contig = info['contig']
            start, end = self.align(reference, contig)
            u, v, _ = edge
            in_degree = self.db_graph.in_degree(u)
            out_degree = self.db_graph.out_degree(v)
            k = self.db_graph.k_mer_len
            if in_degree == 0 and out_degree == 0:
                # |0--->0|
                pass
            elif out_degree == 0:
                # ...0--->0|
                start += k // 2
            elif in_degree == 0:
                # |0--->0...
                end -= (k // 2 + k % 2)
            else:
                # |...0--->0...|
                start += k // 2
                end -= (k // 2 + k % 2)
            self.edge_alignment[edge] = (start, end)

    def align_reads(self):
        coordinates = []
        for e, (start, end) in self.edge_alignment.items():
            coordinates.append(start)
            coordinates.append(end)
        positions = sorted(set(coordinates))
        self.__positions = positions
        self.__buckets = list(zip(positions[:-1], positions[1:]))
        self.read_alignment = {k: defaultdict(int) for k in self.__buckets}
        # split read by buckets
        file_with_reads = SeqIO.parse(
            self.db_graph.file_path,
            self.db_graph.format
        )
        for read in tqdm(file_with_reads, desc='reads alignment'):
            read = str(read.seq)
            start, end = self.align(self.reference, read)
            first_bucket_idx = bisect_right(positions, start) - 1
            last_bucket_idx = bisect_left(positions, end)
            for bucket_id in range(first_bucket_idx, last_bucket_idx):
                bucket = self.__buckets[bucket_id]
                start_substring = max(bucket[0], start) - start
                end_substring = min(bucket[1], end) - start
                sub_read = read[start_substring: end_substring]
                self.read_alignment[bucket][sub_read] += 1

    def split_db_graph(self):
        positions = self.__positions
        self.aligned_db_graph = AlignedDBNetwork()
        encoded_nodes = {kmer: str(n) for n, kmer in enumerate(self.db_graph)}
        self.encoded_nodes = encoded_nodes
        self.__edges_in_bucket = defaultdict(set)
        # create topologically same as De Bruijn graph
        for e in self.db_graph.edges:
            u, v, idx = e
            new_u, new_v = encoded_nodes[u], encoded_nodes[v]
            self.aligned_db_graph.add_edge(new_u, new_v, idx)
            new_e = (new_u, new_v, idx)
            contig = self.db_graph.get_edge_substring(e)
            self.aligned_db_graph.edges[new_e]['contig'] = contig

        # Split edges of created graph and divide by buckets
        for e, (start, end) in self.edge_alignment.items():
            first_bucket_idx = bisect_right(positions, start) - 1
            last_bucket_idx = bisect_left(positions, end)
            u, v, idx = e
            new_u, new_v = encoded_nodes[u], encoded_nodes[v]
            old_edge_in_new_graph = (new_u, new_v, idx)
            contig = self.aligned_db_graph.edges[(new_u, new_v, idx)]['contig']
            self.aligned_db_graph.remove_edge(*old_edge_in_new_graph)
            for bucket_id in range(first_bucket_idx, last_bucket_idx):
                bucket = self.__buckets[bucket_id]
                bucket_start, bucket_end = bucket
                update_new_u, update_new_v = new_u, new_v
                if bucket_start != start:
                    # we should split edge in new graph
                    update_new_u = new_u + '_' + str(bucket_start)
                if bucket_end != end:
                    # first vertex again
                    update_new_v = new_u + '_' + str(bucket_end)
                new_contig = contig[bucket_start - start: bucket_end - start]

                new_idx = self.aligned_db_graph.add_edge(
                    update_new_u, update_new_v
                )
                new_e = (update_new_u, update_new_v, new_idx)
                self.aligned_db_graph.edges[new_e]['contig'] = new_contig
                self.__edges_in_bucket[bucket].add(new_e)

    def unite_same_edges_in_buckets(self):
        bucket_by_edge = dict()
        for bucket, edges in self.__edges_in_bucket.items():
            for e in edges:
                bucket_by_edge[e] = bucket
        id_edge = {e: idx for idx, e in enumerate(self.aligned_db_graph.edges)}
        edge_by_id = {idx: e for e, idx in id_edge.items()}

        for bucket in self.__buckets:
            # find groups of edges in bucket with same sub-sting
            e_groups = defaultdict(set)
            for e in self.__edges_in_bucket[bucket]:
                e_groups[self.aligned_db_graph.get_edge_substring(e)].add(
                    id_edge[e]
                )
            # for each group join duplicated edges and update
            # bucket_by_edge, code_by_edge and edge_by_code
            for substring, edges_codes in e_groups.items():
                if len(edges_codes) < 2:
                    continue
                all_u, all_v = set(), set()
                for code in edges_codes:
                    e = edge_by_id[code]
                    u, v, idx = e
                    all_u.add(u)
                    all_v.add(v)
                united_u = '_' + '_'.join(all_u)
                united_v = '_' + '_'.join(all_v)

                # Work with first node (united_u)
                # Input edges
                u_in = []
                for u in all_u:
                    u_in.extend(self.aligned_db_graph.in_edges(u, keys=True))
                for u_in_edge in u_in:
                    e_id = id_edge[u_in_edge]
                    from_node, u, idx = u_in_edge
                    substring = self.aligned_db_graph.get_edge_substring(
                        u_in_edge)
                    idx = self.aligned_db_graph.add_edge(from_node, united_u)
                    new_edge = (from_node, united_u, idx)
                    self.aligned_db_graph.edges[new_edge]['contig'] = substring
                    # update all maps and buckets
                    bucket = bucket_by_edge[u_in_edge]
                    bucket_by_edge[new_edge] = bucket
                    self.__edges_in_bucket[bucket] -= {u_in_edge}
                    self.__edges_in_bucket[bucket].add(new_edge)
                    id_edge[new_edge] = e_id
                    edge_by_id[e_id] = new_edge

                # Output edges
                u_out = []
                for u in all_u:
                    u_out.extend(self.aligned_db_graph.out_edges(u, keys=True))
                for u_out_edge in u_out:
                    e_id = id_edge[u_out_edge]
                    u, to_node, idx = u_out_edge
                    substring = self.aligned_db_graph.get_edge_substring(
                        u_out_edge)
                    idx = self.aligned_db_graph.add_edge(united_u, to_node)
                    new_edge = (united_u, to_node, idx)
                    self.aligned_db_graph.edges[new_edge]['contig'] = substring
                    # update all maps and buckets
                    bucket = bucket_by_edge[u_out_edge]
                    bucket_by_edge[new_edge] = bucket
                    self.__edges_in_bucket[bucket] -= {u_out_edge}
                    self.__edges_in_bucket[bucket].add(new_edge)
                    id_edge[new_edge] = e_id
                    edge_by_id[e_id] = new_edge

                # Work with second node (united_v)
                # Input edges
                v_in = []
                for v in all_v:
                    v_in.extend(self.aligned_db_graph.in_edges(v, keys=True))
                for v_in_edge in v_in:
                    e_id = id_edge[v_in_edge]
                    from_node, v, idx = v_in_edge
                    if from_node in all_u:
                        continue
                    substring = self.aligned_db_graph.get_edge_substring(
                        v_in_edge)
                    idx = self.aligned_db_graph.add_edge(from_node, united_v)
                    new_edge = (from_node, united_v, idx)
                    self.aligned_db_graph.edges[new_edge]['contig'] = substring
                    # update all maps and buckets
                    bucket = bucket_by_edge[v_in_edge]
                    bucket_by_edge[new_edge] = bucket
                    self.__edges_in_bucket[bucket] -= {v_in_edge}
                    self.__edges_in_bucket[bucket].add(new_edge)
                    id_edge[new_edge] = e_id
                    edge_by_id[e_id] = new_edge

                # Output edges
                v_out = []
                for v in all_v:
                    v_out.extend(self.aligned_db_graph.out_edges(v, keys=True))
                for v_out_edge in v_out:
                    e_id = id_edge[v_out_edge]
                    v, to_node, idx = v_out_edge
                    if to_node in all_u:
                        continue
                    substring = self.aligned_db_graph.get_edge_substring(
                        v_out_edge)
                    idx = self.aligned_db_graph.add_edge(united_v, to_node)
                    new_edge = (united_v, to_node, idx)
                    self.aligned_db_graph.edges[new_edge]['contig'] = substring
                    # update all maps and buckets
                    bucket = bucket_by_edge[v_out_edge]
                    bucket_by_edge[new_edge] = bucket
                    self.__edges_in_bucket[bucket] -= {v_out_edge}
                    self.__edges_in_bucket[bucket].add(new_edge)
                    id_edge[new_edge] = e_id
                    edge_by_id[e_id] = new_edge

                for v in all_v:
                    self.aligned_db_graph.remove_node(v)
                for u in all_u:
                    self.aligned_db_graph.remove_node(u)
                edges_codes.pop()
                for e_id in edges_codes:
                    removed_edge = edge_by_id[e_id]
                    self.aligned_db_graph.remove_edge(*removed_edge)

    def calculate_coverage(self):
        for e in self.aligned_db_graph.edges:
            self.aligned_db_graph.edges[e]['coverage'] = 0
        for bucket, edges in self.__edges_in_bucket.items():
            for read, amount in self.read_alignment[bucket].items():
                edge_inqlude_read = set()
                for e in edges:
                    if e in self.aligned_db_graph.edges:
                        sub = self.aligned_db_graph.get_edge_substring(e)
                        start, end = self.align(sub, read)
                        if sub[start: end] == read:
                            edge_inqlude_read.add(e)
                    # needed improve
                    if len(edge_inqlude_read) > 1:
                        break
                # lol
                if len(edge_inqlude_read) > 1:
                    continue
                for e in edge_inqlude_read:
                    sub = self.aligned_db_graph.get_edge_substring(e)
                    coverage = amount * len(read) / len(sub)
                    self.aligned_db_graph.edges[e]['coverage'] += coverage

        for bucket, edges in self.__edges_in_bucket.items():
            normalize = 0
            for e in edges:
                if e in self.aligned_db_graph.edges:
                    normalize += self.aligned_db_graph.edges[e]['coverage']
            for e in edges:
                if e in self.aligned_db_graph.edges:
                    self.aligned_db_graph.edges[e]['coverage'] /= normalize

        # multi-edge index normalization
        updated_aligned_db_graph = AlignedDBNetwork()
        for e, info in self.aligned_db_graph.edges.items():
            u, v, key = e
            idx = updated_aligned_db_graph.add_edge(u, v)
            contig = info['contig']
            cover = info['coverage']
            updated_aligned_db_graph.edges[(u, v, idx)]['contig'] = contig
            updated_aligned_db_graph.edges[(u, v, idx)]['coverage'] = cover
        self.aligned_db_graph = updated_aligned_db_graph

    def get_aligned_db_graph(self) -> DeBruijnBuildNetwork.DBGraph:
        self.align_db_graph()
        self.align_reads()
        self.split_db_graph()
        self.calculate_coverage()

        return self.aligned_db_graph

