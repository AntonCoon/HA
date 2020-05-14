from src import AlignedDB
from src import Util
from Bio import SeqIO
from math import log
from tqdm import tqdm


class AlignedDBPreprocessor(object):
    def __init__(
            self,
            aligned_db: AlignedDB.AlignedDB,
            probability: float):
        self.aligned_db = aligned_db
        self.probability = probability
        self.reference_size = len(self.aligned_db.ref)
        self.read_len = 0
        self.reads_amount = 0
        self.eriksson_threshold = None

    def normalize_parallel(self):
        for basket, edges in self.aligned_db.baskets.items():
            basket_norm = sum(
                [self.aligned_db.get_edge_data(*e)["coverage"] for e in edges]
            )
            for e in edges:
                self.aligned_db.edges[e]["coverage"] /= basket_norm
                # if e not in self.aligned_db.ref_edges:
                #     self.aligned_db.edges[e]["coverage"] /= basket_norm
                # else:
                #     self.aligned_db.edges[e]["coverage"] = 1
            basket_norm = sum(
                [self.aligned_db.get_edge_data(*e)["coverage"] for e in edges]
            )
            for e in edges:
                self.aligned_db.edges[e]["coverage"] /= basket_norm

    def mean_by_path_parallel(self):
        paths_decomposition = Util.split_graph_by_paths(self.aligned_db)
        for path in paths_decomposition:
            mean = sum(
                self.aligned_db.get_edge_data(*e)["coverage"] for e in path
            ) / len(path)
            for e in path:
                self.aligned_db.edges[e]["coverage"] = mean

    def __calculate_parameters(self):
        whole_len = 0
        for path in self.aligned_db.path_to_reads:
            file_with_reads = SeqIO.parse(path, self.aligned_db.format)
            for read in file_with_reads:
                whole_len += len(read)
                self.reads_amount += 1
        self.read_len = whole_len / self.reads_amount

    def eriksson_clear(self):
        self.__calculate_parameters()
        # Eriksson threshold on distinguishable haplotype
        n = self.reference_size
        p = self.probability
        L = self.read_len
        N = self.reads_amount
        self.eriksson_threshold = - n * log(1 - p ** (1 / n)) / (L * N)
        targeted_edges = set()
        for e, info in tqdm(self.aligned_db.edges.items()):
            if info["coverage"] < self.eriksson_threshold:
                targeted_edges.add(e)
        for e in targeted_edges:
            self.aligned_db.remove_edge(*e)
            basket_key = (e[0].pos, e[1].pos)
            self.aligned_db.baskets[basket_key].remove(e)
        self.normalize_parallel()
        self.mean_by_path_parallel()
        removable = []
        for v in self.aligned_db:
            if self.aligned_db.degree(v) == 0:
                removable.append(v)
        for v in removable:
            self.aligned_db.remove_node(v)
