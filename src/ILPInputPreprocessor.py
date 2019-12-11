from src import Util
from collections import defaultdict
from itertools import product
import networkx as nx


# noinspection PyCallingNonCallable
class DataPreprocessor(object):
    def __init__(self, db_graph):
        self.db_graph = db_graph
        self.haplotypes = list()
        self.haplotypes_edges = dict()

    @staticmethod
    def __multiply_paths(edge_repeats: defaultdict, paths: list) -> list:
        original_path = paths[0]
        edge_repeats_subset = dict()
        for edge in original_path:
            edge_repeats_subset[edge] = edge_repeats[edge]

        all_masks = []
        for edge in original_path:
            all_masks.append(list(range(edge_repeats_subset[edge])))

        all_masks = product(*all_masks)

        modified_paths = []
        for mask, path in zip(all_masks, paths):
            new_path = []
            for edge, idx in zip(path, mask):
                new_path.append((*edge, idx))
            modified_paths.append(new_path)
        return modified_paths

    def find_haplotypes(self) -> tuple:
        self.haplotypes = []
        srcs = []
        dsts = []
        for vertex in self.db_graph.nodes:
            indeg = self.db_graph.in_degree(vertex)
            outdeg = self.db_graph.out_degree(vertex)
            if indeg == 0:
                srcs.append(vertex)
            elif outdeg == 0:
                dsts.append(vertex)

        edges_amount = defaultdict(int)
        for e in self.db_graph.edges:
            edges_amount[e[:2]] += 1

        all_paths = defaultdict(list)
        for src, dst in product(srcs, dsts):
            for path in nx.all_simple_paths(self.db_graph, src, dst):
                edge_path = list(zip(path[:-1], path[1:]))
                all_paths[tuple(path)].append(edge_path)

        for _, paths in all_paths.items():
            paths = self.__multiply_paths(edges_amount, paths)
            for path in paths:
                haplotype = Util.get_haplotype_by_path(self.db_graph, path)
                self.haplotypes.append(haplotype)
                self.haplotypes_edges[haplotype] = path
        return self.haplotypes, self.haplotypes_edges
