from src import AlignedDB
from src import Util


class AlignedDBPreprocessor(object):
    def __init__(
            self,
            aligned_db: AlignedDB.AlignedDB,
            probability: float):
        self.aligned_db = aligned_db
        self.probability = probability

    def normalize_parallel(self):
        for basket, edges in self.aligned_db.baskets.items():
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

    def calculate_parameters(self):
        pass
