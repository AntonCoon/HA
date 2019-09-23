import unittest
from src import DeBruijnBuild
# from src import DeBruijnPaths


def read_edges(path_to_file: str) -> dict:
    real_adj = {}
    with open(path_to_file, 'r') as edges_file:
        for line in edges_file:
            line = line.strip().split()
            real_adj[line[0]] = {k_mer for k_mer in line[1:]}

    return real_adj


def get_clear_adj(dbgraph: DeBruijnBuild.DBGraph) -> dict:
    assembled_adj = {}
    for src, dst in dbgraph.adj.items():
        dsts = set()
        for v in dst.keys():
            dsts.add(v)
        assembled_adj[src] = dsts

    return assembled_adj


# Build De Bruijn graph tests
class InitTest(unittest.TestCase):

    def test(self):
        de_bruijn = DeBruijnBuild.DBGraph('./test/test.fa', 'fasta', 10)
        real_adj = read_edges('./test/test_edges.txt')
        assembled_adj = get_clear_adj(de_bruijn)
        self.assertEqual(assembled_adj, real_adj)


class CoverTest(unittest.TestCase):

    def test(self):
        de_bruijn = DeBruijnBuild.DBGraph('./test/test.fa', 'fasta', 10)
        self.assertEqual(round(de_bruijn.get_mean_cover(), 1), 48.5)


class CompressionTest(unittest.TestCase):

    def test(self):

        de_bruijn = DeBruijnBuild.DBGraph('./test/test.fa', 'fasta', 10)
        de_bruijn.compression()

        real_adj = read_edges('./test/test_edges_compressed.txt')
        assembled_adj = get_clear_adj(de_bruijn)
        self.assertEqual(assembled_adj, real_adj)
