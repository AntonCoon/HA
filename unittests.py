import unittest
import networkx as nx
from src import AlignedDB
from src import AlignedDBPreprocessor
from Bio import SeqIO
from src import Util


# Build Aligned De Bruijn graph tests
class BuildRefTest(unittest.TestCase):

    def test_nodes_amount(self):
        self.path_to_ref = "test/ref.fa"

        db = AlignedDB.AlignedDB(
            ["test/read1.fq", "test/read2.fq"],
            self.path_to_ref,
            "fastq",
            k_mer_len=5
        )

        db.build_ref()

        self.ref = SeqIO.parse(self.path_to_ref, "fasta")
        self.ref = str(next(self.ref).seq)

        self.assertEqual(len(self.ref) - db.k + 1, db.number_of_nodes())

    def test_ref_reconstruction(self):
        self.path_to_ref = "test/ref.fa"

        db = AlignedDB.AlignedDB(
            ["test/read1.fq", "test/read2.fq"],
            self.path_to_ref,
            "fastq",
            k_mer_len=5
        )

        db.build_ref()

        self.ref = SeqIO.parse(self.path_to_ref, "fasta")
        self.ref = str(next(self.ref).seq)

        kmer = list(enumerate(Util.k_mer_pairs(self.ref, db.k)))
        graph_path = list(
            (Util.KMer(s1, i), Util.KMer(s2, i + 1)) for i, (s1, s2) in kmer
        )
        self.assertEqual(self.ref, Util.get_haplotype_by_path(db, graph_path))


class BuildTest(unittest.TestCase):

    def test_graph(self):
        self.path_to_ref = "test/minimal_test/ref.fa"

        db = AlignedDB.AlignedDB(
            ["test/minimal_test/reads.fa"],
            self.path_to_ref,
            "fasta",
            k_mer_len=51
        )

        db.build_ref()
        db.build()

        gt = set()
        with open("test/minimal_test/gt.txt", "r") as file:
            for line in file:
                [h, _] = line.strip().split()
                gt.add(h)

        res = set(db.find_haplotypes())
        self.assertTrue(gt.issubset(res))


# Util function tests
class SplitTest(unittest.TestCase):

    def test_split_graph(self):
        graph = nx.DiGraph()
        graph.add_edge(1, 2)
        graph.add_edge(2, 3)
        graph.add_edge(3, 5)
        graph.add_edge(2, 4)
        graph.add_edge(4, 5)
        graph.add_edge(4, 6)
        graph.add_edge(6, 7)
        graph.add_edge(7, 8)
        graph.add_edge(8, 4)
        graph.add_edge(8, 5)
        graph.add_edge(8, 2)
        graph.add_edge(1, 8)

        paths_decomposition = [
            [(1, 2)],
            [(2, 3), (3, 5)],
            [(2, 4)],
            [(4, 5)],
            [(4, 6), (6, 7), (7, 8)],
            [(8, 4)],
            [(8, 5)],
            [(8, 2)],
            [(1, 8)]]

        self.assertEqual(
            sorted(Util.split_graph_by_paths(graph)),
            sorted(paths_decomposition)
        )


# Preprocessor test
class PreprocessorTest(unittest.TestCase):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        path_to_ref = "test/minimal_test/ref.fa"
        db = AlignedDB.AlignedDB(
            ["test/minimal_test/reads.fa"],
            path_to_ref,
            "fasta",
            k_mer_len=51
        )
        db.build_ref()
        db.build()
        self.prep = AlignedDBPreprocessor.AlignedDBPreprocessor(db, .9)

    def test_normalization(self):

        self.prep.normalize_parallel()

        sums = []
        for edges in self.prep.aligned_db.baskets.values():
            sums.append(
                sum(
                    [self.prep.aligned_db.edges[e]["coverage"] for e in edges]
                )
            )
        for s in sums:
            self.assertAlmostEqual(s, 1, delta=10**-6)

    def test_mean(self):
        # Check that coverage in paths from decomposition is equal
        self.prep.mean_by_path_parallel()
        paths_decomposition = Util.split_graph_by_paths(self.prep.aligned_db)
        for path in paths_decomposition:
            if path:
                cov = self.prep.aligned_db.edges[path[0]]["coverage"]
                for e in path:
                    self.assertEqual(
                        cov, self.prep.aligned_db.edges[e]["coverage"]
                    )


if __name__ == "__main__":
    unittest.main()
