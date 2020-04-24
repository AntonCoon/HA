import unittest
from src import AlignedDB
from Bio import SeqIO
from src import Util
import pickle


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

    def test_nodes_amount(self):
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
        self.assertEqual(gt, res)


if __name__ == "__main__":
    unittest.main()
