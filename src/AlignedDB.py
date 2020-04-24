from Bio import SeqIO
import networkx as nx
from src import Util
from itertools import product
from collections import defaultdict


class AlignedDB(nx.DiGraph):
    def __init__(
            self,
            path_to_reads: list,
            path_to_reference: str,
            read_format: str,
            k_mer_len: int = 61,
            **attr):
        super().__init__(**attr)
        self.path_to_reads = path_to_reads
        self.path_to_reference = path_to_reference
        self.format = read_format
        self.k = k_mer_len
        self.read_len = 0
        self.ref = SeqIO.parse(self.path_to_reference, "fasta")
        self.ref = str(next(self.ref).seq)
        self.baskets = dict()
        self.haplotypes = []

    def get_edge_substring(self, edge: tuple) -> str:
        u, v = edge
        in_degree = self.in_degree(u)
        out_degree = self.out_degree(v)
        k = self.k
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

    def build_ref(self):
        for i, (kmer1, kmer2) in enumerate(Util.k_mer_pairs(self.ref, self.k)):
            kmer1, kmer2 = Util.KMer(kmer1, i), Util.KMer(kmer2, i + 1)
            contig = kmer1.seq + kmer2.seq[-1]
            self.add_edge(kmer1, kmer2)
            self.edges[(kmer1, kmer2)]['coverage'] = 1
            self.edges[(kmer1, kmer2)]['contig'] = contig
            self.baskets[i] = {kmer1}
            self.baskets[i + 1] = {kmer2}

    def build(self):
        self.build_ref()
        aligned_reads = defaultdict(list)
        with Util.BWAContextManager(self.path_to_reads, self.ref) as bwa:
            for read_object in bwa.sam_file:
                read = read_object.seq
                start, end = read_object.pos, read_object.aend
                # take just whole aligned data
                if end - start != len(read):
                    continue
                aligned_reads[read].append(start)
                posed_kmers = [
                    (km1, km2, start + i) for i, (km1, km2) in enumerate(
                        Util.k_mer_pairs(read, self.k))
                ]
                new_edges = [
                    (Util.KMer(km1, i), Util.KMer(km2, i + 1))
                    for km1, km2, i in posed_kmers
                ]
                for u, v in new_edges:
                    contig = u.seq + v.seq[-1]
                    if (u, v) in self.edges:
                        self.edges[(u, v)]['coverage'] += 1
                        self.baskets[u.pos].add(u.seq)
                        self.baskets[v.pos].add(v.seq)
                    else:
                        self.add_edge(u, v)
                        self.edges[(u, v)]['coverage'] = 1
                        self.edges[(u, v)]['contig'] = contig
                        self.baskets[u.pos].add(u.seq)
                        self.baskets[v.pos].add(v.seq)

    def find_haplotypes(self) -> list:
        self.haplotypes = []
        srcs = []
        dsts = []
        for vertex in self.nodes:
            indeg = self.in_degree(vertex)
            outdeg = self.out_degree(vertex)
            if indeg == 0:
                srcs.append(vertex)
            elif outdeg == 0:
                dsts.append(vertex)

        for src, dst in product(srcs, dsts):
            for path in nx.all_simple_paths(self, src, dst):
                edge_path = list(zip(path[:-1], path[1:]))
                haplotype = Util.get_haplotype_by_path(self, edge_path)
                self.haplotypes.append(haplotype)
        return self.haplotypes