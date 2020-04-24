from src import DeBruijnBuildNetwork
from src import ILPInputPreprocessor
from src import ILPMinimizer
from src import DeBruijnNetworkAligner
from src import Util
from os import path
from os import mkdir
import numpy as np


case_name = 'gattaca'
input_path = path.join('./input', case_name)
output_path = path.join('./my_output', case_name)
# mkdir(output_path)
file_path = path.join(input_path, 'reads.fastq')

handmade_alpha = 0.5
k_mer_len = 61

file_extension = 'fastq'
db_graph = DeBruijnBuildNetwork.DBGraph(
    file_path,
    file_extension,
    k_mer_len)
db_graph.build()
db_graph.compression()

heaviest_path = db_graph.get_heaviest_path()
hap = Util.get_haplotype_by_path(db_graph, heaviest_path)

# # ----Plot graph
# import networkx as nx
# from matplotlib import pyplot as plt
# pos = nx.layout.kamada_kawai_layout(db_graph)
# plt.figure(figsize=(8, 6))
# nx.draw_networkx_edges(db_graph, pos, alpha=0.4)
# nx.draw_networkx_nodes(db_graph, pos, node_size=60)
# plt.axis('off')
# plt.show()

aligner = DeBruijnNetworkAligner.NetworkAligner(db_graph)
aligner.align_db_graph()

print()
print(hap)
for e, (start, end) in aligner.edge_alignment.items():
    print(
        '_' * start + db_graph.get_edge_substring(e) + '_' * (len(hap) - end),
        start,
        end
    )


aligner.align_reads_with_bwa()
aligner.split_db_graph()
aligner.unite_same_edges_in_buckets()
aligner.calculate_coverage()

trial_preproc = ILPInputPreprocessor.DataPreprocessor(db_graph)
print('initial path amount')
print(len(trial_preproc.find_haplotypes()[0]))
print()

preproc = ILPInputPreprocessor.DataPreprocessor(aligner.aligned_db_graph)
haps, _ = preproc.find_haplotypes()
print('haplotype amount', len(set(haps)))

minimizer = ILPMinimizer.ILPMinimizer(
    aligner.aligned_db_graph, preproc.haplotypes_edges)
hps_thr_e = minimizer.edges_haplotypes
h_ids = {h: idx for idx, h in enumerate(haps)}

for e, v in hps_thr_e.items():
    coverage = aligner.aligned_db_graph.edges[e]['coverage']
    if len(v) != len(set(haps)):
        print(
            str(
                round(coverage, 5)
            ) + ' -',
            ' - '.join(['F_' + str(h_ids[h]) for h in v])
        )

minimizer.find_alpha(handmade_alpha)
big_val, freqs = minimizer.find_frequencies_square()
non_zero = sum(np.array(list(freqs.values())) != 0)
print(
    'huge alpha  = {}\nnonzero frequencies amount = {}\nphi = {}'.format(
        handmade_alpha,
        non_zero,
        big_val
    )
)

# with open(path.join(output_path, 'haplotypes' + '.txt'), 'w') as my_out:
for h, p in freqs.items():
    if p:
        print('{} {}\n'.format(h, p))
        # my_out.write('{} {}\n'.format(h, p))
