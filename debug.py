from src import DeBruijnBuildNetwork
from src import ILPInputPreprocessor
from src import ILPMinimizer
import matplotlib.pyplot as plt
import numpy as np
import networkx as nx
from collections import defaultdict
from itertools import product
from tqdm import tqdm
from src import DeBruijnNetworkAligner
from src import Util
from bisect import bisect_left
from bisect import bisect_right

# file_path = './example/test_60_30_10.fasta'
# file_path = './example/simulation_60_30_10.fasta'
file_path = './example/simulation_heaviest.fasta'
# file_path = './example/simulation_2e8_path_70_30.fasta'
# file_path = './example/simulation_2e8_path_60_30_10.fasta'
# file_path = './example/silly.fasta'

k_mer_len = 11
file_extension = 'fasta'
db_graph = DeBruijnBuildNetwork.DBGraph(
    file_path,
    file_extension,
    k_mer_len)
db_graph.build()
db_graph.compression()

path = db_graph.get_heaviest_path()
hap = Util.get_haplotype_by_path(db_graph, path)

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


aligner.align_reads()
aligner.split_db_graph()
aligner.unite_same_edges_in_buckets()
aligner.calculate_coverage()

preproc = ILPInputPreprocessor.DataPreprocessor(aligner.aligned_db_graph)
haps, _ = preproc.find_haplotypes()
print('haplotype amount', len(set(haps)))

minimizer = ILPMinimizer.ILPMinimizer(
    aligner.aligned_db_graph, preproc.haplotypes_edges)
hps_thr_e = minimizer.edges_haplotypes
h_ids = {h: idx for idx, h in enumerate(haps)}

# # ----Plot graph
# pos = nx.layout.kamada_kawai_layout(db_graph)
# plt.figure(figsize=(8, 6))
# nx.draw_networkx_edges(db_graph, pos, alpha=0.4)
# nx.draw_networkx_nodes(db_graph, pos, node_size=60)
# plt.axis('off')
# plt.show()
# # ----Plot graph
# pos = nx.layout.kamada_kawai_layout(aligner.aligned_db_graph)
# plt.figure(figsize=(8, 6))
# nx.draw_networkx_edges(aligner.aligned_db_graph, pos, alpha=0.4)
# nx.draw_networkx_nodes(aligner.aligned_db_graph, pos, node_size=60)
# plt.axis('off')
# plt.show()
# #
# ----Print equation
for e, v in hps_thr_e.items():
    # length = len(db_graph.edges[e]['contig'])
    coverage = aligner.aligned_db_graph.edges[e]['coverage']
    print(
        # length,
        "*" if len(v) == len(set(haps)) else "-",
        str(
            round(coverage, 5)
        ) + ' -',
        ' - '.join(['F_' + str(h_ids[h]) for h in v])
    )
#
# # # ----Test huge alpha
# # big_alpha = 3
# # minimizer.find_alpha(big_alpha)
# # big_val, freqs = minimizer.find_frequencies()
# # non_zero = sum(np.array(list(freqs.values())) != 0)
# # big_val -= big_alpha * non_zero
# # print(
# #     'huge alpha  = {}\nnonzero frequencies amount = {}\nphi = {}'.format(
# #         big_alpha,
# #         non_zero,
# #         big_val
# #     )
# # )
#
#
# ----Different alphas
lmbds = np.linspace(0, .85, 85)
freqs = []
targets = []
phis = []
reconstructed = []
for lmbd in tqdm(lmbds):
    minimizer.find_alpha(lmbd)
    val, freq = minimizer.find_frequencies()
    reconstructed.append([(k, v) for k, v in freq.items() if v > 0])
    freqs.append(np.array(list(freq.values())))
    targets.append(val)
    non_zero = sum(freqs[-1] != 0)
    phis.append(val - lmbd * non_zero)

# ----Plot nonzero frequencies
p = plt.plot(
    lmbds,
    np.array(
        [len([x_rez for x_rez in rez if x_rez > 0]) for rez in freqs]
    ),
    '-o'
)
plt.xlabel('importance of zeros')
plt.ylabel('nonzero amount')
plt.show()

# ----Plot target value
p_tv = plt.plot(
    lmbds,
    np.array(
        targets
    ),
    '-o'
)
# plt.plot(
#     lmbds,
#     lmbds + big_val,
#     '--r'
# )
plt.xlabel('importance of zeros')
plt.ylabel('objective function')
plt.show()

# ----Plot error value
p_phi = plt.plot(
    lmbds,
    np.array(
        phis
    ),
    '-o'
)
plt.xlabel('importance of zeros')
plt.ylabel('objective function')
plt.show()

gt = []
with open('./example/heaviest_GT.txt', 'r') as gt_file:
    for line in gt_file:
        h, p = line.strip().split()
        gt.append((h, float(p)))

emd = [Util.earth_mover_distance(ansver, gt) for ansver in reconstructed]
# ----Plot error value
p_phi = plt.plot(
    lmbds,
    emd,
    '-o'
)
plt.xlabel('importance of zeros')
plt.ylabel('EMD')
plt.show()


# # ----Test something
# example = nx.MultiDiGraph()
# example.add_edge(1, 2)
# example.edges[(1, 2, 0)]['w'] = 1
# example.add_edge(2, 3)
# example.edges[(2, 3, 0)]['w'] = .9
# example.add_edge(2, 3)
# example.edges[(2, 3, 1)]['w'] = .1
# example.add_edge(3, 4)
# example.edges[(3, 4, 0)]['w'] = 1
# example.add_edge(4, 5)
# example.edges[(4, 5, 0)]['w'] = .9
# example.add_edge(4, 5)
# example.edges[(4, 5, 1)]['w'] = .1
# example.add_edge(5, 6)
# example.edges[(5, 6, 0)]['w'] = 1
# example.add_edge(6, 7)
# example.edges[(6, 7, 0)]['w'] = .9
# example.add_edge(6, 7)
# example.edges[(6, 7, 1)]['w'] = .1
# example.add_edge(7, 8)
# example.edges[(7, 8, 0)]['w'] = 1
#
# print(example.out_edges(2, keys=True))
# print(example.in_edges(3, keys=True))
#
# paths = nx.all_simple_paths(example, 1, 8)
# for path in map(nx.utils.pairwise, paths):
#     print(list(path))

# ----Plot graph
# pos = nx.layout.kamada_kawai_layout(example)
# plt.figure(figsize=(8, 6))
# nx.draw_networkx_edges(example, pos, alpha=0.4)
# nx.draw_networkx_nodes(example, pos, node_size=60)
# plt.axis('off')
# plt.show()

# print(Util.get_in_edges(example, 5))
# example.remove_node(5)
