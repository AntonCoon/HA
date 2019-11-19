from src import DeBruijnBuildNetwork
from src import ILPInputPreprocessor
from src import ILPMinimizer
import matplotlib.pyplot as plt
import numpy as np
import networkx as nx
from collections import defaultdict
from itertools import product
from tqdm import tqdm

# file_path = './example/test_60_30_10.fasta'
# file_path = './example/simulation_60_30_10.fasta'
# file_path = './example/simulation_2e8_path_70_30.fasta'
file_path = './example/simulation_2e8_path_60_30_10.fasta'
# file_path = './example/silly.fasta'

k_mer_len = 9
file_extension = 'fasta'
db_graph = DeBruijnBuildNetwork.DBGraph(
    file_path,
    file_extension,
    k_mer_len)
db_graph.build()
# for e, info in db_graph.edges.items():
#     print(e, info)
db_graph.compression()
# for e, info in db_graph.edges.items():
#     print(e, info)

db_graph.normalize_coverage()

preproc = ILPInputPreprocessor.DataPreprocessor(db_graph)
haps, _ = preproc.find_haplotypes()

print('haplotype amount', len(set(haps)))

minimizer = ILPMinimizer.ILPMinimizer(db_graph, preproc.haplotypes_edges)


hps_thr_e = minimizer.edges_haplotypes
h_ids = {h: idx for idx, h in enumerate(haps)}

# ----Plot graph
pos = nx.layout.kamada_kawai_layout(db_graph)
plt.figure(figsize=(8, 6))
nx.draw_networkx_edges(db_graph, pos, alpha=0.4)
nx.draw_networkx_nodes(db_graph, pos, node_size=60)
plt.axis('off')
plt.show()

# ----Print equation
for e, v in hps_thr_e.items():
    # length = len(db_graph.edges[e]['contig'])
    coverage = db_graph.edges[e]['coverage']
    print(
        # length,
        "*" if len(v) == len(set(haps)) else "-",
        str(
            round(coverage, 5)
        ) + ' -',
        ' - '.join(['F_' + str(h_ids[h]) for h in v])
    )

# ----Test huge alpha
big_alpha = 3
minimizer.find_alpha(big_alpha)
big_val, freqs = minimizer.find_frequencies()
non_zero = sum(np.array(list(freqs.values())) != 0)
big_val -= big_alpha * non_zero
print(
    'huge alpha  = {}\nnonzero frequencies amount = {}\nphi = {}'.format(
        big_alpha,
        non_zero,
        big_val
    )
)
print(freqs.values())

minimizer.find_alpha(0)
_, freqs = minimizer.find_frequencies()
print(freqs.values())


# ----Different alphas
lmbds = np.linspace(0, 5, 20)
freqs = []
targets = []
phis = []
for lmbd in tqdm(lmbds):
    minimizer.find_alpha(lmbd)
    val, freq = minimizer.find_frequencies()
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
plt.plot(
    lmbds,
    lmbds + big_val,
    '--r'
)
plt.xlabel('importance of zeros')
plt.ylabel('objective function')
plt.show()

# ----Plot target value
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

# # ----Test something
# # example = nx.MultiDiGraph()
# # example.add_edge(1, 2)
# # example.add_edge(2, 3)
# # example.add_edge(2, 3)
# # example.add_edge(3, 4)
# # example.add_edge(4, 5)
# # example.add_edge(4, 5)
# # example.add_edge(5, 6)
# # example.add_edge(6, 7)
# # example.add_edge(6, 7)
# # example.add_edge(7, 8)
# #
# # all_paths = defaultdict(list)
# # for p in nx.all_simple_paths(example, 1, 8):
# #     edge_path = list(zip(p[:-1], p[1:]))
# #     all_paths[tuple(p)].append(edge_path[:])
# #
# #
# # edges_amount = defaultdict(int)
# # for e in example.edges:
# #     edges_amount[e[:2]] += 1
# #
# #
# # def multiply_paths(edge_repeats: defaultdict, paths: list) -> list:
# #     original_path = paths[0]
# #     edge_repeats_subset = dict()
# #     for edge in original_path:
# #         edge_repeats_subset[edge] = edge_repeats[edge]
# #
# #     all_masks = []
# #     for edge in original_path:
# #         all_masks.append(list(range(edge_repeats_subset[edge])))
# #
# #     all_masks = product(*all_masks)
# #
# #     modified_paths = []
# #     for mask, path in zip(all_masks, paths):
# #         new_path = []
# #         for edge, idx in zip(path, mask):
# #             new_path.append((*edge, idx))
# #         modified_paths.append(new_path)
# #     return modified_paths
#
#
# # k = list(all_paths.keys())[0]
# # for path in multiply_paths(edges_amount, all_paths[k]):
# #     print(path)
