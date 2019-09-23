from src import DeBruijnBuild
from itertools import product
import networkx
import scipy
from copy import deepcopy
from scipy.optimize import minimize


def build_network(db_graph: DeBruijnBuild.DBGraph) -> networkx.DiGraph:
    result = networkx.DiGraph()
    for k_mer_1 in db_graph.adj:
        for k_mer_2 in db_graph.adj[k_mer_1]:
            e_id = db_graph.adj[k_mer_1][k_mer_2]
            result.add_edge(k_mer_1, k_mer_2, edge_id=e_id)
    return result


def get_haplotype_by_path(path: list,
                          db_graph: DeBruijnBuild.DBGraph) -> tuple:
    edge_ids_set = set()
    k_mer_len = db_graph.k_mer_len
    first_edge_id = db_graph.adj[path[0]][path[1]]
    edge_ids_set.add(first_edge_id)
    hap = db_graph.edges[first_edge_id][0]
    if len(path) < 3:
        return hap, edge_ids_set
    for u, v in zip(path[1:-1], path[2:]):
        e_id = db_graph.adj[u][v]
        edge_ids_set.add(e_id)
        hap += db_graph.edges[e_id][0][k_mer_len:]
    return hap, edge_ids_set


def find_all_haplos(db_graph: DeBruijnBuild.DBGraph,
                    db_network) -> list:
    result = list()
    all_src = []
    all_dst = []
    for u in db_network.nodes:
        in_deg = db_network.in_degree(u)
        out_deg = db_network.out_degree(u)
        if in_deg == 0:
            all_src.append(u)
        if out_deg == 0:
            all_dst.append(u)

    for u, v in product(all_src, all_dst):
        for path in networkx.all_simple_paths(db_network, u, v):
            result.append(
                get_haplotype_by_path(path, db_graph)
            )

    return result


def get_mean_bridges_cover(
        db_graph: DeBruijnBuild.DBGraph,
        db_network) -> float:
    db_cpy = deepcopy(db_network)
    all_src = []
    all_dst = []
    for u in db_network.nodes:
        in_deg = db_network.in_degree(u)
        out_deg = db_network.out_degree(u)
        if in_deg == 0:
            all_src.append(u)
        if out_deg == 0:
            all_dst.append(u)
    for u in all_src:
        db_cpy.add_edge('src', u)
    for v in all_dst:
        db_cpy.add_edge(v, 'dst')

    db_cpy = db_cpy.to_undirected()

    bridges = list(networkx.bridges(db_cpy))
    sum_bridge_cover, sum_lens, amount = 0, 0, 0
    for e in bridges:
        if e[1] in db_graph.adj[e[0]]:
            e_id = db_graph.adj[e[0]][e[1]]
        else:
            e_id = db_graph.adj[e[1]][e[0]]
        sum_bridge_cover += db_graph.edges[e_id][1]
        sum_lens += len(db_graph.edges[e_id][0])

    return sum_bridge_cover / sum_lens


def get_haps_through_edge(db_graph: DeBruijnBuild.DBGraph,
                          all_haps: list) -> dict:
    haps_through_edges = {e: set() for e in db_graph.edges}
    for idx, hap_and_edges in enumerate(all_haps):
        edge_set = hap_and_edges[1]
        for e in edge_set:
            haps_through_edges[e].add(idx)
    return haps_through_edges


def score_function_factory(haps_through_edges: dict,
                           db_graph: DeBruijnBuild.DBGraph,
                           db_network: networkx.DiGraph,
                           importance_of_zeros=None):

    if importance_of_zeros is None:
        importance_of_zeros = 100 * len(haps_through_edges)

    norm_lens = db_graph.get_true_length()
    b_cover = get_mean_bridges_cover(db_graph, db_network)

    def score(x: scipy.ndarray):
        result = 0
        for k, v in haps_through_edges.items():
            if v:
                edge_result = db_graph.edges[k][1] / norm_lens[k] / b_cover
                for hap_id in v:
                    edge_result -= x[hap_id]
                result += edge_result ** 2
        result += importance_of_zeros * scipy.sum(x != 0)
        return result

    return score


def calc_frequencies(
        haps_through_edges: dict,
        db_graph: DeBruijnBuild.DBGraph,
        db_network: networkx.DiGraph,
        all_haps: list) -> list:
    score = score_function_factory(
        haps_through_edges,
        db_graph,
        db_network
    )
    init_coef = db_graph.get_mean_cover() / 10
    x = (scipy.rand(len(all_haps)) * init_coef).astype(int)
    bounds = [(0, None) for _ in range(len(all_haps))]

    res = minimize(
        score,
        x,
        options={'maxiter': 500000},
        bounds=bounds
    )

    return [(hap[0], freq) for hap, freq in zip(all_haps, res.x.astype(int))]


def get_normalize_pair_list(list_of_pair: list) -> list:
    total_count = sum([count for _, count in list_of_pair])
    return [(k, count / total_count) for k, count in list_of_pair]
