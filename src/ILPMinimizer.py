from src import DeBruijnBuildNetwork
import gurobipy as gb
import numpy as np


class ILPMinimizer(object):
    def __init__(
            self,
            db_graph: DeBruijnBuildNetwork.DBGraph,
            haplotypes_edges: dict):
        self.db_graph = db_graph
        self.haplotypes_edges = haplotypes_edges
        self.edges_haplotypes = {e: [] for e in db_graph.edges}
        for haplotype, edges in haplotypes_edges.items():
            for e in edges:
                self.edges_haplotypes[e].append(haplotype)
        self.alpha = None

    def find_alpha(self, alpha):
        self.alpha = alpha

    def find_frequencies(self) -> tuple:
        alpha = self.alpha
        ord_haplotypes = list(self.haplotypes_edges.keys())
        haplotypes_id = {
            haplotype: idx for idx, haplotype in enumerate(ord_haplotypes)
        }

        m = gb.Model("mip1")
        u = list()
        for e in self.db_graph.edges:
            e_id = '_'.join(map(str, e))
            ui = m.addVar(
                vtype=gb.GRB.CONTINUOUS, name='u_{}'.format(str(e_id)))
            u.append((e, ui))

        f, b = list(), list()
        for idx, _ in enumerate(ord_haplotypes):
            fi = m.addVar(
                vtype=gb.GRB.CONTINUOUS, name='F_{}'.format(str(idx)))
            bi = m.addVar(
                vtype=gb.GRB.BINARY, name='B_{}'.format(str(idx)))
            f.append(fi)
            b.append(bi)

        # Set ILP objective function
        m.setObjective(
            sum([ui for _, ui in u]) + alpha * sum(b), gb.GRB.MINIMIZE
        )

        # Set ILP restrictions
        for e, ui in u:
            idxs = [haplotypes_id[hap] for hap in self.edges_haplotypes[e]]
            freqs = [f[idx] for idx in idxs]
            e_id = '_'.join(map(str, e))
            m.addConstr(
                ui >= self.db_graph.edges[e]['coverage'] - sum(freqs),
                name='c_pos_{}'.format(str(e_id)))
            m.addConstr(
                ui >= sum(freqs) - self.db_graph.edges[e]['coverage'],
                name='c_neg_{}'.format(str(e_id)))

        m.addConstr(sum(f) == 1, name='normalization')

        for fi, bi in zip(f, b):
            m.addConstr(bi >= fi)

        m.params.LogToConsole = False
        m.optimize()

        final_freq = np.array([fi.x for fi in f])

        rez = (
            m.objVal,
            {hap: final_freq[idx] for idx, hap in enumerate(ord_haplotypes)}
        )

        return rez
