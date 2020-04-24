import editdistance
import subprocess
import tempfile
import shutil
import pysam
from os import path
from collections import defaultdict
from ortools.linear_solver import pywraplp
from networkx import MultiDiGraph


def get_haplotype_by_path(db_graph, edge_path: list) -> str:
    return ''.join([db_graph.get_edge_substring(e) for e in edge_path])


def get_in_edges(graph: MultiDiGraph, node) -> list:
    return graph.in_edges(node, keys=True)


def get_out_edges(graph: MultiDiGraph, node) -> list:
    return graph.out_edges(node, keys=True)


def flatten(nested_list: list) -> list:
    return [elem for sublist in nested_list for elem in sublist]


def k_mer_pairs(read, k):
    size = len(read)
    if size < k:
        return None, None
    for i in range(size - k):
        yield read[i:(i + k)], read[(i + 1):(i + k + 1)]


def earth_mover_distance(dist1: list, dist2: list) -> float:
    solver = pywraplp.Solver(
        'earth_mover_distance',
        pywraplp.Solver.GLOP_LINEAR_PROGRAMMING)

    variables = dict()

    dirt_leaving_constraints = defaultdict(lambda: 0)
    dirt_filling_constraints = defaultdict(lambda: 0)

    objective = solver.Objective()
    objective.SetMinimization()

    for x, dirt_at_x in dist1:
        for y, capacity_of_y in dist2:
            amount_to_move_x_y = solver.NumVar(0, solver.infinity(),
                                               'z_{%s, %s}' % (x, y))
            variables[(x, y)] = amount_to_move_x_y
            dirt_leaving_constraints[x] += amount_to_move_x_y
            dirt_filling_constraints[y] += amount_to_move_x_y
            objective.SetCoefficient(amount_to_move_x_y,
                                     editdistance.eval(x, y))

    dist1 = {k: v for k, v in dist1}
    dist2 = {k: v for k, v in dist2}
    for x, linear_combination in dirt_leaving_constraints.items():
        solver.Add(linear_combination == dist1[x])

    for y, linear_combination in dirt_filling_constraints.items():
        solver.Add(linear_combination == dist2[y])

    status = solver.Solve()
    if status not in [solver.OPTIMAL, solver.FEASIBLE]:
        raise Exception('Unable to find feasible solution')

    return objective.Value()


def get_normalize_pair_list(list_of_pair: list) -> list:
    total_count = sum([count for _, count in list_of_pair])

    return [(k, count / total_count) for k, count in list_of_pair]


def read_ground_truth(path_to_file: str) -> list:
    haplos_repr = []
    with open(path_to_file, 'r') as haplos:
        for line in haplos:
            haplo, count = line.strip().split()
            haplos_repr.append((haplo, float(count)))

    return get_normalize_pair_list(haplos_repr)


def read_vgflow(path_to_file: str) -> list:
    sh_frecs = []
    sh_hapls = []
    with open(path_to_file, 'r') as reconstructed:
        for line in reconstructed:
            line = line.strip()
            if line[0] == '>':
                sh_frecs.append(float(line.split('=')[-1]))
            else:
                sh_hapls.append(line)
    return get_normalize_pair_list(list(zip(sh_hapls, sh_frecs)))


class BWAContextManager(object):
    def __init__(self, path_to_reads: list, ref: str):
        self.ref = ref
        self.path_to_reads = path_to_reads
        self.tmp_dir = tempfile.mkdtemp(dir='.')
        self.sam_file = None

    def __run_bwa_mem(self):
        command = [
            'bwa',
            'mem',
            path.join(self.tmp_dir, 'ref.fasta'),
            *self.path_to_reads,
            '-o',
            path.join(self.tmp_dir, 'align.sam')
        ]
        subprocess.run(
            ['bwa', 'index', path.join(self.tmp_dir, 'ref.fasta')],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL
        )
        subprocess.run(
            command,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL
        )

    def __enter__(self):
        with open(path.join(self.tmp_dir, 'ref.fasta'), 'w') as ref_fie:
            ref_fie.write('>ref\n')
            ref_fie.write(self.ref)
        self.__run_bwa_mem()
        self.sam_file = pysam.AlignmentFile(
            path.join(self.tmp_dir, 'align.sam'))

        return self

    def __exit__(self, exc_type, exc_value, exc_traceback):
        shutil.rmtree(self.tmp_dir)


class KMer(object):
    def __init__(self, seq, position):
        self.seq = seq
        self.pos = position

    def __hash__(self):
        return hash((self.seq, str(self.pos)))

    def __eq__(self, other):
        return self.__hash__() == other.__hash__()
