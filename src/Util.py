from networkx import MultiDiGraph
import editdistance
from collections import defaultdict
from ortools.linear_solver import pywraplp


def get_haplotype_by_path(db_graph, path: list) -> str:
    haplotype = ''.join([db_graph.get_edge_substring(e) for e in path])
    return haplotype


def get_in_edges(graph: MultiDiGraph, node) -> list:
    return graph.in_edges(node, keys=True)


def get_out_edges(graph: MultiDiGraph, node) -> list:
    return graph.out_edges(node, keys=True)


def flatten(nested_list: list) -> list:
    return [elem for sublist in nested_list for elem in sublist]


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
