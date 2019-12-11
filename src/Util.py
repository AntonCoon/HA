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
