import argparse
from src import DeBruijnBuild
from src import DeBruijnPaths

parser = argparse.ArgumentParser(
    description='Simple maximum parsimony de novo Haplotype Assembler')

parser.add_argument(
    '--reads',
    help='path to origin read',
    type=str
)
parser.add_argument(
    '--format',
    help='reads file format',
    type=str,
    default='fasta'
)
parser.add_argument(
    '--k',
    help='length of k-mer for De Bruijn graph',
    default=60,
    type=int
)
parser.add_argument(
    '--o',
    help='path to output file',
    type=str,
    default='./haplotypes'
)


if __name__ == '__main__':
    args = vars(parser.parse_args())

    file_path = args['reads']
    file_extension = args['format']
    k_mer_len = args['k']
    out_path = args['o']

    db_graph = DeBruijnBuild.DBGraph(file_path, file_extension, k_mer_len)
    db_graph.compression()
    db_graph.simplify_bad_cover()
    db_graph.update_coverage()

    db_network = DeBruijnPaths.build_network(db_graph)

    all_haplotypes = DeBruijnPaths.find_all_haplos(db_graph, db_network)
    haplotypes_through_edge = DeBruijnPaths.get_haps_through_edge(
        db_graph, all_haplotypes)

    with open('{}.fs'.format(out_path), 'w') as haplotypes_file:
        for idx, hap_freq in enumerate(
                DeBruijnPaths.get_normalize_pair_list(
                    DeBruijnPaths.calc_frequencies(
                        haplotypes_through_edge,
                        db_graph,
                        db_network,
                        all_haplotypes
                    )
                )
        ):
            if hap_freq[1] > 0:
                haplotypes_file.write(
                    '>{} freq={}\n{}\n'.format(
                        str(idx),
                        hap_freq[1],
                        hap_freq[0]
                    )
                )
