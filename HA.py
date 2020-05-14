import argparse
from src import AlignedDB
from src import AlignedDBPreprocessor
from src import ILPInputPreprocessor
from src import ILPMinimizer
from src import Util


parser = argparse.ArgumentParser(
    description="Simple maximum parsimony Haplotype Assembler")

parser.add_argument(
    "--ref",
    help="path to reference genome",
    type=str
)

parser.add_argument(
    "--reads",
    nargs="+",
    help="path to origin read",
    type=str
)

parser.add_argument(
    "--format",
    help="reads file format",
    type=str,
    default="fasta"
)

parser.add_argument(
    "--k",
    help="length of k-mer for De Bruijn graph",
    default=61,
    type=int
)

parser.add_argument(
    "--o",
    help="path to my_output file",
    type=str,
    default="./haplotypes"
)


if __name__ == '__main__':
    args = vars(parser.parse_args())

    path_to_reference = args["ref"]
    paths_to_reads = args["reads"]
    file_extension = args["format"]
    k = args['k']
    out_path = args['o']

    # Build aligned De Bruijn graph by reference genome and reads
    aligned_db = AlignedDB.AlignedDB(
        path_to_reference=path_to_reference,
        path_to_reads=paths_to_reads,
        read_format=file_extension,
        k_mer_len=k
    )
    aligned_db.build()

    # Make normalization of aligned De Bruijn graph
    prep = AlignedDBPreprocessor.AlignedDBPreprocessor(aligned_db, 1 - 10**-4)
    prep.normalize_parallel()
    prep.mean_by_path_parallel()
    prep.eriksson_clear()
    print(prep.aligned_db.number_of_edges())
    aligned_db = prep.aligned_db

    # Find optimized system
    ilp_prep = ILPInputPreprocessor.DataPreprocessor(aligned_db)
    haplotypes, _ = ilp_prep.find_haplotypes()

    print(len(ilp_prep.haplotypes_edges))

    # Find optimal frequencies
    minimizer = ILPMinimizer.ILPMinimizer(
        aligned_db, ilp_prep.haplotypes_edges)
    print(prep.eriksson_threshold)
    minimizer.find_alpha(prep.eriksson_threshold)
    _, result = minimizer.find_frequencies()

    indexes = dict()
    for k, v in ilp_prep.haplotypes_edges.items():
        indexes[k] = (v[0][0].pos, v[-1][-1].pos + aligned_db.k)

    assembled = []
    for h, f in result.items():
        if f > prep.eriksson_threshold / 10:
            left, right = indexes[h][0], indexes[h][1]
            h = aligned_db.ref[:left] + h + aligned_db.ref[right:]
            assembled.append((h, f))

    result = Util.get_normalize_pair_list(assembled)

    with open('{}.fa'.format(out_path), 'w') as haplotypes_file:
        for idx, (haplotype, frequency) in enumerate(result):
            haplotypes_file.write(
                '>{} freq={}\n{}\n'.format(
                    str(idx),
                    frequency,
                    haplotype
                )
            )
