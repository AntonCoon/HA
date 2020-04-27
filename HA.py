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
    aligned_db.build_ref()
    aligned_db.build()

    # Make normalization of aligned De Bruijn graph
    prep = AlignedDBPreprocessor.AlignedDBPreprocessor(aligned_db, 0.99)
    prep.normalize_parallel()
    prep.mean_by_path_parallel()
    prep.eriksson_clear()
    aligned_db = prep.aligned_db

    # Find optimized system
    ilp_prep = ILPInputPreprocessor.DataPreprocessor(aligned_db)
    haplotypes, _ = ilp_prep.find_haplotypes()

    # Find optimal frequencies
    minimizer = ILPMinimizer.ILPMinimizer(
        aligned_db, ilp_prep.haplotypes_edges)
    minimizer.find_alpha(prep.eriksson_threshold / 10)
    _, result = minimizer.find_frequencies()

    result = [(h, f) for h, f in result.items() if f > prep.eriksson_threshold]
    result = Util.get_normalize_pair_list(result)

    with open('{}.fs'.format(out_path), 'w') as haplotypes_file:
        for idx, (haplotype, frequency) in enumerate(result):
            haplotypes_file.write(
                '>{} freq={}\n{}\n'.format(
                    str(idx),
                    frequency,
                    haplotype
                )
            )
