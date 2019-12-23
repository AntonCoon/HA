from numpy.random import uniform
# from numpy.random import normal
from numpy.random import dirichlet
from numpy.random import choice as np_choice
from typing import List
from typing import Callable
from matplotlib import pyplot as plt
from os import path, mkdir
from itertools import product


class Haplotype(object):
    def __init__(self, hap: str, freq: float):
        self.hap = hap
        self.freq = freq


Haplotypes = List[Haplotype]


class ReadMaker(object):
    def __init__(self, distribution: Callable, read_len=150, save_path='.'):
        self.read_len = read_len
        self.distribution = distribution
        self.read_positions = None
        self.path = save_path

    def make_reads(self, name: str, gt: Haplotypes, amount=1000) -> None:
        gt_path = path.join(self.path, '{}_gt.txt'.format(name))
        with open(gt_path, 'w') as gt_file:
            for hap in gt:
                gt_file.write('{} {}\n'.format(hap.hap, hap.freq))
        frequencies = [h.freq for h in gt]
        gt_len = len(gt[0].hap)
        read_file = path.join(self.path, '{}.fastq'.format(name))
        with open(read_file, 'w') as reads:
            low = 0
            high = gt_len - self.read_len + 1
            self.read_positions = self.distribution(low, high, amount)
            for idx, pos in enumerate(self.read_positions):
                pos = int(pos)
                haplotype = np_choice(gt, 1, p=frequencies)[0]
                read = haplotype.hap[pos: pos + self.read_len]
                reads.write('@{}\n{}\n'.format(str(idx), read))
                reads.write('+\n{}\n'.format('I' * self.read_len))


def hap_gen(ref: List[str], pos: List[int], nucls: List[list]) -> List[str]:
    haps = []
    for substitutions in product(*nucls):
        haplotype = ref[:]
        for idx, nucl in zip(pos, substitutions):
            haplotype[idx] = nucl
        haps.append(''.join(haplotype))

    return haps


def special_hap_gen(ref: List[str], pos: List[int], nucls: List[list]) -> List:
    haps = []
    for substitutions in nucls:
        haplotype = ref[:]
        for idx, nucl in zip(pos, substitutions):
            haplotype[idx] = nucl
        haps.append(''.join(haplotype))

    return haps


def read_ref(ref_path: str) -> List[str]:
    ref = []
    with open(ref_path, 'r') as ref_file:
        ref.extend(next(ref_file).strip())

    return ref


def bone_uniform_example():
    bone_name = 'bone_uniform_7_3'
    # classic bone example
    ref = read_ref('./input/ref.txt')
    first, second = len(ref) // 10, len(ref) - len(ref) // 10
    h1 = ''.join(
        ref[:first] + ['T'] + ref[first:second] + ['T'] + ref[second:]
    )
    h2 = ''.join(
        ref[:first] + ['A'] + ref[first:second] + ['A'] + ref[second:]
    )

    freqs = [.7, .3]
    grt = [Haplotype(h, f) for h, f in zip([h1, h2], freqs)]

    uniform_read_maker = ReadMaker(uniform, read_len=200, save_path='./input')
    uniform_read_maker.make_reads(bone_name, grt, amount=1000)

    plt.hist(uniform_read_maker.read_positions)
    plt.xlabel('read start position')
    plt.savefig(path.join('./input', bone_name + '.jpg'))


def pigtail_uniform_example():
    bone_name = 'pigtail_32'
    # pigtail example with 2^5 possible paths and 8 in ground truth
    #    'GTTTGCCAGGA_GATGGAAACCA_AAATGATAGGG_GAATTGGAGGT_TTATCAAAGTA_GACAGTAT'
    haps = [
        'GTTTGCCAGGATGATGGAAACCATAAATGATAGGGAGAATTGGAGGTATTATCAAAGTAAGACAGTAT',
        'GTTTGCCAGGAAGATGGAAACCATAAATGATAGGGAGAATTGGAGGTATTATCAAAGTATGACAGTAT',
        'GTTTGCCAGGATGATGGAAACCATAAATGATAGGGAGAATTGGAGGTTTTATCAAAGTAAGACAGTAT',
        'GTTTGCCAGGAAGATGGAAACCAAAAATGATAGGGAGAATTGGAGGTTTTATCAAAGTATGACAGTAT',
        'GTTTGCCAGGAAGATGGAAACCAAAAATGATAGGGTGAATTGGAGGTATTATCAAAGTAAGACAGTAT',
        'GTTTGCCAGGATGATGGAAACCAAAAATGATAGGGTGAATTGGAGGTATTATCAAAGTATGACAGTAT',
        'GTTTGCCAGGAAGATGGAAACCAAAAATGATAGGGTGAATTGGAGGTTTTATCAAAGTAAGACAGTAT',
        'GTTTGCCAGGATGATGGAAACCAAAAATGATAGGGTGAATTGGAGGTTTTATCAAAGTATGACAGTAT'
    ]
    freqs = dirichlet([1] * len(haps), 1)[0]

    grt = [Haplotype(h, f) for h, f in zip(haps, freqs)]

    uniform_read_maker = ReadMaker(
        uniform, read_len=20, save_path='../example'
    )
    uniform_read_maker.make_reads(bone_name, grt, amount=30000)

    plt.hist(uniform_read_maker.read_positions)
    plt.xlabel('read start position')
    plt.savefig(path.join(uniform_read_maker.path, bone_name + '.jpg'))


def big_pigtail_uniform_example():
    bone_name = 'pigtail_uniform_6_3_1'
    # pigtail example with 2x2x3 possible paths
    ref = read_ref('./input/ref.txt')
    first, second = len(ref) // 10, len(ref) - len(ref) // 10
    third = len(ref) // 2
    h1 = ''.join(
        ref[:first] + ['T'] +
        ref[first:second] + ['T'] +
        ref[second:third] + ['T'] +
        ref[third:]
    )
    h2 = ''.join(
        ref[:first] + ['T'] +
        ref[first:second] + ['A'] +
        ref[second:third] + ['A'] +
        ref[third:]
    )

    h3 = ''.join(
        ref[:first] + ['G'] +
        ref[first:second] + ['T'] +
        ref[second:third] + ['G'] +
        ref[third:]
    )

    freqs = [.6, .3, .1]
    grt = [Haplotype(h, f) for h, f in zip([h1, h2, h3], freqs)]

    uniform_read_maker = ReadMaker(uniform, read_len=200, save_path='./input')
    uniform_read_maker.make_reads(bone_name, grt, amount=3000)

    plt.hist(uniform_read_maker.read_positions)
    plt.xlabel('read start position')
    plt.savefig(path.join('./input', bone_name + '.jpg'))


def full_case_generator(case_name: str, pos_n: int, subs: List[list]) -> None:
    save_path = path.join('./input', case_name)
    mkdir(save_path)

    ref = read_ref('./input/ref.txt')
    pos = [len(ref) // (pos_n + 1) * idx for idx in range(1, pos_n + 1)]
    haps = hap_gen(ref, pos, subs)

    freqs = dirichlet([1] * len(haps), 1)[0]

    grt = [Haplotype(h, f) for h, f in zip(haps, freqs)]

    uniform_read_maker = ReadMaker(uniform, read_len=200, save_path=save_path)
    uniform_read_maker.make_reads('reads', grt, amount=3000)

    plt.hist(uniform_read_maker.read_positions)
    plt.xlabel('read start position')
    plt.savefig(path.join(save_path, 'hist.jpg'))


def part_case_generator(name: str, subs: List) -> None:
    save_path = path.join('./input', name)
    mkdir(save_path)

    pos_n = len(subs[0])
    ref = read_ref('./input/ref.txt')
    pos = [len(ref) // (pos_n + 1) * idx for idx in range(1, pos_n + 1)]
    haps = special_hap_gen(ref, pos, subs)
    freqs = dirichlet([1] * len(haps), 1)[0]

    grt = [Haplotype(h, f) for h, f in zip(haps, freqs)]

    uniform_read_maker = ReadMaker(uniform, read_len=200, save_path=save_path)
    uniform_read_maker.make_reads('reads', grt, amount=10000)

    plt.hist(uniform_read_maker.read_positions)
    plt.xlabel('read start position')
    plt.savefig(path.join(save_path, 'hist.jpg'))


def case_generator_with_freqs(name: str, subs: List, freqs: List) -> None:
    save_path = path.join('./input', name)
    mkdir(save_path)

    pos_n = len(subs[0])
    ref = read_ref('./input/ref.txt')
    pos = [len(ref) // (pos_n + 1) * idx for idx in range(1, pos_n + 1)]
    haps = special_hap_gen(ref, pos, subs)

    grt = [Haplotype(h, f) for h, f in zip(haps, freqs)]

    uniform_read_maker = ReadMaker(uniform, read_len=200, save_path=save_path)
    uniform_read_maker.make_reads('reads', grt, amount=3000)

    plt.hist(uniform_read_maker.read_positions)
    plt.xlabel('read start position')
    plt.savefig(path.join(save_path, 'hist.jpg'))


def base_case_1_1():
    case_name = 'base_case_1_1'
    full_case_generator(case_name, 1, [['A']])


def base_case_2_2():
    case_name = 'base_case_2_2'
    full_case_generator(case_name, 1, [['A', 'T']])


def base_case_3_3():
    case_name = 'base_case_3_3'
    full_case_generator(case_name, 1, [['A', 'T', 'G']])


def base_case_4_2():
    case_name = 'base_case_4_2'
    subs = [
        ['A', 'A'],
        ['T', 'T']
    ]
    part_case_generator(case_name, subs)


def base_case_4_4():
    case_name = 'base_case_4_4'
    subs = [['A', 'T']] * 2
    full_case_generator(case_name, 2, subs)


def base_case_8_2():
    case_name = 'base_case_8_2'
    subs = [
        ['A', 'A', 'A'],
        ['T', 'T', 'T']
    ]
    part_case_generator(case_name, subs)


def base_case_8_4():
    case_name = 'base_case_8_4'
    subs = [
        ['A', 'A', 'A'],
        ['A', 'A', 'T'],
        ['A', 'T', 'T'],
        ['T', 'T', 'T']
    ]
    part_case_generator(case_name, subs)


def base_case_81_3():
    # maximum of nonzero frequency amount = 4 * (3 - 1) + 1
    case_name = 'base_case_81_3'
    subs = [
        ['A', 'A', 'A', 'A'],
        ['T', 'T', 'T', 'T'],
        ['G', 'G', 'G', 'G']
    ]
    part_case_generator(case_name, subs)


def base_case_81_6():
    # maximum of nonzero frequency amount = 4 * (3 - 1) + 1
    case_name = 'base_case_81_6'
    subs = [
        ['ATG', 'ATG', 'ATG', 'ATG'],
        ['ATG', 'ATG', 'TGC', 'TGC'],
        ['TGC', 'TGC', 'TGC', 'TGC'],
        ['TGC', 'TGC', 'GCT', 'GCT'],
        ['GCT', 'GCT', 'GCT', 'GCT'],
        ['GCT', 'GCT', 'ATG', 'ATG']
    ]
    part_case_generator(case_name, subs)


def case_gattaca():
    # maximum of nonzero frequency amount = 4 * (3 - 1) + 1
    case_name = 'gattaca'
    subs = [
        ['G', 'A', 'T', 'T', 'A', 'C', 'A'],
        ['G', 'T', 'T', 'A', 'C', 'A', 'T'],
        ['C', 'C', 'C', 'A', 'G', 'A', 'T']
    ]
    case_generator_with_freqs(case_name, subs, [.4, .35, .25])


if __name__ == '__main__':

    # base_case_1_1()
    # base_case_2_2()
    # base_case_3_3()
    # base_case_4_4()
    # base_case_4_4()
    # base_case_8_2()
    # base_case_8_4()
    # base_case_81_3()
    # base_case_81_6()
    case_gattaca()
