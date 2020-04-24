from src import Util
from os import path


names = [
    'base_case_1_1',
    'base_case_2_2',
    'base_case_3_3',
    'base_case_4_2',
    'base_case_8_2',
    'base_case_8_4',
    'base_case_81_3',
    'base_case_81_6',
    'gattaca'
]


for name in names:
    print(name)
    gt = Util.read_ground_truth(path.join('./input/', name, 'reads_gt.txt'))
    my = Util.read_ground_truth(path.join(
        './my_output', name, 'haplotypes.txt'))
    savage = Util.read_vgflow(path.join(
            './savage+vg-flow', name, 'haps.final.fasta'))
    print('my score: {}'.format(Util.earth_mover_distance(gt, my)))
    print('savage score: {}'.format(Util.earth_mover_distance(gt, savage)))
    gt_set = set([h for h, _ in gt])
    get_set = set([h for h, _ in my])
    savage_set = set([h for h, _ in savage])
    print('whole reconstructed', len(gt_set) - len(gt_set - get_set))
    print('whole reconstructed by savage', len(gt_set) - len(gt_set - savage_set))
    print()
