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

name = names[-1]

gt = Util.read_ground_truth(path.join('./input/', name, 'reads_gt.txt'))
my = Util.read_ground_truth(path.join('./my_output', name, 'haplotypes.txt'))
print('my score: {}'.format(Util.earth_mover_distance(gt, my)))

gt_set = set([h for h, _ in gt])
get_set = set([h for h, _ in my])

for g_h, f in gt:
    print(g_h in get_set, f)

# print(len(set([h for h, _ in gt])))

# # name = 'bone_uniform_7_3'
# name = 'pigtail_uniform_6_3_1'
#
# gt = Util.read_ground_truth('./input/{}_gt.txt'.format(name))
# my = Util.read_ground_truth('./my_output/{}.txt'.format(name))
# vg_flow = Util.read_vgflow('./savage+vg-flow/{}/haps.final.fasta'.format(name))
#
# print('my score: {}'.format(Util.earth_mover_distance(gt, my)))
# print('vg_flow score: {}'.format(Util.earth_mover_distance(gt, vg_flow)))

# name = 'bone_uniform_7_3'