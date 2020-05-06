from random import choices

nuc = "ATGC"

r_lean = 100


def rst(length):
    return "".join(choices(nuc, k=length))


ref = rst(500)
with open("ref.fa", 'w') as ref_file:
    ref_file.write(">ref\n{}".format(ref))

line_indel = ref[:100] + ref[110:]
point_indels = "".join(
    [ch for i, ch in enumerate(ref) if i not in {50, 150, 250, 400}]
)

haps = [
    (ref, 5),
    (line_indel, 2),
    (point_indels, 3)
]


with open("reads.fa", 'w') as reads_file:
    idx = 0
    for h, count in haps:
        for _ in range(count):
            for i in range(len(h) - r_lean + 1):
                idx += 1
                reads_file.write(">{}\n".format(idx))
                reads_file.write("{}\n".format(h[i: i + r_lean]))

with open("gt.txt", 'w') as gt:
    for h, count in haps:
        gt.write("{} {}\n".format(h, count))
