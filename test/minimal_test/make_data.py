from random import choices

nuc = "ATGC"

r_lean = 100


def rst():
    return "".join(choices(nuc, k=120))


left, right = rst(), rst()

haps = [
    ("ATGTCT", 1),
    ("ATGCAT", 1),
    ("ATACAT", 3)
]
haps = [(left + h + right, c) for h, c in haps]

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

with open("ref.fa", 'w') as gt:
    gt.write(">ref\n{}".format(haps[-1][0]))
