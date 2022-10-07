import sys
import itertools
import re

# hard-coded k-1
K = int(sys.argv[4]) 
KM = K - 1 

gfa = dict()
with open(sys.argv[1], 'r') as f:
    for line in f:
        line = line.rstrip().split()
        if line[0] == 'S':
            gfa[int(line[1])] = line[2]
        elif line[0] == 'L':
            KM = int(line[-1][:-1])

with open(sys.argv[2], 'r') as f:
    headers = [(header.rstrip()[1:],seq.rstrip()) for header,seq in itertools.zip_longest(*[f]*2)]

outlines = [[f'>{header}::0:{len(seq)}S:N'] for header,seq in headers]
header_map = dict((header,i) for i,(header,seq) in enumerate(headers))
scores = [[-1000]] * len(outlines)

def seq_rc(seq):
    old_chars = "ACGT"
    replace_chars = "TGCA"
    tab = str.maketrans(old_chars,replace_chars)
    return seq.translate(tab)[::-1]


with open(sys.argv[3], 'r') as f:
    for line in f:
        line = line.rstrip().split()
        name = line[0]
        qlen = int(line[1])
        qbegin = int(line[2])
        qend = int(line[3])
        assert(line[4] == '+')
        path = line[5]
        plen = int(line[6])
        pbegin = int(line[7])
        pend = int(line[8])
        orientation = "".join(re.findall("[<>]", path))
        path = path.replace("<", ">").split(">")[1:]
        unitigs = [gfa[int(a)] for a in path]
        unitigs = [u if o == ">" else seq_rc(u) for u,o in zip(unitigs, orientation)]
        path = unitigs[0] + "".join([u[KM:] for u in unitigs[1:]])
        assert(len(path) == plen)
        path = path[pbegin:pend]
        if line[16][:2] == "cg":
            cigar = line[16].split(":")[2:][0]
        else:
            cigar = f'{qend-qbegin}='
        if qbegin:
            cigar = f'{qbegin}S' + cigar
        if qend != qlen:
            cigar += f'{qlen - qend}S'
        score = int(float(line[13].split(":")[-1]))
        i = header_map[name]
        if score > scores[i][0]:
            outlines[i].append(f'>{name}:{path}:{score}:{cigar}:{path}')
            scores[i].append(score)

for line in outlines:
    if len(line) == 1:
        print(line[0])
    else:
        print("\n".join(set(line[1:])))
