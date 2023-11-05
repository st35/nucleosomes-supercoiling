L = []
with open('Guo_data/gene_list.log', 'r') as f:
    for line in f:
        l = line.strip().split('\t')
        start = int(l[2].strip())
        E = int(l[3].strip())
        if start > E:
            print('Here.')
        L.append((E - start + 1))
print(float(sum(L)) / float(len(L)))
