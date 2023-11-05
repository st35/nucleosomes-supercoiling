import random

for i in range(32):
    G = []
    with open('configs/config_' + str(i) + '.config', 'r') as f:
        index = 0
        for line in f:
            l = line.strip().split('\t')
            G.append(float(l[4].strip()))
    I = [j for j in range(len(G)) if G[j] < 1.0]
    if len(I) == 0:
        I = [j for j in range(len(G)) if G[j] < 10.0]
    chosen_index = random.choice(I)

    with open('configs/config_' + str(i) + '.config', 'r') as f, open('configs_OE/config_' + str(i) + '.config', 'w') as g, open('configs_OE/OE_index_' + str(i) + '.config', 'w') as h:
        index = 0
        for line in f:
            l = line.strip().split('\t')
            g.write(l[0].strip() + '\t' + l[1].strip() + '\t' + l[2].strip() + '\t' + l[3].strip() + '\t')
            if index == chosen_index:
                g.write('100.0' + '\n')
                h.write(line.strip() + '\t' + str(chosen_index) + '\n')
            else:
                g.write(l[4].strip() + '\n')
            index += 1
