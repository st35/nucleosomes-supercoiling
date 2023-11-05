def Get_Data_Subset(start_gene, end_gene, index):
        Chr = []
        Direction = []

        flag = 0

        Starts = []
        Ends = []
        Names = []
        with open('Guo_data/gene_list.log', 'r') as f:
                for line in f:
                        l = line.strip().split('\t')
                        gene_name = l[0].strip()
                        if start_gene == gene_name and flag == 0:
                                flag = 1
                        if flag == 0:
                                continue
                        Chr.append(l[1].strip())
                        startp = 0
                        endp = 0
                        if l[4].strip() == '1':
                                startp = int(l[2].strip())
                                endp = int(l[3].strip())
                                Direction.append(1)
                        else:
                                startp = int(l[3].strip())
                                endp = int(l[2].strip())
                                Direction.append(-1)
                        Starts.append(startp)
                        Ends.append(endp)
                        Names.append(l[0].strip())
                        if end_gene == gene_name and flag == 1:
                                break
        Seg_Start = 0
        Seg_End = 0

        if Direction[0] == 1:
                Seg_Start = Starts[0]
        else:
                Seg_Start = Ends[0]

        if Direction[-1] == 1:
                Seg_End = Ends[-1]
        else:
                Seg_End = Starts[-1]

        chromosome = Chr[0]
        flag = 0
        with open('Guo_data/GapR_Flag.log', 'r') as f, open('configs_guo_data/' + 'config_' + str(index) + '_GapR_Flag.log', 'w') as g:
                count = 0
                for line in f:
                        if count == 0:
                                count += 1
                                continue
                        l = line.strip().split(',')
                        if l[1].strip() == chromosome:
                                pos = int(l[2].strip())
                                if pos >= Seg_Start and pos <= Seg_End:
                                        g.write(str(pos) + ' ' + l[3].strip() + ' ' + l[4].strip() + '\n')
                                if pos > Seg_End:
                                        break
        with open('Guo_data/GapR_Ctrl.log', 'r') as f, open('configs_guo_data/' + 'config_' + str(index) + '_GapR_Ctrl.log', 'w') as g:
                count = 0
                for line in f:
                        if count == 0:
                                count += 1
                                continue
                        l = line.strip().split(',')
                        if l[1].strip() == chromosome:
                                pos = int(l[2].strip())
                                if pos >= Seg_Start and pos <= Seg_End:
                                        g.write(str(pos) + ' ' + l[3].strip() + ' ' + l[4].strip() + '\n')
                                if pos > Seg_End:
                                        break
        with open('Guo_data/RNA_Plus.log', 'r') as f, open('configs_guo_data/' + 'config_' + str(index) + '_RNA_Plus.log', 'w') as g:
                count = 0
                for line in f:
                        if count == 0:
                                count += 1
                                continue
                        l = line.strip().split(',')
                        if l[1].strip() == chromosome:
                                pos = int(l[2].strip())
                                if pos >= Seg_Start and pos <= Seg_End:
                                        g.write(str(pos) + ' ' + l[3].strip() + ' ' + l[4].strip() + '\n')
                                if pos > Seg_End:
                                        break
        with open('Guo_data/RNA_Minus.log', 'r') as f, open('configs_guo_data/' + 'config_' + str(index) + '_RNA_Minus.log', 'w') as g:
                count = 0
                for line in f:
                        if count == 0:
                                count += 1
                                continue
                        l = line.strip().split(',')
                        if l[1].strip() == chromosome:
                                pos = int(l[2].strip())
                                if pos >= Seg_Start and pos <= Seg_End:
                                        g.write(str(pos) + ' ' + l[3].strip() + ' ' + l[4].strip() + '\n')
                                if pos > Seg_End:
                                        break

def Get_Level(level):
        if level < 0.0:
                return(0.01)
        if level < 1.0:
                return(0.1)
        if level < 2.0:
                return(1.0)
        if level < 3.0:
                return(10.0)
        return(100.0)

def Generate_Config_File(Chr, start_gene, end_gene, E, index):
        base_chr = Chr[start_gene]

        buffer = 10000
        TSS = []
        gene_length = []
        Direction = []
        flag = 0
        Starts = []
        Ends = []
        Names = []
        with open('Guo_data/gene_list.log', 'r') as f:
                for line in f:
                        l0 = line.strip().split('\t')
                        gene_name = l0[0].strip()
                        if start_gene == gene_name and flag == 0:
                                flag = 1
                        if flag == 0:
                                continue
                        l = line.strip().split('\t')
                        startp = 0
                        endp = 0
                        if l[4].strip() == '1':
                                startp = int(l[2].strip())
                                endp = int(l[3].strip())
                                Direction.append(1)
                        else:
                                startp = int(l[3].strip())
                                endp = int(l[2].strip())
                                Direction.append(-1)
                        TSS.append(startp)
                        gene_length.append(abs(endp - startp) + 1)
                        Starts.append(startp)
                        Ends.append(endp)
                        Names.append(l[0].strip())
                        if end_gene == gene_name:
                                break
        First_TSS = float("inf")
        for i in range(len(TSS)):
                if TSS[i] < First_TSS:
                        First_TSS = TSS[i]
        gene_count = 0
        with open('configs/' + 'config_' + str(index) + '.config', 'w') as f:
                Chosen_Names = []
                Chosen_Starts = []
                Chosen_Ends = []
                Chosen_Lengths = []
                Chosen_Directions = []
                Chosen_Levels = []
                for i in reversed(range(len(Starts))):
                        if i > 0:
                                x1 = 0
                                x2 = 0
                                if Direction[i] == 1:
                                        x1 = TSS[i]
                                        x2 = Ends[i]
                                else:
                                        x1 = Ends[i]
                                        x2 = TSS[i]
                                y1 = 0
                                y2 = 0
                                if Direction[i - 1] == 1:
                                        y1 = TSS[i - 1]
                                        y2 = Ends[i - 1]
                                else:
                                        y1 = Ends[i - 1]
                                        y2 = TSS[i - 1]
                                if x1 <= y2 and y1 <= x2:
                                        continue
                                overlap = min(x2, y2) - max(x1, y1)
                                if overlap > -10:
                                        continue
                        if i < len(Starts) - 1:
                                x1 = 0
                                x2 = 0
                                if Direction[i] == 1:
                                        x1 = TSS[i]
                                        x2 = Ends[i]
                                else:
                                        x1 = Ends[i]
                                        x2 = TSS[i]
                                y1 = 0
                                y2 = 0
                                if Direction[i + 1] == 1:
                                        y1 = TSS[i + 1]
                                        y2 = Ends[i + 1]
                                else:
                                        y1 = Ends[i + 1]
                                        y2 = TSS[i + 1]
                                if x1 <= y2 and y1 <= x2:
                                        continue
                                overlap = min(x2, y2) - max(x1, y1)
                                if overlap > -10:
                                        continue
                        if Chr[Names[i]] != base_chr:
                                print('Chromosome mismatch')
                                for gn in Names:
                                        print(gn + '\t' + Chr[gn])
                        Chosen_Names.append(Names[i])
                        Chosen_Starts.append(Starts[i])
                        Chosen_Ends.append(Ends[i])
                        Chosen_Lengths.append(gene_length[i])
                        Chosen_Directions.append(Direction[i])
                        Chosen_Levels.append(Get_Level(E[Names[i]]))
                Shift = Chosen_Starts[-1]
                if Chosen_Directions[-1] == -1:
                        Shift = Chosen_Starts[-1] - Chosen_Lengths[-1]
                if Shift < 0.0:
                        print('Calculated coordinates shift is invalid.')
                for i in range(len(Chosen_Names)):
                        f.write(Chosen_Names[i] + '\t' + str(Chosen_Starts[i] - Shift + buffer) + '\t' + str(Chosen_Lengths[i]) + '\t' + str(Chosen_Directions[i]) + '\t' + str(Chosen_Levels[i]) + '\n')

                return (Chosen_Names[-1], Chosen_Names[0])

from math import log10

chromosome = {}
with open('Guo_data/gene_list.log', 'r') as f:
        for line in f:
                l = line.strip().split('\t')
                chromosome[l[0].strip()] = l[1].strip()

E = {}
genes = []
with open('Guo_data/temp_expression_data.log', 'r') as f:
        for line in f:
                l = line.strip().split('\t')
                if float(l[1].strip()) == 0.0:
                        E[l[0].strip()] = 0.0
                else:
                        E[l[0].strip()] = log10(float(l[1].strip()))
                genes.append(l[0].strip())

from random import randint

for i in range(32):
        print(i)

        start_index = randint(0, len(genes) - 1)
        end_index = start_index + randint(5, 25)

        while end_index > len(genes) - 1:
                start_index = randint(0, len(genes) - 1)
                end_index = start_index + randint(5, 25)

        start_gene = genes[start_index]
        end_gene = genes[end_index]

        while (chromosome[start_gene] != chromosome[end_gene]) or (end_index > len(genes) - 1):
                start_index = randint(0, len(genes) - 1)
                end_index = start_index + randint(5, 25)

                start_gene = genes[start_index]
                end_gene = genes[end_index]

        real_start, real_end = Generate_Config_File(chromosome, start_gene, end_gene, E, i)
        Get_Data_Subset(real_start, real_end, i)
