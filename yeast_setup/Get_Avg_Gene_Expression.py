outputfile = open('Guo_data/expression_data.log', 'w')

with open('Guo_data/gene_list.log', 'r') as f:
	for line in f:
		l = line.strip().split('\t')
		gene_name = l[0].strip()
		chr_name = l[1].strip()
		start_gene = int(l[2].strip())
		end_gene = int(l[3].strip())
		dir_gene = int(l[4].strip())

		if chr_name == 'chrmt':
			break

		sig_sum = 0.0
		sig_norm_sum = 0.0
		l_sum = 0.0

		if dir_gene == 1:
			filename = 'Guo_data/RNA_Plus.log'
		else:
			filename = 'Guo_data/RNA_Minus.log'

		with open(filename, 'r') as g:
			count = 0
			for line2 in g:
				if count == 0:
					count += 1
					continue
				l2 = line2.strip().split(',')
				if chr_name not in l2[1].strip():
					continue
				pos = int(l2[2].strip())
				if pos > end_gene:
					break
				if pos >= start_gene and pos <= end_gene:
					sig_sum += float(l2[3].strip())
					sig_norm_sum += float(l2[4].strip())
					l_sum += 1.0
		print(gene_name)
		outputfile.write(gene_name + '\t' + str(sig_sum / l_sum) + '\t' + str(sig_norm_sum / l_sum) + '\n')

outputfile.close()
