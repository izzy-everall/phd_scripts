import argparse
import sys
parser=argparse.ArgumentParser()

parser.add_argument('-i','--input',type=argparse.FileType('r'), help='input de gene file')
parser.add_argument('-a','--annotation',type=argparse.FileType('r'),help='gff annotation of reference genome')
parser.add_argument('-o','--output',type=argparse.FileType('w'),help='output file')
args=parser.parse_args()

################################
locus2product = {}
locus2gene = {}
for line in args.annotation:
	if line.startswith('_'):
		info = line.strip().split('\t')
		annotation_info = info[-1].split(';')
		loc = []
		gene = []
		for x in annotation_info:
			if x.startswith('ID'):
				wanted_loc =str(x[3:])
				loc.append(wanted_loc)
			else:
				pass

			if x.startswith('product='):
				product = x[8:]
			else:
				pass

			if x.startswith('gene='):
				gene_name = x[5:]
				gene.append(gene_name)
			else:
				pass
		locus2product[loc[0]]=(product)
		try:
			locus2gene[loc[0]]=(gene[0])
		except IndexError:
			locus2gene[loc[0]]=('NA')
	else:
		pass


count_out = 0
for line in args.input:
	count_out += 1
	if count_out == 1:
		header = line.split('\t')
		header.insert(1,'gene')
		header.insert(2, 'product')
		head_out = 0
		for x in header:
			head_out += 1
			if head_out == len(header):
				args.output.write(str(x))
			else:
				args.output.write(str(x)+'\t')
	else:
		count_info = line.split('\t')
		gene_name = locus2gene.get(count_info[0])
		product = locus2product.get(count_info[0])
		count_info.insert(1,gene_name)
		count_info.insert(2,product)
		hit_out = 0
		for x in count_info:
			hit_out += 1
			if hit_out == len(count_info):
				args.output.write(str(x))
			else:
				args.output.write(str(x)+'\t')





