import argparse
from collections import defaultdict
import operator
import sys

parser=argparse.ArgumentParser()
parser.add_argument('-i','--input',type=argparse.FileType('rU'),nargs='+',help='list of replicates read count files to merge')
parser.add_argument('-o','--output',type=argparse.FileType('w'),help='combined read count file for DeSEQ1 input file')
parser.add_argument('-d','--deseq',action='store_true',help='if you want just seperate read count file for DESEQ2')
parser.add_argument('-p','--prefix',type=str, help='output file prefix or suffix to output file name')
parser.add_argument('-x','--exclude_intergenic',action='store_true',help='do not include int regions')
args=parser.parse_args()


all_isols2isol2read_counts = defaultdict(dict)
all_regions = []
files = []
file_count = 0


for f in args.input:
	gene2count = []
	file_name = str(f).split(' ')
	isol_name = file_name[2].split('.')

	try:
		final_isol_name = isol_name[0][1:]
		print final_isol_name
	except IndexError:
		print 'Do your input files contain a "_" seperating the name from the replicate? or are the replicate names unique?'
		sys.exit()

	files.append(final_isol_name)
	file_count += 1

	if args.deseq:
		out_file_name=(final_isol_name+'_'+str(args.prefix)+'.txt')
		sep_out_file = open(out_file_name,'w')

	count_line = 0
	for line in f:
		count_line += 1
		if count_line == 1:
			pass
		else:
			info = line.strip().split(',')
			if len(info) == 1:
				pass
			else:
				if info[3] == 'CDS':
					all_isols2isol2read_counts[final_isol_name][info[2]]=(info[4])
					if args.deseq:
						gene_count = (info[2],info[4])
						gene2count.append(gene_count)

					if file_count == 1:
						all_regions.append(info[2])

				else:
					all_isols2isol2read_counts[final_isol_name][info[1]]=(info[4])
					if args.deseq:
						gene_count = (info[1],info[4])
						gene2count.append(gene_count)

					if file_count == 1:
						all_regions.append(info[1])


	sorted_gene2count = sorted(gene2count, key=lambda tup: (tup[0]))
	for x in sorted_gene2count:
		sep_out_file.write(x[0]+'\t'+x[1]+'\n')

if args.output:
	args.output.write('Gene'+'\t'+str(files[0])+'\t'+str(files[1])+'\t'+str(files[2])+'\t'+str(files[3])+'\t'+str(files[4])+'\t'+str(files[5])+'\n')
	for g in all_regions:
		# don't want intergenic 
		if args.exclude_intergenic:
			if 'intergenic' in g:
				pass
			else:
				first_rep = all_isols2isol2read_counts.get(files[0])
				first_rep_gene = first_rep.get(g)
				second_rep = all_isols2isol2read_counts.get(files[1])
				second_rep_gene = second_rep.get(g)
				third_rep = all_isols2isol2read_counts.get(files[2])
				third_rep_gene = third_rep.get(g)
				fourth_rep = all_isols2isol2read_counts.get(files[3])
				fourth_rep_gene = fourth_rep.get(g)
				fifth_rep = all_isols2isol2read_counts.get(files[4])
				fifth_rep_gene = fifth_rep.get(g)
				sixth_rep = all_isols2isol2read_counts.get(files[5])
				sixth_rep_gene = sixth_rep.get(g)
				args.output.write(g+'\t'+first_rep_gene+'\t'+second_rep_gene+'\t'+third_rep_gene+'\t'+fourth_rep_gene+'\t'+fifth_rep_gene+'\t'+sixth_rep_gene+'\n')
		else:
			first_rep = all_isols2isol2read_counts.get(files[0])
			first_rep_gene = first_rep.get(g)
			second_rep = all_isols2isol2read_counts.get(files[1])
			second_rep_gene = second_rep.get(g)
			third_rep = all_isols2isol2read_counts.get(files[2])
			third_rep_gene = third_rep.get(g)
			fourth_rep = all_isols2isol2read_counts.get(files[3])
			fourth_rep_gene = fourth_rep.get(g)
			fifth_rep = all_isols2isol2read_counts.get(files[4])
			fifth_rep_gene = fifth_rep.get(g)
			sixth_rep = all_isols2isol2read_counts.get(files[5])
			sixth_rep_gene = sixth_rep.get(g)
			args.output.write(g+'\t'+first_rep_gene+'\t'+second_rep_gene+'\t'+third_rep_gene+'\t'+fourth_rep_gene+'\t'+fifth_rep_gene+'\t'+sixth_rep_gene+'\n')











