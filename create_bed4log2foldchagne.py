import argparse
import sys
from collections import defaultdict
parser=argparse.ArgumentParser()

parser.add_argument('-i','--input',type=argparse.FileType('r'), help='input de gene file')
parser.add_argument('-a','--annotation',type=argparse.FileType('r'),help='gff annotation of reference genome')
parser.add_argument('-o','--output',type=argparse.FileType('w'),help='output file')
#parser.add_argument('-v','--verbose', type=argparse.FileType('w'),help='genes not counted')
parser.add_argument('-p','--pangenome', type=argparse.FileType('r'),help='pangenome_gene_pres_absence')
parser.add_argument('-r','--rectangle',type=argparse.FileType('w'),help='rectangle.bed')
args=parser.parse_args()
###############################################
args.output.write('chrm'+'\t'+'start'+'\t'+'end'+'\t'+'log2fold'+'\t'+'padj'+'\n')

#########################################################
######### reading in annotation information #############
#########################################################
locus2start_end = {}

for line in args.annotation:
	info = line.strip().split('\t')
	coords = info[1]+'..'+info[2]
	locus2start_end[info[0]]=(coords)

#################################################
############ reading in pangenome file #########
################################################
#now want to read in pangenome file to identify corresponding illumina locus tag.
pacbio_WT_locus2other_loci = {}
pacbio_KO_locus2other_loci = {}
illumina_locus2other_loci = {}
mabs_locus2other_loci = {}

locus2other_loci = defaultdict(list)

pan_locus2annotation = {}
count_pan = 0
for line in args.pangenome:
	count_pan += 1
	if count_pan == 1:
		pass
	else:
		info = line.strip().split('\t')

		for x in info:
			for l in info:
				if x == l:
					pass
				else:
					locus2other_loci[x].append(l)

rm_header = 0
for line in args.input:
	info = line.strip().split('\t')
	rm_header += 1
	if rm_header == 1:
		pass
	else:
		if line.startswith('_'):
			locus_needed = info[1].split('..')
			start_locus = locus_needed[0][16:]
			end_locus = locus_needed[1][16:]

			start_bir1049_wt = locus2other_loci.get(start_locus)
			end_bir1049_wt = locus2other_loci.get(end_locus)

			start_wanted = []
			end_wanted = []
			for x in start_bir1049_wt:
				if x.startswith('BIR1049_WT'):
					start_wanted.append(x)
				else:
					pass

			for x in end_bir1049_wt:
				if x.startswith('BIR1049_WT'):
					end_wanted.append(x)
				else:
					pass


			get_start_coords = locus2start_end.get(start_wanted[0])
			get_end_coords = locus2start_end.get(end_wanted[0])

			#print get_start_coords, get_end_coords
			start_int = get_start_coords.split('..')[1]
			end_int = get_end_coords.split('..')[0]

			end_minus_start = int(end_int) - int(start_int)
			mid = int(end_minus_start)/2

			plot_point = int(start_int) + int(mid)

			args.output.write('chrm1'+'\t'+str(plot_point)+'\t'+str(plot_point)+'\t'+str(info[-2])+'\t'+str(info[-1])+'\n')

			if args.rectangle:
				args.rectangle.write('chrm1'+'\t'+str(start_int)+'\t'+str(end_int)+'\t'+str(info[-2])+'\t'+str(info[-1])+'\n')


		

		else:
			other_loci = locus2other_loci.get(info[0])
			loci_wanted = []
			for x in other_loci:
				if x.startswith('BIR1049_WT'):
					loci_wanted.append(x)
				else:
					pass


			coords = locus2start_end.get(loci_wanted[0])

			start_coord = coords.split('..')[0]
			end_coord = coords.split('..')[1]

			end_minus_start = int(end_coord) - int(start_coord)
			mid = int(end_minus_start)/2

			plot_point = int(start_coord) + int(mid)

			args.output.write('chrm1'+'\t'+str(plot_point)+'\t'+str(plot_point)+'\t'+str(info[-2])+'\t'+str(info[-1])+'\n')

			if args.rectangle:
				args.rectangle.write('chrm1'+'\t'+str(start_coord)+'\t'+str(end_coord)+'\t'+str(info[-2])+'\t'+str(info[-1])+'\n')



