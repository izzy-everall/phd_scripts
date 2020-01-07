import argparse
import sys
from collections import defaultdict
parser=argparse.ArgumentParser()

parser.add_argument('-i','--input',type=argparse.FileType('r'), help='input de gene file')
parser.add_argument('-a','--annotation',type=argparse.FileType('r'),help='gff annotation of reference genome')
parser.add_argument('-o','--output',type=argparse.FileType('w'),help='output file')
parser.add_argument('-v','--verbose', type=argparse.FileType('w'),help='genes not counted')
parser.add_argument('-p','--pangenome', type=argparse.FileType('r'),help='pangenome_gene_pres_absence')
parser.add_argument('-l','--locus_product',type=argparse.FileType('r'),help='locus2product')
args=parser.parse_args()
###############################################
#########################################################
######### reading in annotation information #############
#########################################################
locus2start_end={}
locus2orientation = {}
locus2gene_length = {}
locus2product = {}
locus2gene = {}
locus2all_info = {}
int_end2start = {}
previous_line = []

def organise_annotation_data(annotation):
	for line in annotation:
		if line.startswith('_12163'):
			info = line.strip().split('\t')
			# so don't mess up coords with tRNAs etc
			if info[2] == 'CDS' or info[2] == 'ncRNA':

				contig = info[0]
				start = info[3]
				end = info[4]
				start_end = str(start)+'..'+str(end)
				orientation = info[6]
				annotation_info = info[-1].split(';')
				loc = []
				for x in annotation_info:
					if x.startswith('ID'):
						wanted_loc = str(contig)+'_'+str(x[3:])
						loc.append(wanted_loc)
					else:
						pass

					#if x.startswith('gene'):
					#	gene = x[5:]
					#else:
					#	pass

					if x.startswith('product='):
						product = x[8:]
					else:
						pass


				if len(previous_line) == 0:
					previous_line.append(loc[0])
					previous_line.append(end)
				else:
					# do something with previous line
					gap = float(start) - float(previous_line[1])
					if float(gap) > 0:
						region = previous_line[0]+'..'+loc[0]
						coords = str(previous_line[1])+'..'+str(start)
						int_end2start[region]=(coords)
					else:
						pass
					# clear previous line
					del previous_line[:]
					# put current line info in
					previous_line.append(loc[0])
					previous_line.append(end)


				locus2product[loc[0]]=(product)
				locus2orientation[loc[0]]=(orientation)
				locus2all_info[loc[0]]=(info)
				locus2start_end[loc[0]]=(start_end)
				#locus2gene[loc[0]]=(gene)

			else:
				if args.verbose:
					args.verbose.write(line)
				else:
					pass


organise_annotation_data(args.annotation)

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


##### annotations #######
for line in args.locus_product:
	info = line.strip().split('\t')
	pan_locus2annotation[info[0]]=(info[1])

################################################
#locus2product = {}
#locus2gene = {}
#for line in args.annotation:
#	if line.startswith('_'):
#		info = line.strip().split('\t')
#		annotation_info = info[-1].split(';')
#		loc = []
#		gene = []
#		for x in annotation_info:
#			if x.startswith('ID'):
#				wanted_loc =str(x[3:])
#				loc.append(wanted_loc)
#			else:
#				pass
#
#			if x.startswith('product='):
#				product = x[8:]
#			else:
#				pass
#
#			if x.startswith('gene='):
#				gene_name = x[5:]
#				gene.append(gene_name)
#			else:
#				pass
#		locus2product[loc[0]]=(product)
#		try:
#			locus2gene[loc[0]]=(gene[0])
#		except IndexError:
#			locus2gene[loc[0]]=('NA')
#	else:
#		pass


count_out = 0
for line in args.input:
	count_out += 1
	if count_out == 1:
		header = line.split('\t')
		#header.insert(1,'gene')
		header.insert(1, 'product')
		header.insert(2,'intergenic_flanking_loci')
		header.insert(3,'MABS_locus')
		header.insert(4,'MABS_product')
		head_out = 0
		for x in header:
			head_out += 1
			if head_out == len(header):
				args.output.write(str(x))
			else:
				args.output.write(str(x)+'\t')
	else:
		count_info = line.split('\t')
		out_info = []

		
		for l, p in locus2product.iteritems():
			if count_info[0].startswith('Mycobacterium') and l.endswith(count_info[0]): 
				#gene_name = g
				product = locus2product.get(l)
				out_info.append('NA')
				out_info.append(product)
				# not needed any more - actually is
				#adj_loci = count_info[0]+'"'
				#print count_info[0]
				adj_loci = count_info[0]
				get_mabs_loc = locus2other_loci.get(adj_loci)

				# loop through to see if mabs loc there
				mabs_loc = []
				for x in get_mabs_loc:
					if x.startswith('MAB'):
						mabs_loc.append(x)
					else:
						pass
	
				if len(mabs_loc) == 0:
				#if mabs_loc == '':
					out_info.append('NA')
					out_info.append('NA')
				else:
					mabs_prod = pan_locus2annotation.get(mabs_loc[0])
					out_info.append(mabs_loc[0])
					out_info.append(mabs_prod)
			else:
				pass

		#count_info.insert(1,gene_name)
		#count_info.insert(2,product)

		for i, coord in int_end2start.iteritems():
			#intergenic coordinages
			get_contig = count_info[0].split('_')
			start_contig = get_contig[0]+'_'+get_contig[1]+'_'+get_contig[2]+'_'+get_contig[3]+'_'+get_contig[4]
			if count_info[0].startswith('_') and i.startswith(start_contig):
				# deseq intergenic coord
				de_start_coord = count_info[0].split('_')[-2]
				de_end_coord = count_info[0].split('_')[-1]
				an_start_coord = coord.split('..')[0]
				an_end_coord = coord.split('..')[1]
				if float(de_start_coord) >= float(an_start_coord) and float(de_end_coord) <= float(an_end_coord):
					out_info.append(i)

				else:
					pass

		# intergenic region info
		if len(out_info) == 1:
			#count_info.insert(1,'NA')
			count_info.insert(1,'NA')
			count_info.insert(2,out_info[0])
			count_info.insert(3,'NA')
			count_info.insert(4,'NA')

		# 
		else:
			#count_info.insert(1,gene_name)
			count_info.insert(1,out_info[1])
			count_info.insert(2,'NA')
			count_info.insert(3,out_info[-2])
			count_info.insert(4,out_info[-1])

		
		hit_out = 0
		for x in count_info:
			hit_out += 1
			if hit_out == len(count_info):
				args.output.write(str(x))
			else:
				args.output.write(str(x)+'\t')





