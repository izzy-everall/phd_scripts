import argparse
parser=argparse.ArgumentParser()
parser.add_argument('-i','--input',type=argparse.FileType('r'),help='modifications.csv')
parser.add_argument('-g','--gff',type=argparse.FileType('r'),help='modifications.gff')
parser.add_argument('-m','--motifs',type=argparse.FileType('r'),help='motifs.gff')
parser.add_argument('-l','--lists',type=argparse.FileType('r'),help='file with motifs of interest')
parser.add_argument('-o','--output',type=argparse.FileType('w'),help='output file')
args=parser.parse_args()

###############
#output 
args.output.write('modified'+'\t'+'modification_type'+'\t'+'ipdRatio'+'\t'+'base'+'\t'+'score'+'\t'+'motif'+'\n')
#################

#########################
contig2pos2ipdratio = {}
contig2pos2type = {}
for line in args.gff:
	if line.startswith('#'):
		pass
	else:
		info = line.strip().split('\t')
		contig = info[0].split('|')[0]
		contig2pos = contig+'..'+info[3]
		contig2pos2type[contig2pos]=(info[2])
		ratio = info[-1].split(';')
		for x in ratio:
			if x.startswith('IPDRatio'):
				contig2pos2ipdratio[contig2pos]=(x[9:])
			else:
				pass




########################################################
############# file with motifs of interest in ##########
########################################################
motif_list = []
for line in args.lists:
	info = line.strip()
	motif_list.append(info)


#############################################
######## read in motif gff file #############
#############################################
motif_contig2pos2strand_mod_type = {}
motif_contig2pos2strand_score = {}
motif_contig2pos2strand_motif = {}
motif_contig2pos2strand_mod_nonsense = {}
for line in args.motifs:
	if line.startswith('#'):
		pass
	else:
		info = line.strip().split('\t')
		contig = info[0].split('|')[0]
		contig2pos2strand = contig+'..'+info[3]+'..'+info[6]
		ratio = info[-1].split(';')
		for x in ratio:
			if x.startswith('motif'):
				motif = x[6:]
				all_motifs = motif.split(',')
				for i in all_motifs:
					if i in motif_list:
						motif_contig2pos2strand_mod_type[contig2pos2strand]=(info[2])
						motif_contig2pos2strand_score[contig2pos2strand]=(info[5])
						motif_contig2pos2strand_motif[contig2pos2strand]=(i)
					else:
						# if modified base predicted
						if info[2] != '.':
							motif_contig2pos2strand_mod_nonsense[contig2pos2strand]=(i)
						else:
							pass
			else:
				pass


##############################################################
############### READ IN MODIFICATIONS.CSV FILE ###############
##############################################################
csv_contig2pos_mod_base = {}
csv_contig2pos_ipdR = {}
csv_contig2pos_score = {}
header=0
for line in args.input:
	header +=1
	if header == 1:
		pass
	else:
		info = line.strip().split(',')
		contig = info[0].split('|')[0][1:]
		if info[2] == str(0):
			contig2pos = contig+'..'+info[1]+'..+'
		else:
			contig2pos = contig+'..'+info[1]+'..-'

		csv_contig2pos_mod_base[contig2pos]=(info[3])
		csv_contig2pos_ipdR[contig2pos]=(info[8])
		csv_contig2pos_score[contig2pos]=(info[4])



for k, v in csv_contig2pos_mod_base.iteritems():
	mod_type = motif_contig2pos2strand_mod_type.get(k)
	mod_nonsense = motif_contig2pos2strand_mod_nonsense.get(k)
	modified_ipdR = csv_contig2pos_ipdR.get(k)
	score = csv_contig2pos_score.get(k)
	motif = motif_contig2pos2strand_motif.get(k)
	# modification type was not a base in your list of motifs you trust
	if mod_type == None or mod_type == '.':
		# remove modifications that are just complete nonsense
		if mod_nonsense == None:
			print 'true_unmodified', k, v, mod_type
			un_modified_idpR = csv_contig2pos_ipdR.get(k)
			args.output.write('unmodified'+'\t'+str(v)+'\t'+str(un_modified_idpR)+'\t'+str(v)+'\t'+str(score)+'\t'+'NA'+'\n')

		else:
			pass

	else:
		print k, v, mod_type
		args.output.write('modified'+'\t'+str(mod_type)+'\t'+str(modified_ipdR)+'\t'+str(v)+'\t'+str(score)+'\t'+str(motif)+'\n')





