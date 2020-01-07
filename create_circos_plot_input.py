import argparse
import sys
parser=argparse.ArgumentParser()
parser.add_argument('-i','--input',type=argparse.FileType('r'),help='motifs of interest')
parser.add_argument('-m','--modifications',type=argparse.FileType('r'),help='modifications.csv')
parser.add_argument('-d','--motifs',type=argparse.FileType('r'),help='motifs.gff')
parser.add_argument('-f','--forward',type=argparse.FileType('w'),help='forward.bed')
parser.add_argument('-r','--rect',type=argparse.FileType('w'),help='rectangular output plot')
parser.add_argument('-o','--output',type=argparse.FileType('w'),help='output.bed')
parser.add_argument('-c','--candidates',type=argparse.FileType('r'),help='motif flanking output')
#parser.add_argument('-w','--wanted',tyep=argparse.FileType('r'))
args=parser.parse_args()

pos_strand2idpR = {}
modcsv_contig2base = {}

header = 0
for line in args.modifications:
	header += 1
	if header == 1:
		args.output.write('chrm'+'\t'+'start_pos'+'\t'+'start_end'+'\t'+'value'+'\t'+'strand'+'\n')
	else:
		info = line.strip().split(',')
		if info[0].startswith('"unitig_80'):
			if int(info[2]) == 0:
				args.output.write('chrm1'+'\t'+str(info[1])+'\t'+str(info[1])+'\t'+str(info[8])+'\t'+str('+')+'\n')
				contig_pos_strand = 'chrm1'+'_'+str(info[1])+'_+'
				pos_strand2idpR[contig_pos_strand]=(str(info[8]))
				modcsv_contig2base[contig_pos_strand]=(info[3])
			else:
				args.output.write('chrm1'+'\t'+str(info[1])+'\t'+str(info[1])+'\t'+str(info[8])+'\t'+str('-')+'\n')
				contig_pos_strand = 'chrm1'+'_'+str(info[1])+'_-'
				new_val = '-'+str(info[8])
				pos_strand2idpR[contig_pos_strand]=(new_val)
				modcsv_contig2base[contig_pos_strand]=(info[3])

		else:
			if int(info[2]) == 0:
				args.output.write('chrm2'+'\t'+str(info[1])+'\t'+str(info[1])+'\t'+str(info[8])+'\t'+str('+')+'\n')
				contig_pos_strand = 'chrm2'+'_'+str(info[1])+'_+'
				pos_strand2idpR[contig_pos_strand]=(str(info[8]))
				modcsv_contig2base[contig_pos_strand]=(info[3])
			else:
				args.output.write('chrm2'+'\t'+str(info[1])+'\t'+str(info[1])+'\t'+str(info[8])+'\t'+str('-')+'\n')
				contig_pos_strand = 'chrm2'+'_'+str(info[1])+'_-'
				new_val = '-'+str(info[8])
				pos_strand2idpR[contig_pos_strand]=(new_val)
				modcsv_contig2base[contig_pos_strand]=(info[3])


####################################################
########### Motifs want to investigate #############
####################################################
motifs2investigate = []
for line in args.input:
	info = line.strip()
	motifs2investigate.append(info)


###############################
### go through motifs.gff #####
###############################
total_number_of_motif_of_interest = []
motif_of_interest_and_score_threshold_reached = []
predicted_motif_modification_threshold_reached = []
predicted_motif_modification_threshold_iden_thresh_reached = []
count_test = 0
#intially loop through motifs.gff and idenity all motifs of interest.
contig_position2info = {}
for line in args.motifs:
	if line.startswith('#'):
		pass 
	else:
		# predicted motifs 
		get_info = line.strip().split('\t')
		contig = get_info[0]
		score = get_info[5]
		modification = get_info[2]
		position = get_info[3]
		strand = get_info[6]
		motif = get_info[-1].split(';')
		wanted_motif = []
		coverage = []
		identification_qv = []
		IPDRatio = []
		for x in motif:
			if x.startswith('motif='):
				wanted_motif.append(x[6:])
			else:
				pass
			if x.startswith('identification'):
				identification_qv.append(x[17:])
			else:
				pass
			if x.startswith('coverage'):
				coverage.append(str(x[9:]))
			else:
				pass
			if x.startswith('IPDRatio='):
				IPDRatio.append(str(x[9:]))
			else:
				pass

	##### GATHERING INFORMATION FOR MOTIF OF INTEREST ########
		if len(wanted_motif) > 0:
			if ',' in wanted_motif[0]:
				if wanted_motif[0].split(',')[0] in motifs2investigate or wanted_motif[0].split(',')[1] in motifs2investigate:
					total_number_of_motif_of_interest.append(wanted_motif)
					#contig_position = str(contig)+'_'+str(position)+'_'+str(strand)
					if str(contig).startswith('unitig_80'):
						contig_position = 'chrm1_'+str(position)+'_'+str(strand)
					else:
						contig_position = 'chrm2_'+str(position)+'_'+str(strand)

					modification_base = modcsv_contig2base.get(contig_position)
					contig_position2info[contig_position]=(score,modification,position,strand,wanted_motif[0],coverage,identification_qv, IPDRatio, modification_base)
				else:
					pass
			else:
				if wanted_motif[0] in motifs2investigate:
					total_number_of_motif_of_interest.append(wanted_motif)
					#contig_position = str(contig)+'_'+str(position)+'_'+str(strand)
					if str(contig).startswith('unitig_80'):
						contig_position = 'chrm1_'+str(position)+'_'+str(strand)
					else:
						contig_position = 'chrm2_'+str(position)+'_'+str(strand)

					modification_base = modcsv_contig2base.get(contig_position)
					contig_position2info[contig_position]=(score,modification,position,strand,wanted_motif[0],coverage,identification_qv, IPDRatio, modification_base)
				else:
					pass

		# clear lists so unique for each line.
		wanted_motif[:]
		identification_qv[:]
		coverage[:]
		IPDRatio[:]



##############################################################
######### versions of motif i've identified as modified ######
##############################################################
contig2pos2strand2modified = {}
rm_head = 0
for line in args.candidates:
	rm_head += 1
	if rm_head == 1:
		pass
	else:
		info = line.strip().split('\t')
		if info[0].startswith('unitig_80'):
			contig2pos2strand = 'chrm1_'+str(info[5])+'_'+str(info[7])
			contig2pos2strand2modified[contig2pos2strand]=(str(info[9]))
		
		else:
			contig2pos2strand = 'chrm2_'+str(info[5])+'_'+str(info[7])
			contig2pos2strand2modified[contig2pos2strand]=(str(info[9]))


###### create
mot2be_able2create_positions = {'VAGGATCC':['-4','-3','-2','-1','0','1','2','3'],'BRGATCC':['-3','-2','-1','0','1','2','3'],'BRGATCTSV':['-3','-2','-1','0','1','2','3','4','5']}


########
if args.rect:
	args.rect.write('chrm'+'\t'+'start'+'\t'+'end'+'\t'+'value'+'\t'+'motif'+'\t'+'modified'+'\n')

if args.forward:
	args.forward.write('chrm'+'\t'+'start'+'\t'+'end'+'\t'+'value'+'\t'+'motif'+'\t'+'modified'+'\n')


# loop through each copy of our motifs of interest
for k, v in contig_position2info.iteritems():
	mod = contig2pos2strand2modified.get(k)
	for m in v[4].split(','):
		to_get_all_pos = mot2be_able2create_positions.get(m)
		if to_get_all_pos == None:
			print k, v, m
		else:
			#print k, v, m, to_get_all_pos
			all_positions = []
			for x in to_get_all_pos:
				if x.startswith('-'):
					just_num = x[1:]
					needed_pos = int(v[2]) - int(just_num)
					need_key = k.split('_')
					key = need_key[0]+'_'+str(needed_pos)+'_'+str(need_key[2])
					all_positions.append(key)
				else:
					needed_pos = int(v[2]) + int(x)
					need_key = k.split('_')
					key = need_key[0]+'_'+str(needed_pos)+'_'+str(need_key[2])
					all_positions.append(key)

			for w in all_positions:
				ratio = pos_strand2idpR.get(w)
				stuff = w.split('_')
				if mod == None:
					args.rect.write(str(stuff[0])+'\t'+str(stuff[1])+'\t'+str(stuff[1])+'\t'+str(ratio)+'\t'+str(m)+'\t'+'darkmagenta'+'\n')

				else:
					args.rect.write(str(stuff[0])+'\t'+str(stuff[1])+'\t'+str(stuff[1])+'\t'+str(ratio)+'\t'+str(m)+'\t'+'forestgreen'+'\n')


			# output file with just ipdR of modified adenines in motifs of interest
			if mod == None:
				pass

			else:
				# plot just the ipdRs of the adenines in the 
				if m == 'BRGATCC':
					stuff = k.split('_')
					ratio = pos_strand2idpR.get(k)
					args.forward.write(str(stuff[0])+'\t'+str(stuff[1])+'\t'+str(stuff[1])+'\t'+str(ratio)+'\t'+str(m)+'\t'+'darkmagenta'+'\n')

				else:
					if m == 'VAGGATCC':
						stuff = k.split('_')
						ratio = pos_strand2idpR.get(k)
						args.forward.write(str(stuff[0])+'\t'+str(stuff[1])+'\t'+str(stuff[1])+'\t'+str(ratio)+'\t'+str(m)+'\t'+'forestgreen'+'\n')

					else:
						stuff = k.split('_')
						ratio = pos_strand2idpR.get(k)
						args.forward.write(str(stuff[0])+'\t'+str(stuff[1])+'\t'+str(stuff[1])+'\t'+str(ratio)+'\t'+str(m)+'\t'+'red'+'\n')

				
		#print k, v[4], w, ratio

	# determine motif coordinates of bases represented by motif

	# get ipdR ratios from modifications.csv dictionary




