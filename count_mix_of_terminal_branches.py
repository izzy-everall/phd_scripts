import argparse
import sys
from collections import Counter
from collections import defaultdict
parser=argparse.ArgumentParser()
parser.add_argument('-i','--input',type=argparse.FileType('r'),help='snps per gene output file from track_changes_across_phylogeny_v3.py')
parser.add_argument('-l','--lanes',type=argparse.FileType('r'),help='lanes of terminal branches to keep all others will be discounted')
parser.add_argument('-s','--snp_internal',action='store_true',help='do you want to keep terminal branches if SNP has been observed on internal node')
parser.add_argument('-o','--output',type=argparse.FileType('w'),help='output file name')
args = parser.parse_args()


terminal_branches2keep = []
if args.lanes:
	for line in args.lanes:
		lane = line.strip()
		terminal_branches2keep.append(lane)


terminal_branch_locus2snp_info = defaultdict(list)
final_snps2keep = defaultdict(list)
header = 0
# for every locus in each cluster
for line in args.input:
	header += 1
	if header == 1:
		if args.snp_internal:
			pass
		else:
			info = line.strip().split('\t')
			del info[13:19]
			info.insert(13, 'SNP info')

			out_count = 0
			for x in info:
				out_count += 1
				if out_count == len(info):
					args.output.write(str(x)+'\n')
				else:
					args.output.write(str(x)+'\t')

	else:
		info2keep = []
		info = line.strip().split('\t')
		snp_info = info[-1].split('__')
		total_num_snps = float(info[11])+float(info[12])+float(info[7])
		# which branches and which isolates and what effect did he snps in this locus occur in...
		for individ_snp in snp_info:
			breakdown = individ_snp.split('..')
			if args.lanes:
				# identify terminal branch lanes we want to keep
				for l in terminal_branches2keep:
					if l == breakdown[0].split('>')[-1]:
						info2keep.append(breakdown)

					else:
						pass
			# terminal branch snps will have list of length 4.
			if len(breakdown) == 4:
				if args.snp_internal:
					breakdown.append(info[0])
					print 'WHY', info[3], breakdown
					terminal_branch_locus2snp_info[info[3]].append(breakdown)
				else:
					pass
			else:
				info2keep.append(breakdown)


		
		# Does the SNP count info change from original
		if len(info2keep) == total_num_snps:
			#remove uncessary columns
			del info[13:19]
			if args.snp_internal:
				final_snps2keep[info[3]].append(info)

			else:
				# all SNPs on branch are of interest
				#either not on terminal branch or on terminal branch of interest and there's no change from original
				count_out = 0
				for o in info:
					count_out += 1
					if count_out == len(info):
						args.output.write(str(o)+'\n')
					else:
						args.output.write(str(o)+'\t')


		else:
			# if some of the SNPs in this gene are of interest
			if len(info2keep) >= 1:
				out_snps = []
				snp_type = []
				syn = 0
				nonsyn = 0
				intg = 0
				SNOPs = 0
				STOPs = 0
				STIPs = 0
				total_syn = 0
				total_nonsyn = 0
				# what type of SNPs remain
				for x in info2keep:
					snp_type.append(x[-1])
					#print info[1], info[3], x
					out = '..'.join(x)
					out_snps.append(out)

				for k, v in Counter(snp_type).iteritems():
					if k == 'Synonymous':
						syn = v
					else:
						if k == 'Nonsynonymous':
							nonsyn = v
						else:
							if k == 'Intergenic':
								intg = v
							else:
								if k.startswith('SNOP'):
									SNOPs = v
								else:
									if k.startswith('STOP'):
										STOPs = v
									else:
										if k.startswith('STIP'):
											STIPs = v
										else:
											pass

				new_total_syn = float(syn) + float(STOPs)
				new_total_nonsyn = float(nonsyn)+float(SNOPs)+float(STIPs)

				
				# reasign values
				info.pop(5)
				info.insert(5, syn)
				info.pop(6)
				info.insert(6, nonsyn)
				info.pop(7)
				info.insert(7, intg)
				info.pop(8)
				info.insert(8, SNOPs)
				info.pop(9)
				info.insert(9, STOPs)
				info.pop(10)
				info.insert(10, STIPs)
				info.pop(11)
				info.insert(11, new_total_syn)
				info.pop(12)
				info.insert(12, new_total_nonsyn)

				#remove unecessary columns 
				del info[13:20]

				#re-insert now correct info
				snp_stuff = '__'.join(out_snps)
				info.insert(13, snp_stuff)
				#print info[0], info[3], out_snps, snp_stuff

				if args.snp_internal:
					final_snps2keep[info[3]].append(info)

				else:
					count_out = 0
					for o in info:
						count_out += 1
						if count_out == len(info):
							args.output.write(str(o)+'\n')
						else:
							args.output.write(str(o)+'\t')

			# SNPs are not of interest --> terminal branch SNPs not within 100 days of sampling.
			else:
				pass



if args.snp_internal:
	args.output.write('locus'+'\t'+'product'+'\t'+'syn'+'\t'+'nonsyn'+'\t'+'intergenic'+'\t'+'SNOPs'+'\t'+'STOPs'+'\t'+'STIPs'+'\t'+'total_syn'+'\t'+'total_non_syn'+'\t'+'snp_info'+'\n')
	#final SNPs to keep contains all SNPs occurring on internal branches
	for l, s in final_snps2keep.iteritems():
		print l
		# extract terminal branch SNPs at this loci
		snps_on_t = terminal_branch_locus2snp_info.get(l)
		syn = []
		nonsyn = []
		intg = []
		SNOPs = []
		STOPs = []
		STIPs = []
		snp_info = []

		# for the SNPs on internal nodes add total numbers together.
		for snp_on_int_node in s:
			syn.append(float(snp_on_int_node[5]))
			nonsyn.append(float(snp_on_int_node[6]))
			intg.append(float(snp_on_int_node[7]))
			SNOPs.append(float(snp_on_int_node[8]))
			STOPs.append(float(snp_on_int_node[9]))
			STIPs.append(float(snp_on_int_node[10]))
			add_subsp_cluster = snp_on_int_node[-1]+'..'+snp_on_int_node[0]
			snp_info.append(add_subsp_cluster)


		if snps_on_t == None:
			print 'Just on internal?', l, s

		else:
			# add terminal branch SNPs
			for t_branch_snps in snps_on_t:
				if t_branch_snps[-2] == 'Synonymous':
					syn.append(1)
				else:
					if t_branch_snps[-2] == 'Nonsynonymous':
						nonsyn.append(1)

					else:
						if t_branch_snps[-2] == 'Intergenic':
							intg.append(1)
						else:
							if t_branch_snps[-2].startswith('SNOP'):
								SNOPs.append(1)
							else:
								if t_branch_snps[-2].startswith('STOP'):
									STOPs.append(1)
								else:
									if t_branch_snps[-2].startswith('STIP'):
										STIPs.append(1)
									else:
										print 'SNP TYPE NOT IDENTIFIED', t_branch_snps
										

				snp_stuff = '..'.join(t_branch_snps)
				snp_info.append(snp_stuff)

		total_syn = sum(syn) + sum(STOPs)
		total_nonsyn = sum(nonsyn) + sum(SNOPs) + sum(STIPs)
		all_snps = '__'.join(snp_info)
		args.output.write(str(l)+'\t'+str(s[0][4])+'\t'+str(sum(syn))+'\t'+str(sum(nonsyn))+'\t'+str(sum(intg))+'\t'+str(sum(SNOPs))+'\t'+str(sum(STOPs))+'\t'+str(sum(STIPs))+'\t'+str(total_syn)+'\t'+str(total_nonsyn)+'\t'+str(all_snps)+'\n')
