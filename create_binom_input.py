import argparse
parser=argparse.ArgumentParser()
parser.add_argument('-i','--input',type=argparse.FileType('r'),help='observed snps bp')
parser.add_argument('-e','--expected',type=argparse.FileType('r'),help='expected snps bp')
parser.add_argument('-o','--output',type=argparse.FileType('w'),help='output for r')
parser.add_argument('-s','--snps',type=str,help='type of snps of interest')
parser.add_argument('-l','--locus2prod',type=argparse.FileType('r'),help='locus product')
parser.add_argument('-r','--rna_locus',type=argparse.FileType('r'),help='rna locus tags')
args=parser.parse_args()


locus2expected = {}
for line in args.expected:
	info = line.strip().split('\t')
	locus2expected[info[0]]=(info[1])

locus2line = {}
for line in args.input:
	rm_crap = line.replace('"','')
	rm_crap = line.replace("'",'')
	info = rm_crap.strip().split('\t')
	if info[0] == args.snps:
		locus2line[info[1]]=(info)
	else:
		pass

locus2product = {}
for line in args.locus2prod:
	rm_crap = line.replace('"','')
	rm_crap = line.replace("'",'')
	info = rm_crap.strip().split('\t')
	locus2product[info[0]]=(info[1])

rna_loci = []
for line in args.rna_locus:
	info = line.strip().split(' ')
	rna_loci.append(info[0])


args.output.write('locus'+'\t'+'product'+'\t'+'observed'+'\t'+'snps_per_bp'+'\t'+'positions'+'\t'+'expected'+'\n')

for k, v in locus2product.iteritems():
	expec = locus2expected.get(k)
	observ = locus2line.get(k)
	if k in rna_loci:
		pass

	else:
		if observ == None:
			args.output.write(str(k)+'\t'+str(v)+'\t'+str(0)+'\t'+str(0)+'\t'+'NaN'+'\t'+str(expec)+'\n')
		else:
			observ.pop(0)
			observ.insert(len(observ), expec)
			count_out = 0
			for x in observ:
				count_out +=1 
				if count_out == len(observ):
					args.output.write(x+'\n')
				else:
					args.output.write(str(x)+'\t')

