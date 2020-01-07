import argparse
import random
from collections import Counter
from collections import defaultdict
parser=argparse.ArgumentParser()
parser.add_argument('-i','--input',type=argparse.FileType('r'),help='base per line domains marked')
parser.add_argument('-n','--number',type=int,help='number of SNPs to randomly introduce')
parser.add_argument('-p','--permutations',type=int,help='number of permutations to perform')
parser.add_argument('-d','--domain',type=str,help='domain of interest')
args=parser.parse_args()


def permutate():
    bases_hit = []
    for x in range(args.number):
        bases_hit.append(random.randint(1,len(bases)))

    return bases_hit
        

base2domain = {}
bases = []
count = 0
for line in args.input:
    count += 1
    bases.append(count)
    info = line.strip()
    base2domain[count]=(info)


domain2total_snps = defaultdict(list)

# run the permutations
for x in range(args.permutations):
    distribution = permutate()
    # loop through to get where snp fell:
    domain_hit = []
    for b in distribution:
        domain_hit.append(base2domain.get(b))

    # place domain distribution count in global variable for this permutation
    for k, v in Counter(domain_hit).iteritems():
        domain2total_snps[k].append(v)



# calculate the average #
for k, v in domain2total_snps.iteritems():
    average = float(sum(v))/float(args.permutations)
    print k, average
        
