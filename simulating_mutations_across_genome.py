import argparse
import pandas as pd
import sys
from collections import Counter
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio.Alphabet import IUPAC
from Bio.Data import CodonTable
from collections import defaultdict


parser=argparse.ArgumentParser()
parser.add_argument('-a','--annotation_info',type=argparse.FileType('r'),help='gene2locus2orientation2start2end file')
parser.add_argument('-f','--fasta',type=argparse.FileType('r'),help='fasta file of reference')
parser.add_argument('-s','--synonymous',type=int,help='number of synonymous snps observed')
parser.add_argument('-c','--cds',type=int,help='number of CDS bases in reference')
parser.add_argument('-g','--genes',action='store_true',help='calcuate expected values per gene')
parser.add_argument('-z','--syn_snp_file',type=argparse.FileType('r'),help='tab delimited file with gene to total number of synonymous SNPs')
parser.add_argument('-o','--output',type=argparse.FileType('w'),help='output file name')
args=parser.parse_args()


### calculating the synonymous SNP frequency (total number of synonymous SNPs/number of cd bases)
def synonymous_mutation_rate(obs_syn,cds_bases):
    psn = float(obs_syn)/float(cds_bases)
    return psn

def non_syn_mutation_rate(syn_mut_rate,WGS_NS, WGS_SN, gene_len):
    ### calcuate the whole gneome non-synonymous background rate ###
    NS = sum(WGS_NS)
    SN = sum(WGS_SN)
    R = float(NS)/float(SN)
    pns = float(syn_mut_rate)*float(R)
    print 'per base of gene', pns
    if args.genes:
        pns_per_gene = float(pns)*float(gene_len)
        print 'per gene', pns_per_gene
        return pns_per_gene
    
    else:
        return pns

### reading in fasta ###
chrs = {}
seq_name = []
all_coords = []
def read_fasta(fasta_file):
    for seq in SeqIO.parse(fasta_file,'fasta'):
        seq_name.append(seq.id)
        chrs[seq.id]=seq.seq

### reading in annotation information ###
locus2start_end={}
locus2orientation = {}
locus2gene_length = {}
locus2all_info = {}
def organise_annotation_data(annotation):
    for line in annotation:
        info = line.strip().split('\t')
        start_end = (info[2],info[3])
        locus2start_end[info[0]]=(start_end)
        locus2orientation[info[0]]=(info[1])
        locus2gene_length[info[0]]=(info[-1])
        all_info = [info[1],info[2],info[3],info[-1]]
        locus2all_info[info[0]]=(all_info)

### reading in synonymous SNP info ###
locus2syn_snp = {}
def read_metadata(metadata):
    for line in metadata:
        info = line.strip().split('\t')
        locus2syn_snp[info[0]]=(info[-1])

### determine aa change ###
def aa_change(codon_seq):
    str_codon = ''.join(codon_seq)
    str_codon = Seq(str_codon)
    ref_aa = str_codon.translate()
    return ref_aa

### expecting number of non_syn sites and syn sites per codon ###
def expected_number_of_non_syn_and_syn_sites_per_codon(codon):
    bases = ['A','T','C','G']
    ref_aa = aa_change(codon)
    ref_codon = {'ref':codon}
    nc2sc = []
    changes = []
    for n, i in enumerate(codon):
        for b in bases:
            if i == b:
                pass
            else:
                alt_codon_seq = ref_codon.get('ref')
                alt_codon_seq[n] = b
                new_aa = aa_change(alt_codon_seq)
                if new_aa == ref_aa:
                    changes.append(0)
                else:
                    val = float(1)/float(3)
                    changes.append(val)
                alt_codon_seq[n]=i

    nc = sum(changes)
    sc = float(3) - float(nc)
    nc2sc.append(nc)
    nc2sc.append(sc)
    return nc2sc

### split gene into codons ###
def split_into_codons(gene):
    return [ gene[start:start+3] for start in range(0, len(gene), 3) ]

### calculate the ratio of non-synonymous to synonymous sites across the genome ###
def ratio_of_nonsyn_to_syn_sites(locus,start,end,orientation,length):
    
    expected_syn2expected_non_syn = []
    if locus.startswith('MAB_r') or locus.startswith('MAB_t'):
        pass
        #gene_length = locus2gene_length.get(locus)
        #gene_expected_NS_vs_expected_S = float(3)*float(gene_length)
        #expected_syn2expected_non_syn.append([gene_expected_NS_vs_expexted_S,0])

    else:
        if orientation == 'forward':
            cds = chrs[seq_name[0]][int(start)-1:int(end)]
            
            ### check if complete gene ###
            prot_seq = cds.translate()
            prot_seq = prot_seq.tomutable()
            if prot_seq[-1]=='*':
                pass
            else:
                print 'ERROR: incomplete coding sequence', locus
                sys.exit()

            cds = str(cds).upper()
            cds = list(cds)
            codons = split_into_codons(cds)
            for c in codons:
                expec_non_syn_syn_sites = expected_number_of_non_syn_and_syn_sites_per_codon(c)
                expected_syn2expected_non_syn.append(expec_non_syn_syn_sites)
            

        else:
            cds = chrs[seq_name[0]][int(start)-1:int(end)]
            rev_comp = cds.reverse_complement()
            
            ### check if complete gene ###
            prot_seq = rev_comp.translate()
            prot_seq = prot_seq.tomutable()
            if prot_seq[-1] == '*':
                pass
            else:
                print 'ERROR: incomplete coding sequence', locus
                sys.exit()

            rev_comp = str(rev_comp).upper()
            rev_comp = list(rev_comp)
            codons = split_into_codons(rev_comp)
            for c in codons:
                expec_non_syn_syn_sites = expected_number_of_non_syn_and_syn_sites_per_codon(c)
                #print locus, expec_non_syn_syn_sites
                expected_syn2expected_non_syn.append(expec_non_syn_syn_sites)
                
    all_NS_per_gene = []
    all_SN_per_gene = []
    total_NS_SN_ratio = []
    for x in expected_syn2expected_non_syn:
        all_NS_per_gene.append(x[0])
        all_SN_per_gene.append(x[1])

    total_NS = sum(all_NS_per_gene)
    total_SN = sum(all_SN_per_gene)
    total_NS_SN_ratio.append(total_NS)
    total_NS_SN_ratio.append(total_SN)


    print locus, total_NS, total_SN

    return total_NS_SN_ratio

# calculating the expected non-synonymous mutation frequency across the genome and per gene
### calling functions ###
organise_annotation_data(args.annotation_info)
read_fasta(args.fasta)
read_metadata(args.syn_snp_file)

### for every condon in every locus calculated the expected number NS SN sites ###
locus2expected_NS_SN_sites = {}
for k, v in locus2all_info.iteritems():
    expected_NS_SN_sites = ratio_of_nonsyn_to_syn_sites(k,v[1],v[2],v[0],v[-1])
    locus2expected_NS_SN_sites[k]=(expected_NS_SN_sites)

### to calculate the whole genome non-synonymous background rate ###
WGS_SN_frequency = synonymous_mutation_rate(args.synonymous,args.cds)
WGS_NS = []
WGS_SN = []
for k, v in locus2expected_NS_SN_sites.iteritems():
    WGS_NS.append(v[0])
    WGS_SN.append(v[1])

WGS_NS_frequency = non_syn_mutation_rate(WGS_SN_frequency,WGS_NS,WGS_SN,args.cds)
print "overall non_synonymous background frequency:", WGS_NS_frequency

### to calculate the per gene non-synonymous background rate ###
if args.genes:
    locus2per_gene_background_NS_frequency = {}
    for k, v in locus2expected_NS_SN_sites.iteritems():
        num_syn_snps = locus2syn_snp.get(k)
        gene_length = locus2gene_length.get(k)
        if num_syn_snps == None:
            if k.startswith('MAB_r') or k.startswith('MAB_t'):
                pass
            else:
            ### no observed synonymous SNPs therefore take global synonymous SNP rate as observed number of syn SNPs and multiply by ratio ##
                per_gene_background_non_syn_rate = non_syn_mutation_rate(WGS_SN_frequency,[v[0]],[v[1]],gene_length)
                locus2per_gene_background_NS_frequency[k]=(per_gene_background_non_syn_rate)
        
        else:
            per_gene_background_syn_rate = synonymous_mutation_rate(num_syn_snps,gene_length) 
            per_gene_background_non_syn_rate = non_syn_mutation_rate(per_gene_background_syn_rate,[v[0]],[v[1]],gene_length)
            locus2per_gene_background_NS_frequency[k]=(per_gene_background_non_syn_rate)

    
    for k, v in locus2per_gene_background_NS_frequency.iteritems():
        args.output.write(str(k)+'\t'+str(v)+'\n')




    








#THE Simple approach (LOL):
#does not attempt to correct for mutation biases, all non-synonymous mutations were assumed to occur at the same frequency
#frequency (non-synonymous mutation frequency) derived by multiplying the genome-wide non-synonymous/synonymous ration by the synonymous mutation frequency.
#for each gene, identify the probability of its non-syn mutation count given the expected non-synonymous mutation frequency and its length.  
