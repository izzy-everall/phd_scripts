# produces file with aa changes from list of summarise_snp.out files

import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio.Alphabet import IUPAC
from Bio.Data import CodonTable
import pandas as pd
import sys
from collections import Counter

parser = argparse.ArgumentParser()
parser.add_argument('-f','--fasta_file',type=str,help='fasta file name')
parser.add_argument('-l','--locus',type=argparse.FileType('r'),help='locus2orientation2start2end_positions')
parser.add_argument('-p','--product',type=argparse.FileType('r'),help='locus2product file')
parser.add_argument('-s','--snps',type=argparse.FileType('r'),help='snp2pos.txt',nargs='+')
parser.add_argument('-z','--per_pat_snps',action='store_true',help='write output file per patient')
parser.add_argument('-i','--intergenic',action='store_true',help='include intergenic snps')
parser.add_argument('-o','--output',type=argparse.FileType('w'),help='aa changes output file')
parser.add_argument('-c','--count',action='store_true',help='calculate intergenic region lengths and output in file')
args=parser.parse_args()
########

locus2start_end = {}
locus2orientation = {}
locus2product = {}

########
standard_table = CodonTable.unambiguous_dna_by_name["Standard"]
base_pair_dict = {'A':'T','T':'A','C':'G','G':'C'}
########
def ribosomal_change(start,end,snp_pos,orientation,alt):
    aa_change_info = []
    gene = chrs[seq_name[0]][int(start)-1:int(end)]
    change_position = int(snp_pos)-(int(start))
    get_snp_change = gene.tomutable()
    # zero base error correction
    art_snp_pos = change_position + 1
    final_snp_change = get_snp_change[change_position].upper()+str(art_snp_pos)+alt
    aa_change_info.append(final_snp_change)
    aa_change_info.append('Nonsynonymous')
    #print get_snp_change[change_position], change_position, alt
    return aa_change_info

### get aa change ###
def aa_change(start,end,snp_pos,orientation,alt):

    aa_change_info = []
    if orientation == 'forward':
        gene = chrs[seq_name[0]][int(start)-1:int(end)]

        # where snp occurs
        change_position = int(snp_pos)-(int(start))
        aa_positions = change_position/3

        # reference AA
        actual_prot_seq = gene.translate()
        find_aa = actual_prot_seq.tomutable()

        # alterantive aa
        insert_alt_snp = gene.tomutable()
        insert_alt_snp[change_position]=alt
        alt_prot_seq = insert_alt_snp.toseq().translate()

        # because python is 0 based and artemis/sequence viewers start at 1 need to add one to amino acid pos
        final_aa_change_pos =  aa_positions + 1
        final_aa_change = find_aa[aa_positions]+str(final_aa_change_pos)+alt_prot_seq[aa_positions]

        if find_aa[aa_positions] == alt_prot_seq[aa_positions]:
            effect = 'Synonymous'
        else:
            if find_aa[aa_positions] == '*' and alt_prot_seq[aa_positions] != '*':
                effect = 'STIP (stop codon to non-stop codon)'
            else:
                if alt_prot_seq[aa_positions] == '*' and find_aa[aa_positions] != '*':
                    effect = 'SNOP (non-stop codon to stop codon)'
                else:
                    effect = 'Nonsynonymous'

        aa_change_info.append(final_aa_change)
        aa_change_info.append(effect)
        
        
    else:
        # get gene and reverse complement
        gene = chrs[seq_name[0]][int(start)-1:int(end)]
        rev_comp = gene.reverse_complement()
        
        # where snp occurs
        change_position = int(end)-(int(snp_pos))
        aa_positions = change_position/3

        # reference AA
        actual_prot_seq = rev_comp.translate()
        find_aa = actual_prot_seq.tomutable()

        # alternative 
        insert_alt_snp = rev_comp.tomutable()
        rev_alt_pos = base_pair_dict.get(alt)

        insert_alt_snp[change_position]=rev_alt_pos
        alt_prot_seq = insert_alt_snp.toseq().translate()

        # because python is 0 based and artemis sequence viewers start at 1 need to add one to amino acid pos
        final_aa_change_pos =  aa_positions + 1
        final_aa_change = find_aa[aa_positions]+str(final_aa_change_pos)+alt_prot_seq[aa_positions]

        if find_aa[aa_positions] == alt_prot_seq[aa_positions]:
            effect = 'Synonymous'
        else:
            if find_aa[aa_positions] == '*' and alt_prot_seq[aa_positions] != '*':
                effect = 'STIP (stop codon to non-stop codon)'
            else:
                if alt_prot_seq[aa_positions] == '*' and find_aa[aa_positions] != '*':
                    effect = 'SNOP (non-stop codon to stop codon)'
                else:
                    effect = 'Nonsynonymous'

        aa_change_info.append(final_aa_change)
        aa_change_info.append(effect)


    return aa_change_info


#### function to determine locus ###
def which_locus(snp_pos,ref_base,alt_base):
    locus_found = []
    for k, v in locus2start_end.iteritems():
        if int(snp_pos) >= int(v[0]) and int(snp_pos) <= int(v[1]):
            orientation = locus2orientation.get(k)
            product_info = locus2product.get(k)
            if k.startswith('MAB_r') or '16S' in product_info or '23S' in product_info or '5S' in product_info or 'tRNA' in product_info:
                get_ribosomal_change=ribosomal_change(v[0],v[1],snp_pos,orientation,alt_base)
                locus_found.append(k)
                locus_found.append(get_ribosomal_change[0])
                locus_found.append(get_ribosomal_change[1])
                locus_found.append(product_info)
            else:
                get_aa_change = aa_change(v[0],v[1],snp_pos,orientation,alt_base)
                locus_found.append(k)
                locus_found.append(get_aa_change[0])
                locus_found.append(get_aa_change[1])
                locus_found.append(product_info)
        else:
            pass

    if args.intergenic:
        if len(locus_found) == 0:
            for k, v in intergenic_region2coords.iteritems():
                if int(snp_pos) >= int(v[0]) and int(snp_pos) <= int(v[1]):
                    keys = k.split('..')
                    start_product = locus2product.get(keys[0])
                    end_product = locus2product.get(keys[1])
                    # append name of intergenic region
                    locus_found.append(k)
                    # join ref and alt base
                    ref2pos2alt = str(ref_base)+str(snp_pos)+str(alt_base)
                    locus_found.append(ref2pos2alt)
                    # type of effect
                    locus_found.append('intergenic')
                    # products of flanking genes
                    start_prod2end_prod = start_product+'..'+end_product
                    locus_found.append(start_prod2end_prod)
                    
                else:
                    pass
        else:
            pass

    return locus_found

#### read in fasta file ####
chrs = {}
seq_name = []
all_coords = []
for seq in SeqIO.parse(open(args.fasta_file),'fasta'):
    seq_name.append(seq.id)
    chrs[seq.id]=seq.seq

### read in gene2orientation2positions ###
all_lines = []
all_lines.append(['locus','orientation','start','end'])
for line in args.locus:
    info = line.strip().split('\t')
    start_end = (info[2],info[3])
    locus2start_end[info[0]]=(start_end)
    locus2orientation[info[0]]=(info[1])
    if args.intergenic:
        all_lines.append(info)
        all_coords.append(int(info[2]))
        all_coords.append(int(info[3]))

### to investigate intergenic snsp ###
intergenic_region2coords = {}
if args.intergenic:
    if args.count:
        count_file=open('intergenic_region2length.txt','w')
    l2O2S2E_df = pd.DataFrame(all_lines[1:],columns=all_lines[0])
    all_coords.pop(0)
    all_coords.pop(-1)
    end2next_start_coords = [all_coords[i:i+2] for i in range(0, len(all_coords), 2)]
    # extracting intergenic regions
    for x in end2next_start_coords:
        # overlapping genes no intergenic region 
        if int(x[1]) <= int(x[0]):
            pass
        else: 
            start_int = l2O2S2E_df.loc[l2O2S2E_df['end']==str(x[0])].values.tolist()
            end_int = l2O2S2E_df.loc[l2O2S2E_df['start']==str(x[1])].values.tolist()
            start_int2end_int = start_int[0][0]+'..'+end_int[0][0]
            start_coord2end_coord = (x[0],x[1])
            intergenic_region2coords[start_int2end_int]=(start_coord2end_coord)
            if args.count:
                length = int(x[1])-int(x[0])+1
                count_file.write(start_int2end_int+'\t'+str(length)+'\n')



for line in args.product:
    info = line.strip().split('\t')
    locus2product[info[0]]=(info[1])


#### read snp files ####
for f in args.snps:
    count_line = 0
    isolate = str(f).split(' ')
    isolate = isolate[2].split('/')
    isolate = isolate[-1][:-11]

    # if yout want snps per pat
    if args.per_pat_snps:
        file_name = isolate+'_pat_snps'+".txt"
        file2write = open(file_name,'w')
    
    for line in f:
        count_line += 1
        if count_line == 1:
            info = line.strip().split('\t')
            number_isols = info[5:]
        else:

            info = line.strip().split('\t')
            try:
                #print info
                # all snps called and structural identified
                if int(info[4]) == len(number_isols):
                    pass

                else:
                    # remove snps that are more than likely structural
                    alt = info[3]
                    snps = info[5:]
                    still_structural = Counter(snps)

                    num_alt = still_structural.get(alt)
                    num_not_called = still_structural.get('N')

                    try:
                        if int(num_alt)+int(num_not_called) == len(number_isols):
                            #removes cases where a SNP has been called but all other sites in this position
                            #are an N put most likely a snp so would make it structural
                            pass
                        else:
                            snp_type = which_locus(info[0],info[2],info[3])
                            if len(snp_type) == 0:
                                pass
                            else:
                                 #pass
                                args.output.write(isolate+'\t'+info[0]+'\t'+snp_type[0]+'\t'+info[2]+'\t'+info[3]+'\t'+snp_type[1]+'\t'+snp_type[2]+'\t'+snp_type[3]+'\n')
                                if args.per_pat_snps:
                                     file2write.write(isolate+'\t'+info[0]+'\t'+snp_type[0]+'\t'+info[2]+'\t'+info[3]+'\t'+snp_type[1]+'\t'+snp_type[2]+'\t'+snp_type[3]+'\n')
                                 #print "within pat snps", isolate, info[0], snp_type[0], info[2], info[3], snp_type[1], snp_type[2], snp_type[3]

                    except TypeError:
                        snp_type = which_locus(info[0],info[2],info[3])
                        if len(snp_type) == 0:
                            pass
                        else:
                            #pass
                            args.output.write(isolate+'\t'+info[0]+'\t'+snp_type[0]+'\t'+info[2]+'\t'+info[3]+'\t'+snp_type[1]+'\t'+snp_type[2]+'\t'+snp_type[3]+'\n')
                            if args.per_pat_snps:
                                file2write.write(isolate+'\t'+info[0]+'\t'+snp_type[0]+'\t'+info[2]+'\t'+info[3]+'\t'+snp_type[1]+'\t'+snp_type[2]+'\t'+snp_type[3]+'\n')
                            #print "within pat snps", isolate, info[0], snp_type[0], info[2], info[3], snp_type[1], snp_type[2], snp_type[3]



            # catching heterogeneic sites two alternative bases called
            except ValueError:
                heterogenic_site = info[4].split(',')
                all_sites_at_pos = info[5:]
                count_bases = Counter(all_sites_at_pos)
                het_sites = []
                for k, v in count_bases.iteritems():
                    het_sites.append(int(v))
                alt = info[3].split(',')
                # structural hetero sites
                if sum(het_sites) == len(number_isols):
                    pass
                else:
                    snp_type_1 = which_locus(info[0],info[2],alt[0])
                    snp_type_2 = which_locus(info[0],info[2],alt[1])

                    if len(snp_type_1) == 0:
                        pass

                    else:
                        args.output.write(isolate+'\t'+info[0]+'\t'+snp_type_1[0]+'\t'+info[2]+'\t'+alt[0]+'\t'+snp_type_1[1]+'\t'+snp_type_1[2]+'\t'+snp_type_1[3]+'\n')
                        if args.per_pat_snps:
                            file2write.write(isolate+'\t'+info[0]+'\t'+snp_type_1[0]+'\t'+info[2]+'\t'+alt[0]+'\t'+snp_type_1[1]+'\t'+snp_type_1[2]+'\t'+snp_type_1[3]+'\n')

                    #print "within pat snps", isolate, info[0], snp_type_1[0], info[2], alt[0], snp_type_1[1], snp_type_1[2], snp_type_1[3]

                    if len(snp_type_2) == 0:
                        pass
                    else:
                        args.output.write(isolate+'\t'+info[0]+'\t'+snp_type_2[0]+'\t'+info[2]+'\t'+alt[1]+'\t'+snp_type_2[1]+'\t'+snp_type_2[2]+'\t'+snp_type_2[3]+'\n')

                        if args.per_pat_snps:
                            file2write.write(isolate+'\t'+info[0]+'\t'+snp_type_2[0]+'\t'+info[2]+'\t'+alt[1]+'\t'+snp_type_2[1]+'\t'+snp_type_2[2]+'\t'+snp_type_2[3]+'\n')
                    #print "within pat snps", isolate, info[0], snp_type_2[0], info[2], alt[1], snp_type_2[1], snp_type_2[2], snp_type_2[3]

