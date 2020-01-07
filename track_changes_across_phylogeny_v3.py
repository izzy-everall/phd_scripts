### Tracking changes on particular branches july 2017 ###
import argparse
from collections import defaultdict
import sys
import pandas as pd
from collections import Counter
### arguments ###
parser=argparse.ArgumentParser()
parser.add_argument('-i','--input',type=argparse.FileType('r'),help='reconstruct snps on tree .tab file')
parser.add_argument('-l','--lane2cluster',type=argparse.FileType('r'),help='lane2cluster file')
parser.add_argument('-s','--stats',type=argparse.FileType('w'), help='statistics output file')
parser.add_argument('-b','--branches_to',action='store_true', help='output snps on branches leading to cluster')
parser.add_argument('-g','--locus_file',type=argparse.FileType('r'),help='locus2orientation2start2end_positions2length')
parser.add_argument('-p','--product',type=argparse.FileType('r'),help='locus2product file')
parser.add_argument('-r','--recombination_window_size',type=int,help='window size to detect possible recombinant snps')
parser.add_argument('-t','--terminal',action='store_true',help='DO NOT INCLUDE terminal branches which may contain within host evolution')
parser.add_argument('-o','--output',type=argparse.FileType('w'),help='output file')
parser.add_argument('-c','--concat_output',type=argparse.FileType('w'),help='concatenated output')
parser.add_argument('-d','--recombination_data',type=argparse.FileType('w'),help='recombinant snps')
args=parser.parse_args()
########################################
#     remove recombinant SNPs          #
########################################
def detect_recombination(cluster, branch2snps):
    all_snp_positions = []
    recombinant_snps = []
    for snp in branch2snps:
        all_snp_positions.append(int(snp[0]))

    ordered_snps = sorted(all_snp_positions, key=int)
    for i in range(len(ordered_snps)-1):
        try:
            if (ordered_snps[i+2] - ordered_snps[i]) <= args.recombination_window_size:
                recombinant_snps.append(ordered_snps[i])
                recombinant_snps.append(ordered_snps[i+1])
                recombinant_snps.append(ordered_snps[i+2])
            else:
                pass
        except IndexError:
            pass

    if len(recombinant_snps) == 0:
        return branch2snps
    else:
        to_keep_snps = []
        for keep_snp in branch2snps:
            if int(keep_snp[0]) in recombinant_snps:
                print keep_snp
                isols = '..'.join(keep_snp[-1])
                get_locus = which_locus(keep_snp[0])
                args.recombination_data.write(str(keep_snp[2])+'\t'+str(keep_snp[0])+'\t'+str(keep_snp[5])+'\t'+str(get_locus[0])+'\t'+str(get_locus[1])+'\t'+str(isols)+'\n')
            else:
                to_keep_snps.append(keep_snp)

        return to_keep_snps

#########################################
#                get locus             #
########################################
def which_locus(snp_pos):
    locus_found = []
    for k, v in locus2start_end.iteritems():
        if int(snp_pos) >= int(v[0]) and int(snp_pos) <= int(v[1]):
            orientation = locus2orientation.get(k)
            product_info = locus2product.get(k)
            #print snp_pos, orientation, int(v[0]), int(v[1])
            if orientation == 'forward':
                gene_snp_pos = int(snp_pos) - int(v[0])
                #print k, snp_pos, gene_snp_pos, orientation
            else:
                gene_snp_pos = int(v[1]) - int(snp_pos)
                #print k, snp_pos, gene_snp_pos, orientation

            locus_found.append(k)
            locus_found.append(product_info)
            locus_found.append(gene_snp_pos)
        else:
            pass

    
    if len(locus_found) == 0:
        for k, v in intergenic_region2coords.iteritems():
            if int(snp_pos) >= int(v[0]) and int(snp_pos) <= int(v[1]):
                keys = k.split('..')
                start_product = locus2product.get(keys[0])
                end_product = locus2product.get(keys[1])
                # append name of intergenic region
                locus_found.append(k)
                # products of flanking genes
                start_prod2end_prod = start_product+'..'+end_product
                locus2product[k]=(start_prod2end_prod)
                locus_found.append(start_prod2end_prod)
                    
            else:
                pass
    else:
        pass

    return locus_found

############################################
# read in locus2orientation2start2end file #
############################################
locus2start_end = {}
locus2orientation = {}
locus2length = {}
intergenic_region2coords = {}
locus2product = {}
all_lines = []
all_coords = []
all_lines.append(['locus','orientation','start','end','length'])
for line in args.locus_file:
    info = line.strip().split('\t')
    start_end = (info[2],info[3])
    locus2start_end[info[0]]=(start_end)
    locus2orientation[info[0]]=(info[1])
    if len(info) == 6:
        info.pop(-1)
    all_lines.append(info)
    all_coords.append(int(info[2]))
    all_coords.append(int(info[3]))
    locus2length[info[0]]=(info[-1])
###########################################
#   to get intergenic coords              #
###########################################
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
        int_length = (float(x[1])-float(x[0]))+1
        locus2length[start_int2end_int]=(int_length)
        intergenic_region2coords[start_int2end_int]=(start_coord2end_coord)
#########################################
#      read in locus2product file       #
#########################################
for line in args.product:
    info = line.strip().split('\t')
    locus2product[info[0]]=(info[1])
################################
lane2cluster = {}
cluster2lanes = defaultdict(list)
### read in lane2cluster file ###
header = 0
for line in args.lane2cluster:
    header += 1
    if header == 1:
        pass
    else:
        info = line.strip().split('\t')
        if len(info) == 1:
            info = line.strip().split(',')
            lane2cluster[info[0]]=(info[2])
            cluster2lanes[info[2]].append(info[0])
        else:
            # discounting isolates that are clustered from one CF centre.
            lane2cluster[info[0]]=(info[2])
            cluster2lanes[info[2]].append(info[0])

##################################
pos2strand2node2snp2homoplasy2codon_type2taxa = defaultdict(list)
####### read in tab file #########
cds_block = False
snp_block = ['NA']*7
for line in args.input:
    # determine block
    ### cds blocks ###
    if line.startswith('FT   CDS'):
        position = line.strip().split(' ')[-1]
        snp_block.pop(0)
        snp_block.insert(0,position)
        cds_block = True
    else:
        ## snp blocks ##
        if cds_block == False:
            if line.startswith('FT                   /taxa='):
                taxa = line.strip().replace(' ','')
                taxa = taxa.replace('"','')
                taxa = taxa.split('=')[-1]
                taxa = taxa.split(',')
                snp_block.pop(-1)
                snp_block.insert(len(snp_block), taxa)
                #print snp_block
                pos2strand2node2snp2homoplasy2codon_type2taxa[snp_block[0]].append(snp_block)
                snp_block = ['NA']*7
    
            else:
                if line.startswith('FT   SNP'):
                    pos = line.strip().split(' ')[-1]
                    snp_block.pop(0)
                    snp_block.insert(0, pos)
                else:
                    if line.startswith('FT                   /strand'):
                        strand = line.strip().split('="')[-1][:-1]
                        snp_block.pop(1)
                        snp_block.insert(1,strand)

                    elif line.startswith('FT                   /node='):
                        node = line.strip().replace('"','')
                        node = node.strip().split('=')[-1]
                        snp_block.pop(2)
                        snp_block.insert(2,node)

                    elif line.startswith('FT                   /SNP='):
                        snp = line.strip().replace('"','')
                        snp = snp.split('=')[-1]
                        snp_block.pop(3)
                        snp_block.insert(3, snp)

                    elif line.startswith('FT                   /homoplasy='):
                        homoplasy = line.strip().replace('"','')
                        homoplasy = homoplasy.split('=')[-1]
                        homoplasy = homoplasy.split(',')
                        snp_block.pop(4)
                        snp_block.insert(4,homoplasy)

                    elif line.startswith('FT                   /codon_type'):
                        codon_type = line.strip().replace('"','')
                        codon_type = codon_type.split('=')[-1]
                        snp_block.pop(5)
                        snp_block.insert(5,codon_type)

                    else:
                        pass
                        
        
        ### complete cds block ###
        else:
            if line.startswith('FT   CDS'):
                position = line.strip().split(' ')[-1]
                snp_block.pop(0)
                snp_block.insert(0, position)
            elif line.startswith('FT                   /primary_name='):
                name = line.strip().replace('"','')
                name=name.split('=')[-1]
                snp_block.pop(3)
                snp_block.insert(3,name)

            elif line.startswith('FT                   /node='):
                node = line.strip().replace('"','')
                node = node.strip().split('=')[-1]
                taxa = node.split('->')[-1]
                snp_block.pop(2)
                snp_block.insert(2,node)
                snp_block.pop(-1)
                snp_block.insert(len(snp_block),taxa)

            elif line.startswith('FT                   /colour'):
                pos2strand2node2snp2homoplasy2codon_type2taxa[snp_block[0]].append(snp_block)
                snp_block = ['NA']*7
                cds_block=False

#####################################################################
##### write total number of unique SNPs and total SNPs to file ######
#####################################################################
args.stats.write('total number of unique changes (SNPs and gene pres abs):'+'\t'+str(len(pos2strand2node2snp2homoplasy2codon_type2taxa))+'\n')
############ write total number of SNPs to file #####################
total_number_of_snps = []
for k, v in pos2strand2node2snp2homoplasy2codon_type2taxa.iteritems():
    total_number_of_snps.append(len(v))
args.stats.write('total number of changes:'+'\t'+str(sum(total_number_of_snps))+'\n')
######################################################################
#     identify SNPs that occurred on branches within each cluster    #
######################################################################
cluster2node2snps = defaultdict(dict)
### loop through clusters ###
for k, v in cluster2lanes.iteritems():
    # create dictionary of per branch snps #
    node2snps = defaultdict(list)

    ### uclustered isolates ###
    if k == 'unclustered':
        pass
        #for pos, info in pos2strand2node2snp2homoplasy2codon_type2taxa.iteritems():
        #    for a_snp_at_pos in info:
                # Evolution changes only associated with unclustered isolates #
        #        if set(a_snp_at_pos[-1]).issubset(set(v)):
        #            node2snps[a_snp_at_pos[2]].append(a_snp_at_pos)
        #        else:
        #            pass
    else:
        ### for each SNP how many SNPs occur on each brach in each cluster ###
        for pos, info in pos2strand2node2snp2homoplasy2codon_type2taxa.iteritems():
            for a_snp_at_pos in info:
                # snps or changes only involved in the evolution of a cluster
                if set(a_snp_at_pos[-1]).issubset(set(v)):
                    # remove snps on branches leading to clusters #
                    if len(a_snp_at_pos[-1]) == len(v) and args.branches_to:
                            print 'you want the snps on the branches to the clusters'
                    else:
                        if len(a_snp_at_pos[-1]) == len(v):
                            pass
                             
                        else:
                            node2snps[a_snp_at_pos[2]].append(a_snp_at_pos)
                else:
                    pass
    
    for node, snps in node2snps.iteritems():
        cluster2node2snps[k][node]=(snps)

######################################################################
#                  remove recombinant snps                           #
######################################################################
args.output.write('cluster'+'\t'+'total_number_of_branches_in_cluster'+'\t'+'number_of_isols_in_cluster'+'\t'+'locus'+'\t'+'product'+'\t'+'synonymous'+'\t'+'non_synonymous'+'\t'+'intergenic'+'\t'+'SNOPs'+'\t'+'STOPs'+'\t'+'STIPs'+'\t'+'total_syn'+'\t'+'total_non_syn'+'\t'+'pos_non_syn_snps'+'\t'+'pos_syn_snps'+'\t'+'pos_intergenic_snps'+'\t'+'branches_snps_occurred_on'+'\t'+'isols'+'\t'+'genome_positions'+'\n')

locus2non_syn = defaultdict(list)
locus2syn = defaultdict(list)
locus2non_syn_positions = defaultdict(list)
locus2syn_positions = defaultdict(list)
locus2int = defaultdict(list)
locus2int_positions = defaultdict(list)
for cluster, node in cluster2node2snps.iteritems():
    print 'Analysing within cluster snps for cluster', cluster

    num_in_cluster = cluster2lanes.get(cluster)
    within_cluster_locus2snp = defaultdict(list)
    within_cluster_locus2product = {}
    within_cluster_locus2position2effect = defaultdict(list)
    within_cluster_locus2isol_with_snp = defaultdict(list)
    within_cluster_locus2branch = defaultdict(list)
    within_cluster_locus2genome_position = defaultdict(list)
    within_cluster_locus2original_snp_positions = defaultdict(list)
    within_cluster_locus_node2pos2isol = defaultdict(list)
    syn = defaultdict(list)
    non_syn = defaultdict(list)
    inter = defaultdict(list)
    SNOP = defaultdict(list)
    STIP = defaultdict(list)
    STOP = defaultdict(list)

    # for each branch within a cluster
    for k, v in node.iteritems():
        if args.recombination_window_size:
            # remove recombinant snps
            recombination_removed_from_branch = detect_recombination(cluster, v)
            if len(recombination_removed_from_branch) == 0:
                print 'All SNPs on branch removed due to recombination', cluster, k
            
            else:
                # to include or not include terminal branches
                if len(recombination_removed_from_branch[0][-1]) == 1 and args.terminal:
                    print 'Discounting snps as on terminal branch', k
                else:
                    # determine which locus or intergenic region SNPs fall in #
                    for snp in recombination_removed_from_branch:
                        locus_identified = which_locus(snp[0])
                        if len(locus_identified) == 0:
                            locus_identified.append('end_chrm_intergenic')
                            locus_identified.append('NA')
                        within_cluster_locus2snp[locus_identified[0]].append(snp[5])
                        within_cluster_locus2product[locus_identified[0]]=(locus_identified[1])
                        within_cluster_locus2branch[locus_identified[0]].append(snp[2])
                        needed_pats = '..'.join(snp[6])
                        # original line
                        all_stuff  = snp[2]+'..'+snp[0]+'..'+needed_pats
                        within_cluster_locus_node2pos2isol[locus_identified[0]].append(all_stuff)
                        for p in snp[6]:
                            within_cluster_locus2isol_with_snp[locus_identified[0]].append(p)
                        within_cluster_locus2original_snp_positions[locus_identified[0]].append(snp[0])
                        
                        if len(locus_identified) == 3:
                            pos2effect = (str(locus_identified[2]), snp[5])
                        else:
                            pos2effect = (snp[0], snp[5])
                        within_cluster_locus2position2effect[locus_identified[0]].append(pos2effect)

        #including recombination
        else:
            if len(v[0][-1]) == 1 and args.terminal:
                print 'Discounting snps as on terminal branch', k
            else:
                # determine which locus or intergenic region SNPs fall in #
                for snp in v:
                    locus_identified = which_locus(snp[0])
                    if len(locus_identified) == 0:
                        locus_identified.append('end_chrm_intergenic')
                        locus_identified.append('NA')
                    within_cluster_locus2snp[locus_identified[0]].append(snp[5])
                    within_cluster_locus2product[locus_identified[0]]=(locus_identified[1])
                    within_cluster_locus2branch[locus_identified[0]].append(snp[2])
                    needed_pats = '..'.join(snp[6])
                    all_stuff  = snp[2]+'..'+snp[0]+'..'+needed_pats
                    within_cluster_locus_node2pos2isol[locus_identified[0]].append(all_stuff)
                    for p in snp[6]:
                        within_cluster_locus2isol_with_snp[locus_identified[0]].append(p)

                    within_cluster_locus2original_snp_positions[locus_identified[0]].append(snp[0])
                    if len(locus_identified) == 3:
                        pos2effect = (str(locus_identified[2]), snp[5])
                    else:
                        pos2effect = (snp[0], snp[5])
                    within_cluster_locus2position2effect[locus_identified[0]].append(pos2effect)


    # determine the count of snps per locus within each cluster #
    for l, s in within_cluster_locus2snp.iteritems():
        output_line = []
        total_non_synonymous = []
        total_synonymous = []
        total_intergenic = []
        product_wanted = within_cluster_locus2product.get(l)
        branches = within_cluster_locus2branch.get(l)
        isolates = within_cluster_locus2isol_with_snp.get(l)
        original_positions = within_cluster_locus2original_snp_positions.get(l)

        #####################################
        
        

        summarise_snps = Counter(s)
        syns = summarise_snps.get('Synonymous')
        non_syns = summarise_snps.get('Nonsynonymous')
        inters = summarise_snps.get('Intergenic')
        SNOPs = summarise_snps.get('SNOP (non-stop codon to stop condon)')
        STOPs = summarise_snps.get('STOP (stop codon to stop codon)')
        STIPs = summarise_snps.get('STIP (stop codon to non-stop codon)')

        output_line.append(cluster)
        output_line.append(len(node))
        output_line.append(len(num_in_cluster))
        output_line.append(l)
        output_line.append(product_wanted)
        
        if syns == None:
            output_line.append(0)
        else:
            output_line.append(syns)
            total_synonymous.append(syns)
        if non_syns == None:
            output_line.append(0)
        else:
            output_line.append(non_syns)
            total_non_synonymous.append(non_syns)
        if inters == None:
            output_line.append(0)
        else:
            output_line.append(inters)
        if SNOPs == None:
            output_line.append(0)
        else:
            output_line.append(SNOPs)
            total_non_synonymous.append(SNOPs)
        if STOPs == None:
            output_line.append(0)
        else:
            output_line.append(STOPs)
            total_synonymous.append(STOPs)
        if STIPs == None:
            output_line.append(0)
        else:
            output_line.append(STIPs)
            total_non_synonymous.append(STIPs)

        total_syns = sum(total_synonymous)
        total_non_syns = sum(total_non_synonymous)
        output_line.append(total_syns)
        output_line.append(total_non_syns)

        non_syn_positions = []
        syn_positions = []
        int_positions = []
        
        all_positions = within_cluster_locus2position2effect.get(l)
        for x in all_positions:
            if x[1] == 'Nonsynonymous' or x[1].startswith('SNOP') or x[1].startswith('STIP'):
                non_syn_positions.append(x[0])
            else:
                if x[1] == 'Synonymous' or x[1].startswith('STOP'):
                    syn_positions.append(x[0])
                else:
                    int_positions.append(x[0])

        if len(non_syn_positions) == 0:
            output_line.append('NA')
        else:
            output_line.append('..'.join(non_syn_positions))
        if len(syn_positions) == 0:
            output_line.append('NA')
        else:
            output_line.append('..'.join(syn_positions))
        if len(int_positions) == 0:
            output_line.append('NA')
        else:
            output_line.append('..'.join(int_positions))

        output_line.append('..'.join(branches))
        output_line.append('..'.join(isolates))
        output_line.append('..'.join(original_positions))

        tab_info = within_cluster_locus_node2pos2isol.get(l)
        annoying = '__'.join(tab_info)
        output_line.append(annoying)
            

        count_out = 0
        for out in output_line:
            count_out += 1
            if count_out == len(output_line):
                args.output.write(str(out)+'\n')
            else:
                args.output.write(str(out)+'\t')



#########################
        if output_line[0] == 'unclustered':
            pass
        else:
            # synonymous snps
            if int(output_line[11]) != 0:
                locus2syn[output_line[3]].append(int(output_line[11]))
                locus2syn_positions[output_line[3]].append(output_line[14])
            # non_synonymous snps
            if int(output_line[12]) != 0:
                locus2non_syn[output_line[3]].append(int(output_line[12]))
                locus2non_syn_positions[output_line[3]].append(output_line[13])
            # intergenic snps
            if int(output_line[7]) != 0:
                locus2int[output_line[3]].append(int(output_line[7]))
                locus2int_positions[output_line[3]].append(output_line[-2])
                
# sum all clustered isolates #
for k, v in locus2syn.iteritems():
    positions = locus2syn_positions.get(k)
    pos = '..'.join(positions)
    total = sum(v)
    product = locus2product.get(k)
    length = locus2length.get(k)
    obs_snps_bp = float(total)/float(length)
    args.concat_output.write('syn'+'\t'+str(k)+'\t'+str(product)+'\t'+str(total)+'\t'+str(obs_snps_bp)+'\t'+str(pos)+'\n')

for k, v in locus2non_syn.iteritems():
    positions = locus2non_syn_positions.get(k)
    pos = '..'.join(positions)
    total = sum(v)
    product = locus2product.get(k)
    length = locus2length.get(k)
    obs_snps_bp = float(total)/float(length)
    args.concat_output.write('non_syn'+'\t'+str(k)+'\t'+str(product)+'\t'+str(total)+'\t'+str(obs_snps_bp)+'\t'+str(pos)+'\n')

for k, v in locus2int.iteritems():
    positions = locus2int_positions.get(k)
    pos = '..'.join(positions)
    total = sum(v)
    product = locus2product.get(k)
    length = locus2length.get(k)
    obs_snps_bp = float(total)/float(length)
    args.concat_output.write('int'+'\t'+str(k)+'\t'+str(product)+'\t'+str(total)+'\t'+str(obs_snps_bp)+'\t'+str(pos)+'\n')
    













                
