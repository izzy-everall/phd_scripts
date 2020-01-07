### analyse pangenome script march 2016 ####

### import stuff ###
import argparse
from collections import defaultdict 
from collections import Counter

### arguments ###
parser = argparse.ArgumentParser()

parser.add_argument("-i","--input",type=argparse.FileType("r"),help="gene_presence_abscence.csv")
parser.add_argument("-l","--lanes",type=argparse.FileType("r"),help="file with tab delimited lane 2 cluster")
parser.add_argument("-n", "--number_not_in_group", type = int, help="percentage of non-DCC isolates that it is ok for a gene of interest in DCCs to be in")
parser.add_argument("-b", "--background_group_name", type = str, action = 'store', help="name of group that contains the background isolates (non_DCC, unclustered etc) if providing more than one don't use a comma seperated list!", nargs='+')
parser.add_argument("-d", "--number_of_isols_in_dataset", type=int, help="number of isolates used in pan genome dataset (to define core)")
parser.add_argument("-p", "--narrow_candidate_list", type=int, help="if you want to say look only at genes in greater than half of group")
parser.add_argument("-g","--group_of_interest", type = str, action = 'store', help="if you want information on just one specific group")
parser.add_argument("-c", "--counts", action='store_true', help="if you want dataset stats")
parser.add_argument("-s", "--subspecies", type=argparse.FileType("r"), help="file with tab delimited lane 2 subspecies")
parser.add_argument("-lc", "--lane2cluster", type=argparse.FileType("r"), help="file with tab delimited lane 2 cluster")
parser.add_argument("-o", "--output", type=argparse.FileType("w"), help="output file for specific group information MUST USE -g flag for this")


args=parser.parse_args()

print "WARNING: SCRIPT CURRENTLY REMOVES DUPLICATES and GENES WHICH ARE QC'D as HYPOTHETICAL OR MARKED INVESTIGATE"
#### global dictionaries and lists####
lane2cluster_non_cluster = {}
gene2information = {}
accessory_gene2lanes = {}
gene2only_present_in_clusters = {}
clusters= []
subspecies_total = []
lane2subspecies = {}
subspecies_total_in_pgdataset = []
clusters_total_in_pg_dataset = []
gene2subspecies_presence_absence = {}
abss = []
mass = []
boll = []
abss_mass = []
abss_boll = []
mass_boll = []
pangenome_division_stats = []
clustered = []
unclustered = []
all_clusters = []
lane2_all_clusters = {}
all_clusters_pgdataset = []
all_clusters_stats = []
DCC1 = []
DCC2 = []
DCC3 = []
non_DCC = []
non_DCC_DCC1 = []
non_DCC_DCC2 = []
non_DCC_DCC3 = []
DCC1_DCC2 = []
DCC2_DCC3 = []
DCC1_DCC3 = []
all_DCCs = []
lane2treegubbins={}
count_clusters=[]

#### function to partition Clusters data ####
for line in args.lanes:
	line = line.strip().split("\t")
	cluster_membership = (line[0]+"_"+line[1])
	#lane2cluster_non_cluster[line[0]] = (cluster_membership)
	lane2cluster_non_cluster[line[0]] = (line[1])
	clusters.append(line[1])
	count_clusters.append(line[1])
	if len(line) > 2:
		lane2treegubbins[line[0]]=(line[2])
	else:
		pass

#print lane2treegubbins

# enables stats on DCC gene presence absence
if len(set(clusters)) == 4:
	for k,v in lane2cluster_non_cluster.iteritems():
		if v == "DCC1":
			DCC1.append(k)
			DCC1_DCC3.append(k)
			DCC1_DCC2.append(k)
			all_DCCs.append(k)
			non_DCC_DCC1.append(k)
		elif v == "DCC2":
			DCC2.append(k)
			DCC2_DCC3.append(k)
			DCC1_DCC2.append(k)
			all_DCCs.append(k)
			non_DCC_DCC2.append(k)
		elif v == "DCC3":
			DCC3.append(k)
			DCC1_DCC3.append(k)
			DCC2_DCC3.append(k)
			all_DCCs.append(k)
			non_DCC_DCC3.append(k)
		else:
			non_DCC_DCC1.append(k)
			non_DCC_DCC3.append(k)
			non_DCC_DCC2.append(k)
			non_DCC.append(k)
else:
	pass


#print lane2cluster_non_cluster.get("10665_4#54")
#print Counter(clusters)

if args.subspecies:
	for line in args.subspecies:
		line = line.strip().split("\t")
		lane2subspecies[line[0]] = (line[1])
		subspecies_total.append(line[1])
		if line[1] == "ABSS":
			abss.append(line[0])
			abss_mass.append(line[0])
			abss_boll.append(line[0])
		elif line[1] == "BOLL":
			boll.append(line[0])
			abss_boll.append(line[0])
			mass_boll.append(line[0])
		elif line[1] == "MASS":
			mass.append(line[0])
			abss_mass.append(line[0])
			mass_boll.append(line[0])
	print Counter(subspecies_total)


if args.lane2cluster:
	for line in args.lane2cluster:
		line = line.strip().split("\t")
		lane2_all_clusters[line[0]] = (line[1])
		all_clusters.append(line[1])
		if line[1] == "clustered":
			clustered.append(line[0])
		elif line[1] == "unclustered":
			unclustered.append(line[0])



### number of different groups under analysis###
clusters = list(set(clusters))
### dictionary with cluster as key and empty list as value
cluster2gene2lanes = {key:[] for key in clusters}
#### if want to create a candidate list of genes over certain threshold filling group to cluster dictionary to get isolates ###
if args.narrow_candidate_list:
	cluster2list_of_lanes = {key:[] for key in clusters}
	for k, l in cluster2list_of_lanes.iteritems():
		for lane, val in lane2cluster_non_cluster.iteritems():
			if val == k:
				l.append(lane)
			else:
				pass


count_multi_copy_number = 0
#### looping through gene presence abscence file resulting in dictionary of non-core genes2lanes ####
line_num = 0
for line in args.input:
	line_num += 1
	if line_num == 1:
		line_check = line.strip().split("\t")
		for i in line_check:
			if "#" in i:
				subsp = lane2subspecies.get(i)
				subspecies_total_in_pgdataset.append(subsp)
				clus = lane2cluster_non_cluster.get(i)
				clusters_total_in_pg_dataset.append(clus)
				#all_clus = lane2_all_clusters.get(i)
				count_clus = lane2cluster_non_cluster.get(i)
				all_clusters_pgdataset.append(count_clus)

			else:
				pass
	else:
		#pan-genome into subspecies groups.
		if args.subspecies:
			lanes = []
			get_lanes = line.strip().split("\t")
			for i in get_lanes:
				if "#" in i:
					id_needed = i.split("_")
					#print id_needed
					lane_id = (id_needed[0]+"_"+id_needed[1])
					lanes.append(lane_id)
				else:
					pass
			if set(lanes).issubset(set(abss)):
				gene2subspecies_presence_absence[get_lanes[0]]=("ABSS")
				if set(lanes).issubset(set(DCC1)):
					pangenome_division_stats.append("DCC1")
				elif set(lanes).issubset(set(DCC2)):
					pangenome_division_stats.append("DCC2")
				elif set(lanes).issubset(set(DCC1_DCC2)):
					pangenome_division_stats.append("DCC1_DCC2")
				elif set(lanes).issubset(set(non_DCC_DCC1)):
					pangenome_division_stats.append("non_DCC_DCC1")
				elif set(lanes).issubset(set(non_DCC_DCC2)):
					pangenome_division_stats.append("non_DCC_DCC2")
				else:
					pangenome_division_stats.append("ABSS")
					pangenome_division_stats.append("non-DCC_unique")
			elif set(lanes).issubset(set(boll)):
				gene2subspecies_presence_absence[get_lanes[0]]=("BOLL")
				pangenome_division_stats.append("BOLL")
				pangenome_division_stats.append("non-DCC_unique")
			elif set(lanes).issubset(set(mass)):
				gene2subspecies_presence_absence[get_lanes[0]]=("MASS")
				if set(lanes).issubset(set(DCC3)):
					pangenome_division_stats.append("DCC3")
				elif set(lanes).issubset(set(non_DCC_DCC3)):
					pangenome_division_stats.append("non_DCC_DCC3")
				else:
					pangenome_division_stats.append("MASS")
					pangenome_division_stats.append("non-DCC_unique")
			elif set(lanes).issubset(set(abss_mass)):
				gene2subspecies_presence_absence[get_lanes[0]]=("ABSS_MASS")
				if set(lanes).issubset(set(DCC1_DCC3)):
					pangenome_division_stats.append("DCC1_DCC3")
				elif set(lanes).issubset(set(DCC2_DCC3)):
					pangenome_division_stats.append("DCC2_DCC3")
				elif  set(lanes).issubset(set(all_DCCs)):
					pangenome_division_stats.append("all_DCCs")
				else:
					pangenome_division_stats.append("ABSS_MASS")
					pangenome_division_stats.append("non-DCC_unique")
			elif set(lanes).issubset(set(abss_boll)):
				gene2subspecies_presence_absence[get_lanes[0]]=("ABSS_BOLL")
				pangenome_division_stats.append("ABSS_BOLL")
				pangenome_division_stats.append("non-DCC_unique")
			elif set(lanes).issubset(set(mass_boll)):
				gene2subspecies_presence_absence[get_lanes[0]]=("MASS_BOLL")
				pangenome_division_stats.append("MASS_BOLL")
			else:
				gene2subspecies_presence_absence[get_lanes[0]]=("core")
				pangenome_division_stats.append("core")

		if args.lane2cluster:
			lanes = []
			get_lanes = line.strip().split("\t")
			for i in get_lanes:
				if "#" in i:
					id_needed = i.split("_")
					#print id_needed
					lane_id = (id_needed[0]+"_"+id_needed[1])
					lanes.append(lane_id)
				else:
					pass
			if set(lanes).issubset(set(clustered)):
				#gene2subspecies_presence_absence[get_lanes[0]]=("ABSS")
				all_clusters_stats.append("clustered")
			elif set(lanes).issubset(set(unclustered)):
				#gene2subspecies_presence_absence[get_lanes[0]]=("BOLL")
				all_clusters_stats.append("unclustered")		
			else:
				all_clusters_stats.append("core")				


		remove_duplicates = line.split("\t")
		# at the moment cannot handle duplicates
		if int(remove_duplicates[4]) > args.number_of_isols_in_dataset or float(remove_duplicates[5]) > 1.00 or remove_duplicates[10] == "Investigate" or remove_duplicates[10].startswith("Hypothetical"):
			count_multi_copy_number += 1
			#print remove_duplicates
			pass
		else:
        	# how many isolates have this gene
			number_of_isols = line.count("#")
        	# only interested in non core genes 
			core_defined_99 = round((args.number_of_isols_in_dataset/float(100))*float(99))
			if number_of_isols < core_defined_99:
				#print number_of_isols, args.number_of_isols_in_dataset, line_num
				line = line.split("\t")
				gene_of_interest = str(line[0]+"\t"+line[1]+"\t"+line[2])
				gene2information[line[0]]=(gene_of_interest)
				gene_name=line[0]
				lanes = []
				for i in line:
					if "#" in i:
                    	#safely making lane names the same
						i = i.split("_")
						lane = (i[0]+"_"+i[1])
						lanes.append(lane)
					else:
						pass 
					accessory_gene2lanes[gene_name]=(lanes)
			else:
				pass                  


#print accessory_gene2lanes.get("group_152")
#print accessory_gene2lanes
#print count_multi_copy_number
if args.subspecies:
	print Counter(subspecies_total_in_pgdataset)
	print Counter(pangenome_division_stats)

if args.lane2cluster:
	#print Counter(all_clusters_statsclusters_total_in_pg_dataset)
	print Counter(all_clusters_stats)

#print Counter(clusters_total_in_pg_dataset)

#### reducing stringency list of lanes is currently complete single isol dataset optional with -n flag ###
if args.number_not_in_group:
	#print args.number_not_in_group
    #proportion_of_isols_with_gene1 = round((args.number_not_in_group/float(100))*len(lane2cluster_non_cluster))
    #unhash below command!!!! 
    proportion_of_isols_with_gene = round((args.number_of_isols_in_dataset/float(100))*args.number_not_in_group)

    ### should it be the threshold based on number in the group that is not of interest as opposed to the whole dataset.
	#print "TESTING DOING PROPORTIONS BASED ON NO. IN BACKGROUND CATERGORY AS OPPOSED TO WHOLE DATASET"
	#number_in_backgound = 0
	#for k,v in lane2cluster_non_cluster.iteritems():
#		if v == "unclustered":
#			number_in_backgound += 1
#	proportion_of_isols_with_gene = round((number_in_backgound/float(100))*args.number_not_in_group)

    #print proportion_of_isols_with_gene1, proportion_of_isols_with_gene2

#### counting dataset stats ####
background_group_unique_genes = 0
genes_of_interest_no_group_membership_threshold = 0
genes_of_interest_group_membership_threshold = 0


#### identify the accessory genomes of the clusters of interest###
for gene,lanes in accessory_gene2lanes.iteritems():
	membership = []
	accessory_genes_tree_gubbins = []
	for i in lanes:
		group_membership = lane2cluster_non_cluster.get(i)
		# error catcher -- mising data
		if group_membership == None:
			print "WARNING ISOLATE IN PAN_GENOME HAS NOT BEEN GIVEN CATERGORY!", i, group_membership
		#print i, group_membership
		membership.append(group_membership)
		if len(lane2treegubbins) > 1:
			tree_gubbins_group = lane2treegubbins.get(i)
			accessory_genes_tree_gubbins.append(tree_gubbins_group)
		else:
			pass


	# if these are genes just present in group that isn't of interest:
	#print len(set(membership)), membership[0], args.background_group_name
	#if "DCC2" in membership:
		#print membership, len(membership)
	if len(set(membership)) == 1 and membership[0] in args.background_group_name:
		if args.counts == True:
			background_group_unique_genes += 1
		else:
			pass
	else:
		### work out proportions of isolates in each catergory
		proportions_of_memberships = Counter(membership)
		#print proportions_of_memberships
		### specifically identify proportion in background population
		if len(args.background_group_name) > 1:
			total_in_background = []
			for i in args.background_group_name:
				proportion_in_background = proportions_of_memberships.get(i)
				if proportion_in_background == None:
					pass
				else:
					total_in_background.append(proportion_in_background)
			proportion_in_background = sum(total_in_background)
		else:
			proportion_in_background = proportions_of_memberships.get(args.background_group_name[0])
		### is that proportion insignificant enough to still be of interest as a marker in groups of interest
		if proportion_in_background <= proportion_of_isols_with_gene or proportion_in_background == None:
			if args.counts == True:
				genes_of_interest_no_group_membership_threshold += 1
			#print proportions_of_memberships
			information = gene2information.get(gene)
			## is there a sensible number of the representative group to be meaningful
			## is there greater than the user defined threshold of isolates 
			for k,v in proportions_of_memberships.iteritems():
				#how_many_in_cluster = len(cluster2list_of_lanes.get(k))
				get_how_many_in_cluster = Counter(all_clusters_pgdataset)
				#print get_how_many_in_cluster
				how_many_in_cluster = get_how_many_in_cluster.get(k)
				threshold = round((args.narrow_candidate_list/float(100))*how_many_in_cluster)
				#print information, proportions_of_memberships, k, threshold, len(lanes)
				##### just information on a user defined group ######
				if v >= threshold and k == args.group_of_interest:
					#print information, proportions_of_memberships, threshold, len(lanes)
					if args.counts == True:
						genes_of_interest_group_membership_threshold += 1
					#print information, Counter(tree_gubbins)
					args.output.write(k+"\t"+information+"\t")
					info_for_group_of_interest = proportions_of_memberships.get(args.group_of_interest)
					args.output.write(args.group_of_interest+":"+str(info_for_group_of_interest))
					### wrtiting this to an output file
					count_diff_clusters = Counter(accessory_genes_tree_gubbins)
					even_tabs = 0
					total_not_in_group = []
					for k, v in proportions_of_memberships.iteritems():
						if k == args.group_of_interest and len(proportions_of_memberships) == 1:
							args.output.write("\t\t"+str(threshold)+"\t"+str(len(lanes))+"\n")
						elif k != args.group_of_interest:
							if len(proportions_of_memberships) == 2:
								if len(count_diff_clusters) > 0:
									args.output.write(","+k+":"+str(v)+"\t"+str(v)+"\t"+str(threshold)+"\t"+str(len(lanes))+",")
									for tg,n in count_diff_clusters.iteritems():
										args.output.write(tg+":"+str(n)+",")
									args.output.write("\n")
								else:
									args.output.write(","+k+":"+str(v)+"\t"+str(v)+"\t"+str(threshold)+"\t"+str(len(lanes))+"\n")
							elif len(proportions_of_memberships) == 3:
								even_tabs += 1
								if even_tabs == 1:
									total_not_in_group.append(v)
									args.output.write(","+k+":"+str(v)+",")
								elif even_tabs == 2:
									total_not_in_group.append(v)
									args.output.write(k+":"+str(v)+"\t"+str(sum(total_not_in_group))+"\t"+str(threshold)+"\t"+str(len(lanes))+"\n")
							elif len(proportions_of_memberships) == 4:
								even_tabs += 1
								if even_tabs == 1:
									total_not_in_group.append(v)
									args.output.write(","+k+":"+str(v)+",")
								elif even_tabs == 2:
									total_not_in_group.append(v)
									args.output.write(k+":"+str(v)+",")
								elif even_tabs == 3:
									total_not_in_group.append(v)
									args.output.write(k+":"+str(v)+"\t"+str(sum(total_not_in_group))+"\t"+str(threshold)+"\t"+str(len(lanes))+"\n")

				else:
					pass
			else:
				pass



if args.counts == True:
	print "total number of accessory genes:", len(accessory_gene2lanes)
	print "Genes removed due to only being present in the background:", background_group_unique_genes
#	print "Total number of genes of interest that aren't dominant in the background:", genes_of_interest_no_group_membership_threshold
	print "Total number of genes with a significant majority in a group of interest:", genes_of_interest_group_membership_threshold










