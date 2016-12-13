#from basic_functions_old import *
from basic_functions_new_RNAplfold_format import *
import operator
import numpy
import scipy.stats
import random


def one_poss_to_more_degerated_mers(pos, neg, previous_max_AUC, previous_min_P,fix_poss, seq_file_name, model, testnum, RNAplfold_direct):
    
    fix_mer = poss_to_degernate_motif(fix_poss)
    

    seq_f = open(seq_file_name,'r')
    seq_check = []
    for line in seq_f:
        if line.startswith('>'):
            seq_check.append(line[1:])
    seq_f.close()
    
   
    total_pos = 0
    for line in pos:
        if seq_check.count(line) > 0:
            total_pos = total_pos + 1
   
   
    total_neg = 0
    for line in neg:
        if seq_check.count(line) > 0:
            total_neg = total_neg + 1

    
    filtered_access_total = {}
    
    
    # if positive and negative gene_ID not change, the filtered_access wont change
    
    #keep shifting on the same level of substitution
    

    nt_poss = more_degenerate_poss(fix_poss)
    nt_poss.append(fix_poss)
    
    left_shift = shift(fix_poss, 'Left')
    left_poss = []
    for item in left_shift:
        left_poss.append(item)
    right_shift = shift(fix_poss, 'Right')
    right_poss = []
    for item in right_shift:
        right_poss.append(item)
	
    left_less_shift = shift_less(fix_poss, 'Left')
    left_less_poss = []
    for item in left_less_shift:
        left_less_poss.append(item)
    right_less_shift = shift_less(fix_poss, 'Right')
    right_less_poss = []
    for item in right_less_shift:
        right_less_poss.append(item)
    # print (right_less_poss)
    # print (left_less_poss)
    if len(fix_poss) < 10:
        all_poss = [nt_poss, left_poss, right_poss, left_less_poss, right_less_poss]
    else:
        all_poss = [nt_poss, left_less_poss, right_less_poss]
    all_left_mer = {}
    temp_mer = fix_mer
    temp_result = [previous_max_AUC, previous_min_P,0,0]
    testnum = testnum + len(nt_poss) + len(left_poss) + len(right_poss) + len(left_less_poss) + len(right_less_poss)
    for posses in all_poss:
        
        for poss in posses:
           
	    motifs = required_motif(poss)
	    loc_dic = loc_whole_region(motifs, seq_file_name)
	    if filtered_access_total.has_key(len(poss)):
		if model =='access_seq' or model == 'access_only':
		    access = sum_region_accessibility(filtered_access_total[len(poss)], loc_dic)
		elif model == 'seq' or model == 'site_only':
		    access = sum_site_num(filtered_access_total[len(poss)], loc_dic)
		
		pos_access = access[0].values()
		neg_access = access[1].values()
		pos_zero = total_pos - len(pos_access)
		neg_zero = total_neg - len(neg_access)
    
		pos_access_seq = pos_access[:]
		neg_access_seq = neg_access[:]
		
                if model == 'access_seq' or model == 'seq':
                    for i in range(0, pos_zero):
                        pos_access_seq.append(0)
                    for j in range(0, neg_zero):
                        neg_access_seq.append(0)
                    result = AUC(pos_access_seq, neg_access_seq)
                elif model == 'access_only' or model == 'site_only':
                    result = AUC(pos_access_seq, neg_access_seq)
                
	    else:
		site_len = len(poss)
		#RNAplfold_file = RNAplfold_direct + str(site_len) + '_whole.txt'
		filtered_access = filtered_accessibility(RNAplfold_direct, pos, neg, site_len)
		filtered_access_total = {}
		filtered_access_total[site_len] = filtered_access
		if model == 'access_seq' or model == 'access_only':
		    access = sum_region_accessibility(filtered_access, loc_dic)
		elif model == 'seq'or model == 'site_only':
		    access = sum_site_num(filtered_access_total[len(poss)], loc_dic)
	    
		pos_access = access[0].values()
		neg_access = access[1].values()
		pos_zero = total_pos - len(pos_access)
		neg_zero = total_neg - len(neg_access)
    
		pos_access_seq = pos_access[:]
		neg_access_seq = neg_access[:]

                if model == 'access_seq' or model == 'seq':
                    for i in range(0, pos_zero):
                        pos_access_seq.append(0)
                    for j in range(0, neg_zero):
                        neg_access_seq.append(0)
                    result = AUC(pos_access_seq, neg_access_seq)
                elif model == 'access_only' or model == 'site_only':
                    result = AUC(pos_access_seq, neg_access_seq)

            if poss == fix_poss:
		print poss
		print result 
		print pos_zero
		print neg_zero
	    if result[1] * testnum < 0.05 and result[0] > temp_result[0]:
	    
	        k_mer = poss_to_degernate_motif(poss)
		temp_mer = k_mer
		temp_result = result
                   
    if temp_mer != fix_mer:
        all_left_mer[temp_mer] = temp_result
	  
    print testnum        
    return [all_left_mer, testnum]

def one_poss_to_fix_mer(pos, neg,  fix_mer_groups, seq_file_name, model, testnum, RNAplfold_direct):
   
    seq_f = open(seq_file_name,'r')
    seq_check = []
    for line in seq_f:
        if line.startswith('>'):
            seq_check.append(line[1:])
    seq_f.close()
    
    
    total_pos = 0
    for line in pos:
        if seq_check.count(line) > 0:
            total_pos = total_pos + 1
   
    total_neg = 0
    for line in neg:
        if seq_check.count(line) > 0:
            total_neg = total_neg + 1

    
    filtered_access_total = {}

    all_left_mer = {}
    
    result_AUCs = []
    result_Ps = []
    for temp_mer in fix_mer_groups:
        
	motifs = required_motif(consensus_to_motifs(temp_mer))
	loc_dic = loc_whole_region(motifs, seq_file_name)
	if filtered_access_total.has_key(len(motifs[0])):
	    if model =='access_seq' or model == 'access_only':
		access = sum_region_accessibility(filtered_access_total[len(motifs[0])], loc_dic)
	    elif model == 'seq' or model == 'site_only':
		
		access = sum_site_num(filtered_access_total[len(motifs[0])], loc_dic)
		
	    pos_access = access[0].values()
	    neg_access = access[1].values()
	    pos_zero = total_pos - len(pos_access)
	    neg_zero = total_neg - len(neg_access)
    
	    pos_access_seq = pos_access[:]
	    neg_access_seq = neg_access[:]
    
	    for i in range(0, pos_zero):
		pos_access_seq.append(0)
	    for j in range(0, neg_zero):
		neg_access_seq.append(0)
	    result = AUC(pos_access_seq, neg_access_seq)
	else:
	    site_len = len(motifs[0])
	    #RNAplfold_file = RNAplfold_direct + str(site_len) + '_whole.txt'
	    filtered_access = filtered_accessibility(RNAplfold_direct, pos, neg, site_len)
	    filtered_access_total = {}
	    filtered_access_total[site_len] = filtered_access
	    if model == 'access_seq' or model == 'access_only':
		access = sum_region_accessibility(filtered_access, loc_dic)
	    elif model == 'seq'or model == 'site_only':
		access = sum_site_num(filtered_access_total[len(motifs[0])], loc_dic)
	    
	    pos_access = access[0].values()
	    neg_access = access[1].values()
	    pos_zero = total_pos - len(pos_access)
	    neg_zero = total_neg - len(neg_access)
    
	    pos_access_seq = pos_access[:]
	    neg_access_seq = neg_access[:]

	    if model == 'access_seq' or model == 'seq':
		for i in range(0, pos_zero):
		    pos_access_seq.append(0)
		for j in range(0, neg_zero):
		    neg_access_seq.append(0)
		result = AUC(pos_access_seq, neg_access_seq)
               
	    elif model == 'access_only' or model == 'site_only':
		result = AUC(pos_access_seq, neg_access_seq)
 	print temp_mer
	print result             
	result_AUCs.append(result[0])
	result_Ps.append(result[1])
	    
	
    return [max(result_AUCs), min(result_Ps)*testnum]



def mers_group_to_mer_group(pos, neg, previous_max_AUC,previous_min_P, fix_mer_groups, seq_file_name, model, testnum, RNAplfold_direct):
   
    result = []
    AUCs = []
    Ps = []
    
    #file_in = open(output,'a')
    #file_in.write('new' + '\n')
    for fix_mer in fix_mer_groups:
        fix_poss = consensus_to_motifs(fix_mer)
	
	all_mers_d_all = one_poss_to_more_degerated_mers(pos, neg, previous_max_AUC,previous_min_P,fix_poss, seq_file_name, model, testnum, RNAplfold_direct)
	all_mers_d = all_mers_d_all[0]
        testnum = all_mers_d_all[1]
      

        
	all_mers = all_mers_d.keys()
	for key in all_mers:
	    if result.count(key) == 0:
		#file_in.write(key + '\t' + str(all_mers_d[key]) + '\n')
		AUCs.append(all_mers_d[key][0])
		Ps.append(all_mers_d[key][1])
		result.append(key)
        
    #file_in.close()
    if len(AUCs) > 0:
	max_AUC = max(AUCs)
	min_P = min(Ps)*testnum
    else:
	max_AUC = previous_max_AUC
	min_P = previous_min_P
    return [result, max_AUC,testnum,min_P]

def run2(pos, neg,fix_mer_groups, seq_file_name, model, testnum, RNAplfold_direct):
    '''accessibility_seq model'''
    
    
    [start_max_AUC, start_min_P] = one_poss_to_fix_mer(pos, neg,fix_mer_groups, seq_file_name, model, testnum, RNAplfold_direct)
	
    #file_in = open(output,'w')
    #file_in.write(str(start_max_AUC) + '\t' + str(start_min_P) + '\t' +  str(fix_mer_groups) + '\t'  + '\t' + str(model)+ '\n')
    #file_in.close()
    print 'test'
    print start_max_AUC
    print start_min_P
    print 'test'
    if len(fix_mer_groups) >= 1:
        print fix_mer_groups
        results = mers_group_to_mer_group(pos, neg,start_max_AUC, start_min_P, fix_mer_groups, seq_file_name, model,testnum, RNAplfold_direct)
	
	previous_max_AUC = results[1]
        testnum = results[2]
	previous_min_P = results[3]
				
	if len(results[0]) < 1:
	    previous_result = [fix_mer_groups, start_max_AUC, 4096, start_min_P]
	    results2 = results
	    
	else:
	    i = 0
            previous_result = [fix_mer_groups, start_max_AUC, 4096, start_min_P]
            results2 = results
	    while i <20 and len(results[0]) >0:
		
		results2 = mers_group_to_mer_group(pos, neg,previous_max_AUC, previous_min_P, results[0], seq_file_name, model,testnum, RNAplfold_direct)
                testnum = results2[2]
		if previous_max_AUC < results2[1]: 
               
                    previous_result = results2
                    previous_max_AUC = results2[1]
                    previous_min_P = results2[3]
                    results = results2
                    i = i + 1
		else:
                
		    break
        
    return [previous_result, results2]



def Motif_discovery(final_out_name, detailed_final_out_name,model, pos_file_name,neg_file_name, seq_file_name, RNAplfold_direct):
    
    #testnum = 4096
    #top_num_mer = 4
    #final_out = open(final_out_name,'w')
    
    #RNAplfold_file = RNAplfold_direct + '6_whole.txt'
    genes = []
    
    
    seq_f = open(seq_file_name,'r')
    seqs = []
    for line in seq_f:
	if line.startswith('>'):
	    genes.append(line[1:])
    seq_f.close()
    
    mer_file = open('mer_6.txt','r')
    #mer_file = open('/home/xiaoli/gene_mer_fly/random_mer_pool/mer_6.txt','r')
    mer_6 = []
    for line in mer_file:
	mer_6.append([line.strip()])

    testnum = len(mer_6)
    
    filtered_accesses = filtered_accessibility(RNAplfold_direct,pos_file_name, neg_file_name, 6)
    result_access = {}
    result_site = {}
    row_label = []
    
    filtered_access = filtered_accesses[0]
    for a in filtered_accesses[1].keys():
	filtered_access[a] = filtered_accesses[1][a] 
    
    for gene in genes:
	if filtered_access.has_key(gene.strip()):
	    result_access[gene.strip()] = []
	    result_site[gene.strip()] = []
	    row_label.append(gene.strip())
	    
    for kmer in mer_6:
	loc_dic = loc_whole_region(kmer, seq_file_name)
	
	access_site = sum_region_accessibility_site(filtered_access, loc_dic)
	access = access_site[0]
	site = access_site[1]


	for gene in row_label:
	    if access.has_key(gene):
		result_access[gene].append(access[gene])
		result_site[gene].append(site[gene])
	    else:
		result_access[gene].append('0')
		result_site[gene].append('0')


    #file_access = open('/home/xiaoli/yeast_Hogen/motif_discovery/mer_6/access_6.txt','w')
    #file_site = open('/home/xiaoli/yeast_Hogen/motif_discovery/mer_6/site_6.txt','w')
    pos_f = open(pos_file_name,'r')
    pos = pos_f.readlines()
    pos_f.close()
    neg_f = open(neg_file_name,'r')
    neg = neg_f.readlines()
    neg_f.close()
    
    pos_d = {}
    neg_d = {}
    
    if model == 'access_seq':
	result_model = result_access
    elif model == 'seq':
	result_model = result_site
    else:
	print 'Type the right model!!'
	
    for gene_line in pos:
	gene = gene_line.strip()
	if result_model.has_key(gene):
	    temp = result_model[gene]
	    for i in range(len(mer_6)):
		if pos_d.has_key(mer_6[i][0]):
                    pos_d[mer_6[i][0]].append(float(temp[i]))
                else:
                    pos_d[mer_6[i][0]] = [float(temp[i])]
    for gene_line in neg:
	gene = gene_line.strip()
	if result_access.has_key(gene):
	    temp = result_model[gene]
	    for i in range(len(mer_6)):
		if neg_d.has_key(mer_6[i][0]):
                    neg_d[mer_6[i][0]].append(float(temp[i]))
                else:
                    neg_d[mer_6[i][0]] = [float(temp[i])]

    #print len(neg_d)
    mer_AUC_coll = []
    mer_result_coll = []
    mer_P_coll = []
    result_temp = []#
    for each_mer in mer_6:
	if len(pos_d[each_mer[0]]) == 0 or len(neg_d[each_mer[0]]) == 0:
	    print each_mer
	else:
	    result = AUC(pos_d[each_mer[0]], neg_d[each_mer[0]])
            result_temp.append(float(result[0])) #
	#print result[0]
        if result[0] > 0.50:
            mer_AUC_coll.append(-(result[0]))
	    mer_P_coll.append(result[1])
            mer_result_coll.append([each_mer[0], result])
    IDs=sorted(enumerate(mer_AUC_coll), key = operator.itemgetter(1))
    fix_mer_groups = []
    
    #print len(final_result)
    results = []
    AUCs = []
    detailed_result_file = open(detailed_final_out_name,'w')
    print len(mer_AUC_coll)
    print max(result_temp)#
    print mer_AUC_coll[0]
    for top_num_mer in range(5):
	if top_num_mer < len(mer_AUC_coll):
	    i = top_num_mer
	    print i
	    
	    fix_mer_groups = [mer_result_coll[IDs[i][0]][0]]
	    fix_mer_AUC = IDs[i][1]
	    fix_mer_P = mer_P_coll[IDs[i][0]]
	    if fix_mer_AUC == mer_AUC_coll[IDs[i][0]]:
		print 'yes'
		print fix_mer_AUC
		print fix_mer_P
		
	    print fix_mer_groups
	    
	    a = run2(pos, neg,fix_mer_groups, seq_file_name, model, testnum, RNAplfold_direct)
	    #print a 
	    if len(a[1][0]) == 0:
		all_result = a[0][:]
	    else:
		all_result = a[1][:]
	    
	    #For detailed information
	    detailed_result_file.write('Initial seed: ' + fix_mer_groups[0] + '\t' + \
				       str(-(fix_mer_AUC)) + '\t' + str(fix_mer_P) + '\n')
	    detailed_result_file.write('Motif is: ' + all_result[0][0] + '\n')
	    detailed_result_file.write('AUROC is: ' + str(all_result[1]) + '\n')
	    detailed_result_file.write('Ranksum P value is: ' + str(all_result[3]) + '\n\n')
	    #detailed_result_file.close()
	    #
	    
	    results.append(all_result)
	    AUCs.append(all_result[1])
        
    if len(AUCs)>0:        
        final_result = results[AUCs.index(max(AUCs))]
	    
        final_file = open(final_out_name,'w')
	#print len(final_result)
	#print final_result
	    
        final_file.write('Motif is: ' + final_result[0][0] + '\n')
        final_file.write('AUROC is: ' + str(final_result[1]) + '\n')
        final_file.write('Ranksum P value is: ' + str(final_result[3]) + '\n')
        final_file.close()
    detailed_result_file.close()
			       
			       
        #final_out.write(str(output) + '\t' + str(a) + '\n')
    #final_out.close()
    
    return all_result 

def Motif_discovery_for_cross_validation(model, pos,neg, seq_file_name, RNAplfold_direct):

    #testnum = 4096
    #top_num_mer = 4
    #final_out = open(final_out_name,'w')

    #RNAplfold_file = RNAplfold_direct + '6_whole.txt'
    genes = []


    seq_f = open(seq_file_name,'r')
    seqs = []
    for line in seq_f:
        if line.startswith('>'):
            genes.append(line[1:])
    seq_f.close()

    mer_file = open('mer_6.txt','r')
    #mer_file = open('/home/xiaoli/gene_mer_fly/random_mer_pool/mer_6.txt','r')
    mer_6 = []
    for line in mer_file:
        mer_6.append([line.strip()])

    testnum = len(mer_6)

    filtered_accesses = filtered_accessibility(RNAplfold_direct,pos, neg, 6)
    result_access = {}
    result_site = {}
    row_label = []

    filtered_access = filtered_accesses[0]
    for a in filtered_accesses[1].keys():
        filtered_access[a] = filtered_accesses[1][a]

    for gene in genes:
        if filtered_access.has_key(gene.strip()):
            result_access[gene.strip()] = []
            result_site[gene.strip()] = []
            row_label.append(gene.strip())

    for kmer in mer_6:
        loc_dic = loc_whole_region(kmer, seq_file_name)

        access_site = sum_region_accessibility_site(filtered_access, loc_dic)
        access = access_site[0]
        site = access_site[1]


        for gene in row_label:
            if access.has_key(gene):
                result_access[gene].append(access[gene])
                result_site[gene].append(site[gene])
            else:
                result_access[gene].append('0')
                result_site[gene].append('0')


    #file_access = open('/home/xiaoli/yeast_Hogen/motif_discovery/mer_6/access_6.txt','w')
    #file_site = open('/home/xiaoli/yeast_Hogen/motif_discovery/mer_6/site_6.txt','w')
    pos_d = {}
    neg_d = {}

    if model == 'access_seq':
        result_model = result_access
    elif model == 'seq':
        result_model = result_site
    else:
        print 'Type the right model!!'

    for gene_line in pos:
        gene = gene_line.strip()
        if result_model.has_key(gene):
            temp = result_model[gene]
            for i in range(len(mer_6)):
                if pos_d.has_key(mer_6[i][0]):
                    pos_d[mer_6[i][0]].append(float(temp[i]))
                else:
                    pos_d[mer_6[i][0]] = [float(temp[i])]
    for gene_line in neg:
        gene = gene_line.strip()
        if result_access.has_key(gene):
            temp = result_model[gene]
            for i in range(len(mer_6)):
                if neg_d.has_key(mer_6[i][0]):
                    neg_d[mer_6[i][0]].append(float(temp[i]))
                else:
                    neg_d[mer_6[i][0]] = [float(temp[i])]

    #print len(neg_d)
    mer_AUC_coll = []
    mer_result_coll = []
    mer_P_coll = []
    result_temp = []#
    for each_mer in mer_6:
        if len(pos_d[each_mer[0]]) == 0 or len(neg_d[each_mer[0]]) == 0:
            print each_mer
        else:
            result = AUC(pos_d[each_mer[0]], neg_d[each_mer[0]])
            result_temp.append(float(result[0])) #
        #print result[0]
        if result[0] > 0.50:
            mer_AUC_coll.append(-(result[0]))
            mer_P_coll.append(result[1])
            mer_result_coll.append([each_mer[0], result])
    IDs=sorted(enumerate(mer_AUC_coll), key = operator.itemgetter(1))
    fix_mer_groups = []

    #print len(final_result)
    results = []
    AUCs = []

    for top_num_mer in range(5):
        if top_num_mer < len(mer_AUC_coll):
            i = top_num_mer
            print i

            fix_mer_groups = [mer_result_coll[IDs[i][0]][0]]
            fix_mer_AUC = IDs[i][1]
            fix_mer_P = mer_P_coll[IDs[i][0]]
            if fix_mer_AUC == mer_AUC_coll[IDs[i][0]]:
                print 'yes'
                print fix_mer_AUC
                print fix_mer_P

            print fix_mer_groups

            a = run2(pos, neg,fix_mer_groups, seq_file_name, model, testnum, RNAplfold_direct)
            #print a 
            if len(a[1][0]) == 0:
                all_result = a[0][:]
            else:
                all_result = a[1][:]

            results.append(all_result)
            AUCs.append(all_result[1])

    if len(AUCs)>0:
        final_result = results[AUCs.index(max(AUCs))]
    return final_result[0][0] # the motif identified

def generate_cross_validation_fold(pos_file_name,neg_file_name, fold_num):
    results = []
    pos_f = open(pos_file_name,'r')
    total_pos_ori = pos_f.readlines()
    pos_f.close()
    neg_f = open(neg_file_name, 'r')
    total_neg_ori = neg_f.readlines()
    neg_f.close() 
    #randomize the positive and negative sets
    total_pos = random.sample(total_pos_ori, len(total_pos_ori))
    total_neg = random.sample(total_neg_ori, len(total_neg_ori))
    
    # split the positive and negative sets according to the fold_number
    aa = len(total_pos)/(fold_num)
    bb = len(total_neg)/(fold_num)
    start_ids_pos = range(0, aa*fold_num,(aa))
    start_ids_neg = range(0, bb*fold_num,(bb))

    pos_groups=[] # small groups of positives 
    neg_groups=[] # small groups of negatives
    for i in range(len(start_ids_pos)):
        start_id = start_ids_pos[i]
        if i != len(start_ids_pos)-1:
                end_id = start_ids_pos[i+1]
        else:
                end_id = start_ids_pos[i]+ aa
        pos_groups.append(total_pos[start_id: end_id]) # evenly distributed
    #adding the remaining to the list:
    r_id = 0
    for item in total_pos[end_id:len(total_pos)]:
        pos_groups[r_id].append(item)
        r_id = r_id + 1

    for i in range(len(start_ids_neg)):
        start_id = start_ids_neg[i]
        if i != len(start_ids_neg)-1:
                end_id = start_ids_neg[i+1]
        else:
                end_id = start_ids_neg[i] + bb
        neg_groups.append(total_neg[start_id: end_id])
    r_id = 0
    for item in total_neg[end_id:len(total_neg)]:
        neg_groups[r_id].append(item)
        r_id = r_id + 1
    
    # to run the train and test sets
    if len(pos_groups) !=len(neg_groups):
        print 'error!!!'
    else:
        for i in range(len(pos_groups)):
        	all_ids = range(len(pos_groups)) # len(pos_groups) == len(neg_groups)
                test_id = i
                all_ids.remove(test_id)
                trainning_ids = all_ids[:]
                test_pos = pos_groups[test_id]
                test_neg = neg_groups[test_id]
                trainning_pos = []
                trainning_neg = []
                for item in trainning_ids:
                        trainning_pos = trainning_pos + pos_groups[item]
                        trainning_neg = trainning_neg + neg_groups[item]
		results.append([trainning_pos, trainning_neg, test_pos, test_neg])
    return results

def run_motif_discovery_for_each_cv_fold(final_out_name, model, trainning_pos, trainning_neg, seq_file_name, RNAplfold_direct, test_pos, test_neg):
	file_out = open(final_out_name,'w')
	file_out.write('Train motif\tTest AUC\tTest wmwP\n')

	train_motif = Motif_discovery_for_cross_validation(model,trainning_pos,trainning_neg, seq_file_name, RNAplfold_direct)
        # to do motif scan on the test set              
        if model == 'access_seq':
		result_test = ATS_with_zero(train_motif, test_pos,test_neg, seq_file_name, RNAplfold_direct)
                file_out.write(train_motif + '\t' + str(result_test[0][0]) + '\t' + str(result_test[0][1]) + '\n')
        else:
                print 'Type the right model!!'
	file_out.close()


def ATS(con, pos_file_name,neg_file_name, seq_file_name, RNAplfold_direct):
    motifs = required_motif(consensus_to_motifs(con))
    site_len = len(motifs[0])
    #RNAplfold_file = RNAplfold_direct  + str(site_len) + '_whole.txt'
    loc_dic = loc_whole_region(motifs, seq_file_name)
    filtered_access = filtered_accessibility(RNAplfold_direct, pos_file_name, neg_file_name, site_len)
    access = sum_region_accessibility(filtered_access, loc_dic)
    result = AUC(access[0].values(), access[1].values())
    
    print('AUCROC: ' + str(result[0]))
    print ('Ranksum P value: ' + str(result[1]))
    print ('positive number: ' + str(len(access[0])))
    print ('negative number: ' + str(len(access[1])))
    #return [access[0].values(), access[1].values()]
    return result

def ATS_with_zero(con, pos_file_name,neg_file_name, seq_file_name, RNAplfold_direct):

    
    motifs = required_motif(consensus_to_motifs(con))
    
    site_len = len(motifs[0])
    #RNAplfold_file = RNAplfold_direct + str(site_len) + '_whole.txt'
    loc_dic = loc_whole_region(motifs, seq_file_name)
    filtered_access = filtered_accessibility(RNAplfold_direct, pos_file_name, neg_file_name, site_len)
    access = sum_region_accessibility(filtered_access, loc_dic)
    
    
    seq_ID = []
    seq_f = open(seq_file_name,'r')
    for line in seq_f:
	if line.startswith('>'):
	    seq_ID.append(line[1:].strip())
    if type(pos_file_name)==str:
    	pos_ID_o = []
    	pos_f = open(pos_file_name,'r')
    	for line in pos_f:
		if seq_ID.count(line.strip()) > 0 and pos_ID_o.count(line.strip())==0:
	    		pos_ID_o.append(line.strip())
    else:
	pos_ID_o = pos_file_name
    if type(neg_file_name)==str:
    	neg_ID_o = []
    	neg_f = open(neg_file_name,'r')
    	for line in neg_f:
		if seq_ID.count(line.strip()) > 0 and neg_ID_o.count(line.strip())==0:
	    		neg_ID_o.append(line.strip())
    else:
	neg_ID_o = neg_file_name
    pos_ID = []
    neg_ID = []
    for gene in pos_ID_o:
	if neg_ID_o.count(gene)==0:
		pos_ID.append(gene)
    for gene in neg_ID_o:
	if pos_ID_o.count(gene)==0:
		neg_ID.append(gene)
    
    pos_keys = access[0].keys()
    neg_keys = access[1].keys()
     
    total_access_pos = {}
    total_access_neg = {}
    record_pos = 0
    record_neg = 0
    for gene in pos_ID:
	gene = gene.strip()
	if pos_keys.count(gene)>0:
		total_access_pos[gene] = access[0][gene]
	else:
		total_access_pos[gene] = 0
		record_pos = record_pos + 1
    for gene in neg_ID:
	gene = gene.strip()
	if neg_keys.count(gene)>0:
		total_access_neg[gene] = access[1][gene]
	else:
		total_access_neg[gene] = 0
		record_neg = record_neg +1 

    result = AUC(total_access_pos.values(), total_access_neg.values())
    print "pos not have the sites" + str(len(pos_ID) - len(access[0]))
    print "neg not have the sites" + str(len(neg_ID) - len(access[1]))
    
    for item in set(pos_ID + neg_ID):
	if pos_ID.count(item) + neg_ID.count(item) > 1:
	    print item
    print('AUCROC: ' + str(result[0]))
    print ('Ranksum P value: ' + str(result[1]))
    #return access
    #return [access[0].values(), access[1].values()]
    return [result,[total_access_pos,total_access_neg]]

def ATS_with_zero_region(con, pos_file_name,neg_file_name, seq_file_name, RNAplfold_direct, index_d):
    print 'here'
    motifs = required_motif(consensus_to_motifs(con))
    site_len = len(motifs[0])
    #RNAplfold_file = RNAplfold_direct + str(site_len) + '_whole.txt'
    loc_dic_o = loc_whole_region(motifs, seq_file_name)
    
    loc_dic = {}
    genes = loc_dic_o.keys() 
    for gene in genes:
	locs = []
	if index_d.has_key(gene):
		locs_o = loc_dic_o[gene]
		start_index = index_d[gene][0]
		end_index = index_d[gene][1]
		for loc in locs_o:
			if loc>=start_index and loc+len(con)<=end_index:
				locs.append(loc)
	else:
		gene
	loc_dic[gene] = locs
    
    filtered_access = filtered_accessibility(RNAplfold_direct, pos_file_name, neg_file_name, site_len)
    access = sum_region_accessibility(filtered_access, loc_dic)


    seq_ID = []
    pos_ID_o = []
    neg_ID_o = []
    seq_f = open(seq_file_name,'r')
    for line in seq_f:
        if line.startswith('>'):
            seq_ID.append(line[1:].strip())
    pos_f = open(pos_file_name,'r')
    for line in pos_f:
        if seq_ID.count(line.strip()) > 0 and pos_ID_o.count(line.strip())==0:
            pos_ID_o.append(line.strip())
    neg_f = open(neg_file_name,'r')
    for line in neg_f:
        if seq_ID.count(line.strip()) > 0 and neg_ID_o.count(line.strip())==0:
            neg_ID_o.append(line.strip())

    pos_ID = []
    neg_ID = []
    for gene in pos_ID_o:
        if neg_ID_o.count(gene)==0:
                pos_ID.append(gene)
    for gene in neg_ID_o:
        if pos_ID_o.count(gene)==0:
                neg_ID.append(gene)
    
    pos_keys = access[0].keys()
    neg_keys = access[1].keys()
   
    total_access_pos = {}
    total_access_neg = {}
    record_pos = 0
    record_neg = 0
    for gene in pos_ID:
        gene = gene.strip()
        if pos_keys.count(gene)>0:
                total_access_pos[gene] = access[0][gene]
        else:
                total_access_pos[gene] = 0
                record_pos = record_pos + 1
    for gene in neg_ID:
        gene = gene.strip()
        if neg_keys.count(gene)>0:
                total_access_neg[gene] = access[1][gene]
        else:
                total_access_neg[gene] = 0
                record_neg = record_neg +1
   
    result = AUC(total_access_pos.values(), total_access_neg.values())
    print "pos not have the sites" + str(len(pos_ID) - len(access[0]))
    print "neg not have the sites" + str(len(neg_ID) - len(access[1]))
	
  
    for item in set(pos_ID + neg_ID):
        if pos_ID.count(item) + neg_ID.count(item) > 1:
            print item
    print('AUCROC: ' + str(result[0]))
    print ('Ranksum P value: ' + str(result[1]))
    #return access
    #return [access[0].values(), access[1].values()]
    return [result,[total_access_pos,total_access_neg]]




def PATS_with_zero(con, pos_file_name,neg_file_name, seq_file_name, RNAplfold_direct, RNAplfold_file_p):


    motifs = required_motif(consensus_to_motifs(con))

    site_len = len(motifs[0])
    RNAplfold_file = RNAplfold_direct  + str(site_len) + '_whole.txt'
    loc_dic = loc_whole_region(motifs, seq_file_name)
    filtered_access = filtered_accessibility(RNAplfold_file, pos_file_name, neg_file_name)
    #access = sum_region_accessibility(filtered_access, loc_dic)
    
    #RNAplfold_file_p = '/home/morrislab/xiaoli/posttrans_network/fly/RNAplfold_result/W80L40u1_whole.txt'
   
    filtered_access_p = filtered_accessibility(RNAplfold_file_p, pos_file_name, neg_file_name)
   	
    pos_access = {}
    pos_gene_pool = filtered_access[0].keys()
    for gene in pos_gene_pool:
        sum = 0 # sum for all accessibilty
        num = 0
        if loc_dic.has_key(gene) and len(loc_dic[gene]) > 0:
            for loc in loc_dic[gene]:
                sum = sum + filtered_access[0][gene][loc]*(1-filtered_access_p[0][gene][loc-1])
                num = num + 1
            pos_access[gene] = sum

    neg_access = {}
    neg_gene_pool = filtered_access[1].keys()
    for gene in neg_gene_pool:
        sum = 0
        num = 0
        if loc_dic.has_key(gene) and len(loc_dic[gene]) > 0:
            for loc in loc_dic[gene]:
                sum = sum + filtered_access[1][gene][loc]*(1-filtered_access_p[1][gene][loc-1])
                num = num + 1
            neg_access[gene] = sum
 
   
    seq_ID = []
    pos_ID = []
    neg_ID = []
    seq_f = open(seq_file_name,'r')
    for line in seq_f:
        if line.startswith('>'):
            seq_ID.append(line[1:].strip())
    pos_f = open(pos_file_name,'r')
    for line in pos_f:
        if seq_ID.count(line.strip()) > 0 and pos_ID.count(line.strip())==0:
            pos_ID.append(line.strip())
    neg_f = open(neg_file_name,'r')
    for line in neg_f:
        if seq_ID.count(line.strip()) > 0 and neg_ID.count(line.strip())==0:
            neg_ID.append(line.strip())

    #total_access_pos = access[0].values()
    #total_access_neg = access[1].values()
    pos_keys = pos_access.keys()
    neg_keys = neg_access.keys()
    pos_keys.sort()
    neg_keys.sort()

    #total_access_pos = access[0].values()
    #total_access_neg = access[1].values()
    total_access_pos = []
    total_access_neg = []
    for gene in pos_keys:
        total_access_pos.append(pos_access[gene])
    for gene in neg_keys:
        total_access_neg.append(neg_access[gene])

    for zero_pos in range(len(pos_ID) - len(pos_access)):
        total_access_pos.append(0)
    for zero_neg in range(len(neg_ID) - len(neg_access)):
        total_access_neg.append(0)

    result = AUC(total_access_pos, total_access_neg)

    for item in set(pos_ID + neg_ID):
        if pos_ID.count(item) + neg_ID.count(item) > 1:
            print item
    print('AUCROC: ' + str(result[0]))
    print ('Ranksum P value: ' + str(result[1]))
    #return access
    #return [access[0].values(), access[1].values()]
    return [result,[total_access_pos,total_access_neg]]





def TS(con, pos_file_name,neg_file_name, seq_file_name, RNAplfold_direct):
    motifs = required_motif(consensus_to_motifs(con))
    site_len = len(motifs[0])
    #RNAplfold_file = RNAplfold_direct + str(site_len) + '_whole.txt'
    loc_dic = loc_whole_region(motifs, seq_file_name)
    filtered_access = filtered_accessibility(RNAplfold_direct, pos_file_name, neg_file_name, site_len)
    access = sum_site_num(filtered_access, loc_dic)
    result = AUC(access[0].values(), access[1].values())
    
    print('AUCROC: ' + str(result[0]))
    print ('Ranksum P value: ' + str(result[1]))
    print ('positive number: ' + str(len(access[0])))
    print ('negative number: ' + str(len(access[1])))
        
    return result
    #return [access[0].values(), access[1].values()]
    
def TS_with_zero(con, pos_file_name,neg_file_name, seq_file_name, RNAplfold_direct):
    motifs = required_motif(consensus_to_motifs(con))
    site_len = len(motifs[0])
    #RNAplfold_file = RNAplfold_direct  + str(site_len) + '_whole.txt'
    loc_dic = loc_whole_region(motifs, seq_file_name)
    filtered_access = filtered_accessibility(RNAplfold_direct, pos_file_name, neg_file_name, site_len)
    access = sum_site_num(filtered_access, loc_dic)
    
    
    seq_ID = []
    pos_ID_o = []
    neg_ID_o = []
    seq_f = open(seq_file_name,'r')
    for line in seq_f:
	if line.startswith('>'):
	    seq_ID.append(line[1:].strip())
    pos_f = open(pos_file_name,'r')
    for line in pos_f:
	if seq_ID.count(line.strip()) > 0 and pos_ID_o.count(line.strip()) == 0: 
	    pos_ID_o.append(line.strip())
    neg_f = open(neg_file_name,'r')
    for line in neg_f:
	if seq_ID.count(line.strip()) > 0 and neg_ID_o.count(line.strip()) == 0:
	    neg_ID_o.append(line.strip())
    pos_ID = []
    neg_ID = []
    for gene in pos_ID_o:
	if neg_ID_o.count(gene) ==0:
		pos_ID.append(gene)
    for gene in neg_ID_o:
	if pos_ID_o.count(gene) ==0:
		neg_ID.append(gene)    

    pos_keys = access[0].keys()
    neg_keys = access[1].keys()
    pos_keys.sort()
    neg_keys.sort()
  
    #total_access_pos = access[0].values()
    #total_access_neg = access[1].values()
    total_access_pos = []
    total_access_neg = []
    for gene in pos_keys:
	total_access_pos.append(access[0][gene])
    for gene in neg_keys:
        total_access_neg.append(access[1][gene])


    for zero_pos in range(len(pos_ID) - len(access[0])):
	total_access_pos.append(0)
    for zero_neg in range(len(neg_ID) - len(access[1])):
	total_access_neg.append(0)
    
    result = AUC(total_access_pos, total_access_neg)
    
    for item in set(pos_ID + neg_ID):
	if pos_ID.count(item) + neg_ID.count(item) > 1:
	    print item
    print('AUCROC: ' + str(result[0]))
    print ('Ranksum P value: ' + str(result[1]))
    #return [access[-1].values(), access[1].values()]
    #return [result,[total_access_pos,total_access_neg]]
    return [result,[total_access_pos, total_access_neg],[pos_keys, neg_keys]] 

    
'''final_out_name = '/Users/xiaoli/Desktop/SSCR_run/result/access/final_result.txt'
detailed_final_out_name = '/Users/xiaoli/Desktop/SSCR_run/result/access/access_detailed_result.txt'
model = 'access_seq'
#final_out_name= '/Users/xiao/Desktop/RNA_stability_name/result/access.txt'
pos_file_name= '/Users/xiaoli/Desktop/SSCR_run/pos_gene2.txt'
neg_file_name = '/Users/xiaoli/Desktop/SSCR_run/neg_gene2.txt'
seq_file_name = '/Users/xiaoli/Desktop/SSCR_run/seq.txt'
RNAplfold_direct = '/Users/xiaoli/Desktop/SSCR_run/RNAplfold_result/'
'''

#model = 'access_seq'
#final_out_name: name of the final output
# detailed_final_out_name: name of the final output (the detailed version)
# pos_file_name, neg_file_name: file names of the positive and negative files. Make sure each line has one gene name, also make sure it is consitent in the RNAplfold result
# RNAplfold_direct: directory name of the RNAplfold results

Motif_discovery(final_out_name, detailed_final_out_name, model, pos_file_name, neg_file_name, seq_file_name, RNAplfold_direct)   


