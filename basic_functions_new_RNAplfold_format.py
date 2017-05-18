import random
import numpy
import scipy.stats
import itertools  # knie: added for required_motifs()


# no N shifting adding on right/left
# Considering shorter motif (compared to the input canadiates)

def consensus_to_motifs(consensus):
    '''input = 'UGUH', output = ['U','G','U','ACU']'''

    seq_dic = {}
    seq_dic['A'] = 'A'
    seq_dic['C'] = 'C'
    seq_dic['U'] = 'U'
    seq_dic['G'] = 'G'
    seq_dic['H'] = 'ACU'
    seq_dic['W'] = 'AU'
    seq_dic['N'] = 'ACGU'
    seq_dic['M'] = 'AC'
    seq_dic['Y'] = 'UC'
    seq_dic['K'] = 'GU'
    seq_dic['B'] = 'CGU'
    seq_dic['D'] = 'AGU'
    seq_dic['R'] = 'GA'
    seq_dic['V'] = 'ACG'
    seq_dic['S'] = 'GC'
    motif = []
    for nt in consensus:
        motif.append(seq_dic[nt])
    return motif


def trans_str_to_number(ACGU_str):
    
    result_number = 0
    for i in range(len(ACGU_str)):
        if ACGU_str[i] == 'A':
            result_number = result_number + pow(4,1)
        elif ACGU_str[i] == 'C':
            result_number = result_number + pow(4,2)
        elif ACGU_str[i] == 'G':
            result_number = result_number + pow(4,3)
        elif ACGU_str[i] == 'U':
            result_number = result_number + pow(4,4)
    return result_number


def poss_to_degernate_motif(poss):
    
    '''Input = ['U','ACU'],output = 'UH' '''
    
    number_to_str = {}
    number_to_str[pow(4,1)] = 'A'
    number_to_str[pow(4,2)] = 'C'
    number_to_str[pow(4,3)] = 'G'
    number_to_str[pow(4,4)] = 'U'
    number_to_str[pow(4,1)+pow(4,2) + pow(4,4)] = 'H'
    number_to_str[pow(4,1) + pow(4,4)] = 'W'
    number_to_str[pow(4,1)+pow(4,2) + pow(4,3) + pow(4,4)] = 'N'
    number_to_str[pow(4,1)+pow(4,2)] = 'M'
    number_to_str[pow(4,2)+pow(4,4)] = 'Y'
    number_to_str[pow(4,3)+pow(4,4)] = 'K'
    number_to_str[pow(4,4)+pow(4,2) + pow(4,3)] = 'B'
    number_to_str[pow(4,1)+pow(4,4) + pow(4,3)] = 'D'
    number_to_str[pow(4,1)+pow(4,3)] = 'R'
    number_to_str[pow(4,1)+pow(4,2) + pow(4,3)] = 'V'
    number_to_str[pow(4,2) + pow(4,3)] = 'S'
    
    motif_str = ''
    for item in poss:
        tran_num = trans_str_to_number(item)
        motif_str = motif_str + number_to_str[tran_num]
    return motif_str


# def required_motif(poss):
#
#     '''input = ['U','ACU'], output = ['UU','UA','UC']'''
#     motifs = []
#     motif = ''
#     poss_num = 1
#     for item in poss:
#         poss_num = len(item) * poss_num
#     while len(motifs) < poss_num:
#         motif = ''
#         for i in range(len(poss)):
#             motif = motif + random.choice(poss[i])
#         if motifs.count(motif) == 0:
#             motifs.append(motif)
#     return motifs


def required_motif(poss):

    # knie: modified 2017/05 for faster calculation, imported itertools module

    '''input = ['U','ACU'], output = ['UU','UA','UC']'''
    poss = [p.strip() for p in poss if p.strip()]
    motifs = []
    if poss:
        motifs = [nt for nt in poss[0]]
        for p in poss[1:]:
            motifs = [''.join(i) for i in itertools.product(motifs, p)]
    return motifs


def more_degenerate_nt(input_nt):
    
    all_nt = ['ACU', 'AU', 'ACGU', 'AC', 'UC', 'GU', \
              'CGU', 'AGU', 'GA', 'ACG', 'GC', 'A','C','G','U']
    result = []
    for item in all_nt:
        if len(item) > len(input_nt):
            total = 0
            for i in range(len(input_nt)):
                if item.count(input_nt[i]) > 0:
                    total = total + 1
            if total == len(input_nt):
                result.append(item)
            
    return result


def more_degenerate_poss(input_poss):
    
    result = []
    for i in range(len(input_poss)):
        
        poss_nt = more_degenerate_nt(input_poss[i])
        for item in poss_nt:
            temp = input_poss[:i]
            temp.append(item)
            for follow in input_poss[(i+1):]:
                temp.append(follow)
            result.append(temp)
            #print temp
    return result


def shift(input_motif_poss, leftORright):
    
    result = []
    all_nt = ['ACU', 'AU', 'AC', 'UC', 'GU', \
              'CGU', 'AGU', 'GA', 'ACG', 'GC', 'A','C','G','U']
    if leftORright == 'Left':
        for nt in all_nt:
            temp = []
            
            temp.append(nt)
            for item in input_motif_poss:
                temp.append(item)
            result.append(temp)
       
    elif leftORright == 'Right':
        for nt in all_nt:
            temp = input_motif_poss[:]
            temp.append(nt)
            result.append(temp)
            
    else:
        print 'Wrong shift instruction!'
    return result


def shift_less(input_motif_poss, leftORright):
    
    result = []
    all_nt = ['ACU', 'AU', 'ACGU', 'AC', 'UC', 'GU', \
              'CGU', 'AGU', 'GA', 'ACG', 'GC', 'A','C','G','U']
    if leftORright == 'Left':
        result.append(input_motif_poss[1:])

    elif leftORright == 'Right':
        result.append(input_motif_poss[:-1])
            
    else:
        print 'Wrong shift instruction!'
    return result


def loc_whole_region(motifs, seq_file_name):
    seq_file = open(seq_file_name,'r')
    seq = seq_file.readlines()
    seq_file.close()
    loc_dic = {}
    for line in seq:
        if line.startswith('>'):
            if line.find('|') != -1:
                name = line[1:line.index('|')]
            else:
                name = line[1:].strip()
        else:
            total = 0
            loc = []
            for motif in motifs:
                start = 0
                current = 0
                i = 0
                for iii in range(line.count(motif)):
                #while line.find(motif, start) != -1:
                    current = line.index(motif, start)
                    loc.append(current)
                    start = current + 1
                    i = i + 1
            loc.sort()
            loc_dic[name] = []
            for item in loc:
                loc_dic[name].append(item)
    return loc_dic


def utr3_loc_to_whole_transcript_loc(utr3_loc, whole_seq_file, utr3_seq_file):

    seq_whole_utr3 = {}
    utr3_pos_in_whole = {}
    file_whole = open(whole_seq_file,'r')
    for line in file_whole:
        if line.startswith('>'):
            name = line[1:].strip()
        else:
            seq_whole_utr3[name] = line.strip()
    file_whole.close()

    file_utr3 = open(utr3_seq_file,'r')
    for line in file_utr3:
        if line.startswith('>'):
            name = line[1:].strip()
        else:
            seq = line.strip()
            if seq_whole_utr3.has_key(name):
                position = seq_whole_utr3[name].index(seq)
                utr3_pos_in_whole[name] = position
    file_utr3.close()

    genes = utr3_loc.keys()
    on_whole_loc = {}
    for gene in genes:
        locations = utr3_loc[gene]
        on_whole_loc[gene] = []
        position = utr3_pos_in_whole[gene]
        for loc in locations:
            new_pos = position + int(loc)
            on_whole_loc[gene].append(new_pos)
    return on_whole_loc
                                                                

# taken from previous codes

def loc_in_utr_orf(motifs, seq_file_name, loc_file_name):
    
    '''Again, motifs is the list with all possible motifs, seq_file_name is the 
    name of file containing fasta transcripts sequence. loc_file_name is the
    name of output_file. Basically this function will tell you the location of 
    all possible motifs for each gene and also the location will be classified 
    as 5'UTR, ORF or 3'UTR. The thing need to keep in mind is that this 
    transcripts sequence for yeast is composed of 200nt 5'UTR, ORF and 200nt 
    3'UTRs'''

    '''Again, there is an option that not write the result into file, just use
    '' to replace loc_file_name'''
    
    seq_file = open(seq_file_name, 'r')
    seq = seq_file.readlines()
    seq_file.close()
    
    loc_dic = {}
    for line in seq:
        if line.startswith('>'):
            name = line.strip()
           
        else:
            total = 0
            orf_len = len(line.strip()) - 400
            loc = []
            for motif in motifs:
                start = 0
                current = 0
                for i in range(line.count(motif)):
                    current = line.index(motif, start)
                    loc.append(current)
                    start = current + 1
                    
            loc.sort()
            utr_5_coll = []
            orf_coll = []
            utr_3_coll = []
            for item in loc:
                if item <= 200:
                    utr_5_coll.append(item)
                elif 200 < item and item <= (orf_len + 200):
                    orf_coll.append(item)
                elif item > (orf_len + 200):
                    utr_3_coll.append(item)
            loc_dic[name[1:]] = [utr_5_coll, orf_coll, utr_3_coll]
    
    if loc_file_name != '':
        loc_file = open(loc_file_name, 'w')
        gene_pool = loc_dic.keys()
        for gene in gene_pool:
            if len(loc_dic[gene][0]) + len(loc_dic[gene][1]) + len(loc_dic[gene][2]) > 0:
                loc_file.write(gene + '\t\t')
                for region in loc_dic[gene]:
                    for loc in region:
                        loc_file.write(str(loc) + '\t')
                    loc_file.write('\t')
                loc_file.write('\n')
    
        loc_file.close()
    return loc_dic


def loc_in_utr_orf_region(motifs, seq_file_name,  utrORorf):
     
     '''same as the last function except output file will focus on the region users 
     defined. Here utrORorf can be 'utr_5', 'utr_3' or 'orf'only.'''
     
     loc_dic = loc_in_utr_orf(motifs, seq_file_name, '')
     result = {}
     gene_pool = loc_dic.keys()
     
     if utrORorf == 'UTR_5':
         for gene in gene_pool:
             result[gene] = loc_dic[gene][0]
     elif utrORorf == 'ORF':
         for gene in gene_pool:
             result[gene] = loc_dic[gene][1]
     elif utrORorf == 'UTR_3':
         for gene in gene_pool:
             result[gene] = loc_dic[gene][2]
     elif utrORorf == 'ALL':
         for gene in gene_pool:
             result[gene] = loc_dic[gene][0] + loc_dic[gene][1] + loc_dic[gene][2]
     elif utrORorf == 'ORFandUTR_3':
         for gene in gene_pool:
             result[gene] = loc_dic[gene][1] + loc_dic[gene][2]
     elif utrORorf == 'UTR_5andUTR_3':
         for gene in gene_pool:
             result[gene] = loc_dic[gene][0] + loc_dic[gene][2]
     elif utrORorf == 'UTR_5andORF':
         for gene in gene_pool:
             result[gene] = loc_dic[gene][0] + loc_dic[gene][1]
     else:
         print 'Wrong input of utrORorf!'
             
     return result


def filtered_accessibility_old(RNAplfold_file, pos_list, neg_list):

    '''RNAplfold_file is the name of proper RNAplfold result, Basically, 
    this gives a dict for accessibility for all the possible genes 
    with the required motifs'''
    #print pos
    #print type(pos)
    #print type([])

    if type(pos_list) != type([]):
        file_pos = open(pos_list, 'r')
        pos = file_pos.readlines() # there is '\n' at the end of each name
        file_pos.close()
    
        file_neg = open(neg_list, 'r')
        neg = file_neg.readlines()
        file_neg.close()
    else:
        pos = pos_list[:]
        neg = neg_list[:]
        
    RNAplfold_result = open(RNAplfold_file,'r')
    RNAplfold = RNAplfold_result.readlines()
    RNAplfold_result.close()
    
    pos_access = {}
    neg_access = {}
    exist = 'No'
    for line in RNAplfold:
        if line.startswith('>'):
            
            name = line[1:-1]
            if pos.count(name + '\n') > 0:
                exist = 'pos'
                pos_access[name.strip()] = []
            elif neg.count(name + '\n') > 0:
                exist = 'neg'
                neg_access[name.strip()] = []
            else:
                exist = 'No'
            
        elif not line.startswith('#'):
            if exist == 'pos':
                temp = line.split()
                pos_access[name.strip()].append(float(temp[1]))
            elif exist == 'neg':
                temp = line.split()
                neg_access[name.strip()].append(float(temp[1]))
    filtered_access = [pos_access, neg_access, len(pos), len(neg)]

    return filtered_access


def filtered_accessibility(RNAplfold_dir, pos_list, neg_list, mer_length):

    # modified 10/25/16 to accept the RNAplfold format change
    if type(pos_list) != type([]):
        file_pos = open(pos_list, 'r')
        pos = file_pos.readlines() # there is '\n' at the end of each name
        file_pos.close()

        file_neg = open(neg_list, 'r')
        neg = file_neg.readlines()
        file_neg.close()
    else:
        pos = pos_list[:]
        neg = neg_list[:]
    
    #all_RNAplfold_files = os.listdir(RNAplfold_dir)
    #pos_overlap_list = set(pos) 
    pos_access = {}
    for pos_gene in pos:
        pos_gene = pos_gene.strip()
        # RNAplfold_result = open(RNAplfold_file, 'r')
        access_list = []
        try:
            acc_file_in = open(RNAplfold_dir + pos_gene + '_lunp', 'r')
        except:
            print pos_gene + ' does not exist!'
            continue

        for line in acc_file_in:
            if line.count('#')==0:
                line = line.strip()
                temp = line.split('\t')
                #acc_s = temp[mer_length]
                if len(temp)> mer_length:
                    acc_s = temp[mer_length]
                    if acc_s!='NA':
                        acc = float(acc_s)
                        access_list.append(acc)
                else:
                    print 'longer k mer from RNAplfold settting is required!! Error here!!'
                    access_list = ['NA']
        pos_access[pos_gene] = access_list

    neg_access = {}
    for neg_gene in neg:
        neg_gene = neg_gene.strip()
        #RNAplfold_result = open(RNAplfold_file, 'r')
        access_list = []
        try: 
            acc_file_in = open(RNAplfold_dir + neg_gene +'_lunp', 'r')
        except:
            print neg_gene + ' does not exist!'
            continue

        for line in acc_file_in:
            if line.count('#')==0:
                line = line.strip()
                temp = line.split('\t')
                #acc_s = temp[mer_length]
                if len(temp)> mer_length:
                    acc_s = temp[mer_length]
                    if acc_s!='NA':
                        acc = float(acc_s)
                        access_list.append(acc)
                else:
                    print 'longer k mer from RNAplfold settting is required!! Error here!!'
                    access_list = ['NA']
        neg_access[neg_gene] = access_list

    filtered_access =  [pos_access, neg_access, len(pos), len(neg)]
    # if 'FBtr0335497' in filtered_access[0].keys():
    #     print('Filtered access FBtr0335497:', len(filtered_access[0]['FBtr0335497']))
        # print(filtered_access[0]['FBtr0335497'])

    return filtered_access


def sum_region_accessibility_site(filtered_access, loc_dic):

    '''filtered_access is the result of filtered_accessibility. loc_dic is the                                                                                                                                                                     
    dictionary with all the locations'''

    gene_access = {}
    gene_sites = {}
    gene_pool = filtered_access.keys()
    for gene in gene_pool:
        sum = 0 # sum for all accessibilty                                                                                                                                                                                                         
        num = 0
        if loc_dic.has_key(gene) and len(loc_dic[gene]) > 0:
            for loc in loc_dic[gene]:
                # print(gene, loc, len(loc_dic[gene]), len(filtered_access[gene]))
                sum = sum + filtered_access[gene][loc]
                num = num + 1
            gene_access[gene] = sum
            gene_sites[gene] = num

    return [gene_access, gene_sites]


def sum_region_accessibility(filtered_access, loc_dic):
    
    '''filtered_access is the result of filtered_accessibility. loc_dic is the
    dictionary with all the locations'''
    
    pos_access = {}
    pos_gene_pool = filtered_access[0].keys()
    for gene in pos_gene_pool:
        sum = 0 # sum for all accessibilty
        num = 0
        if loc_dic.has_key(gene) and len(loc_dic[gene]) > 0:
            for loc in loc_dic[gene]:
                # print(gene, loc, len(filtered_access[0][gene]))
                sum = sum + filtered_access[0][gene][loc]
                num = num + 1
            pos_access[gene] = sum
            
    neg_access = {}
    neg_gene_pool = filtered_access[1].keys()
    for gene in neg_gene_pool:
        sum = 0
        num = 0
        if loc_dic.has_key(gene) and len(loc_dic[gene]) > 0:
            for loc in loc_dic[gene]:
                sum = sum + filtered_access[1][gene][loc]
                num = num + 1
            neg_access[gene] = sum
           
    return [pos_access, neg_access, len(pos_gene_pool), len(neg_gene_pool)]


def sum_site_num(filtered_access, loc_dic):

    pos_site_num = {}
    pos_gene_pool = filtered_access[0].keys()
    for gene in pos_gene_pool:
        num = 0
        if loc_dic.has_key(gene) and len(loc_dic[gene]) > 0:
            for loc in loc_dic[gene]:
                num = num + 1
            pos_site_num[gene] = num
            
    neg_site_num = {}
    neg_gene_pool = filtered_access[1].keys()
    for gene in neg_gene_pool:
        num = 0
        if loc_dic.has_key(gene) and len(loc_dic[gene]) > 0:
            for loc in loc_dic[gene]:
                num = num + 1
            neg_site_num[gene] = num
           
    return [pos_site_num, neg_site_num, len(pos_gene_pool), len(neg_gene_pool)]


def AUC(x, y):
        
    """Calculates the AUC, rank sums statistic on the provided scores and
    returns the result.
    x --> positives
    
    y --> negatives
    
    
    Returns: AUC, z-statistic, two-tailed p-value, 
    """
    x,y = map(numpy.asarray, (x, y))
    n1 = len(x)
    n2 = len(y)
    alldata = numpy.concatenate((x,y))
    ranked = scipy.stats.rankdata(alldata)
    x = ranked[:n1]
    y = ranked[n1:]
    s = numpy.sum(x,axis=0) # sum of positive ranks

    # AUC = 1 - s/(n1*n2) + (n1+1)/float(2*n2)
    if n1*n2 != 0:
        AUC = s/(n1*n2) - (n1+1)/float(2*n2) # AUC
    else:
        AUC = 'NaN'
    expected = n1*(n1+n2+1) / 2.0

    # z = (s - expected)/numpy.sqrt(n1*n2*(n1+n2+1)/12.0)
    '''modification starts'''
    z = (s - expected) / numpy.sqrt(n1*n2*(n1+n2+1)/12.0) # Zscore
    T = scipy.stats.tiecorrect(ranked)
    if T == 0:
        # raise ValueError, 'All numbers are identical in ranksums'
        prob = 'tie'
    z /= numpy.sqrt(T)
    '''modification ends'''

    # prob = 2*(1.0 - scipy.stats.zprob(abs(z))) # ranksum p value
    # prob = 2*(min(scipy.stats.zprob(z),scipy.stats.zprob(-z)))
    prob = 2*(min(scipy.stats.stats.special.ndtr(z), scipy.stats.stats.special.ndtr(-z)))
    return [AUC, prob, n1, n2]
