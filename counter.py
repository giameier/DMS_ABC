"""
@author: gmeier
"""


import pysam
import numpy as np
import matplotlib.pyplot as plt
import copy
import mate_generator

CODON_TABLE = {
        'TTT':'F', 'TCT':'S', 'TAT':'Y', 'TGT':'C',
        'TTC':'F', 'TCC':'S', 'TAC':'Y', 'TGC':'C',
        'TTA':'L', 'TCA':'S', 'TAA':'*', 'TGA':'*',
        'TTG':'L', 'TCG':'S', 'TAG':'*', 'TGG':'W',
        'CTT':'L', 'CCT':'P', 'CAT':'H', 'CGT':'R',
        'CTC':'L', 'CCC':'P', 'CAC':'H', 'CGC':'R',
        'CTA':'L', 'CCA':'P', 'CAA':'Q', 'CGA':'R',
        'CTG':'L', 'CCG':'P', 'CAG':'Q', 'CGG':'R',
        'ATT':'I', 'ACT':'T', 'AAT':'N', 'AGT':'S',
        'ATC':'I', 'ACC':'T', 'AAC':'N', 'AGC':'S',
        'ATA':'I', 'ACA':'T', 'AAA':'K', 'AGA':'R',
        'ATG':'M', 'ACG':'T', 'AAG':'K', 'AGG':'R',
        'GTT':'V', 'GCT':'A', 'GAT':'D', 'GGT':'G',
        'GTC':'V', 'GCC':'A', 'GAC':'D', 'GGC':'G',
        'GTA':'V', 'GCA':'A', 'GAA':'E', 'GGA':'G',
        'GTG':'V', 'GCG':'A', 'GAG':'E', 'GGG':'G'
    }
COUNT_TABLE = {
        'TTT':0, 'TCT':0, 'TAT':0, 'TGT':0,
        'TTC':0, 'TCC':0, 'TAC':0, 'TGC':0,
        'TTA':0, 'TCA':0, 'TAA':0, 'TGA':0,
        'TTG':0, 'TCG':0, 'TAG':0, 'TGG':0,
        'CTT':0, 'CCT':0, 'CAT':0, 'CGT':0,
        'CTC':0, 'CCC':0, 'CAC':0, 'CGC':0,
        'CTA':0, 'CCA':0, 'CAA':0, 'CGA':0,
        'CTG':0, 'CCG':0, 'CAG':0, 'CGG':0,
        'ATT':0, 'ACT':0, 'AAT':0, 'AGT':0,
        'ATC':0, 'ACC':0, 'AAC':0, 'AGC':0,
        'ATA':0, 'ACA':0, 'AAA':0, 'AGA':0,
        'ATG':0, 'ACG':0, 'AAG':0, 'AGG':0,
        'GTT':0, 'GCT':0, 'GAT':0, 'GGT':0,
        'GTC':0, 'GCC':0, 'GAC':0, 'GGC':0,
        'GTA':0, 'GCA':0, 'GAA':0, 'GGA':0,
        'GTG':0, 'GCG':0, 'GAG':0, 'GGG':0,
        'reads_filtered':0, 'Triplett_not_identical_in_mates':0,'reads_with_additional_mutation_1':0,'reads_with_additional_mutation_2':0,'reads_with_additional_mutation_3':0,'reads_with_additional_mutation_1_rev':0,'reads_with_additional_mutation_2_rev':0,'reads_with_additional_mutation_3_rev':0, 'reads_with_index_error':0
    }
AA_CODES = {
        'Ala' : 'A', 'A' : 'Ala',
        'Arg' : 'R', 'R' : 'Arg',
        'Asn' : 'N', 'N' : 'Asn',
        'Asp' : 'D', 'D' : 'Asp',
        'Cys' : 'C', 'C' : 'Cys',
        'Glu' : 'E', 'E' : 'Glu',
        'Gln' : 'Q', 'Q' : 'Gln',
        'Gly' : 'G', 'G' : 'Gly',
        'His' : 'H', 'H' : 'His',
        'Ile' : 'I', 'I' : 'Ile',
        'Leu' : 'L', 'L' : 'Leu',
        'Lys' : 'K', 'K' : 'Lys',
        'Met' : 'M', 'M' : 'Met',
        'Phe' : 'F', 'F' : 'Phe',
        'Pro' : 'P', 'P' : 'Pro',
        'Ser' : 'S', 'S' : 'Ser',
        'Thr' : 'T', 'T' : 'Thr',
        'Trp' : 'W', 'W' : 'Trp',
        'Tyr' : 'Y', 'Y' : 'Tyr',
        'Val' : 'V', 'V' : 'Val',
        'Ter' : '*', '*' : 'Ter',
        '???' : '?', '?' : '???'
}



#store dicts of counts for each positions in a dict with position as key and dict as value


def pileup_counter(bam_file,starts,quality_filter,reference_seq):

    pileup_dict=copy.copy(COUNT_TABLE)
    for read1,read2 in mate_generator.read_pair_generator(bam_file,start=starts, end=starts+1):
        if not 'I' in read1.cigarstring and read2.cigarstring and not 'D' in read1.cigarstring and read2.cigarstring:
            triplett1=read1.query_sequence[starts-read1.reference_start:starts-read1.reference_start+3]
            triplett2=read2.query_sequence[starts-read2.reference_start:starts-read2.reference_start+3]
            triplett_qualities1=read1.query_qualities[starts-read1.reference_start:starts-read1.reference_start+3]
            triplett_qualities2=read2.query_qualities[starts-read2.reference_start:starts-read2.reference_start+3]
            if any(x<quality_filter for x in triplett_qualities1 or triplett_qualities2):
                pileup_dict['reads_filtered']= pileup_dict['reads_filtered']+1

            elif triplett1==triplett2:
                
# with this part active one would allow secondary mutations in wt reads                
#                if triplett1==reference_seq[starts:starts+3]:
#                    if triplett1 in pileup_dict.keys():
#                        pileup_dict[triplett1]=pileup_dict[triplett1]+1
#                    else:
#                        pileup_dict[triplett1]=1
#                    samfile_counted.write(read1)
#                    samfile_counted.write(read2) 

#check overlapping sequence for additional mismatches. iterates over all triplett on read and checks for missmatches to query sequence
            
                    break_check=True
                    
                    
                    for i in range(starts+3,read1.reference_end-3):

                        wt_base=reference_seq[i]
                        
                        try:
                            r1_base=read1.query_sequence[i-read1.reference_start]
                        except IndexError:
                            pileup_dict['reads_with_index_error']=pileup_dict['reads_with_index_error']+1
                        try:
                            r2_base=read2.query_sequence[i-read2.reference_start]
                        except IndexError:
                            pileup_dict['reads_with_index_error']=pileup_dict['reads_with_index_error']+1
                            break 
                        q_r1_base=read1.query_qualities[i-read1.reference_start]
                        q_r2_base=read2.query_qualities[i-read2.reference_start]
#                        if i< read1.reference_end():
#                        print(r1_base)
#                        print(wt_base)
                            
                            #check if read1 is wt and read2 is wt
                        if r1_base==r2_base==wt_base:
                                #if so continue with next iteration
                                continue
                            
                        elif r1_base==wt_base and r2_base!=wt_base:#read1 is wt but read2 is not wt
                                #check if low quality in read2 at position i
                            if q_r2_base < 30:
                                continue
                            else:
                                pileup_dict['reads_with_additional_mutation_1']= pileup_dict['reads_with_additional_mutation_1']+1
                                break_check=False
                                break
                                       
                        elif r1_base!=wt_base and r2_base==wt_base: #read1 is not wt but read 2 is wt
                        #check if low quality in read 1
                            if q_r1_base < 30:
                                continue
                            else:
                                pileup_dict['reads_with_additional_mutation_2']= pileup_dict['reads_with_additional_mutation_2']+1
                                break_check=False
                                break
                                       
                        elif r1_base!=wt_base and r2_base!=wt_base: #read1 is not wt and read 2 is not wt
                            if any(x>=30 for x in [q_r1_base,q_r2_base]):
                                pileup_dict['reads_with_additional_mutation_3']= pileup_dict['reads_with_additional_mutation_3']+1
                                break_check=False
       
                                break        
                        else:print('error_druing_base_checking')
                    
                    #if no additional mutation (means no breakcheck) check if other side of read also has no additional mutation
                    if break_check:
                        
                        for i in range(starts-1,read1.reference_start+3,-1):
                            wt_base=reference_seq[i]
                            try:
                                r1_base=read1.query_sequence[i-read1.reference_start]
                            except:
                                pileup_dict['reads_with_index_error']=pileup_dict['reads_with_index_error']+1
                            try:
                                r2_base=read2.query_sequence[i-read2.reference_start]
                            except:
                                pileup_dict['reads_with_index_error']=pileup_dict['reads_with_index_error']+1                            
                            q_r1_base=read1.query_qualities[i-read1.reference_start]
                            q_r2_base=read2.query_qualities[i-read2.reference_start]
                            if r1_base==r2_base==wt_base:
                                    #if so continue with next iteration
                                    continue
                                
                            elif r1_base==wt_base and r2_base!=wt_base:#read1 is wt but read2 is not wt
                                    #check if low quality in read2 at position i
                                if q_r2_base < 30:
                                    continue
                                else:
                                    pileup_dict['reads_with_additional_mutation_1_rev']= pileup_dict['reads_with_additional_mutation_1_rev']+1
                                    break_check=False
                                    break
                                           
                            elif r1_base!=wt_base and r2_base==wt_base: #read1 is not wt but read 2 is wt
                            #check if low quality in read 1
                                if q_r1_base < 30:
                                    continue
                                else:
                                    pileup_dict['reads_with_additional_mutation_2_rev']= pileup_dict['reads_with_additional_mutation_2_rev']+1
                                    break_check=False
                                    break
                                           
                            elif r1_base!=wt_base and r2_base!=wt_base: #read1 is not wt and read 2 is not wt
                                if any(x>=30 for x in [q_r1_base,q_r2_base]):
                                    pileup_dict['reads_with_additional_mutation_3_rev']= pileup_dict['reads_with_additional_mutation_3_rev']+1
                                    break_check=False
                                    break

                            else:print('error_druing_base_checking')
                            #if forward and reverse check ended with no break->add count to pileupdict
                    if break_check:
                        if triplett1 in pileup_dict.keys():
                            pileup_dict[triplett1]=pileup_dict[triplett1]+1
                        else:
                            pileup_dict[triplett1]=1                                                        
                    
            elif triplett1!=triplett2:
                pileup_dict['Triplett_not_identical_in_mates']= pileup_dict['Triplett_not_identical_in_mates']+1

    return pileup_dict             


def create_triplett_count(samfile,quality_filter,positions,reference_seq):
    count_dict={}
    for starts in positions:
        count_dict[int(starts)]=pileup_counter(samfile,int(starts),quality_filter,reference_seq)
        print('processed position'+str(starts))
        
    return count_dict
    
def count_mutants(data_dict,curent_file):
    input_file_directory=data_dict['outputdir']+'/'+curent_file[:-4]+'_codontruncated.bam'
    output_file_directory=input_file_directory+'triplet_count.txt'   
    positions=data_dict['position_list']
    reference_seq=data_dict['reference_sequence']
    print(positions)
    samfile = pysam.AlignmentFile(input_file_directory,'rb')
    output_txt=open(output_file_directory,'w')
    count_dict=create_triplett_count(samfile,30,positions,reference_seq)     
    
    output_txt.write('position'+'\t')             
    keylist=list(count_dict.keys())
    keylist.sort()
    keylist_codon=list(COUNT_TABLE.keys())
    keylist_codon.sort()
    sorted_count_dict=[]
    for key in keylist:
        sorted_count_dict.append((key,count_dict[key]))
        
    for key in keylist_codon:    
        output_txt.write(str(key)+'\t')
    
        
    for i in sorted_count_dict:
        output_txt.write('\n'+str(i[0])+'\t')
        dict_i=i[1]
        for x in keylist_codon:
            output_txt.write(str(dict_i[x])+'\t')    
    
    samfile.close()     
    output_txt.close()