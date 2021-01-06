"""
@author: gmeier

#takes raw counts from a tsv file and returns a tsv count file in HGVS5 format
"""
import pandas as pd
import numpy as np
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


#: Conversions between single- and three-letter amino acid codes
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

def protein_HGSV(reference_position,reference_codon,codon):
    
    wt_AA=threeletter_AA(translate(reference_codon))
    mu_AA=threeletter_AA(translate(codon))
    if wt_AA==mu_AA:
        return ' (p.=)'
    else:
        return ' (p.'+wt_AA+str(int((reference_position-1)/3+2))+mu_AA+')'


def write_HGVS5_format_position(start_shift,reference_position,reference_codon,codon):
    hg5_string=''
    EfrC_shift=start_shift
    if reference_codon[0]!=codon[0]:
        hg5_string=hg5_string+'c.'+str(reference_position+EfrC_shift)+reference_codon[0]+'>'+codon[0]+protein_HGSV(reference_position,reference_codon,codon)
        if reference_codon[1]!=codon[1]:
            hg5_string=hg5_string+', c.'+str(reference_position+1+EfrC_shift)+reference_codon[1]+'>'+codon[1]+protein_HGSV(reference_position,reference_codon,codon)
            if reference_codon[2]!=codon[2]:
                hg5_string=hg5_string+', c.'+str(reference_position+2+EfrC_shift)+reference_codon[2]+'>'+codon[2]+protein_HGSV(reference_position,reference_codon,codon)
        elif reference_codon[2]!=codon[2]:
            hg5_string=hg5_string+', c.'+str(reference_position+2+EfrC_shift)+reference_codon[2]+'>'+codon[2]+protein_HGSV(reference_position,reference_codon,codon)
    elif reference_codon[1]!=codon[1]:
         hg5_string=hg5_string+'c.'+str(reference_position+1+EfrC_shift)+reference_codon[1]+'>'+codon[1]+protein_HGSV(reference_position,reference_codon,codon)       
         if reference_codon[2]!=codon[2]:
             hg5_string=hg5_string+', c.'+str(reference_position+2+EfrC_shift)+reference_codon[2]+'>'+codon[2]+protein_HGSV(reference_position,reference_codon,codon) 
    elif reference_codon[2]!=codon[2]:
        hg5_string=hg5_string+'c.'+str(reference_position+2+EfrC_shift)+reference_codon[2]+'>'+codon[2]+protein_HGSV(reference_position,reference_codon,codon)
    return hg5_string
    
def translate(codon):
    return CODON_TABLE[codon]

def threeletter_AA(one_letter_AA):
    return AA_CODES[one_letter_AA]

def read_csv_pandas(path_to_tsv):
    return pd.read_table(path_to_tsv,header=0)

def make_HDF5(data_dict,curent_file):
    reference_sequence=data_dict['reference_sequence']
    input_directory=data_dict['outputdir']+'/'+curent_file[:-4]+'_codontruncated.bam'+'triplet_count.txt'
    data_frame=read_csv_pandas(input_directory)
    output_HDF5_frame1=open(input_directory+'_readingframe_1_HDF5.csv','w')
    output_HDF5_frame2=open(input_directory+'_readingframe_2_HDF5.csv','w')
    column_list=list(data_frame)[1:65]
    frame_shift_position=data_dict['frameshift_position']
    frame_shift=data_dict['frameshift_offset']
    #create a list of reference positions from data frame
    reference_positions=list(data_frame.iloc[:,0])
    
    
  
#creates a header for all output files
    output_HDF5_frame1.write('\t'+'count'+'\n')
    output_HDF5_frame2.write('\t'+'count'+'\n')

    wt_counter1=0
    wt_counter2=0
#for each reference position in the reference position list determine referenceposition index, codon and 
    for ref in reference_positions:
#define parameters for one position to be analyzed

        ref_ind=reference_positions.index(ref) 
        
        ref_pos=data_frame.iloc[ref_ind,0]
        ref_codon=reference_sequence[ref_pos:ref_pos+3]
#check if reference psition is in frame 1 or 2
        if ref_pos<frame_shift_position:
            for i in column_list:
                if data_frame.loc[ref_ind,i] !=0:
                    if ref_codon == i:
                        wt_counter1=wt_counter1+data_frame.loc[ref_ind,i]
                    else:    
                        output_HDF5_frame1.write(write_HGVS5_format_position(3,ref_pos-2,ref_codon,i)+'\t'+str(data_frame.loc[ref_ind,i])+'\n')
    
        elif ref_pos>frame_shift_position:
            for i in column_list:
                if data_frame.loc[ref_ind,i] !=0:
                    if ref_codon == i:
                        wt_counter2=wt_counter2+data_frame.loc[ref_ind,i]
                    else:    
                        output_HDF5_frame2.write(write_HGVS5_format_position(1,ref_pos-frame_shift_position+frame_shift,ref_codon,i)+'\t'+str(data_frame.loc[ref_ind,i])+'\n')
                 
    output_HDF5_frame1.write('_wt'+'\t'+str(wt_counter1))        
    output_HDF5_frame2.write('_wt'+'\t'+str(wt_counter2))        
  
    output_HDF5_frame1.close()
    output_HDF5_frame1.close()