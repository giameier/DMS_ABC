"""

@author: gmeier

Mate sequence trimmer for DMS of ABC tranporter

Cutting the sequences of a pair of reads to the overlapping region with starts and ends only in frame of the AA sequence.

Trims forwar and reverse reads to overlapp only and removes non overlapp sequences

"""



def modify_cigar(read1,read1_bases_totrim_for):
    """
    Function to trim a read from the front
    First checks if reads are softclipped and removes these sequences
    Then a trimms reads and modifies their cigartuples/cigarstrings/query_qualities and query_sequence
    
    
    read1: read to be trimmed
    read1_bases_totrim_for: number of bases which have to be trimmed from the front
    """
    
    remaining_bases_to_trim=read1_bases_totrim_for
    read1.reference_start=read1.reference_start+read1_bases_totrim_for

#check if first bases are soft_clipped
    if read1.cigartuples[0][0]==4:
        read1_tuple=read1.cigartuples
        insert_lenght=read1.cigartuples[0][1]
        q1=read1.query_qualities
            
        read1.query_sequence=read1.query_sequence[insert_lenght:]
        read1.query_qualities=q1[insert_lenght:]
                      
            
        if read1.cigartuples[0][1]<10:
            read1.cigarstring=read1.cigarstring[2:]
        elif read1.cigartuples[0][1]<100:
            read1.cigarstring=read1.cigarstring[3:]
        elif read1.cigartuples[0][1]<1000:
            read1.cigarstring=read1.cigarstring[4:]    
            
        read1_tuple=read1_tuple[1:] 
        read1.cigartuples=read1_tuple

            
        bases_in_first_string=0


    #loop executed as long as there are still remaining bases to be trimmed
    while remaining_bases_to_trim>0:
        read1_tuple=read1.cigartuples

        #check if first bases are insertion and remove these from read   
        if read1.cigartuples[0][0]==1:
            insert_lenght=read1.cigartuples[0][1]
            q1=read1.query_qualities
            
            read1.query_sequence=read1.query_sequence[insert_lenght:]
            read1.query_qualities=q1[insert_lenght:]
                      
            
            if read1.cigartuples[0][1]<10:
                read1.cigarstring=read1.cigarstring[2:]
            elif read1.cigartuples[0][1]<100:
                read1.cigarstring=read1.cigarstring[3:]
            elif read1.cigartuples[0][1]<1000:
                read1.cigarstring=read1.cigarstring[4:]    
            
            read1_tuple=read1_tuple[1:] 
            read1.cigartuples=read1_tuple

            
            bases_in_first_string=0
            
        #check if first bases are deletions and remove these from read
        elif read1.cigartuples[0][0]==2:
            
            bases_in_first_string=read1.cigartuples[0][1]
            
            if read1.cigartuples[0][1]<10:
                read1.cigarstring=read1.cigarstring[2:]
            elif read1.cigartuples[0][1]<100:
                read1.cigarstring=read1.cigarstring[3:]
            elif read1.cigartuples[0][1]<1000:
                read1.cigarstring=read1.cigarstring[4:]    
            
            read1_tuple=read1_tuple[1:] 
            read1.cigartuples=read1_tuple
            
        #check if first tuple is smaller than the remaining bases to trim->the tuple ist removed completely           
        elif read1.cigartuples[0][1]<= remaining_bases_to_trim:
            bases_in_first_string=read1.cigartuples[0][1]
            q1=read1.query_qualities
            read1.query_sequence=read1.query_sequence[bases_in_first_string:]
            read1.query_qualities=q1[bases_in_first_string:]

            if read1.cigartuples[0][1]<10:
                read1.cigarstring=read1.cigarstring[2:]
                
            elif read1.cigartuples[0][1]<100:
                read1.cigarstring=read1.cigarstring[3:]
            elif read1.cigartuples[0][1]<1000:
                read1.cigarstring=read1.cigarstring[4:]    
            
            read1_tuple=read1_tuple[1:]            
            read1.cigartuples=read1_tuple
            
            
        else:
            bases_in_first_string=read1.cigartuples[0][1]
            q1=read1.query_qualities
            read1.query_sequence=read1.query_sequence[remaining_bases_to_trim:]
            read1.query_qualities=q1[remaining_bases_to_trim:]
            
            
            if read1.cigartuples[0][1]<10:
                read1.cigarstring=str(read1.cigartuples[0][1]-remaining_bases_to_trim)+read1.cigarstring[1:]
            elif read1.cigartuples[0][1]<100:
                read1.cigarstring=str(read1.cigartuples[0][1]-remaining_bases_to_trim)+read1.cigarstring[2:]
            elif read1.cigartuples[0][1]<1000:
                read1.cigarstring=str(read1.cigartuples[0][1]-remaining_bases_to_trim)+read1.cigarstring[3:]
            
            cigar_code=read1_tuple[0][0]
            bases_in_tuple=read1_tuple[0][1]
            
            
            read1_tuple=read1_tuple[1:]
            
            read1_tuple.insert(0,(cigar_code,bases_in_tuple-remaining_bases_to_trim))
                       
            read1.cigartuples=read1_tuple
            

        remaining_bases_to_trim=remaining_bases_to_trim-bases_in_first_string
    else:
        return read1


#Function to trim a read from the back                
def modify_cigar_rev(read1,read1_bases_totrim_rev):
    """
    Function to trim a read from the front
    First checks if reads are softclipped and removes these sequences
    Then a trimms reads and modifies their cigartuples/cigarstrings/query_qualities and query_sequence
    
    
    read1: read to be trimmed
    read1_bases_totrim_for: number of bases which have to be trimmed from the front
    
    """
    
    remaining_bases_to_trim=read1_bases_totrim_rev
    read1_reference_end=read1.reference_end-read1_bases_totrim_rev
    
    #check if first bases are soft_clipped
    if read1.cigartuples[-1][0]==4:
        read1_tuple=read1.cigartuples
        insert_lenght=read1.cigartuples[-1][1]
        q1=read1.query_qualities
            
        read1.query_sequence=read1.query_sequence[:-insert_lenght]
        read1.query_qualities=q1[:-insert_lenght]
                      
            
        if read1.cigartuples[-1][1]<10:
            read1.cigarstring=read1.cigarstring[:-2]
        elif read1.cigartuples[-1][1]<100:
            read1.cigarstring=read1.cigarstring[:-3]
        elif read1.cigartuples[-1][1]<1000:
            read1.cigarstring=read1.cigarstring[:-4]    
            
        read1_tuple=read1_tuple[:-1] 
        read1.cigartuples=read1_tuple

            
        bases_in_first_string=0

            
    #loop executed as long as there are still remaining bases to be trimmed
    while remaining_bases_to_trim>0:

        read1_tuple=read1.cigartuples


        #check if first bases are insertion
        if read1.cigartuples[-1][0]==1:
            insert_lenght=read1.cigartuples[-1][1]
            q1=read1.query_qualities
            
            read1.query_sequence=read1.query_sequence[:-insert_lenght]
            read1.query_qualities=q1[:-insert_lenght]
                      
            
            if read1.cigartuples[-1][1]<10:
                read1.cigarstring=read1.cigarstring[:-2]
            elif read1.cigartuples[-1][1]<100:
                read1.cigarstring=read1.cigarstring[:-3]
            elif read1.cigartuples[-1][1]<1000:
                read1.cigarstring=read1.cigarstring[:-4]    
            
            read1_tuple=read1_tuple[:-1] 
            read1.cigartuples=read1_tuple

            
            bases_in_first_string=0
            
        #check if first bases are deletions
        elif read1.cigartuples[-1][0]==2:
            bases_in_first_string=read1.cigartuples[-1][1]
            if read1.cigartuples[-1][1]<10:
                read1.cigarstring=read1.cigarstring[:-2]
            elif read1.cigartuples[-1][1]<100:
                read1.cigarstring=read1.cigarstring[:-3]
            elif read1.cigartuples[-1][1]<1000:
                read1.cigarstring=read1.cigarstring[:-4]    
            
            read1_tuple=read1_tuple[:-1] 
            read1.cigartuples=read1_tuple
            
        #check if first tule is smaller than the remaining bases to trim->the the tuple ist removed completely           
        elif read1.cigartuples[-1][1]<= remaining_bases_to_trim:
            bases_in_first_string=read1.cigartuples[-1][1]
            q1=read1.query_qualities
            read1.query_sequence=read1.query_sequence[:-bases_in_first_string]
            read1.query_qualities=q1[:-bases_in_first_string]
            
            
            
          
            if read1.cigartuples[-1][1]<10:
                read1.cigarstring=read1.cigarstring[:-2]
                
            elif read1.cigartuples[0][1]<100:
                read1.cigarstring=read1.cigarstring[:-3]
            elif read1.cigartuples[0][1]<1000:
                read1.cigarstring=read1.cigarstring[:-4]    
            
            read1_tuple=read1_tuple[:-1]
            
            read1.cigartuples=read1_tuple
            
            
        else:
            bases_in_first_string=read1.cigartuples[-1][1]
            q1=read1.query_qualities
            read1.query_sequence=read1.query_sequence[:-remaining_bases_to_trim]
            read1.query_qualities=q1[:-remaining_bases_to_trim]
            
            
            if read1.cigartuples[-1][1]<10:
                read1.cigarstring=read1.cigarstring[-2]+str(read1.cigartuples[-1][1]-remaining_bases_to_trim)+read1.cigarstring[-1]
            elif read1.cigartuples[0][1]<100:
                read1.cigarstring=read1.cigarstring[-3]+str(read1.cigartuples[-1][1]-remaining_bases_to_trim)+read1.cigarstring[-1]
            elif read1.cigartuples[0][1]<1000:
                read1.cigarstring=read1.cigarstring[-4]+str(read1.cigartuples[-1][1]-remaining_bases_to_trim)+read1.cigarstring[-1]
            
            cigar_code=read1_tuple[-1][0]
            bases_in_tuple=read1_tuple[-1][1]
            
            
            read1_tuple=read1_tuple[:-1]
            
            read1_tuple.append((cigar_code,bases_in_tuple-remaining_bases_to_trim))
                       
            read1.cigartuples=read1_tuple
            

        remaining_bases_to_trim=remaining_bases_to_trim-bases_in_first_string

    else:
        return read1       
        
        
        
#    else:
        
#        print('gone to endtrim')
#        print(remaining_bases_to_trim)
#        print(read1.cigartuples)
#        print(str(read1.cigarstring)+ '\n\n')
def trim_front_rev(read1,read2,for_reference_cut,rev_reference_cut):
    """
    defines amount of bases to trim
        
    """
    
    read1_bases_totrim_for=for_reference_cut-read1.reference_start
    read2_bases_totrim_for=for_reference_cut-read2.reference_start
    read1_bases_totrim_rev=read1.reference_end-rev_reference_cut
    read2_bases_totrim_rev=read2.reference_end-rev_reference_cut
    
    if read1_bases_totrim_for<147 and read2_bases_totrim_for<147 and read1_bases_totrim_rev<147 and read2_bases_totrim_rev<147:
        #entry point for trimming function
        read1_trimmed_f=modify_cigar(read1,read1_bases_totrim_for)
        read2_trimmed_f=modify_cigar(read2,read2_bases_totrim_for)
        read1_trimmed_f_r=modify_cigar_rev(read1_trimmed_f,read1_bases_totrim_rev)
        read2_trimmed_f_r=modify_cigar_rev(read2_trimmed_f,read2_bases_totrim_rev)
        return read1_trimmed_f_r, read2_trimmed_f_r
    else:

        return([None,None])

def trim(read1,read2,frameshift_position,frameshift_offset):
    """
    This function defines necessery arguments for trimming functions.
    Overlapp start: position on gene were the overlapp of read1 and read2 starts
    Overlapp end: position on gene were the overlapp of read1 and read2 ends
    Initiats forward and reverse trimming of reads
    
    Returns: Trimmed reads
    """
    #Define start of overlapp
    if read1.reference_start < read2.reference_start:
        overlapp_start=read2.reference_start
    elif read1.reference_start >= read2.reference_start:
        overlapp_start=read1.reference_start
    #Define End of overlapp
    if read1.reference_end < read2.reference_end:
        overlapp_end=read1.reference_end
    elif read1.reference_end >= read2.reference_end:
        overlapp_end=read2.reference_end
    #Takes into account the frame shift and frameshif offset if read is on EfrD chain
    
    if overlapp_start > frameshift_position:
        shift_for=frameshift_offset-2
    else:
        shift_for=0
        
    if overlapp_end > frameshift_position:
        shift_rev=frameshift_offset-2
    else:
        shift_rev=0       
        
    #checks if read_overlap_start is in frame 
    if overlapp_start in range(0+shift_for,4000,3):
        overlapp_start=overlapp_start
    #checks if reads have to be trimmed one base to get to next codon start
    elif overlapp_start in range(2+shift_for,4000,3):
        overlapp_start=overlapp_start +1
    #checks if reads have to be trimmed two bases to get to next codon start
    elif overlapp_start in range(1+shift_for,4000,3):
        overlapp_start=overlapp_start +2
        
    #checks if read_overlap_end is in frame
    if overlapp_end in range(0+shift_rev,4000,3):
        overlapp_end=overlapp_end
    #checks if reads have to be trimmed one base to get to last codon end
    elif overlapp_end in range(2+shift_rev,4000,3):
        overlapp_end=overlapp_end-2
    #checks if reads have to be trimmed two bases to get to last codon end
    elif overlapp_end in range(1+shift_rev,4000,3):
        overlapp_end=overlapp_end-1
    
    #Check if read overlapp is at least 6 bp long
    # then start trimming function
    if overlapp_end-overlapp_start>6:
        read1_t,read2_t=trim_front_rev(read1,read2,overlapp_start,overlapp_end)
        
        return read1_t, read2_t
    else:
        return None, None