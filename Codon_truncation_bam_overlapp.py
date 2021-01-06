"""

@author: gmeier
"""
import pysam
import os
import mate_generator
import mate_sequence_trimmer    


def codon_truncation(data_dict,curent_file):
    """
    sequence trimmer function
    trims all sequences in curent file to codons_starts.
    
    
    """

    input_file_path=data_dict['inputdir']+'/'+curent_file
    output_directory=data_dict['outputdir']+'/'+curent_file[:-4]+'_codontruncated.bam'
    freameshift_position=data_dict['frameshift_position']
    frameshift_offset=data_dict['frameshift_offset']      

    #indexing inputfile
  
    os.system('samtools index '+ input_file_path.replace(' ','\ ')+' >'+ input_file_path.replace(' ','\ ')+'.bai')
    #open file for reading and new file for writing

    samfile = pysam.AlignmentFile(input_file_path,'rb')
    samfile_trimmed=pysam.AlignmentFile(output_directory, 'wb',template=samfile)
    
    


    
    unmapped_read_count=0
    softclip_count=0
    hardclip_count=0
    short_overlap_or_non=0
    
    #iterate over all read pairs
    for read1,read2 in mate_generator.read_pair_generator(samfile,start=None, end=None):
        
        #remove all unmapped reads
        if read1.flag==4 or read2.flag ==4 or read1.cigarstring==None or read2.cigarstring==None:
            unmapped_read_count=unmapped_read_count+1
        #remove all hardclipped reads
        elif 'H' in read1.cigarstring or 'H' in read2.cigarstring:
            hardclip_count=hardclip_count+1
        else:
            trimmed_read1,trimmed_read2=mate_sequence_trimmer.trim(read1,read2,freameshift_position,frameshift_offset)
            


            if trimmed_read1==None or trimmed_read2==None:
                short_overlap_or_non=short_overlap_or_non+1

            elif trimmed_read1!=None and trimmed_read2!=None:
                samfile_trimmed.write(trimmed_read1)
                samfile_trimmed.write(trimmed_read2)
    print('reads processed')
    
    samfile.close() 
    samfile_trimmed.close()
    
    #indexing output file
    os.system('samtools index '+ output_directory.replace(' ','\ ')+' >'+ output_directory.replace(' ','\ ')+'.bai')
    
    print(str(unmapped_read_count)+'unmapped reads, '+str(hardclip_count)+' hardclipped reads and ' +str(short_overlap_or_non)+' short overlap reads were removed')

