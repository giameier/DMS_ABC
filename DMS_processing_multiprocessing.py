"""
@author: gmeier

DMS processing module

"""
import json
import Codon_truncation_bam_overlapp
import counter
import create_count_file
import multiprocessing
from functools import partial
from contextlib import contextmanager


@contextmanager
def poolcontext(*args, **kwargs):
    pool = multiprocessing.Pool(*args, **kwargs)
    yield pool
    pool.terminate()

def DMS_processing(i,data_dict):
    """
    Analyses curent file i from data_dict
    
    """
    
    print('analyzing files')
    Codon_truncation_bam_overlapp.codon_truncation(data_dict,i)
    counter.count_mutants(data_dict,i)
    create_count_file.make_HDF5(data_dict,i)    


#running the dms data analysis using config data from jsonfile
def run_analysis(json_file_directory):
    with open(json_file_directory,'r') as jsonfile:
        data_dict=json.load(jsonfile)
    
    #processes defines the amount of processes beeing executed in paralell
    #Nr of processes can be increased if CPUs are available
    with poolcontext(processes=3) as pool:
        pool.map(partial(DMS_processing, data_dict=data_dict), data_dict['data_files'])
    
    
    print('process finished')
    