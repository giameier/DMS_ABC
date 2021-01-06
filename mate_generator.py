#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 19 07:26:41 2019

@author: gmeier

mate generator

"""

from collections import defaultdict
import pysam


def read_pair_generator(bam, start, end):
    """
    Find reads, store them in dict and return them once a pair is found.
    """
    read_dict = defaultdict(lambda: [None, None])
    for read in bam.fetch(reference='EfrCD_wt_sequence',start=start,end=end):
        if not read.is_proper_pair or read.is_secondary or read.is_supplementary:
            continue
        qname = read.query_name
        if qname not in read_dict:
            if read.is_read1:
                read_dict[qname][0] = read
            else:
                read_dict[qname][1] = read
        else:
            if read.is_read1:
                yield read, read_dict[qname][1]
            else:
                yield read_dict[qname][0], read
            del read_dict[qname]
