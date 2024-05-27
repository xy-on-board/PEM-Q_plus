#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#Author: Xinyi Liu
#Date: 2024

"""PEM-Q plus
    
    PEM-Q plus serves as an updated version of PEM-Q, which was originally launched by Mengzhu Liu in 2019.
    PEM-Q plus anaylzes PEM-seq data, which is a high-throughput sequencing method for detecting DNA double-strand breaks (DSBs) and their repair products.

    Copyright (C) 2024  Xinyi Liu

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

Author: Xinyi Liu
Contact: xinyi.liu@stu.pku.edu.cn

Usage:
    DSB_filter   <basename> <yourfile>

Options:
-h --help               Show this screen.
-v --version            Show version.
<de_method>             1:sequence; 2:bait start/end & prey start/end; 3:junction

A more stringent deduplication process is applied to translocation events.

Input file: transloc tab file / Output file: transloc clean tab file

"""

import os
from time import time
from docopt import docopt
import pandas as pd
import numpy as np
from difflib import get_close_matches
import Levenshtein

def find_similar_rows_highest(group):

    similar_rows = []

    for i in range(len(group) - 1):
        row1 = group.iloc[i]
        thresh1 = row1['edit_dist_thresh']
        
        for j in range(i + 1, len(group)):
            row2 = group.iloc[j]
            thresh2 = row2['edit_dist_thresh']

            thresh = max(thresh1, thresh2)
            
            if Levenshtein.distance(row1['random_barcode'], row2['random_barcode']) <= thresh:
                similar_rows.append(row2['Qname'])
                break
    
    return similar_rows

def DSB_filter(basename, yourfile):
    ### column to dedup : Bait_rname,Bait_strand,Bait_start,Bait_end,Prey_rname,Prey_strand,Prey_start,Prey_end,Rname,Strand,Junction,Sequence
    #total dataset
    transloc = pd.read_csv("results/" + basename + yourfile + ".tab",sep = '\t')
    transloc = transloc[~transloc["Rname"].str.contains("_")]
    print("Translocation input: " + str(transloc['Qname'].count()))
    Full = pd.read_csv("results/" + basename + yourfile + ".tab",sep = '\t')
    if 'Offset' not in transloc.columns:
        transloc.loc[:,'Offset'] = np.where(transloc['Strand'] == "+", transloc.Junction - transloc.Qstart, transloc.Junction + transloc.Qlen - transloc.Qend + 1)
    if 'Bait_junction' not in transloc.columns:
        transloc.loc[:,'Bait_junction'] = np.where(transloc['Bait_strand'] == "+", transloc.Bait_end, transloc.Bait_start)
    transloc['Seq_len']=transloc['Sequence'].map(lambda x: 0 if x is np.nan else len(x))

    #reads to be judge 
    de_list = ['Rname','Strand','Offset','Prey_start','Prey_end','R2_barcode']
    
    #mark fist read of duplicated as true
    dup = transloc[transloc.duplicated(de_list, keep= False)]
    dup = dup.sort_values(by=['Barcode_freq','Barcode_len','edit_dist_thresh'], ascending=[False, False, True])
    reads = dup.drop_duplicates(subset=de_list, keep='first')
    print("Reads: " + str(reads['Qname'].count()))
    
    #start dedup
    finalset = transloc.copy()
    add = []
    for i in range(0,reads['Qname'].count()):
        
        Bait_junc = reads['Bait_junction'].iloc[i]
        Rname = reads['Rname'].iloc[i]
        Strand = reads['Strand'].iloc[i]
        Offset = reads['Offset'].iloc[i]
        Junction = reads['Junction'].iloc[i]
        Barcode = reads['Barcode'].iloc[i]
        R2_barcode = reads['R2_barcode'].iloc[i]

        de_condition1 = (transloc['Rname']==Rname) & \
                    (transloc['Strand']==Strand) & \
                    (transloc['Offset']==Offset) & \
                    (transloc['Bait_junction']>=Bait_junc-2) & \
                    (transloc['Bait_junction']<=Bait_junc+2)
        dataset = transloc.drop(transloc[de_condition1].index)
        count = dataset[(dataset['Rname'] == Rname) & (dataset['Junction'] >= Junction-5) & (dataset['Junction'] <= Junction+5)]
        condition = count['Qname'].count() <= 20

        if condition:
            de_condition2 = (finalset['Rname'] == Rname) & \
                        (finalset['Strand'] == Strand) & \
                        (finalset['Offset'] == Offset) & \
                        (finalset['Bait_junction'] >= Bait_junc - 2) & \
                        (finalset['Bait_junction'] <= Bait_junc + 2) & \
                        (finalset['R2_barcode'] == R2_barcode)
        
            finalset.drop(finalset[(de_condition2)].index, inplace=True)
            add.append(i)

    add_reads = reads.iloc[add]
    finalset = pd.concat([finalset, add_reads])
    finalset.reset_index(level=finalset.index.names, inplace=True)
    print("after dedup: " + str(finalset['Qname'].count()))
    
    #Barcode filter
    finalset2 = finalset.copy()
    i_list = []
    f_list = []
    for i in range(0,finalset2['Qname'].count()):
        if i in i_list:
            continue
        
        Rname = finalset2['Rname'].iloc[i]
        Strand = finalset2['Strand'].iloc[i]
        Junction = finalset2['Junction'].iloc[i]
        Barcode = finalset2['Barcode'].iloc[i]

        #add barcode dedup method
        de_barcode_method = (finalset2['Barcode']==Barcode) & \
                            (finalset2['Rname']==Rname) & \
                            (finalset2['Strand']==Strand) & \
                            (finalset2['Junction']>=Junction-5) & \
                            (finalset2['Junction']<=Junction+5)

        add_list = finalset2.index[de_barcode_method].tolist()
        f_list = list(set(add_list) - set([i])) # 去除当前i
        i_list = i_list + f_list
        
    finalset.drop(index=i_list, inplace=True)
    print("after barcode dedup: " + str(finalset['Qname'].count()))
    
    # dedup with RMB

    feature_list = ['Rname','Strand','Prey_start','Prey_end']

    dup = finalset[finalset.duplicated(feature_list, keep=False)]
    other = finalset[~finalset['Qname'].isin(dup['Qname'])]
    other = other.drop_duplicates(subset=feature_list + ['Barcode'])
    dup = dup.drop_duplicates(subset=feature_list + ['Barcode'])
    dup = dup.reset_index(drop=True)

    groups = dup.groupby(feature_list)

    result = []
 
    for _, group in groups:
            if len(group) == 1:
                continue
            else:
                group = group.sort_values(by=['Barcode_freq','Barcode_len','edit_dist_thresh'], ascending=[False, False, True])
                similar_rows = find_similar_rows_highest(group)
                result.extend(similar_rows)
    
    Qname_rm = pd.DataFrame({'Qname': result})


    dup = dup[~dup['Qname'].isin(Qname_rm['Qname'])]

    finalset = pd.concat([other, dup])
    print("after edit dist dedup: " + str(finalset['Qname'].count()))

    removed = Full[~Full['Qname'].isin(finalset['Qname'])]
    removed.to_csv("results/" + basename + yourfile + "_dsbremove.tab", header = True, sep = '\t', index=False)
    
    finalset = finalset.drop_duplicates(['Qname'],keep= 'first')
    finalset = finalset.sort_values(by=['Rname','Junction'], ascending=False)
    finalset.to_csv("results/" + basename + yourfile + "_dsb.tab", header = True, sep = '\t', index=False)
    print("output: " + str(finalset['Qname'].count()))
    

def main():
    start_time = time()
    
    args = docopt(__doc__,version='DSB_filter 1.0')
    
    kwargs = {'basename':args['<basename>'], 'yourfile':args['<yourfile>']}
    print('[PEM-Q] basename: ' + str(kwargs['basename']))
    
    DSB_filter(**kwargs)
    
    print("\nDSB_filter.py done in {}s".format(round(time()-start_time, 3)))
    
if __name__ == '__main__':
    main()