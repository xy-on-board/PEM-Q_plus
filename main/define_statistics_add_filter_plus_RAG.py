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
    define_statistics    <basename> <genome> <cutsite> <transloc_range> <primer> <primer_chr> <primer_strand> <adapter>

Options:
-h --help               Show this screen.
-v --version            Show version.

Final results are generated.

"""


import os
from time import time
from docopt import docopt
import pandas as pd
    
def statistics_add_filter(basename=None,genome=None,cutsite=None,transloc_range=None,primer=None,primer_chr=None,primer_strand=None,adapter=None):
    
    
    print("Processing junction filter and statistics...")
    os.system('mkdir results')
    
    # ~~~~ Germline ~~~~ #
    
    germ_file = "unique/" + basename + "_Germline_final.tab"
    germ = pd.read_csv(germ_file, sep = '\t', index_col=False, low_memory=False)
    germ.to_csv("results/" + basename + "_Germline.tab", header = True, sep = '\t', index=False)
    germ_count = germ["Qname"].count()
    
    # ~~~~ Deletion  ~~~~ #
    
    deletion_file = "unique/" + basename+"_Deletion.tab"
    deletion_all = pd.read_csv(deletion_file, sep = '\t', index_col=False, low_memory=False)
    # print(deletion_all['Qname'].count())
    # Deletion misprime filter
    condition_misprime = deletion_all['Bait_end'] - deletion_all['Bait_start'] > len(primer) + 10
    deletion_all = deletion_all[condition_misprime]
    # print(deletion_all['Qname'].count())
    
    # Deletion (within transloc_range) 
    deletion_all['Microhomolog_len']=deletion_all['Microhomolog'].map(lambda x: 0 if pd.isna(x) else len(str(x)))
    
    if primer_strand == "+":
        #deletion length
        deletion_all['deletion_length'] = deletion_all['Prey_start'] - deletion_all['Bait_end'] + deletion_all['Microhomolog_len'] - 1
        condition_del = deletion_all['deletion_length'] <= int(transloc_range)
    else:
        deletion_all['deletion_length'] = deletion_all['Bait_start'] - deletion_all['Prey_end'] + deletion_all['Microhomolog_len'] - 1
        condition_del = deletion_all['deletion_length'] <= int(transloc_range)
    deletion = deletion_all[condition_del]
    deletion.to_csv("results/" + basename+"_Deletion.tab", header = True, sep = '\t', index=False)
    deletion_count = deletion["Qname"].count()
    
    small_del = deletion[deletion['deletion_length']<=100]['Qname'].count()
    large_del = deletion[deletion['deletion_length']>100]['Qname'].count()
    # print(small_del,large_del)
    del_stats_file = open("results/" + basename+"_del_len_statistics.txt","w")
    del_stats_file.write("small_del(<=100bp)"+"\t"+str(small_del)+"\n")
    del_stats_file.write("large_del(>100bp)"+"\t"+str(large_del)+"\n")
    del_stats_file.close()
    
    deletion_length_freq = deletion['deletion_length'].value_counts()
    deletion_length_freq = deletion_length_freq.sort_index()
    deletion_length_freq.to_csv("results/" + basename+"_deletion_length.txt", header = True, sep = '\t', index=True)
    
    
    # ~~~~ Insertion ~~~~ #
    
    insertion_file = "unique/" + basename+"_Insertion.tab"
    insertion_all = pd.read_csv(insertion_file, sep = '\t', index_col=False, low_memory=False)
    # insertion misprime filter
    condition_misprime = insertion_all['Bait_end'] - insertion_all['Bait_start'] > len(primer) + 10
    insertion_all = insertion_all[condition_misprime]

    insertion_all['Microhomolog_len'] = insertion_all['Microhomolog'].map(lambda x: 0 if pd.isna(x) else len(x))
    if primer_strand == "+":
        #deletion length
        insertion_all['deletion_length'] = insertion_all['Prey_start'] - insertion_all['Bait_end'] + insertion_all['Microhomolog_len'] - 1
        condition_del = insertion_all['deletion_length'] <= int(transloc_range)
    else:
        insertion_all['deletion_length'] = insertion_all['Bait_start'] - insertion_all['Prey_end'] + insertion_all['Microhomolog_len'] - 1
        condition_del = insertion_all['deletion_length'] <= int(transloc_range)
    insertion = insertion_all[condition_del]
    insertion.to_csv("results/" + basename + "_Insertion.tab", header = True, sep = '\t', index=False)
    insertion_count = insertion["Qname"].count()
    
    insertion['insertion_length']=insertion['Insertion'].map(lambda x: 0 if pd.isna(x) else len(str(x)))
    insertion_length_freq = insertion['insertion_length'].value_counts() 
    insertion_length_freq = insertion_length_freq.sort_index()
    insertion_length_freq.to_csv("results/" + basename+"_insertion_length.txt", header = True, sep = '\t', index=True)
    
    one_bp_inser = insertion[insertion['insertion_length']==1]['Qname'].count()
    small_inser = insertion[insertion['insertion_length']<20]['Qname'].count()
    large_inser = insertion[insertion['insertion_length']>=20]['Qname'].count()
    # print(small_inser,large_inser)
    inser_stats_file = open("results/" + basename+"_inser_len_statistics.txt","w")
    inser_stats_file.write("1_bp"+"\t"+str(one_bp_inser)+"\n")
    inser_stats_file.write("small_inser(<20bp)"+"\t"+str(small_inser)+"\n")
    inser_stats_file.write("large_inser(>=20bp)"+"\t"+str(large_inser)+"\n")
    inser_stats_file.close()
    
    # ~~~~ Inversion ~~~~ #
    
    inversion_file = "unique/" + basename+"_inversion.tab"
    inversion = pd.read_csv(inversion_file, sep = '\t', index_col=False, low_memory=False)
    # inversion misprime filter
    condition_misprime = inversion['Bait_end'] - inversion['Bait_start'] > len(primer) + 10
    inversion = inversion[condition_misprime]
    
    inversion['Microhomolog_len'] = inversion['Microhomolog'].map(lambda x: 0 if pd.isna(x) else len(x))
    inversion.to_csv("results/" + basename + "_Inversion.tab", header = True, sep = '\t', index=False)
    inversion_count = inversion["Qname"].count()
    if primer_strand == "+":
        condition_close_inver = abs(inversion['Prey_start'] - inversion['Bait_end']) + inversion['Microhomolog_len'] - 1 < int(transloc_range)
    else:
        condition_close_inver = abs(inversion['Bait_start'] - inversion['Prey_end']) + inversion['Microhomolog_len'] - 1 < int(transloc_range)
    close_inver = inversion[condition_close_inver]
    close_inver.to_csv("results/" + basename+"_close.inver.tab", header = True, sep = '\t', index=False)
    close_inver_count = close_inver["Qname"].count()
    
    # ~~~~ Translocation ~~~~ #
    
    # Part1. Intra_translocations_from_deletions
    
    if primer_strand == "+":
        condition_transloc = (deletion_all['Prey_start'] - deletion_all['Bait_end'] + deletion_all['Microhomolog_len'] - 1) > int(transloc_range)
    else:
        condition_transloc = (deletion_all['Bait_start'] - deletion_all['Prey_end'] + deletion_all['Microhomolog_len'] - 1) > int(transloc_range)
    transloc_del = deletion_all[condition_transloc]

    # transloc_del.to_csv("results/" + basename+"_intra_Translocation.del.tab", header = True, sep = '\t', index=False)
    # transloc_del_count = transloc_del["Qname"].count()
    
    # Part2. Intra_transloctions_from_inversion
    
    if primer_strand == "+":
        condition_inver = abs(inversion['Prey_start'] - inversion['Bait_end']) + inversion['Microhomolog_len'] - 1 > int(transloc_range)
    else:
        condition_inver = abs(inversion['Bait_start'] - inversion['Prey_end']) + inversion['Microhomolog_len'] - 1 > int(transloc_range)
    transloc_inver = inversion[condition_inver]
    # transloc_inver.to_csv("results/" + basename+"_intra_Translocation.inver.tab", header = True, sep = '\t', index=False)
    # transloc_inver_count = transloc_inver["Qname"].count()

    # Part3. Intra_transloctions_from_insertion # new in plus
    
    if primer_strand == "+":
        condition_inser = abs(insertion['Prey_start'] - insertion['Bait_end']) + insertion['Microhomolog_len'] - 1 > int(transloc_range)
    else:
        condition_inser = abs(insertion['Bait_start'] - insertion['Prey_end']) + insertion['Microhomolog_len'] - 1 > int(transloc_range)
    transloc_inser = insertion[condition_inser]
    transloc_inser.to_csv("results/" + basename+"_intra_Translocation.inser.tab", header = True, sep = '\t', index=False)
    transloc_inser_count = transloc_inser["Qname"].count()
    
    # Part4. Inter_translocations
    
    intersv_file = "unique/" + basename+"_interSV.tab"
    intersv = pd.read_csv(intersv_file, sep = '\t', index_col=False, low_memory=False)
    
    # intersv misprime filter
    condition_misprime = intersv['Bait_end'] - intersv['Bait_start'] > len(primer) + 10
    intersv = intersv[condition_misprime]

    # intersv.to_csv("results/" + basename + "_inter_Translocation.tab", header = True, sep = '\t', index=False)
    # intersv_count = intersv["Qname"].count()
    
    # Total translocations
    translocations = pd.concat([transloc_del, transloc_inver])
    translocations = pd.concat([translocations, transloc_inser])
    translocations = pd.concat([translocations, intersv])

    translocation_count = translocations["Qname"].count()
    translocations.to_csv("results/" + basename + "_Translocation.tab", header = True, sep = '\t', index=False)
    
    # Translocation filter
    cmd = "DSB_filter_woCondition_plus_RAG.py {} {}".format(basename, "_Translocation")
    print(cmd)
    os.system(cmd)

    translocations_dsb = pd.read_csv("results/" + basename + "_Translocation_dsb.tab", sep = '\t', index_col=False, low_memory=False)
    dsb_count = translocations_dsb["Qname"].count()
    
    transloc_del = pd.merge(transloc_del, translocations_dsb, on='Qname', how='inner')
    transloc_del.to_csv("results/" + basename+"_intra_Translocation.del.tab", header = True, sep = '\t', index=False)
    transloc_del_count = transloc_del["Qname"].count()

    transloc_inser = pd.merge(transloc_inser, translocations_dsb, on='Qname', how='inner')
    transloc_inser_count = transloc_inser["Qname"].count()

    transloc_inver = pd.merge(transloc_inver, translocations_dsb, on='Qname', how='inner')
    transloc_inver.to_csv("results/" + basename+"_intra_Translocation.inver.tab", header = True, sep = '\t', index=False)
    transloc_inver_count = transloc_inver["Qname"].count()

    inversv = pd.merge(intersv, translocations_dsb, on='Qname', how='inner')
    intersv.to_csv("results/" + basename + "_inter_Translocation.tab", header = True, sep = '\t', index=False)
    intersv_count = inversv["Qname"].count()

    total_events = germ_count + deletion_count + insertion_count + close_inver_count + dsb_count
    edit_events = deletion_count + insertion_count + dsb_count
    edit_effi = edit_events/total_events
    
    stats_file = open("results/" + basename+"_statistics.txt","w")
    stats_file.write("NoJunction"+"\t"+str(germ_count)+"\n")
    stats_file.write("Deletion"+"\t"+str(deletion_count)+"\n")
    stats_file.write("Insertion"+"\t"+str(insertion_count)+"\n")
    stats_file.write("Close_inversion.del"+"\t"+str(close_inver_count)+"\n")
    stats_file.write("Intra_Translocation.del"+"\t"+str(transloc_del_count)+"\n")
    stats_file.write("Intra_Translocation.inser"+"\t"+str(transloc_inser_count)+"\n")
    stats_file.write("Intra_Translocation.inver"+"\t"+str(transloc_inver_count)+"\n")
    stats_file.write("Inter_Translocation"+"\t"+str(intersv_count)+"\n")
    stats_file.write("Translocation dsb"+"\t"+str(dsb_count)+"\n")
    stats_file.write("Editing Events"+"\t"+str(edit_events)+"\n")
    stats_file.write("Total Events"+"\t"+str(total_events)+"\n")
    stats_file.write("Editing Efficiency"+"\t"+str(edit_effi)+"\n")
    stats_file.close()
    
    # Merge editing events
    deletion.drop(columns=['Microhomolog_len','deletion_length'], inplace=True)
    insertion.drop(columns=['Microhomolog_len','insertion_length'], inplace=True)
    close_inver.drop(columns=['Microhomolog_len'], inplace=True)
    #concatenate dataframes
    frames = [deletion, insertion, translocations_dsb]
    merge_df = pd.concat(frames, sort=False)
    merge_df = merge_df.sort_values(by=['Rname','Junction','Insertion','Barcode'], ascending=True)
    merge_df.reset_index(drop=True, inplace=True)
    merge_df.to_csv("results/" + basename + "_Editing_events.tab", header = True, sep = '\t', index=False)
    
    #Plot dot plot
    cmd = "TranslocPlot.R results/{}_Editing_events.tab results/{}_Editing_events_dot_plot.pdf \
           binsize=2000000 assembly={} plotshape=octogon strand=0".format(basename,basename,genome)
    print(cmd)
    os.system(cmd)
    
    #generate html file
    cmd = "TranslocHTMLReads_PEMQ.pl results/{}_Editing_events.tab results/{}_Editing_events.html \
           --primer {} -adapter {}".format(basename,basename,primer,adapter)
    print(cmd)
    os.system(cmd)
    
def main():

    start_time = time()

    args = docopt(__doc__,version='define_statistics_add_filter 1.0')

    kwargs = {'basename':args['<basename>'],'genome':args['<genome>'],'cutsite':args['<cutsite>'],'transloc_range':args['<transloc_range>'],\
              'primer':args['<primer>'],'primer_chr':args['<primer_chr>'],'primer_strand':args['<primer_strand>'],'adapter':args['<adapter>']}
    print('[PEM-Q] basename: ' + str(kwargs['basename']))
    print('[PEM-Q] genome: ' + str(kwargs['genome']))
    print('[PEM-Q] cutsite: ' + str(kwargs['cutsite']))
    print('[PEM-Q] transloc_range: ' + str(kwargs['transloc_range']))
    print('[PEM-Q] primer: ' + str(kwargs['primer']))
    print('[PEM-Q] primer_chr: ' + str(kwargs['primer_chr']))
    print('[PEM-Q] primer_strand: ' + str(kwargs['primer_strand']))
    print('[PEM-Q] adapter: ' + str(kwargs['adapter']))
    
    
    statistics_add_filter(**kwargs)

    print("\Define_statistics_add_filter.py done in {}s".format(round(time()-start_time, 3)))

if __name__ == '__main__':
    main()
