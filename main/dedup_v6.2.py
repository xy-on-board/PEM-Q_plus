#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#Author: Xinyi Liu
#Date: 2024

"""rmb_dedup

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
    dedup   <basename> <barcode_length>

Options:
-h --help               Show this screen.
-v --version            Show version.

Reads are deduplicated based on the barcode and R2 barcode.

Input files: raw files, barcode lists / Output files: deduplicated files

"""
import os
import pysam
from time import time
from docopt import docopt
import pandas as pd
import numpy as np

def dedup(basename=None, bl=None):

    barcode = pd.read_csv("barcode/"+basename+"_barcode_list.txt",sep = '\t',names = ["Qname", "Barcode", "R2_barcode"])
    sid = pd.read_csv("raw/"+basename+"_SID_all.tab",sep = '\t')
    sv = pd.read_csv("raw/"+basename+"_SV.tab",sep = '\t')
    substitution = pd.read_csv("raw/"+basename+"_Substitution.tab",sep = '\t')
    insertion = pd.read_csv("raw/"+basename+"_Insertion.tab",sep = '\t')
    deletion = pd.read_csv("raw/"+basename+"_Deletion.tab",sep = '\t')
    germline = pd.read_csv("raw/"+basename+"_Germline.tab",sep = '\t')
    transloc = pd.read_csv("transloc/"+basename+"_mut.tab",sep = '\t')
    indel = pd.read_csv("indel/"+basename+"_indel_mut.tab",sep = '\t')
    # all_insertion = pd.read_csv(basename+"_All_Insertions.tab",sep = '\t')
    germline_addup = pd.read_table("raw/" + basename + "_Germline_addup_indel.tab",sep = '\t')
    indel_cutoff = pd.read_table("raw/"+basename + "_indel_cutoff.tab",sep = '\t')
    insertion_sv = pd.read_csv("raw/"+basename+"_Insertions_inSV.tab",sep = '\t')

    if "Barcode" not in sid.columns or "R2_barcode" not in sid.columns:
        sid_merge = pd.merge(sid, barcode, on='Qname', how='inner')
        sid_merge.to_csv("raw/"+basename + "_SID_all.tab", header = True, sep = '\t', index=False)
    if "Barcode" not in sv.columns or "R2_barcode" not in sv.columns:
        sv_merge = pd.merge(sv, barcode, on='Qname', how='inner')
        sv_merge.to_csv("raw/"+basename + "_SV.tab", header = True, sep = '\t', index=False)
    if "Barcode" not in substitution.columns or "R2_barcode" not in substitution.columns:
        substitution_merge = pd.merge(substitution, barcode, on='Qname', how='inner')
        substitution_merge.to_csv("raw/"+basename + "_Substitution.tab", header = True, sep = '\t', index=False)
    else:
        substitution_merge = pd.read_csv("raw/"+basename+"_Substitution.tab",sep = '\t')
    if "Barcode" not in insertion.columns or "R2_barcode" not in insertion.columns:
        insertion_merge = pd.merge(insertion, barcode, on='Qname', how='inner')
        insertion_merge.to_csv("raw/"+basename + "_Insertion.tab", header = True, sep = '\t', index=False)
    if "Barcode" not in deletion.columns or "R2_barcode" not in deletion.columns:
        deletion_merge = pd.merge(deletion, barcode, on='Qname', how='inner')
        deletion_merge.to_csv("raw/"+basename + "_Deletion.tab", header = True, sep = '\t', index=False)
    if "Barcode" not in germline.columns or "R2_barcode" not in germline.columns:
        germline_merge = pd.merge(germline, barcode, on='Qname', how='inner')
        germline_merge.to_csv("raw/"+basename + "_Germline.tab", header = True, sep = '\t', index=False)
    else:
        germline_merge = pd.read_csv("raw/"+basename+"_Germline.tab",sep = '\t')
    if "Barcode" not in transloc.columns or "R2_barcode" not in transloc.columns:
        transloc_merge = pd.merge(transloc, barcode, on='Qname', how='inner')
        transloc_merge.to_csv("raw/"+basename + "_transloc.tab", header = True, sep = '\t', index=False)
    if "Barcode" not in indel.columns or "R2_barcode" not in indel.columns:
        indel_merge = pd.merge(indel, barcode, on='Qname', how='inner')
        indel_merge.to_csv("raw/"+basename + "_smallindel.tab", header = True, sep = '\t', index=False)
    if "Barcode" not in germline_addup.columns or "R2_barcode" not in germline_addup.columns:
        germline_addup_merge = pd.merge(germline_addup, barcode, on='Qname', how='inner')
        germline_addup_merge.to_csv("raw/"+basename + "_Germline_addup_indel.tab", header = True, sep = '\t', index=False)
    else:
        germline_addup_merge = pd.read_table("raw/" + basename + "_Germline_addup_indel.tab",sep = '\t')
    if "Barcode" not in indel_cutoff.columns or "R2_barcode" not in indel_cutoff.columns:
        indel_cutoff_merge = pd.merge(indel_cutoff, barcode, on='Qname', how='inner')
        indel_cutoff_merge.to_csv("raw/"+basename + "_smallindel_cutoff.tab", header = True, sep = '\t', index=False)
    if "Barcode" not in insertion_sv.columns or "R2_barcode" not in insertion_sv.columns:
        insertion_sv_merge = pd.merge(insertion_sv, barcode, on='Qname', how='inner')
        insertion_sv_merge.to_csv("raw/"+basename + "_Insertion_sv.tab", header = True, sep = '\t', index=False)

    os.system("mkdir unique")

    germline_dedup = germline_merge.drop_duplicates(['Cigar','Barcode','Position','R2_barcode'],keep= 'first')
    germline_dedup.to_csv("unique/"+basename + "_Germline.tab", header = True, sep = '\t', index=False)

    germline_addup_merge.loc[:,'Position'] = np.where(germline_addup_merge['Bait_strand'] == "+", germline_addup_merge.Bait_start, germline_addup_merge.Bait_end)
    germline_addup_dedup = germline_addup_merge.drop_duplicates(['Rname','Strand','Position','Barcode','R2_barcode'],keep='first')
    germline_addup_dedup.to_csv("unique/"+basename + "_Germline_addup_indel.tab", header = True, sep = '\t', index=False)

    substition_dedup = substitution_merge.drop_duplicates(['MDstring','Barcode','R2_barcode'],keep='first')
    substition_dedup.to_csv("unique/"+basename + "_Substitution.tab", header = True, sep = '\t', index=False)

    cmd = "repeats_dedup_e.dist.py {}_SV.tab {} -f 'Rname,Strand,Bait_junction,Junction,R2_barcode'".format(basename, bl)
    print(cmd)
    os.system(cmd)

    cmd = "repeats_dedup_e.dist.2.py {}_Insertion.tab {} -f 'Strand,Bait_junction,Junction,Insertion'".format(basename, bl)
    print(cmd)
    os.system(cmd)

    cmd = "repeats_dedup_e.dist.2.py {}_Deletion.tab {} -f 'Strand,Bait_junction,Junction'".format(basename, bl)
    print(cmd)
    os.system(cmd)
    
    
def main():
    
    start_time = time()
    
    args = docopt(__doc__,version='dedup 1.0')
    
    kwargs = {'basename':args['<basename>'], 'bl':args['<barcode_length>']}
    print('[PEM-Q] basename: ' + str(kwargs['basename']))
    print('[PEM-Q] RMB length: ' + str(kwargs['bl']))
    
    dedup(**kwargs)

    print("\nDedup.py done in {}s".format(round(time()-start_time, 3)))
    
if __name__ == '__main__':
    main()

    
