#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#Author: Mengzhu
#Date:2019.10.18

"""define_substitution

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
    define_substitution    <basename>  <cutsite>  <cutoff> <primer_strand>

Options:
-h --help               Show this screen.
-v --version            Show version.

The definition of substitution is the same as in PEM-Q. This script is identical to the one in PEM-Q.

"""


import re
from time import time
from docopt import docopt
import pandas as pd

def cal_soft_clipping_number(cigar):
    
    cigar_list = re.split('(\d+)',cigar)
    if cigar_list[2] == 'S':
        number = int(cigar_list[1])
    else:
        number = 0
    return(number)
    
def define_substitution(basename=None,cutsite=None,cutoff=None,primer_strand=None):
    
    substitution_file = "unique/" + basename + "_Substitution.tab"
    data = pd.read_csv(substitution_file, sep = '\t', index_col=False, low_memory=False)
    data_cutoff = pd.DataFrame()
    germline_addup_st = pd.DataFrame()
    
    mdstring = data["MDstring"]
    count = data["MDstring"].count()
    n = 0

    ## new, use list instead of dataframe
    data_cutoff_rows = []

    for i in range(0, count):
        md_string = data["MDstring"][i]
        md_list = re.split('(\d+)',md_string)
        md_string_list = md_list[1:(len(md_list)-1)]
        len_list = []
        seq_length = 0
        sub_length = 0
        for j in range(0, len(md_string_list)):
            if md_string_list[j].isdigit():
                len_list.append(int(md_string_list[j]))
                seq_length = sum(len_list)
            else:
                sub_length = sub_length + 1
                seq_length = seq_length + sub_length
                bed_start = (int(data["Position"][i]) + seq_length - 1) - 1
                bed_end = bed_start + 1 # real position

                if bed_start >= (int(cutsite) - int(cutoff)) and bed_end <= (int(cutsite) + int(cutoff) - 1):
                    n = n + 1
                    data_cutoff_rows.append(data.iloc[i])
                    break
    
    data_cutoff = pd.DataFrame(data_cutoff_rows)
    data_cutoff.to_csv("unique/" + basename + "_Substitution_cutoff.tab", header = True, sep = '\t', index=False)
    print("Substitutions in cutsite +- cutoff:",n)
    if n == 0:
        germline_addup_st = data
    else:
        germline_addup_st = data[~data['Qname'].isin(data_cutoff['Qname'])]
    # germline_addup_st = data.append([data_cutoff])
    germline_addup_st = germline_addup_st.drop_duplicates(keep=False)
    germline_addup_st.to_csv("unique/" + basename + "_Germline_addup_st.tab", header = True, sep = '\t', index=False)
    print(germline_addup_st["Qname"].count())
    
    
    sv_file = "unique/" + basename+"_SV.tab"
    data = pd.read_csv(sv_file, sep = '\t', index_col=False, low_memory=False)
    
    # if not data.empty:
    #inversions
    chrom = data["Bait_rname"][0]
    inversions_con = (data["Prey_rname"] == chrom)
    inversions = data[inversions_con]
    inversions.to_csv("unique/" + basename+"_inversion.tab", header = True, sep = '\t', index=False)
    #translocations(inter)
    interTransloc_con = (data["Prey_rname"] != chrom)
    interTransloc = data[interTransloc_con]
    interTransloc.to_csv("unique/" + basename+"_interSV.tab", header = True, sep = '\t', index=False)
    # else:
    #     print("Warning: SV file is empty!")
    
    germ_file = "unique/" + basename+"_Germline.tab"
    germ_addup_st = pd.read_csv("unique/" + basename+"_Germline_addup_st.tab", sep = '\t', index_col=False, low_memory=False)
    germ_addup_indel = pd.read_csv("unique/" + basename+"_Germline_addup_indel.tab", sep = '\t', index_col=False, low_memory=False)
    data = pd.read_csv(germ_file, sep = '\t', index_col=False, low_memory=False)
    data = pd.concat([data, germ_addup_st])
    data = pd.concat([data, germ_addup_indel])
    data.to_csv("unique/" + basename + "_Germline_final.tab", header = True, sep = '\t', index=False)
    
def main():

    start_time = time()

    args = docopt(__doc__,version='define_substitution 1.0')

    kwargs = {'basename':args['<basename>'],'cutsite':args['<cutsite>'],'cutoff':args['<cutoff>'],'primer_strand':args['<primer_strand>']}
    print('[PEM-Q] basename: ' + str(kwargs['basename']))
    print('[PEM-Q] cutsite: ' + str(kwargs['cutsite']))
    print('[PEM-Q] cutoff: ' + str(kwargs['cutoff']))
    print('[PEM-Q] primer_strand: ' + str(kwargs['primer_strand']))
    
    
    define_substitution(**kwargs)

    print("\ndefine_substitution.py Done in {}s".format(round(time()-start_time, 3)))

if __name__ == '__main__':
    main()
