#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#Author: Xinyi Liu
#Date: 2024

"""define_transloc

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
    define_transloc   <basename> <genome> <cutsite> <primer_strand> <MQ_threshold>

Options:
-h --help               Show this screen.
-v --version            Show version.

This script find bait break site and translocation site from bwa align file.
And also filter reads for indel define.

Input file: reads alignment bam file / Output file: a informative tab file

Author: Mengzhu LIU
Last Update:2019.5.8

"""

    ######################################
    ## i)software dependencies:
    ##    None
    ## ii)data needed:
    ##    1.adapter alignment bam file
    ## iii)python package:
    ##    1.os
    ##    2.pysam
    ##    3.time
    ##    4.docopt
    ##    5.pandas
    ##    6.re
    #######################################
    
from heapq import merge
import os
import pysam
import re
from time import time
from docopt import docopt
import pandas as pd
from Bio import SearchIO
from Bio.Seq import Seq
from Bio import SeqIO
import gzip

class Define_transloc(object):
    
    def __init__(self, basename=None, genome=None, cutsite=None, primer_strand=None, MQ_threshold=None):
        
        self.basename = basename
        self.genome = genome
        bam_sort = "barcode/" + basename + "_sti.filter.sort.bam"
        nosti_r1_bam_sort = "barcode/" + basename + "_nosti.filter.sort.bam"

        nosti_R1_fq = 'flash_out/'+self.basename+'.notCombined_1.fastq.gz'
        nosti_R2_fq = 'flash_out/'+self.basename+'.notCombined_2.fastq.gz'
        
        nosti_R2_chimeric_sam = "transloc/" + self.basename + "_nosti_R2_chimeric.sam"
        nosti_R2_chimeric_bam = "transloc/" + self.basename + "_nosti_R2_chimeric.bam"
        nosti_R2_chimeric_bam_sort = "transloc/" + self.basename + "_nosti_R2_chimeric.sort.bam"

        self.bam_sort = bam_sort
        self.nosti_r1_bam_sort = nosti_r1_bam_sort
        self.nosti_R1_fq = nosti_R1_fq
        self.nosti_R2_fq = nosti_R2_fq
        self.nosti_R2_chimeric_sam = nosti_R2_chimeric_sam
        self.nosti_R2_chimeric_bam = nosti_R2_chimeric_bam
        self.nosti_R2_chimeric_bam_sort = nosti_R2_chimeric_bam_sort

        self.cutsite = int(cutsite)
        self.primer_strand = primer_strand
        self.MQ_threshold = int(MQ_threshold)
    
    def load(self, f_file):
        
        if not os.path.exists(f_file):
            raise ValueError('[PEM-Q] The {} file does not exist.'.format(f_file))   
        
    def reverse_cigar_value(self,cigar):

        letter = re.findall('\D', cigar)
        number = re.findall('\d+', cigar)

        reverse_letter = letter[::-1]
        reverse_number = number[::-1]


        merge_cigar = [None]*(len(reverse_number)+len(reverse_letter))
        merge_cigar[::2] = reverse_number
        merge_cigar[1::2] = reverse_letter

        reverse_cigar = ''.join(merge_cigar)

        return(reverse_cigar)

    def cigar_map_len(self,map_type,cigar):

        #map_type contain: seq/referece
        #reference map len = total - insertion
        #sequence map len = total - deletion
        
        letter = re.findall('\D', cigar)
        number = re.findall('\d+', cigar)
        
        number = list(map(int, number))#convert string to int
        
        #here only want to know mapping length
        #so don't consider strand info and reverse cigar 
        #map_len = (start M + .. + end M) - D
        
        match_list = []
        dele_list = []
        inser_list = []
        for i in range(0,len(letter)):
            if letter[i] == 'M':
                match_list.append(i)
            if letter[i] == 'D':
                dele_list.append(number[i])
            if letter[i] == 'I':
                inser_list.append(number[i])
                
        index_min = min(match_list)
        index_max = max(match_list)

        map_sum = sum(number[index_min:index_max+1])
        del_len = sum(dele_list)
        inser_len = sum(inser_list)
        
        if map_type == 'reference':
            map_len = map_sum - inser_len
        if map_type == 'sequence':
            map_len = map_sum - del_len
        
        return(map_len)

    def cigar_map_seq_start(self,cigar):

        letter = re.findall('\D', cigar)
        number = re.findall('\d+', cigar)
        
        number = list(map(int, number))#convert string to int

        for i in range(0,len(letter)):
            if letter[i] == 'M':
                map_seq_start = sum(number[0:i]) + 1
                break
        return(map_seq_start)

    def cigar_map_seq_end(self,cigar):

        map_seq_start = self.cigar_map_seq_start(cigar)
        map_seq_len = self.cigar_map_len('sequence',cigar)
        map_seq_end = map_seq_start + map_seq_len - 1

        return(map_seq_end)

    def transloc_microhomo(self, rep_sequence, rep_strand, rep_cigarstring, consistant_sup_cigarstring):

        # be careful when use this function because
        # rep_cigarstring have consistant sequential order with sup_cigarstring when have the same trand
        # and vice-versa

        rep_map_seq_start = self.cigar_map_seq_start(rep_cigarstring)
        rep_map_seq_end = self.cigar_map_seq_end(rep_cigarstring)

        sup_map_seq_start = self.cigar_map_seq_start(consistant_sup_cigarstring)
        sup_map_seq_end = self.cigar_map_seq_end(consistant_sup_cigarstring)

        microhomo = ''
        if rep_strand == '-':
            start = rep_map_seq_start
            end = sup_map_seq_end
        else:
            start = sup_map_seq_start
            end = rep_map_seq_end

        if end >= start:
            microhomo = rep_sequence[start-1 : end]
        
        return(microhomo)


    def transloc_find_insertion(self, rep_sequence, rep_strand, rep_cigarstring, consistant_sup_cigarstring):

        rep_map_seq_start = self.cigar_map_seq_start(rep_cigarstring)
        rep_map_seq_end = self.cigar_map_seq_end(rep_cigarstring)

        sup_map_seq_start = self.cigar_map_seq_start(consistant_sup_cigarstring)
        sup_map_seq_end = self.cigar_map_seq_end(consistant_sup_cigarstring)

        insertion = ''
        if rep_strand == '-':
            start = rep_map_seq_start
            end = sup_map_seq_end
        else:
            start = sup_map_seq_start
            end = rep_map_seq_end

        if (start - end) > 1:
            insertion = rep_sequence[end : start-1]
        
        return(insertion)

    def generate_transloc_tab(self): # from stitched reads

        os.system("mkdir transloc")
        os.system("mkdir indel")
        
        bam_file = pysam.AlignmentFile(self.bam_sort, 'rb')
        transloc_tab = open("transloc/" + self.basename + "_mut_sti.tab", "w")
        discard_tab = open("transloc/" + self.basename + "_discard.tab", "w")
        indel_bam = pysam.AlignmentFile("indel/" + self.basename + "_indel.bam", "wb", template=bam_file)
        multi_aln = pd.read_csv("bwa_align/"+self.basename+ "_multi_aln_list.txt", sep='\t', header=None)
        multi_aln.columns = ['Qname','XA_tag']

        transloc_tab.write('Qname'+"\t"+\
                            'Bait_rname'+"\t"+\
                            'Bait_strand'+"\t"+\
                            'Bait_start'+"\t"+\
                            'Bait_end'+"\t"+\
                            'Prey_rname'+"\t"+\
                            'Prey_strand'+"\t"+\
                            'Prey_start'+"\t"+\
                            'Prey_end'+"\t"+\
                            'Rname'+"\t"+\
                            'Strand'+"\t"+\
                            'Junction'+"\t"+\
                            'Sequence'+"\t"+\
                            'B_Qstart'+"\t"+\
                            'B_Qend'+"\t"+\
                            'Qstart'+"\t"+\
                            'Qend'+"\t"+\
                            'Qlen'+"\t"+\
                            'Insertion'+"\t"+\
                            'Microhomolog'+"\t"+\
                            'Prey_MQ'+"\t"+\
                            "Multi_aln"+"\t"+\
                            "annotation"+"\n")

        for read in bam_file:
            multi_aln_check = False

            # because in align_make, reads are already filtered by primer
            # here we only need to classify translocation reads/indel reads/germline reads
            
            # define strand
            
            reference_start = read.reference_start + 1
            reference_end = read.reference_end
            uncut = 0

            if read.is_reverse:
                read_strand = '-'
                if reference_start < self.cutsite:
                    uncut = self.cutsite - reference_start
            else:
                read_strand = '+'
                if read.reference_end > self.cutsite:
                    uncut = read.reference_end - self.cutsite
            
            primer_length = 20   


            #reads contain 'SA' tag have chimeric alignment are potential translocation reads
            
            transloc_check = abs(reference_start - read.reference_end)
            

            condition1 = any('SA' == tg[0] for tg in read.get_tags())
            condition2 = transloc_check > primer_length + 10

            if condition2:

                if condition1:
                    # filter 'uncut' reads
                    if uncut >= 10:
                        continue

                    sa_tag = read.get_tag('SA')
                    # multiple transloclation filter
                    sa_tag_list = sa_tag.split(';')
                
                    sa_tag = ';'.join(sa_tag_list)

                    sa_number = len(sa_tag.split(';'))

                    if sa_number > 2:
                        # c = c + 1
                        continue
                    
                    sa_tag_list = sa_tag.split(',')
                    read_query_name = read.query_name
                    rep_rname = read.reference_name
                    rep_reference_start = reference_start
                    rep_reference_end = read.reference_end
                    rep_strand = read_strand
                    rep_cigarstring = read.cigarstring


                    sup_cigarstring = sa_tag_list[3]
                    sup_map_ref_len = self.cigar_map_len('reference',sup_cigarstring)
            
                    sup_rname = sa_tag_list[0]
                    sup_reference_start = sa_tag_list[1]
                    sup_strand = sa_tag_list[2]                    
                    sup_mapqual = sa_tag_list[4]
                    sup_nm = sa_tag_list[5]
                    sup_reference_end = int(sup_reference_start) + sup_map_ref_len - 1
                    
                    #define transloc junction
                    if sup_strand == '+':
                        transloc_junction = sup_reference_start
                    else:
                        transloc_junction = sup_reference_end
                    # prey mapqual check
                    if int(sup_mapqual) < self.MQ_threshold:
                        continue

                    # multi_aln check
                    XA_tag = None
                    if read_query_name in multi_aln['Qname'].values:
                        XA_tag = multi_aln[multi_aln['Qname'] == read_query_name]['XA_tag'].values[0]
                        XA_list = XA_tag.split(";")
                        XA_list = [item for item in XA_list if item != '']
                        
                        for i in range(len(XA_list)):
                            XA = XA_list[i]
                            XA_detail = XA.split(",")

                            XA_chr = XA_detail[0]
                            if "_" in XA_chr:
                                continue

                            XA_NM = XA_detail[3]
                            if XA_NM <= sup_nm:
                                multi_aln_check = True
                                break
                    
                    if multi_aln_check:
                        continue
                    
                    #define microhomo sequence

                    # rep_cigarstring have consistant sequential order with sup_cigarstring when have the same trand
                    # and vice-versa
            
                    consistant_sup_cigarstring = sup_cigarstring
                    if sup_strand != rep_strand:
                        consistant_sup_cigarstring = self.reverse_cigar_value(sup_cigarstring)#reverse cigarstring

                    rep_sequence = read.query_sequence
                    microhomo = self.transloc_microhomo(rep_sequence,rep_strand,rep_cigarstring,consistant_sup_cigarstring)

                    # define insertion sequence
                    insertion = self.transloc_find_insertion(rep_sequence,rep_strand,rep_cigarstring,consistant_sup_cigarstring)
                    
                    #sequence info
                    if read.is_reverse:
                        read_strand = '-'
                        sequence = Seq(read.query_sequence)
                        sequence = str(sequence.complement())
                        sequence = sequence.upper()
                    else:
                        read_strand = '+'
                        sequence = read.query_sequence
                        sequence = sequence.upper()
                        
                    B_Qstart = self.cigar_map_seq_start(rep_cigarstring)
                    B_Qend = self.cigar_map_seq_end(rep_cigarstring)
                    Qstart = self.cigar_map_seq_start(consistant_sup_cigarstring)
                    Qend = self.cigar_map_seq_end(consistant_sup_cigarstring)
                    Qlen = len(rep_sequence)

                    
                        
                    transloc_tab.write(read_query_name+"\t"+\
                                        rep_rname+"\t"+\
                                        rep_strand+"\t"+\
                                        str(rep_reference_start)+"\t"+\
                                        str(rep_reference_end)+"\t"+\
                                        sup_rname+"\t"+\
                                        sup_strand+"\t"+\
                                        str(sup_reference_start)+"\t"+\
                                        str(sup_reference_end)+"\t"+\
                                        sup_rname+"\t"+\
                                        sup_strand+"\t"+\
                                        str(transloc_junction)+"\t"+\
                                        sequence+"\t"+\
                                        str(B_Qstart)+"\t"+\
                                        str(B_Qend)+"\t"+\
                                        str(Qstart)+"\t"+\
                                        str(Qend)+"\t"+\
                                        str(Qlen)+"\t"+\
                                        insertion+"\t"+\
                                        microhomo+"\t"+\
                                        str(sup_mapqual)+"\t"+\
                                        str(XA_tag)+"\t"+\
                                        "stitched"+"\n"
                                    )
                    #reads do not contain 'SA' tag are potential indel/germline reads
                else:
                    indel_bam.write(read)
            else:
                discard_tab.write(read.query_name+"\n")
        
        
        
        bam_file.close()
        transloc_tab.close()
        indel_bam.close()
        discard_tab.close()

    def find_chimeric_nosti_R1(self): # from non-stitched R1
        bam_file = pysam.AlignmentFile(self.nosti_r1_bam_sort, 'rb')
        R1_transloc_tab = open("transloc/" + self.basename + "_mut_r1.tab", "a")
        discard_tab = open("transloc/" + self.basename + "_discard.tab", "a")
        nosti_indel_bam = pysam.AlignmentFile("indel/" + self.basename + "_nosti_indel.bam", "wb", template=bam_file)
        chimeric_nosti_r1_name = open("transloc/" + self.basename + "_chimeric_nosti_r1_name.txt", "w")
        multi_aln = pd.read_csv("bwa_align/"+self.basename+ "_multi_aln_list.txt", sep='\t', header=None)
        multi_aln.columns = ['Qname','XA_tag']

        R1_transloc_tab.write('Qname'+"\t"+\
                            'Bait_rname'+"\t"+\
                            'Bait_strand'+"\t"+\
                            'Bait_start'+"\t"+\
                            'Bait_end'+"\t"+\
                            'Prey_rname'+"\t"+\
                            'Prey_strand'+"\t"+\
                            'Prey_start'+"\t"+\
                            'Prey_end'+"\t"+\
                            'Rname'+"\t"+\
                            'Strand'+"\t"+\
                            'Junction'+"\t"+\
                            'Sequence'+"\t"+\
                            'B_Qstart'+"\t"+\
                            'B_Qend'+"\t"+\
                            'Qstart'+"\t"+\
                            'Qend'+"\t"+\
                            'Qlen'+"\t"+\
                            'Insertion'+"\t"+\
                            'Microhomolog'+"\t"+\
                            'Prey_MQ'+"\t"+\
                            "Multi_aln"+"\t"+\
                            "annotation"+"\n")

        c = 0
        uncut_reads = 0
        sa_empty = 0        
        mapq_fail = 0
        multi_aln_reads = 0

        for read in bam_file:
            multi_aln_check = False
            
            reference_start = read.reference_start + 1
            reference_end = read.reference_end

            uncut = 0

            if read.is_reverse:
                read_strand = '-'
                if reference_start < self.cutsite:
                    uncut = self.cutsite - reference_start
            else:
                read_strand = '+'
                if read.reference_end > self.cutsite:
                    uncut = read.reference_end - self.cutsite
            
            primer_length = 20   
            #reads contain 'SA' tag have chimeric alignment are potential translocation reads
            
            transloc_check = abs(reference_start - read.reference_end)

            condition1 = any('SA' == tg[0] for tg in read.get_tags())
            condition2 = transloc_check > primer_length + 10
            if condition2:
                if condition1:
                    # filter 'uncut' reads
                    if uncut >= 10:
                        uncut_reads = uncut_reads + 1
                        continue


                    sa_tag = read.get_tag('SA')
                    # multiple transloclation filter
                    sa_tag_list = sa_tag.split(';')
                    sa_tag = ';'.join(sa_tag_list)
                    sa_number = len(sa_tag.split(';'))

                    if sa_number > 2:
                        print(sa_tag)
                        c = c + 1
                        continue
                    
                    sa_tag_list = sa_tag.split(',')
                    read_query_name = read.query_name
                    rep_rname = read.reference_name
                    rep_reference_start = reference_start
                    rep_reference_end = read.reference_end
                    rep_strand = read_strand
                    rep_cigarstring = read.cigarstring

                    sup_cigarstring = sa_tag_list[3]
                    sup_map_ref_len = self.cigar_map_len('reference',sup_cigarstring)
            
                    sup_rname = sa_tag_list[0]
                    sup_reference_start = sa_tag_list[1]
                    sup_strand = sa_tag_list[2]
                    sup_reference_end = int(sup_reference_start) + sup_map_ref_len - 1
                    sup_mapqual = sa_tag_list[4]
                    sup_nm = sa_tag_list[5]
                    
                    #define transloc junction
                    if sup_strand == '+':
                        transloc_junction = sup_reference_start
                    else:
                        transloc_junction = sup_reference_end

                    if int(sup_mapqual) < self.MQ_threshold:
                        continue

                    #define microhomo sequence

                    # rep_cigarstring have consistant sequential order with sup_cigarstring when have the same trand
                    # and vice-versa
            
                    consistant_sup_cigarstring = sup_cigarstring
                    if sup_strand != rep_strand:
                        consistant_sup_cigarstring = self.reverse_cigar_value(sup_cigarstring)#reverse cigarstring

                    rep_sequence = read.query_sequence
                    microhomo = self.transloc_microhomo(rep_sequence,rep_strand,rep_cigarstring,consistant_sup_cigarstring)

                    # define insertion sequence
                    insertion = self.transloc_find_insertion(rep_sequence,rep_strand,rep_cigarstring,consistant_sup_cigarstring)
                    
                    #sequence info
                    if read.is_reverse:
                        read_strand = '-'
                        sequence = Seq(read.query_sequence)
                        sequence = str(sequence.complement())
                    else:
                        read_strand = '+'
                        sequence = read.query_sequence
                        
                    B_Qstart = self.cigar_map_seq_start(rep_cigarstring)
                    B_Qend = self.cigar_map_seq_end(rep_cigarstring)
                    Qstart = self.cigar_map_seq_start(consistant_sup_cigarstring)
                    Qend = self.cigar_map_seq_end(consistant_sup_cigarstring)
                    Qlen = len(rep_sequence)
                    sup_mapqual = sa_tag_list[4]

                    # multi_aln check
                    XA_tag = None
                    if read_query_name in multi_aln['Qname'].values:
                        XA_tag = multi_aln[multi_aln['Qname'] == read_query_name]['XA_tag'].values[0]
                        XA_list = XA_tag.split(";") 
                        XA_list = [item for item in XA_list if item != '']
                      
                        for i in range(len(XA_list)):
                            XA = XA_list[i]
                            XA_detail = XA.split(",")
                            XA_chr = XA_detail[0]
                            XA_loc = int(XA_detail[1])
                            XA_strand = '-' if XA_loc < 0 else '+'
                            XA_start = abs(XA_loc)
                            XA_cigar = XA_detail[2]
                            XA_NM = XA_detail[3]
                            consistant_XA_cigar = XA_cigar
                            if XA_strand != rep_strand:
                                consistant_XA_cigar = self.reverse_cigar_value(XA_cigar) # reverse cigarstring


                            if XA_NM <= sup_nm:
                                XA_len = self.cigar_map_len('reference',XA_cigar)
                                XA_end = int(XA_start) + XA_len - 1
                                junction = XA_start if XA_strand == '+' else XA_end
                                Qstart1 = self.cigar_map_seq_start(consistant_XA_cigar)
                                Qend1 = self.cigar_map_seq_end(consistant_XA_cigar)
                                insertion1 = self.transloc_find_insertion(rep_sequence,rep_strand,rep_cigarstring,consistant_XA_cigar)
                                microhomo1 = self.transloc_microhomo(rep_sequence,rep_strand,rep_cigarstring,consistant_XA_cigar)
                                sup_pos = sup_strand + sup_reference_start
                                XA_tag1 = ",".join([str(sup_rname), str(sup_pos), str(consistant_sup_cigarstring), str(sup_nm)])

                                R1_transloc_tab.write(read_query_name+"\t"+\
                                        rep_rname+"\t"+\
                                        rep_strand+"\t"+\
                                        str(rep_reference_start)+"\t"+\
                                        str(rep_reference_end)+"\t"+\
                                        XA_chr+"\t"+\
                                        XA_strand+"\t"+\
                                        str(XA_start)+"\t"+\
                                        str(XA_end)+"\t"+\
                                        XA_chr+"\t"+\
                                        XA_strand+"\t"+\
                                        str(junction)+"\t"+\
                                        sequence+"\t"+\
                                        str(B_Qstart)+"\t"+\
                                        str(B_Qend)+"\t"+\
                                        str(Qstart1)+"\t"+\
                                        str(Qend1)+"\t"+\
                                        str(Qlen)+"\t"+\
                                        insertion1+"\t"+\
                                        microhomo1+"\t"+\
                                        str(sup_mapqual)+"\t"+\
                                        str(XA_tag1)+"\t"+\
                                        "R1"+"\n")
                        
                    R1_transloc_tab.write(read_query_name+"\t"+\
                                        rep_rname+"\t"+\
                                        rep_strand+"\t"+\
                                        str(rep_reference_start)+"\t"+\
                                        str(rep_reference_end)+"\t"+\
                                        sup_rname+"\t"+\
                                        sup_strand+"\t"+\
                                        str(sup_reference_start)+"\t"+\
                                        str(sup_reference_end)+"\t"+\
                                        sup_rname+"\t"+\
                                        sup_strand+"\t"+\
                                        str(transloc_junction)+"\t"+\
                                        sequence+"\t"+\
                                        str(B_Qstart)+"\t"+\
                                        str(B_Qend)+"\t"+\
                                        str(Qstart)+"\t"+\
                                        str(Qend)+"\t"+\
                                        str(Qlen)+"\t"+\
                                        insertion+"\t"+\
                                        microhomo+"\t"+\
                                        str(sup_mapqual)+"\t"+\
                                        str(XA_tag)+"\t"+\
                                        "R1"+"\n")
                    
                    chimeric_nosti_r1_name.write(read_query_name+"\n")
                else:
                    nosti_indel_bam.write(read)
            else:
                discard_tab.write(read.query_name+"\n")
        bam_file.close()
        R1_transloc_tab.close()
        nosti_indel_bam.close()
        discard_tab.close()
        chimeric_nosti_r1_name.close()

        print("multi_aln_reads",multi_aln_reads)
        print("uncut_reads",uncut_reads)
        print("mapq_fail",mapq_fail)
        print("multi sa",c)

        # merge sti and nosti indel bam
        cmd = "samtools merge -f indel/" + self.basename + "_indel_merge.bam indel/" + self.basename + "_indel.bam indel/" + self.basename + "_nosti_indel.bam"
        os.system(cmd)
        cmd = "mv indel/" + self.basename + "_indel_merge.bam indel/" + self.basename + "_indel.bam"
        os.system(cmd)
        pysam.sort("-o", "indel/" + self.basename + "_indel.sort.bam", "indel/" + self.basename + "_indel.bam")

    def align_chimeric_nosti_R2(self):
        chimeric_nosti_r1_name = "transloc/" + self.basename + "_chimeric_nosti_r1_name.txt"
        nosti_R2_chimeric = "flash_out/" + self.basename + "_nosti_R2_chimeric.fastq"

        with open(chimeric_nosti_r1_name, 'r') as read_name_f:
            read_names = set(line.strip() for line in read_name_f)

        with gzip.open(self.nosti_R2_fq, 'rt', encoding='utf-8') as fastq_f, open(nosti_R2_chimeric, 'w') as output_f:
            for record in SeqIO.parse(fastq_f, 'fastq'):
                read_name = record.id.split("/")[0]
                if read_name in read_names:
                    SeqIO.write(record, output_f, 'fastq')

        bwa_index = os.getenv('BWA_INDEX')
        bwa_index_path = "{}/{}/{}".format(bwa_index, self.genome, self.genome)

        cmd = "bwa mem -Y -t 8 -O 4 -E 1 {} {} > {}".format(bwa_index_path, nosti_R2_chimeric, self.nosti_R2_chimeric_sam)
        os.system(cmd)

        cmd = "samtools view -S -b -h {} > {} \
                && samtools sort {} > {} \
                && samtools index {}".format(self.nosti_R2_chimeric_sam,
                                             self.nosti_R2_chimeric_bam,
                                             self.nosti_R2_chimeric_bam,
                                             self.nosti_R2_chimeric_bam_sort,
                                             self.nosti_R2_chimeric_bam_sort)
        os.system(cmd)

        cmd = "rm {} {}".format(self.nosti_R2_chimeric_sam, self.nosti_R2_chimeric_bam)
        os.system(cmd)
    
    def generate_R2_prey_info(self):
        nosti_R2_bam = pysam.AlignmentFile(self.nosti_R2_chimeric_bam_sort, 'rb')
        R2_transloc_tab = open("transloc/" + self.basename + "_mut_r2.tab", "w")

        R2_transloc_tab.write('Qname'+"\t"+\
                    'R2_rname'+"\t"+\
                    'R2_strand'+"\t"+\
                    'R2_start'+"\t"+\
                    'R2_end'+"\t"+\
                    'R2_MQ'+"\n")

        for read in nosti_R2_bam:
            if read.is_supplementary:
                continue
            # R2 info
            r2_qname = str(read.query_name)
            r2_rname = str(read.reference_name)
            if self.primer_strand == '+':
                r2_reference_start = str(read.reference_start + 1)
                r2_reference_end = str(read.reference_end)
            else:
                r2_reference_start = str(read.reference_end)
                r2_reference_end = str(read.reference_start + 1)
            r2_strand = '-' if read.is_reverse else '+'
            r2_mapqual = str(read.mapping_quality)

            # R2 mapqual control
            if int(r2_mapqual) < self.MQ_threshold:
                continue

            R2_transloc_tab.write(r2_qname+"\t"+\
                                  r2_rname+"\t"+\
                                  r2_strand+"\t"+\
                                  str(r2_reference_start)+"\t"+\
                                  str(r2_reference_end)+"\t"+\
                                  str(r2_mapqual)+"\n")
        
        nosti_R2_bam.close()
        R2_transloc_tab.close()

    def merge_r1_r2_transloc(self):
        dtypes_1 = {
            'Qname': str,
            'Bait_rname': str,
            'Bait_strand': str,
            'Bait_start': str,
            'Bait_end': str,
            'Prey_rname': str,
            'Prey_strand': str,
            'Prey_start': str,
            'Prey_end': str,
            'Rname': str,
            'Strand': str,
            'Junction': str,
            'Sequence': str,
            'B_Qstart': str,
            'B_Qend': str,
            'Qstart': str,
            'Qend': str,
            'Qlen': str,
            'Insertion': str,
            'Microhomolog': str,
            'Prey_MQ': str,
            'Multi_aln': str,
            "annotation": str
        }
        dtypes_2 = {
            'Qname': str,
            'R2_rname': str,
            'R2_strand': str,
            'R2_start': str,
            'R2_end': str,
            'R2_MQ': str
        }
        R1_transloc_tab = pd.read_csv("transloc/" + self.basename + "_mut_r1.tab", sep='\t', dtype=dtypes_1, header=0)
        R2_transloc_tab = pd.read_csv("transloc/" + self.basename + "_mut_r2.tab", sep='\t', dtype=dtypes_2, header=0)
        merge_transloc_tab = pd.merge(R1_transloc_tab, R2_transloc_tab, on='Qname', how='left')

        # filter R2 info and pairing
        ## filter R2 not aligned
        mask = (merge_transloc_tab['R2_rname'] == "None") 
        use_r1_1 = merge_transloc_tab[mask]
        merge_transloc_tab = merge_transloc_tab[~mask]

        ## filter R2 mismatch with R1 prey
        mask = (merge_transloc_tab['R2_strand'] == merge_transloc_tab['Prey_strand']) | (merge_transloc_tab['R2_rname'] != merge_transloc_tab['Prey_rname'])
        use_r1_2 = merge_transloc_tab[mask]
        if use_r1_2.shape[0] > 0:
            merge_transloc_tab = merge_transloc_tab[~mask]

        # filter R2 covered by R1; gap too large;
        merge_transloc_tab['Gap'] = None

        c1 = merge_transloc_tab['Prey_strand'] == merge_transloc_tab['Bait_strand']
        same_strand = merge_transloc_tab[c1]
        diff_strand = merge_transloc_tab[~c1]

        # same strand
        use_r1_41 = pd.DataFrame()
        use_r1_31 = pd.DataFrame()
        if self.primer_strand == '+': # R1 bait and prey both +

            ## filter R2 covered by R1 # fragment length < 150
            mask = (same_strand['Prey_end'] >= same_strand['R2_end'])
            use_r1_41 = same_strand[mask]
            same_strand = same_strand[~mask]

            ## filter gap            
            same_strand['Gap'] = pd.to_numeric(same_strand['R2_start']) - pd.to_numeric(same_strand['Prey_end'])
            mask = (same_strand['Gap'] > 1000) | (same_strand['Gap'] < -150)
            use_r1_31 = same_strand[mask]
            if use_r1_31.shape[0] > 0:
                same_strand = same_strand[~mask]

            # pairing and renew prey info
            mask = (same_strand['Prey_end'] < same_strand['R2_end'])
            same_strand.loc[mask, 'Prey_end'] = same_strand.loc[mask, 'R2_end']
            same_strand.loc[mask, 'Prey_MQ'] = same_strand.loc[mask, 'R2_MQ']            

        if self.primer_strand == '-': # R1 bait and prey both -

            ## filter R2 covered by R1
            mask = (same_strand['Prey_start'] <= same_strand['R2_start'])
            use_r1_41 = same_strand[mask]
            same_strand = same_strand[~mask]

            ## filter gap            
            same_strand['Gap'] = pd.to_numeric(same_strand['Prey_start']) - pd.to_numeric(same_strand['R2_end'])
            mask = (same_strand['Gap'] > 1000) | (same_strand['Gap'] < -150)
            use_r1_31 = same_strand[mask]
            if use_r1_31.shape[0] > 0:
                same_strand = same_strand[~mask]

            # pairing and renew prey info
            mask = (same_strand['Prey_start'] > same_strand['R2_start'])
            same_strand.loc[mask, 'Prey_start'] = same_strand.loc[mask, 'R2_start']
            same_strand.loc[mask, 'Prey_MQ'] = same_strand.loc[mask, 'R2_MQ']

        # diff strand # inversion
        use_r1_42 = pd.DataFrame()
        use_r1_32 = pd.DataFrame()
        if self.primer_strand == '+': # R1 bait +, prey -

            ## filter R2 covered by R1
            mask = (diff_strand['Prey_end'] <= diff_strand['R2_end'])
            use_r1_42 = diff_strand[mask]
            diff_strand = diff_strand[~mask]

            ## filter gap            
            diff_strand['Gap'] = pd.to_numeric(diff_strand['Prey_end']) - pd.to_numeric(diff_strand['R2_start'])
            mask = (diff_strand['Gap'] > 1000) | (diff_strand['Gap'] < -150)
            use_r1_32 = diff_strand[mask]
            if use_r1_32.shape[0] > 0:
                diff_strand = diff_strand[~mask]

            # pairing and renew prey info
            mask = (diff_strand['Prey_end'] > diff_strand['R2_end'])
            diff_strand.loc[mask, 'Prey_end'] = diff_strand.loc[mask, 'R2_end']
            diff_strand.loc[mask, 'Prey_MQ'] = diff_strand.loc[mask, 'R2_MQ']

        if self.primer_strand == '-': # R1 bait -, prey +

            ## filter R2 covered by R1
            mask = (diff_strand['Prey_start'] >= diff_strand['R2_start'])
            use_r1_42 = diff_strand[mask]
            diff_strand = diff_strand[~mask]

            ## filter gap            
            diff_strand['Gap'] = pd.to_numeric(diff_strand['R2_end']) - pd.to_numeric(diff_strand['Prey_start'])
            mask = (diff_strand['Gap'] > 1000) | (diff_strand['Gap'] < -150)
            use_r1_32 = diff_strand[mask]
            if use_r1_32.shape[0] > 0:
                diff_strand = diff_strand[~mask]

            # pairing and renew prey info
            mask = (diff_strand['Prey_start'] < diff_strand['R2_start'])
            diff_strand.loc[mask, 'Prey_start'] = diff_strand.loc[mask, 'R2_start']
            diff_strand.loc[mask, 'Prey_MQ'] = diff_strand.loc[mask, 'R2_MQ']

        use_r1_4 = pd.concat([use_r1_41, use_r1_42], ignore_index=True)
        use_r1_3 = pd.concat([use_r1_31, use_r1_32], ignore_index=True)

        merge_transloc_tab = pd.concat([same_strand, diff_strand], ignore_index=True)
        
        
        # multi_aln check for use_R1_only
        use_r1 = pd.concat([use_r1_1, use_r1_4], ignore_index=True) # drop not cherished R1
        use_r1 = use_r1.drop_duplicates(['Qname'],keep=False) # keep=False = XA tag filter
        c = use_r1.shape[0]

        # multi_aln check for merge_transloc_tab
        merge_clean = merge_transloc_tab.drop_duplicates(['Qname'],keep = False)
        merge_dup = merge_transloc_tab[merge_transloc_tab.duplicated(['Qname'],keep= False)]
        merge_dup = merge_dup.loc[merge_dup.groupby('Qname')['Gap'].idxmin()]
        merge_transloc_tab = pd.concat([merge_clean, merge_dup], ignore_index=True)
        
        # generate output
        merge_transloc_tab = pd.concat([merge_transloc_tab, use_r1], ignore_index=True)
        merge_transloc_tab = merge_transloc_tab[['Qname','Bait_rname','Bait_strand','Bait_start','Bait_end','Prey_rname','Prey_strand','Prey_start',\
                                                 'Prey_end','Rname','Strand','Junction','Sequence','B_Qstart','B_Qend','Qstart','Qend','Qlen','Insertion','Microhomolog','Prey_MQ','Multi_aln']]

        transloc_tab = pd.read_csv("transloc/" + self.basename + "_mut_sti.tab", sep='\t')
        transloc_tab = pd.concat([transloc_tab, merge_transloc_tab], ignore_index=True)
        transloc_tab.to_csv("transloc/" + self.basename + "_mut.tab", sep='\t', index=False)


        

def main():
    
    start_time = time()
    
    args = docopt(__doc__,version='denfine_transloc 1.0')
    
    kwargs = {'basename':args['<basename>'], 'genome':args['<genome>'], 'cutsite':args['<cutsite>'], 'primer_strand':args['<primer_strand>'], 'MQ_threshold':args['<MQ_threshold>']}
    print('[PEM-Q] basename: ' + str(kwargs['basename']))
    print('[PEM-Q] genome: ' + str(kwargs['genome']))
    print('[PEM-Q] cutsite: ' + str(kwargs['cutsite']))
    print('[PEM-Q] primer_strand: ' + str(kwargs['primer_strand']))
    print('[PEM-Q] MQ_threshold: ' + str(kwargs['MQ_threshold']))
    
    define_transloc = Define_transloc(**kwargs)
    define_transloc.generate_transloc_tab()
    define_transloc.find_chimeric_nosti_R1()
    define_transloc.align_chimeric_nosti_R2()
    define_transloc.generate_R2_prey_info()
    define_transloc.merge_r1_r2_transloc()
    

    print("\ndefine_transloc.py Done in {}s".format(round(time()-start_time, 3)))
    
if __name__ == '__main__':
    main()
