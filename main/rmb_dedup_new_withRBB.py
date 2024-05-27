#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#Author: Mengzhu
#Date:2019.4.23
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
    rmb_dedup   <basename> <barcode_length> <adapter>

Options:
-h --help               Show this screen.
-v --version            Show version.

Random barcodes and R2 blue barcodes are extracted from the adapter alignment bam file, and reads are deduplicated according to the barcode attached.

Input file: adapter alignment bam file / Output file: unique reads's query name list

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
    #######################################
    
import os
import pysam
from time import time
from docopt import docopt
import pandas as pd

class Dedup(object):
    
    def __init__(self, basename=None, barcode_length=None, adapter=None):
        
        self.basename = basename
        self.barcode_length = int(barcode_length)
        
        # init adapter 
        adapter_f = "adapter/adapter.fa" 
        self.load(adapter_f)  
        adapt_f = pysam.FastaFile(adapter_f)
        adapter = adapt_f.fetch(reference = 'adapter')
        adapt_f.close()
        print("[PEM-Q] adapter_seqeunce: " + adapter)
        self.adapter = adapter
        
        # init bam file
        fastq_check = "flash_out/" + basename + ".extendedFrags.fastq.gz"
        adapter_sam_check = "barcode/" + basename + "_check.adpt.sam"
        adapter_bam_check = "barcode/" + basename + "_check.adpt.bam"
        adapter_bam_check_sort = "barcode/" + basename + "_check.adpt.sort.bam"
        
        nosti_r2_fastq = "flash_out/" + basename + ".notCombined_2.fastq.gz"
        nosti_adapter_sam_check = "barcode/" + basename + "_nosti_check.adpt.sam"
        nosti_adapter_bam_check = "barcode/" + basename + "_nosti_check.adpt.bam"
        nosti_adapter_bam_check_sort = "barcode/" + basename + "_nosti_check.adpt.sort.bam"
        
        adapter_bam = "bwa_align/" + basename + "_sti.adpt.sort.bam" # get adapter bam file # use raw R2 fq for alignment, nosti included
        primer_bam = "primer/" + basename + "_sti.p.sort.bam" # get primer filter stitch bam file
        nosti_r1_primer_bam = "primer/" + basename + "_nosti_r1.p.sort.bam" # get primer filter stitch bam file

        dedup_bam = "barcode/" + basename + "_sti.dedup.bam"
        dedup_bam_sort = "barcode/" + basename + "_sti.dedup.sort.bam"
        dup_bam = "barcode/" + basename + "_sti.dup.bam"
        dup_bam_sort = "barcode/" + basename + "_sti.dup.sort.bam"
        filter_bam = "barcode/" + basename + "_sti.filter.bam"
        filter_bam_sort = "barcode/" + basename + "_sti.filter.sort.bam"
        nosti_dedup_bam = "barcode/" + basename + "_nosti.dedup.bam"
        nosti_dedup_bam_sort = "barcode/" + basename + "_nosti.dedup.sort.bam"
        nosti_dup_bam = "barcode/" + basename + "_nosti.dup.bam"
        nosti_dup_bam_sort = "barcode/" + basename + "_nosti.dup.sort.bam"
        nosti_filter_bam = "barcode/" + basename + "_nosti.filter.bam"
        nosti_filter_bam_sort = "barcode/" + basename + "_nosti.filter.sort.bam"
        
        
        self.fastq_check = fastq_check
        self.adapter_sam_check = adapter_sam_check
        self.adapter_bam_check = adapter_bam_check
        self.adapter_bam_check_sort = adapter_bam_check_sort
        
        self.nosti_r2_fastq = nosti_r2_fastq
        self.nosti_adapter_sam_check = nosti_adapter_sam_check
        self.nosti_adapter_bam_check = nosti_adapter_bam_check
        self.nosti_adapter_bam_check_sort = nosti_adapter_bam_check_sort
        
        self.adapter_bam = adapter_bam
        self.primer_bam = primer_bam
        self.nosti_r1_primer_bam = nosti_r1_primer_bam

        self.dedup_bam = dedup_bam
        self.dedup_bam_sort = dedup_bam_sort
        self.dup_bam = dup_bam
        self.dup_bam_sort = dup_bam_sort
        self.filter_bam = filter_bam
        self.filter_bam_sort = filter_bam_sort
        
        self.nosti_dedup_bam = nosti_dedup_bam
        self.nosti_dedup_bam_sort = nosti_dedup_bam_sort
        self.nosti_dup_bam = nosti_dup_bam
        self.nosti_dup_bam_sort = nosti_dup_bam_sort
        self.nosti_filter_bam = nosti_filter_bam
        self.nosti_filter_bam_sort = nosti_filter_bam_sort


    def load(self, f_file):
        
        if not os.path.exists(f_file):
            raise ValueError('[PEM-Q] The {} file does not exist.'.format(f_file))
    

    def extract_barcode(self):                
            
        # only keep query with barcode length >= barcode_length-2
        os.system("mkdir barcode")
        bam_file = pysam.AlignmentFile(self.adapter_bam, "rb")
        barcode_list = open("barcode/"+self.basename+"_barcode_list.txt", "w")
        
        for read in bam_file.fetch():
            
            if read.query_alignment_length >= (len(self.adapter) - 2):#aligned adapter length allow 2 deletion
                # RMB
                bar_start = read.query_alignment_end + 1
                bar_end = read.query_length
                bar_length = bar_end - bar_start + 1
                # R2 barcode
                r2_barcode_length = 10
                r2_bar_start = read.query_alignment_start-r2_barcode_length
                r2_bar_end = read.query_alignment_start-1
                r2_barcode = str(read.query_sequence[r2_bar_start:r2_bar_end+1])
                if bar_length <= self.barcode_length and bar_length >=(self.barcode_length - 2) and r2_barcode != '': # RMB allow 2 deletion
                    barcode_seq = read.query_sequence[(bar_start-1) : bar_end]
                    barcode_list.write(str(read.query_name) + "\t" + str(barcode_seq) + "\t" + str(r2_barcode) + "\n")
        barcode_list.close()
    
    def barcode_dedup(self):  
        
        data = pd.read_csv("barcode/"+self.basename+"_barcode_list.txt",sep = '\t',names = ["Qname", "Barcode", "R2_barcode"])
        # get barcode length
        length = data['Barcode'].astype('str')
        length = data['Barcode'].str.len()
        data['Length'] = length
        # get barcode frequency
        data['Freq'] = data.groupby('Barcode')['Barcode'].transform('count')
        data = data.sort_values(by=['Freq','Barcode','Length'], ascending=False)
        data.to_csv("barcode/"+self.basename+"_barcode_sort.txt", header = True, sep = '\t', index=False,
           columns = [u'Qname',u'Barcode',u'Freq',u'Length',u'R2_barcode'])
           
        #also generate duplicated barcode
        datb = data[data.duplicated(['Barcode'])]
        datb = datb.reset_index(drop=True)
        #unique barcode
        data = data.drop_duplicates(['Barcode'], keep = 'first' )#rough dedup
        data = data.reset_index(drop=True)#reset index
        
        data.to_csv("barcode/"+self.basename+"_barcode_uniq.txt", header = True, sep = '\t', index=False,
           columns = [u'Qname',u'Barcode',u'Freq',u'Length',u'R2_barcode'])
        datb.to_csv("barcode/"+self.basename+"_barcode_dup.txt", header = True, sep = '\t', index=False,
           columns = [u'Qname',u'Barcode',u'Freq',u'Length',u'R2_barcode'])

    
    def filter_multiple_adapter(self):
        
        # stitch reads
        chek_file = self.adapter_bam_check_sort
        
        if not os.path.exists(chek_file):
            
            if not self.fastq_check:
                raise ValueError('Stitched fastq is needed for adapter check alignment.')
            if not os.path.exists("adapter/adapter.fa"):
                raise ValueError('adapter/adapter.fa is needed for adapter check alignment.')
                    
            #alignment
            print("[PEM-Q] align to check adapter...")
                              
            cmd = "bwa mem -t 8 adapter/adapter -k 5 -L 0 -T 14 {} > {} 2>barcode/bwa_align_adapter.log".format(self.fastq_check, 
                                                              self.adapter_sam_check)
            os.system(cmd)
            print("[PEM-Q] "+cmd)
            cmd = "samtools view -S -b -h {} > {} \
                   && samtools sort {} > {} \
                   && samtools index {}".format(self.adapter_sam_check, self.adapter_bam_check, self.adapter_bam_check, self.adapter_bam_check_sort, self.adapter_bam_check_sort)
            print("[PEM-Q] sort and index bam...")
            os.system(cmd)
        else:
            print("[PEM-Q] adapter check alignment file exist, jump...")  
        
        #keep record of multiple adapters
        multiple_adapt = open("barcode/"+self.basename+"_multiple_adapt.txt", "w")
        clean_adapt = open("barcode/"+self.basename+"_clean_adapt.txt", "w")
        bam_file = pysam.AlignmentFile(self.adapter_bam_check_sort, "rb")
        multiple_adapt_list = []
        clean_adapt_list = []
        for read in bam_file:
            condition1 = any('SA' == tg[0] for tg in read.get_tags())
            if condition1:
                multiple_adapt_list.append(read.query_name)
                multiple_adapt.write(read.query_name+"\n")
            else:
                clean_adapt_list.append(read.query_name)
                clean_adapt.write(read.query_name+"\n")
        multiple_adapt.close()
        clean_adapt.close()
        bam_file.close()
        #remove reads with multiple adapters
        primer_bam = pysam.AlignmentFile(self.primer_bam, 'rb')
        dedup_bam_sort = primer_bam
        filter_bam = pysam.AlignmentFile(self.filter_bam, "wb", template=dedup_bam_sort)
        name_indexed = pysam.IndexedReads(dedup_bam_sort)
        name_indexed.build()
        for name in clean_adapt_list:
                try:
                    name_indexed.find(name)
                except KeyError:
                    pass
                else:
                    iterator = name_indexed.find(name)
                    for x in iterator:
                        filter_bam.write(x)
        dedup_bam_sort.close()
        primer_bam.close()
        filter_bam.close()
        pysam.sort("-o", self.filter_bam_sort, self.filter_bam)
        

        # non-stitch reads
        cmd = "bwa mem -t 8 adapter/adapter -k 5 -L 0 -T 14 {} > {} 2>barcode/bwa_align_adapter.log".format(self.nosti_r2_fastq, 
                                                            self.nosti_adapter_sam_check)
        os.system(cmd)
        print("[PEM-Q] "+cmd)
        cmd = "samtools view -S -b -h {} > {} \
                && samtools sort {} > {} \
                && samtools index {}".format(self.nosti_adapter_sam_check, self.nosti_adapter_bam_check, self.nosti_adapter_bam_check, self.nosti_adapter_bam_check_sort, self.nosti_adapter_bam_check_sort)
        print("[PEM-Q] sort and index bam...")
        os.system(cmd)

        bam_file = pysam.AlignmentFile(self.nosti_adapter_bam_check_sort, "rb")
        multiple_adapt_list = []
        clean_adapt_list = []
        for read in bam_file:
            condition1 = any('SA' == tg[0] for tg in read.get_tags())
            if condition1:
                multiple_adapt_list.append(read.query_name)
            else:
                clean_adapt_list.append(read.query_name)
        bam_file.close()
        #remove reads with multiple adapters
        primer_bam = pysam.AlignmentFile(self.nosti_r1_primer_bam, 'rb')
        dedup_bam_sort = primer_bam
        filter_bam = pysam.AlignmentFile(self.nosti_filter_bam, "wb", template=dedup_bam_sort)
        name_indexed = pysam.IndexedReads(dedup_bam_sort)
        name_indexed.build()
        for name in clean_adapt_list:
                try:
                    name_indexed.find(name)
                except KeyError:
                    pass
                else:
                    iterator = name_indexed.find(name)
                    for x in iterator:
                        filter_bam.write(x)
        dedup_bam_sort.close()
        primer_bam.close()
        filter_bam.close()
        pysam.sort("-o", self.nosti_filter_bam_sort, self.nosti_filter_bam)
        

        return()

def main():
    start_time = time()
    
    args = docopt(__doc__,version='rmb_dedup 1.0')
    
    kwargs = {'basename':args['<basename>'], 'barcode_length':args['<barcode_length>']}
    print('[PEM-Q] basename: ' + str(kwargs['basename']))
    print('[PEM-Q] barcode_length: ' + str(kwargs['barcode_length']))
    
    dedup = Dedup(**kwargs)
    dedup.extract_barcode()
    dedup.barcode_dedup()
    dedup.filter_multiple_adapter()

    print("\nRmb_dedup.py done in {}s".format(round(time()-start_time, 3)))
    
if __name__ == '__main__':
    main()