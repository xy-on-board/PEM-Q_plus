#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#Author: Xinyi Liu
#Date: 2024

"""Align make

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
    Align_make  [-r <primerChrom>] [-s <primerStart>] [-e <primerEnd>] [-d <primerStrand>] [-p <sequence>] [-a <adapter>] <genome> <fastq_r1>
    Align_make  [-r <primerChrom>] [-s <primerStart>] [-e <primerEnd>] [-d <primerStrand>] [-p <sequence>] [-a <adapter>] <genome> <fastq_r1> <fastq_r2>

Options:
-h --help               Show this screen.
-v --version            Show version.
-p <sequence>           filter reads without primer sequence
-a <adapter>            do adapter alignment
-r <primerChrom>        primer chromosome
-s <primerChrom>        primer start
-e <primerChrom>        primer end
-d <primerChrom>        primer strand('+':1;'-':-1)

In this step, fastq files are aligned to the reference genome and reads are filtered based on primer sequence.

Input file: fastq file / Output file: sorted bam file

"""

    
    ######################################
    ## i)software dependencies:
    ##    1.FLASH
    ##    2.bwa mem
    ##    2.samtools
    ##    3.blat
    ## ii)data needed:
    ##    1.genome data(fa file)
    ##    2.bwa index file of genome
    ## iii)python package:
    ##    1.os
    ##    2.pysam
    ##    3.biopython
    ##    4.docopt
    ##    5.threading
    #######################################
    
    
import os
import sys
import pysam
import threading
from time import time
from docopt import docopt
from Bio import SearchIO


class Align_make(object):
    
    def __init__(self, primerChrom=None, primer=None, primer_start=None, primer_end=None, primer_strand=None, adapter=None, genome=None, fastq_r1=None, fastq_r2=None):
        
        # print("[PEM-Q] initiating fastq file...")
        
        # super().__init__()
        self.primer = primer
        self.primer_chr = primerChrom
        self.primer_start = primer_start
        self.primer_end = primer_end
        self.primer_strand = primer_strand
        self.adapter = adapter
        self.genome = genome
        self.fastq_r1 = fastq_r1
        self.fastq_r2 = fastq_r2
        self.basename = fastq_r1.rpartition('/')[2].partition('.')[0].rpartition('_')[0]
        self.outputdir = self.basename
        
        if fastq_r1:
            self.load(fastq_r1)
        if fastq_r2:
            self.load(fastq_r2)
        else:
            self.fastq_r2=""
        
        if self.fastq_r1 and self.fastq_r2:
            
            self.pairFlag = 1
            
            self.sam = self.basename+'_sti.sam'
            self.bam = self.basename+'_sti.bam'
            self.bam_sort = self.basename+'_sti.sort.bam'
            
            self.p_bam = self.basename+'_sti.p.bam'
            self.p_bam_sort = self.basename+'_sti.p.sort.bam'
            
            self.no_stitch_r1_sam = self.basename+'_nosti_r1.sam'
            self.no_stitch_r1_bam = self.basename+'_nosti.bam'
            self.no_stitch_r1_bam_sort = self.basename+'_nosti_r1.sort.bam'
            
            self.p_no_stich_r1_bam = self.basename+'_nosti_r1.p.bam'
            self.p_no_stich_r1_bam_sort = self.basename+'_nosti_r1.p.sort.bam'
            
            self.adpt_sam = self.basename+'_sti.adpt.sam'
            self.adpt_bam = self.basename+'_sti.adpt.bam'
            self.adpt_bam_sort = self.basename+'_sti.adpt.sort.bam'
            
        else:
            self.pairFlag = 0
            self.sam = self.basename+'_s.sam'
            self.bam = self.basename+'_s.bam'
            self.bam_sort = self.basename+'_s.sort.bam'
        
    def load(self, fastq):
        
        if not os.path.exists(fastq):
            raise ValueError('[PEM-Q] The {} file does not exist.'.format(fastq))
            # exit()
    
    def stitch(self):
        
        chek_file = 'flash_out/'+self.basename+'.extendedFrags.fastq.gz'
        
        if not os.path.exists(chek_file):

            cmd = " flash {} {} -o {} -d flash_out > flash.log".format(self.fastq_r1,self.fastq_r2,self.outputdir)
            
            if self.pairFlag == 1:
                os.system(cmd)
                print("[PEM-Q] stitching reads using FLASH...")
                
                sti_fq = 'flash_out/'+self.basename+'.extendedFrags.fastq'
                nosti_R1_fq = 'flash_out/'+self.basename+'.notCombined_1.fastq'
                nosti_R2_fq = 'flash_out/'+self.basename+'.notCombined_2.fastq'
                
                cmd = "gzip {}".format(sti_fq)
                os.system(cmd)
                cmd = "gzip {}".format(nosti_R1_fq)
                os.system(cmd)
                cmd = "gzip {}".format(nosti_R2_fq)
                os.system(cmd)

            else:
                print("[PEM-Q] Only one fastq file. Stop stitcing...")
        else:
            print("[PEM-Q] flash merge file exist, jump...")
        return()
               
    def align_genome(self): 
             
        chek_file = 'bwa_align/'+self.bam_sort
        if not os.path.exists(chek_file):

            # define if input is pair end data
            if self.pairFlag == 1:
                self.stitch()
                bwa_fq_file = 'flash_out/'+self.basename+'.extendedFrags.fastq.gz'
                no_stich_fq_r1 = 'flash_out/'+self.basename+'.notCombined_1.fastq.gz'
                no_stich_fq_r2 = 'flash_out/'+self.basename+'.notCombined_2.fastq.gz'
            else:
                bwa_fq_file = self.fastq_r1
            
            
            #bwa align
            cmd = "mkdir bwa_align"
            os.system(cmd)
            
            #get bwa index path
            bwa_index = os.getenv('BWA_INDEX')
            bwa_index_path = "{}/{}/{}".format(bwa_index, self.genome, self.genome)
            print("[PEM-Q] index file used " + bwa_index_path)
            print("[PEM-Q] align to genome...")
            
            #-A 2 -B 8 -O 8 -E 4 -T 50 -k 20 -p -M 
            cmd = "bwa mem -Y -t 8 -O 4 -E 1 {} {} > bwa_align/{} 2>bwa_align/bwa_alignstich.log".format(bwa_index_path, 
                                                bwa_fq_file,
                                                self.sam)
            os.system(cmd)
            print("[PEM-Q] "+cmd)
            
            cmd = "samtools view -S -b -h bwa_align/{} > bwa_align/{} \
                   && samtools sort bwa_align/{} > bwa_align/{} \
                   && samtools index bwa_align/{}".format(self.sam, 
                                                self.bam, 
                                                self.bam, 
                                                self.bam_sort, 
                                                self.bam_sort)
            os.system(cmd)

            #no stitch reads
            if self.pairFlag == 1:
                cmd = "bwa mem -Y -t 8 {} {} > bwa_align/{} 2>bwa_align/bwa_alignnostich.log".format(bwa_index_path,
                                                                                                    no_stich_fq_r1,
                                                                                                    self.no_stitch_r1_sam)
                os.system(cmd)
                print("[PEM-Q] "+cmd)

                cmd = "samtools view -S -b -h bwa_align/{} > bwa_align/{} \
                    && samtools sort bwa_align/{} > bwa_align/{} \
                        && samtools index bwa_align/{}".format(self.no_stitch_r1_sam,
                                                               self.no_stitch_r1_bam,
                                                               self.no_stitch_r1_bam,
                                                               self.no_stitch_r1_bam_sort,
                                                               self.no_stitch_r1_bam_sort)
                os.system(cmd)
        else:
            print("[PEM-Q] genome alignment file exist, jump...")
                   
        return()

    def align_adapter(self):
        
        chek_file = 'bwa_align/'+self.adpt_bam_sort
        if not os.path.exists(chek_file):

            print("[PEM-Q]  your adapter sequence: "+self.adapter)
        
            if not self.fastq_r2:
                raise ValueError('R2 fastq file is needed for adapter alignment.')
                                    
            os.system("mkdir adapter")
            os.system("mkdir bwa_align")
            adapter_f = open("adapter/adapter.fa", "w")
            adapter_f.write(">adapter\n")
            adapter_f.write(self.adapter)
            adapter_f.close()
        
            #build bwa index of adapter
            cmd = "samtools faidx adapter/adapter.fa"
            os.system(cmd)
            cmd = "bwa index -a bwtsw -p adapter/adapter adapter/adapter.fa 1>adapter/build_index.o 2>adapter/build_index.e"
            os.system(cmd)
            
            #alignment
            print("[PEM-Q]  align to adapter...")
                              
            cmd = "bwa mem -t 8 adapter/adapter -k 10 -L 0 -T 10 {} > bwa_align/{} 2>bwa_align/bwa_align_adapter.log".format(self.fastq_r2, 
                                                              self.adpt_sam)
            os.system(cmd)
            print("[PEM-Q] "+cmd)
            cmd = "samtools view -S -b -h bwa_align/{} > bwa_align/{} \
                   && samtools sort bwa_align/{} > bwa_align/{} \
                   && samtools index bwa_align/{}".format(self.adpt_sam, self.adpt_bam, self.adpt_bam, self.adpt_bam_sort, self.adpt_bam_sort)
            print("[PEM-Q]  sort and index bam...")
            os.system(cmd)
        else:
            print("[PEM-Q]  adapter alignment file exist, jump...")
        
        return()
        
    def primer_check(self):

        print("[PEM-Q] filter no_primer reads...")
        print("[PEM-Q] Your primer sequence: "+self.primer)

        os.system("mkdir primer")
        primer_f = open("primer/primer.fa", "w")
        primer_f.write(">primer\n")
        primer_f.write(self.primer)
        primer_f.close()

        primer_chr = self.primer_chr
        primer_end = int(self.primer_end)
        primer_start = int(self.primer_start)
        if self.primer_strand == '+':
            primer_strand = 1
        if self.primer_strand == '-':
            primer_strand = -1

        print("[PEM-Q] primer position: " + primer_chr, primer_strand, primer_start, primer_end)
            
        #keep reads with correct primer location
        
        #stitch reads
        bam_file = pysam.AlignmentFile('bwa_align/'+self.bam_sort, 'rb')
        primer_bam = pysam.AlignmentFile('primer/'+self.p_bam, "wb", template=bam_file)
        bam_list = open("primer/bamlist_stitch.txt", "w")
        multi_aln_list = open("bwa_align/"+self.basename+ "_multi_aln_list.txt", "w")
        count = 0
        
        for read in bam_file:
            if primer_strand == 1:
                read_check = read.reference_start+1
                primer_check = primer_start
            else:
                read_check = read.reference_end
                primer_check = primer_end
            
            if read_check is None:
                continue
            
            if read.is_reverse:
                read_strand = -1
            else:
                read_strand = 1
            
            condition1 = bam_file.getrname(read.reference_id) == primer_chr
            condition2 = read_strand == primer_strand
            condition3 = read_check <= primer_check + 4#(index number)
            condition4 = read_check >= primer_check - 4
            
            if  condition1 and condition2 and condition3 and condition4:
                # reads with correct primer location
                count += 1
                bam_list.write(str(read.query_name) + " ")
                bam_list.write(str(read.reference_start + 1) + " ")
                bam_list.write(str(read.reference_end) + "\n")
                primer_bam.write(read)
            else:
                if read.has_tag("XA"):
                    XA_tag = read.get_tag("XA")
                    multi_aln_list.write(str(read.query_name) + "\t" + str(XA_tag) + "\n")

        
        bam_file.close()
        bam_list.close()
        primer_bam.close()
        
        print("[PEM-Q] pass primer filter stitch reads: " + str(count))
        
        pysam.sort("-o", 'primer/'+self.p_bam_sort, 'primer/'+self.p_bam)
        cmd = "samtools index {}".format('primer/'+self.p_bam_sort)
        os.system(cmd)
        
        # no stitch reads        
        bam_file = pysam.AlignmentFile('bwa_align/'+self.no_stitch_r1_bam_sort, 'rb')
        primer_bam = pysam.AlignmentFile('primer/'+self.p_no_stich_r1_bam, "wb", template=bam_file)
        bam_list = open("primer/bamlist_nostitch.txt", "w")
        count = 0
        
        for read in bam_file:
            if primer_strand == 1:
                read_check = read.reference_start+1
                primer_check = primer_start
            else:
                read_check = read.reference_end
                primer_check = primer_end
            
            if read_check is None:
                continue
            
            if read.is_reverse:
                read_strand = -1
            else:
                read_strand = 1
            
            condition1 = read.reference_name == primer_chr
            condition2 = read_strand == primer_strand
            condition3 = read_check <= primer_check + 4 #(index number)
            condition4 = read_check >= primer_check - 4
            
            if  condition1 and condition2 and condition3 and condition4:
                count += 1
                bam_list.write(str(read.query_name) + " ")
                bam_list.write(str(read.reference_start + 1) + " ")
                bam_list.write(str(read.reference_end) + "\n")
                primer_bam.write(read)
            else:
                if read.has_tag("XA"):
                    XA_tag = read.get_tag("XA")
                    multi_aln_list.write(str(read.query_name) + "\t" + str(XA_tag) + "\n")

        bam_list.close()
        print("[PEM-Q] pass primer filter non-stitch reads: " + str(count))
        primer_bam.close()
        bam_file.close()
        multi_aln_list.close()

        pysam.sort("-o", 'primer/'+self.p_no_stich_r1_bam_sort, 'primer/'+self.p_no_stich_r1_bam)
        cmd = "samtools index {}".format('primer/'+self.p_no_stich_r1_bam_sort)
        os.system(cmd)
        
        return()
        
    def thread_1(self):
        # print("threading 1")
        self.align_genome()
        return()

    def thread_2(self):
        # print("threading 2")
        if self.adapter:
            
            try :
                self.align_adapter()
                self.exitcode1 = 1
            except:
                print("[Error] R2 fastq file is needed for adapter alignment.")
                self.exitcode1 = 0
                
        else:
            print("[main] Do NOT process adapter alignment!")
        return()
    
def main():
    
    start_time = time()
    
    args = docopt(__doc__,version='Align_make 1.0')
    
    kwargs = {'primerChrom':args['-r'], 'primer':args['-p'], 'adapter':args['-a'], 'genome':args['<genome>'], \
              'primer_start':args['-s'],'primer_end':args['-e'],'primer_strand':args['-d'],\
              'fastq_r1': args['<fastq_r1>'],'fastq_r2': args['<fastq_r2>']}
              
    print('[PEM-Q] primerChrom: ' + str(kwargs['primerChrom']))
    print('[PEM-Q] primer_start: ' + str(kwargs['primer_start']))
    print('[PEM-Q] primer_end: ' + str(kwargs['primer_end']))
    print('[PEM-Q] primer_strand: ' + str(kwargs['primer_strand']))
    print('[PEM-Q] primer: ' + str(kwargs['primer']))
    print('[PEM-Q] adapter: ' + str(kwargs['adapter']))
    print('[PEM-Q] genome: ' + str(kwargs['genome']))
    print('[PEM-Q] fastq_r1: ' + str(kwargs['fastq_r1']))
    print('[PEM-Q] fastq_r2: ' + str(kwargs['fastq_r2']))
    
    alignment = Align_make(**kwargs)
        
    func = [alignment.thread_1, alignment.thread_2]
    threads = []
    nthreads = range(len(func))
    
    for i in nthreads:
        t = threading.Thread(target=func[i])
        threads.append(t)
    for i in nthreads:
        threads[i].start()
    for i in nthreads:
        threads[i].join()
        
    if alignment.exitcode1 == 0:
        print("adapter alignment error! EXIT")
        exit()
        
    if kwargs['primer']:
        alignment.primer_check()
    else:
        print("[PEM-Q] Do NOT process primer filter!")
    
    print("\nAlign_make.py done in {}s".format(round(time()-start_time, 3)))
        
if __name__ == '__main__':
    main()
