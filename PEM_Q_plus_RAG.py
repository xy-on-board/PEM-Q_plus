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
    PEM-Q.py <genome> <sample> <cutsite> <primer_chr> <primer_start> <primer_end> <primer_strand> <primer>

Options:
-h --help               Show this screen.
-v --version            Show version.
<genome>                Reference genome(hg19/hg38/mm9/mm10).
<sample>                Sample name of input file (Note: input files should be named like: <sample>_R1.fq.gz and <sample>_R2.fq.gz).
<cutsite>               Position of cutsite (3' break end of positive strand).
<primer_chr>            Chromosome of red primer (eg:chr1).
<primer_start>          Start of red primer.
<primer_end>            End of red primer.
<primer_strand>         Strand of red primer (+/-).
<primer>                Sequence of red primer.

"""

import os
import sys
import threading
from time import time
from docopt import docopt

def run_script(sample=None, cutsite=None, genome=None, primer=None, primer_chr=None, primer_start=None, primer_end=None, primer_strand=None):
    
    start_time = time()
    basename = sample
    mq_threshold = 0

    print("######## 01 Reads alignment... ########")
    cmd = "align_make_plus_RAG.py {} {}_R1.fq.gz {}_R2.fq.gz -a CCACGCGTGCTCTACA  -p {} -r {} -s {} -e {} -d {}".format(genome, basename, basename, primer, primer_chr, primer_start, primer_end, primer_strand)
    print(cmd)
    os.system(cmd)

    cmd = "rm bwa_align/*sam bwa_align/*adpt.bam bwa_align/*sti.bam"
    os.system(cmd)
    
    print("######## 02 Barcode Extract... ########")
    cmd = "rmb_dedup_new_withRBB.py {} 17 CCACGCGTGCTCTACA ".format(basename)
    print(cmd)
    os.system(cmd)

    print("######## 03 Define transloc... ########")
    cmd = "define_transloc_3.0_newfil_RAG.py {} {} {} {} {}".format(basename, genome, cutsite, primer_strand, mq_threshold)
    print(cmd)
    os.system(cmd)

    cmd = "rm barcode/*sam barcode/*adpt.bam barcode/*filter.bam"
    os.system(cmd)

    print("######## 04 Define indels... ########")
    cmd = "define_indel_v5.1_mpf.py {} {} {}".format(basename, cutsite, primer_strand)
    print(cmd)
    os.system(cmd)

    print("######## 05 DEDUP... ########")
    cmd = "dedup_v6.2.py {} 17".format(basename)
    print(cmd)
    os.system(cmd)

    print("######## 06 substitutions... ########")
    cmd = "define_substitution.py {} {} {} {}".format(basename,cutsite,10,primer_strand)
    print(cmd)
    os.system(cmd)

    print("######## 07 filter and Statistics... ########")
    cmd = "define_statistics_add_filter_newDSBfilt_RAG.py {} {} {} {} {} {} {} CCACGCGTGCTCTACA".format(basename, genome, cutsite, 500000, primer, primer_chr, primer_strand)
    print(cmd)
    os.system(cmd)
    
    print("PEM-Q_plus_RAG done in {}s".format(round(time()-start_time, 3)))
    
def main():
    args = docopt(__doc__,version='PEM-Q v5.1s')
    
    kwargs = {'sample':args['<sample>'], 'cutsite':args['<cutsite>'],'genome':args['<genome>'],\
    'primer':args['<primer>'],'primer_chr':args['<primer_chr>'],'primer_start':args['<primer_start>'],\
    'primer_end':args['<primer_end>'],'primer_strand':args['<primer_strand>']}
    
    print('[PEM-Q] genome: ' + str(kwargs['genome']))
    print('[PEM-Q] sample: ' + str(kwargs['sample']))
    print('[PEM-Q] cutsite: ' + str(kwargs['cutsite']))
    print('[PEM-Q] primer: ' + str(kwargs['primer']))
    print('[PEM-Q] primer_chr: ' + str(kwargs['primer_chr']))
    print('[PEM-Q] primer_start: ' + str(kwargs['primer_start']))
    print('[PEM-Q] primer_end: ' + str(kwargs['primer_end']))
    print('[PEM-Q] primer_strand: ' + str(kwargs['primer_strand']))
    
    try:
        run_script(**kwargs)
    except KeyboardInterrupt:
        sys.stderr.write("See you~ :)\n")
        sys.exit(0)
    
if __name__ == "__main__":
    main()
