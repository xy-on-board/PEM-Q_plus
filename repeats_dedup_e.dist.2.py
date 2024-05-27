#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#Author: Xinyi Liu

"""repeats_dedup

Usage:
    repeats_dedup [options] <file> <barcode_length>

Options:
-f feat, --feature_list=feat        example: -f 'Rname,Junction,Strand'
-h --help                           Show this screen.
-v --version                        Show version.

"""

from time import time
from docopt import docopts
import pandas as pd
import numpy as np
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

def repeats_dedup(file, bl, feature_list):
    
    file_tab = pd.read_csv("raw/"+file,sep = '\t', header=0)
    bl = int(bl)
    
    if file_tab['Qname'].count() > 0 :

        #~~~ step 1 ~~~#
        #pre-process#

        barcode_length = file_tab['Barcode'].str.len()
        file_tab['Barcode_len'] = barcode_length

        # drop reads with RMB shorter than set length -2
        file_tab = file_tab[barcode_length >=(bl-2)]

        # split barcode into single base
        barcode_list = list(map(str, range(1,(bl+1))))
        if file_tab.shape[0] > 1:
            file_tab[barcode_list] = file_tab['Barcode'].apply(lambda x: pd.Series(list(x)))
        if file_tab.shape[0] == 1:
            if len(file_tab['Barcode'].values[0]) < bl:
                file_tab['Barcode'] = file_tab['Barcode'].values[0] + 'N' * (bl - len(file_tab['Barcode'].values[0]))
                file_tab[barcode_list] = file_tab['Barcode'].apply(lambda x: pd.Series(list(x)))
            if len(file_tab['Barcode'].values[0]) == bl:
                file_tab[barcode_list] = file_tab['Barcode'].apply(lambda x: pd.Series(list(x)))

        # generate stable barcode and random barcode
        file_tab['stable_barcode'] = file_tab[['5', '9', '13']].astype(str).agg(''.join, axis=1)
        columns_to_merge = [str(i) for i in range(1, 18) if str(i) not in ['5', '9', '13']]
        file_tab['random_barcode'] = file_tab[columns_to_merge].apply(lambda x: ''.join(x.dropna().astype(str)), axis=1)
        file_tab['edit_dist_thresh'] = file_tab['stable_barcode'].apply(lambda x: 2 if x == 'TAT' else len([1 for a, b in zip(x, 'TAT') if a != b]) + 2)
        file_tab.loc[file_tab['Barcode_len'] < bl, 'edit_dist_thresh'] += bl - file_tab['Barcode_len']

        # generate 'Offset' and 'Bait_junction'
        file_tab.loc[:,'Offset'] = np.where(file_tab['Strand'] == "+", file_tab.Junction - file_tab.Qstart, file_tab.Junction + file_tab.Qlen - file_tab.Qend + 1)
        file_tab.loc[:,'Bait_junction'] = np.where(file_tab['Bait_strand'] == "+", file_tab.Bait_end, file_tab.Bait_start)

        # get barcode frequency
        counts = pd.value_counts(file_tab[u'Barcode'])
        counts = counts.to_frame()
        counts.reset_index(level=counts.index.names, inplace=True)
        counts.columns = ['barcode','counts']
        # add freq into the original file_tabframe
        file_tab['Barcode_freq'] = file_tab.groupby('Barcode')['Barcode'].transform('count')

        #~~~ step 2 ~~~#
        #dedup#

        # find dup
        dup = file_tab[file_tab.duplicated(feature_list, keep=False)]
        other = file_tab[~file_tab['Qname'].isin(dup['Qname'])]
        other = other.drop(columns=['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17'])

        dup = dup.drop_duplicates(feature_list + ['Barcode'], keep="first")

        total_columns = 17
        all_columns = [str(i) for i in range(1, total_columns + 1)]

        for i in range(total_columns - 2):
            for j in range(i + 1, total_columns - 1):
                current_drop_con = [col for col in all_columns if col not in [all_columns[i], all_columns[j]]]
                drop_con = feature_list + current_drop_con
                dup = dup.drop_duplicates(drop_con, keep="first")
                
        dup = dup.drop(columns=['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17'])
        

        feature_list = feature_list + ['R2_barcode']
        
        dup = dup.sort_values(by=feature_list)
        dup = dup.reset_index(drop=True)

        start_time = time()

        # dedup
        groups = dup.groupby(feature_list)
        print("before dedup by edit dist.:", len(dup))
        print("length of groups: " + str(len(groups)))

        result = []

        for name, group in groups:
                if len(group) == 1:
                    continue
                else:
                    group = group.sort_values(by=['Barcode_freq','Barcode_len','edit_dist_thresh'], ascending=[False, False, True])
                    similar_rows = find_similar_rows_highest(group)
                    result.extend(similar_rows)

        Qname_rm = pd.DataFrame({'Qname': result})

        print("after dedup by edit dist.:", len(dup))
        print ("\nRandom barcode dedup done in {}s".format(round(time()-start_time, 3)))


        dup = dup[~dup['Qname'].isin(Qname_rm['Qname'])]
        dup = dup.drop_duplicates(subset=feature_list + ['Barcode'])

        # clean reads (other + dup)
        output = pd.concat([other, dup])
        output = output.drop_duplicates(['Qname'])

        # reads removed
        # duplicates = rawdata[~rawdata['Qname'].isin(output['Qname'])]

        output.to_csv("unique/" + file, header = True, sep = '\t', index=False)
    else:
        file_tab.to_csv("unique/" + file, header = True, sep = '\t', index=False)
        
def main():
    start_time = time()
    
    args = docopt(__doc__,version='repeats_dedup 1.0')
    feature_list = args['--feature_list'].rsplit(sep=',')
    kwargs = {'file': args['<file>'], 'bl': args['<barcode_length>'], 'feature_list':feature_list}
    print('file: ' + str(kwargs['file']))
    print('feature_list: ' + str(kwargs['feature_list']))
    
    repeats_dedup(**kwargs)
    
    print("\nrepeats_dedup.py Done in {}s".format(round(time()-start_time, 3)))
    
if __name__ == '__main__':
    main()