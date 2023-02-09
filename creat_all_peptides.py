# -*- coding: utf-8 -*-
"""
Created on Mon Feb 14 22:15:43 2022

@author: user
"""

import pandas as pd
import csv
from functools import reduce
#exact path need r
with open (r'C:\\Users\\user\\OneDrive\\Desktop\\final_NGS_analysis\raw_data\barcode_table_final.csv','r',newline='') as f:
    reader=csv.reader(f)
    data = list(reader)
lista=data[0]
listb=data[1]
r='MTTC.....B.*'
path='C:\\\\Users\\\\user\\\\OneDrive\\\\Desktop\\\\final_NGS_analysis\\\\raw_data\\\\'
merge_df=[]
for i in range(0,29):
    exec('df_'+listb[i]+'='+'pd.read_table('+'\"'+path+lista[i]+'\"'+',header=None,delim_whitespace=True'+')')
    exec('refine_'+listb[i]+'='+'df_'+listb[i]+'.loc['+'df_'+listb[i]+'[0].str.contains('+'\''+r+'\''+')]')
    exec('r'+str(i)+'='+'\''+'({})'+'\''+'.format('+'\''+'|'+'\''+'.join('+'refine_'+listb[i]+'[0]'+'))')
    exec('refine_'+listb[i]+'[0]'+'='+'refine_'+listb[i]+'[0].str.extract('+'r'+str(i)+',expand=False).fillna('+'refine_'+listb[i]+'[0])')
    exec('sum_'+listb[i]+'='+'refine_'+listb[i]+'.groupby(0,as_index=False,sort=False).sum()')
    exec('sum_'+listb[i]+'='+'sum_'+listb[i]+'.set_axis(['+'\''+'Peptides'+'\''+','+'\''+listb[i]+'\''+'],axis=1)')
    exec('merge_df'+'='+'merge_df'+'+'+'['+'sum_'+listb[i]+']')

df_merged = reduce(lambda  left,right: pd.merge(left,right,on='Peptides',how='outer'), merge_df).fillna(float(0))
df_merged.to_csv('C:\\\\Users\\\\user\\\\OneDrive\\\\Desktop\\\\final_NGS_analysis\\\\allpeptides.txt',sep='\t',encoding="utf-8")
    
df_merged=pd.read_table('C:\\\\Users\\\\user\\\\OneDrive\\\\Desktop\\\\final_NGS_analysis\\\\allpeptides.txt',index_col=0)

'''
For truncate peptide longer than MTTCxxxxxBxxxxxCSWD to MTTCxxxxxBxxxxxCSWD,i.e., MTTCxxxxxBxxxxxCSWDxxxx truncated to MTTCxxxxxBxxxxxCSWD and sum the number
'''
df_merged['Peptides']=df_merged['Peptides'].str.extract(r'(MTTC.....B.....CSWD)',expand=False)
df_merged=df_merged.dropna()
df_merged=df_merged.groupby("Peptides",as_index=False,sort=False).sum()
df_merged.to_csv('C:\\\\Users\\\\user\\\\OneDrive\\\\Desktop\\\\final_NGS_analysis\\\\allpeptides_only_MTTC_B_CSWD.txt',sep='\t',encoding="utf-8")
