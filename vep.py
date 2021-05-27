from subprocess import call
import subprocess
import datetime
import time

import glob
import pandas as pd


# load data
print('Loading data...\n\n\n')
f=open('hg19_v2/COVID19_HGI_2021.bed','r')
fc=f.read().split()

# make a dictionary with hg19 coordinates as keys and hg38 coordinates as items
hg_dic={}
for i in range(len(fc)):
    if (i+1)%4==0:
        hg_dic[fc[i-2]]=fc[i]

# make a dataframe using the GWAS txt file
df = pd.read_csv('GWASdata/COVID19_HGI_C2_ALL_eur_leave_ukbb_23andme_20210107.txt.gz_1.0E-5.txt',sep='\t')
df1=df[['#CHR','POS','REF','ALT','SNP','rsid']] # make a subset dataframe

# make the VEP input strings
df1['input']=0
for i in range(len(df1)):
    if len(df1['REF'][i])>1 or df1['ALT'][i]=='-':    # Deletion
        pos_dif=len(df1['REF'][i])-1
        df1['input'][i]=str(df1['#CHR'][i]) + ' ' + str(df1['POS'][i]) + ' ' + str(df1['POS'][i]+pos_dif) + ' ' + df1['REF'][i] + '/' + df1['ALT'][i]
    elif len(df1['ALT'][i])>1 or df1['REF'][i]=='-':  # Insertion
        pos_dif=len(df1['ALT'][i])-1
        df1['input'][i]=str(df1['#CHR'][i]) + ' ' + str(df1['POS'][i]) + ' ' + str(df1['POS'][i]-pos_dif) + ' ' + df1['REF'][i] + '/' + df1['ALT'][i]
    else:
        df1['input'][i]=str(df1['#CHR'][i]) + ' ' + str(df1['POS'][i]) + ' ' + str(df1['POS'][i]) + ' ' + df1['REF'][i] + '/' + df1['ALT'][i]

# put the concatenated input strings of variants that match to COVID19_HGI_2021 hg38 coordinates into a list
vep_input=[]
df1['hg38']=df1['#CHR'].map(str) + ':' + df1['POS'].map(str)
for i in hg_dic:
    for j in range(len(df1)):
        if hg_dic[i]==df1['hg38'][j]:
            vep_input.append(df1['input'][j])

# write the list into a txt file
output=open('vep_input.txt','w')
for element in vep_input:
     output.write(element)
     output.write('\n')
output.close()


# def Diff(li1, li2):
#     li_dif = [i for i in li1 + li2 if i not in li1 or i not in li2]
#     return li_dif
# Diff(list(df1['input']),vep_input)
# df1[df1['input']=='9 133270497 133270498 GA/G']
# df1[df1['input']=='9 133270637 133270638 AT/A']

# df1['input'] =df1['#CHR'].map(str) + ' ' + df1['POS'].map(str) + ' ' + df1['REF'].map(str) + '/' + df1['ALT'].map(str)
#
# for i in range(len(df1)):
#     if len(df1['REF'][i])>1 or len(df1['ALT'][i])>1:
#         print(i,df1.loc[[i]])
#
# df1['input'][471]=str(df1['#CHR'][471]) + ' ' + str(df1['POS'][471]) + ' ' + str(df1['POS'][471]+1) + ' ' + df1['REF'][471] + '/' + df1['ALT'][471]
