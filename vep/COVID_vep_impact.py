from subprocess import call
import subprocess
import datetime
import time

import glob
import pandas as pd

# read the liftOver file with hg19 and hg38 coordiantes
df=pd.read_csv('hg19_v2/COVID19_HGI_2021.bed',sep='\t',names=['chr','start','end','hg38_pos'])
# read the vep output file for A2_ALL_eur_leave_23andme
df1=pd.read_csv('vep_output_impact.txt',sep=' ',names=['vep_info','hg38_pos','ALT','impact_summary'])

# make a colum for impact and a column with hg38 first coordinate
# make reference to hg19 coordinate through the liftOver file
df1['hg38_pos_new']=0
df1['impact']=0
df1['hg19_pos']=0
for i in range(len(df1)):
    df1['impact'][i]=df1['impact_summary'][i].split(';')[0].split('=')[1]
    df1['hg38_pos_new'][i]=df1['hg38_pos'][i].split('-')[0]
    df1['hg19_pos'][i]=df[df['hg38_pos']==df1['hg38_pos_new'][i]]['start'].astype(int)

# remove duplicates to make counting easier
df1=df1.drop_duplicates(subset=['hg19_pos','impact'],keep=False,ignore_index=True)

# read all bed.gz files
msa_mapping=glob.glob('hg19_v2/*bed.gz')
msa_mapping.sort()

# read each file as a dataframe and store into a list
pd_list=[]
for i in range(len(msa_mapping)):
    df=pd.read_csv(msa_mapping[i],sep='\t',names=['chr','start','end','hg19_pos','score','strand'])
    species=msa_mapping[i].split('_')[2].split('.')[0]
    df['species']=species
    #df1=df[['chr','start','end','hg19_loc','species']]
    df=df.drop_duplicates()
    pd_list.append(df)

# concatenate all dataframe together
# reindex the dataframe
all=pd.concat(pd_list,ignore_index=True)

# use only the hg19 start position
all['hg19_pos_new']=0
for i in range(len(all)):
    all['hg19_pos_new'][i]=all['hg19_pos'][i].split('-')[0].split(':')[1]

# a dictionary for variant frequencies {IMPACT:{variant:freq}}
counts={'HIGH':{}, 'LOW':{}, 'MODERATE':{}, 'MODIFIER':{}}
for i in range(len(df1)):
    if df1['impact'][i]=='HIGH':
        count_value=all[all['hg19_pos_new']==df1['hg19_pos'][i]].shape[0]
        counts['HIGH'][df1['hg19_pos'][i]]=count_value
    if df1['impact'][i]=='LOW':
        count_value=all[all['hg19_pos_new']==df1['hg19_pos'][i]].shape[0]
        counts['LOW'][df1['hg19_pos'][i]]=count_value
    if df1['impact'][i]=='MODERATE':
        count_value=all[all['hg19_pos_new']==df1['hg19_pos'][i]].shape[0]
        counts['MODERATE'][df1['hg19_pos'][i]]=count_value
    elif df1['impact'][i]=='MODIFIER':
        count_value=all[all['hg19_pos_new']==df1['hg19_pos'][i]].shape[0]
        counts['MODIFIER'][df1['hg19_pos'][i]]=count_value

# write the dictionary into a file
with open("vep_impact_count.txt", 'w') as f:
    for key,value in counts.items():
        for variant,freq in value.items():
#        print(key,variant,freq,sep='\t')
            f.write('%s\t%s\t%s\n' % (key, variant,freq))
