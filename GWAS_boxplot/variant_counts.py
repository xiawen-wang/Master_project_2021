from subprocess import call
import subprocess
import datetime
import time

#load data
print('Loading data...\n\n\n')
f=open('cord.txt','r')
fc=f.read().split()
hg19=[]
hg38=[]
for i in range(len(fc)):
    if i % 2 == 0:
        hg19.append(fc[i])
    else:
        hg38.append(fc[i])
print('Data loaded.\n\n-----------------------\n')

#check hg19 bed.gz files for variant frequencies
print('Checking variant frequencies...\n\n\nThe time now is ' + str(datetime.datetime.now()))
variant_freq=[]
for i in range(len(hg19)):
    zgrep = 'zgrep \'' + hg19[i] + '\' hg19/* | wc -l'
    variant_freq.append(int(subprocess.check_output(zgrep,shell=True))-1)
print('\n\nFinish checking variant frequencies.\n\n-----------------------\n ')

#prepare a dictionary for output
print('Making a dictionary for GWAS and frequency counts\n\n-----------------------\n\n ')
counts = {"A2_ALL_eur_leave_23andme": 0 , "A2_ALL_eur_leave_ukbb_23andme": 0 , "A2_ALL_leave_23andme": 0 , "A2_ALL_leave_UKBB_23andme":0 , "B1_ALL_eur_leave_23andme" : 0 , "B1_ALL_eur_leave_ukbb_23andme" : 0 , "B1_ALL_leave_23andme" : 0 , "B1_ALL_leave_UKBB_23andme" : 0 , "B2_ALL_eur_leave_23andme" : 0 , "B2_ALL_eur_leave_ukbb_23andme" : 0 , "B2_ALL_leave_23andme" : 0 , "B2_ALL_leave_UKBB_23andme" : 0 , "C2_ALL_eur_leave_23andme" : 0 , "C2_ALL_eur_leave_ukbb_23andme" : 0 ,"C2_ALL_leave_23andme" : 0 , "C2_ALL_leave_UKBB_23andme" : 0}
gwas=list(counts.keys())
gwas.sort()

# write frequency counts into the dictionary
print('Writing variant frequencies to corresponding GWAS...\n\n\nThe time now is ' + str(datetime.datetime.now()))
for i in range(len(hg38)):
    for j in range(len(gwas)):
        grep='grep \'' + hg38[i] + '\' GWASdata/* | grep \'' + gwas[j] + '\' | wc -l'
        a=int(subprocess.check_output(grep,shell=True))
        if a!=0:
            counts[gwas[j]]= counts[gwas[j]] + variant_freq[i]
print('\n\nFinish writing variant frequencies.\nFinish time: ' + str(datetime.datetime.now()))

# write the dictionary into a file
print('Writing output...\n\n')
with open("output_new.txt", 'w') as f:
    for key, value in counts.items():
        f.write('%s:%s\n' % (key, value))
print('Output written in output.txt.')
