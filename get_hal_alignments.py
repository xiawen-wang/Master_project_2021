# Extract alignments from halformat

import os
from subprocess import call
import subprocess
import random
import sys
import datetime
import time
import re
sys.path.insert(0, '/exports/cmvm/datastore/sbms/groups/young-lab/rob/scripts/')
sys.path.insert(0, '/exports/cmvm/eddie/sbms/groups/young-lab/rob/scripts/')
from modules import *

#from Bio.Seq import Seq

job_id = sys.argv[1] # This is the job id, now can accept jobs of 1 to 10

# Quick test
hal_genome = 'hg19'

species = ['hg19', 'C57B6J', 'Rattus', 'micOch1', 'jacJac1', 'oryCun2', 'panTro4', 'gorGor3', 'ponAbe2', 'rheMac3', 'oviBos', 'oviAri3', 'bosTau8', 'canFam3', 'felCat8', 'loxAfr3']
species_names = ['human', 'mouse', 'rat', 'prairie_vole', 'egyptian_jerboa', 'rabbit', 'chimpanzee', 'gorilla', 'orangutan', 'macaque', 'musk_ox', 'sheep', 'cow', 'dog', 'cat', 'elephant']

# Rerun to add ffd.bed and second.bed
#promoter_infile = ls_list('*bed')[(int(job_id) - 1)]
#promoter_infile = ls_list('ANDNOT.bed')[0]
promoter_infile = ls_list('hg19/COVID19_HGI_2021.bed')[0]
genome = 'hg19'
genome_dir = genome
output_dir = genome

print job_id
print promoter_infile
print 'GENOME = ' + genome
print 'GENOMEDIR = ' + genome_dir
print 'HALGENOME = ' + hal_genome
print 'OUTDIR = ' + output_dir

replicate = str((int(job_id) % 10))
print replicate

# Make a temp directory if required
#if not os.path.exists('temp/'):
#	os.makedirs('temp')

progress = 0
with open(promoter_infile) as infile, open('temp/' + str(job_id) + '.bed', 'w') as outfile:
	for line in infile:

		if 'ID' not in line:
			progress += 1
			if ((progress % 10) == (int(replicate))):
				chrom = line.split()[0]
				start = line.split()[1]
				end = line.split()[2]

				midpoint = (int(start) + int(end))/2

				id = chrom + ':' + start + '-' + end

				outfile.write(bed_write(chrom = chrom, start = start, score = 0, end = end, name = id) + '\n')

list_infile = 'temp/' + job_id + '.bed'
total_wc = int(subprocess.check_output('wc -l ' + list_infile + ' | awk \'{print $1}\'', shell=True).rstrip())
print total_wc

# Test liftOver -> this takes ~ 4mins for 750 regions -> just over 1 hr in total
# Needs to be shrunk down to ~ 500 regions at a time
for target_genome in species:

	print target_genome

	complete_liftOver = 0

	if os.path.exists(output_dir + '/' + str(job_id) + '_' + target_genome + '.bed.gz'):
		complete_liftOver += 1

	if os.path.exists(output_dir + '/' + str(job_id) + '_' + target_genome + '_temp.bed'):
		complete_liftOver = 0
#		call('rm -f ' + output_dir + '/' + str(job_id) + '_' + target_genome + '_temp.bed',shell=True)

	if complete_liftOver == 0:

		for window in range(500,(total_wc+500),500):
			head_line = window

			get_temp = 'head -n ' + str(head_line) + ' temp/' + str(job_id) + '.bed | tail -n500 > temp/' + str(job_id) + '_temp.bed'
			print get_temp
			call(get_temp, shell=True)


			halLiftover = '/exports/cmvm/eddie/sbms/groups/young-lab/rob/scripts/hal/halLiftover /exports/cmvm/eddie/sbms/groups/young-lab/rob/genomes/1509_outgroups.hal ' + hal_genome + ' temp/' + str(job_id) + '_temp.bed ' + target_genome + ' ' + output_dir + '/' + str(job_id) + '_' + target_genome + '_temp.bed'
			print str(datetime.datetime.now())
			print halLiftover

			call(halLiftover, shell=True)
#			print str(datetime.datetime.now())

			if os.path.exists(output_dir + '/' + str(job_id) + '_' + target_genome + '.bed'):
				cat = 'cat ' + output_dir + '/' + str(job_id) + '_' + target_genome + '.bed ' + output_dir + '/' + str(job_id) + '_' + target_genome + '_temp.bed > ' + output_dir + '/' + str(job_id) + '_' + target_genome + '.bed2'
				print cat
				call(cat, shell=True)
				mv = 'mv ' + output_dir + '/' + str(job_id) + '_' + target_genome + '.bed2 ' + output_dir + '/' + str(job_id) + '_' + target_genome + '.bed'
				print mv
				call(mv, shell=True)

			else:
				mv = 'mv ' + output_dir + '/' + str(job_id) + '_' + target_genome + '_temp.bed ' + output_dir + '/' + str(job_id) + '_' + target_genome + '.bed'
				print mv
				call(mv, shell=True)

		call('rm -f ' + output_dir + '/' + str(job_id) + '_' + target_genome + '_temp.bed',shell=True)
		call('gzip ' + output_dir + '/' + str(job_id) + '_' + target_genome + '.bed', shell=True)

print 'Completed liftover'
exit()

##########################################

##### Create actual projections files ####

##########################################

# What about making a dictionary of projected coordinates?
print str(datetime.datetime.now())
print list_infile

if 'ffd' in promoter_infile:
	projection_outfile = output_dir + '/' + replicate + '_ffd_projections.bed'
	fasta_outfile = output_dir + '/' + replicate + '_ffd_projections.fa'

elif 'second' in promoter_infile:
	projection_outfile = output_dir + '/' + replicate + '_second_projections.bed'
	fasta_outfile = output_dir + '/' + replicate + '_second_projections.fa'
else:
	projection_outfile = output_dir + '/' + replicate + '_projections.bed'
	fasta_outfile = output_dir + '/' + replicate + '_projections.fa'

#else:
#	print promoter_infile
#	print 'ERROR: promoter infile not in expected format'
#	exit()

print projection_outfile + '.gz'

correct_compressed = 0

if os.path.exists(projection_outfile + '.gz'):

	line_count = subprocess.check_output('zcat ' + projection_outfile + '.gz | awk \'{print $4}\' | sort | uniq -c | awk \'{print $1}\' | sort | uniq -c | wc -l', shell=True).rstrip()
	if int(line_count) == 1:
		species_line_count = subprocess.check_output('zcat ' + projection_outfile + '.gz | awk \'{print $4}\' | sort | uniq -c | awk \'{print $1}\' | sort | uniq -c | awk \'{print $2}\'', shell=True).rstrip()
		if int(species_line_count) == 16:
			print 'compressed'
			correct_compressed = 1
			exit()

if correct_compressed == 0:

	if os.path.exists(projection_outfile + '.gz'):
		call('gunzip -f ' + projection_outfile + '.gz', shell=True)

	if os.path.exists(fasta_outfile + '.gz'):
		call('gunzip -f ' + fasta_outfile + '.gz', shell=True)

	if os.path.exists(projection_outfile):
		print projection_outfile
		print str(datetime.datetime.now())

		found_ids = bed_read('head -n -500 ' + projection_outfile + ' | grep human | awk \'{print $4}\'')

		count = len(found_ids)
		print count

		temp_fasta = open(fasta_outfile + '.temp', 'w')
		id = 'NA'
		infile = bed_read('head -n -10 ' + fasta_outfile) # don't lose the last line
		for line in infile:
			if '>' in line:
#				print line
				id = line.split(':')[1] + ':' + line.split(':')[2].split('-')[0] + '-' + line.split(':')[2].split('-')[1]

			if id in found_ids:
				temp_fasta.write(line + '\n')
#			else:
#				print line
#				exit()

		temp_fasta.close()
		print str(datetime.datetime.now())
		print 'line 241'
#		exit()

		temp_projection = open(projection_outfile + '.temp', 'w')
		infile = bed_read('head -n -10 ' + projection_outfile) # don't lost the last line
		for line in infile:
			id = line.split()[3]
			if id in found_ids:
				temp_projection.write(line + '\n')

		temp_projection.close()

		call('mv ' + projection_outfile + '.temp ' + projection_outfile, shell=True)
		call('mv ' + fasta_outfile + '.temp ' + fasta_outfile, shell=True)

		projection_output = open(projection_outfile, 'a')
		fasta_output = open(fasta_outfile, 'a')

	else:
		found_ids = ()
		projection_output = open(projection_outfile, 'w')
		fasta_output = open(fasta_outfile, 'w')

#found_ids = ()
#projection_output = open(projection_outfile, 'w')
#fasta_output = open(fasta_outfile, 'w')

print total_wc
print list_infile

#total_wc = 55

#chr7:16899318..16899336
#for step_index in range(52,(total_wc+1),50):
for step_index in range(1,(total_wc+1),50):
	projections = {}

	for line_index in range(step_index,step_index+50):

		if line_index <= total_wc:

			line = subprocess.check_output('head -n' + str(line_index) + ' ' + list_infile + ' | tail -n1', shell=True)
			id = line.split()[3]

			if id in found_ids:
				found = 1
				print 'Found ' + id + '!'
			else:

				for target_genome in species:

					liftover_file = output_dir + '/' + job_id + '_' + target_genome + '.bed.gz'
#			print liftover_file

#			liftover_file = 'temp/liftover/' + str(job_id) + '_' + target_genome + '.bed'
					check_liftover = int(subprocess.check_output('zcat ' + liftover_file + ' | grep -w ' + id + ' | sortBed | mergeBed -d 1 | wc -l', shell=True).rstrip())
					if check_liftover == 0:
						projection_result = 'GAP:GAP-GAP'
					else:
						top_results = subprocess.check_output('zcat ' + liftover_file + ' | grep -w ' + id + ' | sortBed | mergeBed -d 1 | head -n1', shell=True).rstrip()
						last_results = subprocess.check_output('zcat ' + liftover_file + ' | grep -w ' + id + ' | sortBed | mergeBed -d 1 | tail -n1', shell=True).rstrip()

						proj_chrom = top_results.split()[0]
						proj_start = top_results.split()[1]
						proj_end = last_results.split()[2]
						last_chrom = last_results.split()[0]
						projection_result = proj_chrom + ':' + proj_start + '-' + proj_end

						if proj_chrom != last_chrom:
							proj_chrom = 'MULTIPLE'
							proj_start = 'MULTIPLE'
							proj_end = 'MULTIPLE'
							projection_result = 'MULTIPLE:MULTIPLE-MULTIPLE'

					project_result = target_genome + ',' + projection_result

					if id in projections:
						if not re.search(target_genome, projections[id]): # there was previously a check for duplicates but these were allowed to make sure i got to the end of the file
							projections[id] = projections[id] + '\t' + project_result
					else:
						projections[id] = project_result

	count = str(len(projections))
	print 'Projections Count = ' + count

	tally = 0
	for id in projections:

		print id
		print projections[id]

		if 'ffd' in promoter_infile or 'second' in promoter_infile or 'H3' in promoter_infile:
			id_start = int(id.split(':')[1].split('-')[0])
			id_end = int(id.split(':')[1].split('-')[1])
		else:

			if re.search(hal_genome, projections[id]):
				source_projection = re.search(r'' + hal_genome + ',(\S+)', projections[id]).group(1)

				id_start = int(source_projection.split(':')[1].split('-')[0])
				id_end = int(source_projection.split(':')[1].split('-')[1])

		id_length = id_end - id_start
		projection_result = projections[id]
		tally += 1

		print projection_result

		if re.search(hal_genome, projection_result):
			found_projection = re.search(r'' + hal_genome + ',(\S+)', projection_result).group(1)

			print found_projection

			found_chrom = found_projection.split(':')[0]
			found_start = found_projection.split(':')[1].split('-')[0]
			found_end = found_projection.split(':')[1].split('-')[1]

			if found_start != 'GAP' or found_start != 'MULTIPLE':

				found_start = int(found_projection.split(':')[1].split('-')[0]) - 1
				found_end = int(found_projection.split(':')[1].split('-')[1]) - 1

			genuine_chrom = found_chrom # save these for comparisons at the end
			genuine_start = found_start
			genuine_end = found_end

		else:
			print 'line 341 error'
			exit()

		outfile = open('temp/' + job_id + '_maf.bed', 'w')
		print(bed_write(chrom = found_chrom, start = found_start, end = found_end, name = id))

# To hopefully recover the new way that hal2maf works
		outfile.write(str(found_chrom) + '\t' + str(found_start) + '\t' + str(found_end) + '\t' + str(id) + '\n')

#	outfile.write(bed_write(chrom = found_chrom, start = found_start, end = found_end, name = id) + '\n')
		outfile.close()

#	exit()
		length = int(found_end) - int(found_start)

		sys.stderr.write(id + '\n')
		sys.stderr.write(str(datetime.datetime.now()) + '\n')

		hal2maf = '/exports/igmm/eddie/taylor-lab/rob/scripts/hal/hal2maf /exports/igmm/eddie/taylor-lab/rob/genomes/1509_outgroups.hal temp/' + job_id + '.maf --maxRefGap ' + str(length) + ' --refGenome ' + hal_genome + ' --targetGenomes ' + str(','.join(species)) + ' --refTargets temp/' + job_id + '_maf.bed'

		print hal2maf
#	sys.stderr.write(hal2maf + '\n')
		call(hal2maf, shell=True)

	# Now check the coordinates - of the maf file, not the fasta

		correct_block = 0 # final check to make sure everything is ok
		with open ('temp/' + job_id + '.maf', 'r') as maf_infile, open('temp/' + job_id + '.maf2', 'w') as maf_outfile:
			for line in maf_infile:
				if '#' in line or line.isspace():
					if 'hal' not in line:
						maf_outfile.write(line)
				else:
					length = len(line.split())

					if length == 1 and line.split()[0] == 'a':
						correct_block = 0

					if length > 1:
						target_genome = line.split()[1].split('.')[0]
						maf_chrom = '.'.join(line.split()[1].split('.')[1:])
						maf_start = line.split()[2]
						maf_length = line.split()[3]
						maf_strand = line.split()[4]
						maf_srcLength = line.split()[5]
						maf_seq = line.split()[6]
						maf_line = line.rstrip()

						if maf_strand == '-':
							maf_start = str(int(maf_srcLength) - int(maf_start) - int(maf_length))

						maf_end = int(maf_start) + int(maf_length)

						if target_genome == hal_genome:

							if maf_strand == '-':
								ref_strand = '-'
								maf_seq = Seq(maf_seq)
								maf_seq = maf_seq.reverse_complement()

								if maf_strand == '+':
									maf_strand = '-'
								else:
									maf_strand = '+'
								maf_line = 's\t' + target_genome + '.' + maf_chrom + '\t' + maf_start + '\t' + maf_length + '\t' + maf_strand + '\t' + maf_srcLength + '\t' + str(maf_seq)

							else:
								ref_strand = '+'

							if re.search(target_genome, projection_result):
								found_projection = re.search(r'' + target_genome + ',(\S+)', projection_result).group(1)

#							print found_projection

								projection_chrom = found_projection.split(':')[0]
								projection_start = int(found_projection.split(':')[1].split('-')[0]) - 1
								projection_end = int(found_projection.split(':')[1].split('-')[1]) - 1

								if found_start == 'GAP' or found_start == 'MULTIPLE':
									found_length = id_length
								else:
									found_length = int(found_end) - int(found_start)

								if projection_chrom == maf_chrom:

#								print 'Line 190 = ' + str(line).rstrip()
#								print str(maf_start) + '\t' + str(projection_start) + '\t' + str(maf_end) + '\t' + str(projection_end)
#								print (str(int(maf_start) - int(projection_start))) + '\t' + str((int(projection_end) - int(maf_end)))
									if (int(maf_start) >= int(projection_start)) and (int(maf_end) <= int(projection_end)):

										maf_outfile.write('\na\n')
										maf_outfile.write(maf_line + '\n')
										correct_block += 1

							else:
								print line.rstrip()
								print target_genome
								print project_result
								print 'ERROR line 199'
								exit()

						else:

							if correct_block == 1:
								if ref_strand == '-':

									maf_seq = Seq(maf_seq)
									maf_seq = maf_seq.reverse_complement()
									maf_line = 's\t' + target_genome + '.' + maf_chrom + '\t' + maf_start + '\t' + maf_length + '\t+\t' + maf_srcLength + '\t' + str(maf_seq)

								if re.search(target_genome, projection_result):
									found_projection = re.search(r'' + target_genome + ',(\S+)', projection_result).group(1)

									projection_chrom = found_projection.split(':')[0]
									projection_start = found_projection.split(':')[1].split('-')[0]
									projection_end = found_projection.split(':')[1].split('-')[1]

									if projection_chrom != 'MULTIPLE' and projection_chrom != 'GAP':
										projection_start = int(found_projection.split(':')[1].split('-')[0]) - 1
										projection_end = int(found_projection.split(':')[1].split('-')[1]) - 1

									if found_start == 'GAP' or found_start == 'MULTIPLE':
										found_length = id_length
									else:
										found_length = int(found_end) - int(found_start)

									if projection_chrom == maf_chrom:

										if (int(maf_start) >= int(projection_start)) and (int(maf_start) <= int(projection_end)):

											maf_outfile.write(maf_line + '\n')

										else:
											if (int(maf_start) <= int(projection_start)) and (int(maf_end) >= int(projection_start)):

												maf_outfile.write(maf_line + '\n')

								else:
									print line.rstrip()
									print target_genome
									print project_result
									print 'ERROR line 225'
									exit()

	# Sort maf file
		maf_sort = 'maf-sort temp/' + job_id + '.maf2 | tail -n +2 > temp/' + job_id + '.maf'
		print maf_sort
		call(maf_sort, shell=True)

		maf2fa = '/exports/igmm/software/pkg/el7/apps/phast/1.3/usr/bin/msa_view --out-format FASTA temp/' + job_id + '.maf > temp/' + job_id + '.fa'

		print maf2fa
		call(maf2fa, shell=True)
#	call('rm -f temp/' + job_id + '.maf2', shell=True)

		check_fa_length_cmd = 'cat temp/' + job_id + '.fa | tr \'\\n\' \'\\t\' | tr \'>\' \'\\n\' | tail -n +2 | cut -f2- | tr -d \'\\t\' | awk -F \'\' \'{print NF}\' | sort | uniq | wc -l'
		check_fa_length = subprocess.check_output(check_fa_length_cmd, shell=True).rstrip()
		fa_alignment_lengths = len(check_fa_length)
		print 'fa_length\t' + id + '\t' + str(fa_alignment_lengths)

		fa_data = {}
		with open('temp/' + job_id + '.fa') as fa_infile:
			for line in fa_infile:
				if '>' in line:
					fa_species = line.split('> ')[1].split(':')[0].rstrip()

					if fa_species in fa_data:
						print 'FA DATA = ' + fa_data[fa_species]
						print 'LINE = ' + line
						print 'ERROR: folling fa_daa with ' + fa_species
						exit()
				else:

			# Remove lines with only *
					if len(line) > 1:
						line = line.replace('*', '-')

						if fa_species in species:
							if fa_species in fa_data:
								fa_data[fa_species] = fa_data[fa_species] + line.rstrip()
							else:
								fa_data[fa_species] = line.rstrip()
#								print fa_species


		species_progress = 0

		for target_genome in species:

			target_genome_name = species_names[species_progress]

			if re.search(target_genome, projection_result):
				found_projection = re.search(r'' + target_genome + ',(\S+)', projection_result).group(1)

				found_chrom = found_projection.split(':')[0]
				found_start = found_projection.split(':')[1].split('-')[0]
				found_end = found_projection.split(':')[1].split('-')[1]

				if found_start == 'GAP' or found_start == 'MULTIPLE':
					found_length = id_length
				else:
					found_start = int(found_start) - 1
					found_end = int(found_end) - 1
					found_length = int(found_end) - int(found_start)

			else:
				found_chrom = 'GAP'
				found_start = 'GAP'
				found_end = 'GAP'
				found_length = id_length

			if target_genome != hal_genome:

				if ((found_length > (2*id_length)) or (found_length < (id_length/2))):
					found_chrom = 'GAP'
					found_start = 'GAP'
					found_end = 'GAP'

			if found_start == 'GAP' or found_start == 'MULTIPLE':
				projection_output.write(bed_write(chrom = found_chrom, start = found_start, end = found_end, name = id, score = target_genome + ':' + target_genome_name) + '\n')
			else:
				found_alignment = 0

				if target_genome in species:
					if target_genome in fa_data:

#					print 'Found genome = ' + target_genome

				# pull out fa sequence for correct region, but what if there are multiple regions
						alignment_blocks = grep_file(target_genome, 'temp/' + job_id + '.maf')

						for line in alignment_blocks:
							if '#' not in line:

								alignment_chrom = '.'.join(line.split()[1].split('.')[1:])
								alignment_start = line.split()[2]
								alignment_length = line.split()[3]
								alignment_strand = line.split()[4]
								alignment_srcLength = line.split()[5]
								alignment_seq = line.split()[6]

								if alignment_strand == '-':
									alignment_start = str(int(alignment_srcLength) - int(alignment_start) - int(alignment_length))

								alignment_end = int(alignment_start) + int(alignment_length)

								if found_chrom == alignment_chrom:

									if int(alignment_start) >= int(found_start) and int(alignment_start) <= int(found_end):
											found_alignment += 1

									if int(alignment_start) <= int(found_start) and int(alignment_end) >= int(found_start):
											found_alignment += 1

						if found_alignment == 0:

							print line
							print 'TARGET = ' + target_genome
							print found_chrom + '\t' + alignment_chrom
							print str(alignment_start) + '\t' + str(found_start) + '\t' + str(alignment_end) + '\t' + str(found_end)
							print alignment_length + '\t' + str((int(alignment_end) - int(found_end)))
							print found_alignment
							print 'PROBLEM 1'
							exit()

						projection_output.write(bed_write(chrom = found_chrom, start = found_start, end = found_end, name = id, score = target_genome + ':' + target_genome_name) + '\n')

						fasta_output.write('>' + target_genome + ':' + id + '\n')
						fasta_output.write(fa_data[target_genome] + '\n')

					else: # this is corresponding to not in fa_data, assuming alignment is all gaps
						projection_output.write(bed_write(chrom = 'GAP', start = 'GAP', end = 'GAP', name = id, score = target_genome + ':' + target_genome_name) + '\n')


			species_progress += 1
		fasta_output.write('\n')

fasta_output.close()
projection_output.close()

if 'ffd' in promoter_infile:
	call('gzip -f ' + output_dir + '/' + replicate + '_ffd_projections.bed', shell=True)
	call('gzip -f ' + output_dir + '/' + replicate + '_ffd_projections.fa', shell=True)

elif 'second' in promoter_infile:
	call('gzip -f ' + output_dir + '/' + replicate + '_second_projections.bed', shell=True)
	call('gzip -f ' + output_dir + '/' + replicate + '_second_projections.fa', shell=True)

else:
	call('gzip -f ' + output_dir + '/' + replicate + '_projections.bed', shell=True)
	call('gzip -f ' + output_dir + '/' + replicate + '_projections.fa', shell=True)

print 'complete'
exit()
