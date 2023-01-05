############
# Script runs in the circulator conda environment with an install of seqkit and cd-hit installed.
############

#### import packages ####
import subprocess
import os
import operator
import csv
import concurrent.futures
import sys
import tempfile
import shutil
import argparse


def command_line_interface(args):
	parser = argparse.ArgumentParser(description='Comming soon!',
                                     add_help=False)
	parser.add_argument('-g',
	'--genomes',
	help='a path to input genomes',
	metavar='file.fasta',
	required=True,
	dest='input_genomes',
	nargs='+')
	
	parser.add_argument('-o',
	'--output',
	help='Path to output folder',
	required=True,
	dest='output_folder',
	type=str)

	if len(args) < 1:
		parser.print_help()
		exit(1)
	elif '-help' in args:
		parser.print_help()
		exit(1)

	args = parser.parse_args(args)

	return args


def split_genome_parts(input_genome_path, tmp_folder):
	first_genome_part = os.path.join(tmp_folder, os.path.basename(input_genome_path).split('.f')[0])

	input_args_stats = ['seqkit', 'split', '-s', '1', '-O', first_genome_part, input_genome_path]

	subprocess.run(input_args_stats, capture_output=True, check=True)

	first_genome_part_return = os.path.join(first_genome_part, os.path.basename(input_genome_path).split('.f')[0]) + '.part_001.fna'

	return(first_genome_part_return)


### fix start of genome to be DnaA
def circlator_fix_start(genome_file, tmp_folder, output_folder):
	# Initialise input list
	input_args = ['circlator', 'fixstart']

	# Add the input genome file to input list
	input_args.append(genome_file)
	
	# Construct output file name
	genome_basename = os.path.basename(genome_file)

	genome_basename = genome_basename.split('.')[0]

	output_file = os.path.join(tmp_folder, genome_basename) + '_fixedStart'

	if (os.path.isfile(os.path.join(output_folder, genome_basename) + '_fixedStart.fasta')):
		shutil.copyfile(os.path.join(output_folder, genome_basename) + '_fixedStart.fasta', output_file+'.fasta')
		return output_file

	# Add the output name to the input list
	input_args.append(output_file)

	subprocess.run(input_args)

	# TODO - Clean up files in tmp dir
	# Move fixed genome to output folder
	shutil.copyfile(output_file+'.fasta', os.path.join(output_folder, genome_basename) + '_fixedStart.fasta')
	
	return output_file

### Split genome
def split_genome(fixed_start_genome):
	# Initialise input list
	input_args_stats = ['seqkit', 'stats', '-b', fixed_start_genome+'.fasta']

	seqkit_split_return = subprocess.run(input_args_stats, capture_output=True, text=True, check=True)

	outline = seqkit_split_return.stdout
	outline = int(outline.split(' ')[-1].replace('\n', '').replace(',', ''))

	halfway = outline//2

	with open(fixed_start_genome+'.fasta', 'r') as input_file:
		# Get rid of header line
		input_file.readline()

		# Write first part of genome
		with open(fixed_start_genome+'_part_1.fasta', 'w') as file_1:
			file_1.write(f'>{os.path.basename(fixed_start_genome)}_part_1\n')
			for i in range(0,halfway//60):
				file_1.write(input_file.readline())
			
			remaining_line = input_file.readline()
			file_1.write(remaining_line[0:halfway%60]+'\n')

			# Write second part of genome
			with open(fixed_start_genome+'_part_2.fasta', 'w') as file_2:
				file_2.write(f'>{os.path.basename(fixed_start_genome)}_part_2\n')

				# Read and adjust the remaining genome to fit on 60 chr lines.
				remaing_genome = [remaining_line[halfway%60:]]
				remaing_genome = remaing_genome + list(input_file.readlines())
				remaing_genome = ''.join(remaing_genome)
				remaing_genome = remaing_genome.replace('\n', '')
				remaing_genome = [remaing_genome[i:i+60] for i in range(0, len(remaing_genome), 60)]

				for line in remaing_genome:
					file_2.write(line+'\n')
	
	# Reverse complement secound part of genome
	
	input_args_sliding_1 = ['seqkit', 'seq', '-t', 'DNA', '-r', '-p', fixed_start_genome+'_part_2.fasta', '-o', fixed_start_genome+'_part_2.fasta', '-o', fixed_start_genome+'_rev_part_2.fasta']
	
	try:
		subprocess.run(input_args_sliding_1, check=True, capture_output=True)
	except subprocess.CalledProcessError:
		print(f'Error in {fixed_start_genome}')
		exit(1)
	
	return(fixed_start_genome+'_part_1.fasta', fixed_start_genome+'_rev_part_2.fasta', halfway)
	

### Construct micro windows
def create_windows(file_1, file_2):
	input_args_sliding_1 = ['seqkit', 'sliding', '-g', '-s', '500', '-W', '3000', '-o', f'{file_1}_windows', file_1]
	input_args_sliding_2 = ['seqkit', 'sliding', '-g', '-s', '500', '-W', '3000', '-o', f'{file_2}_windows', file_2]

	subprocess.run(input_args_sliding_1)
	subprocess.run(input_args_sliding_2)

	part_1_names = subprocess.run(["seqkit", "seq", "-n", f"{file_1}_windows"], capture_output=True, text=True, check=True)
	part_2_names = subprocess.run(["seqkit", "seq", "-n", f"{file_2}_windows"], capture_output=True, text=True, check=True)

	part_1_names = part_1_names.stdout.split('\n')
	part_2_names = part_2_names.stdout.split('\n')

	part_1_names = list(filter(None, part_1_names))
	part_2_names = list(filter(None, part_2_names))

	return part_1_names, part_2_names


def seqkit_grep_sequences(seqeunce_list, output_folder, genome_name, window_file, set_num):
	grep_list_name = os.path.join(output_folder, genome_name) + '_grep.txt'
	with open(grep_list_name, 'w', newline='\n') as tmp_list:
		for id in seqeunce_list:
			tmp_list.write(id+'\n')

	read_sub_set = os.path.join(output_folder, genome_name) + '_sub_set_' + str(set_num) + '.fasta'

	seqkit_grep = ['seqkit', 'grep', '-f', grep_list_name, '-n', '-o', read_sub_set, window_file]

	subprocess.run(seqkit_grep, capture_output=True, check=True)

	return read_sub_set

### Cluster all seqeunces in macro windows
def cluster_windows(window_file_1, window_file_2, name_list_1, name_list_2, output_folder):
	# Calculate the number of seqeunces needed to cover a window of a set size:
	num_seq_to_clust = 30000 // 500

	# Find genome name
	genome_name = name_list_1[0].split('_fixedStart')[0]

	cluster_list_1 = [name_list_1[i:i+num_seq_to_clust] for i in range(0, len(name_list_1), int(num_seq_to_clust/3))]
	cluster_list_2 = [name_list_2[i:i+num_seq_to_clust] for i in range(0, len(name_list_2), int(num_seq_to_clust/3))]
	
	# Construct return list with intervals connected across replichores
	connected_intervals = []
	for index, sequence_set_1 in enumerate(cluster_list_1):

		read_sub_set_1 = seqkit_grep_sequences(sequence_set_1, output_folder, genome_name, window_file_1, 1)
		for sequence_set_2 in cluster_list_2:
			read_sub_set_2 = seqkit_grep_sequences(sequence_set_2, output_folder, genome_name, window_file_2, 2)

			cluster_out_name = os.path.join(output_folder, genome_name) + '_clst'

			# cd-hit to find similarities
			cdhit_args = ['cd-hit-est-2d', '-i', read_sub_set_1, '-i2', read_sub_set_2, '-o', cluster_out_name, '-d', '0', '-s', '0.8', '-c', '0.8', '-g', '0', '-n', '5']
			subprocess.run(cdhit_args, capture_output=True, check=True)

			# Find if any clusters that contain more than one seqeunce
			any_clustering = subprocess.run(['clstr2txt.pl', cluster_out_name+'.clstr'], check = True, capture_output=True, text=True)
			clstr_txt = any_clustering.stdout.split('\n')[1:-1]
			clstr_txt = [clstr_line.split('\t') for clstr_line in clstr_txt]

			clusters_oi = [int(clstr[1]) for clstr in clstr_txt if int(clstr[2]) > 1]

			# If any clusters with more than two seqeunces, find position for these
			if clusters_oi:
				# Extract all clusters with more than one seqeunce
				extracted_clusters = {key: [] for key in clusters_oi}
				for line in clstr_txt:
					if int(line[1]) in clusters_oi:
						extracted_clusters[int(line[1])].append([line[0], line[1]])

				# construct new list of lists with [Genome, right-coor, left-coor]
				for key in extracted_clusters:
					coordinates_part1 = extracted_clusters[key][0][0].split(':')[1].split('-')
					coordinates_part1 = [int(coor) for coor in coordinates_part1]

					coordinates_part2 = extracted_clusters[key][1][0].split(':')[1].split('-')
					coordinates_part2 = [int(coor) for coor in coordinates_part2]


					connected_intervals.append([genome_name, *coordinates_part1, *coordinates_part2])

	return connected_intervals



### Cluster overlaps and calculate relative positions
def overlap_n_relative(connected_intervals, half_genome_size):
	connected_intervals = sorted(connected_intervals, key=operator.itemgetter(1,3))

	index = 0
	merge_intervals = [connected_intervals[index]]
	for i in range(1, len(connected_intervals)):
        # If this is not first Interval and overlaps
        # with the previous one, Merge previous and
        # current Intervals
		# if (merge_intervals[index][2]+500 >= connected_intervals[i][1] and merge_intervals[index][4]+500 >= connected_intervals[i][3]): # TRY RUNNING WITH THE 500 removed!
		if (merge_intervals[index][2] >= connected_intervals[i][1] and merge_intervals[index][4] >= connected_intervals[i][3]): # TRY RUNNING WITH THE 500 removed!
			merge_intervals[index][2] = max(merge_intervals[index][2], connected_intervals[i][2])
			connected_intervals[index][2] = max(merge_intervals[index][2], connected_intervals[i][2])
			merge_intervals[index][4] = max(merge_intervals[index][4], connected_intervals[i][4])
			connected_intervals[index][4] = max(merge_intervals[index][4], connected_intervals[i][4])
		else:
			index = index + 1
			merge_intervals.append(connected_intervals[i].copy())


	merge_intervals_relative = []

	for interval in merge_intervals:
		genome = interval[0]

		relative_coords = [int(coor)/half_genome_size for coor in interval[1:]]
		# Round relative coordinates
		relative_coords = [round(coor, 3) for coor in relative_coords]

		genome_coor = interval[1:3]
		sec_replichore_coor = interval[3:][::-1]
		sec_replichore_coor = [half_genome_size*2 -coor for coor in sec_replichore_coor]
		genome_coor.extend(sec_replichore_coor)


		merge_intervals_relative.append([genome, *relative_coords, *genome_coor])

	# Claculate relative position for the merged intervals - loop through each merged interval and the coordinate within
	# merge_intervals_relative = [[line[0], *[round(coor/half_genome_size, 3) for coor in line[1:5]], *line[1:3], *[half_genome_size*2 -coor for coor in line[3:5][::-1]]] for line in merge_intervals]
	return(merge_intervals_relative)


# def write_output(output_list, output_folder):
# 	with open(os.path.join(output_folder, 'random_genome_connections_window_mode.tsv'), 'w') as output_file:
# 		output_file_writer = csv.writer(output_file, delimiter='\t')

# 		# Construct and write header
# 		header = ['Genome', 'relative_replichore_1_start', 'relative_replichore_1_end', 'relative_replichore_2_start', 'relative_replichore_2_end', 'replichore_1_start', 'replichore_1_end', 'replichore_2_start', 'replichore_2_end']
# 		output_file_writer.writerow(header)

# 		# Write result rows
# 		output_file_writer.writerows(output_list)


def extend_output(return_list, output_folder):
	new_file_created = False
	# Check if file exists
	if not os.path.isfile(os.path.join(output_folder, 'random_genome_connections_window_mode_extend.tsv')):
		new_file_created = True
		new_file = open(os.path.join(output_folder, 'random_genome_connections_window_mode_extend.tsv'), 'x')
		new_file.close()

	# Write new lines
	with open(os.path.join(output_folder, 'random_genome_connections_window_mode_extend.tsv'), 'a') as output_file:
		output_file_writer = csv.writer(output_file, delimiter='\t')

		if new_file_created:
			# Construct and write header
			header = ['Genome', 'relative_replichore_1_start', 'relative_replichore_1_end', 'relative_replichore_2_start', 'relative_replichore_2_end', 'replichore_1_start', 'replichore_1_end', 'replichore_2_start', 'replichore_2_end']
			output_file_writer.writerow(header)
			new_file_created = False

		# Write result rows
		output_file_writer.writerows(return_list)

def check_output(genomes, output_folder):
	print(f'Input genomes: {len(genomes) = }')
	genome_basename = [os.path.basename(genome_file) for genome_file in genomes]

	genome_basename = [genome.split('.')[0] for genome in genome_basename]


	if os.path.isfile(os.path.join(output_folder, 'random_genome_connections_window_mode_extend.tsv')):
		with open(os.path.join(output_folder, 'random_genome_connections_window_mode_extend.tsv'), 'r') as output_file:
			# Skip first line
			output_file.readline()
			for line in output_file:
				genome_to_exclude = line.split('\t')[0]
				genomes = [genome[0] for genome in zip(genomes, genome_basename) if genome[1].find(genome_to_exclude) == -1]

	print(f'Filtered genomes genomes: {len(genomes) = }')
	return(genomes)


def process_genome(input_genome_path, output_folder, tmp_folder):
	First_genome_part = split_genome_parts(input_genome_path, tmp_folder)

	fixed_start_genome = circlator_fix_start(First_genome_part, tmp_folder, output_folder)

	split_genome_1, split_genome_2, half_genome_size = split_genome(fixed_start_genome)

	part_1_names, part_2_names = create_windows(split_genome_1, split_genome_2)

	connected_intervals = cluster_windows(f"{split_genome_1}_windows", f"{split_genome_2}_windows", part_1_names, part_2_names, tmp_folder)
	
	merge_overlaps = []

	if connected_intervals:
		merge_overlaps = overlap_n_relative(connected_intervals, half_genome_size)

	return(merge_overlaps)


if __name__ == '__main__':
	cmd_args = command_line_interface(sys.argv[1:])

	tmp_folder_path = tempfile.TemporaryDirectory()

	genomes = check_output(cmd_args.input_genomes, cmd_args.output_folder)
	
	if (genomes):
		with concurrent.futures.ProcessPoolExecutor(max_workers=8) as executor:
			results = [executor.submit(process_genome, genome, cmd_args.output_folder, tmp_folder_path.name) for genome in genomes]
			
			res_index = 1
			for f in concurrent.futures.as_completed(results):
				if res_index % 50 == 0 or res_index == 1:
					print(f'\tGenome: {res_index} processed')
				return_list = f.result()
				if len(return_list):
					extend_output(return_list, cmd_args.output_folder)
				res_index += 1
	else:
		print("All input genomes in output file")
