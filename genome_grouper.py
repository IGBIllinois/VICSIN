#! /usr/bin/python3 -E

import argparse
import os
import shutil
import subprocess
import sys
import time

from multiprocessing import Pool

TEMPOUT = 'tempspineout'
TEMPIN = 'tempspinein'

def run_spine(ij):
    global input_files
    global input_dir
    global verbose
    global spine_path

    i, j = ij
    tempdir = '%s-%s' % (i, j)
    os.mkdir(tempdir)
    os.chdir(tempdir)

    iden1 = input_files[i]
    iden2 = input_files[j]
    
    if verbose:
        print('Running pair %s - %s' % (iden1, iden2))
    
    with open(TEMPIN, 'w') as f:
        f.write('%s\t%s\tfasta\n%s\t%s\tfasta' % (os.path.join(input_dir, input_files[i]),
            input_files[i], os.path.join(input_dir, input_files[j]), input_files[j]))
    
    with open('/dev/null', 'w') as f:
        proc = subprocess.Popen([spine_path, '-f', TEMPIN, '-o', TEMPOUT, '-t', '1'], stdout=f, stderr=f)
        proc.communicate()

    with open('%s.statistics.txt' % (TEMPOUT), 'r') as f:
        genome, genomesize = [], []

        for line in f:
            if line[0].isdigit():
                genome.append(line.split('\t')[1])
                genomesize.append(line.split('\t')[2])
            if line[0] == '-':
                backbone = line.split('\t')[4]

        data1 = float(backbone)/float(genomesize[0])
        data2 = float(backbone)/float(genomesize[1])

    os.chdir('..')
    shutil.rmtree(tempdir)

    return '%s-%s:%s\n%s-%s:%s' % (genome[0], genome[1], str(data1), genome[1], genome[0], str(data2))

if __name__ == '__main__':
    # track how long program takes to run
    start_time = time.time()

    # perform argument parsing
    parser = argparse.ArgumentParser(description='Write all-by-all table of pairwise Spine percent length aligned (PLA) values.', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('input_dir', help='directory containing FASTA formatted files for grouping')
    parser.add_argument('-v', '--verbose', help='increase output verbosity', action='store_true')
    parser.add_argument('-s', '--spine_path', help='spine binary path', default='spine')
    parser.add_argument('-n', '--num_cpus', type=int, help='number of cpus to use', default=1)
    parser.add_argument('-o', '--output', help='output path', default='genome_grouper_output.txt')
    args = parser.parse_args()

    # perform some input validation
    if not os.path.isdir(args.input_dir):
        print('ERROR: Input directory (%s) does not exist.' % args.input_dir)
        sys.exit(1)

    if args.spine_path != 'spine' and not os.path.exists(args.spine_path):
        print('ERROR: Invalid spine binary (%s).' % args.spine_path)
        sys.exit(1)

    if args.num_cpus < 1:
        print('ERROR: Must specify at least one CPU.')
        sys.exit(1)

    # setup variables that will be needed from within run_spine
    global input_files
    global input_dir
    global verbose
    global spine_path
    
    input_dir = args.input_dir
    verbose = args.verbose
    spine_path = args.spine_path
    input_files = os.listdir(input_dir)

    if len(input_files) < 2:
        print('ERROR: Need at least two FASTA files in input directory.')
        sys.exit(1)

    # create data matrix and data structures for mapping between index and file/genome name
    data = [[1.0] * len(input_files) for _ in range(len(input_files))]
    name_to_index = {}
    index_to_name = []

    for f in input_files:
        index_to_name.append(f)
        name_to_index[f] = len(index_to_name) - 1

    # create a list of all pairs of input files, without repeats
    inputs = []
    for i in range(len(input_files)):
        for j in range(i + 1, len(input_files)):
            inputs.append((i, j))

    # use multiprocessing to run multiple instances of spine on the input pairs depending
    # on the number of cpus
    with Pool(processes=args.num_cpus) as pool:
        processed_list = pool.map(run_spine, inputs)

    # extract spine data into data matrix
    for entry in processed_list:
        pair1, pair2 = entry.split('\n')
        key1, val1 = pair1.split(':')
        key2, val2 = pair2.split(':')
        genome1, genome2 = key1.split('-')
        data[name_to_index[genome1]][name_to_index[genome2]] = val1
        data[name_to_index[genome2]][name_to_index[genome1]] = val2
    
    # output the data matrix to a file
    with open(args.output, 'w') as f:
        f.write('\t')
        for genome in index_to_name:
            f.write('%s\t' % genome)
        f.write('\n')

        for i, row in enumerate(data):
            f.write('%s\t' % index_to_name[i])
            for el in row:
                f.write('%s\t' % el)
            f.write('\n')
    
    if verbose:
        print('Program ran for %s seconds' % (time.time() - start_time))

    sys.exit(0)
