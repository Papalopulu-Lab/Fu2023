import multiprocessing as mp
# from multiprocessing import Process as Process
from peakME_functions_2022 import *
from map_build_functions import *
from sys import argv
import os
import pandas as pd

if __name__ == '__main__':

    threads = mp.cpu_count()
    cmf_dir = argv[1]   # close match files
    mir_file = argv[2]  # host file i.e hsa_host output from microME_plus.py
    genome_file = argv[3]   # ensembl genome file needs better description? -> humanGenomInfo.tsv
    outfile = argv[4]   # where to write tfbs_ids
    chip_file = argv[5]     # the chip_file lifted to current genome
    # get_ensembl_ids_2022_ollie(genome_file, chip_file, outfile) #this has been commented out as I already have a completed outfile. Will have to be redone for later database updates - OC
    
    # run first section of make_map

    out_gene_ls = []

    chip_ensmbl = open(outfile)
    out_dir = '{}/map'.format(os.path.split(cmf_dir)[0])
    mir_host = pd.DataFrame(pd.read_csv(mir_file, sep='\t', header=0))

    container_ls = set(mir_host['gene'])

    container_dic = {}
    
    for key in container_ls:
        container_dic[key] = 0

    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)

    # print('\nOutput directory: ', out_dir)

    tf_ensm_dic = {}

    for line in chip_ensmbl:

        if line.startswith('#'):

            continue

        line = line.strip().split('\t')
        tf_ensm_dic[line[0]] = line[1:]

        
    file_ls = os.listdir(cmf_dir)

    file_ls = [item for item in file_ls if not item.startswith('.')]

    test_mode = False
    # print(container_dic)
    # exit()

    if test_mode:

        assign_microrna2019(file_ls, out_dir, cmf_dir, tf_ensm_dic, container_dic, mir_host)
        exit()

    load_out = split_load(threads - 1, file_ls)

    load_out_ls = list(set(list(load_out)))

    # print(load_out_ls)

    '''--------------------------------------------------------------------------------------------------------------'''

    q = None

    jobs = []

    loop_count = 1

    for key in load_out_ls:

        # print(key)

        key_ls = load_out[key]

        p = mp.Process(target=assign_microrna2022_ollie, args=(key_ls, out_dir, cmf_dir, tf_ensm_dic, container_dic, mir_host))

        jobs.append(p)
        loop_count += 1

    for j in jobs:

        j.start()

    for j in jobs:

        j.join()

    print(loop_count)
