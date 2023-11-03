def get_file_paths(control_file):

    """

    :param control_file:
    :return control_dic:

    parses control_file into a dictionary providing the code access to all the paths it needs to run

    control file format is @\tkey=path

    key is used as a dictionary key and the path is assigned to that key

    the dictionary is then returned for use by other functions

    """

    o_control = open(control_file)

    control_dic = {}

    for line in o_control:

        if line.startswith('@'):

            split_line = line.strip().split('\t')

            split1 = split_line[1].split('=')

            control_dic[split1[0]] = split1[1]

    print('LIST ', list(control_dic))

    return control_dic


if __name__ == '__main__':

    import multiprocessing as mp
    # from multiprocessing import Process as Process
    from peakME_functions_2018 import *
    import argparse
    import os

    parser = argparse.ArgumentParser()
    parser.add_argument('control_file')

    args = parser.parse_args()

    control_dic = get_file_paths(args.control_file)

    cwd = os.getcwd()
    threads = mp.cpu_count()
    print(threads)

    # if not os.path.isfile('{}/peakCORE.py'.format(cwd)):
    #
    #     print('\npeakCORE.py is required for peakME.py to run correctly please ensure you have both files in the same working directory')

    genome_file = control_dic['INFILE_DIRECTORY']+control_dic['GENOME_FILE']
    chip_file = control_dic['REMAP_FILE']  # chip_data in bed format
    output_dir = control_dic['NETWORK_DIRECTORY'] # directory to output information to

    # step 2: Make folder system
    work_dir = '{}/work_files'.format(output_dir)
    gen_dir = '{}/genome_split'.format(work_dir)
    chip_dir = '{}/chip_split'.format(work_dir)
    
    make_dir(work_dir, gen_dir, chip_dir)

    dir_ls = make_dir(work_dir, gen_dir, chip_dir)

    # step 3: split up the input files by chromosome

    make_chr_files(genome_file, chip_file, dir_ls, False)

    # we need to split load across the available hardware and the send it to another python file to call
    
    chr_ls, size_ls = match_chr(chip_dir, gen_dir)


    load_out = split_load2018(threads - 1, chr_ls, size_ls)

    load_out_ls = list(set(list(load_out)))

    '''--------------------------------------------------------------------------------------------------------------'''

    q = None

    jobs = []

    loop_count = 1

    for key in load_out_ls:

        key_ls = load_out[key]

        p = mp.Process(target=load_in, args=(key_ls, gen_dir, chip_dir, 50000, 10000))

        jobs.append(p)
        loop_count += 1

    for j in jobs:

        j.start()

    for j in jobs:

        j.join()

    print(loop_count)









