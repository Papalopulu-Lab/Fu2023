def check_args(argv):
    '''This programme should utilise a multicore approach to assigning peaks in a memory efficient way by breaking the
    process up into chromosomes and seperating the files out as such so that only one needs to be loaded into RAM at a time.
    '''

    # step 1: Check we have all the files that we need

    if len(argv) < 4:
        print('\nIncorrect number of arguments provided, please check arguments and try again:\n\n'
              'Usage: python peakME.py genome_file chip_file output_directory')

        exit()


def make_dir(work_dir, gen_dir, chip_dir):

    import os

    if not os.path.isdir(work_dir):
        os.makedirs(work_dir)

    if not os.path.isdir(gen_dir):
        os.makedirs(gen_dir)

    if not os.path.isdir(chip_dir):
        os.makedirs(chip_dir)


    return [work_dir, gen_dir, chip_dir]


def make_chr_files(genome_file, chip_file, dir_ls, anzy):

    gen_dir = dir_ls[1]
    chip_dir = dir_ls[2]
    
    gen_open = open(genome_file)

    line_count = 0
    chr_ls = []
    file_dic = {}

    for line in gen_open:

        if line_count == 0:
            line_count += 1
            print(line)
            head_line = line
            head_split = line.strip().split('\t')

            print(head_split)

            if 'Chromosome Name' in head_split:

                chr_dex = head_split.index('Chromosome Name')
                trans_dex = head_split.index('Transcript type')
            else:

                chr_dex = head_split.index('Chromosome/scaffold name')

                trans_dex = head_split.index('Transcript type')

            continue



        split_line = line.split()
        chr = split_line[chr_dex]

        transcript_type = split_line[trans_dex]


        if chr not in chr_ls:

            chr_ls.append(chr)
            file_dic[chr] = open('{}/chr-{}.txt'.format(gen_dir, chr), 'w')
            file_dic[chr].write(head_line)

        file_dic[chr].write(line)

    gen_open.close()

    for chr in chr_ls:  # close all newly created files

        file_dic[chr].close()

    chr_ls = []
    file_dic = {}

    chip_open = open(chip_file)

    for line in chip_open:

        split_line = line.split()
        chr = split_line[0]

        if chr not in chr_ls:
            chr_ls.append(chr)
            file_dic[chr] = open('{}/chr-{}.txt'.format(chip_dir, chr.replace('chr', '')), 'w')
            # create new files for each choromosome

        file_dic[chr].write(line)

    chip_open.close()

    for chr in chr_ls:  # close all output files

        file_dic[chr].close()

def get_size(path_to_file):

    from os import path

    return path.getsize(path_to_file)


def match_chr(chip_dir, gen_dir):

    import os

    chr_ls = []
    size_ls = []

    chip_chr = [x for x in os.listdir(chip_dir) if 'chr' in x]
    gen_chr = [x for x in os.listdir(gen_dir) if 'chr' in x]

    for chr in chip_chr:

        if chr not in gen_chr:

            continue

        else:

            chr_ls.append(chr)

            # print(chr_ls)
            # print(chip_chr)
            # print(gen_chr)

    for chr in chr_ls:

        size_ls.append(get_size(gen_dir+'/'+chr))

    return chr_ls, size_ls


def split_load(threads, chr_ls):
    # chr_ls = ['chr-1.txt', 'chr-3.txt', 'chr-6.txt', 'chr-12.txt', 'chr-15.txt', 'chr-18.txt', 'chr-20.txt']
    split_dic = {}

    remainder = len(chr_ls) % threads

    division_size = (len(chr_ls) / threads) - (remainder / threads)

    # print(division_size)
    # print(remainder)

    remain_count = remainder
    ls_count = 0

    for x in range(1, threads+1):

        split_dic[x] = []

        for y in range(0, int(division_size)):

            split_dic[x].append(chr_ls[ls_count])
            ls_count += 1

        if remain_count != 0:

            split_dic[x].append(chr_ls[-(remain_count)])
            remain_count -= 1

        else:

            continue

    dic_ls = list(set(list(split_dic)))

    for key in dic_ls:

        print(key, '\t', split_dic[key])

    return split_dic

def dual_sort(list1, list2):
    
    list1, list2 = zip(*sorted(zip(list1, list2), reverse=True))

    return list(list1), list(list2)

def split_load2018(threads, chr_ls, size_ls):
    # chr_ls = ['chr-1.txt', 'chr-3.txt', 'chr-6.txt', 'chr-12.txt', 'chr-15.txt', 'chr-18.txt', 'chr-20.txt']
    split_dic = {}
    size_dic = {}
    pos_counter = 0

    size_ls, chr_ls = dual_sort(size_ls, chr_ls)

    for x in range(0, threads):

        split_dic[x] = []
        size_dic[x] = []

    switch = True

    for i in range(0, len(chr_ls)):

        split_dic[pos_counter].append(chr_ls[i])
        size_dic[pos_counter].append(size_ls[i])

        if pos_counter == threads - 1:

            switch = False

        elif pos_counter == 0:

            switch = True

        if switch:

            pos_counter += 1

        else:

            pos_counter -= 1

    for x in size_dic:

        print(x, sum(size_dic[x]), split_dic[x])

    # exit()

    return split_dic


def tfbs_id_gen(number, file_number):    # to generate a new number for each

    tfbs_id = 0
    tfbs_id_number = number

    tfbs_id_list = ['T', 'F', 'B', 'S', '-', str(file_number), '-', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0']

    tfbs_id_list_len = len(tfbs_id_list)

    tfbs_id_number_str = str(tfbs_id_number)

    # print('Size: ', len(tfbs_id_number_str))

    str_count = 1

    for i in tfbs_id_number_str:

        # print(tfbs_id_number_str[len(tfbs_id_number_str)-str_count])
        tfbs_id_list[tfbs_id_list_len-str_count] = tfbs_id_number_str[len(tfbs_id_number_str)-str_count]

        tfbs_id = ''.join(tfbs_id_list)
        str_count += 1

    return tfbs_id


def load_in(file_ls, gen_dir, chip_dir, up_range, down_range):

    import pandas as pd
    import os
    import numpy as np
    from datetime import datetime

    condition_file = open('{}/conditions.txt'.format(os.path.split(gen_dir)[0]),'w')

    condition_file.write('load_in peakfunctions {}\n'
                         'file_ls: {}\n'
                         'gene_directory: {}\n'
                         'chip_directory: {}\n'
                         'up_range: {}\n'
                         'down_range: {}\n'.format(datetime.now(), file_ls, gen_dir, chip_dir, up_range, down_range))

    close_dir = chip_dir.replace('chip_split', 'chip_close_match')
    all_dir = chip_dir.replace('chip_split', 'chip_all_match')

    if not os.path.isdir(close_dir):

        os.makedirs(close_dir)

    if not os.path.isdir(all_dir):

        os.makedirs(all_dir)



    for file in file_ls:

        # print(file)

        # if file != 'chr-15.txt':
        #
        #     continue

        closest_file = open('{}/{}'.format(close_dir, file.replace('.txt', '-ccm.txt')), 'w')
        all_file = open('{}/{}'.format(all_dir, file.replace('.txt', '-cam.txt')), 'w')

        headings = 'chr\tchip_id\tchip_gene\tchip_pos\tgene_id\ttranscript_id\tgene_name\tstrand\tgene_type\tabs_distance\tact_distance\tTSS\tscore\n'

        closest_file.write(headings)
        all_file.write(headings)

        chip_df = pd.DataFrame(pd.read_csv('{}/{}'.format(chip_dir, file), header=None, sep='\t'))
        gen_df = pd.DataFrame(pd.read_csv('{}/{}'.format(gen_dir, file), sep='\t'))

        # establish coding range

        all_trans_ls = list(gen_df['Gene start (bp)']) + list(gen_df['Gene end (bp)'])

        max_code = max(all_trans_ls) + up_range
        min_code = min(all_trans_ls) - up_range

        # print(chip_df.head())
        # exit()
        # remove rows not in range
        # print(len(chip_df))

        # print(len(chip_df))
        chr = chip_df.iat[0, 0]
        chip_columns = list(chip_df.columns)

        # print(len(chip_columns))

        columns_ls = ['chr', 'start', 'stop', 'TF', 'score', 'strand', 'mid_pos1', 'mid_pos2', 'rgb', 'experiment', 'cell_line']


        chip_df.columns = columns_ls


        del chip_df['mid_pos1']
        del chip_df['mid_pos2']
        del chip_df['rgb']


        # print(chip_df)

        # exit()

        chip_df['mid_point'] = (chip_df['stop'] + chip_df['start']) / 2

        chip_df = pd.DataFrame(chip_df.loc[chip_df['mid_point'] < max_code])
        chip_df = pd.DataFrame(chip_df.loc[chip_df['mid_point'] > min_code])
        # print(chip_df.head())

        # add ranges to genes

        gen_df['up_lim'] = np.where(gen_df['Strand'] == 1,
                                    list(gen_df['Transcription start site (TSS)'] - up_range),
                                    list(gen_df['Transcription start site (TSS)'] - down_range))


        gen_df['down_lim'] = np.where(gen_df['Strand'] == 1,
                                      list(gen_df['Transcription start site (TSS)'] + down_range),
                                      list(gen_df['Transcription start site (TSS)'] + up_range))

        pos_ls = list(chip_df['start'])
        chip_name_ls = list(chip_df['TF'])
        score_ls = list(chip_df['score'])
        experiment_ls = list(chip_df['experiment'])
        cell_ls = list(chip_df['cell_line'])

        # print(gen_df.head())
        # exit()

        for x in range(0, len(chip_df)):

            pos = pos_ls[x]
            tf = chip_name_ls[x]
            score = score_ls[x]
            experiment = experiment_ls[x]
            cellX = cell_ls[x]

            # print(gen_df.head())
            # exit()

            # if pos == 79456775:
            #
            #     print()

            # if pos == 96785795:
            #
            #     test_sub = gen_df.loc[(gen_df['down_lim'] > 96750000) & (gen_df['up_lim'] < 96850000)]
            #
            #     print(test_sub)

            gen_sub = pd.DataFrame(gen_df.loc[(gen_df['up_lim'] <= pos) & (gen_df['down_lim'] >= pos)])

            gen_sub['tss_dis'] = abs(gen_sub['Strand'] * (gen_sub['Transcription start site (TSS)'] - pos))
            gen_sub['tss_dis_real'] = gen_sub['Strand'] * (gen_sub['Transcription start site (TSS)'] - pos)
            #
            # if tf == 'CEBPB' and chr == 'chr15':
            #
            #     match_ls = list(gen_sub['Ensembl Gene ID'])
            #     # print(match_ls)
            #     transcript_ls = list(gen_sub['Ensembl Transcript ID'])
            #     gene_name_ls = list(gen_sub['Associated Gene Name'])
            #     strand_ls = list(gen_sub['Strand'])
            #     distance_ls = list(gen_sub['tss_dis'])
            #
            #     for x1 in range(0, len(match_ls)):
            #
            #         print('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(chr, tf, pos, match_ls[x1], transcript_ls[x1], gene_name_ls[x1], strand_ls[x1], distance_ls[x1]))


            # print(gen_sub.head())

            # print(list(gen_sub['tss_dis']))

            # print(gen_sub)

            gen_sub.sort_values('tss_dis', ascending=True, inplace=True)

            # if pos == 79456775:
            #
            #     pd.DataFrame.to_csv(gen_sub, '{}/HELP.txt'.format(all_dir), sep='\t', index=False)


            match_ls = list(gen_sub['Gene stable ID'])
            transcript_ls = list(gen_sub['Transcript stable ID'])
            gene_name_ls = list(gen_sub['Gene name'])
            strand_ls = list(gen_sub['Strand'])
            transcript_type_ls = list(gen_sub['Gene type'])
            distance_ls = list(gen_sub['tss_dis'])
            act_distance_ls = list(gen_sub['tss_dis_real'])
            tss_ls = list(gen_sub['Transcription start site (TSS)'])



            id_num = tfbs_id_gen(x, chr.replace('chr', ''))

            if len(match_ls) == 0:

                continue
                #
                # closest_file.write('{}\t{}\t{}\t{}\tno_match\tno_match\tno_match\tno_match\tno_match\n'.format(chr, 1, tf, pos))
                # all_file.write('{}\t{}\t{}\t{}\tno_match\tno_match\tno_match\tno_match\tno_match\n'.format(chr, 1, tf, pos))
                #
                # continue

            closest_file.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(chr, id_num, tf, pos, experiment, cellX, match_ls[0], transcript_ls[0], gene_name_ls[0], strand_ls[0], transcript_type_ls[0], distance_ls[0], act_distance_ls[0],tss_ls[0], score))

            # for x2 in range(0, len(distance_ls)):
            #
            #     all_file.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(chr, id_num, tf, pos, experiment, cellX, match_ls[x2], transcript_ls[x2], gene_name_ls[x2], strand_ls[x2], transcript_type_ls[x2], distance_ls[x2], act_distance_ls[0], tss_ls[0], score))

    condition_file.write('finish: {}'.format(datetime.now()))
    condition_file.close()
