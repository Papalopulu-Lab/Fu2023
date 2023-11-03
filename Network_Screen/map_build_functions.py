def assign_microrna(file_ls, out_dir, cmf_dir, tf_ensm_dic, container_ls, mir_host):

    import pandas as pd

    for file in file_ls:

        work_file = open('{}/{}'.format(cmf_dir, file))
        map_file = open('{}/{}-w_mir.txt'.format(out_dir, file), 'w')
        map_file.write(
            'chr\tchip_id\tchip_gene\tchip_gene_id\tgene_type\tchip_pos\texperiment\tcell_type\tgene_id\ttranscript_id\tgene_name\tstrand\tdistance\thost_type\n')

        line_count = 0

        print(container_ls)
        exit()

        for line in work_file:

            if line_count == 0:
                line_count += 1

                continue

            line = line.strip().split('\t')
            # print(line[2])
            # print(tf_ensm_dic[line[2]])
            ensmbl_code = '\t'.join(tf_ensm_dic[line[2]])

            # print(line)

            map_file.write('{}\t{}\t{}\tna\n'.format('\t'.join(line[0:2]), ensmbl_code, '\t'.join(line[3:])))

            # print(line[6])

            try: 
                
                check = container_ls[line[6]]

                print(line[4])

                match_sub = pd.DataFrame((mir_host.loc[(mir_host['gene'] == line[6]) & (mir_host['derives'] != 'self')]))
                mir_set = match_sub.drop_duplicates(subset='mir')
                mir_ls = list(mir_set['mir'])
                alias_ls = list(mir_set['alias'])
                derives_ls = list(mir_set['derives'])
                host_type = list(mir_set['gene_type'])

                for i in range(0, len(mir_ls)):
                    line[4] = derives_ls[i]
                    line[5] = alias_ls[i]
                    line[6] = mir_ls[i]
                    line[8] = "miRNA"

                    print('>>>>', line[6])
                    map_file.write('{}\t{}\t{}\t{}\n'.format('\t'.join(line[0:2]), ensmbl_code, '\t'.join(line[3:])), host_type)

            except KeyError:
                continue

        map_file.close()


def assign_microrna2019(file_ls, out_dir, cmf_dir, tf_ensm_dic, container_dic, mir_host):

    import pandas as pd

    for file in file_ls:

        work_file = open('{}/{}'.format(cmf_dir, file))
        map_file = open('{}/{}-w_mir.txt'.format(out_dir, file), 'w')
        map_file.write(
            'chr\tchip_id\tchip_gene\tchip_gene_id\tgene_type\tchip_pos\texperiment\tcell_type\tgene_id\ttranscript_id\tgene_name\tstrand\tdistance\thost_type\n')

        line_count = 0

        # print(container_dic)
        # exit()

        for line in work_file:

            if line_count == 0:
                line_count += 1
                continue

            line = line.strip().split('\t')
            print(line)
            print(line[2])
            # print(tf_ensm_dic[line[2]])
            ensmbl_code = '\t'.join(tf_ensm_dic[line[2]])

            # print(line)

            map_file.write('{}\t{}\t{}\tna\n'.format('\t'.join(line[0:2]), ensmbl_code, '\t'.join(line[3:])))

            # print(line[6])

            try: 
                
                check = container_dic[line[6]]

            except KeyError:
                continue
                # print(line[4])

            match_sub = pd.DataFrame((mir_host.loc[(mir_host['gene'] == line[6]) & (mir_host['derives'] != 'self')]))
            mir_set = match_sub.drop_duplicates(subset='mir')
            mir_ls = list(mir_set['mir'])
            alias_ls = list(mir_set['alias'])
            derives_ls = list(mir_set['derives'])
            host_type = list(mir_set['gene_type'])

            for i in range(0, len(mir_ls)):
                line[4] = derives_ls[i]
                line[5] = alias_ls[i]
                line[6] = mir_ls[i]
                line[8] = "miRNA"
                host = host_type[i]
                # print('>>>>', line[6])
                map_file.write('{}\t{}\t{}\t{}\n'.format('\t'.join(line[0:2]), ensmbl_code, '\t'.join(line[3:]), host))

            

        map_file.close()


def assign_microrna2022_ollie(file_ls, out_dir, cmf_dir, tf_ensm_dic, container_dic, mir_host):

    import pandas as pd
    import os

    key_error_ls =[]
    key_error_ls_file = os.path.join(os.path.split(cmf_dir)[0], 'excluded_genes.txt')
    error_file = open(key_error_ls_file, 'w')

    for file in file_ls:

        work_file = open('{}/{}'.format(cmf_dir, file))
        map_file = open('{}/{}-w_mir.txt'.format(out_dir, file), 'w')
        map_file.write(
            'chr\tchip_id\tchip_gene\tchip_gene_id\tgene_type\tchip_pos\texperiment\tcell_type\tgene_id\ttranscript_id\tgene_name\tstrand\tdistance\thost_type\n')

        line_count = 0

        # print(container_dic)
        # exit()

        for line in work_file:

            if line_count == 0:
                line_count += 1
                continue

            line = line.strip().split('\t')
            print(line)
            print(line[2])
            # print(tf_ensm_dic[line[2]])

            # Store key errors for chip genes not in ensembl into a list and continue - OC
            try:

                ensmbl_code = '\t'.join(tf_ensm_dic[line[2]])

            except KeyError:

                key_error_ls.append(line[2])
                continue

            # print(line)


            # tom originally had this as the first line, but this seemed to be excluding the gene name, i think extending to 0:3 should change this - OC
            # map_file.write('{}\t{}\t{}\tna\n'.format('\t'.join(line[0:2]), ensmbl_code, '\t'.join(line[3:])))
            map_file.write('{}\t{}\t{}\tna\n'.format('\t'.join(line[0:3]), ensmbl_code, '\t'.join(line[3:])))

            # print(line[6])

            try: 
                
                check = container_dic[line[6]]

            except KeyError:
                continue
                # print(line[4])

            match_sub = pd.DataFrame((mir_host.loc[(mir_host['gene'] == line[6]) & (mir_host['derives'] != 'self')]))
            mir_set = match_sub.drop_duplicates(subset='mir')
            mir_ls = list(mir_set['mir'])
            alias_ls = list(mir_set['alias'])
            derives_ls = list(mir_set['derives'])
            host_type = list(mir_set['gene_type'])

            for i in range(0, len(mir_ls)):
                line[4] = derives_ls[i]
                line[5] = alias_ls[i]
                line[6] = mir_ls[i]
                line[8] = "miRNA"
                host = host_type[i]
                # print('>>>>', line[6])
                map_file.write('{}\t{}\t{}\t{}\n'.format('\t'.join(line[0:2]), ensmbl_code, '\t'.join(line[3:]), host))

        map_file.close()

    # get unique genes from the key error list - OC
    unique_excl_genes = list(set([x for x in key_error_ls]))
    
    # write to error file list of all excluded genes - OC
    for gene in unique_excl_genes:
        error_file.write(f'{gene}\n')


def split_mir(vicious):

    import pandas as pd
    import multiprocessing as mp
    import os

    # split up file
    threads = mp.cpu_count()-1  # leaves thread available for other tasks
    v_file = open(vicious)
    out_dir = "{}_split".format(vicious.replace(".txt", ""))

    if not os.path.isdir(out_dir):

        os.makedirs(out_dir)


    v_count = 0

    for line in v_file:

        v_count += 1

    v_file.close()

    print("Vicious lines: ", v_count)

    remainder = (v_count % threads)
    div_thread = (v_count/threads) - (remainder/threads)

    check = div_thread * threads + remainder

    print("check: ", check)

    count_list = []

    for th in range(0, threads):

        count_list.append(div_thread)

        if remainder > 0:

            count_list[th] += 1

            remainder -= 1

        else:

            continue

    print("check2: ", sum(count_list))

    # print(count_list)

    v_file = open(vicious)

    file_list = []  # to save file locations as they are made

    for y in range(0, len(count_list)):

        counter = count_list[y]

        work_file = open("{}//{}.txt".format(out_dir, y), 'w')
        # print(work_file)
        file_list.append("{}//{}.txt".format(out_dir, y))

        for line in v_file:

            counter -= 1

            # print(counter)

            work_file.write('{}\n'.format(line.strip()))

            if counter == 0:

                break

        work_file.close()

    work_file.close()

    return file_list

def mirna_frame(file_list, genome_file, mir_dir):

    from mirBaseTools import mir_pos

    data_file = mir_pos('hsa')
    gen_file = open(genome_file)

    mir_dic = {}    # for storing microrna data for quick lookup the dic[alias] = [name, derives]

    for line in data_file:

        if 'primary_transcript' in line or line.startswith('#'):

            continue

        line = line.strip().split('\t')
        # print(line)
        line8 = line[8].split(';')

        alias = line8[1].replace('Alias=', '')
        name = line8[2].replace('Name=', '')
        derives = line8[3].replace('Derives_from=', '')

        mir_dic[alias] = [name, derives]

    mir_ls = []
    alias_ls = []
    derives_ls = []
    type_ls = []
    target_name_ls = []
    target_id_ls = []
    target_tranid_ls = []
    quality_ls = []
    target_type = []

    gen_dic = {}

    for line in gen_file:
        line = line.strip().split('\t')

        gen_dic[line[0]] = line[7]

    print('have gene dic')

    count = 0

    print(file_list)

    for filex in file_list:

        v_file = open(filex)

        mir_temp = open('{}/mir_temp{}.txt'.format(mir_dir, count), 'w')

        count+=1

        for line in v_file:  # cycle through the vicious file and create data such that it will fit into chip target frame

            line = line.strip().split('\t')
            # print(line)
            #
            # print(line[4])

            mirinfo = line[2].split(' ')
            target_info = line[0].split('|')

            mir_temp.write('{}\n'.format('\t'.join([mirinfo[0],
                                                    mirinfo[1],
                                                    mir_dic[mirinfo[1]][1],
                                                    line[4],
                                                    'microRNA',
                                                    target_info[2],
                                                    target_info[0],
                                                    target_info[1], gen_dic[target_info[0]]])))

        mir_temp.close()

    mir_temp.close()

    del mir_dic



def assign_mir_mir(file_ls, out_dir, container_ls, mir_host, microrna_dir):

    import pandas as pd

    for file in file_ls:

        work_file = open('{}/{}'.format(microrna_dir, file))
        map_file = open('{}/{}-w_mir.txt'.format(out_dir, file), 'w')
        print('{}/{}-w_mir.txt'.format(out_dir, file))

        line_count = 0

        for line in work_file:

            if line_count == 0:
                line_count += 1

                continue

            line = line.strip().split('\t')
            # print(line)
            ensmbl_code = line[6]
            # print(ensmbl_code)
            # print(line)

            map_file.write('{}\t{}\t{}\n'.format('\t'.join(line[0:3]), ensmbl_code, '\t'.join(line[3:])))

            # print(line[4])

            if line[6] in container_ls:

                # print(line[4])

                match_sub = pd.DataFrame((mir_host.loc[(mir_host['gene'] == line[6]) & (mir_host['derives'] != 'self')]))
                mir_set = match_sub.drop_duplicates(subset='mir')
                mir_ls = list(mir_set['mir'])
                alias_ls = list(mir_set['alias'])
                derives_ls = list(mir_set['derives'])

                for i in range(0, len(mir_ls)):
                    line[6] = derives_ls[i]
                    line[7] = alias_ls[i]
                    line[5] = mir_ls[i]
                    line[8] = "miRNA"
                    # print(line[6])
                    map_file.write('{}\t{}\t{}\n'.format('\t'.join(line[0:3]), ensmbl_code, '\t'.join(line[3:])))
        work_file.close()
        map_file.close()


def load_tf_to_remove(genome_file):

    import os

    work_path = os.path.split(genome_file)[0]

    multi_chip = open(work_path+'/unspec_chip.txt')

    unspecChipLs = []

    for line in multi_chip:

        if line.startswith('#'):

            continue

        unspecChipLs.append(line.strip())

    return unspecChipLs

def get_ensembl_ids(genome_file, tfbs_file, outfile):

        import pandas as pd

        outfile = open(outfile, 'w')
        conversionFile = outfile[0]

        chip_frame = pd.DataFrame(pd.read_csv(tfbs_file, sep='\t', header=None))

        columns_ls = ['chr', 'start', 'stop', 'gene', 'qual', 'strand', 'high_start', 'high_end', 'rgb', 'experiment', 'cell_line']

        chip_frame.columns = columns_ls[0:len(list(chip_frame.columns))]

        gen_ls = list((set(list(chip_frame['gene']))))

        del chip_frame

        gen_frame = pd.DataFrame(pd.read_csv(genome_file, sep='\t', header=0))

        gen_code_ls = []

        sub_len = 0

        for x in gen_ls:

            remove_ls = load_tf_to_remove(genome_file)

            if x in remove_ls:

                continue
        
            gen_sub = gen_frame.loc[gen_frame['Associated Gene Name'] == x]
            sub_ls = list(gen_sub['Ensembl Gene ID'])
            type_ls = list(gen_sub['Gene type'])

            if len(gen_sub) != 0:
                print('{}\t{}\t{}\n'.format(x, sub_ls[0], type_ls[0]))
                outfile.write('{}\t{}\t{}\n'.format(x, sub_ls[0], type_ls[0]))

            else:
                
                print('---------------\n\nNO MATCH\n\n')
                print('OLD GENE: {}\n'.format(x))
                gen_name = input('{}; New Gene Name: '.format(x))
        
                gen_sub = gen_frame.loc[gen_frame['Associated Gene Name'] == x]
                sub_ls = list(gen_sub['Ensembl Gene ID'])
                type_ls = list(gen_sub['Gene type'])

                if gen_name == 'IGNORE':
                    print(f'{x} has been excluded')
                    continue

                if len(gen_sub) != 0:

                    outfile.write('{}\t{}\t{}\n'.format(gen_name, sub_ls[0], type_ls[0]))

                else:

                    gen_code = input('{}; Gene ID: '.format(x))
                    gen_type = input('{}; Gene Type: '.format(x))

                    outfile.write('{}\t{}\t{}\n'.format(gen_name, gen_code, gen_type))

        outfile.close()

        del gen_frame


def get_ensembl_ids_2022_ollie(genome_file, tfbs_file, outfile):

        import pandas as pd

        outfile = open(outfile, 'w')
        conversionFile = outfile[0]

        chip_frame = pd.DataFrame(pd.read_csv(tfbs_file, sep='\t', header=None))

        print(chip_frame.head())

        columns_ls = ['chr', 'start', 'stop', 'gene', 'qual', 'strand', 'high_start', 'high_end', 'rgb', 'experiment', 'cell_line']

        chip_frame.columns = columns_ls[0:len(list(chip_frame.columns))]

        # list of gene names in remap chip file
        gen_ls = list((set(list(chip_frame['gene']))))

        del chip_frame


        gen_frame = pd.DataFrame(pd.read_csv(genome_file, sep='\t', header=0))

        gen_code_ls = []

        sub_len = 0

        for x in gen_ls:

            # think this is just removing TFs specified in unspec_chip.txt file if provided - OC
            remove_ls = load_tf_to_remove(genome_file)

            if x in remove_ls:

                continue

            # Ok can try something here will change Associated Gene Name to just Gene name and Ensembl Gene ID to to just Gene stable ID. As I can't see where these variable names are coming from?
            # Its possible Im loading the wrong file here and I should be using something other than the biomart file but I'm not sure...
            # I think its most likely that ensembl just changed the names they used for things in between me and tom constructing these networks - OC


            # basically if the gene from chip file is in biomart store it and ID and gene type information in temporary lists - OC
            gen_sub = gen_frame.loc[gen_frame['Gene name'] == x]
            # gen_sub = gen_frame.loc[gen_frame['Associated Gene Name'] == x]
            sub_ls = list(gen_sub['Gene stable ID'])
            # sub_ls = list(gen_sub['Ensembl Gene ID'])
            type_ls = list(gen_sub['Gene type'])

            # if gene exists in biomart then write to output file - OC
            if len(gen_sub) != 0:
                print('{}\t{}\t{}\n'.format(x, sub_ls[0], type_ls[0]))
                outfile.write('{}\t{}\t{}\n'.format(x, sub_ls[0], type_ls[0]))

            else:
                
                # No match found either find right match or say to ignore - OC
                print('---------------\n\nNO MATCH\n\n')
                print('\n\nSEARCH FOR GENE IN ENSEMBL AND INPUT ENSEMBL GENE NAME\n')
                print('\n\nTYPE IGNORE IF YOU DO NOT WISH GENE TO BE INCLUDED\n\n------------------\n')
                print('OLD GENE: {}\n'.format(x))
                gen_name = input('{}; New Gene Name: '.format(x))

                # then repeat above process to search biomart for newly provided name and add if matches -OC
                gen_sub = gen_frame.loc[gen_frame['Gene name'] == x]
                # gen_sub = gen_frame.loc[gen_frame['Associated Gene Name'] == x]
                sub_ls = list(gen_sub['Gene stable ID'])
                # sub_ls = list(gen_sub['Ensembl Gene ID'])
                type_ls = list(gen_sub['Gene type'])

                if gen_name == 'IGNORE':
                    print(f'{x} has been excluded')
                    continue

                if len(gen_sub) != 0:

                    outfile.write('{}\t{}\t{}\n'.format(gen_name, sub_ls[0], type_ls[0]))

                else:

                    # if no match then supply gene ID and type manually - OC

                    gen_code = input('{}; Gene ID: '.format(x))
                    gen_type = input('{}; Gene Type: '.format(x))

                    outfile.write('{}\t{}\t{}\n'.format(gen_name, gen_code, gen_type))

        outfile.close()

        del gen_frame


def get_ensembl_ids_jochen(genome_file, tfbs_file, outfile):

    import pandas as pd

    outfile = open(outfile, 'w')

    chip_frame = pd.DataFrame(pd.read_csv(tfbs_file, sep='\t', header=None))

    print(chip_frame.head())

    columns_ls = ['chr', 'start', 'stop', 'gene', 'qual', 'strand', 'high_start', 'high_end', 'rgb', 'experiment']

    chip_frame.columns = columns_ls[0:len(list(chip_frame.columns))]

    gen_ls = list((set(list(chip_frame['gene']))))

    del chip_frame

    gen_frame = pd.DataFrame(pd.read_csv(genome_file, sep='\t', header=0))

    gen_code_ls = []

    sub_len = 0

    for x in gen_ls:

        outfile.write('{}\t{}\t{}\n'.format(x, 'ENSG00000114315', 'protein_coding'))

    del gen_frame

# file_list = split_mir("F:/PhD - Nancy/Vicious/Outputs/ViciousHomo.txt")
# mirna_frame(file_list, 'D:/Dropbox (The University of Manchester)/ENCODE/HumanGenomeInfo.tsv')

