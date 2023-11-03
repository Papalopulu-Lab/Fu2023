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

    return control_dic


def get_ensembl_ids(genome_file):

    """

    :param genome_file:
    :return:
    """
    first_line = True

    o_gene = open(genome_file)

    gen_dic = {}

    for line in o_gene:

        split_line = line.strip().split('\t')

        if first_line:

            id_dex = split_line.index('Gene stable ID')
            name_dex = split_line.index('Gene name')
            type_dex = split_line.index('Gene type')

            first_line = False

            continue

        try:

            gen_dic[split_line[name_dex]][1].append(split_line[type_dex])

            gen_dic[split_line[name_dex]][1] = list(set(gen_dic[split_line[name_dex]][1]))

        except KeyError:

            gen_dic[split_line[name_dex]] = [split_line[id_dex], [split_line[type_dex]]]

    return gen_dic


def get_subbed(subpath):

    import os

    subbed_file = subpath

    subbed_dic = {}

    if os.path.isfile(subbed_file):

        o_sub = open(subbed_file)

        for line in o_sub:

            split_line = line.strip().split('\t')

            if split_line[1] == 'INTERACTION REMOVED':

                continue

            subbed_dic[split_line[0]] = split_line[1]

        return subbed_dic, True

    else:

        return subbed_dic, False


def read_mirtar(mirtar_file, gen_dic):

    from os import path

    o_mirtar = open(mirtar_file)

    outfile = open(path.split(mirtar_file)[0] + '/mirtar_edges.txt', 'w')

    subbed_dic, overwrite = get_subbed(path.split(mirtar_file)[0] + '/mirtar_gene_name_substitutions.txt')

    if overwrite:

        sub_file = open(path.split(mirtar_file)[0] + '/mirtar_gene_name_substitutions.txt', 'w')

    else:

        sub_file = open(path.split(mirtar_file)[0] + '/mirtar_gene_name_substitutions.txt', 'a')

    first_line = True

    for line in o_mirtar:

        split_line = line.strip().split('\t')

        if first_line:

            id_dex = split_line.index('miRTarBase ID')
            miR_dex = split_line.index('miRNA')
            name_dex = split_line.index('Target Gene')
            qual_dex = split_line.index('Support Type')

            first_line = False

            continue

        id_mir = split_line[id_dex]
        mir = split_line[miR_dex]
        target_name  = split_line[name_dex]
        quality = split_line[qual_dex]
        ensembl_id = 0
        old_name = str(target_name)
        subbed = False
        skipped = False

        # try:
        #
        #     ensembl_id = gen_dic[target_name][0]
        #
        #
        #
        # except KeyError:
        #
        #     sub_file.write(old_name + '\tINTERACTION REMOVED\n')
        #
        #     continue

        # ''' This section will be included at a later date, will allow for updating of gene names


        
        while ensembl_id == 0:

            try:

                if target_name in subbed_dic and target_name not in gen_dic:

                    target_name = subbed_dic[target_name]

                ensembl_id = gen_dic[target_name][0]


            except KeyError:

                skipped = True

                break

                target_name = input('{} does not match an ensembl gene name\n\nPlease input a corrected name (type skip to remove): '.format(target_name))

                if target_name != 'skip':

                    subbed = True

                else:

                    subbed = False
                    skipped = True

                    break

                continue

        if subbed and skipped:

            exit('Subbed/skipped ERORR')

        elif subbed:

            sub_file.write(old_name+'\t'+target_name+'\n')
            subbed_dic[old_name] = target_name

        elif skipped:

            sub_file.write(old_name+'\tINTERACTION REMOVED\n')

            continue
        
        # '''

        gen_type = gen_dic[target_name][1]

        if len(gen_type) > 1:

            print(gen_type)

            if 'protein_coding' in gen_type:

                gen_type = 'protein_coding'

            else:

                gen_type = gen_type[0]

        else:

            gen_type = gen_type[0]

        outfile.write('\t'.join([id_mir, mir, 'miRNA', ensembl_id, target_name, gen_type, quality])+'\n')


    outfile.close()

    return 0


if __name__ == "__main__":

    import argparse
    from pprint import pprint as pp
    parser = argparse.ArgumentParser()

    # parser.add_argument('control_file')
    parser.add_argument('biomart_file')
    parser.add_argument('mirtar_file')

    args = parser.parse_args()

    # control_dic = get_file_paths(args.control_file)

    # pp(control_dic)

    gen_dic = get_ensembl_ids(args.biomart_file)

    read_mirtar(args.mirtar_file, gen_dic)

    exit()
