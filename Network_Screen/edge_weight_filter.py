def filter_map_mean(newtwork, filter_dic):

    import os

    # newtwork2022 = newtwork.replace('2018', '2022')
    outfile = open(newtwork.replace('.', 'EdgeFmean.'), 'w')

    type_ls_file = open(os.path.join(os.path.split(newtwork)[0], 'excluded_types_list.txt'), 'w')
    type_ls = []

    

    with open(newtwork) as open_net:

        # no header in network, set to false - OC
        first_line = False

        for n, line in enumerate(open_net):

            # if n > 1140:
            #     break

            if first_line:
                # print(line)
                first_line = False

                outfile.write(line)
                continue

            split_line = line.strip().split('\t')


            # had an issue where some entries were 'Protein coding' instead of 'protein_coding' so just adding this to set it straight - OC
            if 'Protein' in split_line[2] and 'coding' in split_line[2]:
                split_line[2] = 'protein_coding'
    

            #  ok had the issue where in the previous script edge_weight_analysis.py all types which were not protein_coding were excluded from the edge stats
            # this meant genes in the network with different types such as untranscribed_processed_pseudogene etc. were not in edge stats and were thus throwing key errors below
            #  tom said he excluded these genes so I am passing over those lines with a non-protein_coding gene regulator. For records sake I have recorded which types were kept and excluded in an output file
            try:

                if split_line[2] == 'miRNA':
                    outfile.write(line)

                    if split_line[2] not in type_ls:

                        type_ls.append(split_line[2])
                        type_ls_file.write(f'{split_line[2]}\tkept\n')

                    continue

                elif float(split_line[-1]) < filter_dic[split_line[1]]['mean_weight']:

                    continue

                else:

                    outfile.write(line)

            except KeyError:

                if split_line[2] not in type_ls:

                    type_ls.append(split_line[2])
                    type_ls_file.write(f'{split_line[2]}\texcluded\n')

                continue

            if split_line[2] not in type_ls:
                
                    type_ls.append(split_line[2])
                    type_ls_file.write(f'{split_line[2]}\tkept\n')

                
def filter_map_2SDmeanH(newtwork, filter_dic):

    import os


    """

    :param newtwork:
    :param filter_dic:
    :return:

    Filters the network removing anything below 2SD of the mean

    """

    outfile = open(newtwork.replace('.', 'EdgeFmean2SDoM.'), 'w')
 
    type_ls_file = open(os.path.join(os.path.split(newtwork)[0], 'excluded_types_list.txt'), 'w')
    type_ls = []   

    with open(newtwork) as open_net:

        # no header in network, set to false - OC
        first_line = False

        for n, line in enumerate(open_net):

            # if n > 1140:
            #     break

            if first_line:

                first_line = False

                outfile.write(line)

            split_line = line.strip().split('\t')

            # had an issue where some entries were 'Protein coding' instead of 'protein_coding' so just adding this to set it straight - OC
            if 'Protein' in split_line[2] and 'coding' in split_line[2]:
                split_line[2] = 'protein_coding'


            #  ok had the issue where in the previous script edge_weight_analysis.py all types which were not protein_coding were excluded from the edge stats
            # this meant genes in the network with different types such as untranscribed_processed_pseudogene etc. were not in edge stats and were thus throwing key errors below
            #  tom said he excluded these genes so I am passing over those lines with a non-protein_coding gene regulator. For records sake I have recorded which types were kept and excluded in an output file
            # I had originally written this for mean filter function above. Have adapted it for here. Hope it works. - OC May23

            try:

                if split_line[2] == 'miRNA':
                    outfile.write(line)

                    if split_line[2] not in type_ls:

                        type_ls.append(split_line[2])
                        type_ls_file.write(f'{split_line[2]}\tkept\n')

                    continue

                elif float(split_line[-1]) < filter_dic[split_line[1]]['mean_weight']-(2*filter_dic[split_line[1]]['SD_weight']):

                    continue

                else:

                    outfile.write(line)

            except KeyError:

                if split_line[2] not in type_ls:

                    type_ls.append(split_line[2])
                    type_ls_file.write(f'{split_line[2]}\texcluded\n')

                continue

            if split_line[2] not in type_ls:
                
                    type_ls.append(split_line[2])
                    type_ls_file.write(f'{split_line[2]}\tkept\n')


                


def filter_map_median(newtwork, filter_dic):

    outfile = open(newtwork.replace('.', 'EdgeFmed.'), 'w')

    with open(newtwork) as open_net:

        first_line = True

        for line in open_net:

            if first_line:

                first_line = False

                outfile.write(line)
                continue

            split_line = line.strip().split('\t')

            if split_line[2] == 'miRNA':
                outfile.write(line)
                continue

            elif float(split_line[-1]) < filter_dic[split_line[1]]['median_weight']:

                continue

            else:

                outfile.write(line)


def filter_map_1SDmeanH(newtwork, filter_dic):

    """

    :param newtwork:
    :param filter_dic:
    :return:

    Filters the network removing anything below 2SD of the mean

    """

    outfile = open(newtwork.replace('.', 'EdgeF_1D2SM.'), 'w')

    with open(newtwork) as open_net:

        first_line = True

        for line in open_net:

            if first_line:

                first_line = False

                outfile.write(line)

            split_line = line.strip().split('\t')

            if split_line[2] == 'miRNA':
                outfile.write(line)
                continue

            elif float(split_line[-1]) < filter_dic[split_line[1]]['mean_weight']-(filter_dic[split_line[1]]['SD_weight']):

                continue

            else:
                # print(filter_dic[split_line[1]]['mean_weight']-(filter_dic[split_line[1]]['SD_weight']))
                outfile.write(line)


def filter_map_2SDmeanBW(newtwork, filter_dic):

    """

    :param newtwork:
    :param filter_dic:
    :return:

    Filters the network removing anything below or above 2SD of the mean

    """

    outfile = open(newtwork.replace('.', 'EdgeF_2D2SM.'), 'w')

    with open(newtwork) as open_net:

        first_line = True

        for line in open_net:

            if first_line:
                first_line = False

                outfile.write(line)

            split_line = line.strip().split('\t')

            if split_line[2] == 'miRNA':
                outfile.write(line)
                continue

            elif float(split_line[-1]) < filter_dic[split_line[1]]['mean_weight'] - (
                    2 * filter_dic[split_line[1]]['SD_weight']):

                continue

            elif float(split_line[-1]) > filter_dic[split_line[1]]['mean_weight'] + (
                    2 * filter_dic[split_line[1]]['SD_weight']):

                continue

            else:

                outfile.write(line)


def load_filter_dic(edge_stats):

    edge_dic = {}
    head_line = []

    with open(edge_stats) as edge_open:

        first_line = True

        for line in edge_open:

            if first_line:

                first_line = False
                head_line = line.strip().split('\t')
                continue

            split_line = line.strip().split('\t')
            # print(split_line[0])
            edge_dic[split_line[0]] = {}

            for x in range(1, len(head_line)):

                edge_dic[split_line[0]][head_line[x]] = float(split_line[x])



    return edge_dic


if __name__ == "__main__":

    import argparse
    from pprint import pprint as pp

    parser = argparse.ArgumentParser()

    parser.add_argument('network')

    parser.add_argument('edge_stats')

    args = parser.parse_args()

    filter_dic = load_filter_dic(args.edge_stats)

    # pp(filter_dic)


    filter_map_2SDmeanH(args.network, filter_dic)
    # filter_map_2SDmeanBW(args.network, filter_dic)
    # filter_map_1SDmeanH(args.network, filter_dic)
    # filter_map_mean(args.network, filter_dic)
    # filter_map_median(args.network, filter_dic)
