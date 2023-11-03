"""Analyse the edge weights for gene to gene interactions"""
from posixpath import split
from pprint import pprint as pp

def get_edge_per_TR(mapfile):

    from numpy import mean, median, std, percentile
    outfile = open(mapfile.replace('.', '-edge_weight_stats80th.'), 'w')
    open_map = open(mapfile)

    # again, there is no header, so switched to false - OC
    first_line = False

    tr_dic = {}
    tr_dic2 = {}

    type_ls = []

    for n, line in enumerate(open_map):

        # diagnostic purposes - OC
        # if n > 5:
        #     break

        if first_line:
            # print(line)
            first_line = False
            continue

        split_line = line.strip().split('\t')
        
        if split_line[2] not in type_ls:
            type_ls.append(split_line[2])

        # print(split_line[1], split_line[2])

        # had an issue where some entries were 'Protein coding' instead of 'protein_coding' so needed to set them straight so they wouldn't be filtered at next step - OC
        if 'Protein' in split_line[2] and 'coding' in split_line[2]:
            split_line[2] = 'protein_coding'

        if split_line[2] != 'protein_coding':
            continue

        # print(line)

        try:

            tr_dic[split_line[1]].append(float(split_line[-1]))
            tr_dic2[split_line[1]].append(float(split_line[7]))

        except KeyError:

            tr_dic[split_line[1]] = [float(split_line[-1])]
            tr_dic2[split_line[1]] = [float(split_line[7])]

    outfile.write('\t'.join(['TR', 'edge_count',
                             'mean_weight', 'median_weight', 'min_weight',
                             'max_weight', 'SD_weight', 'mean_distance', 'median_distance', 'min_distance', 'max_distance', 'SD_stiacen'])+'\n')

    perc=True

    if perc:

        for x in tr_dic:

            tr_perc = [mean(tr_dic[x])-(2*std(tr_dic[x])), mean(tr_dic[x])+(2*std(tr_dic[x]))]
            tr2_perc = [mean(tr_dic2[x])-(2*std(tr_dic2[x])), mean(tr_dic2[x])+(2*std(tr_dic2[x]))]
            tr_dic[x] = [i for i in tr_dic[x] if tr_perc[0] <= i <= tr_perc[1]]
            tr_dic2[x] = [i for i in tr_dic2[x] if tr2_perc[0] <= i <= tr2_perc[1]]

            if len(tr_dic[x]) == 0:

                print(x, tr_dic[x], tr_perc[0], tr_perc[1], '1')
                exit()

            elif len(tr_dic2[x]) == 0:

                print(x, tr_dic2[x], tr2_perc[0], tr2_perc[1], '2')
                exit()

    for x in tr_dic:

        outfile.write(
            '\t'.join([str(i) for i in [x, len(tr_dic[x]), mean(tr_dic[x]),
                                        median(tr_dic[x]), min(tr_dic[x]), max(tr_dic[x]), std(tr_dic[x]), mean(tr_dic2[x]),
                                        median(tr_dic2[x]), min(tr_dic2[x]), max(tr_dic2[x]), std(tr_dic2[x])]])+'\n')

        print(x, len(tr_dic[x]), mean(tr_dic[x]), median(tr_dic[x]), min(tr_dic[x]), max(tr_dic[x]))

    
    print(type_ls)

if __name__ == "__main__":

    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument('collapsed_network')

    args = parser.parse_args()

    get_edge_per_TR(args.collapsed_network)