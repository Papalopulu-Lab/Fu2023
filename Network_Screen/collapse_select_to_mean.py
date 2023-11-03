
def collapse_to_mean(select_collapse):

    from numpy import mean

    select_open = open(select_collapse)

    select_dic = {}

    outfile = open(select_collapse.replace('.', '_ctm.'), 'w')

    # again there isnt a header so set to false - OC
    first_line = False

    for line in select_open:

        if first_line:

            outfile.write(line)

            first_line = False

            continue

        split_line = line.strip().split('\t')

        key = '\t'.join(split_line[:6])

        try:

            select_dic[key][0].append(int(split_line[7]))
            select_dic[key][1].append(int(split_line[8]))
            select_dic[key][2].append(split_line[9])
            select_dic[key][3].append(split_line[10])
            select_dic[key][4].append(int(split_line[11]))

        except KeyError:

            select_dic[key] = [[int(split_line[7])], [int(split_line[8])], [split_line[9]], [split_line[10]], [int(split_line[11])]]


    for x in select_dic:

        dicfo = select_dic[x]

        min_dis = min(dicfo[1])

        min_dex = dicfo[1].index(min_dis)

        outfile.write(x+'\t'.join(['\t0', str(dicfo[0][min_dex]), str(dicfo[1][min_dex]), ' '.join(dicfo[2]), ' '.join(dicfo[3]), str(mean(dicfo[4]))])+'\n')




if __name__ == "__main__":

    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument("select_collapse")

    args = parser.parse_args()

    collapse_to_mean(select_collapse=args.select_collapse)
