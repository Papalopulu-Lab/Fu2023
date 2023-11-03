from pprint import pprint as pp

def selective_collapse(mapfile):

    collapse_dic = {}
    outfile = open(mapfile.replace('.', '_select-collapsed.'), 'w')
    o_map = open(mapfile)
    head_line = ''

    first_line = True

    for line in o_map:

        if line.startswith('DONT_USE'):

            continue

        # elif first_line:
        #     first_line = False

        #     # dont think we need this line as there is no header anyway and it is just adding to the first data line. could add a header but im in a rush so will just remove it. OC
        #     # outfile.write(line.replace('\n', '\tcount\n'))
        #     continue

        split_line = line.strip().split('\t')

        key = '_'.join([split_line[1], split_line[4], split_line[9], split_line[10]])

        update_dic(collapse_dic, key, line, split_line[5], split_line[8])


    # dont think we need this line as there is no header anyway and it is just adding to the first data line. could add a header but im in a rush so will just remove it. OC
    # outfile.write(head_line.replace('\n', '\tcount\n'))

    # pp(collapse_dic)

    for x in collapse_dic:

        outfile.write(collapse_dic[x][0].replace('\n', '\t{}\n'.format(collapse_dic[x][1])))



def update_dic(dic, key, line, type, abs_distance):

    """

    :param dic: dictionary which to update
    :param key: key for the dictionary
    :param line: line of the file to be stored
    :param type: type of target
    :param abs_distance: peak distance from TSS
    :return dic: update dictionary

    """

    try:

        dic[key][1] += 1

    except KeyError:

        dic[key] = [line, 1, type, abs_distance]

        return dic

    current_type = dic[key][2]
    current_abs_distance = dic[key][3]

    if current_type != "protein_coding" and type == "protein_coding":

        dic[key][2] = type
        dic[key][3] = abs_distance
        dic[key][0] = line

        return dic

    elif current_type == "protein_coding" and type == "protein_coding":

        if abs_distance < current_abs_distance:

            dic[key][2] = type
            dic[key][3] = abs_distance
            dic[key][0] = line

            return dic

        else:

            return dic



if __name__ == "__main__":

    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument('map_file')

    args = parser.parse_args()

    selective_collapse(args.map_file)