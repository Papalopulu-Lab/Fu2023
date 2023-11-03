def check_index(crm_file):
    from pprint import pprint as pp
    o_crm = open(crm_file)
    chr_dic = {}
    first_line = True
    chr_count = 0
    for line in o_crm:

        split_line = line.strip().split('\t')

        start= int(split_line[1])
        chr = split_line[0]

        try:

            chr_dic[chr].append(start)

        except KeyError:

            chr_dic[chr] = [start]

    for key in chr_dic:
        prevX = 0

        for x in chr_dic[key]:

            if x < prevX:

                exit('index inccorect: {} !> {} ()'.format(x, prevX, key))

            else:

                continue

    # print('\nCHRS: ', chr_count)


def get_crm_frame(crm_file):

    import pandas as pd

    crm_frame = pd.DataFrame(pd.read_csv(crm_file, header=None, sep='\t'))

    return crm_frame

def find_in_crm(crm_frame, chip_file):

    import pandas as pd

    o_chip = open(chip_file)

    crm_frame = pd.DataFrame(crm_frame)

    in_count = 0
    out_count = 0

    for line in o_chip:

        split_line = line.strip().split('\t')

        start = split_line[1]
        stop = split_line[2]
        middle = int(split_line[6])
        chr = split_line[0]

        chr_sub = crm_frame.loc[crm_frame[0] == chr]

        middle_sub = chr_sub.loc[(chr_sub[1] < middle) & (chr_sub[2] > middle)]

        if len(middle_sub) > 0:

            in_count += 1

        else:

            out_count += 1

def get_crm_dic(crm_file):

    ocrm = open(crm_file)
    crm_dic = {}

    for line in ocrm:

        split_line = line.strip().split('\t')

        chr = split_line[0]
        strt = int(split_line[1])
        stp = int(split_line[2])

        try:

            crm_dic[chr][0].append(strt)
            crm_dic[chr][1].append(stp)

        except KeyError:

            crm_dic[chr] = [[strt], [stp]]

    print(len(crm_dic))

    return crm_dic


def get_crm_dic_tupled(crm_file):

    ocrm = open(crm_file)
    crm_dic = {}

    for line in ocrm:

        split_line = line.strip().split('\t')

        chr = split_line[0]
        strt = int(split_line[1])
        stp = int(split_line[2])

        try:

            crm_dic[chr].append((strt, stp))


        except KeyError:

            crm_dic[chr] = [(strt, stp)]

    print(len(crm_dic))
    print(list(crm_dic))
    # exit()
    return crm_dic


def in_crm(start, stop, pos):

    if start <= pos <= stop:

        return True

    else:

        return False

def find_best_start(start_ls, pos):

    len_start_ls = len(start_ls)

    if len(start_ls) == 1:

        return start_ls[0]

    elif len_start_ls == 2:

        if start_ls[1][0] < pos:

            return start_ls[1]

        elif start_ls[0][0] < pos:

            return start_ls[0]

        else:

            return False

    mid_pos = int(len_start_ls/2)

    if start_ls[mid_pos][0] < pos:

        return find_best_start(start_ls[mid_pos:], pos)

    else:

        return find_best_start(start_ls[:mid_pos], pos)


def cycle_crms(crm_dic, chip_file):

    chip = open(chip_file)
    outfile = open(chip_file.replace('.', '-crm_filtered.'),'w')

    for line in chip:

        split_line = line.strip().split('\t')

        # print(split_line)

        middle = int(split_line[6])

        chr = split_line[0]

        # The following try except was added by OC, due to KeyErrors arising from chrs present in the chip file which are absent in the crm file.

        try:

            output = find_best_start(crm_dic[chr], middle)

        except KeyError:

            print('WARNING') 
            print('The chromosome ' + chr + ' is not in crm_file. This entry has been ignored.')
            continue

        if output:

            if in_crm(output[0], output[1], middle):

                # print('in crm', output,  line)
                outfile.write(line)
            else:

                # print('removed', line)
                continue
        else:

            # print('removed', line)
            continue



def count_crm_counts(crm):

    o_crm = open(crm)
    counter = 0
    for line in o_crm:

        split_line = line.strip().split('\t')
        # print(split_line)
        counter += int(split_line[4])

        print(counter, end='\r')

    print(counter)

if __name__ == "__main__":

    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument('CRM_file')
    parser.add_argument('chip_file')

    args = parser.parse_args()

    # check_index(args.CRM_file)
    #
    # check_index(args.chip_file)

    crm_dic = get_crm_dic_tupled(args.CRM_file)
    #
    # get_in_crm(crm_dic, args.chip_file)

    cycle_crms(crm_dic, args.chip_file)