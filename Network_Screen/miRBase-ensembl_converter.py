def make_mir_dic(hair_file):

    mir_dic = {}
    o_hair = open(hair_file)

    merged_data = []
    current_sequence = ""

    # Bit of a complicated loop system to stitch together the DNA code which is split over multiple lines. 

    for n, line in enumerate(o_hair):

        # if n > 5:
        #     break

        if line.startswith(">"):

            if current_sequence:

                merged_data.append(current_line + current_sequence)
                current_sequence = ""

            current_line = line + " "

        else:
            current_sequence += line.strip()

    # Append the last sequence
    if current_sequence:
        merged_data.append(current_line + current_sequence)

    # Print or use the merged data
    for n, line in enumerate(merged_data):

        split_line = line.split(' ')

        if not split_line[0].startswith('>hsa'):
            continue

        mir_dic[split_line[0]] = [split_line[1], split_line[-1]] 

    return mir_dic


def make_ensembl_dic(ens_file):

    mir_dic = {}
    o_ens = open(ens_file)

    merged_data = []
    current_sequence = ""

    for n, line in enumerate(o_ens):

        # if n > 5:
        #     break

        # print(line)

        if line.startswith(">"):

            if current_sequence:

                merged_data.append(current_line + current_sequence)
                current_sequence = ""

            current_line = line + " "

        else:
            current_sequence += line.strip()

    # Append the last sequence
    if current_sequence:
        merged_data.append(current_line + current_sequence)

    # Print or use the merged data
    for n, line in enumerate(merged_data):

        split_line = line.strip('>').replace('\n ', '|').split('|')
        # print(split_line)

        mir_dic[split_line[-1]] = [split_line[0], split_line[-2]] 

    return mir_dic


def mature_to_hairpin(mat_mirs, hair_mirs):

    keys_to_remove = []
    count = 0
    dic = {}

    for mat in mat_mirs:
        for hair in hair_mirs:
            if mat_mirs[mat][1] in hair_mirs[hair][1]:
                keys_to_remove.append(mat)
                # count +=1
                # print('MATCH' + str(count))
                # print(mat, hair)
                dic[mat] = [hair, hair_mirs[hair][1]]
                break

    for key in keys_to_remove:
        mat_mirs.pop(key)

    return dic


def mirbase_to_ensembl(mirbase_dic, ens_dic):

    keys_to_remove = []
    count = 0 
    dic = {}

    for n, mir in enumerate(mirbase_dic):

        # if n > 5:
        #     break

        seq = mirbase_dic[mir][-1].replace('U', 'T')

        if seq in ens_dic:
            keys_to_remove.append(mir)
            # count +=1
            # print('MATCH' + str(count))
            # print(mir, mirbase_dic[mir], ens_dic[seq])  
            dic[mir.strip('>')] = [ens_dic[seq][0], ens_dic[seq][1]]    

    for key in keys_to_remove:
        mirbase_dic.pop(key)

    # print('Remaining mirbase_dic')
    # pp(len(mirbase_dic))      
    # pp(mirbase_dic.keys())   
    return dic


def network_converter(network, convert_dic):

    o_net = open(network)
    outfile = open(network.replace('.', '_converted.'), 'w')
    first_line = True

    for n, line in enumerate(o_net):

        if first_line:
            outfile.write(line)
            first_line = False
            continue
        
        # if n >5:
        #     break

        split_line = line.split('\t')

        # print(split_line)

        if split_line[1] in convert_dic:
            # print(split_line[1])
            # print(convert_dic[split_line[1]])

            old_id, old_name = split_line[0], split_line[1]
            new_id, new_name = convert_dic[old_name][0], convert_dic[old_name][1]

            new_line = line.replace(old_id, new_id).replace(old_name, new_name)

            outfile.write(new_line)
        
        else:
            outfile.write(line)


if __name__ == "__main__":

    import argparse
    from pprint import pprint as pp

    parser = argparse.ArgumentParser()

    parser.add_argument('mature_file')
    parser.add_argument('hairpin_file')
    parser.add_argument('ensembl_mirs')
    parser.add_argument('network')

    args = parser.parse_args()

    mat_dic = make_mir_dic(args.mature_file)
    hair_dic = make_mir_dic(args.hairpin_file)

    ensembl_dic = make_ensembl_dic(args.ensembl_mirs)

    mirbase_dic = mature_to_hairpin(mat_dic, hair_dic)

    conversion_dic = mirbase_to_ensembl(mirbase_dic, ensembl_dic)

    network_converter(args.network, conversion_dic)


