"""Build network"""

from posixpath import split


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


def build_network(ccm_folder, mir_edges):

    import os

    ccm_file_ls = [i for i in os.listdir(ccm_folder) if 'ccm' in i]

    finished_map = open(ccm_folder.replace('map/', 'REMAP2022_map.txt'), 'w')

    for ccm in ccm_file_ls:

        o_ccm = open(ccm_folder+'/'+ccm)

        first_line = True

        for n, line in enumerate(o_ccm):

            if first_line:
  
                first_line = False

                continue

            # if n > 5:
            #     break

            split_line = line.strip().split('\t')

            regulator = split_line[3]
            reg_name = split_line[2]
            reg_type = split_line[4]
            target_name = split_line[10]
            target_code = split_line[8]
            target_type = split_line[12]
            experiment = split_line[6]
            cell_type = split_line[7]


            # had to make some slight changes to the slicing indeces here to reflect database changes - OC
            score = split_line[-2]
            act_distance = split_line[-4]
            abs_distance = split_line[-5]


            finished_map.write('\t'.join([regulator, reg_name, reg_type, target_code, target_name, target_type, score, act_distance, abs_distance, experiment, cell_type])+'\n')


        o_ccm.close()


    o_mir = open(mir_edges)

    first_line = True

    for n, line in enumerate(o_mir):

        # commenting out this line because there's no header in the mir_tar_edge file

        # if first_line:
        #     print(line)
        #     first_line = False

        #     continue

        # if n > 5:
        #     break

        split_line = line.strip().split('\t')
        # print(line)
        # print(split_line)
        

        regulator = split_line[0]
        reg_name = split_line[1]
        reg_type = split_line[2]
        target_code = split_line[3]
        target_name = split_line[4]
        target_type = split_line[5]
        score = 'None' # Dont think I need a score for miRTars think thats just a remap thing so setting to None - OC
        act_distance = '0' # These distances were already set to 0. Guess theyre not needed - OC.
        abs_distance = '0' 

        finished_map.write('\t'.join(
            [regulator, reg_name, reg_type, target_code, target_name, target_type, score, act_distance,
             abs_distance, 'None', 'None']) + '\n')

    o_mir.close()

    finished_map.close()


if __name__ == '__main__':

    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument('control_file')

    args = parser.parse_args()

    control_dic = get_file_paths(args.control_file)

    build_network(control_dic['MAP_FOLDER'], control_dic['MIRTAR_edge_file'])

