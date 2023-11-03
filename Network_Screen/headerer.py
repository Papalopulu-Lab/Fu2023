
network = '/home/y91198jk/hdd/Network/NetworkConstruction_2022/work_files/REMAP2022_map.txt'

with open(network) as open_net:

    outfile = open(network.replace('.', '_header.' ), 'w')

    outfile.write('regulator\tregulator_name\tregulator_type\ttarget\ttarget_name\ttarget_type\tscore\tact_distance\tabs_distance\texperiment\tcell_type\tcount\n')

    for n, line in enumerate(open_net):

        # if n > 5:
        #     break

        outfile.write(line)
