def build_mir_dic(mirbase):

    mir_dic = {}

    with open(mirbase) as o_mir:

        for line in o_mir:

            if line.startswith('>'):

                split_line = line.strip().split(' ')

                mir_dic[split_line[0].replace('>', '')] = split_line[1]

    # print(mir_dic)
    return mir_dic



def update_net(network_file, mir_dic):

    out_net = open(network_file.replace('.', '-fixed_mirs.'), 'w')

    with open(network_file) as open_net:

        # set to false, no header - OC
        first_line = False

        for n, line in enumerate(open_net):

            # if n > 15000000:

            #     break

            if first_line:

                out_net.write(line)
                first_line = False
                continue

            split_line = line.split('\t')

            if split_line[2] != "miRNA":

                out_net.write(line)

                continue



            else:
                # print(split_line[0], split_line[1])

                try:

                    split_line[0] = mir_dic[split_line[1]]
                    # print(mir_dic[split_line[1]])
                    out_net.write('\t'.join(split_line))

                except KeyError:

                    split_line[0] = 'unknown'

                    out_net.write('\t'.join(split_line))
                    print(f'unknown {split_line[1]}')

if __name__ == "__main__":

    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument("network")
    parser.add_argument("mirbase")

    args = parser.parse_args()

    update_net(args.network, build_mir_dic(args.mirbase))

    exit()
