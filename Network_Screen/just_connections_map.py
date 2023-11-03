# limits maps to just the TFs which regualte the network

from sys import argv
from pprint import pprint as pp

map_file = open(argv[1])

reg_dic = {}

tf_dic = {}

line_count = 0
head_line = 0

for n, line in enumerate(map_file):

    # if n > 5:
    #     break

    # print(line)

    if line_count == 0:

        line_count += 1
        head_line = line
        continue

    split_line = line.strip().split('\t')
    # print(split_line[4])

    try:

        reg_dic[split_line[4]].append(line)

    except KeyError:

        reg_dic[split_line[4]] = [line]


    tf_dic[split_line[1]] = 0

# pp(reg_dic)
# pp(tf_dic)

outfile = open(argv[1].replace('.', '-just_interactions.'), 'w')
outfile.write(head_line)

for key in list(tf_dic):

    # print(key)

    try:

        for x in reg_dic[key]:

            outfile.write(x)
            # print(x)

    except KeyError:

        continue

map_file.close()
outfile.close()
exit()