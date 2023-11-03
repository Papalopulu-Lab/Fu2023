from sys import argv

def update_map_dic(line, map_dic, reg_nom_dex, regulator_type_dex, targ_nom_dex, targ_type_dex):

    if line[reg_nom_dex] in map_dic[line[regulator_type_dex]]:

        if line[targ_nom_dex] not in map_dic[line[regulator_type_dex]][line[reg_nom_dex]][line[targ_type_dex]]:

            map_dic[line[regulator_type_dex]][line[reg_nom_dex]][line[targ_type_dex]][line[targ_nom_dex]] = 0

    else:

        map_dic[line[regulator_type_dex]][line[reg_nom_dex]] = {'protein_coding': {}, 'miRNA':{}}
        map_dic[line[regulator_type_dex]][line[reg_nom_dex]][line[targ_type_dex]][line[targ_nom_dex]]=0

    return map_dic

def read_map_and_fb(map_file):

    print('read_fb')
    o_map = open(map_file)

    fb_pro = {}
    fb_mir = {}

    map_dic = {'miRNA': {}, 'protein_coding': {}}

    first_line = True

    # count = 1

    for n, line in enumerate(o_map):

        line = line.strip().split('\t')
        # print(line)
        # print(count)
        # count += 1

        # print(n, line)

        if first_line:

            reg_dex, reg_nom_dex, regulator_type_dex, targ_dex, targ_nom_dex, targ_type_dex = line.index('regulator'), line.index('regulator_name'), line.index('regulator_type'), line.index('target'), line.index('target_name'), line.index('target_type')

            first_line = False
            print(reg_dex, reg_nom_dex, regulator_type_dex, targ_dex, targ_nom_dex, targ_type_dex)
            continue
        
        # had to add these two statements to correct for ocurrences of 'Protein coding' instead of 'protein_coding'
        if 'Protein' in line[targ_type_dex] and 'coding' in line[targ_type_dex]:
            line[targ_type_dex] = 'protein_coding'
        if 'Protein' in line[regulator_type_dex] and 'coding' in line[regulator_type_dex]:
            line[regulator_type_dex] = 'protein_coding'

        if line[targ_type_dex] != 'protein_coding' and line[targ_type_dex] != 'miRNA':

            continue

        elif line[reg_nom_dex] == line[targ_nom_dex] and line[regulator_type_dex] == 'protein_coding':

            fb_pro[(line[reg_nom_dex], line[regulator_type_dex])] = 0
            map_dic = update_map_dic(line, map_dic, reg_nom_dex, regulator_type_dex, targ_nom_dex, targ_type_dex)

        elif line[reg_nom_dex] == line[targ_nom_dex] and line[regulator_type_dex] == 'miRNA':

            fb_mir[(line[reg_nom_dex], line[regulator_type_dex])] = 0
            map_dic = update_map_dic(line, map_dic, reg_nom_dex, regulator_type_dex, targ_nom_dex, targ_type_dex)

        else:

            map_dic = update_map_dic(line, map_dic, reg_nom_dex, regulator_type_dex, targ_nom_dex, targ_type_dex)

    # print(map_dic)

    # exit()

    return map_dic, list(fb_pro), list(fb_mir)


def write_out_tup_ls(ls, outfile):

    ls = list(ls)

    o_out = open(outfile, 'w')
    count = {'gene': 0, 'motif': 0, 'gene2': 0}

    for x in ls:

        o_out.write('\t'.join(list(x))+'\n')

        count['gene'] += 1
        count['motif'] += 1

    count['gene2'] = 0

    o_out.close()

    return count

def write_out_dic(dic, outfile):

    o_out = open(outfile, 'w')
    # print(type(dic))

    count = {'gene': 0, 'motif': 0}

    gene2_ls = []

    for x in list(dic):

        # print(x)
        # print(dic[x])
        o_out.write('\t'.join([x, str(len(list(set(list(dic[x]))))), ' '.join(list(set(list(dic[x]))))])+'\n')

        count['gene'] += 1
        count['motif'] += len(list(set(list(dic[x]))))

        gene2_ls += list(set(list(dic[x])))

    count['gene2'] = len(list(set(gene2_ls)))

    o_out.close()

    return count

def tf_mir_tf(map_dic):

    print('tf_mir')
    pros = map_dic['protein_coding']
    mirs = map_dic['miRNA']

    tf_mir_dic = {}

    for tf in pros:
        # print(tf)
        for mirna in pros[tf]['miRNA']:

            if mirna not in mirs:
                print('--------- continue', mirna)
                continue

            if 'protein_coding' in mirs[mirna]:

                if tf in mirs[mirna]['protein_coding']:

                    if tf in tf_mir_dic:

                        if mirna in tf_mir_dic[tf]:
                            # print('--------- continue')
                            continue

                        else:

                            tf_mir_dic[tf].append(mirna)
                            print('---------', mirna)
                    else:

                        tf_mir_dic[tf] = [mirna]
                        print('---------', mirna)
            else:

                continue

    return tf_mir_dic

def tf1_tf2_tf1(map_dic):

    print('tf1_tf2')

    pros = map_dic['protein_coding']
    # print(pros)
    tf_tf2_dic = {}

    for tf in pros:

        for tf2 in pros[tf]['protein_coding']:
            # print(tf2)

            if tf2 not in pros or tf == tf2:

                continue

            if 'protein_coding' in pros[tf2]:

                if tf in pros[tf2]['protein_coding']:

                    if tf in tf_tf2_dic:

                        if tf2 in tf_tf2_dic[tf]:

                            continue

                        else:

                            tf_tf2_dic[tf].append(tf2)

                    else:

                        tf_tf2_dic[tf] = [tf2]

            else:

                continue

    return tf_tf2_dic


def amp_fb(fb_pro, tf_tf2_dic):

    print('amp_fb')

    amp_fb_dic = {}

    for tf_l in fb_pro:

        tf = tf_l[0]

        for x in list(tf_tf2_dic):

            if tf == x:

                # print(tf, 'match')
                amp_fb_dic[tf] = [y for y in tf_tf2_dic[x] if y != tf]

    return amp_fb_dic


def dual_fb(fb_pro, tf_tf2_dic):

    print('dual_fb')

    dual_fb_dic = {}

    tf_ls = [x[0] for x in fb_pro]

    for tf_l in fb_pro:

        tf = tf_l[0]

        if tf not in tf_tf2_dic:

            continue

        dual_fb_dic[tf] = [x for x in tf_tf2_dic[tf] if x in tf_ls and x != tf]

    return dual_fb_dic


def inc_fb_out(pro_fb, dic, outfile):

    o_out = open(outfile, 'w')

    count = {'gene': 0, 'motif': 0, 'microRNA': 0}
    micro_ls = []
    for x in pro_fb:

        try:

            o_out.write('\t'.join([x[0], x[1], str(len(list(set(list(dic[x[0]]))))), ' '.join(list(set(list(dic[x[0]]))))]) + '\n')

            count['gene'] += 1
            count['motif'] += len(list(set(list(dic[x[0]]))))
            micro_ls += list(set(list(dic[x[0]])))

        except KeyError:

            continue

    count['gene2'] = len(list(set(micro_ls)))

    o_out.close()

    return count


def out_summary(dic, outfile):

    o_out = open(outfile, 'w')

    for x in list(dic):

        o_out.write('\t'.join([x, str(dic[x]['gene']), str(dic[x]['motif']), str(dic[x]['gene2'])])+'\n')


def run_single(map_file):

    import os

    motif_dir = os.path.split(map_file)[0]+'/fast_motifs'

    if not os.path.isdir(motif_dir):

        os.makedirs(motif_dir)

    map_dic, fb_pro, fb_mir = read_map_and_fb(map_file)

    fb_pro_count = write_out_tup_ls(fb_pro, motif_dir+'/auto_regulate.txt')

    fb_mir_count = write_out_tup_ls(fb_mir, motif_dir + '/target_host.txt')

    tf_mir_dic = tf_mir_tf(map_dic)

    tf_mir_counts = write_out_dic(tf_mir_dic, motif_dir + '/tf-mir-tf.txt')

    inc_counts = inc_fb_out(fb_pro, tf_mir_dic, motif_dir + '/inc_fb.txt')

    tf_tf2_dic = tf1_tf2_tf1(map_dic)

    amp_dic = amp_fb(fb_pro, tf_tf2_dic)

    amp_counts = write_out_dic(amp_dic, motif_dir + '/amplified_fb.txt')

    dual_dic = dual_fb(fb_pro, tf_tf2_dic)

    dual_counts = write_out_dic(dual_dic, motif_dir + '/dual_fb.txt')

    print('inc_fb: ', inc_counts['gene'], inc_counts['motif'])
    print('pro_fb: ', fb_pro_count['gene'], fb_pro_count['motif'])
    print('mir_fb: ', fb_mir_count['gene'], fb_mir_count['motif'])
    print('amp_fb: ', amp_counts['gene'], amp_counts['motif'])
    print('dual_fb: ', dual_counts['gene'], dual_counts['motif'])
    print('tf_mir: ', tf_mir_counts['gene'], tf_mir_counts['motif'])

    summary_dic = {'incoherent_fb': inc_counts, 'protein_fb': fb_pro_count, 'targets_host': fb_mir_count, 'amplified_feedback': amp_counts, 'dual_feedback': dual_counts, 'tf_mir_loop': tf_mir_counts}


    out_summary(summary_dic, motif_dir + '/motif_summary.txt')
# print(len(argv))

if __name__ == '__main__':

    run_single(argv[1])

    exit()
