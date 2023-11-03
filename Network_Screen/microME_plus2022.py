# code for assigning microRNAs to their primary transcripts
# you will need biomart files for now, may implement sql or perl api subscript
# the appropriate microRNA file


from sys import argv
from sys import version_info
from mirBaseTools import mir_pos
import pandas as pd
import numpy
import os
from pprint import pprint as pp

if version_info[0] < 3:
    print("Must be using Python 3 or above")
    exit()

if len(argv) < 3:

    print('\n\n*******************************************************************************')

    print('\nIncorrect number of arguments entered.\n\nUsage:\tpython microME_plus.py species ensembl_structure_file (optional: miRNase_build_no. default=CURRENT) mirbase_gff3_file\n')

    print('\n\n*******************************************************************************')

    exit()

species = argv[1]
dir = os.path.split(argv[2])[0]
microrna = argv[3]

exon_est = 250

if not os.path.isdir('{}/mirHOST'.format(dir)):

    os.makedirs('{}/mirHOST'.format(dir))

outfile = '{}/mirHOST/{}_hosts.txt'.format(dir, species)
structure_file = '{}/mirHOST/{}_structure.txt'.format(dir, species)

outfile = open(outfile, 'w')
structure_file = open(structure_file, 'w')

# generate headers for intron and exon recording

outfile.write('chr\tmir\talias\tstart\tstop\tderives\tgene_type\tgene_length\tgene\ttranscript\ttranscript_length\tstrand\tstart_feature\tstart_feature_len\tstop_feature\tstop_feature_len\tstart_len\tstop_len\tgene_start\tgene_end\tgene_source\tfull_match\tnum_exons\tnum_introns\texon_mean\tintron_mean\tmax_exon\tmax_intron\tmax_segment\n')
structure_file.write('gene\ttranscript\tfeature\tfeature_number\tnumber_exons\tnumber_introns\tmax_intron\tmax_exon\tutr3_len\tutr5_len\tgene_type\tsize\thost\n')

# microrna = mir_pos(species) this contacts internet to download file, but I keep seeing the connection as refused so have just downloaded hsa.gff3 and provided it as an argument instead on line 30

print('\n\nPlease ensure genome data and mirBase files are on the same build.')

def update_dic(dic, chr):

    dic[chr] = {}
    dic[chr]['name'] = []
    dic[chr]['start'] = []
    dic[chr]['stop'] = []
    dic[chr]['strand'] = []
    dic[chr]['alias'] = []
    dic[chr]['derives'] = []

micro_rna_dic = {}
chr_ls = []

# print(len(microrna))

with open(microrna) as microrna:
    for n, line in enumerate(microrna):

        if line.startswith('#'):

            continue

        line = line.strip()     # removes carriage return
        line = line.split('\t')     # splits line to list

        chr = line[0]       # chr identifier

        if chr not in chr_ls:   # if this is the first time we see a chr it creates a dictionary entry for this chr

            chr_ls.append(chr)
            update_dic(micro_rna_dic, chr)

        info = line[8].split(';')

        micro_rna_dic[chr]['name'].append(info[2].replace('Name=', ''))
        micro_rna_dic[chr]['start'].append(line[3])
        micro_rna_dic[chr]['stop'].append(line[4])
        micro_rna_dic[chr]['strand'].append(line[6])
        micro_rna_dic[chr]['alias'].append(info[1].replace('Alias=', ''))

        if 'derives' in line[8].lower():

            der_ls = line[8].split(';')
            # print(der_ls)
            der = der_ls[3].replace('Derives_from=', '')

        else:

            der = 'self'

        micro_rna_dic[chr]['derives'].append(der)

        length = int(line[3])-int(line[4])

        # print(length)

# pp(micro_rna_dic)
# print((micro_rna_dic['chr1']['alias']))
# print((micro_rna_dic['chr1']['name']))

ensembl = pd.DataFrame(pd.read_csv(argv[2], sep='\t', header=0, dtype={'Chromosome/scaffold name': 'str', 'Gene start (bp)':'int', 'Gene end (bp)':'int'}))
ensembl['Transcript length'] = ensembl['Transcript end (bp)'] - ensembl['Transcript start (bp)']

# print(list(ensembl.columns))
# pp(ensembl)

loop_count = 0
tom = 0
len_list = []
type_list = []
transcript_list = []
len_pro_list = []
gene_id_list = []

for chr in chr_ls:

    loop_count += 1

    chr_sub = ensembl.loc[ensembl['Chromosome/scaffold name'] == chr.replace('chr', '')]

    for x in range(0, len(micro_rna_dic[chr]['alias'])):

        mir = micro_rna_dic[chr]['name'][x]
        alias = micro_rna_dic[chr]['alias'][x]
        start = micro_rna_dic[chr]['start'][x]
        stop = micro_rna_dic[chr]['stop'][x]
        mid = int((int(start)+int(stop))/2)
        derives = micro_rna_dic[chr]['derives'][x]

        strand = micro_rna_dic[chr]['strand'][x]

        if strand == '+':

            strand = 1

        else:

            strand = -1

        hit_sub1 = chr_sub.loc[(chr_sub['Transcript start (bp)'] < int(start)) & (chr_sub['Transcript end (bp)'] > int(start))]
        hit_sub2 = chr_sub.loc[(chr_sub['Transcript start (bp)'] < int(stop)) & (chr_sub['Transcript end (bp)'] > int(stop))]

        hit_sub = pd.DataFrame(pd.concat([hit_sub1, hit_sub2], axis=0))

        # pp(hit_sub.columns)

        hit_sub = hit_sub.loc[hit_sub['Strand'] == strand]

        pp(hit_sub)

        if len(hit_sub) == 0:
            outfile.write('{}\t{}\t{}\t{}\t{}\t{}\tnan\tnan\tnan\tnan'
                          '\tnan\tnan\tnan\tnan\tnan\tnan\tnan\t'
                          'nan\tnan\tnan\tnan\tN\tex\tin\tex\tin\tNA\tNA\tNA\n'.format(chr, mir, alias, start, stop, derives))

        # print(mir, mid, start, stop)
        # print(hit_sub)

        trans_ls = (list(set(list(hit_sub['Transcript stable ID']))))

        pp(trans_ls)

        trans_count = 0
        match_count = 0

        for trans in trans_ls:
            print(hit_sub.columns)
            trans_sub = hit_sub.loc[hit_sub['Transcript stable ID'] == trans]
            gene = list(trans_sub['Gene stable ID'])[0]

            exon_ls = list(set(list(trans_sub['Exon stable ID'])))
            gene_type = list(trans_sub['Gene type'])[0]


            if len(exon_ls) <= 1:

                if trans_count == len(trans_ls) and match_count == 0:

                    print('One Exon', trans_ls[0])

                    outfile.write('{}\t{}\t{}\t{}\t{}\t{}\tnan\tnan\tnan'
                                  '\tnan\tnan\tnan\tnan\tnan\tnan\t'
                                  'nan\tnan\tnan\tnan\tnan\tN\tex\tin\tex\tin\tNA\tNA\tNA\n'.format(chr, mir, alias, start,
                                                                                                 stop, derives))

                    continue

                continue

            trans_count += 1
            exon_rank = list(trans_sub['Exon rank in transcript'])
            exon_start = list(trans_sub['Exon region start (bp)'])
            exon_stop = list(trans_sub['Exon region end (bp)'])
            strand = list(trans_sub['Strand'])[0]
            # gen_source = list(trans_sub['Source of gene name'])[0]

            # print(exon_ls)
            # print(exon_rank)
            # print(exon_start)
            # print(exon_stop)
            # print('Strand: ', strand)

            utr5_start = list(trans_sub['5\' UTR start'])[0]
            utr5_stop = list(trans_sub['5\' UTR end'])[0]

            utr3_start = list(trans_sub['3\' UTR start'])[0]
            utr3_stop = list(trans_sub['3\' UTR end'])[0]

            gene_name = list(trans_sub['Gene name'])[0]
            gene_len = list(trans_sub['Transcript length (including UTRs and CDS)'])[0]
            trans_len = list(trans_sub['Transcript length'])[0]

            intron_start = []
            intron_end = []

            for i in range(0, len(exon_rank)-1):

                if int(strand) == 1:

                    intron_start.append(exon_stop[i])
                    intron_end.append(exon_start[i+1])

                else:

                    intron_start.append(exon_stop[i+1])
                    intron_end.append(exon_start[i])

# '''----------------------------------------------------------------------------------------------------------------'''

            # section added to record information about introns and exons

            number_exon = len(exon_start)
            number_intron = len(intron_start)


            exon_mean_ls = []

            for ex in range(0, len(exon_start)):

                if ex > exon_est - 1:

                    print('exon estimate exceeded:', ex)

                    continue
                size = abs(exon_stop[ex] - exon_start[ex])

                exon_mean_ls.append(size)

            intron_mean_ls = []

            for intron in range(0, len(intron_start)):

                if intron > exon_est - 1:

                    print('intron estimate exceeded:', intron)

                    continue
                size = abs(intron_end[intron] - intron_start[intron])

                intron_mean_ls.append(size)


            exon_mean = sum(exon_mean_ls)/len(exon_mean_ls)
            exon_max = max(exon_mean_ls)

            if len(intron_mean_ls) == 0:

                intron_mean = 0
                intron_max = 0

            else:

                intron_mean = sum(intron_mean_ls)/len(intron_mean_ls)
                intron_max = max(intron_mean_ls)

            segment_max = max([intron_max, exon_max])

            start_loc = 0

            match = False

            for y in range(0, len(exon_ls)):
                # print(y, match)
                # print(len(intron_start))
                # if y < len(intron_start):
                    # print(intron_start[y], intron_end[y])
                if match:
                    # print('here1')


                    continue

                elif int(start) in range(int(exon_start[y]), int(exon_stop[y])):
                    # print('here4')
                    start_loc = (['exon', str(y+1)])
                    start_len = exon_mean_ls[y]

                    match = True

                elif y < len(intron_start):
                    # print('here5')
                    # print('check')
                    if int(start) in range(int(intron_start[y]), int(intron_end[y])):

                        start_loc = (['intron', str(y+1)])
                        start_len = intron_mean_ls[y]
                        # print('---------------------------------------------------------------------------------------')
                        # print(start_loc)
                        # print('---------------------------------------------------------------------------------------')
                        match = True

                    else:

                        continue

                else:
                    # print('here6')
                    start_loc = ['NA', 'NA']
                    start_len=0

                    continue

            stop_loc = 0
            match = False

            for y in range(0, len(exon_ls)):

                if match:

                    continue

                elif int(stop) in range(int(exon_start[y]), int(exon_stop[y])):

                    stop_loc = (['exon', str(y+1)])
                    stop_len = exon_mean_ls[y]

                    match = True

                elif y < len(intron_start):

                    if int(stop) in range(int(intron_start[y]), int(intron_end[y])):

                        stop_loc = (['intron', str(y+1)])
                        stop_len = intron_mean_ls[y]

                        match = True

                    else:

                        continue

                else:

                    stop_loc = ['NA', 'NA']
                    stop_len = 0

                for y in range(0, len(exon_start)):

                    if not numpy.isnan(utr3_start) and not numpy.isnan(utr5_start):

                        if stop_loc == 'NA' and int(utr3_start) <= int(stop) <= int(utr3_stop):

                            stop_loc = '3_utr'
                            tom = input('UTR!!!!!')

                    if not numpy.isnan(utr5_start) and not numpy.isnan(utr3_start):

                        if stop_loc == 'NA' and int(utr3_start) <= int(stop) <= int(utr3_stop):

                            stop_loc = '5_utr'

                            tom = input('UTR!!!!!')



                # if stop_loc == 'stop_no_match':
                #
                #     print('****************************************************************************************')
                #     print('\n\n')
                #     print(match,'\t', y)
                #     print(exon_start)
                #     print(exon_stop)
                #     print('\n')
                #     print(intron_start)
                #     print(intron_end)
                #     print(start)
                #
                #     input()


            if 'NA' in start_loc or 'NA' in stop_loc:

                full_match = 'N'

            else:

                full_match = 'Y'

            outfile.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(chr, mir, alias, start, stop, derives, gene_type, gene_len, gene, trans, trans_len, strand, '\t'.join(start_loc), '\t'.join(stop_loc), start_len, stop_len, list(trans_sub['Transcript start (bp)'])[0], list(trans_sub['Transcript end (bp)'])[0], full_match, number_exon, number_intron, exon_mean, intron_mean, exon_max, intron_max, segment_max))
            len_list.append(int(trans_len))
            len_pro_list.append(int(gene_len))
            type_list.append(gene_type)
            transcript_list.append(trans)
            gene_id_list.append(gene)


# ens_len_ls = list(int(x) for x in list(ensembl['Transcript length (including UTRs and CDS)']))
#
# print('trans: ', sum(ens_len_ls)/len(ens_len_ls), max(ens_len_ls), min(ens_len_ls))
#
print('host: ', sum(len_list)/len(len_list), max(len_list), min(len_list))

outfile.close()

#----------------------------------------------------------------------------------------------------------------------
# get transcript lengths for multi-exonic transcripts

type_short = list(set(type_list))
trans_list_short = list(set(transcript_list))
print(type_short)
#
ensembl_match = ensembl.loc[ensembl['Gene type'].isin(type_short)]
ensembl_match_no_primir = ensembl_match.loc[~ensembl_match['Transcript stable ID'].isin(trans_list_short)]
#
# match_len_ls = list(int(x) for x in list(ensembl_match['Transcript_length']))
# match_len_ls_mirless = list(int(x) for x in list(ensembl_match_no_primir['Transcript_length']))
#
#
# print('matched type: ', sum(match_len_ls)/len(match_len_ls))
# print('match_no_mir: ', sum(match_len_ls_mirless)/len(match_len_ls_mirless))

# take columns from ensembl then overwrite ensembl to reduce memory requirements

other_box = pd.DataFrame(ensembl_match_no_primir[['Gene stable ID', 'Transcript stable ID', 'Transcript length', 'Gene type', 'Transcript length (including UTRs and CDS)']])
other_box['host'] = 'N'
other_box.columns = ['gene_id', 'transcript_id', 'transcript_length', 'transcript_type', 'processed_length', 'host']
other_box['number_mirna'] = 'nan'
other_box['mirnas'] = 'nan'

ensembl_match_no_primir = False
ensembl_match = False

'''------------------------------------------------------------------------------------------------------------------'''

# count microRNAs per host

outframe = pd.DataFrame(pd.read_csv('{}/mirHOST/{}_hosts.txt'.format(dir, species), sep='\t', header=0))

num_mir = []
mirs = []

for trans in transcript_list:

    sub = outframe.loc[(outframe['transcript'] == trans) & (outframe['derives'] == 'self')]
    mirnas = list(set(list(sub['mir'])))
    num_mir.append(len(mirnas))
    mirs.append(' / '.join(mirnas))


# pd.DataFrame.to_csv(ensembl_match_no_primir, '{}/mirHOST/{}_match_trans.txt'.format(dir, species), index=False, sep='\t')

# making table for boxplots



gene_id_list = pd.DataFrame(gene_id_list, columns=['gene_id'])
transcript_list = pd.DataFrame(transcript_list, columns=['transcript_id'])
len_list = pd.DataFrame(len_list, columns=['transcript_length'])
len_pro_list = pd.DataFrame(len_pro_list, columns=['processed_length'])
type_list = pd.DataFrame(type_list, columns=['transcript_type'])
num_mir = pd.DataFrame(num_mir, columns=['number_mirna'])
mirs = pd.DataFrame(mirs, columns=['mirnas'])

mir_box = pd.DataFrame(pd.concat([gene_id_list, transcript_list, type_list, len_list, len_pro_list, num_mir, mirs], axis=1))
mir_box.drop_duplicates(inplace=True)
mir_box.index = list(range(0, len(mir_box)))
other_box.index = list(range(len(mir_box), (len(mir_box)+len(other_box))))
other_box.drop_duplicates(inplace=True)

mir_box['host'] = 'Y'

all_box = pd.DataFrame(pd.concat([mir_box, other_box], axis=0))

pd.DataFrame.to_csv(all_box, '{}/mirHOST/{}_boxdata.txt'.format(dir, species), sep='\t', index=False)

print('host mean: ', sum(mir_box['transcript_length'])/len(mir_box))
print('other mean: ', sum(other_box['transcript_length'])/len(other_box))

mir_host_ls = list(set(mir_box['transcript_id']))
not_host_ls = list(set(other_box['transcript_id']))

del all_box
del mir_box
del other_box

def get_strucutre(tran_list, host):

    for trans in tran_list:

        trans_sub = ensembl.loc[ensembl['Transcript stable ID'] == trans]

        gene = list(trans_sub['Gene stable ID'])[0]

        exon_ls = list(set(list(trans_sub['Exon stable ID'])))
        gene_type = list(trans_sub['Gene type'])[0]

        exon_rank = list(trans_sub['Exon rank in transcript'])
        exon_start = list(trans_sub['Exon region start (bp)'])
        exon_stop = list(trans_sub['Exon region end (bp)'])
        strand = list(trans_sub['Strand'])[0]
        # gen_source = list(trans_sub['Source of gene name'])[0]

        # print(exon_ls)
        # print(exon_rank)
        # print(exon_start)
        # print(exon_stop)
        # print('Strand: ', strand)

        utr5_start = list(trans_sub['5\' UTR start'])[0]
        utr5_stop = list(trans_sub['5\' UTR end'])[0]
        utr5_len = abs(utr5_stop-utr5_start)


        utr3_start = list(trans_sub['3\' UTR start'])[0]
        utr3_stop = list(trans_sub['3\' UTR end'])[0]
        utr3_len = abs(utr3_stop - utr3_start)

        gene_name = list(trans_sub['Gene name'])[0]
        gene_len = list(trans_sub['Transcript length (including UTRs and CDS)'])[0]
        trans_len = list(trans_sub['Transcript length'])[0]

        intron_start = []
        intron_end = []

        for i in range(0, len(exon_rank)-1):

            if int(strand) == 1:

                intron_start.append(exon_stop[i])
                intron_end.append(exon_start[i+1])

            else:

                intron_start.append(exon_stop[i+1])
                intron_end.append(exon_start[i])


        ex_size_ls = []

        for ex in range(0, len(exon_start)):

            size = abs(exon_stop[ex] - exon_start[ex])

            ex_size_ls.append(size)

        in_size_ls = []

        for intron in range(0, len(intron_start)):

            size = abs(intron_end[intron] - intron_start[intron])

            in_size_ls.append(size)

        if len(ex_size_ls) != 0:

            ex_max = max(ex_size_ls)

        else:

            ex_max = 0

        if len(in_size_ls) != 0:

            in_max = max(in_size_ls)

        else:

            in_max = 0

        for ex in range(0, len(exon_start)):

            size = ex_size_ls[ex]

            structure_file.write(
                '{}\t{}\texon\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(gene, trans, ex + 1, len(exon_start),
                                                                                len(in_size_ls), in_max, ex_max
                                                                                , utr3_len, utr5_len,
                                                                                gene_type, size, host))

        for intron in range(0, len(intron_start)):

            size = in_size_ls[intron]

            structure_file.write(
                '{}\t{}\tintron\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(gene, trans, intron + 1,
                                                                                  len(ex_size_ls), len(intron_start),
                                                                                  in_max,
                                                                                  ex_max, utr3_len, utr5_len,
                                                                                  gene_type, size, host))

print('Getting host structure')
get_strucutre(mir_host_ls, 'Y')
print('Getting other structures')
get_strucutre(not_host_ls, 'N')