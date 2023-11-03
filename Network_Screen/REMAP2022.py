def resturcture_file(bedFile, outfile):

    o_bed = open(bedFile)
    o_out = open(outfile, 'w')

    for line in o_bed:

        split_line = line.strip().split('\t')

        peak_info = split_line[3].split('.')

        new_line = '\t'.join(split_line[0:3]+[peak_info[1].upper()]+split_line[4:]+[peak_info[0], peak_info[2]])+'\n'

        # print(new_line)

        o_out.write(new_line)

    o_out.close()
    o_bed.close()

if __name__ == '__main__':

    import argparse

    help_text = 'REMAP2018 reformats the remap bedfile into a format which will work with our mapping tools' \
                '\n\nThe input file and only argument should be the REMAP 2018 single peaks file'

    parser = argparse.ArgumentParser(description=help_text)

    # parser.add_help()

    parser.add_argument('remapBed')

    args = parser.parse_args()

    if '.bed' not in args.remapBed:     # Checks that the input file is a bed file, prevents overwriting of data file

        exit('Unexpected file extension: {}'
             '\n\nExpected .bed'.format(args.remapBed.split('.')[-1]))

    outfile = args.remapBed.replace('.bed', '-REFORMATED.bed')

    resturcture_file(bedFile=args.remapBed, outfile=outfile)

    exit()