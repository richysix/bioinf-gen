#! /usr/bin/env python

''' Script to filter fasta file based on header lines '''

import argparse
import sys
from Bio import SeqIO

def main(args):
    ''' Main body of code '''

    # read in input
    filter_strings = []
    with open(args.filter_file, 'r') as f:
        for line in f:
            line = line.rstrip('\n')
            filter_strings.append(line)
    if args.debug > 0:
        print(filter_strings)

    output_handle = open(args.output_file, "w")

    for record in SeqIO.parse(args.input_file, "fasta"):
        if record.id in filter_strings:
            if args.filter_type == "keep":
                SeqIO.write(record, output_handle, "fasta")
        else:
            if args.filter_type == "remove":
                SeqIO.write(record, output_handle, "fasta")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('input_file', nargs='?', metavar='FILE',
        type=str, default='in.fa', help='Name of the fast file to filter [%(default)s]')
    parser.add_argument('--output_file', metavar='FILE',
        type=str, default="out.fa", help='Name of the output file [%(default)s]')
    parser.add_argument('--filter_file', metavar='FILE',
        type=str, default="ids.txt", help='File of strings to use as the filtering criteria [%(default)s]')
    parser.add_argument('--filter_type', metavar='FILE',
        type=str, choices = ["keep", "remove"], default="remove", 
        help='Whether to keep or remove the sequences based on the filter file [%(default)s]')
    parser.add_argument('--debug', action='count', default=0,
        help='Prints debugging information')
    params = parser.parse_args()
    main(params)
