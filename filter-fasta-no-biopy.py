#! /usr/bin/env python

''' Script to filter fasta file based on header lines '''

import argparse
import sys

def output_sequence(fh, id, seqLines):
    fh.write(">" + id + "\n")
    fh.writelines(seqLines)

def filter(current_id, current_seq, filter_strings, output_handle, args):
    if current_id in filter_strings:
        if args.verbose:
            print(f"Current sequence {current_id}, matches one of the filter strings")
        if args.filter_type == "keep":
            if args.verbose:
                print(f"Keeping...")
            output_sequence(output_handle, current_id, current_seq)
        else:
            if args.verbose:
                print(f"Removing...")
    else:
        if args.verbose:
            print(f"Current sequence {current_id}, does not match one of the filter strings")
        if args.filter_type == "remove":
            if args.verbose:
                print(f"Keeping...")
            output_sequence(output_handle, current_id, current_seq)
        else:
            if args.verbose:
                print(f"Removing...")

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

    output_handle = open(args.output_file, "w", newline='\n')

    current_id = None
    current_seq = []
    with open(args.input_file, 'r') as seq_f:
        for line in seq_f:
            if args.debug:
                print(line)
            if line[0] == ">":
                id = line.replace(">", "")
                id = id.split()[0]
                if args.debug:
                    print('New sequence')
                # new sequence
                if current_id is not None:
                    filter(current_id, current_seq, filter_strings, output_handle, args)
                    current_id = id
                    current_seq = []
                else:
                    current_id = id
            else:
                current_seq.append(line)
        filter(current_id, current_seq, filter_strings, output_handle, args)

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
    parser.add_argument('--verbose', action='store_true', default=False,
        help='Prints some information about which sequences have been kept'),
    parser.add_argument('--debug', action='count', default=0,
        help='Prints debugging information')
    params = parser.parse_args()
    main(params)
