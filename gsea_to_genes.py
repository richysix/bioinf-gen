#!/usr/bin/env python3
''' Script to process a GSEA output file and get the genes that drive the
enrichment'''
import argparse
import sys
import os

def main(args):
    ''' Main body of code '''
    
    # open genes file if it is specified
    genes = {}
    if args.genes_file is not None:
        for line in args.genes_file:
            genes[line.rstrip()] = 1
    
    # constants for indices of columns
    NAME_IDX=0
    NES_IDX=5
    FDR_IDX=7
    FDR_THRESHOLD=args.fdr_threshold
    # open file and read in input
    header = True
    for line in args.input_file:
        if args.debug > 1:
            print(line, sep = "\t")
        if header:
            header = False
            continue
        items = line.rstrip().split("\t")
        if float(items[FDR_IDX]) < FDR_THRESHOLD:
            name = items[NAME_IDX]
            detail_filename = os.path.join(args.base_dir, name + '.xls')
            # open detail file
            # constants for detail file
            GENE_SYMBOL_IDX=2
            GENE_TITLE_IDX=3
            CORE_ENRICHMENT_IDX=7
            detail_header = True
            with open(detail_filename) as detail_file:
                for gene_line in detail_file:
                    if args.debug > 1:
                        print(gene_line, sep = "\t")
                    if detail_header:
                        detail_header = False
                        continue
                    info = gene_line.rstrip().split("\t")
                    if info[CORE_ENRICHMENT_IDX] == "Yes":
                        output = False
                        if args.genes_file is None:
                            output = True
                        else:
                            if info[GENE_SYMBOL_IDX] in genes:
                                output = True
                        if output:
                            gene_info = info[GENE_TITLE_IDX].replace(" | ", "\t").rstrip(" |")
                            print(args.comparison, name, info[GENE_SYMBOL_IDX], gene_info,
                                  sep = "\t", file = args.output_file)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('input_file', nargs='?', metavar='FILE',
        type=argparse.FileType('r'), default=sys.stdin, help='Input file (GSEA .xls file) [default: STDIN]')
    parser.add_argument('output_file', nargs='?', metavar='FILE',
        type=argparse.FileType('w'), default=sys.stdout, help='Output file [default: STDOUT]')
    parser.add_argument('--genes_file', nargs='?', metavar='FILE',
        type=argparse.FileType('r'), default=None, help='File of genes to limit output to')
    parser.add_argument('--base_dir', metavar='Str',
        type=str, default=os.getcwd(), help='Base directory')
    parser.add_argument('--comparison', metavar='Str',
        type=str, default='mut_vs_wt', help='Name of the comparison')
    parser.add_argument('--fdr_threshold', type=float, metavar='FLOAT',
        default=0.05, help='FDR threshold for significant terms')
    parser.add_argument('--debug', action='count', default=0,
        help='Prints debugging information')
    args = parser.parse_args()
    main(args)

# AUTHOR
#
# Richard White <rich@buschlab.org>
#
# COPYRIGHT AND LICENSE
#
# This software is Copyright (c) 2020. University of Cambridge.
#
# This is free software, licensed under:
#
#  The GNU General Public License, Version 3, June 2007
