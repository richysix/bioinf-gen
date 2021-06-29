#!/usr/bin/env python3
''' Script to take haplotype input in region form and collapse adjacent regions
with the same parental haplotypes together '''
import argparse
import sys
from random import randrange

def output_collapsed_region(regions, offspring_cols, args):
    if args.debug:
        print(regions)
    chrom = regions[0][args.chr_col]
    start = regions[0][args.start_col]
    end = regions[len(regions) - 1][args.end_col]
    parent_1_ht = regions[0][args.parent1_col]
    parent_2_ht = regions[0][args.parent2_col]
    
    offspring_haplotypes = assign_offspring_haplotypes(regions, offspring_cols, args)
    output_list = [chrom, start, end, parent_1_ht, parent_2_ht]
    output_list.extend(offspring_haplotypes)
    print("\t".join(output_list))

def assign_offspring_haplotypes(regions, offspring_cols, args):
    offspring_calls = [ {} for x in range(len(offspring_cols)) ]
    for region in regions:
        if args.debug:
            print(region)
        for idx, col_num in enumerate(offspring_cols):
            if args.debug:
                print(col_num, end = "\t")
            if region[col_num] not in offspring_calls[idx]:
                offspring_calls[idx][region[col_num]] = 1
            else:
                offspring_calls[idx][region[col_num]] = offspring_calls[idx][region[col_num]] + 1
    
    offspring = []
    for i in range(len(offspring_calls)):
        # turn counts for haplotypes into haplotypes for counts
        ht_counts = {}
        for ht in offspring_calls[i].keys():
            count = offspring_calls[i][ht]
            if count in ht_counts:
                ht_counts[count].append(ht)
            else:
                ht_counts[count] = [ht]
        # check if there is a haplotype with a majority
        if len(ht_counts[max(ht_counts)]) == 1:
            offspring.append(ht_counts[max(ht_counts)][0])
        else:
            hts = ht_counts[max(ht_counts)]
            offspring.append(hts[randrange(0, len(hts), 1)])
    return(offspring)


def main(args):
    ''' Main body of code '''
    # open file and read in input
    processed_header = False
    current_regions = []
    last_region_idx = 0
    # remove chr, start, end, parent1 and parent2 columns to get the offspring columns
    exclude_cols = [ args.chr_col, args.start_col, args.end_col, args.parent1_col, args.parent2_col ]
    
    for line in args.input_file:
        line = line.rstrip()
        cols = line.split("\t")
        if not processed_header:
            processed_header = True
            # get list of offspring column indexes
            offspring_cols = []
            for i in range(len(cols)):
                if i not in exclude_cols:
                    offspring_cols.append(i)
            if args.header:
                print(line)
                continue
        # if current regions is empty, add region
        if len(current_regions) == 0:
            current_regions.append(cols)
        else:
            # check whether the chromosomes are the same
            if current_regions[last_region_idx][args.chr_col] != cols[args.chr_col]:
                output_collapsed_region(current_regions, offspring_cols, args)
                current_regions = [ cols ]
                last_region_idx = 0
            # check whether the regions are adjacent
            elif int(current_regions[last_region_idx][args.end_col]) + 1 == int(cols[args.start_col]):
                # check whether haplotypes are the same
                p1_ht_1 = current_regions[last_region_idx][args.parent1_col]
                p1_ht_2 = cols[args.parent1_col]
                p2_ht_1 = current_regions[last_region_idx][args.parent2_col]
                p2_ht_2 = cols[args.parent2_col]
                if p1_ht_1 == p1_ht_2 and p2_ht_1 == p2_ht_2:
                    current_regions.append(cols)
                    last_region_idx += 1
                else:
                    # if not output last region
                    output_collapsed_region(current_regions, offspring_cols, args)
                    current_regions = [ cols ]
                    last_region_idx = 0
            else:
                # if not output last region
                output_collapsed_region(current_regions, offspring_cols, args)
                current_regions = [ cols ]
                last_region_idx = 0
    # output final region
    output_collapsed_region(current_regions, offspring_cols, args)

    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('input_file', nargs='?', metavar='FILE',
        type=argparse.FileType('r'), default=sys.stdin, help='')
    parser.add_argument('--header', action='store_true', default=False,
        help='Input has header line')
    parser.add_argument('--parent1_col', action='store', default=3,
        help='Column number for parent 1 column (zero based)')
    parser.add_argument('--parent2_col', action='store', default=4,
        help='Column number for parent 2 column (zero based)')
    parser.add_argument('--chr_col', action='store', default=0,
        help='Column number for chr column (zero based)')
    parser.add_argument('--start_col', action='store', default=1,
        help='Column number for start column (zero based)')
    parser.add_argument('--end_col', action='store', default=2,
        help='Column number for end column (zero based)')
    parser.add_argument('--debug', action='count', default=0,
        help='Prints debugging information')
    args = parser.parse_args()
    main(args)
    