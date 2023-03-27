#! /usr/bin/env python

''' Extract information of Grid Engine jobs 
from the output of the qacct command'''

import argparse
import sys
import re

def main(args):
    ''' Main body of code '''
    
    # read in input
    job_info = {}
    job_id = ''
    for line in args.input_file:
        if re.match(r'^========', line):
            if 'current_job' in job_info.keys():
                if job_info['current_job']['taskid'] == 'undefined':
                    job_id = job_info['current_job']['taskid']
                else:
                    job_id = job_info['current_job']['jobnumber'] + '.' + job_info['current_job']['taskid']
                job_info[job_id] = job_info['current_job']
            job_info['current_job'] = {}
            continue
        if args.debug:
            print(line)
        items = line.rstrip('\n').split()
        k = items[0]
        if len(items) == 2:
            v = items[1]
        else:
            v = " ".join(items[1:])
        job_info['current_job'][k] = v
    
    print("\t".join(args.properties))
    for job in job_info.values():
        items = list(map(lambda k: job[k], args.properties))
        print("\t".join(items))



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('input_file', nargs='?', metavar='FILE',
        type=argparse.FileType('r'), default=sys.stdin, help='')
    parser.add_argument('--properties', action='store', nargs = '*',
        default=['jobname', 'jobnumber', 'taskid', 'exit_status'],
        help='Set job properties to output')
    parser.add_argument('--debug', action='count', default=0,
        help='Prints debugging information')
    args = parser.parse_args()
    if args.debug:
        print(args)
    
    main(args)
