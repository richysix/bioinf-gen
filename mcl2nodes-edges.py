#!/usr/bin/python3

desc = ''' Script to parse an mcl graph to nodes and edges file 

The gene info file should contain the following columns

node_index - index for node corresponding to current network
id - index covering all genes
gene_id - Ensembl gene ID
gene_name - Gene name
description - Description
'''
import argparse
import re
import sys

def main(args):
    ''' Main body of code '''
    # load cluster file first
    cluster_idx_for = {}
    begin = False
    nodelist = []
    with open(args.cluster_file, 'r') as f:
        for line in f:
            info = line.lstrip().rstrip('\n').split()
            if not begin:
                if info[0] == "begin":
                    begin = True
                    continue
            else:
                if info[-1] == "$":
                    # remove $
                    info.pop()
                    complete = True
                else:
                    complete = False
                for x in info:
                    nodelist.append(x)
                if complete:
                    if args.debug > 1:
                        print(nodelist)
                    cluster_idx = None
                    for node_idx in nodelist:
                        if cluster_idx is None:
                            cluster_idx = node_idx
                        else:
                            cluster_idx_for[node_idx] = int(cluster_idx)
                    # reset node list
                    nodelist = []
    if args.debug > 1:
        print(cluster_idx_for)
    
    # load tab file
    with open(args.nodes_file, 'w') as out:
        print('node_idx', 'id', 'gene_id', 'gene_name', 'cluster_id', sep=',', file=out)
        with open(args.gene_info_file, 'r') as f:
            for line in f:
                info = line.rstrip('\n').split('\t') # columns are node_index, id, gene_id, gene_name, description
                node_idx = info[0]
                cluster_idx = cluster_idx_for[node_idx] + 1
                print(str.format('"{0}","{1}","{2}","{3}","{4}"', node_idx, info[1], info[2], info[3], str(cluster_idx)), file=out)
    
    # open file and read in input
    begin = False
    edgelist = []
    edge_idx = 1
    edges = {}
    with open(args.edges_file, 'w') as out:
        # add header
        print('edge_idx', 'source', 'target', 'weight', sep=",", file=out)
        with open(args.mci_file, 'r') as f:
            for line in f:
                info = line.lstrip().rstrip('\n').split()
                if not begin:
                    if info[0] == "begin":
                        begin = True
                        continue
                else:
                    if info[-1] == "$":
                        # remove $
                        info.pop()
                        complete = True
                    else:
                        complete = False
                    for x in info:
                        edgelist.append(x)
                    if complete:
                        # go through list and print edges
                        source_node_idx = None
                        for edge in edgelist:
                            # first idx is sorce node
                            if source_node_idx is None:
                                source_node_idx = re.sub(':[0-9]+', '', edge)
                                continue
                            else:
                                (target_node_idx, weight) = edge.split(":")
                                # remove self edges
                                if source_node_idx == target_node_idx:
                                    continue
                                # check whether this edge already exists in opposite direction
                                edge_string = target_node_idx + '-' + source_node_idx
                                if edge_string in edges:
                                    continue
                                else:
                                    edges[source_node_idx + '-' + target_node_idx] = 1
                            print(edge_idx, source_node_idx, target_node_idx, weight, sep=",", file=out)
                            edge_idx += 1
                        # clear edge list
                        edgelist = []
                        
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('mci_file', metavar='MCIFILE',
        type=str, default='graph.mci', help='mcl graph file')
    parser.add_argument('cluster_file', metavar='CLUSTERFILE',
        type=str, default='graph.mci.I14', help='cluster file from mcl showing which nodes belong to each cluster')
    parser.add_argument('gene_info_file', metavar='GENEFILE',
        type=str, default='graph.tab', help='gene info file containing mcl node index and gene id and name')
    parser.add_argument('--nodes_file', metavar='NODESFILE',
        type=str, default='graph.nodes.csv', help='Name of the output file for the nodes')
    parser.add_argument('--edges_file', metavar='EDGESFILE',
        type=str, default='graph.edges.csv', help='Name of the output file for the edges')
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
