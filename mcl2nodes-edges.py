#!/usr/bin/python3

desc = """Script to parse an mcl graph to nodes and edges file

The gene info file should contain the following columns

node_index - index for node corresponding to current network
id - index covering all genes
gene_id - Ensembl gene ID
gene_name - Gene name
description - Description
"""

import argparse
import re
import sys

def main(args):
    ''' Main body of code '''
    # load cluster file first
    cluster_idx_for = {}
    nodes_for_cluster = {}
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
                nodelist.extend(info)
                if complete:
                    if args.debug > 1:
                        print(nodelist)
                    cluster_idx = None
                    for node_idx in nodelist:
                        if cluster_idx is None:
                            cluster_idx = node_idx
                        else:
                            cluster_idx_for[node_idx] = cluster_idx
                            if cluster_idx not in nodes_for_cluster:
                                nodes_for_cluster[cluster_idx] = [node_idx]
                            else:
                                nodes_for_cluster[cluster_idx].append(node_idx)
                    # reset node list
                    nodelist = []
    if args.debug > 1:
        print(cluster_idx_for, file = sys.stderr)
        print(nodes_for_cluster, file = sys.stderr)

    # open tab file
    gene_id_for = {}
    with open(args.graph_tab_file, 'r') as f:
        for line in f:
            info = line.rstrip('\n').split('\t') # columns are node_index, gene_id
            gene_id_for[info[0]] = info[1]
    gene_info_for = {}
    with open(args.annotation_file, 'r') as f:
        for idx, line in enumerate(f):
            info = line.rstrip('\n').split('\t') # columns are gene_id, chr, start, end, strand, biotype, name, description
            gene_info_for[info[0]] = { 
                "name": info[6], 
                "idx": idx,
                "description": info[7]
            }

    # open nodes output file
    with open(args.nodes_file, 'w') as out:
        print('node_idx', 'id', 'gene_id', 'gene_name', 'description', 'cluster_id', 'singleton', sep=',', file=out)
        for node_idx, gene_id in gene_id_for.items():
            cluster_idx = cluster_idx_for[node_idx]
            singleton = True if len(nodes_for_cluster[cluster_idx]) == 1 else False
            cluster_idx = int(cluster_idx) + 1
            print(str.format('"{0}","{1}","{2}","{3}","{4}","{5}","{6}"', 
                node_idx, gene_info_for[gene_id]["idx"], gene_id, 
                gene_info_for[gene_id]["name"], 
                gene_info_for[gene_id]["description"],
                str(cluster_idx), singleton), file=out)
    
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
                            if args.edge_offset:
                                weight = float(weight) + args.edge_offset
                                weight = f"{weight:.3f}"
                            print(edge_idx, source_node_idx, target_node_idx, weight, sep=",", file=out)
                            edge_idx += 1
                        # clear edge list
                        edgelist = []
                        
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=desc,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('mci_file', metavar='MCI_FILE',
        type=str, default='graph.mci', help='mcl graph file')
    parser.add_argument('cluster_file', metavar='CLUSTER_FILE',
        type=str, default='graph.mci.I14', help='cluster file from mcl showing which nodes belong to each cluster')
    parser.add_argument('graph_tab_file', metavar='GRAPH_TAB_FILE',
        type=str, default='graph.tab', help='graph tab file containing mcl node index and gene id')
    parser.add_argument('annotation_file', metavar='ANNOTATION_FILE',
        type=str, default='annotation.txt', help='annotation file containing gene id and name')
    parser.add_argument('--nodes_file', metavar='NODES FILE',
        type=str, default='graph.nodes.csv', help='Name of the output file for the nodes')
    parser.add_argument('--edges_file', metavar='EDGES FILE',
        type=str, default='graph.edges.csv', help='Name of the output file for the edges')
    parser.add_argument('--edge_offset', metavar='EDGE OFFSET',
        type=float, default=None, help='An amount to add to edge weights')
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
