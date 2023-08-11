#!/usr/bin/env perl

# PODNAME: parse_mcl_output.pl
# ABSTRACT: Script to take the output of MCL clustering and produce an easier
# to read file of which nodes are in which cluster.
# It can also produce graph files of different formats

use warnings;
use strict;
use Getopt::Long;
use autodie;
use Pod::Usage;
use Scalar::Util qw{looks_like_number};
use List::Util qw(all);

# get options
my %options;
get_and_check_options();

if ($options{'debug'}) {
    use Data::Dumper;
}

# open tab file
open my $tab_fh, '<', $ARGV[0];
my %node_idx_for = (); # node_idx_for = { NODE_ID => NODE_IDX }
my %node_id_for = (); # node_id_for = { NODE_IDX => NODE_ID }
while(my $line = <$tab_fh>) {
    chomp $line;
    my ($node_idx, $id) = split /\s+/, $line;
    $node_idx_for{$id} = $node_idx;
    $node_id_for{$node_idx} = $id;
}
if ($options{'debug'} > 1) {
    warn Dumper(%node_idx_for);
}

# open edges file
open my $edge_fh, '<', $ARGV[1];
my %weight_for = (); # %weight_for = { 'SOURCE-TARGET' => WEIGHT }
# In this hash SOURCE is always less than TARGET to avoid duplicating edges
my $begin = 0;
my $edge_list = [];
while(my $line = <$edge_fh>) {
    chomp $line;
    if (!$begin) {
        if ($line eq "begin") {
            $begin = 1;
            next;
        }
    } else {
        if ($line eq ")") {
            last;
        } else {
            # read in cluster numbers
            $line =~ s/\A \s+//xms;
            my @info = split /\s+/, $line;
            my $complete = 0;
            if ($info[ (scalar @info) - 1 ] eq '$') {
                $complete = 1;
                splice(@info, -1);
            }
            push @{$edge_list}, @info;
            if ($complete) {
                if ($options{'debug'} > 1) {
                    warn Dumper($edge_list);
                }
                my $source_node_idx = splice(@{$edge_list}, 0, 1);
                foreach my $edge (@{$edge_list}) {
                    my ($target_node_idx, $weight) = split /:/, $edge;
                    my $edge_name = $source_node_idx < $target_node_idx ?
                        join('-', $source_node_idx, $target_node_idx) :
                        join('-', $target_node_idx, $source_node_idx);
                    $weight_for{$edge_name} = $weight;
                }
                # reset node list
                $edge_list = [];
            }
        }
    }
}
if ($options{'debug'} > 1) {
    warn Dumper(%weight_for);
}

# open cluster file
open my $cluster_fh, '<', $ARGV[2];
$begin = 0;
my $node_list = [];
my %cluster_idx_for; # %cluster_idx_for = { NODE_IDX => CLUSTER_IDX }
my %nodes_for = ();
while(my $line = <$cluster_fh>) {
    chomp $line;
    if (!$begin) {
        if ($line eq "begin") {
            $begin = 1;
            next;
        }
    } else {
        if ($line eq ")") {
            last;
        } else {
            # read in cluster numbers
            $line =~ s/\A \s+//xms;
            my @info = split /\s+/, $line;
            my $complete = 0;
            if ($info[ (scalar @info) - 1 ] eq '$') {
                $complete = 1;
                splice(@info, -1);
            }
            push @{$node_list}, @info;
            if ($complete) {
                my $cluster_idx = undef;
                foreach my $node_idx (@{$node_list}) {
                    if (!defined $cluster_idx) {
                        $cluster_idx = $node_idx;
                    } else {
                        $cluster_idx_for{$node_idx} = $cluster_idx;
                        push @{$nodes_for{$cluster_idx}}, $node_idx;
                    }
                }
                # reset node list
                $node_list = [];
            }
        }
    }
}
if ($options{'debug'} > 1) {
    warn Dumper(%cluster_idx_for);
}

# open info file
my @col_names = ();
my @output_header = ('Node_idx', 'Node_id', 'Cluster_idx');;
my $node_id_col_name = 'Node_id';
my $node_name_col_name;
my @other_cols = ();
my %info_for_idx = (); # $info_for_idx{ IDX => { ID => 'ENSDARG01', NAME => 'sox10', ATTR1 => 'cluster1', ATTR2 => 'square' } }
my %type_for; # %type_for{ COL_NAME1 => [ 'string', 'string' ], COL_NAME1 => [ 'double', 'double ] }
if ($options{'info_file'}){
    open my $info_fh, '<', $options{'info_file'};
    # deal with header
    my %col_id_for = ();
    my $line = <$info_fh>;
    chomp $line;
    @col_names = split /\t/, $line;
    for (my $i = 0; $i < scalar @col_names; $i++) {
        $col_id_for{$col_names[$i]} = $i;
        # if this column is the node id one then change $node_id_col_name
        if ($i == $options{'info_file_node_id_col'}) {
            $node_id_col_name = $col_names[$i];
        } else {
            # check whether a node name column number is defined and
            # whether this column is it
            if ($options{'info_file_node_name_col'} &&
                $i == $options{'info_file_node_name_col'}) {
                $node_name_col_name = $col_names[$i];
            } else {
                push @other_cols, $col_names[$i];
            }
        }
    }
    # replace 'Node_id' in header with node id column name
    $output_header[1] = $node_id_col_name;
    if ($node_name_col_name) {
        push @output_header, $node_name_col_name;
    }
    
    while(my $line = <$info_fh>) {
        chomp $line;
        my @fields = split /\t/, $line;
        my $node_id = $fields[ $options{'info_file_node_id_col'} ];
        if (!exists $node_idx_for{$node_id}) {
            warn "Node with id $node_id doesn't exist in the MCL tab file!\n";
            next;
        }
        my $node_idx = $node_idx_for{$node_id};
        for(my $i = 0; $i < scalar @fields; $i++) {
            my $field_name = $col_names[$i];
            push @{$type_for{$field_name}}, looks_like_number($fields[$i]) ? 'double' : 'string';
            $info_for_idx{$node_idx}{$field_name} = $fields[$i];
        }
    }
}
foreach my $col ( sort keys %type_for ) {
    if (all { $_ eq 'double' } @{$type_for{$col}} ) {
        $type_for{$col} = 'double';
    } else {
        $type_for{$col} = 'string'
    }
}
if ($options{'debug'} > 1) {
    warn Dumper(%type_for);
}

# open cluster output file
open my $output_fh, '>', $ARGV[3];
# print header line
print {$output_fh} join("\t", @output_header, @other_cols, ), "\n";

my $graphml_fh;
my %graph_attr_num_for = ();
if ($options{'graphml_file'}) {
    # open cluster file
    open $graphml_fh, '>', $options{'graphml_file'};
    # print header
    my $graph_id = $options{'graph_id'} ? $options{'graph_id'} : 'G';
    my $graph_header = <<"END_HEADER";
<?xml version="1.0" encoding="UTF-8"?>
<graphml xmlns="http://graphml.graphdrawing.org/xmlns"  
    xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"
    xsi:schemaLocation="http://graphml.graphdrawing.org/xmlns 
     http://graphml.graphdrawing.org/xmlns/1.0/graphml.xsd">
  <key id="d0" for="edge" attr.name="weight" attr.type="double"/>
  <key id="d1" for="node" attr.name="node_id" attr.type="string"/>
  <key id="d2" for="node" attr.name="cluster_idx" attr.type="int"/>
END_HEADER

    print {$graphml_fh} $graph_header;
    my $key_num = 3;
    if ($node_name_col_name) {
        $graph_attr_num_for{$node_name_col_name} = $key_num;
        print {$graphml_fh} qq{  <key id="d${key_num}" for="node" attr.name="${node_name_col_name}" attr.type="$type_for{$node_name_col_name}"/>"}, "\n";
        $key_num++;
    }
    foreach my $attr ( @other_cols ) {
        $graph_attr_num_for{$attr} = $key_num;
        print {$graphml_fh} qq{  <key id="d${key_num}" for="node" attr.name="${attr}" attr.type="$type_for{$attr}"/>"}, "\n";
        $key_num++;
    }
    
    print {$graphml_fh} qq{  <graph id="$graph_id" edgedefault="undirected">}, "\n";
}

#my $csv_graph_node_fh;
#if ($options{'csv_graph_node_file'}) {
#    open $csv_graph_node_fh, '>', $options{'csv_graph_node_file'};
#    
#}

foreach my $cluster_idx ( sort { $a <=> $b } keys %nodes_for ) {
    foreach my $node_idx ( @{$nodes_for{$cluster_idx}} ) {
        my $node_id = $node_id_for{$node_idx};
    
        if ($graphml_fh){
            print {$graphml_fh} qq{    <node id="n${node_idx}">}, "\n";
            print {$graphml_fh} qq{      <data key="d1">${node_id}</data>}, "\n";
            print {$graphml_fh} qq{      <data key="d2">${cluster_idx}</data>}, "\n";
        }
        
        my @output_fields = ($node_idx, $node_id, $cluster_idx);
        if ($info_for_idx{$node_idx}) {
            if ($node_name_col_name){
                my $attr_value = $info_for_idx{$node_idx}{$node_name_col_name};
                push @output_fields, $attr_value;
                if ($graphml_fh){
                    my $key_num = $graph_attr_num_for{$node_name_col_name};
                    print {$graphml_fh} qq{      <data key="d${key_num}">${attr_value}</data>}, "\n";
                }
            }
            foreach my $key ( @other_cols ) {
                my $attr_value = $info_for_idx{$node_idx}{$key};
                push @output_fields, $attr_value;
                if ($graphml_fh){
                    my $key_num = $graph_attr_num_for{$key};
                    print {$graphml_fh} qq{      <data key="d${key_num}">${attr_value}</data>}, "\n";
                }
            }
        }
        print {$output_fh} join("\t", @output_fields, ), "\n";
        if ($graphml_fh){
            print {$graphml_fh} q{    </node>}, "\n";
        }
        
    }
}

if ($graphml_fh){
    my $edge_idx = 0;
    foreach my $edge ( sort keys %weight_for ) {
        my ($source, $target) = split /-/, $edge;
        my $weight = $weight_for{$edge};
        print {$graphml_fh} qq{    <edge id="e${edge_idx}" source="n${source}" target="n${target}">}, "\n";
        print {$graphml_fh} qq{      <data key="d0">${weight}</data>}, "\n";
        print {$graphml_fh} qq{    </edge>}, "\n";
        $edge_idx++;
    }
    print {$graphml_fh} "  </graph>\n</graphml>\n";
}

################################################################################
# SUBROUTINES

#subroutine_name
#
#  Usage       : subroutine_name( arguments )
#  Purpose     : 
#  Returns     : 
#  Parameters  : 
#  Throws      : 
#  Comments    : None

#get_and_check_options
#
#  Usage       : get_and_check_options()
#  Purpose     : parse the options supplied to the script using GetOpt::Long
#  Returns     : None
#  Parameters  : None
#  Throws      : 
#  Comments    : The default option are
#                help which print a SYNOPSIS
#                man which prints the full man page
#                debug
#                verbose

sub get_and_check_options {
    
    GetOptions(
        \%options,
        'info_file=s',
        'info_file_node_id_col=i',
        'info_file_node_name_col=i',
        'graphml_file=s',
        'graph_id=s',
        'csv_graph_node_file=s',
        'csv_graph_node_file=s',
        'no_info_file_header',
        'help',
        'man',
        'debug+',
        'verbose',
    ) or pod2usage(2);
    
    # Documentation
    if( $options{help} ) {
        pod2usage( -verbose => 0, -exitval => 1, );
    }
    elsif( $options{man} ) {
        pod2usage( -verbose => 2 );
    }
    
    # defaults
    $options{'debug'} = $options{'debug'} ? $options{'debug'} : 0;
    $options{'info_file_node_id_col'} = $options{'info_file_node_id_col'} ? $options{'info_file_node_id_col'} - 1: 0;
    $options{'info_file_node_name_col'} = $options{'info_file_node_name_col'} ? $options{'info_file_node_name_col'} - 1: $options{'info_file_node_name_col'};
    if ($options{'csv_graph_node_file'} && !$options{'csv_graph_edge_file'}) {
        if ($options{'csv_graph_node_file'} =~ m/node/xms) {
            $options{'csv_graph_edge_file'} = $options{'csv_graph_node_file'};
            $options{'csv_graph_edge_file'} =~ s/node/edge/xms;
        } else {
            die "csv node file name is set and edge file name is not!"
        }
    }
    
    print "Settings:\n", map { join(' - ', $_, defined $options{$_} ? $options{$_} : 'off'),"\n" } sort keys %options if $options{verbose};
}

__END__

=pod

=head1 NAME

parse_mcl_output.pl

=head1 DESCRIPTION

Script to take the output of MCL clustering and produce an easier to read file
of which nodes are in which cluster. It can also produce graph files of
different formats.

=cut

=head1 SYNOPSIS

    parse_mcl_output.pl [options] tab_file edge_file cluster_file output_file
        --info_file                 Name of a file containing metadata about the graph nodes (Default: undef)
        --info_file_node_id_col     Column in the info file for the node id [1-based]
                                    (Default: 1)
        --info_file_node_name_col   Column in the info file for the node name [1-based]
                                    (Default: undefined)
        --graphml_file              Name of a file to output the graph to in graphml format
        --graph_id                  An id for the graph in graphml file
                                    (Default: G)
        --csv_graph_node_file       Name of a file to output the nodes to in csv format
                                    (Default: undefined)
        --csv_graph_edge_file       Name of a file to output the edges to in csv format
                                    (Default: undefined)
        --no_info_file_header       flag for if the info file doesn't have a header
        --help                      print this help message
        --man                       print the manual page
        --debug                     print debugging information
        --verbose                   turn on verbose output


=head1 ARGUMENTS

=over

=item I<tab_file>

MCL file of node indices to node ids

=item I<edge_file>

MCL file of graph edges (.mci) with weights

=item I<cluster_file>

MCL cluster file showing which cluster membership for each node

=item I<output_file>

Name of the output file. This is tab-separated with columns for
Node Index, Node ID, Cluster Index, and any other columns in the metadata file.

=back

=head1 OPTIONS

=over

=item B<--info_file>

Optional tab-separated file with more metadata for each node.
Must contain a column with unique node ids

=item B<--info_file_node_id_col>

Column in the info file for the node id [1-based]
(Default: 1)

=item B<--info_file_node_name_col>

Column in the info file for the node name [1-based]
(Default: Undefined)

=item B<--graphml_file>

Name of a file to output the graph to in graphml format
Attributes from the info file will be added as node attributes as well as the node id.

=item B<--graph_id>

An id for the graph in graphml file.
(Default: G)
        
=item B<--csv_graph_node_file>

Name of a file to output the nodes to in csv format.
(Default: undefined)

=item B<--csv_graph_edge_file>

Name of a file to output the edges to in csv format.
If --csv_graph_node_file is set and --csv_graph_edge_file is not,
the script will exit, unless the csv_graph_node_file contains 'node'.
In this case the script will substitute 'edge' for 'node'. 

(Default: undefined)

=item B<--debug>

Print debugging information.

=item B<--help>

Print a brief help message and exit.

=item B<--man>

Print this script's manual page and exit.

=back

=head1 DEPENDENCIES

None

=head1 AUTHOR

=over 4

=item *

Richard White <rich@buschlab.org>

=back

=head1 COPYRIGHT AND LICENSE

This software is Copyright (c) 2020. University of Cambridge.

This is free software, licensed under:

  The GNU General Public License, Version 3, June 2007

=cut