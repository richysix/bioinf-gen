#!/usr/bin/env perl

# PODNAME: gene_lists_from_groups_cluego.pl
# ABSTRACT: The script produces a file of gene ids/names for each Function Group
# from a ClueGOo analysis. The input file is the 'Functional Groups With Genes.txt'
# file output from a ClueGO analysis.
# The ClueGO output uses gene names. If a file of gene names to ids is provided
# it will be used to produce lists with ids as well as gene names.

use warnings;
use strict;
use Getopt::Long;
use autodie;
use Pod::Usage;
use File::Slurp;

# get options
my %options;
get_and_check_options();

# open names to ids if it exists
my %gene_id_for; # gene_id_for = { GENE_NAME => GENE_ID }
if ( defined $options{'gene_names2ids'} ) {
    my @lines = read_file( $options{'gene_names2ids'}, chomp => 1, ) ;
    foreach my $line ( @lines ) {
    #foreach my $line ( @lines ) {
      my @fields = split /\t/, $line;
      $gene_id_for{ $fields[ $options{'gene_name_col'} ] } =
        $fields[ $options{'gene_id_col'} ];
    }
}

# open sig genes file if it exists
my %sig_gene_ids; # %sig_gene_ids = { GENE_ID => 1 }
my %sig_gene_names; # %sig_gene_names = { GENE_NAME => 1 }
if ( defined $options{'sig_genes_file'} ) {
    my @lines = read_file( $options{'sig_genes_file'}, chomp => 1, ) ;
    foreach my $line ( @lines ) {
    #foreach my $line ( @lines ) {
      my @fields = split /\t/, $line;
      if (defined $options{'sig_gene_id_col'}) {
        $sig_gene_ids{ $fields[ $options{'sig_gene_id_col'} ] } = 1;
      }
      if (defined $options{'sig_gene_name_col'}) {
        $sig_gene_names{ $fields[ $options{'sig_gene_name_col'} ] } = 1;
      }
    }
}

# open input file
my $input_file = $ARGV[0];
my @groups;
my @lines = read_file( $input_file, chomp => 1 );
my $single_file_fh;
if ($options{'single_file'}) {
  my $filename = $options{'file_base'} . "_all.tsv";
  open $single_file_fh, '>', $filename;
}
foreach my $line ( @lines ) {
  next if $line =~ m/\A Function/xms;
  # get gene names for each group
  my ($group_name, $group_num, $genes) = split /\t/, $line;
  # replace spaces and slashes / with underscores
  $group_name =~ s/\s+/_/xmsg;
  $group_name =~ s/\//_/xmsg;
  
  my @genes = split /\|/, $genes;
  my $fh;
  my $group = join("-", $group_num, $group_name, );
  if (!$single_file_fh) {
    open $fh, '>', $options{'file_base'} . "_" . $group . ".tsv";
  }
  foreach my $gene ( @genes ) {
    # check whether the gene is significant or not if --sig_genes_file is provided
    my $keep = 0;
    if ($options{'sig_genes_file'}){
      if (defined $options{'sig_gene_name_col'} && exists $sig_gene_names{$gene} ) {
        $keep = 1;
      }
      if (defined $options{'sig_gene_id_col'} &&
          exists $gene_id_for{$gene} &&
          exists $sig_gene_ids{ $gene_id_for{$gene} } ) {
        $keep = 1;
      }
    }
    next if !$keep;
    
    my @out_fields = ();
    push @out_fields, $gene;
    if (exists $gene_id_for{$gene}) {
      push @out_fields, $gene_id_for{$gene};
    }
    if ($single_file_fh) {
      push @out_fields, $group;
      print $single_file_fh join("\t", @out_fields, ), "\n";
    } else {
      print $fh join("\t", @out_fields, ), "\n";
    }
  }
  if (!$single_file_fh) {
    close $fh;
  }
}

if ($single_file_fh) {
  close $single_file_fh;
}

################################################################################
# SUBROUTINES
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
        'gene_names2ids=s',
        'gene_name_col=i',
        'gene_id_col=i',
        'sig_genes_file=s',
        'sig_gene_id_col=i',
        'sig_gene_name_col=i',
        'file_base=s',
        'single_file',
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
    $options{'gene_name_col'} = $options{'gene_name_col'} ? $options{'gene_name_col'} - 1 : 0;
    $options{'gene_id_col'} = $options{'gene_id_col'} ? $options{'gene_id_col'} - 1 : 1;
    $options{'sig_gene_name_col'} = $options{'sig_gene_name_col'} ? $options{'sig_gene_name_col'} - 1 : $options{'sig_gene_name_col'};
    $options{'sig_gene_id_col'} = $options{'sig_gene_id_col'} ? $options{'sig_gene_id_col'} - 1 : $options{'sig_gene_id_col'};
    $options{'file_base'} = $options{'file_base'} ? $options{'file_base'} : 'cluego_genes';
    
    # check there is a ids column number of sig_genes_file is defined
    if ($options{'sig_genes_file'}) {
      if (!defined $options{'sig_gene_id_col'} && !defined $options{'sig_gene_name_col'}) {
        pod2usage('At least one of --sig_gene_id_col or --sig_gene_name_col must be specified if --sig_genes_file is');
      }
    }
    
    print "Settings:\n", map { join(' - ', $_, defined $options{$_} ? $options{$_} : 'off'),"\n" } sort keys %options if $options{verbose};
}

__END__

=pod

=head1 NAME

gene_lists_from_groups_cluego.pl

=head1 DESCRIPTION

Description

=cut

=head1 SYNOPSIS

    gene_lists_from_groups_cluego.pl [options] input file | STDIN
        --gene_names2ids        file of gene names to Ensembl IDs
        --gene_name_col         number indicating the column of gene names (1-based)
        --gene_id_col           number indicating the column of gene ids (1-based)
        --sig_genes_file        a file containing the significant genes in the experiment.
        --sig_gene_name_col     number indicating the column of gene names in the sig genes file (1-based)
        --sig_gene_id_col       number indicating the column of gene ids in the sig genes file (1-based)
        --file_base             a base for each file name
        --single_file           produce a single file instead of one for each group
        --help                  print this help message
        --man                   print the manual page
        --debug                 print debugging information
        --verbose               turn on verbose output


=head1 ARGUMENTS

=over

=item input_file

The 'Functional Groups With Genes.txt' output from a ClueGO analysis

=back

=head1 OPTIONS

=over

=item B<--gene_names2ids>

File containing the mapping of gene names to Ensembl gene ids.

=item B<--gene_name_col>

The number of the column containing the gene names (1-based)
default: 1

=item B<--gene_id_col>

The number of the column containing the gene ids (1-based)
default: 2

=item B<--sig_genes_file>

A file containing the significant genes that the ClueGO analysis used.

=item B<--sig_gene_name_col>

The number of the column in the significant genes file containing the gene names (1-based)
One of either --sig_gene_name_col or --sig_gene_id_col must be specified if --sig_genes_file is.

=item B<--sig_gene_id_col>

The number of the column in the significant genes file containing the gene ids (1-based).
One of either --sig_gene_name_col or --sig_gene_id_col must be specified if --sig_genes_file is.

=item B<--file_base>

A base for each file name
default: cluego_genes

=item B<--single_file>

Produce a single file instead of one for each group. In this case, the name
of the group will be in the last column. The default filename is
<FILE_BASE> '_all.tsv'

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

Richard White <richard.white@sanger.ac.uk>

=back

=head1 COPYRIGHT AND LICENSE

This software is Copyright (c) 2018 by Genome Research Ltd.

This is free software, licensed under:

  The GNU General Public License, Version 3, June 2007

=cut