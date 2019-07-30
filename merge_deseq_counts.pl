#!/usr/bin/env perl

# PODNAME: merge_deseq_counts.pl
# ABSTRACT: The script takes a file of the filenames to merge. Optionally,
# a gene list can be supplied to subset the count files with.
# The script can get either normalised or unnormalised counts.

use warnings;
use strict;
use Getopt::Long;
use autodie;
use Pod::Usage;
use File::Slurp;


# get options
my %options;
get_and_check_options();

# open gene list if it exists
# otherwise use all the genes in the first file
my $select_gene; # this hash has all the Ensembl ids to get counts for
if ( defined $options{'gene_list'} ) {
    my @lines = read_file( $options{'gene_list'}, chomp => 1, ) ;
    $select_gene = {
        map { $_ => 1 } @lines
    };
}

# open file of filenames
my $fofn = $ARGV[0];
my @files = read_file( $fofn, chomp => 1 );

#  get count info genes from each file
my $got_gene_info = 0;
my %info_for;
my %counts_for;
my @header;
foreach my $file ( @files ) {
    # open file
    open my $counts_fh, '<', $file;
    my @columns;
    my @info_columns;
    my $line = <$counts_fh>;
    chomp $line;
    my @info = split /\t/, $line;
    for ( my $i = 0; $i < scalar @info; $i++ ) {
        if( $info[$i] =~ m/normalised\scount \z/xms ) {
            if( $options{'normalised'} ) {
                push @columns, $i;
            } else {
                push @info_columns, $i;
            }
        } elsif( $info[$i] =~ m/count \z/xms &&
                    $info[$i] !~ m/normalised/xms ) {
            if( !$options{'normalised'} ) {
                push @columns, $i;
            } else {
                push @info_columns, $i;
            }
        } else {
            push @info_columns, $i;
        }
    }
    warn "$file @columns\n" if $options{'debug'};
    
    if( !$got_gene_info ) {
        push @header, @info[@info_columns];
    }
    push @header, @info[@columns];
    
    # if select_gene is not defined and this is the first file
    # all the ids in this file will be placed in the select_gene hash.
    my $select = !defined $select_gene ? 0 : 1;
    while( my $line = <$counts_fh> ) {
        chomp $line;
        my @info = split /\t/, $line;
        if( !$select ) {
            $select_gene->{$info[0]} = 1;
        } else {
            next if ( !exists $select_gene->{$info[0]} );
        }
        
        my $gene_id = $info[0];
        # get the information about the gene (chr, start, end etc.) if this is
        # the first file
        if( !$got_gene_info ) {
            $info_for{$gene_id} = [ @info[@info_columns] ];
        }
        push @{$counts_for{$gene_id}}, @info[@columns];
    }
    $select = 1;
    $got_gene_info = 1; # set this to 1 so that we don't get the gene information again
}

print join("\t", @header ), "\n";
foreach my $gene_id ( sort keys %counts_for ) {
    print join("\t", @{$info_for{$gene_id}}, @{$counts_for{$gene_id}} ), "\n";
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
        'gene_list=s',
        'normalised',
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
    
    $options{debug} = $options{debug} ? $options{debug} : 0;
    print "Settings:\n", map { join(' - ', $_, defined $options{$_} ? $options{$_} : 'off'),"\n" } sort keys %options if $options{verbose};
}

__END__

=pod

=head1 NAME

merge_deseq_counts.pl

=head1 DESCRIPTION

Description

=cut

=head1 SYNOPSIS

    merge_deseq_counts.pl [options] input file | STDIN
        --gene_list             list of gene ids to limit to
        --normalised            get normalised counts instead of counts
        --help                  print this help message
        --man                   print the manual page
        --debug                 print debugging information
        --verbose               turn on verbose output


=head1 ARGUMENTS

=over

=item file_of_filenames

A file containing the filenames to get the counts from for merging

=back

=head1 OPTIONS

=over

=item B<--gene_list>

File containing a list of Ensembl gene ids to limit to.
If this is not supplied the script will take the gene ids from the first file
in the file of filenames and use that.

=item B<--normalised>

Get normalised counts instead of counts

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