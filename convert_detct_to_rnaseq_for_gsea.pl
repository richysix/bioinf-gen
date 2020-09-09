#!/usr/bin/env perl

# PODNAME: convert_detct_to_rnaseq_for_gsea.pl
# ABSTRACT: Script to convert DeTCT output file (all.tsv) to the same format as
#           RNAseq to allow conversion to GSEA format with
#           convert_deseq2_to_gsea.pl

use warnings;
use strict;
use Getopt::Long;
use autodie;
use Pod::Usage;

# get options
my %options;
get_and_check_options();

# open DeTCT file
open my $in_fh, '<', $ARGV[0];
my $header = <$in_fh>;
my @header = split /\t/, $header;
my $column_number_for = parse_header( @header );
# get the column number for the Log2 fold change column
my $log2fc_idx;
my $gene_id_idx;
my @count_cols = ();
my @norm_count_cols = ();
my @out_cols = ( qw{Gene pval adjp log2fc Chr Start End Strand
                    Biotype Name Description} );
foreach my $col ( keys %{$column_number_for} ){
    my $col_idx = $column_number_for->{$col};
    if ($col =~ m/\A Log\s2\sfold\schange .* \z/xms) {
        $log2fc_idx = $col_idx;
    } elsif ($col =~ m/ Gene\sID /xms) {
        $gene_id_idx = $col_idx;
    } elsif ($col =~ m/ count \z/xms ) {
      if ($col =~ m/ normalised /xms) {
          push @norm_count_cols, $col_idx;
      } else {
          push @count_cols, $col_idx;
      }
    }
}
push @out_cols, @header[@count_cols], @header[@norm_count_cols];

# Hash for keeping track of max absolute log2fc for each gene id
my %max_l2fc_for; # $max_l2fc_for{ GENE_ID => LOG2FC }
my %min_pval_for; # $min_pval_for{ GENE_ID => PVAL }
# Hash for keeping track of the region to use for each gene id
my %region_for; # $region_for{ GENE_ID => REGION }
while( my $line = <$in_fh> ){
    chomp $line;
    my @cols = split /\t/, $line;
    # skip regions where log2fc is NA
    next if $cols[$log2fc_idx] eq "NA";
    # make region name for line
    my $region = region_from_line(\@cols, $column_number_for);
    
    # need to split gene id column by comma
    foreach my $gene_id ( split /,/, $cols[$gene_id_idx] ){
        # check if pvalue is smaller than current smallest
        # and set values and region if applicable
        my $replace = 0;
        my $padj = $cols[ $column_number_for->{"Adjusted p value"} ];
        my $abs_log2fc = abs($cols[$log2fc_idx]);
        if (exists $min_pval_for{ $gene_id }) {
            if ($min_pval_for{ $gene_id } eq "NA") {
                if ($padj ne "NA") {
                    $replace = 1;
                } elsif ($abs_log2fc > $max_l2fc_for{ $gene_id }) {
                    $replace = 1;
                }
            } else {
                if ($padj ne "NA") {
                    if ($padj < $min_pval_for{ $gene_id }) {
                        $replace = 1;
                    }
                }
            }
        } else {
            $replace = 1;
        }
        
        if ($replace) {
            $min_pval_for{ $gene_id } = $padj;
            $max_l2fc_for{ $gene_id } = $abs_log2fc;
            $region_for{ $gene_id } = $region;
        }
        
    }
}

if ($options{"debug"}) {
  use Data::Dumper;
  warn Dumper(%max_l2fc_for, %region_for);
}

# close DeTCT file
close($in_fh);

# print header line
print join("\t", @out_cols), "\n";

# re-open DeTCT file
open $in_fh, '<', $ARGV[0];
$header = <$in_fh>;
while( my $line = <$in_fh> ){
    chomp $line;
    my @cols = split /\t/, $line;
    # make region name for line
    my $region = region_from_line(\@cols, $column_number_for);
    
    # need to split gene id column by comma
    my @gene_ids = split /,/, $cols[$gene_id_idx];
    my @gene_names = split /,/, $cols[$column_number_for->{"Gene name"}];
    my @gene_descriptions = split /,/, $cols[$column_number_for->{"Gene description"}];
    for (my $idx = 0; $idx < scalar @gene_ids; $idx++ ){
        my $gene_id = $gene_ids[$idx];
        my $gene_name = $gene_names[$idx];
        my $gene_description = $gene_descriptions[$idx];
        # check if region appears in the region_for hash for this gene
        if (exists $region_for{ $gene_id } &&
            $region_for{ $gene_id } eq $region) {
            print join("\t", $gene_id,
                       $cols[$column_number_for->{"p value"}],
                       $cols[$column_number_for->{"Adjusted p value"}],
                       $cols[$log2fc_idx],
                       $cols[$column_number_for->{"Chr"}],
                       $cols[$column_number_for->{"Region start"}],
                       $cols[$column_number_for->{"Region end"}],
                       $cols[$column_number_for->{"3' end strand"}],
                       "NA",
                       $gene_name,
                       $gene_description,
                       @cols[@count_cols],
                       @cols[@norm_count_cols],
                       ), "\n";
        }
    }
}


################################################################################
# SUBROUTINES

# parse_header
#
#  Usage       : parse_header( @header )
#  Purpose     : To find the column numbers for the required columns in the
#                output
#  Returns     : HASHREF of column numbers
#  Parameters  : header line split by tab into an array
#  Throws      : 
#  Comments    : None

sub parse_header {
    my @header = @_;
    my %column_number_for;
    
    for (my $i = 0; $i < scalar @header; $i++){
        $column_number_for{ $header[$i] } = $i;
    }
    
    return( \%column_number_for );
}

# region_from_line
#
#  Usage       : region_from_line(\@cols, $column_number_for)
#  Purpose     : Create a region id for a line of input
#  Returns     : STR (region id)
#  Parameters  : ARRAY input line split by tabs
#  Throws      : 
#  Comments    : None

sub region_from_line {
  my ($cols, $column_number_for) = @_;
  my @cols = @{$cols};
  my $region = join(":", $cols[ $column_number_for->{"Chr"} ],
                    $cols[ $column_number_for->{"Region start"} ],
                    $cols[ $column_number_for->{"Region end"} ],
                    $cols[ $column_number_for->{"3' end position"} ],
                    $cols[ $column_number_for->{"3' end strand"} ] );
  return($region);
}

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

convert_detct_to_rnaseq_for_gsea.pl

=head1 DESCRIPTION

This script converts a DeTCT output file to the same format as RNAseq data.
This then allows convert_deseq2_to_gsea.pl to run to produce the files
rquired for GSEA to run successfully. DeTCT output can have multiple
regions per gene id whereas for GSEA there should only be ONE entry
per gene. This script takes the region with the smallest pvalue.
If none of the regions have a p value it takes the one with the largest
log2 fold change.

=cut

=head1 SYNOPSIS

    convert_detct_to_rnaseq_for_gsea.pl [options] input file
        --help                  print this help message
        --man                   print the manual page
        --debug                 print debugging information
        --verbose               turn on verbose output


=head1 ARGUMENTS

=over

input file. This should de a DeTCT output file.

=back

=head1 OPTIONS

=over

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