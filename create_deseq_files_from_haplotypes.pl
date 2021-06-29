#!/usr/bin/env perl

# PODNAME: create_deseq_files_from_haplotypes.pl
# ABSTRACT: Creates the samples files for running DESeq2 based on sample
#           haplotypes provdided in a haplotypes file

use warnings;
use strict;
use Getopt::Long;
use autodie;
use Pod::Usage;
use File::Slurp;
use File::Spec;
use File::Copy;
use File::Path qw(make_path remove_tree);

# get options
my %options;
get_and_check_options();

# make output dir if it doesn't exist
if( ! -e $options{'dir'} ){
    mkdir($options{'dir'}) or
        die('Could not create directory, ', $options{'dir'}, ".\n");
}

# open samples file and get sample names
my @samples_lines = read_file( $options{'samples_file'} );
my @samples = ();
foreach my $line ( @samples_lines ){
    chomp $line;
    my ($sample_name, $condition) = split /\t/, $line;
    next if $condition eq 'condition';
    push @samples, $sample_name;
}
if( $options{'debug'} > 1){
    use Data::Dumper;
    warn Dumper(@samples);
}

# go through haplotypes file
open my $ht_fh, '<', $ARGV[0];
my $header_line = <$ht_fh>;
chomp $header_line;
my $col_for = parse_header( $header_line );
if( $options{'debug'} > 1){
    warn Dumper($col_for);
}

while(<$ht_fh>){
    chomp;
    my @cols = split /\t/, $_;
    my $region = join("_", map { $cols[ $col_for->{$_} ] } qw{chr start end} );
    
    # create new directory for region
    my $dir = File::Spec->catfile($options{'dir'}, 'deseq2-' . $region );
    if( ! -e $dir ){
        mkdir($dir) or die('Could not create directory, ', $dir, ".\n");
    }
    
    # create samples file
    # first get conditions for samples and reorder samples by condition
    # makes plotting easier
    my %samples_for = ();
    foreach my $sample ( @samples ){
        my $condition = $cols[ $col_for->{$sample} ];
        push @{$samples_for{$condition}}, $sample;
    }
    if( $options{'debug'} > 1){
        warn Dumper(%samples_for);
    }
    my @conditions = grep { !/^NA$/ && !/^NC$/ } sort keys %samples_for;
    # check for more than 1 condition
    if( scalar @conditions == 1 ){
        if ($options{randomise_single_conditions}) {
            my @conditions = map { 'group' . $_ } ( 1..$options{randomise_single_conditions} );
            %samples_for = ();
            foreach my $sample ( @samples ){
                my $condition = $conditions[ int( rand @conditions ) ];
                push @{$samples_for{$condition}}, $sample;
            }
        } else {
            warn sprintf
                "Region: %s. Only one condition in embryos. Removing directory ...\n",
                    $region;
            remove_tree($dir, { verbose => 1, });
            next;
        }
    }
    
    # make sample file
    my $samples_file = File::Spec->catfile($dir, 'samples.txt');
    open my $sample_fh, '>', $samples_file;
    print $sample_fh join("\t", '', 'condition'), "\n";
    foreach my $condition ( @conditions ){
        #next if $condition eq "NA" || $condition eq "NC";
        foreach my $sample ( sort sort_by_row @{$samples_for{$condition}} ){
            print $sample_fh join("\t", $sample, $cols[ $col_for->{$sample} ] ),
                "\n";
        }
    }
    close($sample_fh);

    # copy count file to new dir
    my $new_count_file = File::Spec->catfile($dir, 'counts.txt');
    copy($options{'count_file'}, $new_count_file ) or
        die('Could not copy counts file, ', $options{'count_file'}, ".\n");
    
    # create shell script to run deseq for all comparisons
    my $script = File::Spec->catfile($dir, 'run_deseq.sh');
    open my $script_fh, '>', $script;
    my $cmd = <<"CMD";
export R_LIBS_USER=/software/team31/R
perl ~rw4/checkouts/bio-misc/run_deseq2_rnaseq.pl \\
--output_dir %s \\
--counts_file %s \\
--samples_file %s \\
--comparisons %s \\
--nostrip_condition_prefix
CMD
    my @comparisons = ();
    for(my $i = 0; $i < scalar @conditions; $i++ ){
        for(my $j = 0; $j < scalar @conditions; $j++ ){
            next if $j <= $i;
            push @comparisons, join(":", $conditions[$i], $conditions[$j]);
        }
    }
    print $script_fh sprintf $cmd, $dir, $new_count_file, $samples_file,
                                join(q{ }, @comparisons);
    close($script_fh);
    
}

################################################################################
# SUBROUTINES

# parse_header
#
#  Usage       : $col_for = parse_header( header_line )
#  Purpose     : 
#  Returns     : Hashref of column_numbers for column names
#  Parameters  : header_line, Str. Tab-separated string of column names
#  Throws      : 
#  Comments    : None

sub parse_header {
    my ( $header_line ) = @_;
    my @col_names = split /\t/, $header_line;
    my %col_num_for = ();
    for ( my $i = 0; $i < scalar @col_names; $i++ ){
        $col_num_for{$col_names[$i]} = $i;
    }
    return(\%col_num_for);
}

sub sort_by_row {
    my ($row_a, $col_a) = get_row_and_col($a);
    my ($row_b, $col_b) = get_row_and_col($b);
    $row_a cmp $row_b || $col_a <=> $col_b;
}

sub get_row_and_col {
    my ($name, ) = @_;
    my $row = $name;
    $row =~ s/\A .*_([A-P])[0-9]+ \z/$1/xms;
    my $col = $name;
    $col =~ s/\A .*_[A-P]([0-9]+) \z/$1/xms;
    return($row, $col);
}

# get_and_check_options
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
        'samples_file=s',
        'count_file=s',
        'dir=s',
        'randomise_single_conditions+',
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
    $options{samples_file} = $options{samples_file} ?
        $options{samples_file} : 'deseq2/samples.txt';
    $options{count_file} = $options{count_file} ?
        $options{count_file} : 'deseq2/count.txt';
    $options{dir} = $options{dir} ? $options{dir} : './';
    $options{randomise_single_conditions} = $options{randomise_single_conditions}
        ? $options{randomise_single_conditions} + 1: 0;
    $options{debug} = $options{debug} ? $options{debug} : 0;
    
    print "Settings:\n",
        map {join(' - ', $_, defined $options{$_} ? $options{$_} : 'off'),"\n" }
            sort keys %options if $options{verbose};
}

__END__

=pod

=head1 NAME

create_deseq_files_from_haplotypes.pl

=head1 DESCRIPTION

Description

=cut

=head1 SYNOPSIS

    create_deseq_files_from_haplotypes.pl [options] haplotypes_file
        --samples_file          Name of base samples file [default: deseq2/samples.txt]
        --count_file            Name of the counts file [default: deseq2/samples.txt]
        --dir                   Output directory
        --help                  print this help message
        --man                   print the manual page
        --debug                 print debugging information
        --verbose               turn on verbose output


=head1 ARGUMENTS

=over

=item I<haplotypes_file>

File of sample haplotypes. Each line is a region and the sample haplotypes are in columns.
Must contain a header with chr, start, end columns and the sample names must match the samples file.

=back

=head1 OPTIONS

=over

=item B<--samples_file>

Base samples file. Only sample names present in this file are used. The conditions are ignored.

=item B<--count_file>

Count file. This is copied to the new directory to make running DESeq2 easier.

=item B<--dir>

Output directory. Directory within which to create the deseq directories.

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