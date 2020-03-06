#!/usr/bin/env perl

# PODNAME: mutmap_extract_peaks.pl
# ABSTRACT: Extract and rank variation from mutation mapping peaks

## Author         : is1
## Maintainer     : is1
## Created        : 2013-05-05
## Last commit by : $Author$
## Last modified  : $Date$
## Revision       : $Revision$
## Repository URL : $HeadURL$

use warnings;
use strict;
use autodie;
use Carp;
use Try::Tiny;
use Readonly;
use YAML::Tiny;

use Getopt::Long;
use Pod::Usage;
use File::Slurp;

=head1 DESCRIPTION


=head1 EXAMPLES


=cut

# Default options
my $mutant_sample_name = 'MUT';
my $sibling_sample_name = 'SIB';
my $formatting_yaml = 'excel-format.yml';
my ( $debug, $help, $man );

# Get and check command line options
get_and_check_options();

# open input file
open my $in_fh, '<', $ARGV[0];
# find max number of fields
my $max_fields = 0;
while( my $line = <$in_fh> ) {
    my @fields = split /\t/, $line;
    $max_fields = scalar @fields > $max_fields ? scalar @fields : $max_fields;
}
close($in_fh);

# output header
my @header_fields = ("Distance from peak", "Chr", "Pos", "Ref", "Alt", "FILTER",
                     $mutant_sample_name, $sibling_sample_name );
my $num_fields_before_csq = scalar @header_fields;
for(my $i = $num_fields_before_csq; $i < $max_fields; $i++) {
    push @header_fields, "Consequence" . ($i - $num_fields_before_csq + 1)
}
push @header_fields, "FILT";
print join("\t", @header_fields, ), "\n";

# open input file again
open $in_fh, '<', $ARGV[0];
Readonly my $FILTER_FIELD => 5;

if ($debug){ warn $FILTER_FIELD, "\n" }

while( my $line = <$in_fh> ) {
    chomp $line;
    my @fields = split /\t/, $line;
    
    # make extra column to extract FILTER info
    my $filter = "";
    if ($fields[$FILTER_FIELD] =~ m/LOW_QUAL;GT_INCONSISTENT/){
        $filter = "LOW_QUAL;GT_INCONSISTENT"
    } elsif ($fields[$FILTER_FIELD] =~ m/LOW_QUAL/){
        $filter = "LOW_QUAL"
    } elsif ($fields[$FILTER_FIELD] =~ m/GT_INCONSISTENT/){
        $filter = "GT_INCONSISTENT"
    } elsif ($fields[$FILTER_FIELD] =~ m/PASS/){
        $filter = "PASS"
    }
    
    my @output = @fields[0..($num_fields_before_csq - 1)];
    for(my $i = $num_fields_before_csq; $i < $max_fields; $i++) {
        push @output, defined $fields[$i] ? $fields[$i] : "";
    }
    push @output, $filter;
    print join("\t", @output, ), "\n"; 
}

# output YAML file for conditional formatting in Excel
my @colours_for = (
    {
        text => 'PASS',
        fontColour => "#375623",
        bgFill => "#C6E0B4"
    },
    {
        text => 'LOW_QUAL',
        fontColour => "#833C0C",
        bgFill => "#ED7D31"
    },
    {
        text => 'GT_INCONSISTENT',
        fontColour => "#9C6500",
        bgFill => "#FFEB9C"
    },
    {
        text => 'LOW_QUAL;GT_INCONSISTENT',
        fontColour => "#595959",
        bgFill => "#C9C9C9"
    }
);

sub get_Excel_column_from_number {
    my ($max_fields, $debug) = @_;
    my @letters = ( qw{ A B C D E F G H I J K L M N O P Q R S T U V W X Y Z} );
    my $denom = scalar @letters;
    Readonly my @LETTERS => ( "", @letters );
    
    my $column_letter = "";
    my $first_letter_idx = int( ($max_fields - 1) / $denom );
    $column_letter .= $LETTERS[ $first_letter_idx ];
    my $letter_idx = ( ($max_fields - 1 ) % $denom ) + 1;
    $column_letter .= $LETTERS[ $letter_idx ];
    
    if ($debug) {
        warn join(q{ }, $max_fields, $first_letter_idx, $letter_idx, $column_letter, ), "\n";
    }
    return($column_letter);
}

if ($formatting_yaml) {
    my $yaml = YAML::Tiny->new;
    $yaml->[0]->{formatting} = [];
    my $i = 0;
    my $column_letters = get_Excel_column_from_number($max_fields + 1);# add one to max fields because last field will be FILT
    foreach my $info ( @colours_for ) {
        my $formatting_rule = '$' . $column_letters . q{1=="} . $info->{'text'} . q{"};
        $yaml->[0]->{formatting}->[$i] = {
            type => 'expression',
            rule => $formatting_rule,
            fontColour => $info->{'fontColour'},
            bgFill => $info->{'bgFill'},
        };
        $i++;
    }
    $yaml->write( $formatting_yaml );
}

# Get and check command line options
sub get_and_check_options {

    # Get options
    GetOptions(
        'mutant_sample_name=s'   => \$mutant_sample_name,
        'sibling_sample_name=s'  => \$sibling_sample_name,
        'formatting_yaml=s'      => \$formatting_yaml,
        'debug'                  => \$debug,
        'help'                   => \$help,
        'man'                    => \$man,
    ) or pod2usage(2);

    # Documentation
    if ($help) {
        pod2usage(1);
    }
    elsif ($man) {
        pod2usage( -verbose => 2 );
    }

    return;
}

=head1 USAGE

    mutmap_extract_peaks.pl
        [--mutant_sample_name name]
        [--sibling_sample_name name]
        [--formatting_yaml file]
        [--debug]
        [--help]
        [--man]

=head1 OPTIONS

=over 8

=item B<--mutant_sample_name STRING>

Name of mutant sample. Default = MUT

=item B<--sibling_sample_name STRING>

Name of sibling sample. Default = SIB

=item B<--formatting_yaml FILE>

name of YAML file to output for conditional formatting of the final Excel file.

=item B<--debug>

Print debugging information.

=item B<--help>

Print a brief help message and exit.

=item B<--man>

Print this script's manual page and exit.

=back

=cut

