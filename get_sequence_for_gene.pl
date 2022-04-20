#!/usr/bin/env perl

# PODNAME: get_sequence_for_gene.pl
# ABSTRACT: Script to retrieve sequences from the Ensembl database

use warnings;
use strict;
use Getopt::Long;
use autodie;
use Pod::Usage;
use Carp;

use Bio::EnsEMBL::Registry;
use Bio::Seq;
use Bio::SeqIO;

# get options
my %options;
get_and_check_options();

# Connnect to Ensembl database
Bio::EnsEMBL::Registry->load_registry_from_db(
    -host => $options{'ensembl_dbhost'},
    -port => $options{'ensembl_dbport'},
    -user => $options{'ensembl_dbuser'},
    -pass => $options{'ensembl_dbpass'},
);

# Get genebuild version
my $genebuild_version = 'e' . Bio::EnsEMBL::ApiVersion::software_version();
warn 'Genebuild version: ', $genebuild_version, "\n" if $options{'debug'};

# Get Ensembl adaptor
my $gene_adaptor = Bio::EnsEMBL::Registry->get_adaptor($options{'species'}, 'core', 'Gene');

# Ensure database connection isn't lost; Ensembl 64+ can do this more elegantly
## no critic (ProhibitMagicNumbers)
if ( Bio::EnsEMBL::ApiVersion::software_version() < 64 ) {
## use critic
    Bio::EnsEMBL::Registry->set_disconnect_when_inactive();
}
else {
    Bio::EnsEMBL::Registry->set_reconnect_when_lost();
}

# create output files
my $dna_fasta_out = Bio::SeqIO->new(
    -file     => '>' . $options{'dna_fa_filename'},
    -alphabet => 'dna',
    -format   => 'fasta',
);
my $aa_fasta_out;
if ($options{'peptide'}) {
    $aa_fasta_out = Bio::SeqIO->new(
        -file     => '>' . $options{'peptide_fa_filename'},
        -alphabet => 'protein',
        -format   => 'fasta',
    );
}

# go through genes ids and get sequences
while(<>) {
    my $gene_id = $_;
    chomp $gene_id;
    if ($gene_id !~ m/\A ENS [A-Z]* G [0-9]{11} \z/xms) {
        die("Input, $gene_id, isn't a valid Ensembl ID [ENS [A-Z]* G [0-9]{11}]");
    }
    my $gene = $gene_adaptor->fetch_by_stable_id( $gene_id );
    output_transcript_and_protein_sequence( $gene, $options{'species'}, $dna_fasta_out, $aa_fasta_out, $options{'individual_files'});
}

################################################################################
# SUBROUTINES

# output_transcript_and_protein_sequence
#
#  Usage       : output_transcript_and_protein_sequence( $gene, $species, $dna_fasta_out, $aa_fasta_out, $individual_files )
#  Purpose     : Print the dna and aa sequences for the canonical transcript for a gene
#  Returns     : 
#  Parameters  : Gene               EnsEBML::Gene object
#                Species            Str
#                dna_output_file    Bio::SeqIO object
#                aa_output_file     Bio::SeqIO object
#                individual_files   Bool
#  Throws      : If biotype of gene is not 'protein_coding'
#  Comments    : None

sub output_transcript_and_protein_sequence {
    my ($gene, $species, $dna_fasta_out, $aa_fasta_out, $individual_files) = @_;
    
    if ($gene->biotype() ne 'protein_coding') {
        croak "$gene is not protein coding";
    }
    foreach my $transcript ( @{ $gene->get_all_Transcripts } ) {
        next if !$transcript->is_canonical();
        # output individual fasta files if option is set
        my $label = $gene->stable_id . "_" . join("_", $transcript->stable_id, $species);
        $label =~ s/\s/_/xmsg;
        
        my $dna_seq = Bio::Seq->new(
            -id       => $label,
            -alphabet => 'dna',
            -seq      => $transcript->translateable_seq(),
        );
        $dna_fasta_out->write_seq($dna_seq);

        my $pep_seq;
        if ($aa_fasta_out) {
            $pep_seq = Bio::Seq->new(
                -id       => $label,
                -alphabet => 'protein',
                -seq      => $transcript->translate()->seq(),
            );
            $aa_fasta_out->write_seq($pep_seq);
        }
        
        if($individual_files) {
            # create output files
            my $ind_dna = $species . '-' . $gene->stable_id . '.fa';
            my $ind_dna_out = Bio::SeqIO->new(
                -file     => '>' . $ind_dna,
                -alphabet => 'dna',
                -format   => 'fasta',
            );
            $ind_dna_out->write_seq($dna_seq);
            
            if ($aa_fasta_out) {
                my $ind_aa = $species . '-' . $gene->stable_id . '-pep.fa';
                my $ind_aa_out = Bio::SeqIO->new(
                    -file     => '>' . $ind_aa,
                    -alphabet => 'protein',
                    -format   => 'fasta',
                );
                $ind_aa_out->write_seq($pep_seq);
            }
        }
    }
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
        'species=s',
        'peptide',
        'individual_files',
        'dna_fa_filename=s',
        'peptide_fa_filename=s',
        'ensembl_dbhost=s',
        'ensembl_dbport=s',
        'ensembl_dbuser=s',
        'ensembl_dbpass=s',
        'help',
        'man',
        'debug+',
        'verbose',
    ) or pod2usage(2);
    
    # Documentation
    if( $options{'help'} ) {
        pod2usage( -verbose => 0, -exitval => 1, );
    }
    elsif( $options{'man'} ) {
        pod2usage( -verbose => 2 );
    }
    
    $options{'debug'} = $options{'debug'} ? $options{'debug'} : 0;
    $options{'species'} = $options{'species'} ? $options{'species'} : 'danio_rerio';
    $options{'dna_fa_filename'} = $options{'dna_fa_filename'} ? $options{'dna_fa_filename'} : 'ensembl_seqs.fa';
    if ($options{'peptide'}) {
        $options{'peptide_fa_filename'} = $options{'peptide_fa_filename'} ? $options{'peptide_fa_filename'} : 'ensembl_seqs-pep.fa';
    }
    
    # database defaults
    $options{'ensembl_dbhost'} = $options{'ensembl_dbhost'} ? $options{'ensembl_dbhost'} : 'ensembldb.ensembl.org';
    $options{'ensembl_dbuser'} = $options{'ensembl_dbuser'} ? $options{'ensembl_dbuser'} : 'anonymous';
    
    print "Settings:\n", map { join(' - ', $_, defined $options{$_} ? $options{$_} : 'off'),"\n" } sort keys %options if $options{'verbose'};
}

__END__

=pod

=head1 NAME

get_sequence_for_gene.pl

=head1 DESCRIPTION

Takes a list of Ensembl ids and a species and outputs fasta files for ORF DNA and PEPTIDE sequences.
The default is to get the ORF DNA sequence for the canonical/longest transcript.
Gene ids should be supplied on STDIN or as a filename of gene ids.

=cut

=head1 SYNOPSIS

    get_sequence_for_gene.pl [options] input file | STDIN
        --species               species Default: danio_rerio
        --peptide               output the peptide sequence as well
        --individual_files      output an individual file for each sequence
        --dna_fa_filename       name of the overall DNA sequence file
        --peptide_fa_filename   name of the overall peptide sequence file
        --ensembl_dbhost        Ensembl Host name
        --ensembl_dbport        Ensembl Port number
        --ensembl_dbuser        Ensembl username
        --ensembl_dbpass        Ensembl password
        --help                  print this help message
        --man                   print the manual page
        --debug                 print debugging information
        --verbose               turn on verbose output


=head1 ARGUMENTS

=over

Names of files containing Ensembl gene ids. Alternatively, gene ids can be supplied on STDIN.

=back

=head1 OPTIONS

=over

=item B<--species SPECIES>

Species to use (defaults to Danio rerio).

=item B<--peptide>

Whether to output the pepide sequence as well

=item B<--individual_files>

Output a seperate file for each sequence as well as one with all in

=item B<--dna_fa_filename>

Name of the output fasta file for all the DNA sequences.

=item B<--peptide_fa_filename>

Name of the output fasta file for all the peptide sequences.

=item B<--ensembl_dbhost HOST>

Ensembl MySQL database host.

=item B<--ensembl_dbport PORT>

Ensembl MySQL database port.

=item B<--ensembl_dbuser USERNAME>

Ensembl MySQL database username.

=item B<--ensembl_dbpass PASSWORD>

Ensembl MySQL database password.

=item B<--debug>

Print debugging information.

=item B<--help>

Print a brief help message and exit.

=item B<--man>

Print this script's manual page and exit.

=back

=head1 DEPENDENCIES

Ensembl API perl modules

=head1 AUTHOR

=over 4

=item *

Richard White <rich@buschlab.org>

=back

=head1 COPYRIGHT AND LICENSE

This software is Copyright (c) 2022. Queen Mary University of London.

This is free software, licensed under:

  The GNU General Public License, Version 3, June 2007

=cut