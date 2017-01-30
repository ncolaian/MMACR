#! /usr/bin/env perl

# Need to load BLAST and muscle

use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;
use Carp;
use Readonly;
use Path::Class;
use Data::Dumper;
use Log::Log4perl qw(:easy);
use Log::Log4perl::CommandLine qw(:all);

# My Variables
my $help = 0;
my $man = 0;
my $sample_name_file;
my $min_align_length;
my $read_dir_path;
my $output_path;

# Read in the variable from the command line
GetOptions ( 'man'  =>  \$man,
            'help'  =>  \$help,
            'sample_file|sf=s'  =>  \$sample_name_file,
            'min_len|ml=i'      =>  \$min_align_length,
            'read_dir|rd=s'     =>  \$read_dir_path,
            'out_path|o=s'     =>  \$output_path,
            ) || die("There was an error in the command line arguements \n");

# Pod Usage for the manual and help pages
if ( $help ) { pod2usage(0) }
if ( $man ) { pod2usage(-verbose => 3) }

# Setup logging environment
my $logger = get_logger();

#EXPERIMENT SPECIFIC
#abyss assembly users
Readonly::Hash my %ABYSS_USERS => map { $_ => 1 } qw(
    albost_Cambodia
    immigrans
    kep_Brunei
    koh_Sarawak
    taxonF
    taxonG
);
my $abyss_dir = "/proj/cdjones_lab/projects/waddell/assembly/abyss_assemblies";

# MAIN #

$logger->info("Creating an array with sample names");
my $sample_fileo = file($sample_name_file);
my @sample_names = $sample_fileo->slurp( chomp=>1 );

# Make Blast Databases for all the genomes
make_blast_db( \@sample_names );

# Prgrm path that retrieves sequences from a fasta file
my $parser = "/nas02/home/n/c/ncolaian/my_scripts/phylogeny_scripts/pull_out_bres_portion.pl";
my $starter = $sample_names[0];
#Find the conserved regions within all the assemblies
find_conserved_regions( $parser, $starter, \@sample_names );

# Pull out only the conserved sequences that are over min_length
# Will also print out stats on the sequences over the min_length
# Returns an array ref of the region names
my $region_names_aref = pull_out_final_conserved( $starter, @sample_names );

#create blast databases for each individual region
create_individual_gene_bdb( $region_names_aref );

#run blast on each genome for each gene.
run_blasts_on_each_genome( $region_names_aref, \@sample_names );

# Decide which genomes to run the alignments on and update the stats
my ( $used_aref, $excluded_aref ) = decide_and_create_prealignment_file( $region_names_aref, \@sample_names, $parser);

#create used region and excluded region file
create_used_and_exclude_file( $used_aref, $excluded_aref );

#Run muscle alignemnt on each multple_align file
run_muscle( $used_aref );

# Combine the alignment files into one massive alignment file
combine_alignment_files( $used_aref );


# Subroutines #

sub stall {
    my ( $job_name, $out_dir ) = @_;
    $logger->info("Stalling until <$job_name> jobs finish");
    
    my $cmd = "bjobs -J " . $job_name . " > $out_dir/ACTIVE_JOBS";
    system( $cmd );
    
    #If file is empty then stop stalling
    if ( -s "$out_dir/ACTIVE_JOBS" ) {
        return 1;
    }
    
    # Remove ACTIVE_JOBS file
    `rm "$out_dir/ACTIVE_JOBS"`;
    
    return 0; # send the continue to stall signal
}

sub make_blast_db {
    my ( $sample_aref ) = @_;
    $logger->info("Making blast databases for all genomes");
    
    foreach my $name ( @$sample_aref ) {
        my $bdb_cmd = "bsub -o $output_path/lsf_bdb.out -J make_bdb makeblastdb ";
        if ( $ABYSS_USERS{ $name } ) {
            $bdb_cmd .= "-in $abyss_dir/$name.fasta ";
        }
        else {
            $bdb_cmd .= "-in $read_dir_path/$name/contigs.fasta "
        }
        $bdb_cmd .= "-title $name -dbtype nucl -out $output_path/blast_dbs/$name/$name";
        mkdir "$output_path/blast_dbs";
        mkdir "$output_path/blast_dbs/$name";
        $logger->debug($bdb_cmd);
        system( $bdb_cmd );
    }
    while( stall( "make_bdb", $output_path) ) {
        sleep 20;
    }
    return 1;
}

sub find_conserved_regions {
    my ( $parse_pth, $start_file_name, $sample_aref ) = @_;
    my $start_file = "$read_dir_path/$start_file_name/contigs.fasta";
    $logger->info( "Finding the conserved sequences using $start_file as the original file" );
    
    for( my $i = 1; $i < scalar(@$sample_aref); $i++ ) {
        my $nname = $sample_aref->[$i];
        #Blastn
        my $cmd = "bsub -J blast -o $output_path/lsf_blast.out blastn -db $output_path/blast_dbs/$nname/$nname -query $start_file -out $output_path/blast_results/br_$i " . '-outfmt "6 qseqid sseqid qstart qend score length nident"' . " -max_target_seqs 25";
        mkdir "$output_path/blast_results";
        $logger->debug( $cmd );
        system( $cmd );
        while( stall( "blast", $output_path ) ){
            sleep 60;
        }
        
        #Update query file
        my $p_cmd = "bsub -J parser -o $output_path/lsf_parse.out perl $parse_pth -btf $output_path/blast_results/br_$i -af $start_file -o $output_path/updated_query_files/conserved_$i";
        mkdir "$output_path/updated_query_files";
        $logger->debug( $p_cmd );
        system( $p_cmd );
        
        #Update query file name
        $start_file = "$output_path/updated_query_files/conserved_$i";
        while( stall( "parser", $output_path ) ) {
            sleep 30;
        }
    }
    return $start_file;
}

sub pull_out_final_conserved {
    my ( $file_path, @sample_names ) = @_;
    $logger->info( "Getting conserved sequences and producing statistics" );
    
    my $file_num = scalar(@sample_names) - 1;
    my $final_conserved_fo = file("$output_path/updated_query_files/conserved_$file_num");
    my @cons_file = $final_conserved_fo->slurp( chomp=>1 );
    
    #keep track of region names
    my @region_names;
    mkdir "$output_path/indiv_regions/";
    #STATS
    my $max = 0;
    my $total_sequences=0;
    my $total_sequence_length=0;
    
    for( my $j = 0; $j<scalar(@cons_file); $j++ ) {
        next if ( $cons_file[$j] =~ qr/start=/ );
        my $length = length $cons_file[$j];
        if ( $length >= $min_align_length ) {
            my $prev_line = $cons_file[$j-1];
            #handle stats
            $total_sequences++;
            $total_sequence_length += $length;
            if ( $max < $length ) {
                $max = $length;
            }
            #print out file
            my @split_prev = split " ", $prev_line;
            my $reg_name = $split_prev[0];
            $reg_name =~ s/>//;
            mkdir "$output_path/indiv_regions/$reg_name";
            open my $OUT, ">", "$output_path/indiv_regions/$reg_name/orig_gene_file.txt";
            print $OUT "$prev_line\n", $cons_file[$j], "\n";
            push @region_names, $reg_name;
            close($OUT);
        }
    }
    # Print out the stats
    open my $STATS, ">", "$output_path/conserved_regions_stats.txt";
    
    print $STATS "These are the stats for the conserved sequences that are above $min_align_length\n\nTotal Sequences:\t$total_sequences\nTotal Sequence Length:\t$total_sequence_length\nMaximum Conserved Region:\t$max\nAverage Region Length:\t", $total_sequence_length/$total_sequences;
    
    close($STATS);
    
    return \@region_names;
}

sub create_individual_gene_bdb {
    my ( $names_aref ) = @_;
    $logger->info( "Creating individual region databases" );
    
    foreach my $names ( @$names_aref ) {
        my $ind_cmd = "bsub -o $output_path/indiv_regions/lsf_bdb.out -J make_bdb makeblastdb -in $output_path/indiv_regions/$names/orig_gene_file.txt -title $names -dbtype nucl -out $output_path/indiv_regions/$names/blast_db/$names";
        
        mkdir "$output_path/indiv_regions/$names/blast_db/";
        $logger->debug($ind_cmd);
        system($ind_cmd);
    }
    while( stall( "make_bdb", $output_path ) ) {
        sleep 60;
    }
    return 1;
}

sub run_blasts_on_each_genome {
    my ( $region_aref, $sample_aref ) = @_;
    $logger->info( "Running blast on each region" );
    
    foreach my $region ( @$region_aref ){
        mkdir "$output_path/indiv_regions/$region/blast_files";
        foreach my $sample ( @$sample_aref ) {
            my $blast_cmd = "bsub -J blast -o $output_path/indiv_regions/$region/lsf_blast.out ";
            if ( $ABYSS_USERS{$sample} ) {
                $blast_cmd .= "blastn -query $abyss_dir/$sample.fasta "
            }
            else {
                $blast_cmd .= "blastn -query $read_dir_path/$sample/contigs.fasta "
            }
            $blast_cmd .= "-db $output_path/indiv_regions/$region/blast_db/$region -out $output_path/indiv_regions/$region/blast_files/$sample.txt" . ' -outfmt "6 qseqid sseqid qstart qend score length nident"' . " -max_target_seqs 10";
            $logger->debug($blast_cmd);
            system( $blast_cmd );
        }
    }
    
    while( stall( "blast", $output_path) ) {
        sleep 60;
    }
}

sub decide_and_create_prealignment_file {
    my ( $region_aref, $sample_aref, $parse_pth ) = @_;
    $logger->info( "Deciding which regions to keep" );
    #regions used
    my @regions_used;
    my @regions_excluded;
    my $usable_region_length = 0;
    
    foreach my $region ( @$region_aref ) {
        #get the length of the conserved region
        open my $ORIG, "<", "$output_path/indiv_regions/$region/orig_gene_file.txt";
        my $first_line = <$ORIG>;
        my $sequence_line = <$ORIG>;
        chomp $sequence_line;
        my $reg_length = length($sequence_line);
        
        #determine if the region is usable
        my $check = 0;
        #check each samples file for a result
        my $single_line;
        foreach my $sample ( @$sample_aref ) {
            my $file = "$output_path/indiv_regions/$region/blast_files/$sample.txt";
            #check if file is empty
            if ( -z $file ) {
                $check = 1;
                last;
            }
            #check to see if one of the blast results matches length
            my $file_o = file($file);
            my @f_lines = $file_o->slurp( chomp=>1, split=>qr/\t/);
            my $single_line;
            for( my $f = 0; $f <= scalar(@f_lines); $f++) {
                if ( $f == scalar(@f_lines) ) {
                    $check = 1;
                    print "$sample stopped $region from being included in the analysis\n";
                    last;
                }
                if ( $f_lines[$f]->[5] <= ($reg_length + ($reg_length*.05)) &&
                    $f_lines[$f]->[5] >= ($reg_length - ($reg_length*.05)) ) {
                    $single_line = join "\t", @{$f_lines[$f]};
                    last;
                }
            }
            
            #re-write blast result file so only one gene exists
            if ($check == 0) {
                open my $REBLAST, ">", $file;
                print $REBLAST $single_line;
                close($REBLAST);
            }
        }
        if ($check != 0) {
            push @regions_excluded, $region;
            next;
        }
        push @regions_used, $region;
        $usable_region_length += $reg_length;#length;
    }
    
    #get the sequences for each gene
    $logger->info("Getting sequences for each gene then concatenating the file");
    #Running pull_bres
    foreach my $region ( @regions_used ) {
        my $psub = "bsub -J parser -o $output_path/indiv_regions/$region/lsf_parse.out" . ' "';
        foreach my $samp ( @$sample_aref ) {
            $psub .= "perl $parse_pth -btf $output_path/indiv_regions/$region/blast_files/$samp.txt";
            if ( $ABYSS_USERS{$samp} ) {
                $psub .= " -af $abyss_dir/$samp.fasta ";
            }
            else {
                $psub .= " -af $read_dir_path/$samp/contigs.fasta ";
            }
            $psub .= "-o $output_path/indiv_regions/$region/conserved_sequences/$samp.fasta -name $samp; ";
        }
        mkdir "$output_path/indiv_regions/$region/conserved_sequences/";
        $psub .= '"';
        $logger->debug($psub);
        system( $psub );
    }
    
    while( stall( "parser", $output_path ) ) {
        sleep 60;
    }
    
    #concatenate all the files
    foreach my $region ( @regions_used ) {
        my $cat_cmd = "cat $output_path/indiv_regions/$region/conserved_sequences/* > $output_path/indiv_regions/$region/pre_align_file.txt";
        $logger->debug( $cat_cmd );
        system( $cat_cmd );
    }
    
    #make sure all the sequences are on the same strand
    $logger->info("Ensuring the sequences are on the same strand");
    foreach my $region ( @regions_used ) {
        #Have to fix this path
        my $process_cmd = "bsub -J process -o $output_path/indiv_regions/$region/pro.out perl /nas02/home/n/c/ncolaian/my_scripts/utility_scripts/pre_mult_align_seq_flip.pl -af $output_path/indiv_regions/$region/pre_align_file.txt -o $output_path/indiv_regions/$region/rdy2_align_file.fasta";
        $logger->debug( $process_cmd );
        system( $process_cmd );
    }
    
    #update_stats_file
    $logger->info("updating stats file to include automated usable regions info");
    open my $STATS_UP, ">>", "$output_path/conserved_regions_stats.txt";
    my $num_reg = scalar( @regions_used );
    print $STATS_UP "\n\n**** USED REGIONS STATS ****\n\nRegions Used:\t$num_reg\n Used Sequence Length:\t$usable_region_length\n\nCheck the regions_used.txt and regions_excluded.txt for the id's of the regions not automatically added to the analysis";
    close($STATS_UP);
    
    while( stall( "process", $output_path) ) {
        sleep 30;
    }
    
    return ( \@regions_used, \@regions_excluded );
}

sub create_used_and_exclude_file{
    my ( $u_aref, $e_aref ) = @_;
    $logger->info("Creating region information files");
    
    open my $USED, ">", "$output_path/regions_used.txt";
    open my $EXCL, ">", "$output_path/regions_excluded.txt";
    
    print $USED join "\n", @$u_aref;
    print $EXCL join "\n", @$e_aref;
    
    close($USED);
    close($EXCL);
    return 1;
}

sub run_muscle {
    my ($regions_aref) = @_;
    $logger->info("Running muscle alignment on all the usable regions");
    
    foreach my $region ( @$regions_aref ) {
        my $mus_cmd = "bsub -J muscle -o $output_path/lsf_musc.out muscle -in $output_path/indiv_regions/$region/rdy2_align_file.fasta -out $output_path/indiv_regions/$region/aligned.fasta";
        $logger->debug( $mus_cmd );
        system($mus_cmd);
    }
    while( stall( "muscle", $output_path) ) {
        sleep 60;
    }
    return 1;
}

sub combine_alignment_files {
    my ( $used_aref ) = @_;
    $logger->info("Combining all the individual alignment files into one combined alignment file");
    
    my %align_holder_hash;
    
    #Fill in a hash that holds the combined alignment from each regions alignment
    foreach my $region ( @$used_aref ) {
        my $fo = file( "$output_path/indiv_regions/$region/aligned.fasta" );
        my @slurp_fo = $fo->slurp( chomp=>1 );
        my $samp_name;
        foreach my $line ( @slurp_fo ) {
            if ( $line =~ qr/>/ ) {
                my @names = split / /, $line;
                $samp_name = $names[0];
                $samp_name =~ s/>//;
                if ( !exists $align_holder_hash{$samp_name} ) {
                    $align_holder_hash{$samp_name} = "";
                }
            }
            else {
                $align_holder_hash{$samp_name} .= $line;
            }
        }
    }
    #iterate through hash and print out the massive alignment. It should be ordered by region
    open my $FINAL_OUT, ">", "$output_path/combined_alignment.fasta";
    foreach my $key ( keys %align_holder_hash ) {
        print $FINAL_OUT ">$key\n", $align_holder_hash{$key}, "\n";
    }
    close($FINAL_OUT);
    return 1;
}


__END__
=head1 MMACR

This is a tool to create a massive alignment file from alignments coming from conserved regions between a set of genomic assemblies. The resulting file will be a multiple alignment file that has combined all the indiviidual regions alignments.

=head1 VERSION

This documentation refers to MMACR 0.0.1

=head1 INCLUDED MODULES

Getopt::Long;
Pod::Usage;
Carp;
Readonly;
Path::Class;
Data::Dumper;
Log::Log4perl qw(:easy);
Log::Log4perl::CommandLine qw(:all);

=head1 INHERIT

=head1 SYNOPSIS

=head1 PARAMETERS

=head1 CONFIGURATION AND ENVIRONMENT

    Need to load muscle and BLAST

=head1 DEPENDENCIES

    muscle
    BLAST
    Master_aln
    
=head1 INCOMPATIBILITIES

    None reported.

=head1 BUGS AND LIMITATIONS

No bugs have been reported.

Please report any bugs or feature requests	
	
=head1 Author

Nicholas Colaianni
contact via C<< <ncolaian@live.unc.edu> >>

=head1 LICENCE AND COPYRIGHT

Copyright (c) 2016, Nicholas Colaianni
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met: 

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer. 
2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution. 

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

The views and conclusions contained in the software and documentation are those
of the authors and should not be interpreted as representing official policies, 
either expressed or implied, of the FreeBSD Project.

=head1 DISCLAIMER OF WARRANTY

BECAUSE THIS SOFTWARE IS LICENSED FREE OF CHARGE, THERE IS NO WARRANTY
FOR THE SOFTWARE, TO THE EXTENT PERMITTED BY APPLICABLE LAW. EXCEPT WHEN
OTHERWISE STATED IN WRITING THE COPYRIGHT HOLDERS AND/OR OTHER PARTIES
PROVIDE THE SOFTWARE "AS IS" WITHOUT WARRANTY OF ANY KIND, EITHER
EXPRESSED OR IMPLIED, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. THE
ENTIRE RISK AS TO THE QUALITY AND PERFORMANCE OF THE SOFTWARE IS WITH
YOU. SHOULD THE SOFTWARE PROVE DEFECTIVE, YOU ASSUME THE COST OF ALL
NECESSARY SERVICING, REPAIR, OR CORRECTION.

IN NO EVENT UNLESS REQUIRED BY APPLICABLE LAW OR AGREED TO IN WRITING
WILL ANY COPYRIGHT HOLDER, OR ANY OTHER PARTY WHO MAY MODIFY AND/OR
REDISTRIBUTE THE SOFTWARE AS PERMITTED BY THE ABOVE LICENCE, BE
LIABLE TO YOU FOR DAMAGES, INCLUDING ANY GENERAL, SPECIAL, INCIDENTAL,
OR CONSEQUENTIAL DAMAGES ARISING OUT OF THE USE OR INABILITY TO USE
THE SOFTWARE (INCLUDING BUT NOT LIMITED TO LOSS OF DATA OR DATA BEING
RENDERED INACCURATE OR LOSSES SUSTAINED BY YOU OR THIRD PARTIES OR A
FAILURE OF THE SOFTWARE TO OPERATE WITH ANY OTHER SOFTWARE), EVEN IF
SUCH HOLDER OR OTHER PARTY HAS BEEN ADVISED OF THE POSSIBILITY OF
SUCH DAMAGES.

=cut