#!/usr/bin/perl

use strict;

use Getopt::Long;
use LWP::UserAgent;
use Net::SFTP::Foreign;

my $query;
my $search;
my %argOpts = ();
my $uniprotContent;
my $uniprotFlattFile;
my %uniprotHash;
my %cosmicHash;

my %queryField =(
    'a'      =>  'accession:',
    'e'      =>  'mnemonic:',
    'r'      =>  'cluster:',
);


GetOptions (\%argOpts, 'h', 'a=s', 'e=s','r=s','ens=s','o:s','email=s' );
&help if (defined $argOpts{'h'});

die "Pass accession id with '-a' OR Pass entry name with '-e' OR Pass refseq id with '-r' OR Pass ensembly id with '-ens'\nPlease, pass help message with '-h'\n"  unless ( defined $argOpts{'a'} | defined $argOpts{'e'} | defined $argOpts{'r'} | defined $argOpts{'ens'});
die "Please enter your available email adress for COSMIC\n" unless (defined $argOpts{'email'});
die "Please enter output file name with -o\n" unless (defined $argOpts{'o'});

#if ( defined $argOpts{'r'} ) { $query = $queryField{'r'}.$argOpts{'r'}; $search = $argOpts{'r'}; }
if ( defined $argOpts{'e'} ) { $query = $queryField{'e'}.$argOpts{'e'}; $search = $argOpts{'e'}; }
if ( defined $argOpts{'a'} ) { $query = $queryField{'a'}.$argOpts{'a'}; $search = $argOpts{'a'}; }


my $url = "http://www.uniprot.org/uniprot/?query=$query&columns=id,entry name,genes,sequence,length,feature(NATURAL VARIANT),feature(ALTERNATIVE SEQUENCE),feature(MODIFIED RESIDUE)&format=tab&sort=score";

my $agent = LWP::UserAgent->new;
my $response = $agent->get($url);

if($response->is_success) {
    $uniprotContent = $response->decoded_content;
    print "Query file Loaded!\n";
}else{
    print "Unsuccessfully load operation!\n";
}

my @entries;
my $outputFile;
my %ensHash;
my $isoID;
my %variantHash;
my %content;

my @lines = split("\n", $uniprotContent);
foreach my $line (@lines){
    #print $line."\n";
    if ($line =~ m/^P.*/) {
        my @columns = split("\t", $line);
        my $entry = $columns[0];
        $outputFile = $search.".txt";
        unless (-e $outputFile) {
            $url = "http://www.uniprot.org/uniprot/$entry.txt";
            $response = $agent->get($url);
            open(my $file, '>' , $outputFile);
                if($response->is_success) {
                    print $file $response->content;
                    print "$search.txt File Downloaded!\n";
                }else{
                    die "Unsuccesfully downloaded!\n";
                }
            close $file;
        }else{
            print "$outputFile already exists!\n";
        }
        my $entry_name = $columns[1];
        my $genes = $columns[2];
        my @genes = split(" ", $genes);
        my $gene = $genes[0];
        my $sequence = $columns[3];
        my $length = $columns[4];
        my $variantsCol = $columns[5];
        $content{$entry}{$entry_name} = {
            gene => $gene,
            sequence => $sequence,
            len => $length,
            variantCol => $variantsCol
        }
    }
}

if ($isoID = &check_canoncial){
    if(%ensHash = &get_ensembl($outputFile, $isoID)){
        foreach my $enstID (keys %ensHash){
            foreach my $entry (keys %content){
                foreach my $entry_name (keys $content{$entry}){
                    if (%variantHash = &parse_variant($content{$entry}{$entry_name}{'variantCol'})){
                        foreach my $featureID (keys %variantHash){
                            foreach my $mutation (keys $variantHash{$featureID}){
                                my $geneEnsembl = $content{$entry}{$entry_name}{'gene'}."_".$enstID;
                                $uniprotHash{$geneEnsembl}{$mutation} ={
                                        entry => $entry,
                                        entryName => $entry_name,
                                        gene => $content{$entry}{$entry_name}{'gene'},
                                        enstID => $enstID,
                                        isoID => $isoID,
                                        featureID => $featureID,
                                        len => $content{$entry}{$entry_name}{'len'},
                                        sequence => $content{$entry}{$entry_name}{'sequence'},
                                        description => $variantHash{$featureID}{$mutation},
                                        ensp => $ensHash{$enstID}{'ensp'},
                                        ensg => $ensHash{$enstID}{'ensg'}
                                }
                            }
                        }
                    }
                }
            }
        }
        }
}else{
    die "Its not canonical\n";
}

open(my $file,'>',$argOpts{'o'});
foreach my $geneEnsembl (sort keys %uniprotHash){
    if(my %cosmicHash = &cosmic_part($geneEnsembl,$argOpts{'email'})){
        foreach my $mutation (sort keys $uniprotHash{$geneEnsembl}){
            if (exists $cosmicHash{$geneEnsembl}{$mutation}){
                print $file "$geneEnsembl\t$mutation\t";
                foreach my $valueUniprot (sort keys $uniprotHash{$geneEnsembl}{$mutation}){
                    print $file "$uniprotHash{$geneEnsembl}{$mutation}{$valueUniprot}\t";
                }
                foreach my $valueCosmic (sort keys $cosmicHash{$geneEnsembl}{$mutation}){
                    if(ref($cosmicHash{$geneEnsembl}{$mutation}{$valueCosmic}) eq 'ARRAY'){
                        print $file join( ";", @{$cosmicHash{$geneEnsembl}{$mutation}{$valueCosmic}});
                    }else{
                        print $file $cosmicHash{$geneEnsembl}{$mutation}{$valueCosmic}."\t";
                    }
                }
                print $file "\n\n";
            }
        }
    }
}
print $file "\n";
close $file;

sub cosmic_part {
    my @args = @_;
    my $geneEnsemblUniprot = $args[0];
    my $user = $args[1];
            if(my $file  = &cosmic_sftp($geneEnsemblUniprot,$user)){
                my @file = split("\n", $file);
                foreach my $row (@file){
                    my @columns = split("\t", $row);
                    my $geneEnsemblCosmic = $columns[0];
                    my $mutationCosmic = $columns[1];
                    my $primarySite = $columns[2];
                    if (exists $cosmicHash{$geneEnsemblCosmic}{$mutationCosmic}){
                        $cosmicHash{$geneEnsemblCosmic}{$mutationCosmic}{'count'}++;
                        push (@{$cosmicHash{$geneEnsemblCosmic}{$mutationCosmic}{'primarySite'}}, $primarySite);
                    }else{
                        $cosmicHash{$geneEnsemblCosmic}{$mutationCosmic} = {
                            siteSubtype1 => $columns[3],
                            siteSubtype2 => $columns[4],
                            siteSubtype3 => $columns[5],
                            primaryHistology => $columns[6],
                            histologySubtype1 => $columns[7],
                            histologySubtype2 => $columns[8],
                            histologySubtype3 => $columns[9],
                            cosmicID => $columns[10],
                            description => $columns[11],
                            tumorOrigin => $columns[12],
                            primarySite => [$primarySite],
                            count => 1,
                        }
                    }
                }
            }
    return %cosmicHash;
}

sub check_canoncial {
    open(my $file,'<',$outputFile);
    while(my $row = <$file>){
        if (($isoID, my $flag) = $row =~ m/^CC\s+IsoId\=([A-Z]+[0-9]+\-[0-9]+)\;\s+Sequence\=(.*)\;/) {
            if ($flag eq "Displayed"){
                return $isoID;
            }
        }
    }
    close $file;
}


sub get_ensembl {
    my @args = @_;
    $outputFile = $args[0];
    $isoID = $args[1];
    open(my $file,'<',$outputFile);
        while(my $row = <$file>){
            if ( (my $enstID, my $enspID, my $ensgID) = $row =~ /^DR\s+Ensembl\;\s+([A-Z]+[0-9]+)\;\s+([A-Z]+[0-9]+)\;\s+([A-Z]+[0-9]+)\.\s+\[($isoID)\]/){
                $ensHash{$enstID}={
                    ensg => $ensgID,
                    ensp => $enspID
                }
            }
        }
        return %ensHash;
    close $file;
}

sub parse_variant {
    my @args = @_;
    my $variantsCol = $args[0];
    my @variants = split(/\.\; /, $variantsCol);
    foreach my $variant (@variants){
        if( (my $position1,my $position2,my $fromMutation,my $toMutation,my $description,my $featureID) = $variant =~ /^VARIANT\s([0-9]*)\s([0-9]*)\s([A-Z]+) -> ([A-Z]+)\s(.*)\/FTId\=(VAR\_[0-9]+)/){
            my $mutation = "p.".$fromMutation.$position1.$toMutation;
            $variantHash{$featureID}{$mutation} = $description;
        }
    }
    return %variantHash;
}

sub cosmic_sftp {
    my @args = @_;
    my $word = $args[0];
    my $user = $args[1];
    my $cosmicFile = "CosmicMutantExport.tsv.gz";
    my $host = "sftp-cancer.sanger.ac.uk";
    my $cosmicOutput = "CosmicMutant.tsv";
    my %args = (
        user     => $user,
        port      => 22,
        # ...
    );


    unless (-e "$cosmicFile") {
        #check file is downloaded or not before.
        print "Cosmic Mutant file is downloading. This will get some time. Please wait!\n";
        my $sftp = Net::SFTP::Foreign->new($host, %args);
        $sftp->error and die "Something bad happened: " . $sftp->error;
        $sftp->get("/files/grch38/cosmic/v79/$cosmicFile") or die "put failed: " . $sftp->error;
    }

    print "Cosmic Mutant file is extracting...\n";
    my $command = `gunzip -c $cosmicFile | awk -F["\t"] '{print \$1\"\t\"\$19\"\t\"\$8\"\t\"\$9\"\t\"\$10\"\t\"\$11\"\t\"\$12\"\t\"\$13\"\t\"\$14\"\t\"\$15\"\t\"\$17\"\t\"\$20\"\t\"\$34}' | grep ^$word`;

    return $command;

}


sub help {
    my $script = "$0";

    print<<"END";

  Usage: perl $script [Options]

  Options:
    -h                    Output usage message.

    -a accession          Accession ID represents ENTRY in the UniProt.

    -e entryname          Entry name(names & taxonomy) in the UniProt.

    -r refseq             Refseq(Sequence Database).

    -ens ensembl          Ensembl(Genome annotation database).

    -o filename           Output file.

END

    exit;
}
