#!/usr/bin/perl

## ---------------------------------------------------------------

## extract info from gff3 file for genes that are specified in an
##  external list; extract min, max, and mean positions of CDS and
##  strandedness

## additionally outputs a list with alternative transcripts per unique gene

## Assuming following gff3 format definition:
## http://gmod.org/wiki/GFF3#GFF3_Format
## https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md

## Example: gff3Info3.pl -l SP1_goodIDs.txt -g SP1.gff3 -o SP1_map.txt

## ---------------------------------------------------------------

use warnings;
use strict;

use File::Basename;
use Getopt::Long;
##use List::MoreUtils 'any';


## set defaults
my $topclass = 'gene'; ## 'gene' (has to come before mRNA in gff3 file;
##                         its 'ID' defines 'Parent' of mRNA)
my $baseclass = 'mRNA'; ## 'mRNA' (has to come before CDS in gff3 file;
##                         its 'ID' defines 'Parent' of CDS)
my $baseclass2 = 'krkrkrkrkrkrkr'; ## second string to define baseclass 
##                         (partial match)
my $cds = 'CDS';
my $genenamestr = 'ID'; ## ('Name' might not be unique on different scaffs)
my $cdsnamestr = 'Name'; ## 'Name' (might not exist -> use 'ID')
my $idstr = 'ID'; ## 'ID'
my $parentstr = 'Parent'; ## 'Parent'
my $col = 1; ## column in protein list to use, \s+ sep
my $iso = 0; ## split off '.1'-style isoform to find gene matches
my $format = 1;

my ($genes, $gff, $out);
GetOptions(
    'l|L=s'  => \$genes,
    'g|G=s'  => \$gff,
    'o|O=s'  => \$out,
    'u|U=s'  => \$topclass,
    't|T=s'  => \$baseclass,
    'v|V=s'  => \$baseclass2,    
    'm|M=s'  => \$genenamestr,
    'n|N=s'  => \$cdsnamestr,
    'c|C=i'  => \$col,
    'i|I=i'  => \$iso,
    'x|X=i'  => \$format,
    'h|help' => \&usage
);

&usage if (!defined($genes) || !defined($gff) || !defined($out));

if($format == 1){
    $genenamestr = $genenamestr."=";
    $cdsnamestr = $cdsnamestr."=";
    $idstr = $idstr."=";
    $parentstr = $parentstr."=";
}
## print "$cdsnamestr\n";
## ---------------------------------------------------------------------

## open genes file
open (GENES, $genes) or die "Couldn't read $genes file\n";

## open gff3 file
open (GFF, $gff) or die "Couldn't read $gff file\n";

## open outfile
open (OUT, ">", "$out") or die "Couldn't open $out to write to\n";


if($format == 1){
    my $transcripts = $out;
    $transcripts =~ s/\.txt//;
    $transcripts .= "_transcripts.txt";
    open (TRANS, ">", "$transcripts") or die "Couldn't open $transcripts to write to\n";
}

## read in the genes list and prepare hashes to store gene info
my %models;
my @geneorder = ();
my @line = ();
my $nam;
my $fullnam;
my %minpos;
my %cnt;

while(<GENES>){
    chomp $_;
    @line = split(/\s+/,$_);
    $fullnam = $line[$col - 1];
    $nam = $fullnam;
    if($iso == 1){
	$nam =~ s/\.[0-9]+//;
    }
    if(exists($models{$nam})){
	die "Error: gene $nam is not unique\n";
    }
    push (@geneorder, $nam);
    $models{$nam} = $fullnam;
    $minpos{$nam} = undef;
    $cnt{$nam} = 0;
}
close(GENES);

my %maxpos = %minpos;
my %meanpos = %minpos;
my %nseq = %minpos;
my %strand = %minpos;
my %chr = %minpos;

my %uniGenes;
my %uniRnas;
my @allgenes = ();
my $newRnaId = '';
my $newGeneId = '';


if($format == 1){
    ## read elements from gff3 file
    my $cnt = 0;
    my $line = '';
    my @fields = ();
    my @attrib = ();
    my %ids;
    my $cdsParent = '';
    my $rnaParent = '';
    my $unknownGene = '';
    my $unknownRna = '';
    my $first = 0;
    my $last = 0;

    while($line = <GFF>){
	chomp $line;
	$cnt++;
	unless($line =~ /^#/){
	    @fields = split (/\t/,$line);
	    if($fields[2] eq $topclass){ ## assuming this feature is always
		## before mRNA and its 'ID=' provides the 'Parent=' of the mRNA
		$fields[8] =~ s/;+$//; ## remove trailing ';', if any
		@attrib = split (/;/,$fields[8]);
		## get gene name
		my $namefield = &getIndex($cnt,$topclass,$genenamestr,\@attrib);
		$nam = $attrib[$namefield]; ## Name= in attributes field, or
		##                             as specified with -m flag
		$nam =~ s/$genenamestr//;
		## get gene ID
		my $IDfield = &getIndex($cnt,$topclass,$idstr,\@attrib);
		$newGeneId = $attrib[$IDfield];
		$newGeneId =~ s/$idstr//; ## id in first attributes field
		if(exists($uniGenes{$newGeneId})){
		    die "Error: $topclass $newGeneId is not unique\n";
		}
		push (@allgenes, $newGeneId);
		$uniGenes{$newGeneId} = $nam.":";
		$unknownGene = '';
		$newRnaId = '';
		$unknownRna = '';
	    }
	    elsif($fields[2] eq $baseclass){ ## assuming this feature is always
		## before CDS and its 'ID=' provides the 'Parent=' of the CDS
		$fields[8] =~ s/;+$//; ## remove trailing ';', if any
		@attrib = split (/;/,$fields[8]);
		## get mRNA ID
		my $IDfield = &getIndex($cnt,$baseclass,$idstr,\@attrib);
		$newRnaId = $attrib[$IDfield];
		$newRnaId =~ s/$idstr//; ## id in first attributes field
		if(exists($uniRnas{$newRnaId})){
		    die "Error: $baseclass $newRnaId is not unique\n";
		}
		$uniRnas{$newRnaId} = '';
		## get mRNA parent
		my $parentfield = &getIndex($cnt,$baseclass,$parentstr,\@attrib);
		$rnaParent = $attrib[$parentfield];
		$rnaParent =~ s/$parentstr//; ## Parent= in attributes field
		## check names
		if($rnaParent ne $newGeneId){
		    print "$gff is not sorted for $topclass and $baseclass for $rnaParent at line $cnt - check output file\n";
		    $unknownGene = $rnaParent;
		}
		else{
		    $unknownGene = '';
		}
		$unknownRna = '';
	    }
	    elsif($fields[2] =~ m/$baseclass2/){ 
                ## alternative partial match, assuming this feature is always
		## before CDS and its 'ID=' provides the 'Parent=' of the CDS
		$fields[8] =~ s/;+$//; ## remove trailing ';', if any
		@attrib = split (/;/,$fields[8]);
		## get mRNA ID
		my $IDfield = &getIndex($cnt,$baseclass2,$idstr,\@attrib);
		$newRnaId = $attrib[$IDfield];
		$newRnaId =~ s/$idstr//; ## id in first attributes field
		if(exists($uniRnas{$newRnaId})){
		    die "Error: $baseclass2 $newRnaId is not unique\n";
		}
		$uniRnas{$newRnaId} = '';
		## get mRNA parent
		my $parentfield = &getIndex($cnt,$baseclass2,$parentstr,\@attrib);
		$rnaParent = $attrib[$parentfield];
		$rnaParent =~ s/$parentstr//; ## Parent= in attributes field
		## check names
		if($rnaParent ne $newGeneId){
		    print "$gff is not sorted for $topclass and $baseclass2 for $rnaParent at line $cnt - check output file\n";
		    $unknownGene = $rnaParent;
		}
		else{
		    $unknownGene = '';
		}
		$unknownRna = '';
	    }
	    elsif($fields[2] eq $cds){
		@attrib = split (/;/,$fields[8]);
		## get CDS name
		my $namefield = &getIndex($cnt,$cds,$cdsnamestr,\@attrib);
		if($namefield eq ''){
		    next;
		}
		$nam = $attrib[$namefield]; ## Name= in attributes field, or
		##                             as specified with -n flag
		$nam =~ s/$cdsnamestr//;
		## get CDS parent
		my $parentfield = &getIndex($cnt,$cds,$parentstr,\@attrib);
		$cdsParent = $attrib[$parentfield];
		$cdsParent =~ s/$parentstr//; ## Parent= in attributes field
		## check names
		if($cdsParent eq $newGeneId){ ## parent is gene, not mRNA
		    print "Warning: no $baseclass or $baseclass2 for $cds at line $cnt\n";
		    $unknownRna = '';
		}
		elsif($cdsParent ne $newRnaId){
		    print "$gff is not sorted for $baseclass or $baseclass2 and $cds for $cdsParent at line $cnt - check output file\n";
		    $unknownRna = $cdsParent;
		}
		else{
		    $unknownRna = '';
		}

		if (exists($models{$nam})){
		    $first = $fields[3];
		    $last = $fields[4];
		    if ($cnt{$nam}==0){ ## first CDS
			$cnt{$nam}++;
			$ids{$nam} = $cdsParent;
			## set strandedness
			$strand{$nam} = $fields[6];
			## set chromosome/scaffold id
			$chr{$nam} = $fields[0];
			$minpos{$nam} = $first;
			$maxpos{$nam} = $last;
			$meanpos{$nam} = $first + ($last - $first)/2.0;
			$nseq{$nam} = 1;
			## add name to genes
			if (($unknownGene eq '') && ($unknownRna eq '')){
			    $uniGenes{$newGeneId} .= " ".$models{$nam};
			}
		    }
		    elsif ($ids{$nam} eq $cdsParent){ ## next CDS
			if($minpos{$nam} > $first){
			    $minpos{$nam} = $first;
			}
			if($maxpos{$nam} < $last){
			    $maxpos{$nam} = $last;
			}
			$meanpos{$nam} += $first + ($last - $first)/2.0;
			$nseq{$nam}++;
		    }
		}
	    }
	}
    }
}

close(GFF);

## print outfile
print OUT "name\tchr\tstart\tend\tmean\tstrand\n";
for $nam (@geneorder){
    print OUT $models{$nam}."\t";
    if(defined($minpos{$nam})){
	print OUT $chr{$nam}."\t";
	print OUT $minpos{$nam}."\t";
	print OUT $maxpos{$nam}."\t";
	printf OUT ("%.0f",$meanpos{$nam} / $nseq{$nam});
	print OUT "\t".$strand{$nam}."\n";
    }
    else{
	print OUT "NA\tNA\tNA\tNA\tNA\n";
    }
}

close(OUT);

if($format == 1){
    for $newGeneId (@allgenes){
	unless($uniGenes{$newGeneId} =~ /:$/){
	    print TRANS $uniGenes{$newGeneId}."\n";
	}
    }
    close(TRANS);
}

exit;



## Get index in attributes field
sub getIndex{
    my ($cnt,$type,$pattern,$attref) = @_;
    my $idx = '';
    my @attr = @{$attref};
    for(my $j=0;$j<scalar(@attr);$j++){
	if($attr[$j] =~ /$pattern/){
    	    $idx = $j;
    	}
    }
    if($idx eq ''){
    	##die "Error: $pattern does not exist as attribute of $type in line $cnt\n";
	print "Warning: $pattern does not exist as attribute of $type in line $cnt\n";
    }
    return $idx;
}



## Show usage
## -------------------------------------------------------------------------
sub usage{
    print "\n";
    print "  Usage:\n";
    print "    ".basename($0)."\n";
    print "      -l <list of protein names>\n";
    print "      -g <gff3 file>\n";
    print "      -o <outfile>\n";
    print "      -u <type of feature that defines parent of mRNA (default=gene)>\n";
    print "      -t <type of feature that defines parent of CDS (default=mRNA)>\n";
    print "      -v <alternative type of feature that defines parent of CDS\n";
    print "          as a partial match (default=krkrkrkrkrkrkr)>\n";
    print "      -m <identifier in attributes field of gene that will \n";
    print "          be used in output (default=ID)>\n";
    print "      -n <identifier in attributes field of CDS that will \n";
    print "          be used to match name in protein list (default=Name)>\n";
    print "      -c <column in protein list with name (s+ sep, default=1)>\n";
    print "      -i <split off '.1'-style isoform to find gene matches\n";
    print "          (default=0)>\n";
    print "      -x <format (default=1, gff3)>\n";
    print "\n";
    exit;
}
## -------------------------------------------------------------------------
