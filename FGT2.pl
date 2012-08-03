#!/usr/bin/perl

#Based on Federico Abascal's "mitobank.pl" (http://darwin.uvigo.es/)
#Copyright (C) 2012 Evan McCartney-Melstad, based on Federico Abascal's
#mitobank.pl (http://darwin.uvigo.es/)
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>

use strict;
use warnings;
use Bio::Perl;
use Bio::DB::GenBank;
use Bio::SeqIO;
use Bio::AlignIO;
use Bio::Align::Utilities; ## for concatenating alignment files
use Bio::Root::IO;
use Bio::Seq;
use Bio::SimpleAlign;
use Tie::File;
use List::Util qw(first);
use Bio::Align::Utilities qw(:all);


#STEP ONE, PULL ALL RELEVANT mtGENOMES FROM GENBANK AND CREATE NUCLEOTIDE AND AMINO ACID FASTA FILES
#FOR ALL GENES AND rRNAs
#ALWAYS RUN SCRIPT WITH -nuc AND -rna SWITCHES (PROGRAM WILL PRINT BOTH AMINO ACIDS AND NUCLEOTIDES)

my $takeRNA = 0; 	#By default no RNA is processed
my $nuc = 0;		#By default amino acid sequences are taken.
my @taxids;
if(scalar(@ARGV) == 0) {
	&usage("");
}
while($ARGV[0] =~ /^\-/) {
	my $param = shift(@ARGV);
	if($param eq "-nuc") {
		$nuc = 1;
	} elsif ($param eq "-rna") {
		$takeRNA = 1;
	} else {
		&usage("\nParam \"$param\" not recognized.");
	}
} 
if(scalar(@ARGV) == 0) {
	&usage("\nAt least one NCBItaxid is required");
}
@taxids = @ARGV;

sub usage {
	my $error = shift(@_);
	print STDERR "$error\n";
	print STDERR "mitobank.pl, by Federico Abascal (2005), version 2.0 [http://darwin.uvigo.es/]\n";
	print STDERR "  Usage: perl mitobank.pl [options] list_of_NCBItaxids\n";
	print STDERR "         By default (no options specified) it takes the amino acid sequence\n";
	print STDERR "         of coding genes.\n";
	print STDERR "         Option -rna tells the program to obtain the sequence of RNA genes\n";
	print STDERR "         Option -nuc tells the program to take nucleotide sequences\n";
	print STDERR "         If you want to log output (recommended) redirect it with \"> name_of_file\"\n";
	print STDERR "  Example: to get the nucleotide sequence of all genes (coding and\n";
	print STDERR "         RNA genes) of Primates and Glires: \n";
	print STDERR "         perl mitobank.pl -nuc -rna 9443 314147 > genomes.log\n\n";
	exit;
}

#Make folder for storing the raw FASTA sequences and metadata (genbank file, gene order file, gene hits, taxon names vs. NCBI ID)
mkdir("FASTA_Files");
mkdir("FASTA_Files/nuc_pre-alignment");
mkdir("FASTA_Files/prot_pre-alignment");
mkdir("Other_Data");
mkdir("Other_Data/jmodeltest_results");
mkdir("FASTA_Files/prot_post-alignment");
mkdir("FASTA_Files/nuc_post-alignment");

my @genes_check = qw(ATP6 ATP8 COX1 COX2 COX3 CYTB ND1 ND2 ND3 ND4L ND4 ND5 ND6 12SrRNA 16SrRNA tRNA-Ala tRNA-Arg tRNA-Asn tRNA-Asp tRNA-Cys tRNA-Gln tRNA-Glu tRNA-Gly tRNA-His tRNA-Ile tRNA-Leu tRNA-Lys tRNA-Met tRNA-Phe tRNA-Pro tRNA-Ser tRNA-Thr tRNA-Trp tRNA-Tyr tRNA-Val);
foreach my $gene_check ( @genes_check ) {
	if(-e "FASTA_Files/nuc_pre-alignment/".$gene_check.".fas") {
		print STDERR "\nERROR: File $gene_check.fas already exists. Exiting...\n";
		print STDERR "Remove existing files or run the program in a different directory.\n\n";
		exit;
	}
}
if(-e "Other_Data/gene_order.txt") {
		print STDERR "\nERROR: File gene_order.txt already exists. Exiting... \n";
		print STDERR "Remove existing files or run the program in a different directory.\n\n";
		exit;
	}

my $qry_string;
foreach my $tmp (@taxids) {
	$qry_string .= "txid".$tmp." "; 
}
chop($qry_string);
$qry_string =~ s/ /\[Organism:exp\] OR /g;
$qry_string .= "[Organism:exp]";
$qry_string = "($qry_string)";
print "Query string: $qry_string AND mitochondri* AND complete genome NOT plasmid[title] NOT chromosome NOT chloroplast NOT synthetic construct\n";

my $seq;
my $gb = new Bio::DB::GenBank;
my $query = Bio::DB::Query::GenBank->new
#	(-query   =>'txid41705[Organism:exp] AND mitochondri*',   #	ProtacanthopterygiiÊ
#	(-query   =>'txid32443[Organism:exp] AND mitochondri*',   #	TeleosteiÊ
#	(-query   =>'txid9263[Organism:exp] AND mitochondri*',    #	MarsupialesÊ
	(-query   =>$qry_string . ' AND mitochondri* AND "complete genome" NOT plasmid[title] NOT chromosome NOT chloroplast NOT synthetic construct',    
	 -db      => 'nucleotide');


print "Query returned the following ", $query->count, " results:\n";


my(%genes,  %names, %yausados, %already_got); #, %classes);
my($specie, $taxid);
my $count = 0;
my $seqio = $gb->get_Stream_by_query($query);
#my $seqio  = new Bio::SeqIO(-file => "vertebrata.genbank", -format => 'genbank');
my $seqout = new Bio::SeqIO(-file => ">Other_Data/genomes.genbank", -format => 'genbank');
while(defined ($seq = $seqio->next_seq )) {
	$seqout->write_seq($seq);
	$specie = $seq->species->binomial   ;
	$taxid  = $seq->species->ncbi_taxid ;
	#$specie = $seq->species->common_name;
	my $tmp = $specie . ' 'x(30-length($specie)) . "[NCBI_TaxID: " . $seq->species->ncbi_taxid . "]";
	#my @classification = $seq->species->classification;
	my $num = 60-length($tmp);
	$count++;
	open(GENE_ORDER, ">>Other_data/gene_order.txt");
	print GENE_ORDER "\n>$tmp:\n";

###Decomment below to get at feature tags (translations, accession numbers, etc...)	
#	foreach my $feat_object ($seq->get_SeqFeatures) {
#		print "primary tag: ", $feat_object->primary_tag, "\n";
#		foreach my $tag ($feat_object->get_all_tags) {
#			print "  tag: ", $tag, "\n";
#			foreach my $value ($feat_object->get_tag_values($tag)) {
#				print "    value: ", $value, "\n";
#			}
#		}
#	}

	$specie =~ s/ /\_/g;
	my($first, $second) = (split(/\_/, $specie))[0,1];
	#my $name = substr($first,0,1) . "." . $second;
	my $name = $first . "_" . $second;
	while(exists($names{$name})) {
		$name = "_" . $name;
	}
	$names  {$name} = $taxid;
	#$classes{$name} = [@classification];
	foreach my $feat_object ($seq->get_SeqFeatures) {
		if ($feat_object->primary_tag eq "CDS") {
			if ($feat_object->has_tag('gene')) {
				my $val = ($feat_object->get_tag_values('gene'))[0];
				#some rules to correct naming inconsitencies:
				$val = uc($val);
				if($val eq "CO1" || $val eq "COXI" || $val eq "COI" || $val eq "COX-I")
				{	$val = "COX1";		} 
				elsif($val eq "CO2" || $val eq "COXII" || $val eq "COII" || $val eq "COX-II")
				{	$val = "COX2";		} 
				elsif($val eq "CO3" || $val eq "COXIII" || $val eq "COIII"  || $val eq "COX-III")
				{	$val = "COX3";		}
				elsif($val eq "ATP 6" || $val eq "ATPASE 6" || $val eq "ATPASE6")
				{	$val = "ATP6";		}
				elsif($val eq "ATP 8" || $val eq "ATPASE 8" || $val eq "ATPASE8")
				{	$val = "ATP8";		}
				elsif($val eq "CYT B" || $val eq "CYTOCHROME B")
				{	$val = "CYTB";		}				
				#Rules end.
				
				#Stop if there are results in the feature table that don't start with COX, ATP, ND, or CYT
				if($val !~ /^COX/ && $val !~/^ATP/ && $val !~ /^ND/ && $val !~ /^CYT/) {
					next;
				}
				open(GENE_ORDER, ">>Other_data/gene_order.txt");
				print GENE_ORDER "$val\t";
				close(GENE_ORDER);
				$genes{$val}++;
				next if(exists $already_got{$name}{$val});
				$already_got{$name}{$val} = 1;
				if($feat_object->has_tag('translation')) {
					open(O, ">>FASTA_Files/nuc_pre-alignment/$val.fas");   ##CREATE FASTA FILE FOR NUCLEOTIDE SEQUENCES
						#print O ">", $taxid, "\n";
						print O ">", $name, '_*_', $taxid, "\n";
						print O $feat_object->spliced_seq->seq, "\n"; 	##NUCLEOTIDE SEQUENCE
					close(O);
				}
			}
		} elsif($takeRNA == 1 && ($feat_object->primary_tag eq "tRNA" || $feat_object->primary_tag eq "rRNA")) {
				if(!$feat_object->has_tag('product')) {
					print STDERR "No TAG for: ", $name, " (feat_prim_tag: ", 
												 $feat_object->primary_tag, "\n";
					next;
				}
				my $val = ($feat_object->get_tag_values('product'))[0];
				#some rules to correct naming inconsitencies:
				if($val eq "12S ribosomal RNA" || $val eq "s-rRNA" || $val eq "s-RNA" || $val eq "small subunit ribosomal RNA" || $val eq "12 ribosomal RNA" || $val eq "12S ribosomal RNA subunit")
				{	$val = "12SrRNA";		} 
				elsif($val eq "16S ribosomal RNA" || $val eq "l-rRNA" || $val eq "l-RNA" || $val eq "large subunit ribosomal RNA" || $val eq "16S ribosamal RNA" || $val eq "16S ribosomal RNA subunit" || $val eq "16S rivbosomal RNA") 
				{	$val = "16SrRNA";		}
				elsif($val eq "tRNA-Glx") 
				{	$val = "tRNA-Glu";		}
				#rules end.
				print "$val ";
				$genes{$val}++;
				open(O, ">>FASTA_Files/nuc_pre-alignment/$val.fas");  ###create FASTA file for RNAs
					print O ">", $name, '_*_', $taxid, "\n";
					#print O ">", $taxid, "\n";
					print O $feat_object->spliced_seq->seq, "\n";
				close(O);
		}
	}    
	print "\n";
}
print "\n\nThere are: $count mitochondrial genomes for this taxonomy criterion\n";

open(O, ">Other_Data/number_of_results_per_gene.txt");
print O "Genes in dataset:\n";
foreach my $gene (keys %genes) {
	print O "$gene:\t$genes{$gene}\n";
}
close(O);

open(O, ">Other_Data/Species_names_vs_NCBI_taxID.txt");
print O "Names vs NCBI_TaxID:\n";
foreach my $tmp (keys %names) {
	print O "$tmp:\t$names{$tmp}\n";
}
close(O);

my @fasta_files_pre_alignment = ("ATP6.fas","ATP8.fas","COX1.fas","COX2.fas","COX3.fas","CYTB.fas","ND1.fas","ND2.fas","ND3.fas","ND4.fas","ND4L.fas","ND5.fas","ND6.fas");
my @RNA_fastas = ("12SrRNA.fas","16SrRNA.fas","tRNA-Ala.fas","tRNA-Arg.fas","tRNA-Asn.fas","tRNA-Asp.fas","tRNA-Cys.fas","tRNA-Gln.fas","tRNA-Glu.fas","tRNA-Gly.fas","tRNA-His.fas","tRNA-Ile.fas","tRNA-Leu.fas","tRNA-Lys.fas","tRNA-Met.fas","tRNA-Phe.fas","tRNA-Pro.fas","tRNA-Ser.fas","tRNA-Thr.fas","tRNA-Trp.fas","tRNA-Tyr.fas","tRNA-Val.fas");
my @all_fastas = ("ATP6.fas","ATP8.fas","COX1.fas","COX2.fas","COX3.fas","CYTB.fas","ND1.fas","ND2.fas","ND3.fas","ND4.fas","ND4L.fas","ND5.fas","ND6.fas","12SrRNA.fas","16SrRNA.fas","tRNA-Ala.fas","tRNA-Arg.fas","tRNA-Asn.fas","tRNA-Asp.fas","tRNA-Cys.fas","tRNA-Gln.fas","tRNA-Glu.fas","tRNA-Gly.fas","tRNA-His.fas","tRNA-Ile.fas","tRNA-Leu.fas","tRNA-Lys.fas","tRNA-Met.fas","tRNA-Phe.fas","tRNA-Pro.fas","tRNA-Ser.fas","tRNA-Thr.fas","tRNA-Trp.fas","tRNA-Tyr.fas","tRNA-Val.fas");

##GET RID OF DUPLICATE TAXA
for my $file (@all_fastas) {
	#Put all lines of file into an array--each line becomes an element. Array[0] is the first line of the file, etc...
	tie my @tied_array, 'Tie::File', "FASTA_Files/nuc_pre-alignment/$file", recsep => "\r\n" or die "Rats! Couldn't get the lines of the file into an array--is the file name correct?";

	#Go through each array element (each line of the file) to check if the array element contains >_
	#If it does, then get the index of that array element and splice it and the following element out of the array
	#the "following element" mentioned above is the corresponding DNA sequence
	for my $line(@tied_array) {
		my $index = first {$tied_array[$_] =~ m/>_/ } 0 .. $#tied_array;
		if ( length $index ){
		splice(@tied_array,$index,2);
		reset ('$index');
		}
			}
	untie @tied_array;
	}
	
######    #need to fix gene order, number_of_results_per_gene, species_names_vs_ncbi_taxID

#TRANSLATION
#Translate DNA sequences into amino acids and place them in the prot_pre-alignment folder
for my $file ( @fasta_files_pre_alignment ) {  ##Loop through all files in pre-alignment directory
	##Create SeqIO object to read in sequences from file
	my $seqio_obj = Bio::SeqIO->new(-file => "FASTA_Files/nuc_pre-alignment/$file", -format => "fasta" );
	my $seqio_prot_obj = Bio::SeqIO->new(-file => ">FASTA_Files/prot_pre-alignment/aa_$file", -format => 'fasta' ); ##Creates protein sequence SeqIO object
	
	##Translate each sequence and print to file one-by-one
	while (my $seq_obj = $seqio_obj->next_seq) {
		my $translate_prot = $seq_obj->translate(-codontable_id => 2);
		$seqio_prot_obj->write_seq($translate_prot);
		}
}

#MULTIPLE ALIGNMENT
#Muscle is in C:\Users\Evan\Desktop\My Dropbox\Columbia\Fish Genome Tool\perl scripts\FGT\subs
#Must ensure Muscle is in PATH for muscle to run
#output to CLUSTAL format for next step

#Align protein-coding amino acids
for my $file ( @fasta_files_pre_alignment ) {
	system("subs/muscle.exe", "-in", "FASTA_Files/prot_pre-alignment/aa_$file", "-out", "FASTA_Files/prot_post-alignment/aa_aligned_$file");
	}

#Align RNA nucleotides
for my $file ( @RNA_fastas ) {
	system("subs/muscle.exe", "-in", "FASTA_Files/nuc_pre-alignment/$file", "-out", "FASTA_Files/nuc_post-alignment/aligned_$file");
	}
	
#REVERSE TRANSLATE BACK TO ALIGNED DNA SEQUENCES
#use pal2nal.pl, by Mikita Suyama, David Torrents, and Peer Bork (2006). "PAL2NAL: robust conversion of protein sequence alignments into the corresponding codon alignments." Nucleic Acids Res. 34, W609-W612.
for my $file ( @fasta_files_pre_alignment ) {
	open(POSTALIGN, ">FASTA_Files/nuc_post-alignment/aligned_$file") || die "Can't create nucleotide post alignment file";
	close POSTALIGN;
	system("perl subs/pal2nal.pl FASTA_Files/prot_post-alignment/aa_aligned_$file FASTA_Files/nuc_pre-alignment/$file -output fasta -codontable 2 > FASTA_Files/nuc_post-alignment/aligned_$file");
	}

#close genomes.genbank filehandle
$seqout->close();



#Partitions, concatenations
#create in new subdirectory? and pull from subdirectory for future analyses?
#Create Seq objects from FASTA files, and concatenate using Bio::SeqUtils->cat(@seqs);
## Title     : cat
## Usage     : $aln123 = cat($aln1, $aln2, $aln3)
## Function  : Concatenates alignment objects. Sequences are identified by id.
##             An error will be thrown if the sequence ids are not unique in the
##             first alignment. If any ids are not present or not unique in any
##             of the additional alignments then those sequences are omitted from
##             the concatenated alignment, and a warning is issued. An error will
##             be thrown if any of the alignments are not flush, since
##             concatenating such alignments is unlikely to make biological
##             sense.
## Returns   : A new Bio::SimpleAlign object
## Args      : A list of Bio::SimpleAlign objects
####      {
####      my $counter = 0;
####      for my $file ( @all_fastas ) {
####      	my $temp_alignment = Bio::AlignIO->new(-file => $file, -format => 'fasta');
####      	my $temp_output_alignment = Bio::AlignIO->new(-file => ">outputfilename$counter", -format => 'nexus');
####      	while ( my $aln . $counter = $temp_alignment->next_aln ) {
####      	$temp_output_alignment->write_aln($aln . $counter);
####      	}
####      	$counter++;
####      
####      }
####      }


#THIS ONE WORKS (SORTA)
# my $ATP6_alignobject = Bio::AlignIO->new(-file => "FASTA_Files/nuc_post-alignment/aligned_ATP6.fas",-format => 'fasta');
# my $ATP8_alignobject = Bio::AlignIO->new(-file => "FASTA_Files/nuc_post-alignment/aligned_ATP8.fas",-format => 'fasta');
# my $ND1_alignobject = Bio::AlignIO->new(-file => "FASTA_Files/nuc_post-alignment/aligned_ND1.fas",-format => 'fasta');
# my $ND2_alignobject = Bio::AlignIO->new(-file => "FASTA_Files/nuc_post-alignment/aligned_ND2.fas",-format => 'fasta');
# my $ND3_alignobject = Bio::AlignIO->new(-file => "FASTA_Files/nuc_post-alignment/aligned_ND3.fas",-format => 'fasta');
# my $ND4_alignobject = Bio::AlignIO->new(-file => "FASTA_Files/nuc_post-alignment/aligned_ND4.fas",-format => 'fasta');
# my $ND4L_alignobject = Bio::AlignIO->new(-file => "FASTA_Files/nuc_post-alignment/aligned_ND4L.fas",-format => 'fasta');
# my $ND5_alignobject = Bio::AlignIO->new(-file => "FASTA_Files/nuc_post-alignment/aligned_ND5.fas",-format => 'fasta');
# my $ND6_alignobject = Bio::AlignIO->new(-file => "FASTA_Files/nuc_post-alignment/aligned_ND6.fas",-format => 'fasta');
# my $CYTB_alignobject = Bio::AlignIO->new(-file => "FASTA_Files/nuc_post-alignment/aligned_CYTB.fas",-format => 'fasta');
# my $COX1_alignobject = Bio::AlignIO->new(-file => "FASTA_Files/nuc_post-alignment/aligned_COX1.fas",-format => 'fasta');
# my $COX2_alignobject = Bio::AlignIO->new(-file => "FASTA_Files/nuc_post-alignment/aligned_COX2.fas",-format => 'fasta');
# my $COX3_alignobject = Bio::AlignIO->new(-file => "FASTA_Files/nuc_post-alignment/aligned_COX3.fas",-format => 'fasta');

# my $out = Bio::AlignIO->new(-file => ">FASTA_Files/nuc_post-alignment/concatenated.fas",-format=> 'fasta');

# my $aln0 = $ATP6_alignobject->next_aln();
# my $aln1 = $ATP8_alignobject->next_aln();
# my $aln2 = $ND1_alignobject->next_aln();
# my $aln3 = $ND2_alignobject->next_aln();
# my $aln4 = $ND3_alignobject->next_aln();
# my $aln5 = $ND4_alignobject->next_aln();
# my $aln6 = $ND4L_alignobject->next_aln();
# my $aln7 = $ND5_alignobject->next_aln();
# my $aln8 = $ND6_alignobject->next_aln();
# my $aln9 = $CYTB_alignobject->next_aln();
# my $aln10 = $COX1_alignobject->next_aln();
# my $aln11 = $COX2_alignobject->next_aln();
# my $aln12 = $COX3_alignobject->next_aln();

# my $concat_alignment = cat($aln0,$aln1,$aln2,$aln3,$aln4,$aln5,$aln6,$aln7,$aln8,$aln9,$aln10,$aln11,$aln12);

# $out->write_aln($concat_alignment);


#Convert to PhyML for modeltesting
#use seqConverter.pl
# system("perl", "subs/seqConverter.pl", "-dATP6_align.fas", "-ope");
# system("perl", "subs/seqConverter.pl", "-dATP8_align.fas", "-ope");
# system("perl", "subs/seqConverter.pl", "-dCOX1_align.fas", "-ope");
# system("perl", "subs/seqConverter.pl", "-dCOX2_align.fas", "-ope");
# system("perl", "subs/seqConverter.pl", "-dCOX3_align.fas", "-ope");
# system("perl", "subs/seqConverter.pl", "-dCYTB_align.fas", "-ope");
# system("perl", "subs/seqConverter.pl", "-dND1_align.fas", "-ope");
# system("perl", "subs/seqConverter.pl", "-dND2_align.fas", "-ope");
# system("perl", "subs/seqConverter.pl", "-dND3_align.fas", "-ope");
# system("perl", "subs/seqConverter.pl", "-dND4_align.fas", "-ope");
# system("perl", "subs/seqConverter.pl", "-dND4L_align.fas", "-ope");
# system("perl", "subs/seqConverter.pl", "-dND5_align.fas", "-ope");
# system("perl", "subs/seqConverter.pl", "-dND6_align.fas", "-ope");


##REMOVE ALL INVARIANT SITES?


#Model testing
#run model testing for our different partitions of interest
#USE mraic.pl
#need phmyl.exe in the same folder, or phyml.exe in the PATH
#three partitions per protein coding gene + 1 partition per RNA?
#Convert FASTA files to PhyML format using seqConverter.pl
#sleep(15); #make sure files are created first

#system("perl", "mraic.pl", "-modeltest", "ATP6_align.phylip");
#system("perl", "mraic.pl", "-modeltest", "ATP8_align.phylip");
#system("perl", "mraic.pl", "-modeltest", "COX1_align.phylip");
#system("perl", "mraic.pl", "-modeltest", "COX2_align.phylip");
#system("perl", "mraic.pl", "-modeltest", "COX3_align.phylip");
#system("perl", "mraic.pl", "-modeltest", "CYTB_align.phylip");
#system("perl", "mraic.pl", "-modeltest", "ND1_align.phylip");
#system("perl", "mraic.pl", "-modeltest", "ND2_align.phylip");
#system("perl", "mraic.pl", "-modeltest", "ND3_align.phylip");
#system("perl", "mraic.pl", "-modeltest", "ND4_align.phylip");
#system("perl", "mraic.pl", "-modeltest", "ND4L_align.phylip");
#system("perl", "mraic.pl", "-modeltest", "ND5_align.phylip");
#system("perl", "mraic.pl", "-modeltest", "ND6_align.phylip");



#sleep 15 seconds to allow files to be created and closed before starting the model testing
print "brief pause \n";
foreach (1..15) {
	print ".";
	sleep(1);
	}
print "\n";

#change directory to jmodeltest2
chdir "subs/jModelTest2";

##HAVE TO CHANGE THESE TO GIVE JAVA MORE RAM
system("java -jar jModelTest.jar -d ../../FASTA_Files/nuc_post-alignment/aligned_ATP6.fas -g 4 -i -f -AIC -tr 3 > ../../Other_Data/jmodeltest_results/ATP6_jmodeltest_results.txt");
system("java -jar jModelTest.jar -d ../../FASTA_Files/nuc_post-alignment/aligned_ATP8.fas -g 4 -i -f -AIC -tr 3 > ../../Other_Data/jmodeltest_results/ATP8_jmodeltest_results.txt");
system("java -jar jModelTest.jar -d ../../FASTA_Files/nuc_post-alignment/aligned_CYTB.fas -g 4 -i -f -AIC -tr 3 > ../../Other_Data/jmodeltest_results/CYTB_jmodeltest_results.txt");
system("java -jar jModelTest.jar -d ../../FASTA_Files/nuc_post-alignment/aligned_COX1.fas -g 4 -i -f -AIC -tr 3 > ../../Other_Data/jmodeltest_results/COX1_jmodeltest_results.txt");
system("java -jar jModelTest.jar -d ../../FASTA_Files/nuc_post-alignment/aligned_COX2.fas -g 4 -i -f -AIC -tr 3 > ../../Other_Data/jmodeltest_results/COX2_jmodeltest_results.txt");
system("java -jar jModelTest.jar -d ../../FASTA_Files/nuc_post-alignment/aligned_COX3.fas -g 4 -i -f -AIC -tr 3 > ../../Other_Data/jmodeltest_results/COX3_jmodeltest_results.txt");
system("java -jar jModelTest.jar -d ../../FASTA_Files/nuc_post-alignment/aligned_ND1.fas -g 4 -i -f -AIC -tr 3 > ../../Other_Data/jmodeltest_results/ND1_jmodeltest_results.txt");
system("java -jar jModelTest.jar -d ../../FASTA_Files/nuc_post-alignment/aligned_ND2.fas -g 4 -i -f -AIC -tr 3 > ../../Other_Data/jmodeltest_results/ND2_jmodeltest_results.txt");
system("java -jar jModelTest.jar -d ../../FASTA_Files/nuc_post-alignment/aligned_ND3.fas -g 4 -i -f -AIC -tr 3 > ../../Other_Data/jmodeltest_results/ND3_jmodeltest_results.txt");
system("java -jar jModelTest.jar -d ../../FASTA_Files/nuc_post-alignment/aligned_ND4.fas -g 4 -i -f -AIC -tr 3 > ../../Other_Data/jmodeltest_results/ND4_jmodeltest_results.txt");
system("java -jar jModelTest.jar -d ../../FASTA_Files/nuc_post-alignment/aligned_ND4L.fas -g 4 -i -f -AIC -tr 3 > ../../Other_Data/jmodeltest_results/ND4L_jmodeltest_results.txt");
system("java -jar jModelTest.jar -d ../../FASTA_Files/nuc_post-alignment/aligned_ND5.fas -g 4 -i -f -AIC -tr 3 > ../../Other_Data/jmodeltest_results/ND5_jmodeltest_results.txt");
system("java -jar jModelTest.jar -d ../../FASTA_Files/nuc_post-alignment/aligned_ND6.fas -g 4 -i -f -AIC -tr 3 > ../../Other_Data/jmodeltest_results/ND6_jmodeltest_results.txt");
#change working directory back to FGT folder
chdir("../..");


#Parsimony analysis
#fork and run at same time as model testing? or create new instance/pipeline somehow?
#Run parsimony analysis on the big concatenation with TNT



#RAxML likelihood tree
#pull model scores from mraic.pl output file



#prune taxa from tree snd analyze, repeat etc...
#use treePruner.pl