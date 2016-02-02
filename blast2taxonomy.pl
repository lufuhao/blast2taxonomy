#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use Bio::LITE::Taxonomy::NCBI;
use Bio::LITE::Taxonomy::NCBI::Gi2taxid;
use constant USAGE=><<EOH;

SYNOPSIS:

perl $0 --input my.fa [Options]
Version: LUFUHAO20160202

Requirements:
	Perl Modules
		Getopt::Long
		Bio::LITE::Taxonomy::NCBI
		Bio::LITE::Taxonomy::NCBI::Gi2taxid

Descriptions:
	Accept BLASTX tab-delimited output, retrieve taxid, identify seqs
	Not belong to some species/genus/family/order/class/phylum/kingdom,
	and output seqID list which are potential contaminants.

Options:
	--help/-h
		Print this help/usage;
	--input|-i <FILE>
		[necessary] Input fasta file to extract from;
	--col_taxid|-n <INT>
		which column in '--input' is taxid;
	--min_ident|-t <Float>
		Minimum percentage identity in 3rd column of '--input';
		Default: 0
	--min_alnlen|-l <INT>
		Minimum alignment length in 4th column of '--input';
		For protein database, it's length of AA;
		Default: 0
	--max_evalue|-e <Float>
		Maximum E-value in 11th column of '--input';
		Default: 10
	--gi_taxid|-g <File>
		NCBI dictionary file to convert GI to TaxID, NCBI FTP:
		For Protein GI:
		ftp://ftp.ncbi.nih.gov/pub/taxonomy/gi_taxid_prot.dmp.gz
		OR ftp://ftp.ncbi.nih.gov/pub/taxonomy/gi_taxid_prot.zip
		For Nucleotide GI:
		ftp://ftp.ncbi.nih.gov/pub/taxonomy/gi_taxid_nucl.dmp.gz
		OR ftp://ftp.ncbi.nih.gov/pub/taxonomy/gi_taxid_nucl.zip
	--nodesdmp|-d <File>
		NCBI taxonomy nodes.dmp file; NCBI FTP:
		ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz
		OR ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdmp.zip
	--namesdmp|-m <File>
		NCBI taxonomy names.dmp file; NCBI_FTP:
		ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz
		OR ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdmp.zip
	--gi2taxid|-x <File>
		2 column and tab delimited file. GI in 1st column, and taxid in 2nd;
		For self-learning;
	--species_tree|-s <File>
		Species tree file; tab-delimited:
		spcies/Tab/genus/Tab/family/Tab/order/Tab/class/Tab/phylum/Tab/kingdom
		For self-learning
	--output/-o <FILE>
		[Optional] New file for the extracted sequences;
	--verbose
		[Optional] Detailed output for trouble-shooting;
	--version/-v
		[Optional] Print current SCRIPT version;
		

Example:
	perl $0 --i myBlast.oufmt8.fasta --min_ident 80.00 --min_alnlen 50 \
		--gi_taxid Path/to/gi_taxid_prot.dmp \
		--nodesdmp Path/to/nodes.dmp
		--namesdmp Path/to/names.dmp
		--output output.txt
		

Author:
	Fu-Hao Lu
	Post-Doctoral Scientist in Micheal Bevan laboratory
	Cell and Developmental Department, John Innes Centre
	Norwich NR4 7UH, United Kingdom
	E-mail: Fu-Hao.Lu\@jic.ac.uk
EOH
###HELP ends##########################
die USAGE if (scalar(@ARGV)==0);


###Receving parameter#################
our ($help, $input, $col_taxid, $min_ident, $min_aln_len, $max_evalue, $gi_taxid, $nodesdmp, $namesdmp, $gi2taxid_file, $species_tree, $output, $seqtogo, $verbose, $ver);
GetOptions(
	"help|h!" => \$help,
	"input|i=s" => \$input,
	"col_taxid|n:i" => \$col_taxid,
	"min_ident|t:f" => \$min_ident,
	"min_alnlen|l:i" => \$min_aln_len,
	"max_evalue|e:f" => \$max_evalue,
	"gi_taxid|g:s" => \$gi_taxid,
	"nodesdmp|d=s" => \$nodesdmp,
	"namesdmp|m=s" => \$namesdmp,
	"gi2taxid|x:s" => \$gi2taxid_file,
	"species_tree|s:s" => \$species_tree,
	"output|o:s" => \$output,
	"seqtogo|k:s" =>\$seqtogo, 
	"verbose!" => \$verbose,
	"version|v!" => \$ver) or die USAGE;
#die USAGE unless (@ARGV);
($help or $ver) and die USAGE;



###Default and initalization########################
$col_taxid=0 unless (defined $col_taxid);
$min_ident=0 unless (defined $min_ident);
$min_aln_len=0 unless (defined $min_aln_len);
$max_evalue=10 unless (defined $max_evalue);
$species_tree='Temp.species.tree' unless (defined $species_tree);
$gi2taxid_file='Temp.GI2TaxID' unless (defined $gi2taxid_file);
our %gi2taxid=();
our %taxid2taxonomy=();
our %taxid_tax=();
our %GitoKeep=();
our %GiToRemove=();
our ($KPspcies, $KPgenus, $KPfamily, $KPorder, $KPclass, $KPphylum, $KPkingdom)=('Triticum aestivum', 'Triticum', 'Poaceae', 'Poales', 'Liliopsida', 'Streptophyta', 'Viridiplantae');




###File and Forlders###############################
###File check
die "Can not find the --input file specified\n" unless (defined $input and (-e $input));
die "Can not find --gi_taxid file specified by --gi_taxid. Plase double check the path and filename\n" unless (defined $gi_taxid and (-e $gi_taxid));
die "Can not find --namesdmp file specified by --namesdmp. Plase double check the path and filename\n" unless (defined $namesdmp and (-e $namesdmp));
die "Can not find --nodesdmp file specified by --nodesdmp. Plase double check the path and filename\n" unless (defined $nodesdmp and (-e $nodesdmp));
###input
our ($input_basename, $input_base_no_ext); 
$input_basename=&RetrvBase($input);
$input_base_no_ext=&RetrvNoExt($input);
#define output
$output=$input_base_no_ext.".taxID" unless (defined $output);
&FileExistsTest($output);
$seqtogo=$input_base_no_ext.'.SeqToRemove' unless (defined $seqtogo);
&FileExistsTest($seqtogo);



###Initialize the DB
#ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/gi_taxid_nucl.dmp.gz
#ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/gi_taxid_prot.dmp.gz
our $gi_taxid_dir=&RetrvDir($gi_taxid);
our $gi_taxid_basenoext=&RetrvNoExt($gi_taxid);
our $gi_taxid_bin=$gi_taxid_basenoext.'.bin';
use Bio::LITE::Taxonomy::NCBI::Gi2taxid qw/new_dict/;
new_dict (in => "$gi_taxid", out => "$gi_taxid_dir/$gi_taxid_bin") unless ($col_taxid>0 or (-e "$gi_taxid_dir/$gi_taxid_bin"));
our $dictbin='';
our $taxdb='';
###Print all input parameters in --verbose
print "The running parameters are:\nInput: $input\nInput_basename: $input_basename\nInput_base_no_ext: $input_base_no_ext\nOutput: $output\nCol_taxid: $col_taxid\nMin_ident: $min_ident\nMin_aln_len $min_aln_len\nMax_evalue: $max_evalue\ngi_taxid: $gi_taxid\nNodesdmp: $nodesdmp\nNamesdmp: $namesdmp\nspecies_tree: $species_tree\nseqtogo: $seqtogo\n" if (defined $verbose);



###load --gi2taxid into %gi2taxid if any
if (-s $gi2taxid_file) {
	open (GI2TAXIDFILE, "$gi2taxid_file") || die "Can not write to the file specified by '--gi2taxid'\n";
	my $G2T_test=0;
	while (my $G2T_line=<GI2TAXIDFILE>) {
		chomp $G2T_line;
		my @G2T_arr=(); @G2T_arr=split(/\t/, $G2T_line);
		if ((defined $G2T_arr[0]) and (scalar(@G2T_arr) == 2) and ($G2T_arr[0] ne '')) {
			$gi2taxid{$G2T_arr[0]}=$G2T_arr[1];
			$G2T_test=1;
		}
		elsif ($G2T_test==0) {
			&backup($gi2taxid_file);
			last;
		}
	}
	close GI2TAXIDFILE;
}




###Load species tree if provided into %taxid_tax
if (-s $species_tree) {
	open (SPECIESTREE, "$species_tree") || die "Can not open $species_tree specified by --species_tree";
	my $ST_test=0;
	while (my $ST_line=<SPECIESTREE>) {
		chomp $ST_line;
		my @ST_arr=();
		@ST_arr=split(/\t/, $ST_line);
		if ((defined $ST_arr[0]) and (scalar(@ST_arr) == 8) and ($ST_arr[0] ne '')) {
			$taxid_tax{$ST_arr[0]}=[($ST_arr[1], $ST_arr[2], $ST_arr[3], $ST_arr[4], $ST_arr[5], $ST_arr[6], $ST_arr[7])];
			$ST_test=1;
		}
		elsif ($ST_test == 0) {
			&backup($species_tree);
			last;
		}
	}
	close SPECIESTREE;
}



###BLASTreading####################################
open (INPUT, "$input") || die "Can not open $input file\n";
open (OUTPUT, ">>$output") || die "Can not write to $output file\n";
open (GI2TAXIDFILE, ">>$gi2taxid_file") || die "Can not write to the file specified by '--gi2taxid'\n";
open (SPECIESTREE, ">>Temp.species.tree") || die "Can not write to Temporary species tree file\n";
<INPUT>;<INPUT>;<INPUT>;<INPUT>;<INPUT>;
while (my $blast_line = <INPUT>) {
	chomp $blast_line;
	my @blast_arr=split(/\t/, $blast_line);
	my $gi_num=''; 
	my $taxid='';
	our @species_arr=(); 
	unless (exists $GitoKeep{$blast_arr[0]}) {
		if ($blast_arr[2]>=$min_ident and $blast_arr[3]>=$min_aln_len and $blast_arr[10]<=$max_evalue) {
###convert GI to taxID
			($gi_num=$blast_arr[1])=~s/^gi\|(\d+)\|.*$/$1/;
			if ($col_taxid != 0) {
				$taxid=$blast_arr[$col_taxid-1];
			}
			elsif (exists $gi2taxid{$gi_num}) {
				$taxid=$gi2taxid{$gi_num};
			}
			else {
				$taxid=&GI2TaxID($gi_num);
				$gi2taxid{$gi_num}=$taxid;
				print GI2TAXIDFILE $gi_num."\t".$taxid."\n";
			}
###TaxID to Taxonomy
			unless (exists $taxid_tax{$taxid}){
				&TaxID2Taxa($taxid);
				print SPECIESTREE $taxid."\t".join("\t", @{$taxid_tax{$taxid}})."\n";
			}
			@species_arr=@{$taxid_tax{$taxid}};
			my $species_line=join("\t", @species_arr); 
			print OUTPUT $blast_line."\t".$gi_num."\t".$species_line."\n";
			print $gi_num."\t".$taxid."\t"."@species_arr\n" if (defined $verbose);
			if (($species_arr[0] ne ('none'||'')) and (($species_arr[0] =~ /Oryza/) or ($species_arr[2] eq "$KPfamily") or ($species_arr[3] eq "$KPorder") or ($species_arr[4] eq "$KPclass") or ($species_arr[5] eq "$KPphylum"))) {
				$GitoKeep{$blast_arr[0]}++;
			}
			elsif (! (exists $GitoKeep{$blast_arr[0]})) {
				$GiToRemove{$blast_arr[0]}++;
			}
		}
	}
}
close SPECIESTREE;
close GI2TAXIDFILE;
close INPUT;
close OUTPUT;

###Write SeqIDs to be removed to FILE
our @GiToRemove_arr=keys %GiToRemove;
open (SEQKEEP, ">>$seqtogo") || die "Can not write the SeqIDs kept to the file specified by --seqtokeep\n";
foreach my $GiToRemove_ind (@GiToRemove_arr) {
	unless ($GitoKeep{$GiToRemove_ind}) {
		print SEQKEEP $GiToRemove_ind."\n";
	}
}
close SEQKEEP;



###############################################
###sub functions###############################
###############################################

### Backup file
sub backup {
	my $BUori_name=shift @_;
	my $BU_new_name='';
	$BU_new_name=$BUori_name.".bak.".int(rand(10000));
	if (-e "$BU_new_name") {
		&backup($BUori_name);
	}
	rename($BUori_name, $BU_new_name);
	print "The Original file $BUori_name was renamed to $BU_new_name\n";
}

###Return file path
sub RetrvDir {
	my $RD_ori=shift @_;
	chomp $RD_ori;
	my $RD_new='';
	if ($RD_ori=~/\//) {
		($RD_new=$RD_ori)=~ s/(.*)\/.*$/$1/s;
	}
	else {
		$RD_new='.';
	}
	return ($RD_new);
	
	
}

###return file base name
sub RetrvBase {
	my $RB_ori=shift @_;
	chomp $RB_ori;
	my $RB_new='';
	($RB_new=$RB_ori)=~ s/.*\///s; 
	return $RB_new;
}

###Return file name without extension
sub RetrvNoExt {
	my $RNE_ori=shift @_;
	chomp $RNE_ori;
	my ($RNE_new, $RNE_base)=('', '');
	($RNE_base=$RNE_ori)=~ s/.*\///s; 
	($RNE_new=$RNE_base)=~s/^(.*)\.\w+$/$1/;
	return $RNE_new;
}

###determine how to process those files existed
sub FileExistsTest {
	my $FET_fileid=shift @_;
	chomp $FET_fileid;
	if (-e $FET_fileid) {
	if ($verbose) {
		print "The $FET_fileid file specified already existed\nNeed to overwrite [y] or backup [n] or others [AnyKey], default[y]:\n";
		my $FET_test=''; $FET_test=<STDIN>;
		chomp $FET_test;
		if (lc ($FET_test) eq ('y'||'yes')) {
			unlink($FET_fileid);
		}
		elsif (lc ($FET_test) eq ('n'||'no')) {
			&backup($FET_fileid);
		}
		else {
			die "Please specify a new name for output\n";
		}
	}
	else {
		unlink($FET_fileid);
	}
	}
}

###Retrieve taxid by NCBI GI number
sub GI2TaxID{
	my $G2T_gi=shift @_;
	my $G2T_taxid='';
	use Bio::LITE::Taxonomy::NCBI::Gi2taxid;
	$dictbin= Bio::LITE::Taxonomy::NCBI::Gi2taxid->new(dict=>"$gi_taxid_dir/$gi_taxid_bin");
	$G2T_taxid = $dictbin->get_taxid($G2T_gi);
	return $G2T_taxid;
}

###Retrieve taxonomy info by TaxID
sub TaxID2Taxa {
	my $TTtaxid=shift @_;
	use Bio::LITE::Taxonomy::NCBI;
	$taxdb=Bio::LITE::Taxonomy::NCBI->new (
		db=>'NCBI',
		names=>"$namesdmp",
		nodes=>"$nodesdmp",
	);
	my ($TTspecies, $TTgenus, $TTfamily, $TTorder, $TTclass, $TTphylum, $TTkingdom)=('', '', '', '', '', '', ''); 
	$TTspecies=$taxdb->get_term_at_level($TTtaxid, "species");
	$TTspecies='none' unless ((defined $TTspecies) and ($TTspecies ne''));
	chomp $TTspecies;
	$TTgenus=$taxdb->get_term_at_level($TTtaxid, "genus");
	chomp $TTgenus;
	$TTgenus='none' unless ((defined $TTgenus) and ($TTgenus ne ''));
	$TTfamily=$taxdb->get_term_at_level($TTtaxid, "family");
	chomp $TTfamily;
	$TTfamily='none' unless ((defined $TTfamily) and ($TTfamily ne ''));
	$TTorder=$taxdb->get_term_at_level($TTtaxid, "order");
	chomp $TTorder;
	$TTorder='none' unless ((defined $TTorder) and ($TTorder ne ''));
	$TTclass=$taxdb->get_term_at_level($TTtaxid, "class");
	chomp $TTclass;
	$TTclass='none' unless ((defined $TTclass) and ($TTclass ne ''));
	$TTphylum=$taxdb->get_term_at_level($TTtaxid, "phylum");
	chomp $TTphylum;
	$TTphylum='none' unless ((defined $TTphylum) and ($TTphylum ne ''));
	$TTkingdom=$taxdb->get_term_at_level($TTtaxid, "kingdom");
	chomp $TTkingdom;
	$TTkingdom='none' unless ((defined $TTkingdom) and ($TTkingdom ne ''));
	$taxid_tax{$TTtaxid}=[($TTspecies, $TTgenus, $TTfamily, $TTorder, $TTclass, $TTphylum, $TTkingdom)];
	print "subfunction01: ".$TTspecies."\t".$TTgenus."\t".$TTfamily."\t".$TTorder."\t".$TTclass."\t".$TTphylum."\t".$TTkingdom."\n" if (defined $verbose); ###For test
}
