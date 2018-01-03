#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Try::Tiny;
use Bio::TreeIO;
use JSON::Parse 'json_file_to_perl';
use IO::String;
use List::Util qw(min);

my ($pantherdir,$pantherhmm,$fastafile,$outfile,$annotationfile,$raxmlloc,$hmmerloc,$keep, 
    $algo, $cpus,  $auto, $help);

my $directory;
$Getopt::Long::ignorecase=0;

#Todos:
#Make the hmmer location agnostic to whether you are running hmmscan or hmmsearch
#Step the number of cpus
#Auto detect whether hmmscan or hmmsearch should be used
#Look for the press files if running hmmscan


&GetOptions
    (
     "f=s" => \$fastafile, # -f for the input fasta file
     "r=s" => \$raxmlloc, # -r for the location of RAXML if not in PATH
     "s=s" => \$hmmerloc, # -s for the location of hmmscan if not in PATH
     "o=s" => \$outfile, # -o for the output file
     "d=s" => \$directory, # -d for directory of the package
     "algo=s" => \$algo,
     "auto" => \$auto,
     "cpus=i" => \$cpus,
     "k=s" => \$keep, #-k for keeping the tmp files
     "h"   => \$help, #print the usage statement
    ) or die "Invalid option passed.\n";


&usage if($help);

#Make sure we know where we are going to write to.
&usage ("Please specify output file\n") unless ($outfile);

#Check the fasta file is defined, present and has size 
&usage ("Please specify input fasta file\n") unless ($fastafile);
if(!-e $fastafile){
  die "Your fasta file, $fastafile, does not exisit.\n";
}elsif(! -s $fastafile){
  die "Your fasta file, $fastafile, has not size\n";
}

#Check that the PANTHER directory is present, and contains the HMM files
&usage ("Please specify the directory of the pipeline\n") unless ($directory);
if(!-d $directory){
  die "Your directory, $directory, does not exisit.\n";
}

#We expect this directory to have a certain structure
#TODO:Deal with different releases, hardcoded to 12.0
$pantherhmm = "$directory/PANTHER12_famhmm/PANTHER12.0_all_fam.hmm";
if(!-e $pantherhmm){
  die "The PANTHER hmm file, $pantherhmm, does not exist.\n";
}elsif( !-s $pantherhmm ){
  die "The PANTHER hmm file, $pantherhmm,  has no size.\n";
}

#Check that the directory containing the PANTHER alignments is present.
$pantherdir = "$directory/PANTHER12_Tree_MSF";

if(!-d $pantherdir){
  die "The PANTHER alignments directory, $pantherdir, does not exisit.\n";
}
#-------------------------------------------------------------------------------------
$annotationfile = "$directory/PANTHER12_PAINT_Annotations/PANTHER12_PAINT_Annotatations_TOTAL.txt";

if(!-e $annotationfile){
  die "The PANTHER annotation file, $annotationfile, does not exist.\n";
}elsif( !-s $annotationfile ){
  die "The PANTHER hmm file, $annotationfile,  has no size.\n";
}

#Now read in the annotations file to work out which PANTHER entries have annotations.
my %annotations; my %pthrs;
open ANO, "< $annotationfile" or die "cannot open $annotationfile\n";
while(<ANO>){
  chomp;
  my @array = split(/\t/);
  $annotations{$array[0]} = $array[1];
  my @a = split(/:/,$array[0]);
  $pthrs{$a[0]} =1; # store the PANTHER families with annotations
}
close ANO;

#-------------------------------------------------------------------------------------

if($algo and $auto){
  warn "Please specify either -algo <hmmscan|hmmsearch> or -auto, but not both.\n";
  warn "If you are unsure, let the script decide the most efficient approach.\n"; 
}

if(defined($algo)){
  #Check that the algorithm is either hmmscan or hmmsearch
  if($algo !~ /^(hmmscan|hmmsearch)/){
    print "Unknown algorithm $algo\n";
  }
}elsif(!defined($auto)){
  $auto = 1;
}

#If we need to, auto detect the file
if(defined($auto)){
  #
  $algo = autodetect($fastafile, $pantherhmm);
}

#Need to check for the presence of the hmmpress files
if($algo eq "hmmscan"){
  
}

#In the case of hmmer cpu=0 switches off threading. 
#Otherwise hmmer will autodetect the number of CPUs
#and use all of them, which does not work well on HPCs.

if(!defined($cpus)){
  $cpus = 0;
}

open (FINALOUT, '>',  $outfile) or die "cannot output to $outfile\n";

#-------------------------------------------------------------------------------------
#my $hmmscanout = "$fastafile.hmmscanout.$$"; #TODO replace this.
my $hmmerout = "$fastafile.hmmerout"; 

unless (-s $hmmerout){
  my $hmmercommand;
  if ($hmmerloc) {
    $hmmercommand = "$hmmerloc/$algo --notextw --cpu $cpus -o $hmmerout $pantherhmm $fastafile > /dev/null";
  } else {
    $hmmercommand = "$algo --notextw --cpu $cpus -o $hmmerout $pantherhmm $fastafile > /dev/null";
  }
  system($hmmercommand) and die "Error running $hmmercommand";
}


#-------------------------------------------------------------------------------------
#Now parse the hmmer output


my $queryid;my $matchpthr; my @matchalign; my @hmmstart ; my @hmmend; my @codes; my @hmmalign;
my $hmmline =0; my $domainline =0; my $alignmentline =0;
my $matchn =0;
open HMM, "< $hmmerout" or die "cannot open $hmmerout\n";
while(<HMM>){
  $hmmline++;
  if ($_ =~ /^Query:\ +([^ ]+) /){
    $queryid = $1;
  }
  elsif ($_ =~ />> (PTHR[0-9]+)/){ # use the first match
    $matchn++;
    if ($matchn==1){$matchpthr = $1;}
    next;
  }
  if ($matchn ==1){
    if ($hmmline >= $domainline+4){
      if ($_ =~ /!/){
	my @a = split(/ +/);
	push(@hmmstart,$a[7]);
	push(@hmmend,$a[8]);
      }
    }
    if ($_ =~ /== domain/){
      $alignmentline =$hmmline;
    }
    if ($hmmline == $alignmentline+1){
      my @a = split(/ +/,$_);
      push(@hmmalign,$a[3]);
    }
    elsif ($hmmline == $alignmentline+3){
      my @a = split(/ +/,$_);
      push(@matchalign,$a[3]);
    }
    elsif ($hmmline == $alignmentline+4){
      my @a = split(/ +/,$_);
      push(@codes,$a[1]);
    }
  }
  if ($_ =~ /^\/\/$/){
    &pipline($queryid,$matchpthr,\@matchalign,\@hmmstart,\@hmmend,\@codes,\@hmmalign);
    ($queryid,$matchpthr,@matchalign,@hmmstart,@hmmend,@codes,@hmmalign,$matchn,$domainline,$alignmentline) = ();
  }
}
close HMM;
close FINALOUT;

sub pipline{
  my ($queryid,$matchpthr,$matchalign,$hmmstart,$hmmend,$codes,$hmmalign) = @_;
  print "$queryid,$matchpthr,$matchalign,$hmmstart,$hmmend,$codes\n";
#  return;
####################################################################
# add the query sequence to precalculated PANTHER msf using hmmalign
###################################################################

  unless (exists $pthrs{$matchpthr}){
    print "$queryid\t$matchpthr\tAnnotationsNotAvailable\n";
#    return;
  }
  $queryid =~ s/[^\w]/\_/g;
  my $hmmalignmsf = "$pantherdir/$matchpthr.AN.fasta";
  my $queryfasta = "$directory/tmp/$queryid.$matchpthr.fasta";
  
  my $lines; my $lastline;
  open IN, "< $hmmalignmsf" or die;
  my $first = <IN>;
  while(<IN>){
    chomp;
    if ($_ =~ /^>/){
      last;
    }
    else{
      $lines++;
      $lastline = $_;
    }
  }
  close IN;
  my $length_msf = ($lines-1)*80 + length($lastline);
  my $querymsf;

  my $domains = scalar @$matchalign;
  my $start = $hmmstart->[0];
  my $align = $matchalign->[0];
  my $halign = $hmmalign->[0];
  my $code = $codes->[0];
  foreach my $i (1..$start-1){
    $querymsf .= "-";
  }
  my @aligna = split(//,$align);
  my @codea = split(//,$code);
  my @haligna = split(//,$halign);
  my $alignas = length($code);
  foreach my $a (0..$alignas-1){
    my $aa = $aligna[$a]; 
    my $c = $codea[$a];
    my $haa = $haligna[$a];

    if (($aa =~ /[A-Z]/) and ($haa eq ".")){
      next;
    }
    $querymsf .= $aa;
  }
  foreach my $j (1..$domains-1){ 
    my $start = $hmmstart->[$j];
    my $align = $matchalign->[$j];
    my $end = $hmmend->[$j-1];
    foreach my $i ($end+1..$start-1){
      $querymsf .= "-";
    }
    my $halign = $hmmalign->[$j];
    my $code = $codes->[$j];
    my @aligna = split(//,$align);
    my @codea = split(//,$code);
    my @haligna = split(//,$halign);
    my $alignas = length($code);
    foreach my $a (0..$alignas-1){
      my $aa = $aligna[$a]; 
      my $c = $codea[$a];
      my $haa = $haligna[$a];

      if (($aa =~ /[A-Z]/) and ($haa eq ".")){
	next;
      }
      $querymsf .= $aa;
    }
  }
  my $prevend = $hmmend->[$domains-1];
  foreach my $i ($prevend..$length_msf-1){
    $querymsf .= "-";
  }
  $querymsf = uc($querymsf);
  my $l = length($querymsf);

  my @parts = $querymsf =~ /(.{1,80})/g;
  open OUT,"> $queryfasta" or die "cannot open $queryfasta\n";
 # print  ">query_$queryid\n";
  print OUT ">query_$queryid\n";
  foreach my $line (@parts){
    print OUT "$line\n";
  #  print "$line\n";
  }
  close OUT;
  
  unless ($l eq $length_msf){
    print "ERROR MSF of $queryid should have length $length_msf, actual length is $l\n";
    return;
  }

  `cat $hmmalignmsf >> $queryfasta`;
  
  ###############################################################
  ###run RAxML ################################################## 
  ###############################################################
  my $raxmldir = "$directory/tmp/$matchpthr"."_$queryid"."_"."raxml$$";
  my $bifurnewick = "$pantherdir/$matchpthr.bifurcate.newick";
  unless (-e $bifurnewick) {print STDERR "no bifurcate newickfile for $matchpthr\n";return;};

  mkdir($raxmldir); my $raxmlcommand;
  if ($raxmlloc){
  $raxmlcommand = "$raxmlloc -f y -p 12345 -t $bifurnewick -G 0.05 -m PROTGAMMAWAG  -s $queryfasta -n $matchpthr -w $raxmldir "; }
  else{ $raxmlcommand = "raxmlHPC-SSE3 -f y -p 12345 -t $bifurnewick -G 0.05 -m PROTGAMMAWAG  -s $queryfasta -n $matchpthr -w $raxmldir ";}
  try{ system($raxmlcommand); } catch {
    warn "caught error for $queryid $matchpthr: $_";
  };

  ## now find the location of the graft sequence in the tree 
  my $mapANs;
  try { $mapANs = &mapto($raxmldir,$matchpthr,"query_$queryid");} catch{ warn "caught error in mapto for $queryid $matchpthr $_";};
  unless ($mapANs){
    $mapANs = "root";
    my $annotation = $annotations{$matchpthr.":".$mapANs};
    print FINALOUT "$queryid\t$matchpthr\t$annotation\n";
    return;
  }
  my $commonAN = &commonancestor($mapANs,$matchpthr);
  unless ($commonAN) {$commonAN = "root"};
  my $annotation = $annotations{$matchpthr.":".$commonAN};
  print FINALOUT "$queryid\t$matchpthr\t$annotation\n";
  unless($keep){
    my $command = "rm -rf $directory/tmp/*";
    system($command);
  }
}


############find common ancestor of a list of ANs which are leaf genes#########

sub commonancestor{
  my $mapANs = shift;
  my $matchpthr = shift;
  my @a = split(/;/,$mapANs);
  unless ($a[1]){
    return $mapANs;
  }
  my $newick = $pantherdir."/$matchpthr.newick";
  my $treeio = Bio::TreeIO->new(-file => $newick, -format => "newick");
  my $tree = $treeio->next_tree;
  my $size = scalar @a;
  my @anceslist; my %anstore;
  my $current = $tree->find_node(-id=>$a[0]);

  while(1){
    my $mom = $current->ancestor||last;
    unless ($mom->id) {last;}
    push(@anceslist,$mom->id);
    $current = $mom;
  }

  foreach my $i (1..$size-1){
    my $current = $tree->find_node(-id=>$a[$i]);
    while(1){
      my $mom = $current->ancestor||last;
      unless ($mom->id) {last;}
      $anstore{$mom->id}++;
      $current = $mom;
    }
  }
  foreach my $mom (@anceslist){
    if ($anstore{$mom} == $size-1){
      return $mom;
    }
  }
  print STDERR "couldn't find common ancestor for $mapANs $matchpthr??\n";
}

  # this is a subroutine to grasp the location of the query sequence in the raxml output
sub mapto{

  # I found it is very convient to use the jplace output file
  my $dir = shift;
  my $raxml_name = shift;
  my $test =shift;
  
  my $classification;
  $classification = $dir."/RAxML_portableTree.$raxml_name.jplace";

  unless (-e $classification){
    print STDERR "raxml didn't run successfully for $raxml_name, $test\n";
    return;  }
  
  my $jsonperl; my $success;
  try{
    $jsonperl = json_file_to_perl ($classification);
  }  catch { #so
    $success = $_;
    warn "caught error: $_"; # not $@ 
  };
    
  if ($success){
    return;  }
  
  my $treestring = $jsonperl->{"tree"};
  my $locationsref = $jsonperl->{"placements"}->[0]->{"p"};
  my $nloc = scalar @$locationsref; # number of possible locations

  my $allchildren;
  # then modify the treestring

  $treestring =~ s/:[0-9\.]+\{([0-9]+)\}/R$1/g;

  my %AN_label;

  while($treestring =~ /(AN[0-9]+)(R[0-9]+)/g){
    $AN_label{$1} = $2;
    $AN_label{$2} = $1;
    $allchildren .= ";".$1;
  }

  $treestring =~ s/AN[0-9]+//g;
  
  # modification of the tree string to simple newick format 
  my $io = IO::String->new($treestring);

  my $input = Bio::TreeIO->new(-fh => $io, -format => "newick");
  my $rtree = $input->next_tree;

  # find the nodes, and find their common ancestor

  if ($nloc ==1){
    my $map = $locationsref->[0]->[0];
    unless ($map =~ /[0-9]/){
      print STDERR "No RAXML mapping result for longid\n?";
      return $allchildren;
    }
    $map = "R$map";
    if (exists $AN_label{$map}){ # so if map to a leaf node
      return $AN_label{$map};
    }
    else{ # so, map to internal node. IN this case, find the node, and return all the descendants
      my @node = $rtree->find_node(-id => $map);
      my $childrenids;
      foreach my $child ($node[0]->get_all_Descendents){
	next unless ($child->is_Leaf);
	my $cid = $child->id;
	my $can = $AN_label{$cid};
	if ($can){
	  $childrenids .= $can.";";
	}
	else{
	  print "NO AN number for $cid for longid?\n"; 
	}
      }
      return $childrenids;
    }
  }
  elsif ($nloc > 1){
    my %ancestors; my @ancestororder;
    my $indicator;
    foreach my $i (0..$nloc-1){
      if ($i == $nloc-1){
	$indicator = 1;
      }
      my $map = $locationsref->[$i]->[0];
      $map = "R$map";
      my $node = $rtree->find_node(-id => $map);
      my $current = $node;
      unless ($node->is_Leaf){
	$ancestors{$map} ++;
      }
      while(1){
	my $ancestor = $current->ancestor || last;
	my $anid = $ancestor->id;
	if ($anid){
	  $ancestors{$anid}++;
	  if ($indicator){
	    push(@ancestororder,$anid);
	  }
	}
	else{
	  last;
	}
	$current = $ancestor;
      }      
    }
    my $commonancestorid;
    foreach my $anid (@ancestororder){
      if ($ancestors{$anid} >= $nloc){
	$commonancestorid = $anid;
	last;
      }
    }
    unless ($commonancestorid){
      print STDERR "cannot find the common ancestor for multiple placements for longid\n";
      print STDERR Dumper($locationsref);
      return $allchildren;
    }

    my $node = $rtree->find_node(-id => $commonancestorid);
    my $childrenids;
    foreach my $child ($node->get_all_Descendents){
      next unless ($child->is_Leaf);
      my $cid = $child->id;
      my $can = $AN_label{$cid};
      if ($can){
	$childrenids .= $can.";";
      }
      else{
	print "NO AN number for $cid for longid?\n";
      }
    }
    return $childrenids;
  }
}


sub usage{
  my $errorm=shift;
  if ($errorm){
    print "\nError: $errorm\n";
  }

  print qq(
tree grafting pipeline in Perl. This program grafts input sequences in fasta format to best positions in best-matching PANTHER trees

  Where args are:
  -d for directory of the package
  -f for the input fasta file
  -r for the location of RAXML if not in PATH
  -s for the location of hmmscan if not in PATH
  -o for the output file
  -a for the annotation file
  -k for keeping the tmp files
);
  exit;
}


