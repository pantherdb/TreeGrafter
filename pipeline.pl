#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Try::Tiny;
use Bio::TreeIO;
use JSON::Parse 'json_file_to_perl';
use IO::String;
use List::Util qw(min);

my ($pantherdir,$pantherhmm,$fastafile,$outfile,$annotationfile,$raxmlloc,$hmmscanloc,$keep, $help);

$Getopt::Long::ignorecase=0;

&GetOptions
    (
     "p=s" => \$pantherdir, # -p for the PANTHER12 directory of trees and msfs
     "h=s" => \$pantherhmm, # -h for the PANTHER12 hmm directory with family only hmms
     "f=s" => \$fastafile, # -f for the input fasta file
     "r=s" => \$raxmlloc, # -r for the location of RAXML if not in PATH
     "s=s" => \$hmmscanloc, # -s for the location of hmmscan if not in PATH
     "o=s" => \$outfile, # -o for the output file
     "a=s" => \$annotationfile, # -a for the annotation file
     "k=s" => \$keep, #-k for keeping the tmp files
     "h"   => \$help, #print the usage statement
    ) or die "Invalid option passed.\n";

&usage if($help);
&usage ("Please specify input fasta file\n") unless ($fastafile);
&usage ("Please specify output file\n") unless ($outfile);
&usage ("Please specify PANTHER12 hmm directory\n") unless ($pantherhmm);
&usage ("Please specify PANTHER12 directory of trees and msfs\n") unless ($pantherdir);
&usage ("Please specify PANTHER12 GO annotation summary file\n") unless ($annotationfile);

my $hmmscanout = "$fastafile.hmmscanout.$$";

my $hmmscancommand;
if ($hmmscanloc){
  $hmmscancommand = "$hmmscanloc --notextw -o $hmmscanout $pantherhmm/PANTHER12.0_all_fam.hmm $fastafile > /dev/null";}
else{
  $hmmscancommand ="hmmscan --notextw -o $hmmscanout $pantherhmm/PANTHER12.0_all_fam.hmm $fastafile > /dev/null";}

system($hmmscancommand) and die "Error running $hmmscancommand";

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

my $queryid;my $matchpthr;my $matchalign; my $hmmstart;my $hmmend;
open HMM, "< $hmmscanout" or die;
while(<HMM>){
  if ($_ =~ /^Query:\ +([^ ]+) /){
    $queryid = $1;
  }
  elsif ($_ =~ />> (PTHR[0-9]+)/){
    $matchpthr = $1;
  }
  elsif (($matchpthr) and ($_ =~ /$matchpthr[\w\.]+\ +([0-9]+) [\w\.]+ ([0-9]+)/)){
    $hmmstart = $1;$hmmend = $2;
  }
  elsif ($queryid){
    my @a = split(/ +/,$_);
    if( ($a[1] eq $queryid)){
      $matchalign = $a[3];
      &pipline($queryid,$matchpthr,$matchalign,$hmmstart,$hmmend);
      ($queryid,$matchpthr,$matchalign,$hmmstart,$hmmend) = ();
    }
  }
}
close HMM;

sub pipline{
  my ($queryid,$matchpthr,$matchalign,$hmmstart,$hmmend) = @_;
#  print "$queryid,$matchpthr,$matchalign,$hmmstart,$hmmend\n";
####################################################################
# add the query sequence to precalculated PANTHER msf using hmmalign
###################################################################

  unless (exists $pthrs{$matchpthr}){
    print "$queryid\t$matchpthr\tAnnotationsNotAvailable\n";
    return;
  }
  my $hmmalignmsf = "$pantherdir/$matchpthr.AN.fasta";
  my $queryfasta = "./tmp/$queryid.$matchpthr.fasta";
  
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
  my $length_msf = ($lines-1)*80 + lengh($lastline);
  my $querymsf;
  foreach my $i (1..$hmmstart-1){
    $querymsf .= "-";
  }
  $querymsf .= $matchalign;
  foreach my $i ($hmmend+1..$length_msf){
    $querymsf .= "-";
  }
  my @parts = $querymsf =~ /(.{1,80})/g;
  open OUT,"> $queryfasta" or die;
  print OUT ">query_$queryid";
  foreach my $line (@parts){
    print OUT "$line\n";
  }
  close OUT;
  `cat $hmmalignmsf >> $queryfasta`;
  
  ###############################################################
  ###run RAxML ################################################## 
  ###############################################################
  my $raxmldir = "./tmp/$matchpthr"."_$queryid"."_"."raxml$$";
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
    print "$queryid\t$matchpthr\t$annotation\n";
    return;
  }
  my $commonAN = &commonancestor($mapANs,$matchpthr);
  unless ($commonAN) {$commonAN = "root"};
  my $annotation = $annotations{$matchpthr.":".$commonAN};
  print "$queryid\t$matchpthr\t$annotation\n";
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
  my $current = $a[0];
  while(1){
    my $mom = $current->ancestor||last;
    unless ($mom->id) {last;}
    push(@anceslist,$mom->id);
    $current = $mom;
  }

  foreach my $i (1..$size-1){
    my $current = $a[$i];
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
  -p for the PANTHER12 directory of trees and msfs
  -h for the PANTHER12 hmm directory with family only hmms
  -f for the input fasta file
  -r for the location of RAXML if not in PATH
  -s for the location of hmmscan if not in PATH
  -o for the output file
  -a for the annotation file
  -k for keeping the tmp files
);
  exit;
}


