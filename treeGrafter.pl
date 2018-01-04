#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use IO::String;
use Try::Tiny;
use Bio::TreeIO;
use JSON::Parse 'json_file_to_perl';
use Cwd 'abs_path';

$Getopt::Long::ignorecase=0;

my ($fastafile, 
    $outfile, 
    $raxmlloc, 
    $hmmerloc, 
    $directory, 
    $keep, 
    $algo, 
    $cpus,  
    $auto,
    $hmmer,
    $help);

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
     "hmmer=s" => \$hmmer,
     "k" => \$keep, #-k for keeping the tmp files
     "h"   => \$help, #print the usage statement
    ) or die "Invalid option passed.\n";


my $options = {};
processOptions( $options, $help, $outfile, $directory, $fastafile, 
               $raxmlloc, $hmmerloc, $algo, $auto, $cpus, $hmmer, $keep);
#-------------------------------------------------------------------------------------
#This is really the main body of the script


my $matches = {};

runhmmer($options) unless($options->{hmmerprecal});

#Parse the hmmer results
parsehmmer($options, $matches );

#Now that we have parsed the HMMER searc/scan results, graft the matched into the tree.
my $allResults = [];
graftMatches($options, $matches, $allResults);

#Now print out the results.
printResults($options, $allResults);

exit;
#--------------------------------------------------------------------------------------

sub processOptions {
  my ( $options, $help, $outfile, $directory, $fastafile, 
        $raxmlloc, $hmmerloc, $algo, $auto, $cpus, $keep) = @_;

  &usage if($help);

  #Make sure we know where we are going to write to.
  &usage ("Please specify output file\n") unless ($outfile);
  $options->{outfile} = $outfile;

  #Check the fasta file is defined, present and has size 
  &usage ("Please specify input fasta file\n") unless ($fastafile);
  if(!-e $fastafile){
    die "Your fasta file, $fastafile, does not exisit.\n";
  }elsif(! -s $fastafile){
    die "Your fasta file, $fastafile, has no size\n";
  }
  $options->{fastafile} = $fastafile;
  

  #Check that the PANTHER directory is present, and contains the HMM files
  &usage ("Please specify the directory of the pipeline\n") unless ($directory);
  # Need to make sure this directory is an absolute path, RaxML requires this.
  $directory = abs_path($directory);

  if(!-d $directory){
    die "Your directory, $directory, does not exisit.\n";
  }
  $options->{directory} = $directory;

  #We expect this directory to have a certain structure
  #TODO:Deal with different releases, hardcoded to 12.0
  my $pantherhmm = "$directory/PANTHER12_famhmm/PANTHER12.0_all_fam.hmm";
  if(!-e $pantherhmm){
    die "The PANTHER hmm file, $pantherhmm, does not exist.\n";
  }elsif( !-s $pantherhmm ){
    die "The PANTHER hmm file, $pantherhmm,  has no size.\n";
  }
  $options->{pantherhmm} = $pantherhmm;
  #Check that the directory containing the PANTHER alignments is present.
  my $pantherdir = "$directory/PANTHER12_Tree_MSF";

  if(!-d $pantherdir){
    die "The PANTHER alignments directory, $pantherdir, does not exisit.\n";
  }
  $options->{pantherdir} = $pantherdir;
  #-------------------------------------------------------------------------------------
  my $annotationfile = "$directory/PANTHER12_PAINT_Annotations/PANTHER12_PAINT_Annotatations_TOTAL.txt";

  if(!-e $annotationfile){
    die "The PANTHER annotation file, $annotationfile, does not exist.\n";
  }elsif( !-s $annotationfile ){
    die "The PANTHER hmm file, $annotationfile,  has no size.\n";
  }

  #Now read in the annotations file to work out which PANTHER entries have annotations.
  open ANO, "<", $annotationfile or die "cannot open $annotationfile\n";
  while(<ANO>){
    chomp;
    my @array = split(/\t/);
    $options->{annotations}->{$array[0]} = $array[1];
    my @a = split(/:/,$array[0]);
    $options->{pthrAnn}->{$a[0]} = 1; # store the PANTHER families with annotations
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
    $algo = autodetect($options);
  }

  #Need to check for the presence of the hmmpress files
  if($algo eq "hmmscan"){
    my $error = 0;
    foreach my $ext (qw(h3f h3i h3m h3p)){
      if(!-e $pantherhmm.".$ext"){
        $error = 1; 
      }
    }
    if($error){
      my $hmmpress = '';
      if($hmmerloc){
        $hmmpress .= $hmmerloc."/";
      }
      $hmmpress .= "hmmpress";

      system("$hmmpress $pantherhmm") and die "Error running hmmpress\n";
    }
  }
  $options->{algo} = $algo; 
  if($hmmer){
    #Check that the algo is set.
    #Then check that the file is present.
    $options->{hmmerprecal} = 1; 
    $options->{hmmerout} = $hmmer; 
  }else{
    $options->{hmmerout} = "$fastafile.$algo.out"; 
  }

  $options->{hmmerloc} = $hmmerloc if($hmmerloc);
  $options->{raxmlloc} = $raxmlloc if($raxmlloc);
  $options->{keep} = 1 if($keep); 
  #In the case of hmmer cpu=0 switches off threading. 
  #Otherwise hmmer will autodetect the number of CPUs
  # and use all of them, which does not work well on HPCs.

  if(!defined($cpus)){
    $cpus = 0;
  }
  $options->{cpus} = $cpus;
}

#-------------------------------------------------------------------------------------

sub runhmmer {
  my($options) = @_;
  if (!-s $options->{hmmerout}){
    my $hmmercommand = '';
    if (defined($options->{hmmerloc})) {
      $hmmercommand = $options->{hmmerloc}."/";
    }

    $hmmercommand .= $options->{algo}. " --notextw --cpu ". 
                     $options->{cpus}. " -o ".
                     $options->{hmmerout}. " ".
                     $options->{pantherhmm}." ".
                     $options->{fastafile}." > /dev/null";
    
    system($hmmercommand) and die "Error running $hmmercommand";
  }
}


#-------------------------------------------------------------------------------------
#Now parse the hmmer output

sub parsehmmer {
  my ($options, $matches ) = @_;

  my $queryid;my $matchpthr; 
  my $hmmline =0; my $domainline =0; my $alignmentline =0;
  my $matchn =0;

  #TODO - this is written for hmmscan, will need to change for hmmscan.
  open HMM, "<",  $options->{hmmerout} or die "cannot open ".$options->{hmmerout}."\n";
  while(<HMM>){
    #Ignore hmmer header lines.
    next if(/^#/);

    $hmmline++;
    if ($_ =~ /^Query:\s+(\S+) /){
      $queryid = $1;
      next;
    }elsif ($_ =~ /^>> (PTHR[0-9]+)/){ # use the first match
      $matchn++;
      if ($matchn==1){
        $matchpthr = $1;
      }
      next;
    }

    #This simple takes the first alignment (best match)
    #What happens in the case of gene fusion events? 
      # How many times does a sequence match two PANTHER entries?
    if ($matchn ==1){
      if ($hmmline >= 4){
        if ($_ =~ /!/){
	        my @a = split(/ +/);
	        push(@{$matches->{$matchpthr}->{$queryid}->{hmmstart}}, $a[7]);
	        push(@{$matches->{$matchpthr}->{$queryid}->{hmmend}},   $a[8]);
        }
      }

      if ($_ =~ /== domain/){
        $alignmentline = $hmmline;
      }

      if ($hmmline == $alignmentline+1){
        my @a = split(/ +/,$_);
	      push(@{$matches->{$matchpthr}->{$queryid}->{hmmalign}},   $a[3]);
      }elsif ($hmmline == $alignmentline+3){
        my @a = split(/ +/,$_);
	      push(@{$matches->{$matchpthr}->{$queryid}->{matchalign}},   $a[3]);
      } 
    }

    if ($_ =~ /^\/\/$/){
      ($queryid,$matchpthr,$domainline,$alignmentline) = ();
      $hmmline = $domainline = $alignmentline = $matchn = 0;
    }
  }
  close HMM;
}



#-------------------------------------------------------------------------------------

sub graftMatches {
  my ( $options, $matches, $allResults) = @_;
  foreach my $pthr (keys(%$matches)){
  #Make sure that we have annotations for this!
    unless (exists $options->{pthrAnn}->{$pthr}){
      foreach my $queryid (keys %{$matches->{$pthr}}){
        warn "$queryid\t$pthr\tAnnotationsNotAvailable\n";
      }
      next;
    }

    #Now get the length of the file.
    my $pthrAlignLength = _getAlignLength($options, $pthr);
    if($pthrAlignLength < 1){
      warn "Could not find alignment for $pthr\n";
      next;
    }
  
    foreach my $queryid (keys %{$matches->{$pthr}}){
      my $res = _graftPipeline($options, $queryid, $pthr, $pthrAlignLength, $matches->{$pthr}->{$queryid});
      if(defined $res){
        push(@$allResults, $res);
      }
    }
  }
}
#-------------------------------------------------------------------------------------

sub printResults{
  my($options, $allResults) = @_; 
  open (FINALOUT, '>',  $options->{outfile}) or die "cannot open ".$options->{outfile}."\n";
  foreach my $res (@$allResults){
    print FINALOUT $res;
  }
  close FINALOUT;
}

#---------------------------------------------------------------
sub _getAlignLength {
  my($options, $matchpthr) = @_;

  #Deterine the length of the PANTHER alignment. This does assume
  #that every position in the alignment is a match state.
  #This also assumes that the alignment are wrapped at 80 chars

  my $hmmalignmsf = $options->{pantherdir}."/$matchpthr.AN.fasta";
  if(!-e $hmmalignmsf or !-s $hmmalignmsf){
    return 0;
  }


  my $lines; my $lastline;
  open IN, "< $hmmalignmsf" or die "Could not open $hmmalignmsf:[$!]";
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
  
  return ($length_msf);
}  

#---------------------------------------------------------------
sub _graftPipeline{
  my ($options, $queryid, $matchpthr, $length_msf, $matchdata) = @_;
  
  my $querymsf = _querymsf($matchdata, $length_msf, $queryid);
  return unless($querymsf);
  
  my $queryfasta = _generateFasta($options, $querymsf, $matchpthr, $queryid);
  my $resString  = _runRAxMLAndAnnotate( $options, $queryid, $matchpthr, $queryfasta);
  
  unless($options->{keep}){
    my $command = "rm -rf ".$options->{directory}."/tmp/*";
    #TODO - try and replace with perl solution.
    system($command);
  }
  return($resString);
}

#---------------------------------------------------------------

sub _querymsf {
  my ($matchdata, $length_msf, $queryid) = @_;

  #Now build up the query sequence    
  my $querymsf; #aligned sequence placeholder

  my $start = $matchdata->{hmmstart}->[0];
  #N-terminaly padd the sequence
  foreach my $i (1..$start-1){
    $querymsf .= "-";
  }
  
  #For the first element/domain, extract the query string
  my @aligna = split(//,$matchdata->{matchalign}->[0]);
  my @haligna = split(//,$matchdata->{hmmalign}->[0]);
  my $alignas = length($matchdata->{hmmalign}->[0]);
  foreach my $a (0..$alignas-1){
    next if (($haligna[$a] eq ".")); #hmm insert state
    $querymsf .= $aligna[$a];
  }

  #Now we need to add in any additional domains in the next part
  my $domains = scalar @{$matchdata->{matchalign}};
  foreach my $j (1..$domains-1){ 
    my $start = $matchdata->{hmmstart}->[$j];
    my $end = $matchdata->{hmmend}->[$j-1];
    #This bridges the gap between the hits
    foreach my $i ($end+1..$start-1){
      $querymsf .= "-";
    }
    
    my @aligna = split(//,$matchdata->{matchalign}->[$j]);
    my @haligna = split(//,$matchdata->{hmmalign}->[$j]);
    my $alignas = length($matchdata->{hmmalign}->[$j]);
    foreach my $a (0..$alignas-1){
      next if ($haligna[$a] eq "."); #insert state, so skipping.
      $querymsf .= $aligna[$a];
    }
  }

  #Pad out as necessary
  my $prevend = $matchdata->{hmmend}->[$domains-1];
  foreach my $i ($prevend..$length_msf-1){
    $querymsf .= "-";
  }

  my $l = length($querymsf);
  unless ($l eq $length_msf){
    print "ERROR MSF of $queryid should have length $length_msf, actual length is $l\n";
    return 0;
  }

  $querymsf = uc($querymsf);
  return $querymsf;
} 


#---------------------------------------------------------------
# Write out the sequence file and concatenate.

sub _generateFasta {
  my ($options, $querymsf, $matchpthr, $queryid) = @_; 

  my @parts = $querymsf =~ /(.{1,80})/g;
  $queryid =~ s/[^\w]/\_/g;
  my $queryfasta = $options->{directory}."/tmp/$queryid.$matchpthr.fasta";
  open OUT,"> $queryfasta" or die "cannot open $queryfasta\n";
  print OUT ">query_$queryid\n";
  foreach my $line (@parts){
    print OUT "$line\n";
  }

  my $hmmalignmsf = $options->{pantherdir}."/$matchpthr.AN.fasta";
  #Concatenate the alignment. 
  open(A, "<", $hmmalignmsf) or die "Failed to open $hmmalignmsf for reading:[$!]\n";
  while(<A>){
    print OUT $_;
  }
  close(A);
  close OUT;
   
  return $queryfasta;
}
#---------------------------------------------------------------


sub _runRAxMLAndAnnotate {
  my ($options, $queryid, $matchpthr, $queryfasta) = @_;  
  ###############################################################
  ###run RAxML ################################################## 
  ###############################################################
  #TODO:make sure this tmp dir is unique to the process.
  $queryid =~ s/[^\w]/\_/g;
  my $raxmldir = $options->{directory}."/tmp/$matchpthr"."_$queryid"."_"."raxml$$";
  my $bifurnewick = $options->{pantherdir}."/$matchpthr.bifurcate.newick";
  if (! -e $bifurnewick) {
    print STDERR "no bifurcate newickfile for $matchpthr\n";
    return;
  };

  mkdir($raxmldir) or die "failed to make the direcorry $raxmldir:[$!]\n"; 
  
  my $raxmlcommand = '';
  if (defined($options->{raxmlloc})){
    $raxmlcommand = $options->{raxmlloc}."/";
  } 
  $raxmlcommand .= "raxmlHPC-SSE3 -f y -p 12345 -t $bifurnewick -G 0.05 -m PROTGAMMAWAG -T 4 -s $queryfasta -n $matchpthr -w $raxmldir >/dev/null";
    
  try{ 
    system($raxmlcommand); } 
  catch {
    warn "caught error for $queryid $matchpthr: $_";
  };

  ## now find the location of the graft sequence in the tree 
  my $mapANs;
  try { $mapANs = &mapto($raxmldir,$matchpthr,"query_$queryid"); } 
  catch{ warn "caught error in mapto for $queryid $matchpthr $_";};
  
  my $result;
  if(!$mapANs){  
    $mapANs = "root";
    my $annotation = $options->{annotations}->{$matchpthr.":".$mapANs};
    $result = "$queryid\t$matchpthr\t$annotation\n";
  }else{
    my $commonAN = &commonancestor($options, $mapANs, $matchpthr);
    if (!$commonAN) {$commonAN = "root"};
    my $annotation = $options->{annotations}->{$matchpthr.":".$commonAN};
    $result = "$queryid\t$matchpthr\t$annotation\n";
  }
  return $result;
}




#RDF checked up to here!




############find common ancestor of a list of ANs which are leaf genes#########

sub commonancestor{
  my ($options, $mapANs, $matchpthr) = @_;
  
  my @a = split(/;/,$mapANs);
  unless ($a[1]){
    return $mapANs;
  }
  my $newick = $options->{pantherdir}."/$matchpthr.newick";
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
  #my $io = IO::String->new($treestring);

  my $input = Bio::TreeIO->new(-fh => IO::String->new($treestring), -format => "newick");
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

sub autodetect {
  my ($options) = @_;  

  print "Reading HMM file\n";
  my $hmmsize = 0;
  open(HMM, "<", $options->{pantherhmm}) or die "Could not open $options->{pantherhmm}:[$!]\n";
  while(<HMM>){
    if(/^LENG  (\d+)/){
      $hmmsize += $1;
    }
  }
  close(HMM);

  #Now, multiple by 40, as residue for residues, a profile HMM is 40x bigger than a sequence
  $hmmsize = $hmmsize * 40;
  print "hmm database size in memory: $hmmsize\n";
  
  #Now, read the fasta file.  As soon as it is bigger than hmmsize variable, it is more efficient to use hmmsearch
  my $fastasize = 0;
  my $algo = 'hmmscan';

  open(FA, "<", $options->{fastafile}) or die "Could not open $options->{fastafile}:[$!]\n";
  while(<FA>){
    next if(/^>/);
    chomp;
    $fastasize += length;
    if($fastasize > $hmmsize){
      $algo= 'hmmsearch';
      last;
    }
  }
  close(FA);
  print "fasta file size in memory: $fastasize\n";
  print "Best algorithm is $algo\n";
  return $algo;
}

sub usage{
  my $errorm=shift;
  if ($errorm){
    print "\nError: $errorm\n";
  }

#TODO - this need updating.
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


