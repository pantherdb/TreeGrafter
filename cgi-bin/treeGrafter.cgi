#!/usr/bin/perl

#####
#####
##
## myprog.pl
##
#####
#####

# necessary libraries
# use FastaFile;
# use FamLibBuilder;
#use Algs::Hmmer;
use Data::UUID;
use CGI;
use Try::Tiny;

# default values
my $cgi = CGI->new;
#my $library = "PANTHER8.0";
#my $library = "./PANTHER8.0";
my $library = "../Current_Release4";
#my $library = "/opt/panther/pantherScoring/PANTHER8.0";
#my $library = "/var/www/panther_data/PANTHER8.0";
my $test_out = "./test.out";
#open (TO, ">$test_out");
# get command-line arguments
use Getopt::Std;
getopts('ao:i:e:vVh') || &usage();
&usage() if ($opt_h);           # -h for help
$outFile = $opt_o if ($opt_o);  # -o for (o)utput file (redirect STDOUT)
$inFile = $opt_i if ($opt_i);   # -i for (i)Input file (redirect STDIN)
$errFile = $opt_e if ($opt_e);  # -e for (e)rror file (redirect STDERR)
$verbose = 1 if ($opt_v);       # -v for (v)erbose (debug info to STDERR)
$verbose = 2 if ($opt_V);       # -V for (V)ery verbose (debug info STDERR)

# get input file if specified on command-line
&usage() if ($#ARGV > 0);
if ($#ARGV == 0) {
    &usage("Please specify only one input file.") if ($opt_i);
    $inFile = $ARGV[0];
}

# do file redirects
if ($outFile) {
    open (STDOUT, ">$outFile") || &usage("Cannot open $outFile");
}
if ($errFile) {
    open (STDERR, ">>$errFile") || &usage("Cannot open $errFile");
}

&pukeSetup() if ($verbose);

#######################################################
###  MAIN BEGINS HERE
#######################################################

#print TO "$0 starts at ". localtime() . "\n";
# my $flb = new FamLibBuilder($library,"prod");
# die "Cannot read library\n" unless ($flb->exists());

# Uniqify seq in and out files for each session w/ UUIDs
$ug = Data::UUID->new;
$uuid = $ug->create();
$uuid_str = $ug->to_string($uuid);

#### ****DONT ALTER FOLLOWING ***###############################
## Determine the form's REQUEST_METHOD (GET or POST) and split the form   #
## fields up into their name-value pairs.  If the REQUEST_METHOD was      #
## not GET or POST, send an error.                                        #

my %FORM;
if ($ENV{'REQUEST_METHOD'} eq 'GET') {
  @pairs = split(/&/, $ENV{'QUERY_STRING'});
} elsif ($ENV{'REQUEST_METHOD'} eq 'POST') {
  read(STDIN, $buffer, $ENV{'CONTENT_LENGTH'});
  # Split the name-value pairs
  @pairs = split(/&/, $buffer);
}
foreach $pair (@pairs) {
  local($name, $value) = split(/=/, $pair);
  $name =~ tr/+/ /;
  $name =~ s/%([a-fA-F0-9][a-fA-F0-9])/pack("C", hex($1))/eg;
  $value =~ tr/+/ /;
  $value =~ s/%([a-fA-F0-9][a-fA-F0-9])/pack("C", hex($1))/eg;
  $value =~ s/<!--(.|\n)*-->//g;
  $FORM{$name} = $value;
}

my $seq = $FORM{'sequence'};

unless ($seq) {
    print $cgi->header(-type=>'text/plain', -status=> '400 Bad Request');
    print "<Response>\n";
    print "  <Error>Missing input sequence</Error>\n";
    print "</Response>\n";
    die "Missing input sequence\n";
}

unless ($seq=~/^>/) {
  $seq = ">$uuid_str\n" . $seq;
}
# die "Missing input sequence\n" unless ($seq);

my $tmpIn = "../tmp/$uuid_str.fasta";
my $tmpOut = "../tmp/$uuid_str-out.txt";
my $tmpInHmmscan = "$tmpIn.hmmscan.out";

open(SEQIN, '>', $tmpIn) || die "Could not open $tmpIn\n";
print SEQIN $seq;
close(SEQIN);

try {
  # Add RAxML cmd to path for use in treeGrafter.pl
  local $ENV{PATH} = "$ENV{PATH}:/opt/panther/TreeGrafter/RAxML/standard-RAxML-master/";
  # now call treeGrafter script - use -k to keep tmp files (../resources/PANTHER14.1_data/tmp) but delete after use here
  my $cmd = "perl ../treeGrafter.pl -f $tmpIn -d ../resources/PANTHER14.1_data/ -algo hmmscan -o $tmpOut -k > /dev/null";
  print TO "treeGrafter.pl started\n";
  if (system($cmd)) { 
    die "FATAL ERROR in treeGrafter : $cmd \nSystem command returned error status: ($!)\n"; 
  } 

  open(TMP,"$tmpOut") || die "Cannot open output file\n";
  my $graftPoint = "";
  while(my $line = <TMP>) {
    my @array = split(/\t/,$line);
    chomp @array;
    $graftPoint = $array[3];
    $matchFam = $array[1];
    # print "$graftPoint\n";
    last;  # Only print first
  }

  # Get MSA file passed into RAxML - TODO: read and pass into XML <MSA>
  $tmpQueryFasta = "../resources/PANTHER14.1_data/tmp/$uuid_str.$matchFam.fasta";
  $tmpQueryFasta =~ s/-/\_/g;

  $query_seq = "";
  $seq_count = 0;
  open MSA, "<", $tmpQueryFasta or die "cannot open ".$tmpQueryFasta."\n";
  while(<MSA>) {
    if ($_ =~ />/) {
      $seq_count++;
    }
    if ($seq_count == 2) {
      last;
    }
    $query_seq = $query_seq.$_;
  }

  print $cgi->header();
  print "<Response>\n";
  print "  <GraphPointNode>$graftPoint</GraphPointNode>\n";
  print "  <MSA>$query_seq</MSA>\n";
  print "</Response>\n";

  close(TMP);

  unlink $tmpIn;
  unlink $tmpOut;
  unlink $tmpInHmmscan;
}
catch {
  print $cgi->header(-type=>'text/plain', -status=> '500 Internal Server Error');
  print "<Response>\n";
  print "  <Error>Exception: $_</Error>\n";
  print "</Response>\n";
  die "Exception: $_\n";
};
#print $cgi->header();

# open(FH,"$tmpScoreRes") || die "Cannot open output file\n";
# print "<scores>\n";
# my %scores;
# while(my $line = <FH>) {
#   my @array = split(/\t/,$line);
#   chomp @array;
#   my $hmmId = $array[1];
#   my $eval = $array[3];
#   print "  <id>$hmmId</id>\n";
#   print "  <eval>$eval</eval>\n";
#   print "  <align>\n";
#   #print TO "  <id>$hmmId</id>\n";
#   #print TO "  <eval>$eval</eval>\n";
#   #print TO "  <align>\n";
#   #now get alignment
#   my ($fam,$sf) = split(/:/,$hmmId); 
#   my $hmmType = ($sf ? $sf : "mag"); 
#   my $fle = $flb->getLibEntry($fam); 
#   my $hmmFile = $fle->hmmerHmmFile($hmmType); 
#   die "HMM file does not exist\n" unless (-s $hmmFile); 
#   my $tmpAlnRes = "/usr/tmp/scoreOnFly.$$.aln";
#   #my $hmmerObj = run Hmmer($tmpFasta,$hmmFile,undef,$tmpAlnRes);
#   my $run_hmmsearch = "./hmmsearch $hmmFile $tmpFasta > $tmpAlnRes";
#   system("$run_hmmsearch");
#   #print out alignment
#   open(FH,"$tmpAlnRes") || die "Cannot read result file\n"; 
#   my $start =0; 
#   while(my $line = <FH>) { 
#     last if ($line=~/^Internal pipeline statistics summary:/); 
#     if ($line=~/^Domain annotation for each sequence/) { 
#       $start = 1; 
#       next; 
#     } 
#     next unless ($start); 
#     $line=~s/>/&gt;/g;
#     $line=~s/</&lt;/g;
#     print $line;
#     print TO $line;	
#   }
#   close(FH);
#   print "  </align>\n";
#   print TO "  </align>\n";
# #  unlink $tmpAlnRes;
#   last;
# }
# close(FH);
# print "</Response>\n";
# #print TO "</scores>\n";
# unlink $tmpFasta;
# unlink $tmpScoreRes;
# #close (TO);



#######################################################
###  SUBROUTINES BEGIN HERE
#######################################################

sub pukeSetup {
    print STDERR "_"x50, "\n";
    print STDERR "Verbose level is high.\n";
    print STDERR "Input file is: $inFile\n";
    print STDERR "Output file is: $outFile\n";
    print STDERR "Error file is: $errFile\n";
    print STDERR "_"x50, "\n";
}

sub usage {
    my $error = shift;

    print "Error: $error\n\n" if ($error);

    print <<__EOT;

myprog.pl - a program to do something

Usage:
myprog.pl {args} 

Where args are:
\t-h for help (this message)
\t-i (i)nput file (redirect STDIN)
\t-o (o)utput file (redirect STDOUT)
\t-e (e)rror file (redirect STDERR)
\t-v (v)erbose (debug info to STDERR)
\t-V (V)ery verbose (debug info to STDERR)

__EOT

exit(-1);
}
