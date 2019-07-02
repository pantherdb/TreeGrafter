use Getopt::Long;
use FamLibBuilder;

my $lib;

&GetOptions
(
    "l=s" => \$log,         # -l for (l)ib path (ie ./PANTHER_7.0)
    "s=s" => \$source_dir,  # -s for (s)ource path of Panther release files
    "d=s" => \$dest_dir,    # -d for (d)estination path of where to save output
);

my $processor = shift;

# my $dir = "/home/pmd-02/pdt/pdthomas/panther/famlib/rel/PANTHER13.1";
# my $tar = "/home/pmd-02/pdt/pdthomas/panther/debert/TreeGrafter/resources/PANTHER13.1_data/Tree_MSF";
my $dir = $source_dir;
my $tar = $dest_dir;

opendir (DH,$dir) or die "cannot open $dir";
my @files= readdir(DH);
closedir(DH);

$node_num = $ENV{SLURM_PROCID};

open LOG, ">$log$node_num" or die "$log: $!";
print LOG "$0 starts at ". localtime() . "\n";

my $flb = new FamLibBuilder($dir);
# &usage("Cannot find library $dir.\n") unless ($flb->exists());

# loop over all books in the library and divide into groups
my $all;  #arrayref of array
my $i = 0;  
for my $f (sort $flb->bookNames) {  #put one book into one group at a time
    push @{$all->[$i]}, $f;
    $i++;
    $i = 0 if $i == $processor;
}

$node_num--; # array index = which group minus one
my @books = @{$all->[$node_num]};
print LOG "group_members:". join(",", @books) ."\n";

# foreach my $file (@files){
foreach my $file (@books){
  if ($file =~ /PTHR([0-9]+)/){
    # my $number = $1;
    # next unless ($number%99 == $i);


    my $tree = $dir."/books/$file/tree.tree";
    my $ttree = $tar."/$file.tree";
    `cp $tree $ttree`;
    
    my $hmm = $dir."/books/$file/hmmer.hmm";
    my $fasta = $dir."/books/$file/cluster.fasta";
    my $sto = $tar."/$file.sto";
    
    `/home/pmd-02/pdt/pdthomas/panther/scripts/hmmer-3.1b2/binaries/hmmalign -o $sto $hmm $fasta` unless -e $sto;

    print LOG "$file finished\n";
  }
}

close LOG;
