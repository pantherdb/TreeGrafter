#!/usr/bin/perl

print "In script\n";
use Bio::TreeIO;
print "Using folder\n";
use lib q(HAIMING);
print "Using ToNewick pm\n";
use HAIMING::ToNewick;
print "ToNewick used\n";
#use strict;
use warnings;

use Getopt::Std;
getopts('h:s:f:V');
my $source_dir = $opt_s if ($opt_s);  # -s for source directory of Panther and PAINT files - Maybe split this apart to separate monthly PAINT updates
my $destination_base_dir = $opt_f if ($opt_f);   # -f for destination base folder under which Tree_MSF, PAINT_Annotations, famhmm folders live

#my $limit =shift;
#my $annotationfile = "/auto/pmd-02/pdt/pdthomas/panther/famlib/dev/UPL/PANTHER12.0/DBLoad/annotation.dat";
#my $anqualiferfile = "/auto/pmd-02/pdt/pdthomas/panther/famlib/dev/UPL/PANTHER12.0/DBLoad/annotation_qualifier.dat";
#my $genenodefile = "/auto/pmd-02/pdt/pdthomas/panther/famlib/dev/UPL/PANTHER12.0/DBLoad/gene_node.dat";
#my $nodefile = "/auto/pmd-02/pdt/pdthomas/panther/famlib/dev/UPL/PANTHER12.0/DBLoad/node.dat";
# my $annotationfile = "/auto/pmd-02/pdt/pdthomas/panther/xiaosonh/UPL/PANTHER13.1/library_building/DBload/annotation.dat";
# my $anqualiferfile = "/auto/pmd-02/pdt/pdthomas/panther/xiaosonh/UPL/PANTHER13.1/library_building/DBload/annotation_qualifier.dat";
# my $genenodefile = "/auto/pmd-02/pdt/pdthomas/panther/xiaosonh/UPL/PANTHER13.1/library_building/DBload/gene_node.dat";
# my $nodefile = "/auto/pmd-02/pdt/pdthomas/panther/xiaosonh/UPL/PANTHER13.1/library_building/DBload/node.dat";
my $annotationfile = "$source_dir/DBload/annotation.dat";
my $anqualiferfile = "$source_dir/DBload/annotation_qualifier.dat";
my $genenodefile = "$source_dir/DBload/gene_node.dat";
my $nodefile = "$source_dir/DBload/node.dat";

my %gene_node; my %node;
my %annotations;
my %pthrs;

print "Parsing node file\n";
open IN, "< $nodefile" or die;
while(<IN>){
  chomp;
  my @a = split(/\t/);
  $node{$a[0]} = $a[1];
}
close IN;

print "Parsing gene-node file\n";
open IN, "< $genenodefile" or die;
while(<IN>){
  chomp;
  my @a = split(/\t/);
  $gene_node{$a[1]} = $a[0];
}
close IN;

print "Parsing annotation qualifier file\n";
open IN, "< $anqualiferfile" or die;
while(<IN>){
  chomp;
  my @a = split(/\t/);
  my $go = $a[1];
  $go =~ s/PTHR[0-9]+://;
  if ($a[2] eq "NOT"){
    $annotations{$a[0]}{$a[2]}{$go} =1;
  }
}
close IN;

print "Parsing annotation file\n";
open IN, "< $annotationfile" or die;
while(<IN>){
  chomp;
  my @a = split(/\t/);
  next unless ($a[2] eq "GO");
  my $go;
  if ($a[1] =~ /(PTHR[0-9]+):(GO:[0-9]+)/){
    $go = $2;
    $pthrs{$1} =1;
  }
  unless (exists $annotations{$a[0]}{"NOT"}{$go}){
    $annotations{$a[0]}{"GAIN"}{$go} =1;
  }
}
close IN;

# Run this from TreeGrafter/resources/PANTHER##.#_data/ folder
#my $predir = "/home/pmd-02/pdt/pdthomas/users/haiming/PANTHER12_treesandmsf";
#my $todir = "/home/pmd-02/pdt/pdthomas/users/haiming/TREEMETHOD/treegraftingpipeline/PANTHER12_Tree_MSF";
#my $predir = "./PANTHER13.1_treesandmsf";  # Apparently unused
my $treedir = "$source_dir/tree";
my $todir = "$destination_base_dir/Tree_MSF";

#my $leafan = "./PANTHER12_leaf_GO_annotations.tab"; 
#my $internalan = "./PANTHER12_internalnodes_GO_annotations.tab";
my $leafan = "$destination_base_dir/PAINT_Annotations/leaf_GO_annotations.tab";
my $internalan = "$destination_base_dir/PAINT_Annotations/internalnodes_GO_annotations.tab";
open LEAF, ">> $leafan" or die;
open INTER, ">> $internalan" or die;

print "Doing it's thing\n";
my $indicator =0;
foreach my $pthr (sort keys %pthrs){
  #if ($pthr eq $limit) {
  #  $indicator =1;
  #  next;
  #}
  #next unless ($indicator ==1);
  my $newick = $todir."/$pthr.newick";
  print "working on $pthr\n";
  unless (-s $newick){
    my $tree= $treedir."/$pthr.orig.tree";
    tonewickOnlyAN($tree,$newick);
  }

  &propagate($pthr,$newick);
#  last;

  my $bifurcate = $todir."/$pthr.bifurcate.newick";
  unless (-s $bifurcate){
    tobifurcate($newick,$bifurcate);
  }
}
close LEAF;
close INTER;

sub propagate{
  my $pthr = shift;
  my $newick = shift;

  my $treeio = Bio::TreeIO->new(-format =>"newick",-file=>$newick);
  my $tree = $treeio->next_tree;

  my $instance;
  my $root = $tree->get_root_node;
  unless ($root->id){
    my @tmp = $root->each_Descendent;
    $root = $tmp[0];
  }
  my $rootid = "$pthr:".$root->id;
  my $ptn = $node{$rootid};
  foreach my $go (keys %{$annotations{$rootid}{"GAIN"}}){
    $instance .= "GAIN:$go".";";
  }
  foreach my $go (keys %{$annotations{$rootid}{"NOT"}}){
    $instance .= "NOT:$go".";";
  }
  my $toprint = &processinstance($instance);
  print INTER "$pthr:root\t$toprint\t$ptn\n";
  print INTER "$rootid\t$toprint\t$ptn\n";
  &checkchild($root,$instance,$pthr);
}

close INTER;
close LEAF;

sub checkchild{
  my $node = shift;
  my $annot = shift;
  my $pthr = shift;
  return if ($node->is_Leaf);

  foreach my $child ($node->each_Descendent){
    my $childid = "$pthr:".$child->id;
    my $instance = $annot;
    foreach my $go (keys %{$annotations{$childid}{"GAIN"}}){
      $instance .= "GAIN:$go".";";
    }
    foreach my $go (keys %{$annotations{$childid}{"NOT"}}){
      $instance .= "NOT:$go".";";
    }
    my $toprint = &processinstance($instance);
    if ($child->is_Leaf){
      my $longid = $gene_node{$childid};
      print LEAF "$childid\t$toprint\t$longid\n";
    }
    else{
      my $ptn = $node{$childid};
      print INTER "$childid\t$toprint\t$ptn\n";
      &checkchild($child,$instance,$pthr);
    }
  }
}


sub processinstance{
  my $instance = shift;
  my @a = split(/;/,$instance);
  my %store;
  foreach my $a (@a){
    my @s = split(/:/,$a);
    if ($s[0] eq "GAIN"){
      $store{$s[1].$s[2]} ++;
    }
    else{
      $store{$s[1].$s[2]}--;
    }
  }
  my $out;
  foreach my $key (keys %store){
    if ($store{$key} >0){
      $out .= $key.";";
    }
  }
  return $out;
}
