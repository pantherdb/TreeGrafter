Copyright (C) 2017 Paul Thomas
This file may be copied and redistributed freely, without advance permission,
provided that this Copyright statement is reproduced with each copy. 

LIMITATION OF WARRANTY
NOTHING IN THIS AGREEMENT WILL BE CONSTRUED AS A REPRESENTATION MADE OR
WARRANTY GIVEN BY PAUL THOMAS OR ANY THIRD PARTY THAT THE USE OF
DATA PROVIDED HEREUNDER WILL NOT INFRINGE ANY PATENT, COPYRIGHT, TRADEMARK
OR OTHER RIGHTS OF ANY THIRD PARTY. DATA IS PROVIDED "AS IS" WITHOUT
WARRANTY OF ANY KIND WHATSOEVER, EXPRESS OR IMPLIED, INCLUDING IMPLIED
WARRANTIES OF MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE. PAUL
THOMAS MAKES NO WARRANTY THAT ITS DATA DOES NOT CONTAIN ERRORS.

############################################################################

TreeGrafter - Version 1.01
TreeGrafter is a new software tool for annotating uncharacterized protein sequences,
 using annotated phylogenetic trees.

http://www.pantherdb.org/downloads/
10/8/2017

Written by Haiming Tang

#############################################################################

IMPORTANT!!
If you need GO annotations for existing protein sequences in PANTHER, like HUMAN/MOUSE/ECOLI/ etc.
You don't need to run this tool, it is better to directly parse out the results from PANTHER12_PAINT_Annotations/PANTHER12_leaf_GO_IBDannotations.tab  using gene ids of your choice

Description of the file format: seperated by tab 
First column: PANTHER gene label
Second column:  GO anntations
Third column: PANTHER long gene id: 
Format: uniprot_species_Mnemonic|Various_Gene_Database=gene_id|UniprotKB=uniprotKB_id
Example: "MONBE|Gene=28576|UniProtKB=A9V8K6"

For a complete list of existing species in PANTHER12
Check: http://pantherdb.org/panther/summaryStats.jsp

##############################################################################

##########
Requirements:

Perl
Perl modules: 
     Bio::TreeIO  available at http://search.cpan.org/CPAN/authors/id/C/CJ/CJFIELDS/BioPerl-1.6.924.tar.gz
     Try::Tiny  http://search.cpan.org/~ether/Try-Tiny-0.28/lib/Try/Tiny.pm
     JSON::Parse http://search.cpan.org/~bkb/JSON-Parse-0.49/lib/JSON/Parse.pod
     IO::String http://search.cpan.org/~gaas/IO-String-1.08/String.pm
RAxML 8.2.4  available at https://github.com/stamatak/standard-RAxML	
HMMER 3.1b2  available at http://hmmer.org/download.html

The location to  Perl and Perl modules must be defined in your $PATH variable.  If you have
 any questions on how to set up $PATH, please contact your UNIX system administrator.

The locations of RAxML and HMMER  are required as input for this tool if not defined in $PATH.

##########
Usage:

% cd treeGrafter1.01
% perl treeGrafter.pl -f <input fasta file> -o <output file> -d <full directory path to this script> -r <optional: RAxML location> -s <optional: hmmscan location> -k <keep temporary files: 1 for yes>

# An example:
% cd treeGrafter1.01
% pwd
/home/users/haiming/Documents/treeGrafter1.01
% perl treeGrafter.pl -f sample.fasta -o sample.out -d /home/users/haiming/Documents/treeGrafter1.01

For option k, enter 0 if you don't want to keep the tmp files.

the input file must be a list of sequences in fasta format

##########
sample input fasta file:

sample.fasta

>STRCO|Gene=CAB44512|UniProtKB=Q9XAS3
--DRTAYSLVATDLDGTLLRGDDTVSDRSLAALARVAGAGARHLVVTGRPAPRVRPLLDRLGCTGLAVCGQGAQVYDAGH
RMLWSVTLDRELAETALGIEAEVGQVHAAVDQDGVTP-----DYLMPHPTAVRVERRAQLWS-TPISKVLLR-HPELTDD
ELAATARAVVGSLATVTMSGPGTVELQPCGITKATGLALAAEHLGLERRRTIAFGDMPNDIPMFQWAAHGVAMAGAHPEL
KAVADEVTTTNEDDGVAVVLERIF--
>STAA8|EnsemblGenome=SAOUHSC_02831|UniProtKB=Q2FVA2
------VKAIAVDMDGTFLDSKKTYDKLRFEAITELRNRDITFIAASGNQYAKLKSIFGDRD--MYFISENGAVIYNG--
NELYNKSFNRQVFQQVVDLNMKQSIDQLVICGKH-TAF--KEDTRFYYHQLKEIDSLQQLPE-DDYVKIAFNIN-RETHP
NVDEEVATQFSNDIKLVSSGHDSIDIIMPNMTKGQALKRLLDKWEMSPSELMAFGDANNDKDMLAFAKHSYVMENSHDEE
LNIASAVAPSNDKQGVLTIIEQ----

Sample  output info:

sample.ouput

STRCO_Gene_CAB44512_UniProtKB_Q9XAS3	PTHR10000	GO:0016311;GO:0016791;GO:0005737;
STAA8_EnsemblGenome_SAOUHSC_02831_UniProtKB_Q2FVA2	PTHR10000	GO:0016311;GO:00167
91;GO:0005737;


Output format, seperated by tab:  
first column: gene id (non word characters replaced by '_');
second column:  best-matched PANTHER family 
third column: predicted GO anntoations for the gene id

##########

Troubleshooting:

Ultimately, if you have any problems, please contact us at: feedback@pantherdb.org

##########
