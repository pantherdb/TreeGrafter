# TreeGrafter: a new software tool for annotating uncharacterized protein sequences using annotated phylogenetic trees.
 
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

*  *  *  *  * 

# TreeGrafter - Version 1.01

TreeGrafter is a new software tool for annotating uncharacterized protein sequences,
using annotated phylogenetic trees.

The scripts are located in:
https://github.com/haimingt/TreeGrafter

The PANTHER data for this tool is located in:
ftp://ftp.pantherdb.org/downloads/TreeGrafter/

2/28/2018

Software authors: Haiming Tang, Robert D Finn

*  *  *  *  * 

## IMPORTANT!!
If you need GO annotations (and/or PANTHER subfamily classifications) for protein sequences from organisms that are already in the PANTHER trees, you do not need to run this tool.  You can download the annotations directly from:
ftp://ftp.pantherdb.org/sequence_classifications/current_release/

You can also parse them from the annotation file distributed with TreeGrafter:
PANTHER12_PAINT_Annotations/PANTHER12_leaf_GO_IBDannotations.tab (in ftp://ftp.pantherdb.org/downloads/TreeMethod/treeGrafter1.01_supplemental.tar.gz) 

Description of the file format: separated by tab 

First column: PANTHER gene label

Second column:  SF and GO anntations

Third column: PANTHER long gene id: 

Format: `uniprot_species_mnemonic|Various_Gene_Database=gene_id|UniprotKB=uniprotKB_id`

Example: `"MONBE|Gene=28576|UniProtKB=A9V8K6"`

For a complete list of existing species in the current version of PANTHER:
Check: http://pantherdb.org/panther/summaryStats.jsp

*  *  *  *  * 

## Prerequisites:

- RAxML 8.2.4  available at https://github.com/stamatak/standard-RAxML	
- HMMER 3.1b2  available at http://hmmer.org/download.html
- Perl
- Perl modules: 
  - `Bio::TreeIO` available at http://search.cpan.org/CPAN/authors/id/C/CJ/CJFIELDS/BioPerl-1.6.924.tar.gz
  - `Try::Tiny` available at http://search.cpan.org/~ether/Try-Tiny-0.28/lib/Try/Tiny.pm
  - `JSON::Parse` available at http://search.cpan.org/~bkb/JSON-Parse-0.49/lib/JSON/Parse.pod
  - `IO::String` available at http://search.cpan.org/~gaas/IO-String-1.08/String.pm

When using Perl scripts it can be challenging to install all of the dependencies and there
can be "drift" in functionality of the external Perl modules that this script was
originally written for. To migate this, we have generate a standalone executable
that contain all of the dependencies. 

##### OS X
If you are using OSX, you can skip installing the Perl modules, and use the executable 
directly.

##### Windows
If you are using Windows, you can use Cygwin to run the tool. You follow the additional steps below:

1. Install Cygwin and Perl
  - Install Cygwin from https://cygwin.com/install.html,
  - Install “Perl” and “make” in “Select Packages” while running “setup.exe”. 
2.	Open Cygwin and verify Perl installation with command: % perl -v. 
3.	Install Cpanminus http://search.cpan.org/~miyagawa/App-cpanminus-1.7040/lib/App/cpanminus.pm#Installing_to_system_perl. Cpanminus could be installed using command cpan App::cpanmius
4.	Install required modules: `Try::Tiny`, `BioPerl`, `JSON::Parse`, `IO::String`, using commands below:
    ```
    % cpanm Try::Tiny 
    % cpanm Bio::Perl 
    % cpanm JSON::Parse 
    % cpanm IO::String 
    ```
The location to  Perl and Perl modules must be defined in your `$PATH` variable.  If you have
 any questions on how to set up `$PATH`, please contact your UNIX system administrator.

The locations of RAxML and HMMER  are required as input for this tool if not defined in $PATH.

## Usage:

Download the PANTHER data 
```
% wget  ftp://ftp.pantherdb.org/downloads/TreeMethod/treeGrafter1.01_supplemental.tar.gz
```

Uncompress the PANTHER data
```
% tar xvfz treeGrafter1.01_supplemental.tar.gz
```

Download or clone the treeGrafter from github
Go to https://github.com/haimingt/TreeGrafting

```
% cd treeGrafting-master
% perl treeGrafter.pl -f <input fasta file> -o <output file> -d <directory path to treeGrafter1.01_supplemental> -r <optional: RAxML location> -algo <hmmsearch or hmmscan> -auto <Please specify either -algo <hmmscan|hmmsearch> or -auto, but not both> -k <keep temporary files: 1 for yes> -hmmer <for previously stored output of hmmscan or hmmsearch>
```

## Examples:
1. using the Test data, automatically choose hmmscan or hmmsearch:
    ```
    cd treeGrafting-master
    perl treeGrafter.pl -f ./Test/sample.fasta -o ./Test/sample.1.out -d ./Test/PANTHER_mini -auto
    ```

2. using the prevously stored hmmscan output

    **IMPORTANT!** provide -algo hmmscan for hmmscan.out; -algo hmmsearch for hmmsearch.out

    ```
    perl treeGrafter.pl -f ./Test/sample.fasta -o ./Test/sample.2.out -d ./Test/PANTHER_mini -algo hmmscan -hmmer ./Test/sample.fasta.hmmscan.out
    ```

3. using hmmsearch algo

    ```
    perl treeGrafter.pl -f ./Test/sample.fasta -o ./Test/sample.3.out -d ./Test/PANTHER_mini -algo hmmsearch
    ```

4. using the full dataset treeGrafter1.01_supplemental

    ```
    perl treeGrafter.pl -f ./Test/sample.fasta -o ./Test/sample.4.out -d ./path/to/treeGrafter1.01_supplemental -algo hmmsearch
    ```

5. using the executable on mac

    ```
    treeGrafter.osx -f ./Test/sample.fasta -o ./Test/sample.5.out -d ./path/to/treeGrafter1.01_supplemental -algo hmmsearch
    ```

The input file must be a list of sequences in fasta format

*  *  *  *  * 
Sample input fasta file:

`sample.fasta`

```
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
```

Sample output info:

`sample.ouput`

```
STRCO_Gene_CAB44512_UniProtKB_Q9XAS3	PTHR10000	SF:PTHR10000:SF8;NAME:Sugar phosphatase YidA; GO-Annotations:GO:0016311;GO:0016791;GO:0005737;
STRR6_EnsemblGenome_spr1125_UniProtKB_Q8DPK0		PTHR10000		    SF:PTHR10000:SF27;NAME:5-amino-6-(5-phospho-D-ribitylamino)uracil phosphatase YbjI;GO-Annotations:GO:0016311;GO:0016791;GO:0005737;
```

Output format, separated by tab:  
first column: gene id (non word characters replaced by '_');
second column:  best-matched PANTHER family 
third column: predicted SF and GO annotations for the gene id


## Docker container

Alternatively script can be run using Docker container

Latest Docker container is available in https://hub.docker.com/r/ningzhithm/treegrafter/

Docker installation instructions for Windows, OSX and Linux platforms are available on https://docs.docker.com/install/.

##### Run software in container

Pull Docker container and run test
```
$ docker run --rm --name treegrafter -v /path/to/output:/output ningzhithm/treegrafter:1.01 -f ./Test/sample.fasta -o /output/sample.1.out -d ./Test/PANTHER_mini -auto
Unable to find image 'ningzhithm/treegrafter:1.01' locally
1.01: Pulling from ningzhithm/treegrafter
cc1a78bfd46b: Already exists 
d2c05365ee2a: Pull complete 
231cb0e216d3: Pull complete 
3d2aa70286b8: Pull complete 
e80dfb6a4adf: Pull complete 
d0b05f6cad70: Pull complete 
2d2cd662bc46: Pull complete 
6cfaf1c5a3ac: Pull complete 
15fc59d44f95: Pull complete 
765395a78139: Pull complete 
27e5292b1937: Pull complete 
b63a1923c3c0: Pull complete 
268d0765615f: Pull complete 
28f9ca7819ff: Pull complete 
c4f4dd48c9e4: Pull complete 
46c6672b0ab5: Pull complete 
229f69bb16ad: Pull complete 
64a809c054ca: Pull complete 
6ec0ce3584d7: Pull complete 
5cf5741bea6d: Pull complete 
eb44c517d033: Pull complete 
442b151cd5d0: Pull complete 
e7c68b3b30cf: Pull complete 
4a7312bbc726: Pull complete 
2476efc68f5b: Pull complete 
fa7417f5079f: Pull complete 
d42bd620ebf8: Pull complete 
1b11a5c09de7: Pull complete 
Digest: sha256:85fb4be0512140c50838155f70a8ba3fdcc7556e0c1e9cd519cdc5ac0e382610
Status: Downloaded newer image for ningzhithm/treegrafter:1.01
Reading HMM file
hmm database size in memory: 10640
fasta file size in memory: 10816
Best algorithm is hmmsearch
```

PANTHER DB can be mounted as a data volume using
`-v /path/to/treeGrafter1.01_supplemental:/opt/supplemental`

```
$ docker run --rm --name treegrafter -v /path/to/sample:/sample -v /path/to/output:/output -v /path/to/treeGrafter1.01_supplemental:/opt/supplemental ningzhithm/treegrafter:1.01 -f /sample/sample.fasta -o /output/sample.1.out -d /opt/supplemental -auto
```


*  *  *  *  * 

Troubleshooting:

If you have any problems, please contact us at: feedback@pantherdb.org

*  *  *  *  * 
