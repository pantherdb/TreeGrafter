# Copyright 2018 EMBL - European Bioinformatics Institute
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

FROM perl:latest
MAINTAINER Ola Tarkowska (EMBL-EBI) <olat@ebi.ac.uk>

RUN apt-get update
RUN apt-get -y install zip

# RAxML
RUN wget -O /opt/RAxML.zip https://github.com/stamatak/standard-RAxML/archive/master.zip 
RUN unzip /opt/RAxML.zip -d /opt && \
    rm -f /opt/RAxML.zip

RUN mkdir -p /opt/standard-RAxML-master
WORKDIR /opt/standard-RAxML-master

RUN make -f Makefile.gcc
RUN make -f Makefile.SSE3.gcc
RUN make -f Makefile.PTHREADS.gcc
RUN make -f Makefile.SSE3.PTHREADS.gcc

RUN cp raxmlHPC* /usr/local/bin/

ENV PATH=/opt/standard-RAxML-maste:$PATH


# HMMER
RUN mkdir -p /opt/hmmer-3.1b2
WORKDIR /opt/hmmer-3.1b2
RUN wget -O /opt/hmmer-3.1b2.tar.gz http://eddylab.org/software/hmmer3/3.1b2/hmmer-3.1b2.tar.gz
RUN tar -pxvzf /opt/hmmer-3.1b2.tar.gz -C /opt/hmmer-3.1b2 --strip-components=1 && \
  rm -f /opt/hmmer-3.1b2.tar.gz

RUN ./configure
RUN make
RUN make install


# treegrafter
RUN mkdir -p /opt/treegrafter
WORKDIR /opt/treegrafter

# perl modules
RUN cpanm Try::Tiny
RUN cpanm Bio::Perl
RUN cpanm JSON::Parse
RUN cpanm IO::String


RUN git clone https://github.com/haimingt/TreeGrafter.git /opt/TreeGrafter.git

WORKDIR /opt/TreeGrafter.git

ENTRYPOINT ["perl", "treeGrafter.pl"]

# Example CMD
# docker run --rm --name treegrafter -v /path/to/TreeGrafter.git/Test:/tmp treegrafter -f ./Test/sample.fasta -o /tmp/sample.1.out -d /tmp/PANTHER_mini -auto