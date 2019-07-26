import string
import argparse
import os
import time
import datetime
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio.SeqIO import FastaIO

parser = argparse.ArgumentParser()
parser.add_argument('-t', "--tree_file")
parser.add_argument('-s', "--sto_file")
parser.add_argument('-o', "--fasta_file")
parser.add_argument('-n', "--num_tasks")
parser.add_argument('-d', '--base_dir')

def convert_family(tree_file, sto_file):
    seqs = []
    seq_map = {}
    node_map = {}
    bp_seqs = []

    with open(tree_file) as tree_f:
        for l in tree_f.readlines():
            if l.startswith("AN"):
                an_id = l.split(":")[0]
                long_id = l.split(":")[1]
                long_id = long_id.replace(";\n", "")
                node_map[long_id] = an_id

    with open(sto_file) as sto_f:
        for sl in sto_f.readlines():
            # l = l.replace(" ", "")
            # l = l.replace("\n", "")
            # l = l.replace("\r", "")
            if sl != "\n" and sl != "//\n" and not sl.startswith("#"):
                seqs.append(sl)

    trantable = str.maketrans("", "", string.ascii_lowercase)
    for s in seqs:
        seq_bits = s.split(" ")
        seq_bits = list(filter(bool, seq_bits))
        long_id = seq_bits[0]
        try:
            an_id = node_map[long_id]
        except:
            print("Node mapping missing for {} in file {}".format(long_id, tree_file))
            continue
            # an_id = node_map[long_id]
        sequence = seq_bits[1]
        sequence = sequence.replace("\n", "")
        sequence = sequence.replace(".", "")
        sequence = sequence.translate(trantable)
        if an_id in seq_map:
            seq_map[an_id] = seq_map[an_id] + sequence
        else:
            seq_map[an_id] = sequence

    for an_id in seq_map:
        bs = Seq(seq_map[an_id])
        sr = SeqRecord(seq=bs, id=an_id, description="")
        bp_seqs.append(sr)

    return bp_seqs

if __name__ == "__main__":
    args = parser.parse_args()

    node_num = int(os.environ['SLURM_PROCID'])
    num_tasks = int(args.num_tasks)
    counter = 0

    process_books = []
    base_dir = args.base_dir
    
    # books = os.listdir("/auto/pmd-02/pdt/pdthomas/panther/famlib/rel/PANTHER13.1/books/")
    # for b in books:
    for f in os.listdir(base_dir):
        # Use presence of .sto files as qualification for processing
        if not f.startswith("PTHR") or not f.endswith(".sto"):
            continue
        if counter == node_num:
            process_books.append(f.replace(".sto", ""))
        counter += 1
        if counter == num_tasks:
            counter = 0

    log_file = open("log/sto_to_fasta.log.{}".format(node_num), "w+")
    log_file.write("Processing {} books in task {}:\n".format(len(process_books), node_num))
    log_file.write(",".join(process_books) + "\n")

    for pb in process_books:
        log_file.write("Working on {}\n".format(pb))
        tree_file = "{}/Tree_MSF/{}.tree".format(base_dir, pb)
        sto_file = "{}/Tree_MSF/{}.sto".format(base_dir, pb)
        fasta_file = "{}/Tree_MSF/{}.AN.fasta".format(base_dir, pb)
        fasta_seqs = convert_family(tree_file, sto_file)
        # SeqIO.write(fasta_seqs, fasta_file, "fasta")
        with open(fasta_file, "w+") as wf:
            writer = FastaIO.FastaWriter(wf, wrap=80)  # treeGrafter.pl expects line length of 80
            writer.write_file(fasta_seqs)

    # fasta_seqs = convert_family("resources/PTHR31977.tree", "resources/PTHR31977.sto")
    # fasta_seqs = convert_family(args.tree_file, args.sto_file)
    # SeqIO.write(fasta_seqs, args.fasta_file, "fasta")

    log_file.write("Job task {} ended {}\n".format(node_num, datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S')))
    log_file.close()
