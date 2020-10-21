import os
import argparse

### This script simply removes surrounding pair of parentheses.

parser = argparse.ArgumentParser()
parser.add_argument('tree_dir', help="Directory of .bifurcate.newick files needing fixing")

args = parser.parse_args()

fixed_count = 0

tree_dir = args.tree_dir
for f in os.listdir(tree_dir):
    if f.endswith(".bifurcate.newick"):
        tree_text = None
        new_tree_text = None
        with open(tree_dir + f) as rf:
            tree_text = rf.read()
            if tree_text.startswith("(") and tree_text.endswith(");"):
                new_tree_text = tree_text[1:len(tree_text)-2] + ";"
                fixed_count += 1
        # Overwrite tree file if new_tree_text
        if new_tree_text:
            with open(tree_dir + f, "w+") as wf:
                wf.write(new_tree_text)

print("{} .bifurcate.newick files fixed".format(fixed_count))