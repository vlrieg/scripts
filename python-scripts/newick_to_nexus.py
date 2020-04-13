#!/usr/bin/env python3
# run like newick_to_nexus.py < C1.trees > C1-trees.nex

import dendropy
import sys

ds = dendropy.DataSet()
tns = ds.new_taxon_namespace()
ds.attach_taxon_namespace(tns)

ds.read(file=sys.stdin, schema="newick")

ds.write(file=sys.stdout, schema="nexus")


