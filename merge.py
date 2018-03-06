#!/usr/bin/env python
# -*- coding:utf-8 -*-

import sys

from collections import defaultdict

from gzopen import gzopen

gcode = {
   'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
   'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
   'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
   'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
   'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
   'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
   'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
   'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
   'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
   'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
   'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
   'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
   'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
   'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
   'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
   'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
}

def translate(DNA):
   n = len(DNA)
   return ''.join([gcode.get(DNA[i:i+3], '_') for i in range(0,n,3)])


def main(f0, f1, f14):
   # Create one int dictionary per file.
   dict0 = defaultdict(int)
   dict1 = defaultdict(int)
   dict14 = defaultdict(int)

   for line in f0:
      DNA,count = line.split()
      dict0[DNA] = int(count)
   for line in f1:
      DNA,count = line.split()
      dict1[DNA] = int(count)
   for line in f14:
      DNA,count = line.split()
      dict14[DNA] = int(count)

   for DNA in set(dict0).union(dict1).union(dict14):
      print DNA, translate(DNA), dict0[DNA], dict1[DNA], dict14[DNA]


if __name__ == "__main__":
   with gzopen(sys.argv[1]) as f0, gzopen(sys.argv[2]) as f1, \
            gzopen(sys.argv[3]) as f14:
      main(f0, f1, f14)
