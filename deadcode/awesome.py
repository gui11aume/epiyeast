#!/usr/bin/env python
# -*- coding:utf-8 -*-

import seeq
import string
import sys

from itertools import izip

comp = string.maketrans('gatcGATC', 'ctagCTAG')


def revcomp(seq):
   '''Reverse complement a DNA string.'''
   return seq[::-1].translate(comp)

   
def compare(seq1, seq2):
   '''Compute the Hamming match. Total number of matches'''
   len1 = len(seq1)
   len2 = len(seq2)
   overlap = 0
   for pos in range(0,min(len1,len2)):
      if seq1[pos] != 'N' and seq1[pos] == seq2[pos]:
         overlap = overlap + 1         
   return overlap


def merge(SEQ1, SEQ2, s):
   '''Find the consensus between two sequences and print it.
   The Ns of a sequence are replaced by the nucleotides of the
   other one. In case of conflict, the function returns witout
   printing.'''
   len1 = len(SEQ1)
   len2 = len(SEQ2)
   # Initialize consensus with first sequence.
   consensus = [a for a in SEQ1[:min(len1,len2)]]
   for pos in range(min(len1,len2)):
      # Return if nucleotide is unknown.
      if SEQ1[pos] == SEQ2[pos] == 'N':
         return
      # Return if conflict.
      if SEQ1[pos] != SEQ2[pos] and SEQ1[pos] != 'N' and SEQ2[pos] != 'N':
         return
      if SEQ1[pos] == 'N':
         consensus[pos] = SEQ2[pos]
   print '%s: %s' % (s, ''.join(consensus))


f = open ('/data/GF3_new1_11735_CCTTCA_read1.fastq','r')
g = open ('/data/GF3_new1_11735_CCTTCA_read2.fastq','r')

central = seeq.compile("CTCTCTTGCGAGATGATCCCGCATTTT", 7)


samples = {
   (20,22): 's1',#
   (20,23): 's2',
   (20,24): 's3',
   (20,25): 's4',#
   (21,22): 's5',
   (21,23): 's6',#
   (21,24): 's7',
   (21,25): 's8',
   (22,22): 's9',
   (22,23): 's10',
   (22,24): 's11',#
   (22,25): 's12',
   (23,22): 's13',
   (23,23): 's14',#
   (23,24): 's15',
   (23,25): 's16',#
}

linenumber = 0
for (read1,read2) in izip(f,g):
   linenumber = linenumber + 1
   if linenumber % 4 == 2:
      # First figure out the orientation of the read pair.
      if central.match(read1):
         # Read 1 is in forward orientation.
         SEQ1 = read1.rstrip()
         SEQ2 = revcomp(read2.rstrip())
      elif central.match(read2):
         # Read 2 is in forward orientation.
         SEQ1 = read2.rstrip()
         SEQ2 = revcomp(read1.rstrip())
      else:
         # The central portion was not found.
         continue

      # Demultiplex.
      lind = 'ACAGGCCGTACGCAGTTGTCGAA'
      rind = 'AGCAGAATTACCCTCCACGTTGATT'

# s1
# GGCCGTACGCAGTTGTCGAA
# AGCAGAATTACCCTCCACGTTG

      (L,R) = (0,0)
      for shift in range(4):
         SL = compare(SEQ1, lind[shift:])
         if SL > 15:
            L = len(lind)-shift
            break
      SR = compare(SEQ2[-len(rind):], rind)
      if SR > 15: R = len(rind)
      if R == 0:
         for shift in range(1,4):
            SR = compare(SEQ2[-len(rind)+shift:], rind[:-shift])
            if SR > 15:
               R = len(rind)-shift
               break

      S = compare(SEQ1, SEQ2)
      # Doing the stats, we are sure that we found the right frame if
      # we find more than 53 matches.
      THRESHOLD = 53
      if S > THRESHOLD:
         merge(SEQ1, SEQ2, samples.get((L,R), 'sX'))
         continue
      for shift in range(1,10):
         S = compare(SEQ1, SEQ2[shift:])
         if S > THRESHOLD:
            merge(SEQ1, SEQ2[shift:], samples.get((L,R), 'sX'))
            break
         S = compare(SEQ1[shift:], SEQ2)
         if S > THRESHOLD:
            merge(SEQ1[shift:], SEQ2, samples.get((L,R), 'sX'))
            break
