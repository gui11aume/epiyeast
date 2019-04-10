#!/usr/bin/env python
# -*- coding:utf-8 -*-

import datetime
import os
import string
import sys

import seeq

from itertools import izip, combinations
from collections import defaultdict

from gzopen import gzopen


## EXCEPTIONS ##

class BadSpecifications(Exception):
   pass

class BadRead(Exception):
   pass

class NoSample(Exception):
   pass

class SpikeException(Exception):
   def __init__(self, spikename, *args, **kwargs):
      self.spikename = spikename
      super(Exception, self).__init__(*args, **kwargs)


## UTILITIES ##

def dist_less_than(seq1, seq2, cutoff):
   '''Convenience function to tell whether the Hamming distance
   between a primer and a read is less or equal to the given cutoff.
   Note that sequences are cut to the length of the shortest, so
   only the distance at the beginning of the sequences is counted.'''

   # Make sure 'seq1' is the smallest.
   if len(seq1) > len(seq2): (seq1,seq2) = (seq2,seq1)
   
   dist = 0
   for pos in range(len(seq1)):
      dist += seq1[pos] != seq2[pos]
      # This allows to quickly bail out 
      if dist > cutoff: return False

   return True


## CLASSES ##

class LaneInfo:
   '''Just a name space for the informations to record about
   how the lane was preprocessed.'''

   def __init__(self):
      # Initialize the total number of reads and the number of 
      # aberrant reads (those that cannot be recognized).
      self.ntotal    = 0
      self.naberrant = 0
       
      # Initialize the number of counts per sample and the 
      # numbers of counts for the reads that cannot be
      # demultiplexed.
      self.counts    = defaultdict(int)
      self.nosamples = defaultdict(int)

      # Initialize the number of counts for the spikes.
      self.spdict = defaultdict(lambda: defaultdict(int))

      # Sequences errors.
      self.nttot = 0
      self.nterr = 0

   def write_to_file(self, f):
      '''Write a log file.'''

      # Start with a time stamp.
      dt = datetime.datetime # (for short).
      f.write(dt.strftime(dt.now(),
         '%Y-%m-%d %H:%M:%S Preprocessing summary\n'))

      # Total and percent reads lost.
      lost = float(self.naberrant + sum(self.nosamples.values()))
      f.write('Reads lost:\t%d (%.2f %%)\n' % \
            (lost, 100 * lost / self.ntotal))

      # Aberrant reads.
      f.write('Aberrant reads:\t%d (%.2f %%)\n' % \
            (self.naberrant, 100 * self.naberrant / self.ntotal))

      # Reads per sample.
      f.write('Reads per sample:\n')
      for smpl,count in sorted(self.counts.items()):
         f.write('%s\t%d\n' % (smpl,count))

      # Demultiplexing failures.
      f.write('Demultiplexing failures:\n')
      for ((L,R),count) in sorted(self.nosamples.items()):
         f.write('%d/%d:\t%d\n' % (L,R,count))
      if self.spdict:
         f.write('Spikes:\n')
         for smpl,D in sorted(self.spdict.items()):
            f.write('%s:\n' % smpl)
            for sp,cnt in D.items():
               f.write('  %s:\t%d\n' % (sp,cnt))

      # Sequencing error rate.
      f.write('Sequencing errors:\n')
      f.write('Total:\t%d\n' % self.nttot)
      f.write('Mismatches:\t%d\n' % self.nterr)

      # End of report.
      f.write('---\n')

class MultiplexSpecifications:
   '''Stores information about multiplexing. Essentially a name space
   for attributes, together with a parser.'''

   # Attributes.
   cst     = ''
   samples = {}

   def __init__(self, cst, fwd, rev, samples, spikes=dict()):
      '''Simple constructor instantiating attributes.'''

      # Base multiplexing information.
      self.cst = seeq.compile(cst, len(cst)/5)
      self.samples = samples

      # Store the sequences of the primers used
      # in a frozen set for reference.
      self.Lseq = frozenset([fwd[-a:] for (a,b) in self.samples.keys()])
      self.Rseq = frozenset([rev[-b:] for (a,b) in self.samples.keys()])

      # Check that the primers are not too close to each
      # other, otherwise it will be impossible to demultiplex.
      # If all the primers are at a distance greater than
      # twice the cut-off, then no double hit is possible.
      for (a,b) in combinations(self.Lseq, 2):
         if dist_less_than(a,b, 2*len(fwd)/5):
            raise BadSpecifications('primer sequences too close', a,b)
      for (a,b) in combinations(self.Rseq, 2):
         if dist_less_than(a,b, 2*len(rev)/5):
            raise BadSpecifications('primer sequences too close', a,b)

      # Check that the spikes (if present) are not too close to
      # each other to avoid multiple matches. The distance between
      # any two spikes has to be at least 2. This does not guarantee
      # that the matches are unique but we have to accomodate the
      # experiments.
      for (a,b) in combinations(spikes.values(), 2):
         if dist_less_than(a,b, 1):
            raise BadSpecifications('spike sequences too close', a,b)

      # Check that the spikes can be found by looking for the
      # constant part.
      for spseq in spikes.values():
         if self.cst.match(spseq) is None:
            raise BadSpecifications('spike too divergent', spseq)

      # Spikes are not too close. Replace the values of the
      # dictionary by a compiled seeq pattern allowing 2 errors.
      self.spikes = dict()
      for spname in spikes:
         self.spikes[spname] = seeq.compile(spikes[spname], 1)

      # Check that the sample specification corresponds to
      # the primers used for the PCR.
      if len(fwd) != max([a for (a,b) in samples.keys()]):
         raise BadSpecifications('inconsistent sample keys')
      if len(rev) != max([b for (a,b) in samples.keys()]):
         raise BadSpecifications('inconsistent sample keys')



   @classmethod
   def parse(cls, fname):
      '''Factory creating an object from a specifications file. The
      file should be tab-separated and formatted as shown below.

      AGGTTTGGATCAGGATTTGCGCCTTTGGAT
      fwd	 TTGCTCTCGGTCAAGCTTTTAAA
      rev	 TCGAACAGGCCGTACGCAGTTGT
      1t0	 20	20
      1t1	 21	21
      1t14	 22	22
      2t0	 23	23
      2t1	 23	21
      2t14	 20	23
      spike1 ACTGGCGATGCGAGTGAGCTGAGTA
      spike2 ACTGGCGTTGCCAGTGAGCTGAGTA

      The first line is the constant part in forward direction.
      The next two lines are the primers used for the PCR (longest
      forms). The next lines specify the lengths of the primers
      used in each experiment.'''

      with open(fname) as f:
         # The first line is the constant fragment. Instantiate
         # a 'seeq' pattern allowing up to 20% errors.
         cst = next(f).rstrip()

         # Second and third lines are the primers.
         check,fwd = next(f).split()
         assert(check == "fwd")
         check,rev = next(f).split()
         assert(check == "rev")

         # Following lines are the lengths used in different samples
         # and possibly the spikes. Create a dictionary of samples
         # like the following:
         # samples = {
         #    (20,20): '1t0',
         #    (21,21): '1t1',
         #    (22,22): '1t14',
         #    (23,23): '2t0',
         #    (23,21): '2t1',
         #    (20,23): '2t14',
         # }
         # Possibly create a dictionary of spikes like the following:
         # spikes = {
         #    'spike1': 'ACTGGCGATGCGAGTGAGCTGAGTA',
         #    'spike2': 'ACTGGCGTTGCCAGTGAGCTGAGTA',
         # }
         samples = dict()
         spikes = dict()
         for line in f:
            if line.startswith('spike'):
               spikename,spikeseq = line.split()
               spikes[spikename] = spikeseq
            else:
               name, lfwd, lbwd = line.split()
               samples[(int(lfwd),int(lbwd))] = name

         return cls(cst, fwd, rev, samples, spikes)



class PairedEndRead:
   '''Sequence and read quality of paired-end reads.'''

   comp = string.maketrans('gatcGATC', 'ctagCTAG')

   def __init__(self, read1, read2, qual1, qual2):
      self.read1 = read1
      self.read2 = read2
      self.qual1 = qual1
      self.qual2 = qual2
   

   @classmethod
   def parse(cls, fname1, fname2):
      '''Iterator that yields objects from a pair of fastq files.'''

      with gzopen(fname1) as f, gzopen(fname2) as g:
         for lineno,(line1,line2) in enumerate(izip(f,g)):
            if lineno % 4 == 1:
               # Read the sequence.
               read1 = line1.rstrip()
               read2 = line2.rstrip()
            if lineno % 4 == 3:
               # Read the quality and yield the object.
               qual1 = line1.rstrip()
               qual2 = line2.rstrip()
               yield cls(read1, read2, qual1, qual2)
            

   def check_if_spike(self, specs):
      for spname,pattern in specs.spikes.items():
         if pattern.match(self.read1):
            raise SpikeException(spname)


   def orient(self, specs):
      '''Make sure that read1 is forward and read2 is reverse
      (according to specifications). Modifies the object in place.'''

      m = specs.cst.matchBest(self.read1)
      if m is not None:
         # Record the place of the hit for
         # the purposes of merging.
         self.anchor1 = m.matchlist[0]
         return
      m = specs.cst.matchBest(self.read2)
      if m is not None:
         # PE read must be reversed. Swap reads 1 and 2.
         (self.read1, self.read2) = (self.read2, self.read1)
         (self.qual1, self.qual2) = (self.qual2, self.qual1)
         self.anchor1 = m.matchlist[0]
         return

      raise BadRead('constant part not found')


   def identify_sample(self, specs):
      '''Look for similarity between the beginning of the reads
      and the primers used for multiplexing.'''

      (L,R) = (0,0)
      # Look for identity allowing 20% error.
      for seq in specs.Lseq:
         if dist_less_than(seq, self.read1, len(seq)/5):
            L = len(seq)
      for seq in specs.Rseq:
         if dist_less_than(seq, self.read2, len(seq)/5):
            R = len(seq)

      # Record the length of the primers.
      self.L = L
      self.R = R

      # Record the sample.
      try:
         self.sample = specs.samples[(L,R)]
      except KeyError:
         raise NoSample


   def merge(self, specs, info, trim=True):
      '''Merge the two reads (same as FLASH).'''

      # Never try to merge reads that were not oriented first.
      if not hasattr(self, 'anchor1'): self.orient(specs)

      # Reverse complement read 2.
      revread2 = self.read2[::-1].translate(self.comp)

      # Look for the constant part on the reversed read.
      m = specs.cst.matchBest(revread2)
      if m is None:
         raise BadRead('no conensus found')

      # The shift is the difference between the
      # positions of the anchors.
      self.anchor2 = m.matchlist[0]
      shift = self.anchor1[0] - self.anchor2[0]
      if self.anchor1[1] - self.anchor2[1] != shift:
         # The beginning and the end of the constant regions
         # must be at the same position (ruling out the presence
         # of an indel in a read and not in the other.
         raise BadRead('no conensus found')

      revqual2 = self.qual2[::-1]

      # Initialize consensus. We will need to update the string
      # in place so we need a 'bytearray', which is mutable.
      # The left part (if any) is initialized with read1, the
      # right part with read2. There are two situations that look
      # as shown below.
      #
      # Case 1: overlap is in 5' of the reads. Consensus must be
      # the union of the sequences.
      #
      #         ------------------->
      #                  <---------------------
      #         | read1  |       read2        |
      #
      # Case 2: overlap is in 3' of the reads. Consensus must be
      # the intersection of the sequences because the overhangs
      # are the sequencing adapters.
      #
      #                  ------------------->
      #         <---------------------
      #                  |   read2   |

      '''
      def there_is_a_low_quality_nucleotide(qu):
         for nt in qu:
            if nt-33 < 20: return True
         return False

      def there_is_a_low_quality_nucleotide_in_variable_regions(qu):
         #   [ oligo [ variable 1 [ constant [ variable 2 [ oligo ]
         #   ^       ^            ^          ^            ^
         #   0       L       anchor1[0]  anchor1[1]      -R
         if not hasattr(self, 'L') or not hasattr(self, 'R'):
            exit()
         for nt in qu[self.L:self.anchor1[0]]:
            if nt-33 < 20: return True
         for nt in qu[self.anchor1[1]:-self.R]:
            if nt-33 < 20: return True
         return False
      '''


      def pos(x): return  x if x > 0 else 0
      def neg(x): return -x if x < 0 else 0

      nttot = nterr = 0
      cs = bytearray(self.read1[:pos(shift)] + revread2[neg(shift):])
      qu = bytearray(self.qual1[:pos(shift)] + revqual2[neg(shift):])
      for i in range(pos(shift), min(len(self.read1), len(cs))):
         # Record the disagreement in order compute the error
         # rate of the batch.
         there_is_no_N = self.read1[i] != 'N' and cs[i] != 78
         nttot += 1 if there_is_no_N else 0
         if ord(self.read1[i]) != cs[i]:
            nterr += 1 if there_is_no_N else 0
            if self.qual1[i] > revqual2[i-shift]:
               cs[i] = self.read1[i]
               qu[i] = self.qual1[i]

      if 'N' in cs:
         raise BadRead('no consensus found')

      #if there_is_a_low_quality_nucleotide_in_variable_regions(qu):
      #   raise BadRead('low quality nucleotide')

      info.nttot += nttot
      info.nterr += nterr


      # If the primers need to be trimmed, do it now.
      if trim and hasattr(self, 'L') and hasattr(self, 'R'):
         self.cs = str(cs[self.L:-self.R])
      else:
         self.cs = str(cs)



def preprocess(specs, fname1, fname2):
   '''Pre-process paired fastq files.'''

   outfiles = dict()
   for name in specs.samples.values():
      outfiles[name] = open(name + '.txt', 'w')

   info = LaneInfo()

   for PEread in PairedEndRead.parse(fname1, fname2):
      info.ntotal += 1
      try:
         PEread.orient(specs)
         PEread.identify_sample(specs)
         PEread.check_if_spike(specs)
         PEread.merge(specs, info)
      except BadRead:
         info.naberrant += 1
         continue
      except NoSample:
         info.nosamples[(PEread.L, PEread.R)] += 1
         continue
      except SpikeException as sp:
         info.spdict[PEread.sample][sp.spikename] += 1
         info.counts[PEread.sample] += 1
         continue

      outfiles[PEread.sample].write(PEread.cs + '\n')
      info.counts[PEread.sample] += 1

   for outf in outfiles.values():
      outf.close()

   return info


if __name__ == '__main__':

   # argv[1]: specification file
   # argv[2]: fastq file read1
   # argv[3]: fastq file read2
   # argv[4]: output directory

   # Parse run specifications.
   specs = MultiplexSpecifications.parse(sys.argv[1])

   # Change directory for output.
   topdir = os.getcwd()
   if not os.path.exists(sys.argv[4]): os.makedirs(sys.argv[4])
   os.chdir(sys.argv[4])

   # Preprocess reads.
   info = preprocess(specs, sys.argv[2], sys.argv[3])

   # Go back to top directory for log file.
   os.chdir(topdir)

   # Write log file.
   with open('logs.txt', 'a') as logf:
      info.write_to_file(logf)
