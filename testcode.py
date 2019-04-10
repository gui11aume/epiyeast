import sys
import tempfile
import unittest

from StringIO import StringIO
from textwrap import dedent

import preprocess

class TestUtils(unittest.TestCase):

   def test_dist_less_than(self):
      # Use 'dlt' for short.
      dlt = preprocess.dist_less_than

      # Test case 1 (3 differences).
      a = 'CGTAGTCGAGGAGCGCTGAGCT'
      b = 'CGTAGTCgAGGaGCGcTGAGCT'
      self.assertTrue(dlt(a, b, 3))
      self.assertFalse(dlt(a, b, 2))
      # Check symmetry
      self.assertTrue(dlt(b, a, 3))
      self.assertFalse(dlt(b, a, 2))

      # Test case 2 (3 differences).
      a = 'CGTAGTCGAGGAGCGCTGAGCT'
      b = 'CGTAGTCgAGGaGCGcTGAGCTcgatggctagag'
      self.assertTrue(dlt(a, b, 3))
      self.assertFalse(dlt(a, b, 2))
      # Check symmetry
      self.assertTrue(dlt(b, a, 3))
      self.assertFalse(dlt(b, a, 2))

class TestLaneInfo(unittest.TestCase):

   def test_constructor(self):
      info = preprocess.LaneInfo()
      self.assertEqual(info.ntotal, 0)
      self.assertEqual(info.naberrant, 0)
      self.assertEqual(info.counts, dict())
      self.assertEqual(info.nosamples, dict())
      self.assertEqual(info.spdict, dict())

   def test_write_to_file(self):
      info = preprocess.LaneInfo()
      info.ntotal = 1000000
      info.naberrant = 100000
      info.counts = {
         '1t1': 400000,
         '1t2': 400000,
      }
      info.nosamples = {
         (0,0): 100000,
      }
      info.spdict = {
         '1t1': { 'sp1': 100 },
         '1t2': { 'sp1': 100 },
      }
      info.nttot = 10
      info.nterr = 1

      buffer = StringIO()
      info.write_to_file(buffer)

      txt = '''Reads lost:\t200000 (20.00 %)
         Aberrant reads:\t100000 (10.00 %)
         Reads per sample:
         1t1\t400000
         1t2\t400000
         Demultiplexing failures:
         0/0:\t100000
         Spikes:
         1t1:
           sp1:\t100
         1t2:
           sp1:\t100
         Sequencing errors:
         Total:\t10
         Mismatches:\t1
         ---'''

      check = '\n'.join(buffer.getvalue().splitlines()[1:])
      self.assertEqual(check, dedent(' '*9 + txt))
      
      buffer.close()


class TestMultiplexSpecifications(unittest.TestCase):

   def test_constructor(self):
      # Use 'MS' for short.
      MS = preprocess.MultiplexSpecifications

      # Example taken from segment 1.
      cst = 'AGGTTTGGATCAGGATTTGCGCCTTTGGAT'
      fwd = 'TTGCTCTCGGTCAAGCTTTTAAA'
      rev = 'ACAACTGCGTACGGCCTGTTCGA'
      samples = {
         (20,20): '1t0',
         (21,21): '1t1',
         (22,22): '1t14',
         (23,23): '2t0',
         (23,21): '2t1',
         (20,23): '2t14',
      }

      ms = MS(cst, fwd, rev, samples)

      # Check attributes.
      Lseq = frozenset([
         'TTGCTCTCGGTCAAGCTTTTAAA',
          'TGCTCTCGGTCAAGCTTTTAAA',
           'GCTCTCGGTCAAGCTTTTAAA',
            'CTCTCGGTCAAGCTTTTAAA'
      ])
      Rseq = frozenset([
         'ACAACTGCGTACGGCCTGTTCGA',
          'CAACTGCGTACGGCCTGTTCGA',
           'AACTGCGTACGGCCTGTTCGA',
            'ACTGCGTACGGCCTGTTCGA',
      ])


      self.assertIsNotNone(ms.cst.match(cst))
      self.assertEqual(ms.Lseq, Lseq)
      self.assertEqual(ms.Rseq, Rseq)
      self.assertEqual(ms.samples, samples)
      self.assertEqual(ms.spikes, dict())

      # Now add some spikes
      spikes = {
         'sp1': 'AGGTTTGGATCAGGATTTGCGCCTTTGGAT',
         'sp2': 'CAGTTTGGATCAGGATTTGCGGCTTTGGAA',
      }

      ms = MS(cst, fwd, rev, samples, spikes)

      self.assertIsNotNone(ms.cst.match(cst))
      self.assertEqual(ms.Lseq, Lseq)
      self.assertEqual(ms.Rseq, Rseq)
      self.assertEqual(ms.samples, samples)

      self.assertIsNotNone(ms.spikes['sp1'].match(spikes['sp1']))
      self.assertIsNotNone(ms.spikes['sp2'].match(spikes['sp2']))


   def test_triangle_inequality(self):
      # Use 'MS' for short.
      MS = preprocess.MultiplexSpecifications

      cst = 'AGGTTTGGATCAGGATTTGCGCCTTTGGAT'
      # When deleting the first nucleotides of this primer,
      # the sequences will have distance 0 so the condition
      # on the triangle inequality should fail.
      fwd = 'AAAAAAAAAAAAAAAAAAAAAAA'
      rev = 'ACAACTGCGTACGGCCTGTTCGA'
      samples = {
         (20,20): '1t0',
         (21,21): '1t1',
         (22,22): '1t14',
         (23,23): '2t0',
         (23,21): '2t1',
         (20,23): '2t14',
      }

      # Check that error is raised.
      with self.assertRaises(preprocess.BadSpecifications):
         MS(cst, fwd, rev, samples)
      with self.assertRaises(preprocess.BadSpecifications):
         MS(cst, rev, fwd, samples)

      # This time give a better forward primer, but bad spikes.
      fwd = 'TTGCTCTCGGTCAAGCTTTTAAA'
      spikes = {
         'sp1': 'CTGACTGGATCGGCGCGTAG',
         'sp2': 'CTGACTGGATCGGCGCGTAT',
      }
      # Check that error is raised.
      with self.assertRaises(preprocess.BadSpecifications):
         MS(cst, fwd, rev, samples, spikes)


   def test_primer_consistency(self):
      # Use 'MS' for short.
      MS = preprocess.MultiplexSpecifications

      # Example taken from segment 1.
      cst = 'AGGTTTGGATCAGGATTTGCGCCTTTGGAT'
      fwd = 'TTGCTCTCGGTCAAGCTTTTAAA' # Length 23.
      rev = 'ACAACTGCGTACGGCCTGTTCGA' # Length 23.

      # Check that error is raised.
      samples = { (24,23): '2t0', }
      self.assertRaises(preprocess.BadSpecifications, MS,
            cst, fwd, rev, samples)
      samples = { (23,24): '2t0', }
      self.assertRaises(preprocess.BadSpecifications, MS,
            cst, fwd, rev, samples)

   def test_parse(self):
      # Test case taken from segment 1.
      text = '''AGGTTTGGATCAGGATTTGCGCCTTTGGAT
         fwd\tTTGCTCTCGGTCAAGCTTTTAAA
         rev\tACAACTGCGTACGGCCTGTTCGA
         1t0\t20\t20
         1t1\t21\t21
         1t14\t22\t22
         2t0\t23\t23
         2t1\t23\t21
         2t14\t20\t23
         spike1\tAGGTTTGGATCAGGATTTGCGCCTTTGGAT
         spike2\tTGGTTTGGATCAGGATTTGCGCCTTTGGAG'''

      tmpf = tempfile.NamedTemporaryFile()
      tmpf.write(dedent(' '*9 + text))
      tmpf.seek(0)
      
      ms = preprocess.MultiplexSpecifications.parse(tmpf.name)
      tmpf.close()

      # Test attributes.
      Lseq = frozenset([
         'TTGCTCTCGGTCAAGCTTTTAAA',
          'TGCTCTCGGTCAAGCTTTTAAA',
           'GCTCTCGGTCAAGCTTTTAAA',
            'CTCTCGGTCAAGCTTTTAAA'
      ])
      Rseq = frozenset([
         'ACAACTGCGTACGGCCTGTTCGA',
          'CAACTGCGTACGGCCTGTTCGA',
           'AACTGCGTACGGCCTGTTCGA',
            'ACTGCGTACGGCCTGTTCGA',
      ])

      samples = {
         (20,20): '1t0',
         (21,21): '1t1',
         (22,22): '1t14',
         (23,23): '2t0',
         (23,21): '2t1',
         (20,23): '2t14',
      }

      cst = 'AGGTTTGGATCAGGATTTGCGCCTTTGGAT'
      sp1 = 'AGGTTTGGATCAGGATTTGCGCCTTTGGAT'
      sp2 = 'TGGTTTGGATCAGGATTTGCGCCTTTGGAG'
      self.assertIsNotNone(ms.cst.match(cst))
      self.assertEqual(ms.Lseq, Lseq)
      self.assertEqual(ms.Rseq, Rseq)
      self.assertEqual(ms.samples, samples)
      self.assertIsNotNone(ms.spikes['spike1'].match(sp1))
      self.assertIsNotNone(ms.spikes['spike2'].match(sp2))


class TestPairedEndRead(unittest.TestCase):
   def test_constructor(self):
      # Use 'PER' for short.
      PER = preprocess.PairedEndRead

      read1 = 'ACTGCTGAGC'
      read2 = 'CTGGAGAGCTTAGA'
      qual1 = 'BB3B@GEGFG'
      qual2 = 'BBCCCGGGGGGFGD'

      PEread = PER(read1, read2, qual1, qual2)

      self.assertEqual(PEread.read1, read1)
      self.assertEqual(PEread.read2, read2)
      self.assertEqual(PEread.qual1, qual1)
      self.assertEqual(PEread.qual2, qual2)

   def test_parse(self):
      # Use 'PER' for short.
      PER = preprocess.PairedEndRead

      text1 = '''@HWI-D00733:44:C7
         ACAACTGCGTACGGCCTGTT
         +
         BBBBCGG1FGGG/EGGG0FG
         @HWI-D00733:44:C7
         CAACGGCGTACGGCCTGTTC
         +
         A:BB0BCCGGGG@GGGE?F1'''

      text2 = '''@HWI-D00733:44:C7
         TTCTCTCGGTCAAGCTTTTA
         +
         =3=:A>FGGGGGGGGGGGGG
         @HWI-D00733:44:C7
         TGCTCTCGGTAAAGCTTTTA
         +
         BBBBB11=;F0>DEGDEG1D'''

      tmpf1 = tempfile.NamedTemporaryFile()
      tmpf2 = tempfile.NamedTemporaryFile()

      tmpf1.write(dedent(' '*9 + text1))
      tmpf2.write(dedent(' '*9 + text2))

      tmpf1.seek(0)
      tmpf2.seek(0)

      (PE1, PE2) = [a for a in PER.parse(tmpf1.name, tmpf2.name)]

      # Read 1.
      self.assertEqual(PE1.read1, 'ACAACTGCGTACGGCCTGTT')
      self.assertEqual(PE1.read2, 'TTCTCTCGGTCAAGCTTTTA')
      self.assertEqual(PE1.qual1, 'BBBBCGG1FGGG/EGGG0FG')
      self.assertEqual(PE1.qual2, '=3=:A>FGGGGGGGGGGGGG')

      # Read 2.
      self.assertEqual(PE2.read1, 'CAACGGCGTACGGCCTGTTC')
      self.assertEqual(PE2.read2, 'TGCTCTCGGTAAAGCTTTTA')
      self.assertEqual(PE2.qual1, 'A:BB0BCCGGGG@GGGE?F1')
      self.assertEqual(PE2.qual2, 'BBBBB11=;F0>DEGDEG1D')

   def test_orient(self):
      # Setup of the test (example taken from segment 1).
      cst = 'AGGTTTGGATCAGGATTTGCGCCTTTGGAT'
      fwd = 'TTGCTCTCGGTCAAGCTTTTAAA'
      rev = 'ACAACTGCGTACGGCCTGTTCGA'
      samples = {
         (20,20): '1t0',
         (21,21): '1t1',
         (22,22): '1t14',
         (23,23): '2t0',
         (23,21): '2t1',
         (20,23): '2t14',
      }

      ms = preprocess.MultiplexSpecifications(cst, fwd, rev, samples)

      # Use 'PER' for short.
      PER = preprocess.PairedEndRead

      read1 = 'ACTGCGTACGGCCTGTTCGACAAGTCGACTCCTGACCTGGAAAGCGCTT' \
         'CATCCACTGGCGCAAATCCTGATCCAAACATTTCTATCCCTTTGATGGCTCCAA' \
         'TGGCATCTTTAAAAGCTTANCC'
      read2 = 'GCTCTCGGTCAAGCTTTTAAAAAGGCCATAGGGAATGTCAAAGGTATAG' \
         'CAAGGTTTGGATCAGGATTGGCGACTTTGGGTGAAGCACTTTCTAGATCAGTTG' \
         'TTGACTTGTCGAACAGGCCGTA'
      qual1 = ''
      qual2 = ''

      PEread = PER(read1, read2, qual1, qual2)
      PEread.orient(ms)

      self.assertEqual(PEread.read1, read2)
      self.assertEqual(PEread.read2, read1)

   def test_identify_sample(self):
      # Test setup (taken from segment 1).
      cst = 'AGGTTTGGATCAGGATTTGCGCCTTTGGAT'
      fwd = 'TTGCTCTCGGTCAAGCTTTTAAA'
      rev = 'ACAACTGCGTACGGCCTGTTCGA'
      samples = {
         (20,20): '1t0',
         (21,21): '1t1',
         (22,22): '1t14',
         (23,23): '2t0',
         (23,21): '2t1',
         (20,23): '2t14',
      }

      ms = preprocess.MultiplexSpecifications(cst, fwd, rev, samples)

      # Use 'PER' for short.
      PER = preprocess.PairedEndRead

      # The are already in the proper orientation.
      read1 = 'GCTCTCGGTCAAGCTTTTAAAAAGGCCATAGGGAATGTCAAAGGTATAG' \
         'CAAGGTTTGGATCAGGATTGGCGACTTTGGGTGAAGCACTTTCTAGATCAGTTG' \
         'TTGACTTGTCGAACAGGCCGTA'
      read2 = 'AACTGCGTACGGCCTGTTCGACAAGTCGACTCCTGACCTGGAAAGCGCTT' \
         'CATCCACTGGCGCAAATCCTGATCCAAACATTTCTATCCCTTTGATGGCTCCAA' \
         'TGGCATCTTTAAAAGCTTANC'
      qual1 = 'ABBB@GGCCGGGGGFGGGGG<?E1?CFGGD/FFFGEGGGGEGGG1EEGG' \
         'G>GGGF?11:>CF/<<BEF0@F@DGGCG90<<FEFCF>GGGGGGGEF@FGGG0;' \
         ';;;FDD9;FCDGGEGEEGD###'
      qual2 = 'A3:ABGGGGGGGGGEGGGGGGGFGGGGGGGGGGGGGGCFGGFGGGFGGG' \
         'GGGGG1FGGGGGGGGG:BF1>EE8<EF1FC/CFBGEGGGGGFGGGGGEGGGGG>' \
         'CGGGGGGGGGGGGB:EEGGGGG'

      PEread = PER(read1, read2, qual1, qual2)
      PEread.identify_sample(ms)

      self.assertEqual(PEread.sample, '1t1')


   def test_identify_sample_fail(self):
      # Test setup (taken from segment 1).
      cst = 'AGGTTTGGATCAGGATTTGCGCCTTTGGAT'
      fwd = 'TTGCTCTCGGTCAAGCTTTTAAA'
      rev = 'ACAACTGCGTACGGCCTGTTCGA'
      samples = {
         (20,20): '1t0',
         (21,21): '1t1',
         (22,22): '1t14',
         (23,23): '2t0',
         (23,21): '2t1',
         (20,23): '2t14',
      }

      ms = preprocess.MultiplexSpecifications(cst, fwd, rev, samples)

      # Use 'PER' for short.
      PER = preprocess.PairedEndRead

      # The are already in the proper orientation.
      read1 = 'GCTCTCGGTCAAGCTTTTAAAAAGGCCATAGGGAATGTCAAAGGTATAG' \
         'CAAGGTTTGGATCAGGATTGGCGACTTTGGGTGAAGCACTTTCTAGATCAGTTG' \
         'TTGACTTGTCGAACAGGCCGTA'
      read2 = 'ACTGCGTACGGCCTGTTCGACAAGTCGACTCCTGACCTGGAAAGCGCTT' \
         'CATCCACTGGCGCAAATCCTGATCCAAACATTTCTATCCCTTTGATGGCTCCAA' \
         'TGGCATCTTTAAAAGCTTANCC'
      qual1 = 'ABBB@GGCCGGGGGFGGGGG<?E1?CFGGD/FFFGEGGGGEGGG1EEGG' \
         'G>GGGF?11:>CF/<<BEF0@F@DGGCG90<<FEFCF>GGGGGGGEF@FGGG0;' \
         ';;;FDD9;FCDGGEGEEGD###'
      qual2 = 'A3:ABGGGGGGGGGEGGGGGGGFGGGGGGGGGGGGGGCFGGFGGGFGGG' \
         'GGGGG1FGGGGGGGGG:BF1>EE8<EF1FC/CFBGEGGGGGFGGGGGEGGGGG>' \
         'CGGGGGGGGGGGGB:EEGGGGG'

      PEread = PER(read1, read2, qual1, qual2)
      with self.assertRaises(preprocess.NoSample):
         PEread.identify_sample(ms)


   def test_merge(self):
      # Test setup (taken from segment 1).
      cst = 'AGGTTTGGATCAGGATTTGCGCCTTTGGAT'
      fwd = 'TTGCTCTCGGTCAAGCTTTTAAA'
      rev = 'ACAACTGCGTACGGCCTGTTCGA'
      samples = {
         (20,20): '1t0',
         (21,21): '1t1',
         (22,22): '1t14',
         (23,23): '2t0',
         (23,21): '2t1',
         (20,23): '2t14',
         (21,20): 'XXX',
      }

      ms = preprocess.MultiplexSpecifications(cst, fwd, rev, samples)
      info = preprocess.LaneInfo()

      # Use 'PER' for short.
      PER = preprocess.PairedEndRead

      # Standard case without trimming.
      read1 = 'GCTCTCGGTCAAGCTTTTAAAAAGGCCATAGGGAATGTCAAAGGTATAG' \
         'CAAGGTTTGGATCAGGATTTGCGCCTTTGGATGAAGCACTTTCTAGATCAGTTG' \
         'TTGACTTGTCGAACAGGCCGTA'
      read2 = 'ACTGCGTACGGCCTGTTCGACAAGTCGACTCCTGACCTGGAAAGCGCTT' \
         'CATCCACTGGCGCAAATCCTGATCCAAACATTTCTATCCCTTTGATGGCTCCAA' \
         'TGGCATCTTTAAAAGCTTANCC'
      qual1 = 'A3:ABGGGGGGGGGEGGGGGGGFGGGGGGGGGGGGGGCFGGFGGGFGGG' \
         'GGGGG1FGGGGGGGGG:BF1>EE8<EF1FC/CFBGEGGGGGFGGGGGEGGGGG>' \
         'CGGGGGGGGGGGGB:EEG#GGG'
      qual2 = 'ABBB@GGCCGGGGGFGGGGG<?E1?CFGGD/FFFGEGGGGEGGG1EEGG' \
         'G>GGGF?11:>CF/<<BEF0@F@DGGCG90<<FEFCF>GGGGGGGEF@FGGG0;' \
         ';;;FDD9;FCDGGEGEEGG###'

      cs = 'GCTCTCGGTTAAGCTTTTAAAAAGGCCATAGGAAATGTCAAAGGTATAGCAA' \
         'GGTTTGGATCAGGATTTGCGCCTTTGGATGAAGCACTTTCCAGATCAGTTGTTG' \
         'ACTTGTCGAACAGGCCGTACGCAGT'

      PEread = PER(read1, read2, qual1, qual2)
      PEread.merge(ms, info, trim=False)

      self.assertEqual(PEread.cs, cs)
      self.assertEqual(info.nttot, 118)
      self.assertEqual(info.nterr, 20)

      # Swap the reads to check that there is no asymetry (they
      # will be oriented during the call).

      PEread = PER(read2, read1, qual2, qual1)
      PEread.merge(ms, info, trim=False)

      self.assertEqual(PEread.cs, cs)
      # Nucleotides are cumulative.
      self.assertEqual(info.nttot, 236)
      self.assertEqual(info.nterr, 40)

      # Now with trimming.
      cs = 'AAGGCCATAGGAAATGTCAAAGGTATAGCAAGGTTTGGATCAGGATTTGCGC' \
         'CTTTGGATGAAGCACTTTCCAGATCAGTTGTTGACTTG'

      PEread = PER(read1, read2, qual1, qual2)
      PEread.identify_sample(ms)
      self.assertEqual(PEread.sample, 'XXX')

      PEread.merge(ms, info, trim=True)
      self.assertEqual(PEread.cs, cs)
      self.assertEqual(info.nttot, 354)
      self.assertEqual(info.nterr, 60)


      # A case modified by removing 13 nucleotides in each side.
      # This makes a negative shift between the reads (i.e. read1
      # is right of read2 in the alignment).
      read1 = 'CTTTTAAAAAGGCCATAGGGAATGTCAAAGGTATAGCAAGGTTTGGATC' \
         'AGGATTTGCGCCTTTGGATGAAGCACTTTCTAGATCAGTTGTTGACTTGTCGAA' \
         'CAGGCCGTA'
      read2 = 'TGTTCGACAAGTCGACTCCTGACCTGGAAAGCGCTTCATCCACTGGCGC' \
         'AAATCCTGATCCAAACATTTCTATCCCTTTGATGGCTCCAATGGCATCTTTAAA' \
         'AGCTTANCC'
      qual1 = 'GEGGGGGGGFGGGGGGGGGGGGGGCFGGFGGGFGGGGGGGG1FGGGGGG' \
         'GGG:BF1>EE8<EF1FC/CFBGEGGGGGFGGGGGEGGGGG>CGGGGGGGGGGGG' \
         'B:EEGGGGG'
      qual2 = 'GFGGGGG<?E1?CFGGD/FFFGEGGGGEGGG1EEGGG>GGGF?11:>CF' \
         '/<<BEF0@F@DGGCG90<<FEFCF>GGGGGGGEF@FGGG0;;;;FDD9;FCDGG' \
         'EGEEGD###'

      PEread = PER(read1, read2, qual1, qual2)
      PEread.merge(ms, info, trim=False)

      cs = 'CTTTTAAAAAGGCCATAGGAAATGTCAAAGGTATAGCAAGGTTTGGATCAGG' \
         'ATTTGCGCCTTTGGATGAAGCACTTTCCAGATCAGTTGTTGACTTGTCGAACA'

      self.assertEqual(PEread.cs, cs)
      self.assertEqual(info.nttot, 459)
      self.assertEqual(info.nterr, 79)



   def test_merge_fail(self):
      # Test setup (taken from segment 1).
      cst = 'AGGTTTGGATCAGGATTTGCGCCTTTGGAT'
      fwd = 'TTGCTCTCGGTCAAGCTTTTAAA'
      rev = 'ACAACTGCGTACGGCCTGTTCGA'
      samples = {
         (20,20): '1t0',
         (21,21): '1t1',
         (22,22): '1t14',
         (23,23): '2t0',
         (23,21): '2t1',
         (20,23): '2t14',
      }

      ms = preprocess.MultiplexSpecifications(cst, fwd, rev, samples)
      info = preprocess.LaneInfo()

      # Use 'PER' for short.
      PER = preprocess.PairedEndRead

      # The central part is not present on the second read.
      read1 = 'GCTCTCGGTCAAGCTTTTAAAAAGGCCATAGGGAATGTCAAAGGTATAG' \
         'CAAGGTTTGGATCAGGATTGGCGACTTTGGGTGAAGCACTTTCTAGATCAGTTG' \
         'TTGACTTGTCGAACAGGCCGTA'
      read2 = 'ACTGCGTACGGAGTGCGCTAGGTACGTAGCGGGCGCATGAGCTGAGCGT' \
         'CCCCCCGCTAGGACTGAGAGATCGCGATTATATATAGCGCGATCGGCATGCGCG' \
         'CGGCGCGCGCTAAGATGCGGGC'
      qual1 = 'ABBB@GGCCGGGGGFGGGGG<?E1?CFGGD/FFFGEGGGGEGGG1EEGG' \
         'G>GGGF?11:>CF/<<BEF0@F@DGGCG90<<FEFCF>GGGGGGGEF@FGGG0;' \
         ';;;FDD9;FCDGGEGEEGD###'
      qual2 = 'A3:ABGGGGGGGGGEGGGGGGGFGGGGGGGGGGGGGGCFGGFGGGFGGG' \
         'GGGGG1FGGGGGGGGG:BF1>EE8<EF1FC/CFBGEGGGGGFGGGGGEGGGGG>' \
         'CGGGGGGGGGGGGB:EEGGGGG'

      PEread = PER(read1, read2, qual1, qual2)
      with self.assertRaises(preprocess.BadRead):
         PEread.merge(ms, info)

      # The sequences have an unresolved N.
      read1 = 'GCTCTCNGTCAAGCTTTTAAAAAGGCCATAGGGAATGTCAAAGGTATAG' \
         'CAAGGTTTGGATCAGGATTTGCGCCTTTGGATGAAGCACTTTCTAGATCAGTTG' \
         'TTGACTTGTCGAACAGGCCGTA'
      read2 = 'ACTGCGTACGGCCTGTTCGACAAGTCGACTCCTGACCTGGAAAGCGCTT' \
         'CATCCACTGGCGCAAATCCTGATCCAAACATTTCTATCCCTTTGATGGCTCCAA' \
         'TGGCATCTTTAAAAGCTTANCN'
      qual1 = 'A3:ABGGGGGGGGGEGGGGGGGFGGGGGGGGGGGGGGCFGGFGGGFGGG' \
         'GGGGG1FGGGGGGGGG:BF1>EE8<EF1FC/CFBGEGGGGGFGGGGGEGGGGG>' \
         'CGGGGGGGGGGGGB:EEGGGGG'
      qual2 = 'ABBB@GGCCGGGGGFGGGGG<?E1?CFGGD/FFFGEGGGGEGGG1EEGG' \
         'G>GGGF?11:>CF/<<BEF0@F@DGGCG90<<FEFCF>GGGGGGGEF@FGGG0;' \
         ';;;FDD9;FCDGGEGEEGD###'

      PEread = PER(read1, read2, qual1, qual2)
      with self.assertRaises(preprocess.BadRead):
         PEread.merge(ms, info)

if __name__ == '__main__':
   unittest.main()
