#!/usr/bin/python

# codes in "bioinfo" are under the MIT License
# Copyright (c) 2013 Jean-Etienne Morlighem <jem.nvnt@gmail.com>
# https://github.com/jem-gh/bioinfo



from bioinfo import *

import unittest
import random           # for test_countNucleotides, test_transcribe, test_GCcontent
import re               # for test_countNucleotides, test_transcribe



class TestBioinfo(unittest.TestCase):
    
    def test_countNucleotides(self):
        self.assertIs(countNucleotides(""), None)
        self.assertIs(countNucleotides("tu"), False)
        
        self.assertRaises(AssertionError, countNucleotides, (123))
        self.assertRaises(AssertionError, countNucleotides, ("not a correct seq"))
        
        for nucleotides in ["acgt", "ACGT", "acgu", "ACGU"]:
            for length in [1, 10, 100, 1000, 10000]:
                sequence = ""
                for i in range(length):
                    sequence += random.choice(nucleotides)
                a, c, g, t, u = countNucleotides(sequence)
                self.assertEqual(a, len(re.findall('a', sequence.lower())))
                self.assertEqual(c, len(re.findall('c', sequence.lower())))
                self.assertEqual(g, len(re.findall('g', sequence.lower())))
                self.assertEqual(t, len(re.findall('t', sequence.lower())))
                self.assertEqual(u, len(re.findall('u', sequence.lower())))
        
        print "countNucleotides test ... OK"
    
    
    def test_transcribe(self):
        self.assertIs(transcribe(""), None)
        self.assertIs(transcribe("acgu"), False)
        
        self.assertRaises(AssertionError, countNucleotides, (067))
        self.assertRaises(AssertionError, countNucleotides, ("> !! not a seq"))
        
        for nucleotides in ["acgt", "ACGT"]:
            for length in [1, 10, 100, 1000, 10000]:
                dna = ""
                for i in range(length):
                    dna += random.choice(nucleotides)
                rna = transcribe(dna)
                posT = [pos.start() for pos in re.finditer("t", dna.lower())]
                posU = [pos.start() for pos in re.finditer("u", rna)]
                self.assertEqual(posT, posU)
        
        print "transcribe test ... OK"
    
    
    def test_GCcontent(self):
        self.assertIs(GCcontent(""), None)
        self.assertIs(GCcontent("acgtu"), False)
        
        self.assertRaises(AssertionError, GCcontent, (987))
        self.assertRaises(AssertionError, GCcontent, ("> !! not a seq"))
        
        for n in range(0, 1001):
            sequence = [random.choice("cg") for i in range(n)]
            sequence += [random.choice("at") for i in range(1000-n)]
            random.shuffle(sequence)
            sequence = ''.join([x for x in sequence])
            self.assertAlmostEqual(GCcontent(sequence), n / 10.)
        
        print "GCcontent test ... OK"
    
    
    def test_complement(self):
        self.assertIs(complement(""), None)
        self.assertIs(complement("acgtu"), False)
        
        self.assertRaises(AssertionError, complement, (178439246))
        self.assertRaises(AssertionError, complement, ('"@ not a seq'))
        
        for d in [("acgt", "tgca"), ("acgu", "ugca")]:
            seq, comp = "", ""
            for length in [10, 100, 1000, 10000]:
                for n in range(length):
                    nuc = random.randint(0,3)
                    seq  += d[0][nuc]
                    comp += d[1][nuc]
                self.assertEqual(complement(seq), comp)
        
        print "complement test ... OK"
    
    
    def test_reverse_complement(self):
        self.assertIs(reverse_complement(""), None)
        self.assertIs(reverse_complement("acgtu"), False)
        
        self.assertRaises(AssertionError, reverse_complement, (-356))
        self.assertRaises(AssertionError, reverse_complement, ('"@ not a seq'))
        
        for d in [("acgt", "tgca"), ("acgu", "ugca")]:
            seq, revcomp = "", ""
            for length in [10, 100, 1000, 10000]:
                for n in range(length):
                    nuc = random.randint(0,3)
                    seq     += d[0][nuc]
                    revcomp  = d[1][nuc] + revcomp
                self.assertEqual(reverse_complement(seq), revcomp)
        
        print "reverse_complement test ... OK"



if __name__ == "__main__":
    unittest.main()

