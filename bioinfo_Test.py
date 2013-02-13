#!/usr/bin/python

# codes in "bioinfo" are under the MIT License
# Copyright (c) 2013 Jean-Etienne Morlighem <jem.nvnt@gmail.com>
# https://github.com/jem-gh/bioinfo



from bioinfo import *

import unittest
import random           # for test_countNucleotides, test_transcribe
import re               # for test_transcribe


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



if __name__ == "__main__":
    unittest.main()

