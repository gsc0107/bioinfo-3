#!/usr/bin/python



from bioinfo import *

import unittest
import random           # for test_countNucleotides



class TestBioinfo(unittest.TestCase):
    
    def test_countNucleotides(self):
        self.assertIs(countNucleotides(""), None)
        self.assertIs(countNucleotides("tu"), False)
        
        for nucleotides in ["acgt", "ACGT", "acgu", "ACGU"]:
            for length in [1, 10, 100, 1000, 10000]:
                sequence = ""
                for i in range(length):
                    sequence += random.choice(nucleotides)
                a, c, g, t, u = countNucleotides(sequence)
                self.assertEqual(a, sequence.lower().count('a'))
                self.assertEqual(c, sequence.lower().count('c'))
                self.assertEqual(g, sequence.lower().count('g'))
                self.assertEqual(t, sequence.lower().count('t'))
                self.assertEqual(u, sequence.lower().count('u'))
        
        print "countNucleotides test ... OK"



if __name__ == "__main__":
    unittest.main()
