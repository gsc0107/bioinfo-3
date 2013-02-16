#!/usr/bin/python

# codes in "bioinfo" are under the MIT License
# Copyright (c) 2013 Jean-Etienne Morlighem <jem.nvnt@gmail.com>
# https://github.com/jem-gh/bioinfo



from bioinfo import *

import unittest
import random           # for "test_countNucleotides", "test_transcribe", 
                        # "test_GCcontent", "mutateSeq"
import re               # for test_countNucleotides, test_transcribe



def mutateSeq(seq, num_mut, type_nuc="DNA"):
    """ introduce x number of random point mutation(s) (num_mut) in a DNA or 
        RNA sequence. 
        'seq' can be a list or a string and a list/string will be returned 
        respectively. """
    
    assert num_mut >= 0
    
    if not seq: # empty string/list
        print "mutateSeq: no sequence provided"
        return None
    
    seq_length = len(seq)
    
    if seq_length < num_mut:
        print "mutateSeq: number of mutations higher than length of sequence "
        return False
    
    if "t" in seq and "u"in seq:
        print "mutateSeq: T and U found in sequence"
        return False
    
    if isinstance(seq, str):
        isString = True
        seq = [x for x in seq]
    elif isinstance(seq, list):
        isString = False
        seq = list(seq)
    
    pos_mutated = []
    
    if "u" in seq or (type_nuc=="RNA" and "t" not in seq):
        nuc_dic = "acgu"
    else:
        nuc_dic = "acgt"
    
    while num_mut > 0:
        pos = random.randrange(0, seq_length)
        
        if pos not in pos_mutated:
            nuc_ori = nuc_mut = seq[pos]
            
            while nuc_mut == nuc_ori:
                nuc_mut = random.choice(nuc_dic)
            
            seq[pos] = nuc_mut
            
            pos_mutated.append(pos)
            num_mut -= 1
    
    assert seq_length == len(seq)
    
    if isString:
        seq = ''.join(seq)
    
    return seq



class TestBioinfo(unittest.TestCase):
    
    #@unittest.skip("skip ... test_countNucleotides")
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
    
    
    #@unittest.skip("skip ... test_transcribe")
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
    
    
    #@unittest.skip("skip ... test_GCcontent")
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
    
    
    #@unittest.skip("skip ... test_complement")
    def test_complement(self):
        self.assertIs(complement(""), None)
        self.assertIs(complement("acgtu"), False)
        
        self.assertRaises(AssertionError, complement, (178439246))
        self.assertRaises(AssertionError, complement, ('"@ not a seq'))
        
        for d in [("acgt", "tgca"), ("acgu", "ugca")]:
            seq, comp = "", ""
            type_nuc = "DNA" if d == ("acgt", "tgca") else "RNA"
            for length in [10, 100, 1000, 10000]:
                for n in range(length):
                    nuc = random.randint(0,3)
                    seq  += d[0][nuc]
                    comp += d[1][nuc]
                self.assertEqual(complement(seq, type_nuc), comp)
        
        print "complement test ... OK"
    
    
    #@unittest.skip("skip ... test_reverse_complement")
    def test_reverse_complement(self):
        self.assertIs(reverse_complement(""), None)
        self.assertIs(reverse_complement("acgtu"), False)
        
        self.assertRaises(AssertionError, reverse_complement, (-356))
        self.assertRaises(AssertionError, reverse_complement, ('"@ not a seq'))
        
        for d in [("acgt", "tgca"), ("acgu", "ugca")]:
            seq, revcomp = "", ""
            type_nuc = "DNA" if d == ("acgt", "tgca") else "RNA"
            for length in [10, 100, 1000, 10000]:
                for n in range(length):
                    nuc = random.randint(0,3)
                    seq     += d[0][nuc]
                    revcomp  = d[1][nuc] + revcomp
                self.assertEqual(reverse_complement(seq, type_nuc), revcomp)
        
        print "reverse_complement test ... OK"
    
    
    #@unittest.skip("skip ... test_hammingDistance")
    def test_hammingDistance(self):
        for n in [("a",""), ("","a"), ("","")]:
            self.assertIs(hammingDistance(*n), None)
        self.assertRaises(AssertionError, hammingDistance, *("ttuu", "acgt"))
        
        self.assertRaises(AssertionError, hammingDistance, *(-356, "acgt"))
        self.assertRaises(AssertionError, hammingDistance, *('"@ not a seq', 'acg'))
        
        self.assertIs(hammingDistance("aaa", "aa"), False)
        
        for d in ["acgt", "acgu"]:
            seq = ''.join(random.choice(d) for n in range(100))
            
            for errors in range(101):
                seq_mutated = mutateSeq(seq, errors)
                self.assertEqual(hammingDistance(seq, seq_mutated), errors)
        
        print "hammingDistance test ... OK"



if __name__ == "__main__":
    unittest.main()

