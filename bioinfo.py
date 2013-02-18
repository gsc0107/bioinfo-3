#!/usr/bin/python

# codes in "bioinfo" are under the MIT License
# Copyright (c) 2013 Jean-Etienne Morlighem <jem.nvnt@gmail.com>
# https://github.com/jem-gh/bioinfo



from __future__ import division

import re           # for "transcribe" (assertion)
from string import maketrans        # for "complement"

from bioinfoLibrary import RNA_CODON_TABLE



def countNucleotides(seq):
    """ simple nucleotide counting function for DNA or RNA sequences without 
        gap nor nucleotide ambiguity """
    
    assert isinstance(seq, str)
    
    if not seq: # empty string
        print "countNucleotides: no sequence provided"
        return None
    
    seq = seq.lower()
    
    a = seq.count('a')
    c = seq.count('c')
    g = seq.count('g')
    t = seq.count('t')
    u = seq.count('u')
    
    if t and u:
        print "countNucleotides: both T and U found in the sequence"
        return False
    
    assert (a + c + g + t == len(seq)) or (a + c + g + u == len(seq))
    
    print "countNucleotides result:"
    print "{:>5} {:>5} {:>5} {:>5} {:>5}".format("A","C","G","T","U")
    print "{:>5} {:>5} {:>5} {:>5} {:>5}".format(a,c,g,t,u)
    
    return a, c, g, t, u



def transcribe(seq):
    """ Transcribed a DNA sequence into its corresponding RNA sequence """
    
    assert isinstance(seq, str)
    
    if not seq: # empty string
        print "transcribe: no sequence provided"
        return None
    
    assert countNucleotides(seq)
    
    seq = seq.lower()
    
    if "u" in seq:
        print "transcribe: U already in the sequence"
        return False
    
    result = seq.replace("t", "u")
    
    assert ([pos.start() for pos in re.finditer("t", seq)] == 
            [pos.start() for pos in re.finditer("u", result)])
    
    print "transcribe result:"
    print result
    
    return result



def GCcontent(seq):
    """ return the GC content in % of a nucleotide sequence """
    
    assert isinstance(seq, str)
    
    if not seq: # empty string
        print "GCcontent: no sequence provided"
        return None
    
    try:
        a,c,g,t,u = countNucleotides(seq)
    except TypeError:
        print "GCcontent: sequence not valid"
        return False
    
    result = (c + g) / (a + c + g + t + u) * 100
    
    # due to the Python floating-point approximation issue, assertion is made 
    # by comparing results to an acceptable error margin
    assert abs(result + ((a+t+u)/(a+c+g+t+u)*100) - 100) <= 0.00000000001
    
    print "GCcontent result:"
    print result
    
    return result



def complement(seq, type_nuc="DNA"):
    """ return the complement of a DNA or RNA sequence """
    
    assert isinstance(seq, str)
    
    if not seq: # empty string
        print "complement: no sequence provided"
        return None
    
    try:
        a,c,g,t,u = countNucleotides(seq)
    except TypeError:
        print "complement: sequence not valid"
        return False
    
    in_seq, out_seq  = ("acgu", "ugca") if (u or type_nuc=="RNA") else ("acgt", "tgca")
    result = seq.lower().translate(maketrans(in_seq, out_seq))
    
    assert countNucleotides(result) in [(t, g, c, a, u), (u, g, c, t, a)]
    
    print "complement result:"
    print result
    
    return result



def reverse_complement(seq, type_nuc="DNA"):
    """ return the reverse complement of a DNA or RNA sequence """
    
    assert isinstance(seq, str)
    
    if not seq: # empty string
        print "reverse_complement: no sequence provided"
        return None
    
    if "u" in seq.lower():
        type_nuc = "RNA"
    
    result_c = complement(seq, type_nuc)
    if result_c:
        result_rc = result_c[::-1]
    else:
        print "reverse_complement: sequence not valid"
        return False
    
    assert complement(result_rc, type_nuc)[::-1] == seq.lower()
    
    print "reverse_complement result:"
    print result_rc
    
    return result_rc



def hammingDistance(seq1, seq2):
    """ return the Hamming distance (number of differences) between two DNA or 
        RNA sequences of the same length """
    
    assert isinstance(seq1, str) and isinstance(seq2, str)
    
    if not (seq1 and seq2): # empty strings
        print "hammingDistance: two sequences needed"
        return None
    
    assert countNucleotides(seq1) and countNucleotides(seq2)
    
    if len(seq1) != len(seq2):
        print "hammingDistance: sequences do not have the same length"
        return False
    
    result = 0
    
    for i in range(len(seq1)):
        if seq1[i] != seq2[i]:
            result += 1
    
    assert 0 <= result <= len(seq1)
    
    print "hammingDistance result:"
    print result
    
    return result



def findMotif(seq, motif):
    """ return all beginning positions of a motif in a sequence """
    
    assert isinstance(seq, str) and isinstance(motif, str)
    
    if not (seq and motif): # empty strings
        print "findMotif: two sequences needed"
        return None
    
    if len(motif) > len(seq):
        print "findMotif: motif length longer than the sequence length"
        return False
    
    seq, motif = seq.lower(), motif.lower()
    
    positions = []
    current_pos = 0
    
    while current_pos != -1:
        current_pos = seq.find(motif, current_pos)
        if current_pos != -1:
            positions.append(current_pos+1)
            current_pos += 1
    
    assert len(positions) <= (len(seq) - len(motif) + 1)
    
    print "findMotif result:"
    print None if not positions else ' '.join(str(p) for p in positions)
    
    return positions



def translation(seq):
    """ return the protein sequence encoded by an RNA sequence (seq) """
    
    assert isinstance(seq, str)
    
    if not len(seq) >= 3: # sequence shorter than a codon
        print "translation: an RNA sequence of at least 3 nucleotides is required"
        return None
    
    try:
        a,c,g,t,u = countNucleotides(seq)
    except TypeError:
        print "translation: not a valid RNA sequence"
        return False
    if t:
        print "translation: sequence provided is DNA instead of RNA"
        return False
    
    length = len(seq)
    seq = seq.upper()
    
    result = ""
    
    for n in range(0, length, 3):
        try:
            aa = RNA_CODON_TABLE[ seq[n:n+3] ]
        except KeyError:     # occurs if remaining RNA seq is less than 3 nucleotides 
            break
        
        if aa == "*":
            break
        else:
            result += aa
    
    assert len(result) <= length / 3
    
    print "translation result:"
    print result if result else "No protein obtained from the RNA sequence"
    
    return result



