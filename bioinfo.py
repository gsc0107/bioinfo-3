#!/usr/bin/python



def countNucleotides(seq):
    """ simple nucleotide counting function for DNA or RNA sequences without 
        gap nor nucleotide ambiguity """
    
    if not seq:
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


