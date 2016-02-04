# -*- coding: utf-8 -*-
"""
Mini project 1: Gene Finder

@author: Kristyn Walker

"""

import random
from amino_acids import aa, codons, aa_table   # you may find these useful
from load import load_seq


def shuffle_string(s):
    """Shuffles the characters in the input string
        NOTE: this is a helper function, you do not
        have to modify this in any way """
    return ''.join(random.sample(s, len(s)))

# YOU WILL START YOUR IMPLEMENTATION FROM HERE DOWN ###


def get_complement(nucleotide):
    """ Returns the complementary nucleotide

        nucleotide: a nucleotide (A, C, G, or T) represented as a string
        returns: the complementary nucleotide
    >>> get_complement('A')
    'T'
    >>> get_complement('C')
    'G'
    >>> get_complement('T')
    'A'
    >>> get_complement('G')
    'C'
    """
    if nucleotide == 'A':
        return 'T'
    elif nucleotide == 'C':
        return 'G'
    elif nucleotide == 'T':
        return 'A'
    else: #if nucleotide = G
        return 'C'
   


def get_reverse_complement(dna):
    """ Computes the reverse complementary sequence of DNA for the specfied DNA
        sequence

        dna: a DNA sequence represented as a string
        returns: the reverse complementary DNA sequence represented as a string
    >>> get_reverse_complement("ATGCCCGCTTT")
    'AAAGCGGGCAT'
    >>> get_reverse_complement("CCGCGTTCA")
    'TGAACGCGG'
    >>> get_reverse_complement("ATGC")
    'GCAT'
    """
    sequence = ""
    for letter in dna: 
        complement = get_complement(letter)
        sequence = sequence + complement  
    reverse_sequence = sequence[::-1]
    return reverse_sequence


def rest_of_ORF(dna):
    """ Takes a DNA sequence that is assumed to begin with a start
        codon and returns the sequence up to but not including the
        first in frame stop codon.  If there is no in frame stop codon,
        returns the whole string.

        dna: a DNA sequence
        returns: the open reading frame represented as a string
    >>> rest_of_ORF("ATGTGAA")
    'ATG'
    >>> rest_of_ORF("ATGAGATAGG")
    'ATGAGA'
    >>> rest_of_ORF("ATGTGA")
    'ATG'
    >>> rest_of_ORF("ATGTGCCC")
    'ATGTGCCC'
    """
    templist = [] #empty list
    for i in range(0,len(dna), 3): #breaks dna into index of 3
        templist.append(dna[i:i+3]) #slices dna into indexes 0-3,3-6, etc.. 

    codonindex = 0
    for codon in templist: #for each group of 3 in the list
        if codon == 'TAG' or codon =='TGA' or codon == 'TAA': #checks if index == to stop codon
            codonindex = templist.index(codon) #tracks which index is stopcodon
            break #once stop codon found, breaks out of if loop

    if codonindex != 0:
        newlist = dna[0:codonindex * 3] #adds codons to list
    else:
        newlist = dna
    return newlist


def find_all_ORFs_oneframe(dna):
    """ Finds all non-nested open reading frames in the given DNA
        sequence and returns them as a list.  This function should
        only find ORFs that are in the default frame of the sequence
        (i.e. they start on indices that are multiples of 3).
        By non-nested we mean that if an ORF occurs entirely within
        another ORF, it should not be included in the returned list of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_oneframe("ATGCATGAATGTAGATAGATGTGCCC")
    ['ATGCATGAATGTAGA', 'ATGTGCCC']
    """


    i = 0
    results = []
    while i<len(dna)-1:
        if dna[i:i+3] == 'ATG':
            temp = rest_of_ORF(dna[i:])
            results.append(temp)
            i = i+len(temp)+3
        else:
            i = i+3
    return results
   
def find_all_ORFs(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence in
        all 3 possible frames and returns them as a list.  By non-nested we
        mean that if an ORF occurs entirely within another ORF and they are
        both in the same frame, it should not be included in the returned list
        of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs

    >>> find_all_ORFs("ATGCATGAATGTAG")
    ['ATGCATGAATGTAG', 'ATGAATGTAG', 'ATG']
    """
    all_ORFs = []
    for i in range(3): #checking starting at index 0,1,2
        ORF_frame= find_all_ORFs_oneframe(dna[i:])
        all_ORFs.extend(ORF_frame)
    return all_ORFs



def find_all_ORFs_both_strands(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence on both
        strands.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_both_strands("ATGCGAATGTAGCATCAAA")
    ['ATGCGAATG', 'ATGCTACATTCGCAT']
    """
    both_ORFs = []
    forward_ORFs = find_all_ORFs(dna)
    reverse_dna = get_reverse_complement(dna)
    reverse_ORFs = find_all_ORFs(reverse_dna)
    both_ORFs = forward_ORFs + reverse_ORFs
    return both_ORFs 




def longest_ORF(dna):
    """ Finds the longest ORF on both strands of the specified DNA and returns it
        as a string
    >>> longest_ORF("ATGCGAATGTAGCATCAAA")
    'ATGCTACATTCGCAT'
    """
    max_strand_length = 0
    orf_long = ""
    for strand in find_all_ORFs_both_strands(dna):
        strand_length =len(strand)
        if strand_length > max_strand_length:
            max_strand_length = strand_length
            orf_long = strand
    return orf_long


def longest_ORF_noncoding(dna, num_trials):
    #return max length ORF
    """ Computes the maximum length of the longest ORF over num_trials shuffles
        of the specfied DNA sequence

        dna: a DNA sequence
        num_trials: the number of random shuffles
        returns: the maximum length longest ORF 
        >>> longest_ORF_noncoding
        """
    longest_length = 0
    longest_string = ''
    for i in range(num_trials):
        dna_string = shuffle_string(dna)    
        candidate =  longest_ORF(dna_string)
        if len(candidate) > longest_length:
            longest_string = candidate
            longest_length = len(candidate)
    return longest_string


def coding_strand_to_AA(dna):
    """ Computes the Protein encoded by a sequence of DNA.  This function
        does not check for start and stop codons (it assumes that the input
        DNA sequence represents an protein coding region).

        dna: a DNA sequence represented as a string
        returns: a string containing the sequence of amino acids encoded by the
                 the input DNA fragment

        >>> coding_strand_to_AA("ATGCGA")
        'MR'
        >>> coding_strand_to_AA("ATGCCCGCTT")
        'MPA'
    """
    amino = ''
    i=0
    while i+3 < len(dna) + 1:
        amino_acid = aa_table[dna[i:i+3]]
        amino += amino_acid
        i += 3
    return amino
        


def gene_finder(dna):
    """ Returns the amino acid sequences that are likely coded by the specified dna

        dna: a DNA sequence
        returns: a list of all amino acid sequences coded by the sequence dna.
    """
    amino_acid_sequence = []
    threshold = longest_ORF_noncoding(dna, 1500)
    for item in find_all_ORFs_both_strands(dna):
        if len(item) > len(threshold):
            amino_acid_sequence.append(coding_strand_to_AA(item))
    return amino_acid_sequence      

if __name__ == "__main__":
    import doctest
    dna = load_seq("./data/X73525.fa")
    print gene_finder(dna)
    #doctest.testmod()
    #print len(longest_ORF_noncoding(dna,1))
    #doctest.run_docstring_examples(coding_strand_to_AA, globals())
