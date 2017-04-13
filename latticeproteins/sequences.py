#!/usr/bin/python
# Begin sequences.py
#---------------------------------------------------------------------------
"""
Originally written by Jesse Bloom, 2004.

Updated by Zach Sailer, 2017."""
#---------------------------------------------------------------------------
import random, shelve, os
#---------------------------------------------------------------------------
class SequenceError(Exception):
    """Error with a lattice protein sequence."""
    pass
#---------------------------------------------------------------------------
# codes for all residues
_residues = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P',
'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
assert len(_residues) == 20

def hamming_distance(seq1, seq2):
    """Returns the Hamming distance between two sequences."""
    if len(seq1) != len(seq2):
        raise SequenceError("Sequences differ in length.")
    d = 0
    for i in range(len(seq1)):
        if seq1[i] != seq2[i]:
            d += 1
    return d

def find_differences(s1, s2):
    """Return the index of differences between two sequences."""
    indices = list()
    for i in range(len(s1)):
        if s1[i] != s2[i]:
            indices.append(i)
    return indices

def random_sequence(length):
    """Returns a random sequence of the specified length."""
    if not (isinstance(length, int) and length > 0):
        raise SequenceError("Invalid sequence length of %r." % length)
    s = [random.choice(_residues) for i in range(length)]
    return s

def mutate_sequence(seq, mutrate):
    """Mutates a protein sequence.

    Parameters
    ----------
    seq :
        is a protein sequence, specified as either a string or a list.
    mutrate :
        Mutates each residue in 'seq' to some different residue with
        probability 'mutrate'.  So 'mutrate' is the per residue
        mutation rate.

    Returns
    -------
    newseq :
        the new sequence as a list."""
    mutated = False
    for ires in range(len(seq)):
        if random.random() < mutrate:
            if not mutated:
                mutated = True
                newseq = list(seq)
            newres = random.choice(_residues)
            while newres == seq[ires]:
                newres = random.choice(_residues)
            newseq[ires] = newres
    if mutated:
        return newseq
    else:
        return seq

def n_mutants(seq, nmutations, nsequences):
    """Returns sequences with a specified number of mutations.

    Parameters
    ----------
    seq :
        is a string or list specifying the protein we wish to mutate.
    nmutations :
        is the number of mutations each mutant of 'seq' should
        have.  It must be <= 'len(seq)' and > 0.
    nsequences :
        is the number of mutant sequences to make.  It can be
        'ALL', in which case we make all possible mutants with 'nmutations',
        or it can be some positive integer in which case we make this
        many randomly chosen mutants with 'nmutations' mutations.
        'ALL' is only a valid option only when 'nmutations' is 1 or 2.

    Return
    ------
    seqlist : list
        List of mutant sequences n mutations away.
    """
    if not (0 < nmutations <= len(seq)):
        raise SequenceError("Invalid 'nmutations' of %r." % nmutations)
    seqlist = []
    if nsequences == 'ALL':
        if nmutations == 1:
            for ires in range(len(seq)):
                for mutres in _residues:
                    if mutres != seq[ires]:
                        newseq = list(seq)
                        newseq[ires] = mutres
                        seqlist.append(newseq)
        elif nmutations == 2:
            for ires in range(len(seq)):
                for imutres in _residues:
                    if imutres != seq[ires]:
                        for jres in range(ires + 1, len(seq)):
                            for jmutres in _residues:
                                if jmutres != seq[jres]:
                                    newseq = list(seq)
                                    newseq[ires] = imutres
                                    newseq[jres] = jmutres
                                    seqlist.append(newseq)
        else:
            raise SequenceError("'nsequences' cannot be 'ALL' when 'nmutations' is %r." % nmutations)
    elif isinstance(nsequences, int) and nsequences > 0:
        for imutant in range(nsequences):
            newseq = list(seq)
            for imut in range(nmutations):
                ires = random.choice(range(len(seq)))
                while newseq[ires] != seq[ires]:
                    ires = random.choice(range(len(seq)))
                mutres = random.choice(_residues)
                while mutres == seq[ires]:
                    mutres = random.choice(_residues)
                newseq[ires] = mutres
            seqlist.append(newseq)
    else:
        raise SequenceError("Invalid 'nsequences' of %r." % nsequences)
    return seqlist
