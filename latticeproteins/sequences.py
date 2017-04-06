#!/usr/bin/python
# Begin sequences.py
#---------------------------------------------------------------------------
"""Module for lattice protein sequences.

Written by Jesse Bloom, 2004."""
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
#---------------------------------------------------------------------------
def most_abundant(population):
    """Returns the most abundant sequence in a population.

    Call is: 'seq = MostAbundant(population)'
    'population' is a list of protein sequences.
    'seq' is returned as the sequence that is most abundant in 'population'.
        If there are several equally abundant sequences, just one of
        them is returned.
    """
    seq = population[0]
    n = population.count(seq)
    for seq2 in population[1 : ]:
        if seq2 != seq:
            if population.count(seq2) > n:
                n = population.count(seq2)
                seq = seq2
    return seq
#---------------------------------------------------------------------------
def pairwise_hamming_distances(seqlist):
    """Computes the pairwise Hamming distances between many sequences.

    Call is: 'dlist = pairwise_hamming_distances(seqlist)'
    'seqlist' is a list of sequences all of the same length.
    'dlist' is returned as a list of numbers representing the
        Hamming distances between all pairs of sequences."""
    if not (isinstance(seqlist, list) and len(seqlist) > 1):
        raise SequenceError("'seqlist' is not a list of at least 2 entries.")
    length = len(seqlist[0])
    dlist = []
    for i1 in range(len(seqlist)):
        seq1 = seqlist[i1]
        if len(seq1) != length:
            raise SequenceError("Invalid length sequence of %r." % seq1)
        for i2 in range(i1 + 1, len(seqlist)):
            seq2 = seqlist[i2]
            if len(seq2) != length:
                raise SequenceError("Invalid length sequence of %r." % seq1)
            dlist.append(hamming_distance(seq1, seq2))
    if len(dlist) != len(seqlist) * (len(seqlist) - 1) / 2:
        raise SequenceError("Incorrect number of distances.")
    return dlist
#---------------------------------------------------------------------------
def hamming_distance(seq1, seq2):
    """Returns the Hamming distance between two sequences.

    Call is: 'd = hamming_distance(seq1, seq2)'
    'seq1' and 'seq2' are two sequences of the same length.
    'd' is returned as the Hamming distance between these two sequences."""
    if len(seq1) != len(seq2):
        raise SequenceError("Sequences differ in length.")
    d = 0
    for i in range(len(seq1)):
        if seq1[i] != seq2[i]:
            d += 1
    return d
#---------------------------------------------------------------------------
def random_sequence(length):
    """Returns a random sequence of the specified length.

    Call is 's = random_sequence(length)'
    'length' is an integer >= 1.  Returns a sequence of length 'length'
        as a list of randomly chosen residues."""
    if not (isinstance(length, int) and length > 0):
        raise SequenceError("Invalid sequence length of %r." % length)
    s = [random.choice(_residues) for i in range(length)]
    return s
#------------------------------------------------------------------------
def mutate_sequence(seq, mutrate):
    """Mutates a protein sequence.

    Call is: 'seqnew = mutate_sequence(seq, mutrate)'
    'seq' is a protein sequence, specified as either a string or a list.
    Mutates each residue in 'seq' to some different residue with
        probability 'mutrate'.  So 'mutrate' is the per residue
        mutation rate.
    Returns the new sequence as the list 'seqnew'."""
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
#------------------------------------------------------------------------
def n_mutants(seq, nmutations, nsequences):
    """Returns sequences with a specified number of mutations.

    Call is: 'seqlist = n_mutants(seq, nmutations, nsequences)'
    'seq' is a string or list specifying the protein we wish to mutate.
    'nmutations' is the number of mutations each mutant of 'seq' should
        have.  It must be <= 'len(seq)' and > 0.
    'nsequences' is the number of mutant sequences to make.  It can be
        'ALL', in which case we make all possible mutants with 'nmutations',
        or it can be some positive integer in which case we make this
        many randomly chosen mutants with 'nmutations' mutations.
        'ALL' is only a valid option only when 'nmutations' is 1 or 2.
    The mutant sequences are returned in the list 'seqlist'.  Each entry
        in 'seqlist' is a list representing a mutant sequence."""
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
