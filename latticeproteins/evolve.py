import numpy as np
import np2d

from latticeproteins.sequences import find_differences, _residues

def fixation(fitness1, fitness2, N=10e8, *args, **kwargs):
    """ Simple fixation probability between two organism with fitnesses 1 and 2.
    Note that N is the effective population size.
    .. math::
        p_{\\text{fixation}} = \\frac{1 - e^{-N \\frac{f_2-f_1}{f1}}}{1 - e^{-\\frac{f_2-f_1}{f1}}}
    """
    sij = (fitness2 - fitness1)/abs(fitness1)
    # Check the value of denominator
    denominator = 1 - np.exp(-N * sij)
    numerator = 1 - np.exp(- sij)
    # Calculate the fixation probability
    fixation = numerator / denominator
    fixation = np.nan_to_num(fixation)
    return fixation

def monte_carlo_fixation_walk(seq, lattice, selected_trait="fracfolded", fix_threshold=0.5, max_mutations=15, target=None):
    """Use Monte Carlo method to walk

    Parameters
    ----------
    seq : str
        seq
    lattice : LatticeThermodynamics object
        Lattice protein calculator
    selected_trait : str
        The trait to select.
    max_mutations : int (default = 15)
        Max number of mutations to make in the walk.
    target : str
        selected lattice target conformation. If None, the lattice will
        fold to the natural native conformation.
    """
    length = len(seq)
    fitness_method = getattr(lattice, selected_trait)
    fitness0 = fitness_method(seq, target=target)
    finished = False
    mutant = list(seq[:])
    path, fitness, probs = [seq], [fitness0], [0]
    # Monte Carlo move.
    m = 0
    while finished is False and m < max_mutations:
        # Construct grid of all stabilities of all amino acids at all sites
        AA_grid = np.array([_residues]*length)
        fits = np.zeros(AA_grid.shape, dtype=float)
        for (i,j), AA in np.ndenumerate(AA_grid):
            seq1 = mutant[:]
            seq1[i] = AA_grid[i,j]
            fits[i,j] = fitness_method(seq1, target=target)

        # Calculate fitness for all neighbors in sequence space
        fix = fixation(fitness0, fits)

        # Normalize
        denom = fix.sum()
        p = fix / denom

        # Sample moves
        mutation, indices = np2d.random.choice(AA_grid, p=p)
        site = indices[0,0]
        AA = indices[0,1]
        mutant[site] = mutation

        # Update our trajectory
        path.append("".join(mutant))
        fitness.append(fits[site, AA])
        probs.append(fix[site, AA])
        fitness0 = fits[site, AA]

        # Check criteria to kill the trajectory
        # If the total probability of a mutation fixing is < 5%,
        # then kill the loop.
        if denom < fix_threshold:
            finished = True
        m += 1

    return path, fitness, probs
