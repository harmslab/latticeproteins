#!/usr/bin/python
# Begin fitness.py
#---------------------------------------------------------------------
"""Module for calculating fitnesses of lattice protein sequences.

Written by Jesse Bloom, 2004."""
#----------------------------------------------------------------------
import os
import math, sys
import numpy as np

from . import conformations
from .interactions import miyazawa_jernigan
#----------------------------------------------------------------------
class ThermodynamicsError(Exception):
    """Error computing lattice protein thermodynamics."""

def fold_energy(sequence, conformation, interactions=miyazawa_jernigan):
    """Calculate the energy of the sequence with the given conformation.

    Parameters
    ----------
    sequence : str
        Amino acid sequence to fold.
    conformation : str
        Conformation according to latticemodel's conformations format (e.g. 'UDLLDRU')

    Returns
    ------
    energy : float
        energy of the conformation (sum of all contact energies)
    """
    contacts = lattice_contacts(sequence, conformation)
    energy = sum([interactions[c] for c in contacts])
    return energy

def lattice_contacts(sequence, conformation):
    """Find all contacts in conformation.

    Parameters
    ----------
    sequence : str
        Amino acid sequence to fold.
    conformation : str
        Conformation according to latticemodel's conformations format (e.g. 'UDLLDRU')

    Returns
    -------
    contacts : list
        list of contact pairs
    """
    sites = list(sequence)
    length = len(sites)
    try:
        moves = list(conformation)
    except TypeError:
        raise Exception("""Protein conformation is None; is there a native state? """)
    # build a coordinate system, note that odd rotation of intuitive coordinates
    # since we are working in numpy array grid.
    coordinates = {"U": [-1,0], "D":[1,0], "L":[0,-1], "R":[0,1]}
    grid = np.zeros((2*length+1, 2*length+1), dtype=str)
    x = y = round(length/2.0) # initial position on the grid is at the center of the 2d array
    grid[x,y] = sites[0]
    # move on grid, populate with amino acid at that site, and store all contacting neighbors.
    contacts = []
    for i in range(length-1):
        step = coordinates[moves[i]]
        x += step[0]
        y += step[1]
        grid[x,y] = sites[i+1]
        neighbors = [sites[i+1] + grid[x+c[0], y+c[1]] for c in coordinates.values()]
        contacts += [n for n in neighbors if n in miyazawa_jernigan]
    # subtract the contacts that have bonds between them.
    for i in range(1,length):
        try:
            contacts.remove(sequence[i-1:i+1])
        except ValueError:
            contacts.remove(sequence[i] + sequence[i-1])
    return contacts

#----------------------------------------------------------------------
class LatticeThermodynamics(object):
    """Attaches thermodynamic evaluators to a lattice protein conformation database.

    Parameters
    ----------
    temp : float
        the temperature at which the fitness is computed.
    conformations : conformations.Conformations object
        is the 'conformations.Conformations' object
        used to fold the protein sequences.  'conformations.Length()'
        specifies the length of the protein sequences that can be
        folded.
    dGdependence :
        specifies how the fitness depends on the free
        energy of folding of the sequence:
        * if 'dGdependence' is number, then any sequence with a
            free energy of folding <= 'dGdependence' has a fitness
            of one, and any sequence with a free energy of folding
            > 'dGdependence' has a fitness of 'nofitness'.
        * if 'dGdependence' is the string 'fracfolded' then the fitness
            of the sequence is the fraction of the sequences that
            will be folded at temperature 'temp' at equilibrium.
            If 'dG' is the free energy of folding, this fraction is:
            'f = 1 / (1 + exp(dG / temp)'
        * if 'dGdependence' is the string 'negstability', then the
            fitness of the sequence is just negative one times the
            stability of the sequence.  The negative one is so
            that more stable sequences have higher fitnesses.
    ligand :
        is an optional argument that is used if we are looking
        for a protein that binds a ligand.  By default, it is 'None'
        meaning that no ligand binding is considered.  If it is set
        to another value, it should be the 3-tuple '(ligand,
        ligandconf, stabcutoff)' where 'ligand' and 'ligandconf'
        are both strings describing a ligand as detailed in the
        documentation string for the 'conformations.BindLigand'
        method.  'stabcutoff' is a number specifying the
        stability cutoff for the protein to fold.  In this case,
        'dGdependence' no longer has any meaning.  The protein
        is folded according to the parameters 'temp' and
        'targets' as normal.  If the free energy of folding of
        the protein is > 'stabcutoff', then the returned fitness is
        zero.  If the free energy of folding of the protein is <=
        'stabcutoff', then the returned fitness 'exp(-be)' where
        'be' is the binding energy of the ligand to the protein
        in the folded conformation.
    """
    #------------------------------------------------------------------
    def __init__(self, temp, conformations, ligand=None):
        # Assign class instance variables and error check
        self._temp = temp
        if not (isinstance(self._temp, (int, float)) and temp > 0):
            raise ThermodynamicsError("Invalid 'temp' of %r." % temp)
        self._conformations = conformations
        self._ligand = ligand
        if ligand == None:
            pass
        elif isinstance(ligand, tuple) and len(ligand) == 3:
            (lig, ligconf, stabcutoff) = ligand
            if not (isinstance(lig, str) and isinstance(ligconf, str) and len(lig) == len(ligconf) + 1):
                raise ThermodynamicsError("%r does not specify a valid ligand." % ligand)
            if not (isinstance(stabcutoff, (int, float))):
                raise ThermodynamicsError( "Invalid 'stabcutoff' of %r." % stabcutoff)
        else:
            raise ThermodynamicsError("Invalid 'ligand' of %r." % ligand)

    @classmethod
    def from_length(cls, length, temp, ligand=None, database_dir="database/", interactions=miyazawa_jernigan):
        """Create a thermodynamic object for sequences of a given length.
        """
        if not os.path.exists(database_dir):
            print("Creating a conformations database in %s" % database_dir)
            os.makedirs(database_dir)
        confs = conformations.Conformations(length, database_dir=database_dir, interaction_energies=interactions)
        self = cls(temp, confs, ligand)
        return self

    def native_conf(self, seq):
        """Return the native conformation."""
        (minE, conf, partitionsum, numcontacts, folds) = self._nativeE(seq)
        return conf

    def k_lowest_confs(self, seq, k):
        """Get the `k` lowest conformations in the sequence's conformational ensemble.
        """
        length = len(seq)
        # Calculate the kth lowest conformations
        ncontacts = self._conformations.max_contacts()
        confs = np.array(self._conformations.unique_conformations(ncontacts))
        energies = np.empty(len(confs), dtype=float)
        for i, conf in enumerate(confs):
            energies[i] = fold_energy(seq, conf, interactions=self._conformations._interaction_energies)
        sorted_e = np.argsort(energies)
        states = confs[sorted_e[0:k]]
        return states

    def nativeE(self, seq, target=None):
        """Compute the native energy and return it.

        Parameters
        ----------
        seq : str or list
            sequence to fold.

        Returns
        -------
        minE : float
            Energy of the native state.
        """
        if target is None:
            (minE, conf, partitionsum, numcontacts, folds) = self._nativeE(seq)
        else:
            minE = fold_energy(seq, target, self._conformations._interaction_energies)
        return minE

    def _nativeE(self, seq):
        """Compute the lattice native energy and partition sum of a sequence.

        Parameters
        ----------
        seq : str or list
            sequence to fold.

        Returns
        -------
        minE : float
            Energy of the native state.
        conf : str
            Native conformation
        partitionsum : float
            Partition function sum
        numcontacts : int
            Number of contacts in native state
        folds : boolean
            True if a single native structure exists. False is not.
        """
        if len(seq) != self.length():
            raise ThermodynamicsError("Invalid 'seq' of %r." % seq)
        return self._conformations.fold_sequence(seq, self._temp)
    #---------------------------------------------------------------------
    def stability(self, seq, target=None):
        """Computes the stability of a sequence if it is below cutoff.

        Parameters
        ----------
        seq : str or list
            sequence to fold.

        Returns
        -------
        stability : float
            Folding stability of the native state.
        """
        nativeE_results = self._nativeE(seq)
        if target is not None:
            minE = fold_energy(seq, target, self._conformations._interaction_energies)
            nativeE_results = list(nativeE_results)
            nativeE_results[0] = minE
            nativeE_results[1] = target
            nativeE_results[-1] = True
        stability_results = self._stability(*nativeE_results)
        return stability_results[0]

    def _stability(self, minE, conf, partitionsum, numcontacts, folds):
        """Computes a stability from minE and partition function.

        Parameters
        ----------
        seq : str or list
            sequence to fold.

        Returns
        -------
        stability : float
            Folding stability of the native state.
        conf : str
            Native conformation
        partitionsum : float
            Partition function sum
        numcontacts : int
            Number of contacts in native state
        folds : boolean
            True if a single native structure exists. False is not.
        """
        #(minE, conf, partitionsum, numcontacts, folds) = self._NativeE(seq)
        # Calculate a stability... if calculation does not work, stability = 0
        if folds:
            gu = - self._temp * math.log(partitionsum - math.exp(-minE / self._temp))
            dGf = minE - gu
            return (dGf, conf, partitionsum, numcontacts, folds)
        else:
            return (0, conf, partitionsum, numcontacts, folds)
    #---------------------------------------------------------------------
    def fracfolded(self, seq, target=None):
        """Compute the fraction folded of the sequence.

        Parameters
        ----------
        seq : str or list
            sequence to fold.

        Returns
        -------
        fracfolded : float
            fractioned folded.
        """
        nativeE_results = self._nativeE(seq)
        if target != None:
            minE = fold_energy(seq, target, self._conformations._interaction_energies)
            nativeE_results = list(nativeE_results)
            nativeE_results[0] = minE
            nativeE_results[1] = target
            nativeE_results[-1] = True
        stability_results = self._stability(*nativeE_results)
        fracfolded_results = self._fracfolded(*stability_results)
        return fracfolded_results[0]

    def _fracfolded(self, dG, conf, partitionsum, numcontacts, folds):
        """Computes the fitness from a given stability value

        Parameters
        ----------
        seq : str or list
            sequence to fold.

        Returns
        -------
        fracfolded : float
            Energy of the native state.
        conf : str
            Native conformation
        partitionsum : float
            Partition function sum
        numcontacts : int
            Number of contacts in native state
        folds : boolean
            True if a single native structure exists. False is not.
        """
        # folding to a target conformation
        #(dG, conf, partitionsum, numcontacts, folds) = self._Stability(seq)
        if folds is False:
            return (0, conf, partitionsum, numcontacts, folds)
        else:
            if self._ligand:
                if dG > self._ligand[2]:
                    return (0, conf, partitionsum, numcontacts, False) # does not stably fold
                else:
                    be = conformations.bind_ligand(seq, conf, self._ligand[0], self._ligand[1])[0]
                    return (math.exp(-be), conf, partitionsum, numcontacts, folds)  # compute the fitness
            else:
                f = 1.0 / (1.0 + math.exp(dG / self._temp))
                return (f, conf, partitionsum, numcontacts, folds)

    #---------------------------------------------------------------------
    def all_metrics(self, seq):
        """Compute lattice NativeE, Stability, and Fitness of a given sequence.

        Parameters
        ----------
        seq : str
            protein sequence string.

        Returns
        -------
        nativeE : float
            energy of the native state.
        dG : float
            stability of the native state.
        fitness : float
            fitness of the native state.
        """
        metrics = self._all_metrics(seq)
        return metrics[0], metrics[1], metrics[2]

    def _all_metrics(self, seq):
        """Compute lattice NativeE, Stability, and Fitness of a given sequence.

        Parameters
        ----------
        seq : str
            protein sequence string.

        Returns
        -------
        minE : float
            energy of native conformation
        stability : float
            stability of native energy
        fracfolded : float
            fraction folded.
        conf : str
            Native conformation
        partitionsum : float
            Partition function sum
        numcontacts : int
            Number of contacts in native state
        folds : boolean
            True if a single native structure exists. False is not.
        """
        nativeE_results = self._nativeE(seq)
        stability_results = self._stability(*nativeE_results)
        fracfolded_results = self._fracfolded(*stability_results)
        return (nativeE_results[0], stability_results[0],
                fracfolded_results[0], *nativeE_results[1:])
    #---------------------------------------------------------------------
    def length(self):
        """Returns the sequence length for which fitnesses are computed."""
        return self._conformations.length()
