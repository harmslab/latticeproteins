#!/usr/bin/python
# Begin fitness.py
#---------------------------------------------------------------------
"""Module for calculating thermodynamics of lattice protein sequences.

Originally written by Jesse Bloom, 2004.

Updated by Zach Sailer, 2017."""
#----------------------------------------------------------------------
import os
import math, sys
import numpy as np

from . import conformations
from .interactions import miyazawa_jernigan
#----------------------------------------------------------------------
class ThermodynamicsError(Exception):
    """Error computing lattice protein thermodynamics."""

#----------------------------------------------------------------------
class LatticeThermodynamics(object):
    """Attaches thermodynamic evaluators to a lattice protein conformation database.

    Parameters
    ----------
    temp : float
        the temperature at which the fitness is computed.
    confs : conformations.Conformations object
        is the 'conformations.Conformations' object
        used to fold the protein sequences.  'conformations.Length()'
        specifies the length of the protein sequences that can be
        folded.
    """
    #------------------------------------------------------------------
    def __init__(self, temp, confs):
        # Assign class instance variables and error check
        if not (isinstance(temp, (int, float)) and temp > 0):
            raise ThermodynamicsError("Invalid 'temp' of %r." % temp)
        # Check Conformations.
        if not isinstance(confs, conformations.Conformations) and not isinstance(confs, conformations.ConformationList):
            raise TypeError("conformations must be a Conformations or ConformationList object.")
        self._temp = temp
        self._conformations = confs

    @classmethod
    def from_length(cls, length, temp, database_dir="database/", interactions=miyazawa_jernigan):
        """Create a thermodynamic object for sequences of a given length.
        """
        if not os.path.exists(database_dir):
            print("Creating a conformations database in %s" % database_dir)
            os.makedirs(database_dir)
        confs = conformations.Conformations(length, database_dir=database_dir, interaction_energies=interactions)
        self = cls(temp, confs)
        return self

    def native_conf(self, seq):
        """Return the native conformation."""
        (minE, conf, partitionsum, folds) = self._nativeE(seq)
        return conf

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
            (minE, conf, partitionsum, folds) = self._nativeE(seq)
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

    def _stability(self, minE, conf, partitionsum, folds):
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
            return (dGf, conf, partitionsum, folds)
        else:
            return (0, conf, partitionsum, folds)
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

    def _fracfolded(self, dG, conf, partitionsum, folds):
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
            return (0, conf, partitionsum, folds)
        else:
            if self._ligand:
                if dG > self._ligand[2]:
                    return (0, conf, partitionsum, False) # does not stably fold
                else:
                    be = conformations.bind_ligand(seq, conf, self._ligand[0], self._ligand[1])[0]
                    return (math.exp(-be), conf, partitionsum, folds)  # compute the fitness
            else:
                f = 1.0 / (1.0 + math.exp(dG / self._temp))
                return (f, conf, partitionsum, folds)

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
        out = (nativeE_results[0], stability_results[0],
                fracfolded_results[0]) + tuple(nativeE_results[1:])
        return out
    #---------------------------------------------------------------------
    def length(self):
        """Returns the sequence length for which fitnesses are computed."""
        return self._conformations.length()


class GroupThermodynamics(object):
    """Efficiently calculates thermodynamic properties for a list of lattice proteins.

    Parameters
    ----------
    seqlist : list
        List of lattice proteins.
    temp : float
        temperature of the system.
    confs : Conformations or ConformationList object
        Conformation database for lattice with set length
    target : str (optional, default=None)
        target conformation to fold protein list

    Attributes
    ----------
    seqlist : list
        list of sequences.
    temp : float
        temperature of the system.
    nativeEs : array
        native (or target) energy for sequences in seqlist
    stabilities : array
        array of stabilities for sequences in seqlist
    fracfolded : array
        array of fraction folded for sequences in seqlist
    """
    def __init__(self, seqlist, temp, confs, target=None):
        # Assign class instance variables and error check
        if not (isinstance(temp, (int, float)) and temp > 0):
            raise ThermodynamicsError("Invalid 'temp' of %r." % temp)
        self.temp = temp

        # Set conformations after checking
        if not isinstance(confs, conformations.Conformations) and not isinstance(confs, conformations.ConformationList):
            raise TypeError("conformations must be a Conformations or ConformationList object.")
        self._conformations = confs
        self.length = self._conformations.length

        # Check length of target
        if len(target)-1 != self.length:
            raise ThermodynamicsError("Length of target must be length-1")
        self._target = target

        self.seqlist = seqlist
        self.nativeEs = np.empty(len(self.seqlist), dtype=float)
        self.confs = np.empty(len(self.seqlist), dtype="U|" + self.length)
        self._partitionsum = np.empty(len(self.seqlist), dtype=float)

        for i, seq in enumerate(self.seqlist):
            # Check that sequence is valid
            if len(seq) != self.length:
                raise ThermodynamicsError("Length of sequence, %s, does not equal conformation lengths" % seq)

            # Fold sequence and set energy
            out = self._conformations.fold_sequence(seq, self.temp)
            if target is not None:
                self.nativeEs[i] = conformations.fold_energy(seq, target)
            else:
                self.nativeEs[i] = out[0]
            self.confs[i] = out[1]
            self._partitionsum[i] = out[2]

    @property
    def stabilities(self):
        """Folding stability for all sequences in seqlist."""
        gu = - self._temp * np.log(self._partitionsum - np.exp(-self.nativeEs / self.temp))
        dGf = minE - gu
        return dGf

    @property
    def fracfolded(self):
        """Fracfolded folded for all sequences in seqlist."""
        return 1.0 / (1.0 + np.exp(self.stability / self._temp))
