import numpy as np
import pandas as pd

from .conformations import Conformations
from .thermodynamics import LatticeThermodynamics, LatticeGroupThermodynamics

class LatticeProtein(object):
    """A single lattice protein.

    Attributes
    ----------
    native_energy : float
        native energy

    native_conf : str
        native conformation

    partition_sum : float
        partition function.

    folded : bool
        is it folded or not?

    stability : float
        folding stability of nativec conf (or target)

    fracfolded : float
        fraction of protein folded (the probability of the native state.)
    """
    def __init__(self, sequence, target=None, conformations=None,
                 temp=1.0, **kwargs):
        self.sequence = sequence
        self.length = len(sequence)
        self.target = target
        self.temp = temp

        self.conformations = conformations
        if conformations is None:
            self.conformations = Conformations(self.length)

        _ = self.conformations.fold_sequence(self.sequence, self.temp)

        self.native_energy = _[0]
        self.partition_sum = _[2]
        self.folded = _[3]

        if target is  None:
            self.native_conf = _[1]
        else:
            self.native_conf = target

        # Calculate stability
        self.stability = self.native_energy + (
            self.temp * np.log(self.partition_sum -
            np.exp(-self.native_energy / self.temp))
        )

        # Calculate fraction folded,
        self.fracfolded = 1.0 / (1.0 + np.exp(self.stability / self.temp))

    def k_lowest_confs(self, k):
        return self.conformations.k_lowest_confs(
            self.sequence, self.temp, k)


class LatticeProteins(object):
    """A group of lattice proteins.

    Attributes
    ----------
    native_energy : float
        native energy

    native_conf : str
        native conformation

    partition_sum : float
        partition function.

    folded : bool
        is it folded or not?

    stability : float
        folding stability of nativec conf (or target)

    fracfolded : float
        fraction of protein folded (the probability of the native state.)
    """
    def __init__(self, sequence_list, target=None, conformations=None,
                 temp=1.0, **kwargs):
        self.sequence_list = sequence_list
        self.length = len(sequence_list[0])
        self.n = len(self.sequence_list)
        self.target = target
        self.temp = temp

        self.conformations = conformations
        if conformations is None:
            self.conformations = Conformations(self.length)

        self.native_energy = np.empty(self.n, dtype=float)
        self.native_conf = np.empty(self.n, dtype=str)
        self.partition_sum = np.empty(self.n, dtype=float)
        self.folded = np.empty(self.n, dtype=bool)

        for i, seq in enumerate(self.sequence_list):
            _ = self.conformations.fold_sequence(seq, self.temp)
            self.native_energy[i] = _[0]
            self.native_conf[i] = _[1]
            self.partition_sum[i] = _[2]
            self.folded[i] = _[3]

        # Calculate stability
        self.stability = self.native_energy + (
            self.temp * np.log(self.partition_sum -
            np.exp(-self.native_energy / self.temp))
        )

        # Calculate fraction folded,
        self.fracfolded = 1.0 / (1.0 + np.exp(self.stability / self.temp))
