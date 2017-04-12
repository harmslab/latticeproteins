
import os, shutil
import unittest
from nose import with_setup
from nose import tools


from ..conformations import lattice_contacts, Conformations, ConformationList
from ..thermodynamics import LatticeThermodynamics, ThermodynamicsError


class TestLatticeThermodynamics(unittest.TestCase):

    def setUp(self):
        self.database = 'database'
        self.length = 6
        self.seq = "AHKEDP"
        self.temp = 1
        # Initialize a Conformations object
        self.confs = Conformations(self.length, database_dir=self.database)
        # Initialize a ConformationList object with 3 states
        self.states = self.confs.k_lowest_confs(self.seq, self.temp, 3)
        self.conflist = ConformationList(self.length, self.states)

    def tearDown(self):
        # Remove a database if successfully created.
        try:
            shutil.rmtree("database")
        except FileNotFoundError:
            pass

    def test__init__(self):
        # Test initializeation with regular Conformations object
        lattice1 = LatticeThermodynamics(self.temp, self.confs)
        self.assertEquals(lattice1._temp, self.temp)
        self.assertIsInstance(lattice1._conformations, Conformations)

        # Test initialization with ConformationList
        lattice2 = LatticeThermodynamics(self.temp, self.conflist)
        self.assertEquals(lattice2._temp, self.temp)
        self.assertIsInstance(lattice2._conformations, ConformationList)

        # Test that error is throw with bad arguments to initialization.
        self.assertRaises(TypeError, LatticeThermodynamics, self.temp, 4)
        self.assertRaises(ThermodynamicsError, LatticeThermodynamics, self.conflist, self.conflist)
