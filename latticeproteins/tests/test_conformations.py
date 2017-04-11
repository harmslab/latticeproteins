
import os, shutil
import unittest
from nose import with_setup
from ..conformations import Conformations

class TestConformations(unittest.TestCase):

    def setUp(self):
        self.database = 'database'
        self.length = 6
        self.seq = "AHKEDP"
        self.temp = 1

    def tearDown(self):
        # Remove a database if successfully created.
        try:
            shutil.rmtree("database")
        except FileNotFoundError:
            pass

    def test__init__(self):
        c = Conformations(self.length, self.database)
        self.assertEqual(os.path.exists(self.database), True)
        self.assertNotEqual(os.path.exists(self.database), False)
        self.assertEqual(c._length, self.length)

    def test_fold_sequence(self):
        c = Conformations(self.length, self.database)
        minE, conf, partitionsum, folds = c.fold_sequence(self.seq, self.temp)
        self.assertEqual(type(folds), bool)
        self.assertEqual(len(conf), self.length-1)

    def test_k_lowest_confs(self):
        c = Conformations(self.length, self.database)
        out = c.k_lowest_confs(self.seq, self.temp, 4)
        # Number of lowest conformations better be 4
        self.assertEqual(len(out), 4)
        # Length of conformation string is length -1
        self.assertEqual(len(out[0]), self.length - 1)

    def test_max_contacts(self):
        c = Conformations(self.length, self.database)
        n = c.max_contacts()
        # Should only see 2 contacts in 6-site lattice protein
        self.assertEqual(n, 2)

    def test_num_conformations(self):
        c = Conformations(self.length, self.database)
        n = c.num_conformations()
        # Should have 36 conformations for a 6 site lattice protein
        self.assertEqual(n, 36)
