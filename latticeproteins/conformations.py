#!/usr/bin/python
# Begin conformations.py
#---------------------------------------------------------------------------
"""Module for constructing conformation database for sequences of set length.

Originally written by Jesse Bloom, 2004.

Updated by Zach Sailer, 2017."""
#-----------------------------------------------------------------------
import math, sys, os
import numpy as np
from latticeproteins.interactions import miyazawa_jernigan
# Python 3 compatibility
try:
    import cPickle as pickle
except ImportError:
    import pickle
#----------------------------------------------------------------------
class ConformationsError(Exception):
    """Error finding or storing a conformation."""
#----------------------------------------------------------------------
# '_loop_in_C' is a Boolean switch specifying whether we try to speed
# up the energy calculations by doing the looping with the C-extension
# 'contactlooper'.  'True' means we try to do this.
_loop_in_C = True
if _loop_in_C:
    from latticeproteins.contactlooper import ContactLooper

class PickleProtocolError(Exception):
    """Error is pickle version is too old. """

PROTOCOL = pickle.HIGHEST_PROTOCOL
if PROTOCOL < 2:
    raise PickleProtocolError("Version of pickle is outdated for this package.")

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
    x = y = length # initial position on the grid is at the center of the 2d array
    grid[x,y] = sites[0]
    # move on grid, populate with amino acid at that site, and store all contacting neighbors.
    contacts = []
    for i in range(length-1):
        step = coordinates[moves[i]]
        x += step[0]
        y += step[1]
        grid[x,y] = sites[i+1]
        neighbors = [sites[i+1] + grid[x+c[0], y+c[1]] for c in coordinates.values()]
        contacts += [n for n in neighbors if len(n) == 2]
    # subtract the contacts that have bonds between them.
    for i in range(1,length):
        # Try forward
        try:
            contacts.remove(sequence[i-1:i+1])
        # Else try reverse
        except ValueError:
            contacts.remove(sequence[i] + sequence[i-1])
    return contacts

#----------------------------------------------------------------------
# Classes
#----------------------------------------------------------------------

class Conformations(object):
    """Creates a database of conformations for a protein of specified length.

    The created 'Conformations' object 'c' stores the contact
        lists and the number of conformations with these contact sets
        for all self-avoiding walks of length 'length'.  It can then
        be used to compute the free energy of a protein folding to
        the lowest energy conformation.

    Parameters
    ----------
    length :
        is an integer specifying the length of the protein
        for which we are computing the contacts.  It must be >= 2.
    database_dir :
        specifies the name of the database directory storing
        existing conformations.  If the conformation instance
        already exists in this database we return the
        existing data, and if it doesn't we store it in the database.
    interaction_energies :
        specifies the interaction energies between
        residues. By default, this is interactions.miyazawa_jernigan.


    Attributes
    ----------
    _numconformations : dict
        A dictionary mapping the number of contact sets to the number of conformations
        with that contact set.
    _contactsets : list of lists
        'self._contactsets' is a list of contact sets.  'self._contactsets[i]'
        is the contact set for contact i.  It is a list of numbers.
        'x = self._contactsets[i]' describes the residues in contact
        in contact 'i'.  If this contact is between residues 'ires'
        and 'jres', then 'x = self._length * ires + jres' where
        0 <= ires, jres < 'self._length', and ires < jres + 1 contact sets
    _contactsetdegeneracy : list
        'self._contactsetdegeneracy' is a list of integers giving the
        degeneracy of the contact sets (the number of different conformations
        with this contact set).  'self_contactsetdegeneracy[i]' is the
        degeneracy of the contact set 'self._contactsets[i]'
    _contactsetconformation : list
        'self._contactsetconformation' is a list of the conformations
        associated with each contact set.  If contact set 'self._contactsets[i]'
        is degenerate ('self._contactsetdegeneracy[i]' > 1), the value
        'self._contactsetconformation[i]' is 'None'.  Otherwise, it
        is the string representing the conformation that gives rise
        to contact set 'self._contactsets[i]'.  The conformations are given
        such that 'self._contactsetconformation[i][j]'
        gives the conformation of bond 'j' (0 <= j < 'self._length' - 1)
        as 'U' (Up), 'R' (Right), 'D' (Down), or 'L' (Left).
        We require the first bond to be Up, and the first
        non-Up bond to be Right.
    _numcontactsets : dict
        'self._numcontactsets[i]' holds the number of different contact
        sets with 'i' contacts.
    """
    def __init__(self, length, database_dir="database/", interaction_energies=miyazawa_jernigan):
        self._interaction_energies = interaction_energies
        if not os.path.isdir(database_dir):
            os.makedirs(database_dir)
            #raise IOError("Cannot find database directory of %s." % database_dir)
        object_properties = ['_length', '_numconformations', '_contactsets', '_contactsetdegeneracy', '_contactsetconformation', '_numcontactsets']
        foundone = didntfindone = False
        for prop in object_properties:
            f = "%s/%d%s.pickle" % (database_dir, length, prop)
            if os.path.isfile(f):
                try:
                    val = pickle.load(open("%s/%d_length.pickle" % (database_dir, length), 'rb'))
                    foundone = True
                except ValueError:
                    foundone = False
            else:
                didntfindone = True
        if foundone and not didntfindone:
            # return existing values
            self._length = pickle.load(open("%s/%d_length.pickle" % (database_dir, length), 'rb'))
            if self._length != length:
                raise ValueError("Length mismatch.")
            self._numconformations = pickle.load(open("%s/%d_numconformations.pickle" % (database_dir, length), 'rb'))
            self._contactsets = pickle.load(open("%s/%d_contactsets.pickle" % (database_dir, length), 'rb'))
            self._contactsetdegeneracy = pickle.load(open("%s/%d_contactsetdegeneracy.pickle" % (database_dir, length), 'rb'))
            self._contactsetconformation = pickle.load(open("%s/%d_contactsetconformation.pickle" % (database_dir, length), 'rb'))
            self._numcontactsets = pickle.load(open("%s/%d_numcontactsets.pickle" % (database_dir, length), 'rb'))
            # sort the contact set information so that the conformations
            # with the most contacts come first
            n = len(self._contactsets)
            decorated_list = [(len(self._contactsets[i]), self._contactsets[i], self._contactsetdegeneracy[i], self._contactsetconformation[i]) for i in range(n)]
            decorated_list.sort()
            decorated_list.reverse()
            self._contactsets = [decorated_list[i][1] for i in range(n)]
            self._contactsetdegeneracy = [decorated_list[i][2] for i in range(n)]
            self._contactsetconformation = [decorated_list[i][3] for i in range(n)]
            self._foldedsequences = {}
            return
        elif foundone and didntfindone:
            raise ValueError("Found some but not all conformations for length %d in %s." % (length, database_dir))
        # If we have made it here, we are generating new information for this conformation
        # A listing of all class variables
        # 'self._length' is the length of the protein (self-avoiding walk)
        self._length = length
        # there are 'self._numconformations[i]' conformations with 'i' contacts
        self._numconformations = {}
        # 'self._contactsets' is a list of contact sets.  'self._contactsets[i]'
        # is the contact set for contact i.  It is a list of numbers.
        # 'x = self._contactsets[i]' describes the residues in contact
        # in contact 'i'.  If this contact is between residues 'ires'
        # and 'jres', then 'x = self._length * ires + jres' where
        # 0 <= ires, jres < 'self._length', and ires < jres + 1
        self._contactsets = []
        # 'self._contactsetdegeneracy' is a list of integers giving the
        # degeneracy of the contact sets (the number of different conformations
        # with this contact set).  'self_contactsetdegeneracy[i]' is the
        # degeneracy of the contact set 'self._contactsets[i]'
        self._contactsetdegeneracy = []
        # 'self._contactsetconformation' is a list of the conformations
        # associated with each contact set.  If contact set 'self._contactsets[i]'
        # is degenerate ('self._contactsetdegeneracy[i]' > 1), the value
        # 'self._contactsetconformation[i]' is 'None'.  Otherwise, it
        # is the string representing the conformation that gives rise
        # to contact set 'self._contactsets[i]'.  The conformations are given
        # such that 'self._contactsetconformation[i][j]'
        # gives the conformation of bond 'j' (0 <= j < 'self._length' - 1)
        # as 'U' (Up), 'R' (Right), 'D' (Down), or 'L' (Left).
        # We require the first bond to be Up, and the first
        # non-Up bond to be Right.
        self._contactsetconformation = []
        # 'self._numcontactsets[i]' holds the number of different contact
        # sets with 'i' contacts.
        self._numcontactsets = {}
        # 'self._foldedsequences = {}' holds recently folded sequences.
        # The keys are the arguments to 'FoldSequence', and the items
        # are how many times the saved sequence was accessed and the
        # information about the folded sequence.  Sequences in the keys
        # are strings.
        self._foldedsequences = {}
        #
        # Now do some error checking on input variables
        if not (isinstance(self._length, int) and self._length >= 2):
            raise ConformationsError("Invalid 'length' of %r." % self._length)
        #
        # The initial generation of conformations uses the variable
        # 'contactsets' which is a dictionary as described below.  Once
        # we have constructed this dictionary, we use it to create
        # 'self._contactsets', 'self._contactsetdegeneracy', and
        # 'self._contactsetconformation'.  The format of 'contactsets':
        # The keys are string representing the contact sets.  These
        # strings are the numbers that will appear in 'self._contactsets',
        # separated by spaces.
        # The values of 'contactsets' provides information
        # about the conformations encoding the contact sets.  If there is
        # a single conformation encoding a contact set 'cs', then
        # 'contactsets[cs]' is a string.  If the contact set stored as
        # 'contactsets[cs]' is coded for by multiple conformations,
        # then 'self._contactsets[cs]' is an integer > 1 that represents
        # the number of conformations coding for this contact set.
        contactsets = {}
        # Now being looping over all possible conformations
        # The initial conformation is all 'U'
        dx = {'U':0, 'R':1, 'D':0, 'L':-1}
        dy = {'U':1, 'R':0, 'D':-1, 'L':0}
        next = {'U':'R', 'R':'D', 'D':'L', 'L':'U'}
        opposite = {'U':'D', 'D':'U', 'R':'L', 'L':'R'}
        n = self._length - 2 # index of last bond in 'conformation'
        conformation = ['U' for i in range(n + 1)]
        first_R = n # index of the first 'R' in the conformation
        ncount = 0
        while True: # keep finding new conformations
            # See if the current conformation has overlap
            # If no overlap, store the contact set
            x = y = j = 0
            res_positions = {(x, y):j} # keyed by coords, items are residue numbers
            res_coords = [(x, y)] # 'res_coords[j]' is coords of residue 'j'
            for c in conformation:
                x += dx[c]
                y += dy[c]
                coord = (x, y)
                if coord in res_positions: # overlap
                    # increment at the step that gave the problem
                    for k in range(j + 1, n + 1):
                        conformation[k] = 'U'
                    conformation[j] = next[conformation[j]]
                    while conformation[j] == 'U':
                        j -= 1
                        conformation[j] = next[conformation[j]]
                    if j == first_R and conformation[j] not in ['R', 'U']:
                        first_R -= 1
                        conformation[first_R] = 'R'
                        for k in range(j, n + 1):
                            conformation[k] = 'U'
                    break
                j += 1
                res_positions[coord] = j
                res_coords.append(coord)
            else: # loop finishes normally, this is a valid conformation
                # generate the contact set
                cs = []
                numcontacts = 0
                for j in range(len(res_coords)):
                    (x, y) = res_coords[j]
                    partners_list = []
                    for c in ['U', 'R', 'D', 'L']:
                        try:
                            k = res_positions[(x + dx[c], y + dy[c])]
                            if k > j + 1:
                                partners_list.append(k)
                        except KeyError:
                            pass
                    partners_list.sort()
                    jtimeslength = j * self._length
                    for k in partners_list:
                        cs.append("%d " % (jtimeslength + k))
                        numcontacts += 1
                cs = ''.join(cs)
                # add conformation to count
                try:
                    self._numconformations[numcontacts] += 1
                except KeyError:
                    self._numconformations[numcontacts] = 1
                # store contact set
                try:
                    contactsets[cs] += 1
                except KeyError:
                    contactsets[cs] = ''.join(conformation)
                    try:
                        self._numcontactsets[numcontacts] += 1
                    except KeyError:
                        self._numcontactsets[numcontacts] = 1
                except TypeError:
                    contactsets[cs] = 2
                # generate the next conformation
                i = n
                conformation[i] = next[conformation[i]]
                while conformation[i] == 'U':
                    i -= 1
                    conformation[i] = next[conformation[i]]
                # make sure first non-'U' is 'R'
                if i == first_R and conformation[i] not in ['R', 'U']:
                    first_R -= 1
                    conformation[first_R] = 'R'
                    for j in range(i, n + 1):
                        conformation[j] = 'U'
            # see if we are done
            if first_R == 0:
                break
        #
        # Now use 'contactsets' to generate 'self._contactsets',
        # 'self._contactsetdegeneracy', and 'self._contactsetconformation'
        for (cs, n_or_conf) in contactsets.items():
            # convert 'cs' to the appropriate list
            clist = [int(x) for x in cs.split()]
            self._contactsets.append(clist)
            if isinstance(n_or_conf, str):
                self._contactsetdegeneracy.append(1)
                self._contactsetconformation.append(n_or_conf)
            else:
                self._contactsetdegeneracy.append(n_or_conf)
                self._contactsetconformation.append(None)
        #
        del contactsets
        #
        # to make the energy evaluations quicker, sort so that
        # contact sets with more conformations are first:
        decorated_list = [(len(self._contactsets[i]), self._contactsets[i], self._contactsetdegeneracy[i], self._contactsetconformation[i]) for i in range(len(self._contactsets))]
        del self._contactsets, self._contactsetdegeneracy, self._contactsetconformation
        decorated_list.sort()
        decorated_list.reverse()
        self._contactsets = [decorated_list[i][1] for i in range(len(decorated_list))]
        self._contactsetdegeneracy = [decorated_list[i][2] for i in range(len(decorated_list))]
        self._contactsetconformation = [decorated_list[i][3] for i in range(len(decorated_list))]
        #
        # store the conformations data in the database
        pickle.dump(self._length, open("%s/%d_length.pickle" % (database_dir, length), 'wb'), protocol=PROTOCOL)
        pickle.dump(self._numconformations, open("%s/%d_numconformations.pickle" % (database_dir, length), 'wb'), protocol=PROTOCOL)
        pickle.dump(self._contactsets, open("%s/%d_contactsets.pickle" % (database_dir, length), 'wb'), protocol=PROTOCOL)
        pickle.dump(self._contactsetdegeneracy, open("%s/%d_contactsetdegeneracy.pickle" % (database_dir, length), 'wb'), protocol=PROTOCOL)
        pickle.dump(self._contactsetconformation, open("%s/%d_contactsetconformation.pickle" % (database_dir, length), 'wb'), protocol=PROTOCOL)
        pickle.dump(self._numcontactsets, open("%s/%d_numcontactsets.pickle" % (database_dir, length), 'wb'), protocol=PROTOCOL)

    def length(self):
        """Returns the length of the protein these conformations are for."""
        return self._length

    def k_lowest_confs(self, seq, temp, k):
        """Get the `k` lowest conformations in the sequence's conformational ensemble.
        """
        length = len(seq)
        # Calculate the kth lowest conformations
        ncontacts = self.max_contacts()
        confs = []
        nconfs, n = 0, ncontacts
        while nconfs < k and n > 1:
            confs += self.unique_conformations(n)
            nconfs = len(confs)
            n -= 1
        # Loop through conformations and find lowest energy.
        confs = np.array(confs)
        energies = np.empty(len(confs), dtype=float)
        for i, conf in enumerate(confs):
            energies[i] = fold_energy(seq, conf, interactions=self._interaction_energies)
        sorted_e = np.argsort(energies)
        states = confs[sorted_e[0:k]]
        return states

    def fold_sequence(self, seq, temp):
        """Folds a protein sequence; calculates native energy and partition sum.

        Parameters
        ----------
        seq : string
            is the sequence of the protein to be folded as one-letter amino
            acid codes.  It should be a string or list of length 'c.Length()'.
        temp :
            is the temperature at which the protein is to be folded.  It
            must be a number > 0.  It represents a reduced temperature, scaled
            so that a value of 1 represents 273 K.

        Returns
        -------
        minE : float
            The energy of the lowest energy conformation.
        conf : str
            Lowest energy conformation.
        partitionsum : float
            Total partition function sum.
        numcontacts : int
            Number of contacts in the native conformation.
        folds : bool
            True if lattice protein has a single lowest energy.
        """
        # do some error checking on the input variables
        if len(seq) != self._length:
            raise ConformationsError("Invalid 'seq' length: is %r but should be %r." % (len(seq), self._length))
        try:
            temp = float(temp)
            if temp <= 0.0:
                raise ConformationsError("'temp' is <= 0: %r." % temp)
        except KeyError:
            raise ConformationsError("Invalid 'temp' of %r." % temp)

        # create 'res_interactions'.  'res_interactions[j]' holds the energy
        # for the interaction 'j' as specified in 'self._contactsets[i][j]'
        res_interactions = [0.0 for i in range(self._length**2)]
        for ires in range(self._length):
            itimeslength = ires * self._length
            for jres in range(ires + 1, self._length):
                respair = "%s%s" % (seq[ires], seq[jres])
                res_interactions[itimeslength + jres] = self._interaction_energies[respair]

        # Use C function to loop through contacts
        (minE, ibest, partitionsum, folds) = ContactLooper(res_interactions, self._contactsets, self._contactsetdegeneracy, float(temp))

        # Set folds to boolean
        if folds is 1:
            folds = True
        else:
            folds = False

        # Get the native conformation
        conf = self._contactsetconformation[ibest]
        return minE, conf, partitionsum, folds

    def unique_conformations(self, numcontacts):
        """Gets all unique conformations with specified number of contacts.

        Parameters
        ----------
        numcontacts : int
            Number of contacts to include in unique conformations list.

        Returns
        -------
        clist : list
            is of all unique conformations
            with exactly 'numcontacts' contacts.  A conformation
            is "unique" if it is the only conformation that gives
            rise to its particular contact set.  If there are
            no unique conformations with 'numcontacts' contacts,
            'clist' is an empty list.  Conformations are specified
            as strings of 'U', 'R', 'L', and 'D' as described in
            'FoldSequence'."""
        if not (isinstance(numcontacts, int) and numcontacts >= 0):
            raise ConformationsError("Invalid 'numcontacts' of %r." % numcontacts)
        clist = []
        for i in range(len(self._contactsets)):
            if self._contactsetdegeneracy[i] == 1:
                if len(self._contactsets[i]) == numcontacts:
                    clist.append(self._contactsetconformation[i])
        return clist

    def max_contacts(self):
        """Gets the most contacts of any conformation.

        Returns
        -------
        n : int
            is returned as the number of contacts for the conformation
            with the most contacts."""
        n = 0
        for cs in self._contactsets:
            if len(cs) > n:
                n = len(cs)
        return n

    def num_conformations(self, contacts = None):
        """Returns the number of conformations.

        If 'contacts' has its default value of 'None', returns the total
            number of conformations (self-avoiding walks).
        If 'contacts' has an integer value, returns the number of conformations
            with 'contacts' contacts.  If there are no walks with this number
            of contacts, returns 0.
        """
        if contacts == None:
            n = 0
            for x in self._numconformations.values():
                n += x
        else:
            try:
                n = self._numconformations[contacts]
            except KeyError:
                if isinstance(contacts, int) and contacts >= 0:
                    n = 0
                else:
                    raise ConformationsError("Invalid 'contacts' of %r." % contacts)
        return n

    def num_contact_sets(self, contacts = None):
        """Returns the number of unique contact sets.

        If 'contacts' has its default value of 'None', returns the total
            number of unique contact sets (defined as the list of all
            contacts of non-adjacent residues).
        If 'contacts' has an integer value, returns the number of unique
            contact sets with 'contacts' contacts.  If there are no
            contact sets with this number of contacts, returns 0."""
        if contacts == None:
            return len(self._contactsets)
        else:
            try:
                return self._numcontactsets[contacts]
            except KeyError:
                if isinstance(contacts, int) and contacts >= 0:
                    n = 0
                else:
                    raise Conformationserror("Invalid 'contacts' of %r." % contacts)

class ConformationList(object):
    """Build an Conformations like object without the database. Uses a list of
    conformations provided by user to construct a Conformations object.

    **Note** This will likely be much slower at calculating large lists of
    conformations.

    Parameters
    ----------
    length :
        is an integer specifying the length of the protein
        for which we are computing the contacts.  It must be >= 2.
    conflist :
        a list of conformations.
    interaction_energies :
        specifies the interaction energies between
        residues. By default, this is interactions.miyazawa_jernigan.
    """
    def __init__(self, conflist, interaction_energies=miyazawa_jernigan):
        self._interaction_energies = miyazawa_jernigan
        self._conflist = conflist
        self._length = len(self._conflist[0]) + 1

    def length(self):
        """Returns the length of the protein these conformations are for."""
        return self._length

    def fold_sequence(self, seq, temp):
        """Folds a protein sequence, calculate native energy and partition sum.

        Parameters
        ----------
        seq : string
            is the sequence of the protein to be folded as one-letter amino
            acid codes.  It should be a string or list of length 'c.Length()'.
        temp :
            is the temperature at which the protein is to be folded.  It
            must be a number > 0.  It represents a reduced temperature, scaled
            so that a value of 1 represents 273 K.

        Returns
        -------
        minE : float
            The energy of the lowest energy conformation.
        conf : str
            Lowest energy conformation.
        partitionsum : float
            Total partition function sum.
        numcontacts : int
            Number of contacts in the native conformation.
        folds : bool
            True if lattice protein has a single lowest energy.
        """
        Es = []
        minE, partitionsum = 0, 0
        conf = ""
        folds = False
        for c in self._conflist:
            E = fold_energy(seq, c)
            Es.append(E)
            partitionsum += math.exp(-E/temp)
            if E < minE:
                minE = E
                conf = c
                folds = True
            # If the energy equals the current minimum, sequence doesn't fold
            elif E == minE:
                folds = False
        return minE, conf, partitionsum, folds

    def unique_conformations(self, numcontacts):
        """Gets all unique conformations with specified number of contacts.

        Parameters
        ----------
        numcontacts : int
            Number of contacts to include in unique conformations list.

        Returns
        -------
        clist : list
            is of all unique conformations
            with exactly 'numcontacts' contacts.  A conformation
            is "unique" if it is the only conformation that gives
            rise to its particular contact set.  If there are
            no unique conformations with 'numcontacts' contacts,
            'clist' is an empty list.  Conformations are specified
            as strings of 'U', 'R', 'L', and 'D' as described in
            'FoldSequence'."""
        clist = list(np.unique(self._conflist))
        return clist

    def num_conformations(self, contacts = None):
        """Returns the number of conformations.

        If 'contacts' has its default value of 'None', returns the total
            number of conformations (self-avoiding walks).
        If 'contacts' has an integer value, returns the number of conformations
            with 'contacts' contacts.  If there are no walks with this number
            of contacts, returns 0.
        """
        return len(self._conflist)

    def max_contacts(self):
        """Gets the most contacts of any conformation.

        Returns
        -------
        n : int
            is returned as the number of contacts for the conformation
            with the most contacts.
        """
        n = 0
        dummy_seq = "0"*(self._length)
        for c in self._conflist:
            cs = lattice_contacts(dummy_seq, c)
            if len(cs) > n:
                n = len(cs)
        return n

    def num_contact_sets(self, contacts = None):
        """Returns the number of unique contact sets.

        If 'contacts' has its default value of 'None', returns the total
            number of unique contact sets (defined as the list of all
            contacts of non-adjacent residues).
        If 'contacts' has an integer value, returns the number of unique
            contact sets with 'contacts' contacts.  If there are no
            contact sets with this number of contacts, returns 0.
        """
        contacts = []
        dummy_seq = "0"*self.length + 1
        for c in self._conflist:
            contacts.append(lattice_contacts(dummy_seq, c))
        return len(contacts)
