import numpy as np
from .interactions import miyazawa_jernigan

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
    """ Find all contacts in conformation.

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
