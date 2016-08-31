=================================
Lattice protein simulator
=================================

This is a fork of the lattice protein simulator written by `Jesse Bloom`_.

Please publish Jesse's papers if you use this software:

    * `Protein stability promotes evolvability`_

    * `Stability and the evolvability of function in a model protein`_

This maintains the `GNU Public License`_ of the original package.

The program uses a few C extensions, and so compilation requires the ``gcc`` compiler. To install the package from source::

    python setup.py build
    sudo python setup.py install

This package calculates exact thermodynamic stabilities within the model by summing over the entire ensemble of conformations to compute the partition function. When run, it first creates a database that stores all of these conformations, which can take a substantial amount of time. It will run quickly for short proteins <20 residues, but get increasingly slow and memory-intensive after that.


.. _`Jesse Bloom`: http://research.fhcrc.org/bloom/en.html
.. _`Protein stability promotes evolvability`: http://www.ncbi.nlm.nih.gov/pubmed/16581913
.. _`Stability and the evolvability of function in a model protein`: http://www.ncbi.nlm.nih.gov/pubmed/15111394
.. _`GNU Public License`: http://www.gnu.org/licenses/gpl.html
.. _`Zachary Sailer`: https://github.com/Zsailer
.. _`Harms lab`: http://harmslab.uoregon.edu/
.. _`Version 0.1`: https://github.com/jbloom/latticeproteins/tree/v0.1
.. _`Version 0.2`: https://github.com/jbloom/latticeproteins/tree/v0.2
