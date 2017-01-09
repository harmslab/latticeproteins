Tutorials
=========

Getting started
---------------
There are a few things you need to know before you get started. Due to the
poor scaling of lattice protein calculations, the ``latticeproteins`` package takes
a few precautions. First, folding a sequence is done in separate steps (and functions) rather than
through a single call. This forces you to be aware of the magnitude of each call.
You'll notice in these tutorials, we have to import many functions.

Second, the hardest step (most memory and time) in the calculation, by far, is enumerating conformations on
the grid. ``latticeproteins`` tries to reduce the pain by creating a pickled
database of conformations after the first creation. If you delete this directory,
it will have to recreate it next time you run calculations.

Folding a random sequence
-------------------------



.. code-block:: python

    import
