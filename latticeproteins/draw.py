__doc__ = """
Module for creating SVG's of protein lattice configurations.

Originally written by Jesse Bloom, 2004.

Updated by Zach Sailer, 2017.

Example call:

    >>> # Create an instance
    >>> drawing = latticegpm.svg.Configuration(sequence, configuration, filename="drawing1.svg")
    >>> # Save to file
    >>> # drawing.save()
    >>> # Print in Jupyter (IPython) notebook
    >>> drawing.notebook

"""

import svgwrite
import numpy as np
from functools import wraps
try:
    from IPython.display import SVG
except ImportError:
    raise ImportError("This module requires a IPython's display module to create "
                    "an SVG drawing. Please install IPython before using this module.")

ROTATE = {"U":"R", "R":"D", "D":"L", "L":"U"}

COLORS = {
    "r": "red",
    "b": "blue",
    "k": "black",
    "g": "green",
    "m": "magenta",
    "y": "yellow",
    ".": "black"
}

def in_notebook(sequence, conf, **kwargs):
    """Creates a Python SVG configuration object. Automagically displays in notebook."""
    drawing = Configuration(sequence, conf, **kwargs)
    return drawing

def to_file(sequence, conf, filename, **kwargs):
    """"""
    drawing = Configuration(sequence, conf, **kwargs)
    drawing.save(filename)

def configuration_to_array(sequence, configuration):
    """Create a square numpy array with the configuration laid out.
    """
    moves = {"U":[0,-1], "D":[0,1], "R":[1,0], "L":[-1,0]}
    # find boundaries for drawing
    xmoves, ymoves = [0], [0]
    # Figure out the perimeter of the configuration
    for i in range(len(configuration)):
        move = configuration[i]
        xmoves.append(xmoves[i] + moves[move][0])
        ymoves.append(ymoves[i] + moves[move][1])
    # bounds of configuration
    xmax = max(xmoves)
    xmin = min(xmoves)
    ymax = max(ymoves)
    ymin = min(ymoves)
    xorigin = 0
    yorigin = 0
    # add letters to grid
    xmax2 = xmax + xmax + 3
    xmin2 = xmin + xmin - 2
    ymax2 = ymax + ymax + 3
    ymin2 = ymin + ymin - 2
    # Size of the 2d grid
    xgrid = abs(xmax2 - xmin2)
    ygrid = abs(ymax2 - ymin2)
    # Build an empty grid object
    grid = [[" " for i in range(xgrid)] for j in range(ygrid)]
    # Set dots on grid
    for i in range(len(grid)):
        for j in range(len(grid[i])):
            if i%2==0 and j%2 == 0:
                grid[i][j] = "."
    # Realign origin to fit in grid
    grid_xorg = xorigin - xmin2
    grid_yorg = yorigin - ymin2
    # initial position
    xpos = grid_xorg
    ypos = grid_yorg
    grid[ypos][xpos] = sequence[0]
    for i in range(len(configuration)):
        let = sequence[i+1]
        bond = configuration[i].lower()
        # Get direction
        step = moves[configuration[i]]
        # Define change in x and y for both edges and letter
        dx = step[0]
        dy = step[1]
        # Add bond
        xpos += dx
        ypos += dy
        grid[ypos][xpos] = bond
        # Add residue
        xpos += dx
        ypos += dy
        grid[ypos][xpos] = let
    return grid


class Configuration(SVG):
    """ Main class for drawing an SVG of a lattice protein's fold.

    Parameters
    ----------
    sequence : str
        Amino acid sequence
    configuration : str
        sequence of direction letters describing the 2d configuration.
    colors : list of strings
        list of colors for each amino acid in sequence
    rotation : int
        rotate the configuration by 0, 90, 180, or 270 degrees
    font_size : int
        Font size, in pixels, of sequence in configuration. The svg will scale
        with the font size of the letters.

    Examples
    --------
    >>> # Create an instance
    >>> drawing = Configuration(sequence, configuration)
    >>> # Save to file
    >>> # drawing.save()
    >>> # Print in Jupyter (IPython) notebook
    >>> drawing.notebook
    """
    def __init__(self, sequence, configuration,
        color_sequence=None,
        rotation=0,
        font_size=20,
        dot_scale=1.0,
        font_weight="normal"
        ):
        # Rotate configuration if given
        # Set sequence and configuration
        self.sequence = sequence
        self.rotation = 0
        self.font_size = font_size
        self.dot_scale = dot_scale
        self.font_weight = font_weight
        # Set configuration
        self.configuration = configuration
        # set color
        if color_sequence is None:
            self.color_sequence = "k"*len(self.sequence)
        elif len(color_sequence) != len(self.sequence):
            raise Exception("color_sequence must have the same length as sequence.")
        else:
            self.color_sequence = color_sequence
        # Sets rotation and configuration
        self.rotate(rotation)

    @property
    def data(self):
        """Return svg as a string."""
        return self.drawing.tostring()

    @property
    def notebook(self):
        """ Display SVG in Jupyther notebook. """
        try:
            # Import IPython display for notebook
            from IPython.display import SVG as ipython_display
            # Display in notebook
            return ipython_display(self.string)
        except ImportError:
            raise Warning(""" IPython not installed. """)

    def rotate(self, rotation):
        """Rotate the drawing by 90, 180, or 270.
        """
        # Rotate svg if desired.
        n = int(rotation/90)
        self.rotation += n
        for i in range(n):
            self.configuration = "".join([ROTATE[c] for c in self.configuration])
        self._build_drawing()

    def save(self, filename):
        """ save svg """
        self.drawing.saveas(filename)

    def _add_item(self, x, y):
        """Adds item at position (x,y) in array.
        """
        # Try if it's a bond
        try:
            if self.mapping[self.array[x,y]] is not None:
                # Get SVG element to add
                item = self.mapping[self.array[x,y]]
                # Add element
                item(2*self.offset+self.font_size*y,
                    2*self.offset+self.font_size*x)
        # Otherwise its a letter
        except KeyError:
            # Add color to specific letters if given
            color = COLORS[self.color_array[x,y]]
            self.drawing.letter(2*self.offset+self.font_size*y,
                2*self.offset+self.font_size*x,
                self.array[x,y],
                color=color)

    def _build_drawing(self):
        """Build drawing object."""
        self.array = np.array(configuration_to_array(self.sequence, self.configuration))
        self.color_array = np.array(configuration_to_array(self.color_sequence, self.configuration))
        # Build SVG grid object
        self.shape = self.array.shape
        self.height = self.font_size * self.shape[0]
        self.width = self.font_size * self.shape[1]
        self.drawing = Drawing(
            font_size=self.font_size,
            size=(self.width, self.height),
            dot_scale=self.dot_scale
        )
        self.offset = 0.25 * self.font_size
        
        # Object for how to draw a configuration
        self.mapping = {"d":self.drawing.down,
            "u":self.drawing.up,
            "r":self.drawing.right,
            "l":self.drawing.left,
            ".":self.drawing.dot,
            " ":None
        }
        # Draw grid in svg
        for x in range(self.shape[0]):
            for y in range(self.shape[1]):
                self._add_item(x,y)

class Drawing(svgwrite.Drawing):
    """Wrap svgwrite.Drawing object with extra methods that make drawing lattice
    proteins much easier
    """
    def __init__(self, font_size=20, dot_scale=1.0, font_weight="normal", **kwargs):
        self.font_size = font_size
        self.stepsize = 0.25*font_size
        self.linewidth = 0.1*font_size
        self.dot_scale = dot_scale
        self.font_weight = font_weight
        super(Drawing, self).__init__(**kwargs)

    def bond(self, start, end):
        """ Create an svg line object to add to a svgwrite drawing object.
        """
        line = self.line(start, end, stroke=svgwrite.rgb(10,10,16, '%'), style="stroke-width:" + str(self.linewidth))
        self.add(line)

    def right(self, x, y):
        """ Add a right facing line to figure. """
        start = (x-self.stepsize,y)
        end = (x+self.stepsize,y)
        self.bond(start,end)

    def left(self, x, y):
        """ Add a right facing line to figure. """
        start = (x+self.stepsize,y)
        end = (x-self.stepsize,y)
        self.bond(start,end)

    def up(self, x, y):
        """ Add a right facing line to figure. """
        start = (x,y-self.stepsize)
        end = (x,y+self.stepsize)
        self.bond(start, end)

    def down(self, x, y):
        """ Add a right facing line to figure. """
        start = (x,y+self.stepsize)
        end = (x,y-self.stepsize)
        self.bond(start,end)

    def letter(self, x, y, char, color="black"):
        """ Add letter to grid where residues exist. """
        xoffset = -self.stepsize-0.1*self.stepsize
        yoffset = self.stepsize
        letter = self.text(char, insert=(x + xoffset,y + yoffset),
            style="font-size:"+str(self.font_size)+
            "px;font-family:Courier;"
            "font-weight:" + self.font_weight + ";"
            "fill:" + str(color))
        return self.add(letter)

    def dot(self, x, y):
        """ Add a dot on grid where no letter exists"""
        xoffset = (-self.stepsize-0.15*self.stepsize)*self.dot_scale
        yoffset = (self.stepsize/2)
        letter = self.text('.', insert=(x + xoffset, y + yoffset),
            style="font-size:"+str(self.font_size * self.dot_scale)+
            "px;font-family:Courier")
        return self.add(letter)
