import sympy as sp


class Node:
    """A node in a scampy optical network.

    An instance of the `Node` class is created for each node added to a
    model. Nodes track the attached components, and hold the sympy symbols
    associated with the node.

    Attributes
    ----------
    name : str
            Unique name of the node.
    components : list
            List of components attached to node.
    symbols : tuple of sympy symbols
            Tuple of symbols associated with node.
    """

    def __init__(self, name):
        self.name = name
        self.components = []
        self.symbols = (sp.symbols(self.name + '_aP'),
                        sp.symbols(self.name + '_aS'),
                        sp.symbols(self.name + '_bP'),
                        sp.symbols(self.name + '_bS'))
