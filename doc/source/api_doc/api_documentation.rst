API documentation
=================================

`scampy` is structured into four modules. The `Model` module is the main entry
point for using `scampy`, holding the lists of optical components and nodes that
define the optical network, along with functions for building and evaluating the
model.

Optical components are described in the `components` module, this is currently a
monolithic file containing classes for each optical component.

The class describing optical detectors is held in the `Detector` module, all
results from the model are read out from instances of this class.

Systems of optical components in scampy are described by a network of
interconnected nodes, the basic node structure is defined in the `Node` module.
There should be no need for end users to interact with this module directly.

.. toctree::
   :maxdepth: 2
   :caption: Modules:

   model
   detector
   node
   components