.. muLAn documentation master file, created by
   sphinx-quickstart on Sat Apr 23 19:04:35 2016.

muLAn: gravitational MicroLensing Analysis code
===============================================

**muLAn** is an easy-to-use and flexible software to model gravitational
microlensing events. Why using muLAn?

- The code combines the-state-of-the-art numerical methods to compute
  the magnification. These algorithms are direclty plugged into muLAn
  from their original versions.

- The code uses several methods to explore the parameters space: a
  grid search method, a differential MCMC algorithm based on the package
  EMCEE, and a combination of these two methods are implemented.

- The code is designed to allow a fast integration of any additional
  numerical method independently developed. muLAn provides a powerfull
  tool to compare several algorithms on a same light curve.

- The code has been designed to be easy-to-use and easy-to-install by
  any interested astrophysicist. In particular, muLAn may be run
  through a python sript like any python package, or through command
  line instructions in a terminal.

- The code and the documentation are tested and maintained by several
  researchers.

Contents:

.. toctree::
   :maxdepth: 2

   contents/introduction
   contents/tutorial
   contents/configfiles
   contents/plotroutines
   contents/statistics


.. Indices and tables
   ==================

   * :ref:`genindex`
   * :ref:`modindex`
   * :ref:`search`


.. Installation
   ============

   *Create a directory.*
