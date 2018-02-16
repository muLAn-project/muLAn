Getting Started
===============

Description
-----------

*This section will be added soon.*

Installation
------------

First, download the last version of **muLAn** from GitHub: ::

   $ git clone https://github.com/muLAn-project/muLAn.git

Then, you can install **muLAn** localy as any standard python
package by doing: ::

   $ pip install --user muLAn/

It is recommended to use virtual environments to run muLAn. This
way, the code will manage automatically the proper package versions
without changing your current python distribution. There are several
ways to do so; we show how to do it with the widely spread
Anaconda_ python
distribution.

Create a new Anaconda_ virtual environment called *mulanenv*, activate it, and
install **muLAn**: ::

   $ conda create -n mulanenv python=2.7
   $ source activate mulanenv 
   $ pip install muLAn/

Following this procedure, the missing python packages will be installed automatically
within the virtual environment. Once you have finished to work with **muLAn**,
you can deactivate the corresponding environment: ::

   $ source deactivate mulanenv 

Dependencies
------------

This package depends in particular on the following packages:

- *Cython*
- PyAstronomy_ (version >= 0.10.1)
- *GetDist* (version >= 0.2.6)

First example
-------------

The package muLAn can be run in two different ways, depending on
user's preference.

- It can be used as a standard python package from a python script.

- It can also be used through a terminal in command line.

There is no difference between these two modes. We only provide a
flexible ready-to-use python script to allow a command line mode. 

An example is provided in the GitHub directory.

**Python style execution**

*This section will be added soon.*

**Command line style**

1. Open a terminal and move to the event directory: ::

   $ cd path-to-the-event

2. Edit the configuration files.

3. Run **muLAn** from the terminal: ::

   $ python mulan.py

It is possible to overwrite some parameters from the configuration files
by adding one or several options in the terminal. For example, the command::

   $ python mulan.py --stop

will cut short the exploration of the parameters space and properly stop
the code.

Collecting the Markov chains and create a CSV file from the best fitting
to the worse fitting model: ::

   $ python mulan.py --sort

A list of all options is printed by calling the help: ::

   $ python mulan.py -h

Here is the corresponding output: ::

   usage: python mulan.py [-h] [-a ARCHIVE] [-f] [--nchains NCHAINS]
                          [--ncores NCORES] [-o] [-p PLOT] [--resume] [-s] [-sn]
                          [--stop] [-v {0,1,2,3,4,5}]

   The command line options below overwrite the configuration file.

   optional arguments:
     -h, --help            show this help message and exit
     -a ARCHIVE, --archive ARCHIVE
                           Replace <ARCHIVE> by the name of the archive.
     -f, --fit             Ask muLAn to fit the data.
     --nchains NCHAINS     Replace <NCHAINS> by the number of chains.
     --ncores NCORES       Replace <NCORES> by the number of cores.
     -o, --optimize        Optimize inputs/outputs using pickle package when
                           possible.
     -p PLOT, --plot PLOT  Ask muLAn to plot the model number <PLOT> (default:
                           best fitting model). For several models, use coma
                           separator without space (e.g. -p 1,5,7) or a dash
                           without space for a range of models (e.g. -p 1-5).
     --resume              Resume an exploration.
     -s, --sort            Sort the samples and create a table.
     -sn, --sortno         Skip sort stage when running muLAn.
     --stop                Stop the exploration and save properly the results.
     -v {0,1,2,3,4,5}, --verbose {0,1,2,3,4,5}
                           Choose a verbose level.


.. References

.. _Anaconda: https://www.anaconda.com/what-is-anaconda/
.. _PyAstronomy: http://www.hs.uni-hamburg.de/DE/Ins/Per/Czesla/PyA/PyA/index.html
