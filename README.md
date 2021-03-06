muLAn: gravitational MICROlensing Analysis Software
======

[![Build Status](https://travis-ci.org/muLAn-project/muLAn.svg?branch=py3)](https://travis-ci.org/muLAn-project/muLAn)
[![License](https://img.shields.io/github/license/muLAn-project/muLAn)](https://github.com/muLAn-project/muLAn/blob/master/LICENSE)
[![Build status](https://ci.appveyor.com/api/projects/status/kyufdmvy3ske3qlb?svg=true)](https://ci.appveyor.com/project/clementranc/mulan)

Goals
-----

**muLAn** is an easy-to-use and flexile software to model gravitational microlensing events

Getting Started
---------------

### Dependencies

This package depends in particular on the following packages:

- `Cython`
- `pyAstronomy` (version >= 0.10.1,  <a href="http://www.hs.uni-hamburg.de/DE/Ins/Per/Czesla/PyA/PyA/index.html">link</a>)
- `GetDist` (version >= 0.2.6)

### Installation

Install one of the available relesases (recommended), or install the current development version of **muLAn** package from source:
```
git clone https://github.com/muLAn-project/muLAn.git
```

To install **muLAn** in your user python environment: 
```
pip install --user muLAn/
```

We however recommend to use an [Anaconda](https://anaconda.org) user environment:
```
conda create -n muLAn python=3.8
source activate muLAn
pip install muLAn/
```

### Loading

The python script has to be executed in the directory of the microlensing event to be analyzed:

```python
from muLAn import mulan
```

### muLAn basic commands

- Launch **muLAn** (based on `setup.ini` options, see documentation):

```python
mulan.run()
```

- Clean stop (waiting for all current sub-processes to finish):

```python
mulan.stop()
```

- Collect all Markov chains and re-order them with best fit first (necessary after unexpected stop, e.g. wall time exceeded):

```python
mulan.sortmodels()
```

Examples
--------

See the documentation section.

License
-------

This software is licensed under the MIT license. See the [LICENSE](LICENSE) file
for details.

