muLAn: gravitational MICROlensing Analysis code
======

<!-- Commentaire <a href="https://travis-ci.org/muLAn-project/muLAn"><img src="https://travis-ci.org/muLAn-project/muLAn.svg?branch=master"></a> --> 

Goals
-----

<b>muLAn</b> is an easy and flexile software to model gravitational microlensing events

Getting Started
---------------

### Dependencies

This package depends on the following non-standard packages:

- `pyAstronomy` (version >= 0.10.1,  <a href="http://www.hs.uni-hamburg.de/DE/Ins/Per/Czesla/PyA/PyA/index.html">link</a>)
- `GetDist` (version >= 0.2.6)

### Installation

To install the current development version of <b>muLAn</b> package from source:

```
$ git clone https://github.com/muLAn-project/muLAn.git
$ pip install --user muLAn/
```

### Loading

```python
from muLAn import mulan
```

### Launcing muLAn

In the directory of the event to be anayzed, run:
```python
from muLAn import mulan
mulan.run()
```

Examples
--------

See the documentation section.

License
-------

This software is licensed under the MIT license. See the [LICENSE](LICENSE) file
for details.

