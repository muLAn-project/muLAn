# muLAn: gravitational MICROlensing Analysis code

Installation of muLAn:

1) Basic installation:

You will need to download and install separately the _pyAstronomy_ package (version 0.10.1)

<pre><code>git clone https://github.com/muLAn-project/muLAn.git
pip install --user muLAn/</code></pre>

2) Installation in a conda environment:

$ conda create -n yourenv python=2.7

$ source activate yourenv

$ pip install Cython

$ pip install muLAn/

3) Basic use of muLAn to model microlensing events:

In the event directory, run:

$ python test.py

where test.py contains:

  *from muLAn import mulan*

  *mulan.run()*


More inforation and examples to come!

The muLAn project team
