name = 'muLAn'

import sys
import os
import platform

from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
import setuptools

if platform.system() == 'Darwin':
    macosv = platform.mac_ver()[0]
    macosv_full = [int(a) for a in macosv.split('.')]
    if (macosv_full[0] >= 10) & (macosv_full[1] > 9):
        txt = "-mmacosx-version-min={:d}.{:d}".format(macosv_full[0], macosv_full[1])
        txt = "-mmacosx-version-min=10.9"
        extensions = [Extension(name="muLAn.models.vbb.vbb",
                                sources=["muLAn/models/vbb/vbb.pyx"],
                                language="c++",
                                extra_link_args=["-stdlib=libc++", txt])]
    else:
        extensions = [Extension(name="muLAn.models.vbb.vbb",
                                sources=["muLAn/models/vbb/vbb.pyx"],
                                language="c++")]
else:
        extensions = [Extension(name="muLAn.models.vbb.vbb",
                            sources=["muLAn/models/vbb/vbb.pyx"],
                            language="c++")]

setup(
      name = name,
      ext_modules = cythonize(extensions),
      data_files=[('lib/python2.7/site-packages/muLAn/models/espl_inputs', ['muLAn/models/espl_inputs/tab_Bo.p', 'muLAn/models/espl_inputs/tab_B1.p'])]
      )

pjoin = os.path.join
here = os.path.abspath(os.path.dirname(__file__))

packages = []
for d, _, _ in os.walk(pjoin(here, name)):
    if os.path.exists(pjoin(d, '__init__.py')):
        packages.append(d[len(here)+1:].replace(os.path.sep, '.'))

version_ns = {}
with open(pjoin(here, name, '_version.py')) as f:
    exec(f.read(), {}, version_ns)

setup_args = dict(
    name            = name,
    version         = version_ns['__version__'],
    packages        = packages,
    description     = "muLAn: gravitational MICROlensing Analysis code",
    long_description= "muLAn: gravitational MICROlensing Analysis code",
    author          = 'Clement Ranc, Arnaud Cassan',
    author_email    = 'clement.ranc@nasa.gov, cassan@iap.fr',
    url             = 'https://github.com/muLAn-project/muLAn.git',
    license         = 'MIT',
    platforms       = "Linux, Mac OS X, Windows",
    keywords        = ['Astronomy', 'Microlensing', 'Science'],
    classifiers     = [
        'Intended Audience :: Developers',
        'Intended Audience :: System Administrators',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT',
        'Programming Language :: Python',
        'Programming Language :: Python :: 2.7',
                       ],
)

if 'develop' in sys.argv or any(a.startswith('bdist') for a in sys.argv):
    import setuptools

setuptools_args = {}

install_requires = setuptools_args['install_requires'] = [
                                                          'scikit-learn>=0.19.0',
                                                          'argparse>=1.1',
                                                          'astropy>=1.3.2',
                                                          'bokeh>=0.12.4',
                                                          'configparser>=3.5.0',
                                                          'Cython>=0.27.3',
                                                          'emcee>=2.2.1',
                                                          'GetDist==0.2.6',
                                                          'matplotlib>=2.1.2',
                                                          'numpy>=1.12.1',
                                                          'pandas>=0.18.1',
                                                          'PyAstronomy>=0.10.1',
                                                          'scipy>=0.18.1'
]

extras_require = setuptools_args['extras_require'] = {
}

if 'setuptools' in sys.modules:
    setup_args.update(setuptools_args)

if __name__ == '__main__':
    setup(**setup_args)
