QUANTBOX
========

The quantum sandbox ("quantbox") is a repository of Octave/Matlab code in the context of quantum information theory. See the source files for detailed inline documentation.


Installation
------------

Run `install_quantbox` to add all folders to the path. Use `savepath` if you want to automatically add the quantum sandbox folders to your Octave/Matlab search path.

Most functionality should work out of the box.

The semidefinite programming wrapper [`solve_sdp`](sdp/solve_sdp.m) requires either [SDPT3](http://www.math.nus.edu.sg/~mattohkc/sdpt3.html) or [SeDuMi](http://sedumi.ie.lehigh.edu/). We recommend [this distribution](https://github.com/sqlp), which works with both Octave and Matlab.

The self-tests use the [octave-doctest](https://github.com/catch22/octave-doctest) package. They can be run via `make test`.


References
----------

Previous versions of the quantum sandbox have been used in the following works:

- S. Seehars, Symmetric extensions in entanglement theory and quantum cryptography (2011, Master's thesis, ETH Zurich)
- M. Walter et al, Quantum state tomography of 1000 bosons (2012)
- A. Hansen, Swapped bound entanglement (2013, Master's thesis, ETH Zurich)
