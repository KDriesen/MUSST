**MUSST** is a set of modern fortran subroutines dedicated to the resolution of some lubrication problems.

Table of Contents
-----------------

- [Brief description](#brief-description)
- [Installation](#Installation)
- [Data-format](#data-format)
- [Limitations](#limitations)
- [Prerequesites](#prerequesites)
- [Wrappers](#wrappers)
- [Make](#make)
- [License](#license)

Brief description
-----------------

**MUSST** implements the finite element method on Reynolds' equation. It can solve a large variety of problems with standard geometries like sliders, bearings and pockets.
However, the user can provide the surface heights (as a ```.sur``` file for instance), whatever the surface.
The lubricant that can be modelled are perfect gases, incompressible fluids and compressible fluids (treated as gaz-fluid mixtures <a href="#foot01">[1]</a>). The latter allows for cavitation phenomena to appear.

The objective of the multiscale approach is to study flows between large rough surfaces needing very fine meshes while maintaining a reasonable computation time.
For this purpose, the domain is split into a number of subdomains (bottom-scale meshes) connected by a coarse mesh (top-scale).
The pressure distribution at the top-scale is used as boundary conditions for the bottom-scale problems.
This pressure is adjusted to ensure global mass flow balance between the contiguous subdomains.

This multiscale method allows for a significant reduction of the number of operations as well as a satisfactory accuracy of the results if the top-scale mesh is properly fitted to the roughness lateral scale.
Furthermore the present method is well-suited to parallel computation, leading to much more significant computation time reduction <a href="#foot01">[2]</a>.

[top](#table-of-contents)


Installation
--------------

**MUSST** needs to call efficient sparse linear system solvers like [MUMPS](http://mumps.enseeiht.fr/index.php?page=doc), [UMFPACK](http://faculty.cse.tamu.edu/davis/suitesparse.html),
[SuperLU](http://crd-legacy.lbl.gov/~xiaoye/SuperLU/#superlu) or [MA48](http://www.hsl.rl.ac.uk/catalogue/ma48.html). As the use of the latter is restricted, the *MSOLV* ```makefile``` has to be tuned accordingly.
If the user can supply the following files :

* ```MSOLV/src/hsl_common90.f90```
* ```MSOLV/src/hsl_ddeps90.f90```
* ```MSOLV/src/hsl_ma48d.f90```

and

* ```MSOLV/lib/libhsl_ma48.a```

then ```MA48_LIB = 1```, otherwise ```MA48_LIB = 0```.

To build the whole package: ```$cd MUSST``` and ```$make all```. To build **MUSST** core, ```$make```, or ```$make debug```, ```$make gprof```, ```$make clean```.

To launch some test examples:

* ```$./main cfg/TEST_01.dat``` -> slider bearing,  deterministic, incompressible
* ```$./main cfg/TEST_02.dat``` -> journal bearing, deterministic, compressible
* ```$./main cfg/TEST_04.dat``` -> rough surface,   deterministic, compressible
* ```$./main cfg/TEST_05.dat``` -> pockect bearing, deterministic, perfect gas
* ```$./main cfg/TEST_11.dat``` -> slider bearing,  multiscale,    compressible
* ```$./main cfg/TEST_14.dat``` -> rough surface,   multiscale,    compressible


[top](#Installation)

Data format
-----------

[top](#table-of-contents)

Limitations
-----------

[top](#table-of-contents)

Prerequesites
-------------
	
[top](#table-of-contents)

Wrappers
--------

[top](#table-of-contents)

Make
----

[top](#table-of-contents)

License
-------


<p id="foot01">[1]</p> Brunetière N. A General Model for Liquid and Gas Lubrication, Including Cavitation. *ASME. J. Tribol.* **2017**;*140*(2):021702-021702-10. [DOI](https://doi:10.1115/1.4037355)

<p id="foot01">[2]</p> Brunetière, N.; Francisco, A. Multiscale Modeling Applied to the Hydrodynamic Lubrication of Rough Surfaces for Computation Time Reduction. *Lubricants* **2018**, *6*, 83. [DOI](https://doi.org/10.3390/lubricants6030083)

