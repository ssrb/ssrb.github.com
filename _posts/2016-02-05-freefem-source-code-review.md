---
layout: post
title: "FreeFem++ Source Code Review"
description: ""
category: ""
tags: []
---
{% include JB/setup %}

In this post I will go through the source code of a tool I really like, [FreeFem++](http://www.freefem.org/ff++/), and document my findings.

<!-- more -->
Here is the structure of the `src` directory containing most if not all of `ff++` source code (excluding dependencies):

	ff++$ tree -d src
	src
	├── agl
	├── Algo
	├── bamg
	├── bamglib
	├── bin-win32
	├── Eigen
	├── femlib
	├── fflib
	├── glx
	├── Graphics
	├── iml
	├── lglib
	├── libMesh
	├── medit
	├── mpi
	├── nw
	├── solver
	├── std
	└── x11

The directories I'm particularly interested in are:

* `fflib`: lexer, abstract syntax tree data structures & interpreter, high level assembly and resolution of discretized PDE;
* `lglib`: the Bison grammar definition for the `ff++` language;
* `femlib`: low level assembly of discretized PDE, various FE families implementation (P1, P2, [Raviart-Thomas](https://en.wikipedia.org/wiki/Raviart%E2%80%93Thomas_basis_functions) ...) as well as matrix and vector definitions;
* `solver`: mainly wrappers to sparse direct solvers such as [UMFPACK](http://faculty.cse.tamu.edu/davis/suitesparse.html), [MUMPS](http://mumps.enseeiht.fr), [pardiso](http://www.pardiso-project.org/), etc.

Stuff I'm not too intersted in:

* The `bamg` and `bamglib` directories contain all the 2D mesh generation code: bamg stands for Bidimensional Anisotropic Mesh Generator and is the name of a stand alone 2D mesh generation tool written by [Frédéric Hecht](http://www.ann.jussieu.fr/hecht/), the main author of `ff++`. It comes with a library you can re-use in your own project.
The 3D mesh generation relies on [Tetgen](http://wias-berlin.de/software/tetgen);
* The `Algo`, `Eigen` and `iml` directories contain implementation of iterative methods for solving linear systems (CG, GMRES, ...) and find eigen values/vectors;
* The `libMesh` directory contains the code of the [libMesh library](https://www.rocq.inria.fr/gamma/gamma/Membres/CIPD/Loic.Marechal/Research/LM6.html) defining mesh datastructures and the corresponding read/write functions.

## Language

Work in progress !

## FE spaces and assembly

Work in progress !

### Numerical integration

Work in progress !

