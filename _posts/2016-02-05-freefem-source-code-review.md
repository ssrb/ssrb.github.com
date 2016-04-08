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

Given a concise description of a PDE weak formulation / variational form as well as a set of boudary conditions, `ff++` will automagically discretize and solve the PDE.
Thanks to its domain specific language it's also possible to build more advance numerical methods on top of the FE solver itself (domain decomposition, multi-physics, ...).

This is a good tool if you're tired of writting ad-hoc solvers from scratch each time you need to perform a FE analysis.
Of course before using such a tool it's nice to clearly understand what's going on under the hood.

Similar tools are

* [FEnics](http://fenicsproject.org);
* [the MATLAB PDE toolbox](http://au.mathworks.com/products/pde/?requestedDomain=www.mathworks.com)

## Codebase srtucture

Here is the structure of the `src` directory containing most of `ff++` source code (excluding dependencies):

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
* `lglib`: the grammar definition for the `ff++` language;
* `femlib`: low level assembly of discretized PDE, various FE families implementation (P1, P2, [Raviart-Thomas](https://en.wikipedia.org/wiki/Raviart%E2%80%93Thomas_basis_functions) ...) as well as matrix and vector definitions;
* `solver`: mainly wrappers to sparse direct solvers such as [UMFPACK](http://faculty.cse.tamu.edu/davis/suitesparse.html), [MUMPS](http://mumps.enseeiht.fr), [pardiso](http://www.pardiso-project.org/) as well as sparse iterative solvers such as [HYPRE](https://computation.llnl.gov/project/linear_solvers/index.php), [HYPS](http://hips.gforge.inria.fr/), [pARM](http://www-users.cs.umn.edu/~saad/software/pARMS) etc.

Stuff I'm not too intersted in:

* The `bamg` and `bamglib` directories contain all the 2D mesh generation code: bamg stands for Bidimensional Anisotropic Mesh Generator and is the name of a stand alone 2D mesh generation tool written by [Frédéric Hecht](http://www.ann.jussieu.fr/hecht/), the main author of `ff++`. It comes with a library you can re-use in your own project.
The 3D mesh generation relies on [Tetgen](http://wias-berlin.de/software/tetgen);
* The `Algo`, `Eigen` and `iml` directories contain implementation of iterative methods for solving linear systems (CG, GMRES, ...) and find eigen values/vectors;
* The `libMesh` directory contains the code of the [libMesh library](https://www.rocq.inria.fr/gamma/gamma/Membres/CIPD/Loic.Marechal/Research/LM6.html) defining mesh datastructures and the corresponding read/write functions.

## Language

`ff++` language is imperative. It borrows some syntax from `C++`. For example the stream operators (`<<` and `>>`) to perform formated I/Os which is not really the best part of `C++` ...

### Lexer

The hand written lexer (`class mylex`) lives in `fflib/lex.hpp`.

### Parser

The parser is using [GNU Bison](https://www.gnu.org/software/bison). Its grammar is described in the `lglib/lg.ypp` file. I guess the `lg` prefix stands for "language". That prefix is used in other places in the code base. For example: `fflib/lgfem.hpp`, `fflib/lgmesh.hpp`, `fflib/lgmat.hpp`, `fflib/lgmesh.hpp`, etc where the language builtin types and functions are defined. This convention doesn't seem to be enforced.

The generated parser entry point is called `yyparse()` (it's #defined to `lgparse`)

The first frames of a call stack look like this:

	#0  yyparse () at lg.tab.cpp:1851
	#1  Compile () at lg.ypp:777
	#2  mainff (argc=2, argv=0x7fffffffddf8) at lg.ypp:947
	#3  main (argc=2, argv=0x7fffffffde08) at ../Graphics/sansrgraph.cpp:199

Now, even though `yyparse()` is called by a function called `Compile()` and `ff++` outputs messages such as 

	times: compile 0.100339s, execution 89.3461s

the language is **not** compiled: no machine or byte code is generated whatsoever. It's more of an interpreter.

### Interpreter

In `ff++`, the [Bison semantic actions](https://www.gnu.org/software/bison/manual/bison.html#Semantic-Actions) do not usually interprete statements on the fly. 
Instead an abstract syntax / expression tree is built using types defined in `fflib/AFunction.hpp`, `fflib/AFunction2.hpp`, `fflib/Operator.hpp`, etc. 

An exception would be the `load` statement which will trigger a call to `bool load(string ss)` defined in `fflib/load.cpp`, loading shared objects ("plugins") on the fly.

In the top level semantic action, executed once the end of file is reached, the tree is walked and statements are interpreted:

{% highlight C++ %}
start:   input ENDOFFILE {
	[...]	
	$1.eval(stack);
	[...]
}
{% endhighlight %}

In this Bison/C++ code snipet `S1` is the result of the `input` semantic action and is of type `class CListOfInst` defined in `fflib/AFunction.hpp`.

Ultimately the call to `CListOfInst::eval` will execute:

{% highlight C++ %}
AnyType ListOfInst::operator()(Stack s) const {     
    AnyType r; 
    double s0=CPUtime(),s1=s0,ss0=s0;
    StackOfPtr2Free * sptr = WhereStackOfPtr2Free(s);
    try { // modif FH oct 2006 
	for (int i=0;i<n;i++) 
	{
	    TheCurrentLine=linenumber[i];
	    r=(*list[i])(s);
	    sptr->clean(); // modif FH mars 2006  clean Ptr
	    s1=CPUtime();
	    if (showCPU)  
		cout << " CPU: "<< i << " " << s1-s0 << "s" << " " << s1-ss0 << "s" << endl;
	    s0=CPUtime();
	}
    }
    catch( E_exception & e) 	
    {
	if (e.type() != E_exception::e_return)  
	    sptr->clean(); // pour ne pas detruire la valeur retourne  ...  FH  jan 2007
	throw; // rethow  
    }
    catch(...)
    {
	sptr->clean();
	throw; 
    }
    return r;
}
{% endhighlight %}

This method can be seen as the interpreter main loop I guess.

By the way, often the code is mixing (broken) English and French which is not making it easier to dive into the codebase.

After parsing is complete and the interpreter starts doing its job, the stack should look like this:

	#0  ListOfInst::operator() at AFunction2.cpp:793
	#1  eval at ./../fflib/AFunction.hpp:1459
	#2  yyparse () at lg.tab.cpp:1851
	#3  Compile () at lg.ypp:777
	#4  mainff (argc=2, argv=0x7fffffffddf8) at lg.ypp:947
	#5  main (argc=2, argv=0x7fffffffde08) at ../Graphics/sansrgraph.cpp:199

### Activation record

Storage for variables within the current lexical scope is allocated on a stack of type `struct StackType` defined in `fflib/ffstack.hpp` (beware that `typedef StackType & Stack;`).

### Builtins

## FE spaces and assembly

Work in progress !

### Numerical integration

Work in progress !

