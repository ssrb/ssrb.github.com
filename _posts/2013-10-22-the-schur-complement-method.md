---
layout: post
title: "The Schur Complement Method: Part 1"
description: "An easy way to solve finite element problems in parallel"
category: "Numerical methods"
tags: [mathematics, parallel computing]
assets: assets/schur
---

{% include JB/setup %}

This post is going to illustrate the primal Schur complement method briefly described [here](http://en.wikipedia.org/wiki/Schur_complement_method).
The Schur complement method is a strategy one can use to divide a finite element problem into independant sub-problems.
It's not too involved but requires good understanding of block Gaussian elimination, reordering degrees of freedom plus a few "tricks of the trade"
to avoid computing inverse of large sparse matrices.

<!-- more -->

## A (finite element) problem

In that part, we're going to create all the data needed to implement the method:

* a mesh
* a system of linear equations

### A 2D Poisson's equation

First step is to choose a problem.
For example, solve a [Poisson equation](http://en.wikipedia.org/wiki/Poisson's_equation).

Given a domain `\(\Omega\)`, the problem is to find `\(\varphi\)` such that
$$\begin{cases}\left( \frac{\partial^2}{\partial x^2} + \frac{\partial^2}{\partial y^2}\right)\varphi(x,y) & = & 1, & \mbox{on } \Omega \\\ \varphi(x,y) & = & 0, & \mbox{on } \partial\Omega\end{cases}$$

I chose a L-shape for `\(\Omega\)`:

![Alt The Domain]({{ site.url }}/{{ page.assets }}/poisson2D.domain.png)

### FreeFem++

Since this post doesn't aim to illustrate in great details the [finite element method](http://en.wikipedia.org/wiki/Finite_element_method) itself, I'm going to translate the [weak formulation](http://en.wikipedia.org/wiki/Weak_formulation) of the problem into the [FreeFem++](http://www.freefem.org/ff++/index.htm) language.
What FreeFem++ is going to do for us is to create the desired mesh and system of linear equations.

The weak formulation of the problem is to find `\(u\)` such that `\(\int_\Omega \left( \frac{\partial u}{\partial x}\frac{\partial v}{\partial x} + \frac{\partial u}{\partial y}\frac{\partial v}{\partial y} \right)dxdy + \int_\Omega vdxdy = 0, \forall v\)`. 

Note that I will be using Lagrange P1 elements to make reordering of the degrees of freedom (DoF) easier in the next step.

This is how the problem translates into the FreeFem++ langage:
{% highlight c linenos %}
// Describe the geometry of the domain
border red(t=0,1){x=t;y=0;};
border green(t=0,0.5){x=1;y=t;};
border yellow(t=0,0.5){x=1-t;y=0.5;};
border blue(t=0.5,1){x=0.5;y=t;};
border orange(t=0.5,1){x=1-t;y=1;};
border violet(t=0,1){x=0;y=1-t;};

// Build a mesh over that domain
mesh Th = buildmesh (red(6) + green(4) + yellow(4) + blue(4) + orange(4) + violet(6));

// Define the finite element space, using Lagrange P1 element
fespace Vh(Th,P1);

// Define u, the function we are looking for, as well as v, the test function in the weak formulation.
Vh u=0,v;

// Define the value of the Laplacian of phi within the domain
func f= 1;

// Define the value of phi on the border of the domain
func g= 0;

// Describe the problem using the weak formulation
problem Poisson2D(u,v) =
    int2d(Th)(dx(u)*dx(v) + dy(u)*dy(v)) // First integral
  + int2d(Th) (v*f)  // Second integral
  + on(red,green,yellow,blue,orange,violet,u=g); // What happens on the border
{% endhighlight %}

Given that script, FreeFem++ is going to:

* build a mesh;
* build a system of linear equations: a square matrix `\(A\)` and a column vector `\(b\)`. Each DoF (unknown) is the value of the function we are looking for at a given vertex of the mesh (because we choosed Lagrange P1 elements).

Hereafter is the mesh (built by FreeFem++) I am going to work with in the next steps:
![Alt Mesh]({{ site.url }}/{{ page.assets }}/poisson2D.mesh.png)

It is made of 7215 triangles and 3735 vertices.

* [Raw Mesh data]({{ site.url }}/{{ page.assets }}/poisson2D.mesh)

Hereafter is the sparsity pattern of the corresponding `\(A\)` matrix:
![Alt Spy A]({{ site.url }}/{{ page.assets }}/spyA.png)

Its size is 3735x3735. It's a band matrix because of the way FreeFem++ numbered the DoF.

* [Raw A matrix data]({{ site.url }}/{{ page.assets }}/poisson2D.A)
* [Raw b column vector data]({{ site.url }}/{{ page.assets }}/poisson2D.b)

### Conclusion and spoiler

So there we are, we got all the data needed to get started.
What needs to be done next is:

* partition the domain
* reorder the DoF
* implement the method itself
* celebrate

By the way, if you want to have a look at the solution to the problem, click [here]({{ site.url }}/{{ page.assets }}/poisson2D.sol.filled.png)

## Partitioning the domain

The next stage to illustrate the method is to partition the initial domain into two sub-domains.
In order to do so, one can use specialized tools such as [Chaco](http://www.sandia.gov/~bahendr/chaco.html)
 or [Metis](http://glaros.dtc.umn.edu/gkhome/views/metis). I'm arbitrarily going to use Metis.

### Metis

Metis input is a graph and a number of sub-domains. It doesn't use the coordinates of the vertices.
I generated Metis input graph from the mesh using this one liner:
{% highlight bash%}
perl -ne 'push @T, "$1 $2 $3\n" if /(\d+) (\d+) (\d+) (\d+)/; END{print @T."\n"; print for @T}' poisson2D.mesh > poisson2D.metis_graph
{% endhighlight %}

* [Raw Metis input graph data]({{ site.url }}/{{ page.assets }}/poisson2D.metis_graph)

The command `mpmetis poisson2D.metis_graph 2` will output two files:

* [poisson2D.metis_graph.epart.2]({{ site.url }}/{{ page.assets }}/poisson2D.metis_graph.epart.2)
* [poisson2D.metis_graph.npart.2]({{ site.url }}/{{ page.assets }}/poisson2D.metis_graph.npart.2)

`poisson2D.metis_graph.epart.2` (respectively `poisson2D.metis_graph.npart.2`) contains the domain indices (starting from 0) of each
triangle (respectively vertex) of the mesh. For example, line 42 of `poisson2D.metis_graph.epart.2` will contain a 0 if triangle number
42 belongs to the domain number 0.

### Identifying the boundary vertices

Sadly, the Metis output doesn't tell us what are the boundary vertices so that we need to do a little bit of post-processing.
Remember that we need to identify the boundary DoFs to pack them into a single matrix block in a later stage.
The script hereafter is going to fan out vertices ids into separate files,`domain1.vids` (vertices within domain 1 exclusively), `domain2.vids` and `interface.vids` (vertices on the interface).

{% highlight perl linenos %}
#! /usr/bin/perl
use strict;
use warnings;

open my $FDTRIANGLES, "poisson2D.metis_graph";
open my $FDIDS, "poisson2D.metis_graph.epart.2";

my @vertexToDomain;
my %interface;

while (<$FDTRIANGLES>) {
	if (/(?<vid1>\d+) (?<vid2>\d+) (?<vid3>\d+)/) {
		my $domain = <$FDIDS> + 1;
		CheckVertex($+{vid1}, $domain);
		CheckVertex($+{vid2}, $domain);
		CheckVertex($+{vid3}, $domain);
	}
}

my @DOMAINFD;
open $DOMAINFD[1], "> domain1.vids";
open $DOMAINFD[2], "> domain2.vids";
open my $INTERFACEFD, "> interface.vids";
while (my ($vid, $domain) = each @vertexToDomain) {
	next if !$vid;
	my $fd = exists $interface{$vid} ? $INTERFACEFD : $DOMAINFD[$domain];
	print $fd $vid."\n";
}

# Count the number of domain a vertex belongs to, if it's greater than 1
# that vertex is on the interface
sub CheckVertex {
	my ($vid, $domain) = @_;
	if ($vertexToDomain[$vid] && $vertexToDomain[$vid] != $domain) {
		$interface{$vid}++;
	} else {
		$vertexToDomain[$vid] = $domain;
	}
}
{% endhighlight %}

* [domain1.vids]({{ site.url }}/{{ page.assets }}/domain1.vids)
* [domain2.vids]({{ site.url }}/{{ page.assets }}/domain2.vids)
* [interface.vids]({{ site.url }}/{{ page.assets }}/interface.vids])

### Conclusion

Finally, here is what the partition looks like:

![Alt Mesh]({{ site.url }}/{{ page.assets }}/poisson2D.part2.png)

Domains 1 and 2 respectively contain 1839 and 1834 DoFs. The interface contains 62 DoFs.
Tools such as Metis are doing a good job at balancing the size of the sub-domains while minimizing the size of the interfaces.

## Reordering the degrees of freedom

Now that I partitioned the domain, we have everything needed to give the sytem of linear equations the structure needed by the Schur complement method:

$$\begin{equation}\left[\begin{matrix} A_{11} & 0 & A_{1\Gamma} \\\ 0 & A_{22} & A_{2\Gamma} \\\ A_{\Gamma 1} & A_{\Gamma 2} & A_{\Gamma\Gamma}\end{matrix}\right]\left[\begin{matrix} U_1 \\\ U_2 \\\ U_\Gamma\end{matrix}\right] = \left[\begin{matrix} b_1 \\\ b_2 \\\ b_\Gamma\end{matrix}\right]\label{eq:mainsystem}\end{equation}$$

where

* `\(A_{ii}\)` is the block coupling the interior of domain `\(i\)` to itself;
* `\(A_{\Gamma\Gamma}\)` is the block coupling the interface to itself;
* `\(A_{\Gamma i}\)` is the block coupling the interface to the interior of domain `\(i\)`;
* `\(A_{i\Gamma}\)` is the block coupling the interior of domain `\(i\)` to the interface

### Octave

From now on, I will be using [Octave](http://www.gnu.org/software/octave).

To do the reordering I'm going to create a permutation matrix `\(P\)` out of the `domain1.vids`, `domain2.vids` and `interface.vids` and apply it to `\(A\)` (`\(PAP^\top\)`) and `\(b\)` (`\(Pb\)`) from the first step.

{% highlight octave linenos %}
load("-ascii", "domain1.vids");
load("-ascii", "domain2.vids");
load("-ascii", "interface.vids");
q = [domain1; domain2; interface];
A = A(q,q);
b = b(q);
{% endhighlight %}

Hereafter is the spy of the desired block matrix after applying the permutation:
![Alt Spy reordered A]({{ site.url }}/{{ page.assets }}/spy_reordered1.png)

### Matrix pr0n:

`\(A_{11}\)`, `\(A_{22}\)` sparsity patterns are similar to the initial `\(A\)` matrix one since we took care of preserving the relative order
of DoFs in the interior of the domains.

`\(A_{11}\)`: ![Alt Spy A11]({{ site.url }}/{{ page.assets }}/spyA11.png)
`\(A_{22}\)`: ![Alt Spy A22]({{ site.url }}/{{ page.assets }}/spyA22.png)

However `\(A_{TT}\)` sparsity pattern doesn't look great considering that the geometry here is a simple polyline. The initial DoF ordering
does not make sens on the interface. I'm going to address this in the next section.

`\(A_{TT}\)`: ![Alt Spy ATT]({{ site.url }}/{{ page.assets }}/spyATT.png)

`\(A_{T1}\)`: ![Alt Spy A11]({{ site.url }}/{{ page.assets }}/spyAT1.png)

`\(A_{T2}\)`: ![Alt Spy A11]({{ site.url }}/{{ page.assets }}/spyAT2.png)


### `\(A_{TT}\)` reloaded
I'm going to improve `\(A_{TT}\)` sparsity pattern.

The idea is to walk the interface from one end to the other, dumping DoF ids as we go.
This new permutation is going to preserve locality and will greatly improve the sparsity pattern of `\(A_{TT}\)`.

{% highlight perl linenos %}
#! /usr/bin/perl
use strict;
use warnings;

open my $FDTRIANGLES, "poisson2D.metis_graph";
open my $FDIDS, "poisson2D.metis_graph.epart.2";

my %triangleToVertices;
my %vertexToTriangles;
my %interface;
my %interfaceGraph;

open my $INTERFACEFD, "interface.vids";
while (<$INTERFACEFD>) {
	$interface{int($_)}++;
}

# Build connectivity and reverse connectivity tables for 
# domain 1 elements on the interface
my $triangleId = 0;
while (<$FDTRIANGLES>) {
	if (/(?<vid1>\d+) (?<vid2>\d+) (?<vid3>\d+)/) {
		++$triangleId;
		my $domain = <$FDIDS> + 1;
		if ($domain == 1) {
			AddVertices($triangleId, $1, $2, $3);
			AddTriangle($+{vid1}, $triangleId);
			AddTriangle($+{vid2}, $triangleId);
			AddTriangle($+{vid3}, $triangleId);
		}
	}
}

# Go through each edge and check if it's on the interface
while (my ($tid, $vids) = each %triangleToVertices) {
	for (0..$#{$vids}) {
		my ($vid1, $vid2) = @{$vids}[$_, ($_ + 1) % @{$vids}];
		if (IsInterfaceEdge($vid1, $vid2)) {
			AddInterfaceEdge($vid1, $vid2);
		}
	}
}

# Walk the interface
my $start = FindInterfaceEnd();
WalkInterface($start, $start);

sub AddVertices {
	my ($tid, $vid1, $vid2, $vid3) = @_;
	if (exists $interface{$vid1} || exists $interface{$vid2} || exists $interface{$vid3}) {
		$triangleToVertices{$tid} = [$vid1, $vid2, $vid3];
	}
}

sub AddTriangle {
	my ($vid, $tid) = @_;
	if (exists $interface{$vid}) {
		$vertexToTriangles{$vid}{$tid}++;
	}
}

sub AddInterfaceEdge {
	my ($vid1, $vid2) = @_;
	$interfaceGraph{$vid1}{$vid2}++;
	$interfaceGraph{$vid2}{$vid1}++;
}

# An interface edge has its two vertices on the interface 
# and belongs to one and only one domain 1 triangle.
sub IsInterfaceEdge {
	my ($vid1, $vid2) = @_;

	return 0 unless exists $interface{$vid1} && exists $interface{$vid2};

	my @tids1 = keys %{$vertexToTriangles{$vid1}};
	my @tids2 = keys %{$vertexToTriangles{$vid2}};

	my %union;
	my $triangleCount = 0;
	for my $tid (@tids1, @tids2) {
		if (++$union{$tid} == 2) {
			$triangleCount++;
		}
	}

	return $triangleCount == 1;
}

# The interface should have two ends.
# An end is simply a vertex connected to only one other vertex.
sub FindInterfaceEnd {
	while (my ($from, $tos) = each %interfaceGraph) {
		if (keys %{$tos} == 1) {
			return $from;
		}
	}
	exit -1;
}

# Perform a simple DFS to walk the interface 
# and print the vertex ids as we go.
sub WalkInterface {
	my ($u, $v) = @_;
	print $v."\n";
	for (keys %{$interfaceGraph{$v}}) {
		next if $_ == $u;
		return WalkInterface($v, $_);
	}
}
{% endhighlight %}

[interface2.vids]({{ site.url }}/{{ page.assets }}/interface2.vids])

And ... voila:

`\(A_{TT}\)`: ![Alt Spy ATT]({{ site.url }}/{{ page.assets }}/spyATT2.png)

If you're wondering how can one interface DoF be coupled to 4 other DoFs, this happens when
one interface DoF belongs to 2 triangles and each of these triangles has all its
vertices on the interface.

### Conclusion

We're done doing the reordering. Here is the final version:

![Alt Spy A reordered final]({{ site.url }}/{{ page.assets }}/spy_reordered21.png)
![Alt Spy A reordered final detail ATT]({{ site.url }}/{{ page.assets }}/spy_reordered22.png)

We're now ready to define a few more variables in Octave:

{% highlight octave linenos %}
i1 = length(domain1);
i2 = i1 + length(domain2);

A11 = AA(1:i1, 1:i1);
A1T = AA(1:i1, i2+1:end);
A22 = AA(i1+1:i2, i1+1:i2);
A2T = AA(i1+1:i2, i2+1:end);
AT1 = AA(i2+1:end, 1:i1);
AT2 = AA(i2+1:end, i1+1:i2);
ATT = AA(i2+1:end, i2+1:end);
U = zeros(length(A), 1);
clear A;

F1 = F(1:i1);
F2 = F(i1+1:i2);
FT = F(i2+1:end);
clear F;
{% endhighlight %}

## Solving on the interface

### Block Gaussian elimination

Being able to solve first for the interface without even considering what's happening in the interior of the subdomains 
might seem magical: this is nothing more than block Gaussian elimination applied to the matrix we crafted
in the previous stages. In `\(\eqref{eq:mainsystem}\)`, multiply the first line by `\(A_{T1}A_{11}^{-1}\)`, the second line by `\(A_{T2}A_{22}^{-1}\)`,
substract both quantities from the third line and you end up with:

$$(A_{TT} - A_{T1}A_{11}^{-1}A_{1T} - A_{T2}A_{22}^{-1}A_{2T})U_{T} = F_{T} - A_{T1}A_{11}^{-1}F_{1} - A_{T2}A_{22}^{-1}F_{2}$$

In other words, we need to solve:

$$\begin{equation}\Sigma U_{T} = \tilde{b}\label{eq:interfacesystem}\end{equation}$$

where `\(\Sigma = A_{TT} - A_{T1}A_{11}^{-1}A_{1T} - A_{T2}A_{22}^{-1}A_{2T}\)`
and `\(\tilde{b} = F_{T} - A_{T1}A_{11}^{-1}F_{1} - A_{T2}A_{22}^{-1}F_{2}\)`

### Schur complement

In the numerical analysis lingo, `\(\Sigma\)` is known as the Schur complement of `\(A_{TT}\)` in `\(A\)`.

As we can see, both `\(\Sigma\)` and `\(\tilde{b}\)` depend on `\(A_{11}^{-1}\)` and `\(A_{22}^{-1}\)`.
However, `\(A_{11}\)` and `\(A_{22}\)` are large matrices we should try not to invert.

Anyway, let's explicitely compute the Schur complement for our baby problem and have a look:

{% highlight octave linenos %}
% Bad bad bad !
spy(ATT - AT1 * inv(A11) * A1T + AT2 * inv(A22) * A2T);
{% endhighlight %}

It's crazy dense!

![Alt Spy Schur]({{ site.url }}/{{ page.assets }}/spySchur.png)

Furthermore Octave tells us that `\(A_{11}\)` and `\(A_{22}\)` are singular to machine precision:
`warning: inverse: matrix singular to machine precision, rcond = 1.65116e-30`

### Iterative method

To summarize we would like to avoid:


1. explicitely computing `\(\Sigma\)` in order to solve `\(\eqref{eq:interfacesystem}\)`;
2. explicitely inverting `\(A_{11}\)` or `\(A_{22}\)`


If we try to satisfy these two constraints, it is clear that we can't use a direct method to solve `\(\eqref{eq:interfacesystem}\)`.
On the other hand, if we choose an iterative method such as the preconditioned conjugate gradient method (PCG),
we just need to mulitply, each iteration of the method, a current solution `\(u\)` by `\(\Sigma\)`. Multiplying by `\(\Sigma\)` is 
a different story than explicitely computing `\(\Sigma\)`:


$$\Sigma u = A_{TT}u - A_{T1}A_{11}^{-1}A_{1T}u - A_{T2}A_{22}^{-1}A_{2T}u$$


### Decoupled Dirichlet problems

I'm now going to describe a simple tactic to evaluate `\(A_{11}^{-1}A_{1T}u\)` and `\(A_{22}^{-1}A_{2T}u\)` without inverting `\(A_{11}\)` and `\(A_{22}\)`.

Assume you're given a very large invertible matrix `\(A\)` and one (1) vector `\(b\)` and you're tasked to compute the quantity `\(A^{-1}b\)`.
Would you explicitely compute `\(A^{-1}\)` ? 
Probably not.
What you would do instead is to find `\(x\)` such that `\(Ax=b\)`.
Then `\(x=A^{-1}b\)`, which is what you were tasked to compute.

Following this idea, I'm going to use a direct method to compute these quantities: I will compute the LU decomposition of `\(A_{11}\)` (respectively `\(A_{22}\)`) once and for all, and every PCG iteration, compute `\(b = A_{1T} * u\)` (respectively `\(b = A_{2T} * u\)`) and solve `\(A_{11}x=b\)` (respectively `\(A_{22}x=b\)`). The same tactic is used to compute `\(A_{11}^{-1}F_{1}\)` and `\(A_{22}^{-1}F_{2}\)` in `\(\tilde{b}\)`.

Each iteration of the PCG, we end up solving decoupled [Dirichlet problems](http://en.wikipedia.org/wiki/Dirichlet_problem) on each domain.

All together, this looks like this:
{% highlight octave linenos %}
[L11, U11, p11, q11] = lu(A11, 'vector');
q11t(q11) = 1:length(q11);
clear A11;
clear q11;

[L22, U22, p22, q22] = lu(A22, 'vector');
q22t(q22) = 1:length(q22);
clear A22;
clear q22;

FT -= AT1 * (U11 \ (L11 \ F1(p11)))(q11t) ...
	+ AT2 * (U22 \ (L22 \ F2(p22)))(q22t);

U(i2+1:end) = pcg(@(x) ...
	ATT * x ...
	- AT1 * (U11 \ (L11 \ (A1T * x)(p11)))(q11t) ...
	- AT2 * (U22 \ (L22 \ (A2T * x)(p22)))(q22t), FT, 1.e-9, 200);

clear ATT;
clear FT;
{% endhighlight %}

Solution on the interface:

![Alt Solution interface]({{ site.url }}/{{ page.assets }}/poisson2D.solInterface.png)

### Conclusion

I think the method in this post is along the lines of what was described in the wikipedia article:


> The important thing to note is that the computation of any quantities involving `\(A_{11}^{-1}\)` or `\(A_{22}^{-1}\)` involves solving 
> decoupled Dirichlet problems on each domain, and these can be done in parallel. Consequently, we need not store the Schur complement matrix 
> explicitly; it is sufficient to know how to multiply a vector by it.


Note that it is not the only method.

## Solving on the subdomains

Once we computed `\(U_{T}\)`, the solution on the interface, we can quickly solve in parallel for `\(U_{1}\)` and `\(U_{2}\)` as we already computed the LU decomposition of both `\(A_{11}\)` and `\(A_{22}\)`.

### Subdomain 1
First line of `\(\eqref{eq:mainsystem}\)` gives us
$$A_{11} U_{1} = F_{1} - A_{1T} U_{T}$$

{% highlight octave linenos %}
U(1:i1) = (U11 \ (L11 \ (F1 - A1T * U(i2+1:end))(p11)))(q11t);
{% endhighlight %}

![Alt Solution domain 1]({{ site.url }}/{{ page.assets }}/poisson2D.sol1.png)

### Subdomain 2
Similarly, second line of `\(\eqref{eq:mainsystem}\)` gives us
$$A_{22} U_{2} = F_{2} - A_{2T} U_{T}$$

{% highlight octave linenos %}
U(i1+1:i2) = (U22 \ (L22 \ (F2 - A2T * U(i2+1:end))(p22)))(q22t);
{% endhighlight %}

![Alt Solution domain 2]({{ site.url }}/{{ page.assets }}/poisson2D.sol2.png)

Ultimately we need to revert ordering of the DoFs to display the solution:
{% highlight octave linenos %}
% Reorder back DoFs
p(q) = 1:length(q);
U = U(p);
{% endhighlight %}

### Fusiiiiiiiiiiiiion !

![Alt Fusion]({{ site.url }}/{{ page.assets }}/fusion.png)

### Boum \o/

![Alt Solution]({{ site.url }}/{{ page.assets }}/poisson2D.sol.filled.png)

## Conclusion
Hopefully, you now have a good overview of how one can implement the Schur complement method.
I also hope that it wasn't too boring and that you're going to reuse that knowledge for your own projects.
I didn't talk about implementing the method for more than 2 subdomains, preconditioning, parallel implementation, &c.
Another time.

### In the next episode
In the next episode, "The Schur Complement Method: Part 2", I'm going to describe ... the _DUAL_ Schur complement method.
It's going to be a little more involved but not that much.










