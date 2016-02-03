---
layout: post
title: "The alternator simulation project"
description: ""
category: "Numerical methods"
tags: [mathematics, electrotechnics]
assets: assets/alternator

---
{% include JB/setup %}

In this post I'm going to illustrate how to implement a finite element simulation of an [alternator](https://en.wikipedia.org/wiki/Induction_generator) in motion coupled to a linear electrical circuit.

<!-- more -->

The simulation computes the magnetic vector potential within a (2D) radial cross section of the machine I'm about to describe.

Here's the result:

<iframe width="500" height="500" frameborder="0" scrolling="no" marginheight="0" marginwidth="0" src="http://ssrb.github.io/alternator-fem-webapp/alternator.html" style="border: 1px solid black">unwantedtext</iframe><br/><small><a href="http://ssrb.github.io/alternator-fem-webapp/alternator.html">View Larger Simulation</a></small>

Code is [here](https://github.com/ssrb/alternator-fem-webapp)

## Description of the machine

The machine we're going to study is an unusual induction generator of the [squirel-cage](https://en.wikipedia.org/wiki/Squirrel-cage_rotor) type.

This is `\(1/12^{th}\)` of the machine:
![machine]({{ site.url }}/{{ page.assets }}/machine.png)

The rotor has 12 massive coper rods shorted by rings at both ends: the so called squirel cage.
Rods are enclosed in a cyclinder made of iron laminations.

The stator, made of iron laminations too, has 12 winding slots. The winding pattern is as follow: 

![winding1]({{ site.url }}/{{ page.assets }}/winding-1.png)

Notice that the field coil goes up/up, down/down, up/up ...
This will be very important later on when computing the [current density](https://en.wikipedia.org/wiki/Current_density) in the coils.

![winding2]({{ site.url }}/{{ page.assets }}/winding-2.png)

Again, notice that the armature coil goes up/down when the field coil goes up/up, and then down/up, up/down ...


The diameter of the machine is 80mm, its depth/thickness is 50mm.

## Modelization

### Periodicity

Because of the problem periodicity, we can reduce the ammount of computation by modeling only a slice of the machine.

To find the size of the smallest possible slice we compute the [greatest common divisor](https://en.wikipedia.org/wiki/Greatest_common_divisor) (`\(gcd\)`) of the number of rods and slots and only simulate a `\(360 / gcd\)` degrees machine slice.

Next, in order to decide if we need to use [periodic or anti-periodic boundary conditions](https://en.wikipedia.org/wiki/Periodic_boundary_conditions) to join our alternator slices, we simply observe that if we apply [Kirchof current's law](https://en.wikipedia.org/wiki/Kirchhoff%27s_circuit_laws) to the rotor squirel cage topology, the [eddy currents](https://en.wikipedia.org/wiki/Eddy_current) flowing through one rod will flow in the opposite direction in the previous and next rods so that we use periodic (respectively anti-periodic) boundary condition if the slice contains an even (respectively odd) number of rod.

It comes that for this particular endeavour we can simulate `\(1/12^{th}\)` of the machine (30 degrees), so a single rod and a single slot, and must use _anti-periodic_ boundary condition.

Whem implementing the finite element method using first order Lagrange element to solve non periodic problem we have a one to one mapping between mesh vertices and [degrees of freedom](https://en.wikipedia.org/wiki/Degrees_of_freedom_(physics_and_chemistry))(`\(dof\)`) so that we can identify vertices with `\(dof\)`.

For example a single quadrangular first order Lagrange element has 4 vertices and leads to a 4 `\(dof\)` [stiffness matrix](https://en.wikipedia.org/wiki/Stiffness_matrix):

![element]({{ site.url }}/{{ page.assets }}/element-1.png)

However if we consider an (anti-)periodic along one axis element, it leads to a 2 `\(dof\)` stiffness matrix:

![periodic element]({{ site.url }}/{{ page.assets }}/element-2.png)

We can also illustrate further what's going on by writing the equations:

$$\begin{aligned} 
    & 
    \begin{pmatrix}
    a_{0,0} & a_{0,1} & a_{0,2} & a_{0,3}\\\
    a_{1,0} & a_{1,1} & a_{1,2} & a_{1,3}\\\
    a_{2,0} & a_{2,1} & a_{2,2} & a_{2,3}\\\
    a_{3,0} & a_{3,1} & a_{3,2} & a_{3,3}
    \end{pmatrix}
    \begin{pmatrix}
    x_0\\\
    x_1\\\
    x_2\\\
  	x_3
    \end{pmatrix}
    =
    \begin{pmatrix}
    b_0 \\\
    b_1 \\\ 
    b_2 \\\
    b_3
    \end{pmatrix}\\\
	\Leftrightarrow & 
    \begin{pmatrix}
    a_{0,0} & a_{0,1} & a_{0,2} & a_{0,3}\\\
    a_{0,1} & a_{1,1} & a_{1,2} & a_{1,3}\\\
    a_{0,2} & a_{1,2} & a_{2,2} & a_{2,3}\\\
    a_{0,3} & a_{1,3} & a_{2,3} & a_{3,3}
    \end{pmatrix}
    \begin{pmatrix}
    x_0\\\
    \pm x_0\\\
    x_2\\\
  	\pm x_2
    \end{pmatrix}
    =
    \begin{pmatrix}
    b_0 \\\
    b_1 \\\ 
    b_2 \\\
    b_3
    \end{pmatrix}, x_1 = \pm x_0, x_3 = \pm x_2 \text{, as per stiffness matrix symmetry and element periodicity}\\\
	\Leftrightarrow & 
    \begin{pmatrix}
    a_{0,0} \pm a_{0,1} & a_{0,2} \pm a_{0,3}\\\
    a_{0,1} \pm a_{1,1} & a_{1,2} \pm a_{1,3}\\\
    a_{0,2} \pm a_{1,2} & a_{2,2} \pm a_{2,3}\\\
    a_{0,3} \pm a_{1,3} & a_{2,3} \pm a_{3,3}
    \end{pmatrix}
    \begin{pmatrix}
    x_0\\\
    x_2
    \end{pmatrix}
    =
    \begin{pmatrix}
    b_0 \\\
    b_1 \\\ 
    b_2 \\\
    b_3
    \end{pmatrix}, x_1 = \pm x_0, x_3 = \pm x_2\\\
    \Leftrightarrow & 
    \begin{pmatrix}
    a_{0,0} \pm 2 a_{0,1} + a_{1,1} & a_{0,2} \pm a_{0,3} \pm a_{1,2} + a_{1,3}\\\
    a_{0,2} \pm a_{0,3} \pm a_{1,2} + a_{1,3} & a_{2,2} \pm 2 a_{2,3} + a_{3,3}
    \end{pmatrix}
    \begin{pmatrix}
    x_0\\\
    x_2
    \end{pmatrix}
    =
    \begin{pmatrix}
    b_0 \pm b_1\\\
    b_2 \pm b_3
    \end{pmatrix}, x_1 = \pm x_0, x_3 = \pm x_2
    \end{aligned}$$

Notice that we also end up with a symmetric linear system.

So that the key to implement (anti-)peridic boundary conditions is to compute the [non injective](https://en.wikipedia.org/wiki/Injective_function) mapping from mesh vertices to `\(dof\)` as well as introduce a sign/polarity vector in order to remember the sign of each nodal value relative to its `\(dof\)`
(or in short: to deal with all these `\(\pm\)` in the equations).

* This is how I compute the mapping:

{% gist ssrb/44a84951fe6dca229bf3 %}

Relevant project class is in <https://github.com/ssrb/alternator-fem-webapp/blob/master/domain.ts>

* And this is how I assemble stiffness matrix and load vector using the computed mapping:

{% gist ssrb/6ff15cec796b1e1a0116 %}

Relevant project class is in <https://github.com/ssrb/alternator-fem-webapp/blob/master/solver.ts>

### Motion

In order to simulate the entire machine, we need to match rotor and stator `\(dof\)` along their common interface.
One way to get a match is to first discretize the interface and then use these interface vertices to compute the [Delaunay triangulations](https://en. wikipedia.org/wiki/Delaunay_triangulation) of both the stator and rotor domains.

![mesh]({{ site.url }}/{{ page.assets }}/alter012.000.png)

How about motion ? Idealy we would like to solve for *any* rotor/stator angle.
That's doable but a bit too ambitious and we're instead going to only consider rotations multiple of a fixed "rotation step".
Furthermore we choose a "rotation step" equal to a unit fraction of the domain slice angle, for example `\(1/32\)` of 30 degrees, so that, by uniformly dividing the interface into 32 segments we can match rotor and stator `\(dof\)` regardless of the "discrete" rotation angle.

The advantage of this method is that we compute mesh for both domains once and for all, regardless of the rotation angle.

A drawback of this method is that given the angular speed of the rotor we need to choose the time step of the simulation accordingly.
In our example, for a rotor speed of 3000rpm, the time step will be `\(60 / (12 * 32 * 3000)\)` seconds.

As it rotates, the rotor slice will face both a stator and an "anti-stator" slice. 
But we won't explicitely model the anti-stator slice since its nodal values are the opposite of those of the stator slice.
Instead we will model a rotor slice which wraps around and use the sign/polarity vector introduced earlier in order to decide if a rotor nodal value is coupled to a stator or "anti-stator" `\(dof\)`.

One last thing we need to discuss is how to generate the mapping between rotor/stator "local" `\(dof\)` to the machine "global" `\(dof\)`:

1. We first map interior "local" `\(dof\)` for both domains;
2. We then match interface rotor/stator local `\(dof\)` and assign a single "global" `\(dof\)` to each pair.

Here is a sketch:
![motion]({{ site.url }}/{{ page.assets }}/element-3.png)

Here is how I compute the mapping:

{% gist ssrb/c56b14b2e89c331e9a40 %}

Relevant project class is in <https://github.com/ssrb/alternator-fem-webapp/blob/master/domain.ts>

### Maxwell's equations

In this section, I will describe as concisely as I can the derivation of the [PDE](https://en.wikipedia.org/wiki/Partial_differential_equation) we're about to solve.

DISCLAIMER: I'M BY NO MEANS A MAXWELL'S EQUATIONS EXPERT.

Within the machine, the physical quantities we're intested in, such as magnetic vector potential and current density, are governed by [Maxwell's equations](https://en.wikipedia.org/wiki/Maxwell%27s_equations).

* The main equation we're interested in is the differential form of [Ampère's circuital law](https://en.wikipedia.org/wiki/Amp%C3%A8re%27s_circuital_law#Differential_form) which
relates the magnetic vector potential `\(A\)` to the [electric field](https://en.wikipedia.org/wiki/Electric_field) `\(E\)`:

$$\nabla \times \left(\mu \nabla \times A \right) = \sigma E \$$

`\(\mu\)` is the [relative magnetic permeability](https://en.wikipedia.org/wiki/Permeability_(electromagnetism)) and `\(\sigma\)` the [electrical conductivity](https://en.wikipedia.org/wiki/Electrical_resistivity_and_conductivity). In this toy simulation, these quantities are considered constant for a particular medium (coper, iron, air).

* The next equation we use is the [Maxwell–Faraday induction law](https://en.wikipedia.org/wiki/Faraday%27s_law_of_induction#Maxwell.E2.80.93Faraday_equation).
It relates the electric field `\(E\)` to the magnetic vector potential `\(A\)` and the [electric potential](https://en.wikipedia.org/wiki/Electric_potential) `\(V\)` :

$$E = -\nabla V - \frac{\partial A}{\partial t}$$

Putting these two together we got:

$$\nabla \times \left(\mu \nabla \times A \right) + \sigma \frac{\partial A}{\partial t} = - \sigma \nabla V$$

For this 2D problem `\(A\)` and `\(\nabla V\)` are reduced to their z components. That is `\(A \equiv \left(0,0,A\right)\)` and `\(\nabla V \equiv \left(0,0,\partial V / \partial z\right) \)`.

Now we're going to make a gross simplification and say that `\(\partial V / \partial z \neq 0\)` only within the coils. 
Furthermore, within a coil, we're going to estimate this quantity as the difference of electric potential at both ends of a single conductor 
divided by its length: 

$$\frac{\partial V}{\partial z} = \frac{V_{z_b} - V_{z_a}}{L}$$

Since the machine is 50mm thick, we set `\(L = 0.05\)`.

Then, using [Ohm's law](https://en.wikipedia.org/wiki/Ohm%27s_law), we say that:

$$V_{z_a} - V_{z_b} = iR$$

with `\(i\)` the current going through the coil and `\(R\)` the resistance of that 50mm long coper conductor.
Remember that we care about the sign of `\(i\)` which depends on the winding pattern.

We can estimate `\(R\)` using [Pouillet's law](https://en.wikipedia.org/wiki/Electrical_resistivity_and_conductivity):

$$R = \frac{L}{\sigma \mathscr{A}}$$

with `\(\mathscr{A}\)` the cross-sectional area of the conductor.

Ultimately within a single conductor we got:

$$- \sigma \nabla V \equiv \sigma \frac{V_{z_a} - V_{z_b}}{L} = \sigma \frac{iR}{L} = \sigma \frac{i}{L} \cdot \frac{L}{\sigma \mathscr{A}} = \frac{i}{\mathscr{A}}$$

In the code we estimate the cross-sectional area of a single conductor dividing the cross-sectional area of the entire coil by the number of conductors within the coil.

Putting everything back together we want to solve:

$$\begin{aligned} 
\nabla \times \left(\mu \nabla \times A \right) + \sigma \frac{\partial A}{\partial t} = \sum_{k = 1}^{n_{coil}}{\psi_k i_k}\\\
\psi_k(x) = \begin{cases}
\frac{1}{\mathscr{A}_k} \text{ within coil } k\\\
0 \text{ otherwise}
\end{cases}
\end{aligned}$$

for `\(A \equiv \left(0,0,A\right)\)` and `\(i_k, k = 1 \dots n_{coil}\)`.

### Coupling between magnetic & electrical equations

At this stage we got two unknowns:

* the magnetic vector potential;
* the current flowing through the coils

We're going to model a _linear_ electrical circuit the alternator coils are connected to using a [lumped element model](https://en.wikipedia.org/wiki/Lumped_element_model).

The lumped element model of an N buses linear electrical circuit can be concisely characterised  using the [nodal admittance matrix](https://en.wikipedia.org/wiki/Nodal_admittance_matrix):

Work in progress ;-)

<https://en.wikipedia.org/wiki/Lenz%27s_law>

## Conclusion

Work in progress ;-)


