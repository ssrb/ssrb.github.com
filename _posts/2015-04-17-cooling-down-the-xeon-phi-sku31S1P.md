---
layout: post
title: "Cooling Down The Xeon Phi SKU 31S1P"
description: ""
category: "HPC"
assets: assets/cooldownphi
tags: [parallel computing]
---
{% include JB/setup %}

<script src="https://embed.github.com/view/3d/ssrb/ssrb.github.com/master/assets/cooldownphi/phifan.stl"></script>

<!-- more -->

A few months ago, Dr Z who's with the [CAPSL team](http://www.capsl.udel.edu/people.shtml) sent me a [Xeon Phi co-processor](http://en.wikipedia.org/wiki/Xeon_Phi):

![Phi]({{ site.url }}/{{ page.assets }}/2015-04-29 21.07.45.jpg)

This particular board, SKU 31S1P, is not usually installed inside desktop computers as it is meant to be passively cooled in HPC data center environment.

Airflow requirements for the passive cooling solutions are listed in the section 3.3.2 of the [coprocessor datasheet]("{{ site.url }}/{{ page.assets }}/Xeon-Phi-Coprocessor-Datasheet.pdf"):

> In order to ensure adequate cooling of the SE10P/7120P/3120P 300W and 31S1P
> 270W SKUs with a 45oC inlet temperature, the system must be able to provide 33 ft3/
> min of airflow to the card with 7.2 ft3/min on the secondary side and the remainder on
> the primary side. The total pressure drop (assuming a multi-card installation conforming
> to the PCI Express* mechanical specification) is 0.54 in H2O at this flow rate.

## Rapid prototyping

To cool down the Phi inside an ATX desktop case I came up with a custom 8mm fan duct adapter I designed in Blender and 3D printed using my DIY [Prusa Mendel](http://reprap.org/wiki/Prusa_Mendel):

<div class="galleria">
  <a href="{{ site.url }}/{{ page.assets }}/blenderviewport.png">
    <img src="{{ site.url }}/{{ page.assets }}/blenderviewport.png" 
    data-title="" 
    data-description=""
    data-big="{{ site.url }}/{{ page.assets }}/blenderviewport.png"/>
  </a>
  <a href="{{ site.url }}/{{ page.assets }}/2015-04-28 00.00.03.jpg">
    <img src="{{ site.url }}/{{ page.assets }}/2015-04-28 00.00.03.jpg" 
    data-title="" 
    data-description=""
    data-big="{{ site.url }}/{{ page.assets }}/2015-04-28 00.00.03.jpg"/>
  </a>
  <a href="{{ site.url }}/{{ page.assets }}/2015-03-05 00.49.48.jpg">
    <img src="{{ site.url }}/{{ page.assets }}/2015-03-05 00.49.48.jpg" 
    data-title="" 
    data-description=""
    data-big="{{ site.url }}/{{ page.assets }}/2015-03-05 00.49.48.jpg"/>
  </a>
  <a href="{{ site.url }}/{{ page.assets }}/2015-03-05 00.50.05.jpg">
    <img src="{{ site.url }}/{{ page.assets }}/2015-03-05 00.50.05.jpg" 
    data-title="" 
    data-description=""
    data-big="{{ site.url }}/{{ page.assets }}/2015-03-05 00.50.05.jpg"/>
  </a>
</div>

<video width="640" height="480" controls="controls">
<source src="https://www.dropbox.com/s/b9pl26lpjuxc63z/phifantimelapse3.mp4?raw=1" type="video/mp4" />
</video>

Here are the files if your interested in customizing or printing your own:


* [Blend file]("{{ site.url }}/{{ page.assets }}/phifan.blend")
* [STL file]("{{ site.url }}/{{ page.assets }}/phifan.stl")

## Tests

I tested the fan duct with 2 different 12VDC fans:


* a **YATE LOON D80SH-12**, 12V, 0.18A **32ft3/m**, 34 dBA ;
* a **Delta Electronics FFB0812EHE**, 12V, 1.35A, **80.16ft3/m**, 52.5 dBA

<div class="galleria">
  <a href="{{ site.url }}/{{ page.assets }}/2015-04-27 13.47.57.jpg">
    <img src="{{ site.url }}/{{ page.assets }}/2015-04-27 13.47.57.jpg" 
    data-title="YATE LOON VS Delta Electronics" 
    data-description=""
    data-big="{{ site.url }}/{{ page.assets }}/2015-04-27 13.47.57.jpg"/>
  </a>
  <a href="{{ site.url }}/{{ page.assets }}/2015-04-27 11.00.50.jpg">
    <img src="{{ site.url }}/{{ page.assets }}/2015-04-27 11.00.50.jpg" 
    data-title="Phi installed" 
    data-description=""
    data-big="{{ site.url }}/{{ page.assets }}/2015-04-27 11.00.50.jpg"/>
  </a>
</div>

To test, I'm running the [Embree](https://embree.github.io) high performance raytracing kernels examples on the Phi
and monitoring temperature VS power.

### D80SH-12

* 67&deg;C@103W
* 73&deg;C@121W
* 86&deg;C@168W

<div class="galleria">
  <a href="{{ site.url }}/{{ page.assets }}/phi1.png">
    <img src="{{ site.url }}/{{ page.assets }}/phi1.png" 
    data-title="" 
    data-description=""
    data-big="{{ site.url }}/{{ page.assets }}/phi1.png"/>
  </a>
  <a href="{{ site.url }}/{{ page.assets }}/phi2.png">
    <img src="{{ site.url }}/{{ page.assets }}/phi2.png" 
    data-title="" 
    data-description=""
    data-big="{{ site.url }}/{{ page.assets }}/phi2.png"/>
  </a>
  <a href="{{ site.url }}/{{ page.assets }}/phi3.png">
    <img src="{{ site.url }}/{{ page.assets }}/phi3.png" 
    data-title="" 
    data-description=""
    data-big="{{ site.url }}/{{ page.assets }}/phi3.png"/>
  </a>
</div>

According to Intel that's too much:

![Thermal spec]({{ site.url }}/{{ page.assets }}/thermalspec.png)

### Delta Electronics FFB0812EHE

* 39&deg;C@93W
* 42&deg;C@112W
* 48&deg;C@156W

<div class="galleria">
  <a href="{{ site.url }}/{{ page.assets }}/phi4.png">
    <img src="{{ site.url }}/{{ page.assets }}/phi4.png" 
    data-title="" 
    data-description=""
    data-big="{{ site.url }}/{{ page.assets }}/phi4.png"/>
  </a>
  <a href="{{ site.url }}/{{ page.assets }}/phi5.png">
    <img src="{{ site.url }}/{{ page.assets }}/phi5.png" 
    data-title="" 
    data-description=""
    data-big="{{ site.url }}/{{ page.assets }}/phi5.png"/>
  </a>
  <a href="{{ site.url }}/{{ page.assets }}/phi6.png">
    <img src="{{ site.url }}/{{ page.assets }}/phi6.png" 
    data-title="" 
    data-description=""
    data-big="{{ site.url }}/{{ page.assets }}/phi6.png"/>
  </a>
</div>

This is good but ...

it's extremely noisy (rated 52.5 dBA ...) so that I simply can't work near the box and must unplug/remove the Phi from the computer when I'm not using it. Not a big problem as it's just a hobby.

In the next post of this "series", I will try to implement a parallel [Lattice Boltzmann](https://en.wikipedia.org/wiki/Lattice_Boltzmann_methods) fluid simulation targetting the Phi. Stay tuned :-)

<script type="text/javascript" src="{{ site.url }}/rungalleria.js"></script>
