---
layout: post
title: "The Lyttelton port wave penetration project: Part 1"
description: ""
category: "Numerical methods"
tags: [mathematics, parallel computing, numerical hydraulics]
assets: assets/lyttelton
---
{% include JB/setup %}

In this series of mini-posts, I will try to create a wave penetration simulator, a baby version of [HARES](http://www.svasek.com/modelling/hares.html).

For this kickoff post I will focus on retrieving the [Lyttelton port](http://www.lpc.co.nz) geometry, generate a mesh, partition it into more than two subdomains (5 for example)
and modify the Perl scripts we used in [the first Schur complement post]({% post_url 2013-10-22-the-schur-complement-method %}) in order to deal with many subdomains. I will also demonstrate how we can use the 
[Google Maps JavaScript API](https://developers.google.com/maps/documentation/javascript/tutorial), the [OpenStreetMap](http://www.openstreetmap.org/#map=16/-43.6075/172.7155) tiles as well as [WebGL](https://developer.mozilla.org/en-US/docs/Web/WebGL) to display the mesh and the sudomains.

<!-- more -->

## Retrieving the port geometry

### Raw data

I will retrieve the Lyttelton port geometry from the OSM data using the [Overpass API](http://wiki.openstreetmap.org/wiki/Overpass_API) and taking advantage 
of the [coastline tag](http://wiki.openstreetmap.org/wiki/Coastline).

Here is the query I'm using the get the coastline nearby the Lyttelton port:
{% highlight xml linenos %}
<osm-script output="xml">
  <query type="way">
   <has-kv k="natural" v="coastline"/>
   <bbox-query
    w="172.70914950376547" s="-43.611191483302505"
    e="172.72557535177208" n="-43.60308883447627"/>
  </query>
 <print mode="meta"/>
 <recurse type="down"/>
 <print mode="meta" order="quadtile"/>
</osm-script>
{% endhighlight %}

Let's run it in an Overpass sandbox:
<iframe width="850" height="350" frameborder="0" scrolling="yes" marginheight="0" marginwidth="0" src="http://overpass-turbo.eu/s/2O4" style="border: 1px solid black">unwantedtext</iframe><br/><small><a href="http://overpass-turbo.eu/s/2O4">View Larger Map</a></small>

As you can see the query is returning too much coastline.

Here is the exported data in the OSM XML file format:

* [Raw Lyttleton port coastline]({{ site.url }}/{{ page.assets }}/lyttleton_raw.osm)

### Open boundary

The next stage is to cleanup unwanted edges and close the geometry with an open boundary.
An open boundary is not a physical boundary. The condition on the open boundary will be discused at a later stage.

In order to do so, I'm using the [josm](https://josm.openstreetmap.de) editor.

<div class="galleria">
   <a href="{{ site.url }}/{{ page.assets }}/lyttelton_raw.png">
      <img src="{{ site.url }}/{{ page.assets }}/lyttelton_raw.png" 
      data-title="Raw data loaded into josm" 
      data-description=""
      data-big="{{ site.url }}/{{ page.assets }}/lyttelton_raw.png"/>
   </a>
   <a href="{{ site.url }}/{{ page.assets }}/lyttelton_with_open_boundary.png">
      <img src="{{ site.url }}/{{ page.assets }}/lyttelton_with_open_boundary.png" 
      data-title="The final geometry with an open boundary" 
      data-description=""
      data-big="{{ site.url }}/{{ page.assets }}/lyttelton_with_open_boundary.png"/>
   </a>
</div>

Here is the final geometry in the OSM XML file format:

* [Lyttleton port with open boundary]({{ site.url }}/{{ page.assets }}/lyttleton.osm)

It's worth noticying that manually added nodes are flaged with an `action='modify'` tag by the editor:
{% highlight xml %}
<node id='-2498' action='modify' visible='true' lat='-43.613555333973856' lon='172.71590659752025' />
{% endhighlight %}
I will take advantage of this in order to distinguish between closed and open boundaries later on.

### Conclusion
I'm done with this part. Note that the current geometry is currently in "latitude/longitude" space. Next stage is to generate a mesh.

## Mesh generation

In this part I'm going to use [bamg](http://www.ann.jussieu.fr/hecht/ftp/bamg/bamg.pdf), which is shiped with the Metis software I demonstrated in the first post.
Given a geometry, bamg will generate a mesh using sophisticated algorithm.

I must convert the OSM file into the bamg mesh file format.

One thing we have to take care of is to project the geometry from "lat/lon" space to [Mercator](http://en.wikipedia.org/wiki/Mercator_projection) space and then make the projected geometry dimensionless.

Here is a script doing that:
{% highlight perl linenos %}
#! /usr/bin/perl
use strict;
use warnings;
use XML::XPath;
use Math::Trig 'pi';
use List::Util 'min' ,'max';

use Geo::Proj4;
my $proj = Geo::Proj4->new(init => "epsg:3857") or die;

my $osm = XML::XPath->new(filename => 'lyttelton.osm');
my $nodeset = $osm->find('/osm/node');

my $bamgid = 0;
my %nodes;

my $xmin = 'inf';
my $xmax = '-inf';
my $ymin = 'inf';
my $ymax = '-inf';

my $earthRadius = 6378137;

foreach my $node ($nodeset->get_nodelist) {
  ++$bamgid;
    my $osmid = $node->getAttribute("id");
    my $lat = $node->getAttribute("lat");
    my $lon = $node->getAttribute("lon");
    my $vnode = $node->getAttribute("action");

    # EPSG 3857 projection
    my ($x,$y) = $proj->forward($lat, $lon);

    # Google mercator convertion
    $x = 256 * (0.5 + $x / (2 * pi * $earthRadius));
    $y = 256 * (0.5 - $y / (2 * pi * $earthRadius));

    $xmin = min($x,$xmin);
    $xmax = max($x,$xmax);
        
    $ymin = min($y,$ymin);
    $ymax = max($y,$ymax);

    $nodes{$osmid} = [$bamgid, $x, $y, $vnode];
}

# Make the geometry dimensionless
my $d = max($xmax - $xmin, $ymax - $ymin);
foreach my $node (values(%nodes)) {
    $node->[1] = ($node->[1] - $xmin) / $d;
    $node->[2] = ($node->[2] - $ymin) / $d;
}

# Remember these in order to display on top of a Google map later on
print "Translate: ($xmin, $ymin)\n";
print "Scale:$d\n";

# Create bamg input file
open my $out, "> lyttelton_0.mesh";
print $out "MeshVersionFormatted 1\n\n";
print $out "Dimension 2\n\n";
print $out "Vertices " . $bamgid . "\n";

foreach my $node (sort {$a->[0] <=> $b->[0]} values(%nodes)) {
  printf $out  "%.12f %.12f 1\n", $node->[1], $node->[2];
}

my $edgeset = $osm->find('/osm/way/nd');

print $out "\nEdges " . ($edgeset->size - 1) . "\n";

my $start = 0;

foreach my $edge ($edgeset->get_nodelist) {
    my $end = $edge->getAttribute("ref");
    if ($start != 0) {
      my $boundaryId = $nodes{$start}->[3] || $nodes{$end}->[3] ? 2 : 1;
      print $out $nodes{$start}->[0] ." " . $nodes{$end}->[0] . " " . $boundaryId . "\n";
    }
    $start = $end;
}

{% endhighlight %}

Here is the bamg input:

* [Lyttelton geometry for bamg]({{ site.url }}/{{ page.assets }}/lyttelton_0.mesh)

Here is the command line for bamg:

{% highlight bash%}
bamg  -g lyttelton_0.mesh -coef 0.06  -ratio 20 -o lyttelton.mesh
{% endhighlight %}

And here is the generated mesh file:

* [Lyttelton mesh]({{ site.url }}/{{ page.assets }}/lyttelton.mesh)

This is what a Lyttelton port small mesh looks like:
![yttelton mesh pic]({{ site.url }}/{{ page.assets }}/lyttelton_mesh.png)

### Conclusion

It seems quite easy to use the OSM data for baby coastal engineering. 
I'm missing the water depth inside the port though: I will address this in another post (3rd one).
In the next stage we will decompose the generated mesh into 5 subdomains.

## Partitioning

To make this project more challenging, I will implement a parallel solver.
In order to do so I will apply domain decomposition techniques.
I detailed the partioning tools and the methodology in [the first Schur complement post]({% post_url 2013-10-22-the-schur-complement-method %}) 
so that this part is going to be concise.

* Metis input graph : [lyttelton.metis_graph]({{ site.url }}/{{ page.assets }}/lyttelton.metis_graph)
* Metis output element partition : [lyttelton.metis_graph.epart.5]({{ site.url }}/{{ page.assets }}/lyttelton.metis_graph.epart.5)
* Metis output vertex partition : [lyttelton.metis_graph.npart.5]({{ site.url }}/{{ page.assets }}/lyttelton.metis_graph.npart.5)

That's it.

### Interface computation and unknown reordering

One last thing we have to do is to compute the domains interface and reorder the interface unknowns to improve the sparsity pattern of the "interface-interface" stiffness matrix.
Again I addressed this in the first Schur complement post, but only for the two subdomains case. I had to slightly modify the scripts in that post to make them work for many
subdomains.

Here are the modifications on github:

* [Find the interface for many subdomains](https://github.com/ssrb/ssrb.github.com/commit/d780cf949630b57695e050f1359d9434cc268803#diff-1ea8218905c12a730bfa1a0e6f720716)
* [Reorder the unknowns for many subdomains](https://github.com/ssrb/ssrb.github.com/commit/d780cf949630b57695e050f1359d9434cc268803#diff-74e54e3a95ea9e1eb5241ed29445b3df)

Ultimately here is the interface:

* [interface2.vids]({{ site.url }}/{{ page.assets }}/interface2.vids])

## Visualization

So here we are ! Ready to display our domain and subdomains on top of the matching OSM tiles thanks to the Google Maps JavaScript API and WebGL:

<iframe width="850" height="700" frameborder="0" scrolling="no" marginheight="0" marginwidth="0" src="{{ site.url }}/{{ page.assets }}/map.html" style="border: 1px solid black">unwantedtext</iframe><br/><small><a href="{{ site.url }}/{{ page.assets }}/map.html">View Larger Map</a></small>

## Conclusion
In the next mini-post, I will implement several parallel versions of the "trace first" Schur complement method applied to the data I generated using [MPI](http://en.wikipedia.org/wiki/Message_Passing_Interface).

<script type="text/javascript" src="{{ site.url }}/rungalleria.js"></script>