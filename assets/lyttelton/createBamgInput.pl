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
