#! /usr/bin/perl
use warnings;
use strict;

use Geo::ShapeFile;
use Math::Trig 'pi';
use List::Util 'min' ,'max';

my $shapefile = new Geo::ShapeFile("depth-area-polygons-hydro-14k-122k");
$shapefile->build_spatial_index();

my $rtree = $shapefile->get_spatial_index();

my $d = 0.00955042966330666;
my $xmin = 250.816739768889;
my $ymin = 162.523014276134;
my $earthRadius = 6378137;

my @vertices;

open my $mesh, "< ../lyttelton.mesh";
while (<$mesh>) {
	last if /Vertices/;
}
my $nbVertex = <$mesh>;
foreach my $vid (1..$nbVertex) {
	<$mesh> =~ /^([-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?)\s+([-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?)\s+\d+$/;

	my $x = $d * $1 + $xmin;
	my $y = $d * $2 + $ymin;
 		
	# From Google mercator convertion
    $x = 2 * pi * $earthRadius * ($x / 256 - 0.5);
    $y = -2 * pi * $earthRadius * ($y / 256 - 0.5);

    my @shapes;
    $rtree->query_point(($x, $y),\@shapes);

	my $z = 12;
	foreach my $shape (@shapes) {
		if ($shape->contains_point(Geo::ShapeFile::Point->new(X => $x, Y => $y))) {
			my %db = $shapefile->get_dbf_record($shape->shape_id());
			$z = $db{'DRVAL2'};
			last;
		}
	}
	push @vertices, [$x, $y, $z];
}

while (<$mesh>) {
	last if /Triangles/;
}
my $nbTriangle = <$mesh>;
foreach my $vid (1..$nbTriangle) {
	<$mesh> =~ /^(\d+)\s+(\d+)\s+(\d+)\s+\d+$/;

	my $v1 = $1 - 1;
	my $v2 = $2 - 1;
	my $v3 = $3 - 1;

	print min($vertices[$v1]->[2], $vertices[$v2]->[2], $vertices[$v3]->[2]) ."\n";
}
