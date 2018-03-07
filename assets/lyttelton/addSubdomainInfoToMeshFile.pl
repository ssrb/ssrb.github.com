#! /usr/bin/perl

use strict;
use warnings;

if (@ARGV != 2) {
	print("Usage: " . 	(split /\//, $0)[-1] . " domain.mesh domain.metis_graph.epart.3\n");
	exit 1
}

open my $MESH, $ARGV[0];
open my $FDIDS, $ARGV[1];

while (<$MESH>) {
	last if/Triangles/;
	print;
}

print "Triangles\n";
my $nbTriangle = <$MESH>;
print $nbTriangle;

foreach (1..$nbTriangle) {
	my $triangle = <$MESH>;
	$triangle =~ /(\d+) (\d+) (\d+)/;
	my $domain = <$FDIDS>;
	++$domain;
	print "$1 $2 $3 $domain\n";
}

print while <$MESH>;