#! /usr/bin/perl
use strict;
use warnings;

open my $INTERFACEFD, "interface.vids";
open my $MESHFD, "poisson2D.mesh";

my %interface;
my $i = 1;
while (<$INTERFACEFD>) {
	$interface{int($_)} = $i++;
}

while (<$MESHFD>) {
	print;
	last if /Vertices/;
}

$_ = <$MESHFD>;
$_ =~ /(?<nbVertices>\d+)/;
my $nbVertices = $+{nbVertices};

print $nbVertices + keys(%interface) . "\n";

my @duplicated;
foreach my $vid (1..$nbVertices) {
	$_ = <$MESHFD>;
	print;
	push @duplicated, $_ if exists $interface{$vid};
}
print foreach @duplicated;

while (<$MESHFD>) {
	print;
	last if /Triangles/;
}

$_ = <$MESHFD>;
$_ =~ /(?<nbTriangles>\d+)/;
my $nbTriangles = $+{nbTriangles};

print $nbTriangles, "\n";

foreach my $tid (1..$nbTriangles) {
	$_ = <$MESHFD>;
	/(?<vid1>[\d\.]+) (?<vid2>[\d\.]+) (?<vid3>\d+) (?<domain>\d+)/;
	if ($+{domain} == 1) {
		my ($vid1, $vid2, $vid3) = ($+{vid1}, $+{vid2}, $+{vid3});
		$vid1 = $nbVertices + $interface{$vid1} if exists $interface{$vid1};
		$vid2 = $nbVertices + $interface{$vid2} if exists $interface{$vid2};
		$vid3 = $nbVertices + $interface{$vid3} if exists $interface{$vid3};
		print "$vid1 $vid2 $vid3 $+{domain}\n";
	} else {
		print;
	}
}

print while <$MESHFD>;
