#! /usr/bin/perl
use strict;
use warnings;
use Data::Dumper;

if (@ARGV != 2) {
	print("Usage: " . 	(split /\//, $0)[-1] . " domain.metis_graph domain.metis_graph.epart.3\n");
	exit 1
}

my %interfaceGraph;
my %interface;

my %triangleToVertices;
my %vertexToTriangles;

open my $INTERFACEFD, "interface.vids";
while (<$INTERFACEFD>) {
	$interface{int($_)}++;
}

$ARGV[1] =~ /.epart.(\d+)$/ or die "$ARGV[1] filename should end with the number of subdomains";
my $nbDomain = $1;

foreach my $did (2..$nbDomain) {

	open my $FDTRIANGLES, $ARGV[0];
	open my $FDIDS, $ARGV[1];

	%triangleToVertices = ();
	%vertexToTriangles = ();

	# Build connectivity and reverse connectivity tables for 
	# all but one doamins elements on the interface
	my $triangleId = 0;
	while (<$FDTRIANGLES>) {
		if (/(?<vid1>\d+) (?<vid2>\d+) (?<vid3>\d+)/) {
			++$triangleId;
			my $domain = <$FDIDS> + 1;
			if ($domain == $did) 
			{
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
			if (IsBoundaryEdge($vid1, $vid2)) {
				AddInterfaceEdge($vid1, $vid2);
			}
		}
	}

}

# Walk the interface
WalkInterface();

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

# A boundary edge has its two vertices on the interface 
# and belongs to only one domain 1 triangle.
sub IsBoundaryEdge {
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

sub WalkInterface {
	while ((my $start = FindInterfaceEnd()) != -1) {
		WalkInterfaceI($start, $start);
	}
}

# The interface should have two ends.
# An end is simply a vertex connected to only one other vertex.
#=============> NO LONGER TRUE, check medit mesh file edge list <==============
sub FindInterfaceEnd {
	while (my ($from, $tos) = each %interfaceGraph) {
		if (keys %{$tos} == 1) {
			return $from;
		}
	}
	return -1;
}

# Perform a simple DFS to walk the interface 
# and print the vertex ids as we go.
sub WalkInterfaceI {
	my ($u, $v) = @_;
	print $v."\n";
	delete $interfaceGraph{$u}{$v};
	delete $interfaceGraph{$v}{$u};
	for (keys %{$interfaceGraph{$v}}) {
		WalkInterfaceI($v, $_);
	}
	delete $interfaceGraph{$v};
}

