#! /usr/bin/perl
use strict;
use warnings;

if (@ARGV != 2) {
	print("Usage: " . 	(split /\//, $0)[-1] . " domain.metis_graph domain.epart.3\n");
	exit 1
}

open my $FDTRIANGLES, $ARGV[0];
open my $FDIDS, $ARGV[1];

my $nbDomain = 0;
my @vertexToDomain;
my %interface;

while (<$FDTRIANGLES>) {
	# Match triangles
	if (/(?<vid1>\d+) (?<vid2>\d+) (?<vid3>\d+)/) {
		my $domain = <$FDIDS> + 1;

		$nbDomain = $domain if ($domain > $nbDomain);

		CheckVertex($+{vid1}, $domain);
		CheckVertex($+{vid2}, $domain);
		CheckVertex($+{vid3}, $domain);
	}
}

my @DOMAINFD;
foreach (1..$nbDomain) {
	open $DOMAINFD[$_], "> domain" . $_ . ".vids";
}
open my $INTERFACEFD, "> interface.vids";
while (my ($vid, $domain) = each @vertexToDomain) {
	next if !$vid;
	my $fd = exists $interface{$vid} ? $INTERFACEFD : $DOMAINFD[$domain];
	print $fd $vid."\n";
}

# Count the number of domains a vertex belongs to, if it's greater than 1
# that vertex is on the interface
sub CheckVertex {
	my ($vid, $domain) = @_;
	if ($vertexToDomain[$vid] && $vertexToDomain[$vid] != $domain) {
		$interface{$vid}++;
	} else {
		$vertexToDomain[$vid] = $domain;
	}
}