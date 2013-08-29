#!/usr/bin/perl
use strict; use warnings; use fatal;
use DataBrowser;

die "
usage: $0 <hmm file> <emission file> <length> <filename>\n"
	unless @ARGV == 4;
my ($HMMFILE, $EMITFILE, $LENGTH, $FILENAME) = @ARGV;

# read HMM file
open(my $hfh, $HMMFILE) or die;
my %state;
while (<$hfh>) {
	next unless /^STATE:/;
	my $line = <$hfh>;
	my ($name) = $line =~ /NAME:\s+(\S+)/;
	next if $name eq 'INIT';
	while (<$hfh>) {last if /^TRANSITION/}
	my %trans;
	while (<$hfh>) {
		last if /^END/;
		my ($s2, $prob) = split;
		chop $s2;
		$trans{$s2} = $prob;
	}
	my $garbage = <$hfh>;
	my $eline = <$hfh>;
	my ($pdf) = $eline =~ /PDF:\s+(\S+)/;
	
	$state{$name} = {pdf => $pdf, trans => \%trans};
}
close $hfh;

# create cumulative distribution function for all transitions
# also make sure they sum to 1.0 in the process
foreach my $csn (keys %state) {
	my $sum = 0;
	foreach my $s (keys %{$state{$csn}{trans}}) {
		$sum += $state{$csn}{trans}{$s};
	}
	foreach my $s (keys %{$state{$csn}{trans}}) {
		$state{$csn}{trans}{$s} /= $sum;
	}

	my %cdf;
	$sum = 0;
	foreach my $s (keys %{$state{$csn}{trans}}) {
		$sum += $state{$csn}{trans}{$s};
		$cdf{$s} = $sum;
	}
	$state{$csn}{cdf} = \%cdf;
}

# read emission file
open(my $efh, $EMITFILE) or die;
my $header = <$efh>;
chomp $header;
my @model = split(/\s/, $header);
my %emit;
while (<$efh>) {
	chomp;
	my @f = split;
	for (my $i = 0; $i < @f; $i++) {
		push @{$emit{$model[$i]}}, $f[$i];
	}
}
close $efh;

# create cumulative distribution functions
my %cdf;
foreach my $label (keys %emit) {
	my $sum = 0;
	for (my $i = 0; $i < @{$emit{$label}}; $i++) {
		$sum += $emit{$label}[$i];
		$cdf{$label}[$i] = $sum;
	}
}

# emit sequence (always starts in 1x state)
open(my $fasta, ">$FILENAME.fa") or die;
open(my $tsv, ">$FILENAME.tsv") or die;
print $fasta ">$FILENAME\n";
print $tsv "#$FILENAME\n";

my $csn = "1x"; # current state name
my $csb = 0;    # current state begin
for (my $i = 0; $i < $LENGTH; $i++) {
	
	# emit a symbol from the state
	my $r = rand;
	my $emit;
	for (my $i = 0; $i < @{$cdf{$state{$csn}{pdf}}}; $i++) {
		if ($r < $cdf{$state{$csn}{pdf}}[$i]) {
			$emit = $i;
			last;
		}
	}
	print $fasta $emit;
	if ($i == $LENGTH -1) {print $fasta "\n"}
	else                  {print $fasta ","}
	
	# choose a new state
	$r = rand;
	foreach my $s (keys %{$state{$csn}{cdf}}) {
		if ($r < $state{$csn}{cdf}{$s}) {
			if ($s ne $csn) {
			 	printf $tsv "%s %d\t%d\n", $csn, $csb+1, $i+1;
			 	$csb = $i + 1;
			 	$csn = $s;
			}
			last;
		}
	}
}





