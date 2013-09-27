#!/usr/bin/perl
use strict; use warnings; use fatal;
use Getopt::Std;
use vars qw($opt_f $opt_s);
getopts('f:s:');

my $FILENAME = "generated";

die "
usage: $0 [options] <hmm file> <emission file> <length>
options:
  -f  <filename> [$FILENAME]
  -s  <int> static state mode
" unless @ARGV == 3;

my $STATIC_MODE = $opt_s;
$FILENAME = $opt_f if $opt_f;

my ($HMMFILE, $EMITFILE, $LENGTH) = @ARGV;

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

if ($STATIC_MODE) {
	
	if ($STATIC_MODE == 1) {
		my @state = qw(1x 2x 3x 4x);
		my $offset = 0;
		for (my $i = 0; $i < 1000; $i++) {
			for (my $s = 0; $s < @state; $s++) {
				generate(\%cdf, $state[$s], $LENGTH, $offset, $fasta, $tsv);
				$offset += $LENGTH;
			}
		}
	}
	elsif ($STATIC_MODE == 2) {
		my @state = qw(1x 1x-S 1x 2x 2x-S 2x 3x 3x-S 3x 4x 4x-S 4x);
		my $SLEN = 10;
		my $offset = 0;
		for (my $i = 0; $i < 1000; $i++) {
			for (my $s = 0; $s < @state; $s++) {
				if ($state[$s] =~ /S$/) {
					generate(\%cdf, $state[$s], $SLEN, $offset, $fasta, $tsv);
					$offset += $SLEN;
				} else {
					generate(\%cdf, $state[$s], $LENGTH, $offset, $fasta, $tsv);
					$offset += $LENGTH;
				}
			}
		}
	}
	elsif ($STATIC_MODE == 3) {
		my @state = qw(1x 1x-Z 1x 2x 2x-Z 2x 3x 3x-Z 3x 4x 4x-Z 4x);
		my $ZLEN = 10;
		my $offset = 0;
		for (my $i = 0; $i < 1000; $i++) {
			for (my $s = 0; $s < @state; $s++) {
				if ($state[$s] =~ /Z$/) {
					generate(\%cdf, $state[$s], $ZLEN, $offset, $fasta, $tsv);
					$offset += $ZLEN;
				} else {
					generate(\%cdf, $state[$s], $LENGTH, $offset, $fasta, $tsv);
					$offset += $LENGTH;
				}
			}
		}
	}
	elsif ($STATIC_MODE == 4) {
		my @state = qw(1x 1x-Z 1x 1x-S 1x 2x 2x-Z 2x 2x-S 2x 3x 3x-Z 3x 3x-S 3x 4x 4x-Z 4x 4x-S 4x);
		my $ZLEN = 10;
		my $SLEN = 10;
		my $offset = 0;
		for (my $i = 0; $i < 1000; $i++) {
			for (my $s = 0; $s < @state; $s++) {
				if ($state[$s] =~ /Z$/) {
					generate(\%cdf, $state[$s], $ZLEN, $offset, $fasta, $tsv);
					$offset += $ZLEN;
				} elsif ($state[$s] =~ /S$/) {
					generate(\%cdf, $state[$s], $SLEN, $offset, $fasta, $tsv);
					$offset += $SLEN;
				} else {
					generate(\%cdf, $state[$s], $LENGTH, $offset, $fasta, $tsv);
					$offset += $LENGTH;
				}
			}
		}
	} else {
		die;
	}

} else {
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
		print $fasta $emit, ",";
	
		# choose a new state
		$r = rand;
		foreach my $s (keys %{$state{$csn}{cdf}}) {
			if ($r < $state{$csn}{cdf}{$s}) {
				if ($s ne $csn) {
				 	printf $tsv "%s %d\t%d\n", $csn, $csb+1, $i+1;
				 	$csb = $i + 1;
				 	$csn = $s;
				 	print $fasta "\n";
				}
				last;
			}
		}
	}
	printf $tsv "%s %d\t%d\n", $csn, $csb+1, $LENGTH;
}

sub generate {
	my ($cdf, $model, $len, $offset, $fasta, $tsv) = @_;
	
	# write tsv
	print $tsv join("\t", $offset + 1, $offset + $len, $model), "\n";
	
	# write fasta
	my @seq;
	if    ($model =~ /S$/) {$model = 'Spike'}
	elsif ($model =~ /Z$/) {$model = 'Zero'}
	my $m = $cdf->{$model};
	for (my $i = 0; $i < $len; $i++) {
		my $r = rand;
		for (my $j = 0; $j < @$m; $j++) {
			if ($r < $m->[$j]) {
				push @seq, $j;
				last;
			}
		}
	}
	print $fasta join(",", @seq), "\n";
}


