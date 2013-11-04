#!/usr/bin/perl
use strict;
use warnings;
use FileHandle;
use Getopt::Std;
use vars qw($opt_h $opt_q $opt_d $opt_l $opt_m $opt_t
	$opt_z $opt_Z $opt_s $opt_S $opt_w);
getopts('hq:dt:l:m:z:Z:s:S:w:');

my $WDIR   = "dosamatic_scratch";
my $QUALITY = 20;
my $TRANS  = "1e-6";
my $ZERO1  = "1e-2";
my $ZERO2  = "0.5";
my $SPIKE1 = "1e-2";
my $SPIKE2 = "0.9";
my $LOW    = 0;
my $LIMIT  = 256;
my $WSIZE  = 0;

die "
usage: hmm_gen.pl <sam file(s)>
options:
  -d <string> working directory name [$WDIR]
  -q <int>    quality value threshold [$QUALITY]
  -t <float>  CNV transition probability [$TRANS]
  -l <int>    low-count threshold [$LOW]
  -m <int>    maximum value from any emission [$LIMIT]
  -z <float>  transition probability to zero state [$ZERO1]
  -Z <float>  transition probability stay zero state [$ZERO2]
  -s <float>  transition probability to spike state [$SPIKE1]
  -S <float>  transition probability stay spike state [$SPIKE2]
  -w <int>    window size manual override
  -h          help (this usage statement)
" unless @ARGV;

$QUALITY = $opt_q if $opt_q;
$LOW     = $opt_l if $opt_l;
$TRANS   = $opt_t if $opt_t;
$LIMIT   = $opt_m if $opt_m;
$ZERO1   = $opt_z if $opt_z;
$ZERO2   = $opt_Z if $opt_Z;
$SPIKE1  = $opt_s if $opt_s;
$SPIKE2  = $opt_S if $opt_S;
$WSIZE   = $opt_w if $opt_w;

# Set up working directory
run("mkdir $WDIR") unless -d $WDIR;

#
# Extract coordinates from SAM file and save into individual files
#
my %FH;
my @filenames;
if (-s "$WDIR/filenames.txt") {
	open(IN, "$WDIR/filenames.txt") or die;
	while (<IN>) {
		chomp;
		push @filenames, $_;
	}
} else {
	open(OUT, ">$WDIR/filenames.txt") or die;
	foreach my $file (@ARGV) {
		if ($file =~ /\.gz$/) {open(IN, "gunzip -c $file |") or die}
		else                  {open(IN, $file) or die}
		while (<IN>) {
			next if /^@/;
			my ($foo, $bar, $chrom, $coor, $qual) = split;
			next if $qual < $QUALITY;
			if (not defined $FH{$chrom}) {
				$FH{$chrom} = new FileHandle;
				$FH{$chrom}->open(">$WDIR/$chrom.loc");
				print OUT "$WDIR/$chrom.loc\n";
				push @filenames, "$WDIR/$chrom.loc";
			}
			$FH{$chrom}->print("$coor\n");
		}
	}
	close OUT;
}

####################################
# Determine the proper window size #
####################################
my $max_count = 0;
my %zero;
my %logcount;
my $total_windows = 0;
foreach my $file (@filenames) {
	open(IN, "cpw $file 1024 |") or die;
	my $run = 0;
	while (<IN>) {
		$total_windows++;
		chomp;
		$logcount{int (log($_)/log(2))}++ if $_ != 0;
		if ($_ > $max_count) {$max_count = $_}
		if ($_ <= $LOW) {
			$run++;
		} else {
			if ($run) {
				$zero{$run}++;
				$run = 0;
			}
		}
	}
	close IN;
}

my $max_log = 0;
my $max_n;
foreach my $n (keys %logcount) {
	if ($logcount{$n} > $max_log) {
		$max_log = $logcount{$n};
		$max_n = $n;
	}
}

my $window = 1024 * 2**(3 - $max_n); # consider changing

open(my $wfh, ">$WDIR/windowsize.txt") or die;
if ($WSIZE) {
	print $wfh $WSIZE, "\n"; # manual override
} else {
	print $wfh $window, "\n";
}
close $wfh;

############################
# Create the diploid model #
############################
my @count2x;
foreach my $file (@filenames) {
	open(IN, "cpw $file $window |") or die;
	while (<IN>) {
		chomp;
		next if $_ == 0 or $_ > $LIMIT;
		$count2x[$_]++;
	}
	close IN;
}
my $freq2x = count2freq(\@count2x);

#########################
# Create derived models #
#########################
my %dname = (
	'1x' => 0.5,
	'3x' => 1.5,
	'4x' => 2,
);

# sampled models
my %model;
foreach my $name (sort keys %dname) {
	$model{$name} = resample($freq2x, $dname{$name});
}
$model{'2x'} = $freq2x;

# geometric models
for (my $i = 0; $i < $LIMIT; $i++) {
	$model{'Zero'}[$i] = 2 ** -($i+1);
	$model{'Spike'}[$i] = 2 ** -($LIMIT - $i);	
}

###########################
# Write the emission file #
###########################
open(OUT, ">$WDIR/emission.txt") or die;	
foreach my $mn (sort keys %model) {print OUT "$mn\t"}
print OUT "\n";

my %model_idx;
my $idx = 0;
foreach my $mn (sort keys %model) {$model_idx{$mn} = $idx++}

for (my $i = 0; $i < $LIMIT; $i++) {
	foreach my $mn (sort keys %model) {
		if (not defined $model{$mn}[$i]) {print OUT "0\t"}
		else                             {print OUT "$model{$mn}[$i]\t"}
	}
	print OUT "\n";
}
close OUT;

######################################
# Create parameter file for StochHMM #
######################################
my %state;
my $n = 0;
foreach my $mn (sort keys %model) {
	next if $mn eq 'Spike' or $mn eq 'Zero';
	$state{$mn} = $n++;
	$state{"$mn-S"} = $n++;
	$state{"$mn-Z"} = $n++;
}

my %matrix;
foreach my $n1 (keys %state) {
	foreach my $n2 (keys %state) {
		if ($n1 eq $n2) {
			if    ($n1 =~ /S$/) {$matrix{$n1}{$n1} = $SPIKE2}
			elsif ($n1 =~ /Z$/) {$matrix{$n1}{$n1} = $ZERO2}
			else                {$matrix{$n1}{$n1} = 1}
		} else {
			my $s1 = substr($n1, 0, 2);
			my $s2 = substr($n2, 0, 2);
			if ($s1 ne $s2) {
				if ($n1 !~ /[SZ]$/ and $n2 !~ /[SZ]$/) {
					$matrix{$n1}{$n2} = $TRANS;
				} else {
					# $matrix{$n1}{$n2} = 0;
				}
			}
			elsif ($n1 =~ /[SZ]$/ and $n2 =~ /[SZ]$/) {next}
			elsif ($n1 =~ /S$/) {$matrix{$n1}{$n2} = 1 - $SPIKE2}
			elsif ($n1 =~ /Z$/) {$matrix{$n1}{$n2} = 1 - $ZERO2}
			elsif ($n2 =~ /S$/) {$matrix{$n1}{$n2} = $SPIKE1}
			elsif ($n2 =~ /Z$/) {$matrix{$n1}{$n2} = $ZERO1}
		}
	}
}

my $DATE = `date`; chomp $DATE;
open(OUT, ">$WDIR/cnv12.hmm") or die;
print OUT "\
#STOCHHMM MODEL FILE

MODEL INFORMATION
======================================================
MODEL_NAME:	Dosamatic
MODEL_DESCRIPTION:	12 state model for CNV detection
MODEL_CREATION_DATE:	$DATE

TRACK SYMBOL DEFINITIONS
======================================================
SEQ:	REAL_NUMBER
";

print OUT "\
STATE DEFINITIONS
######################################################
STATE:	
	NAME:	INIT
TRANSITION:	STANDARD:	P(X)
";
my $default_init = 1 / (scalar keys %state);
foreach my $mn (sort keys %state) {
	print OUT "\t$mn\t$default_init\n";
}

foreach my $s1 (sort keys %state) {
	my $pdf;
	if    ($s1 =~ /^(\d)x$/) {$pdf = $s1}
	elsif ($s1 =~ /S$/)      {$pdf = 'Spike'}
	elsif ($s1 =~ /Z$/)      {$pdf = 'Zero'}
	else                     {$pdf = $s1}
	
	print OUT "######################################################\n";
	print OUT "STATE:\n";
	print OUT "\tNAME:\t$s1\n";
	print OUT "\tGFF_DESC:\t$s1\n";
	print OUT "\tPATH_LABEL: $state{$s1}\n";
	print OUT "TRANSITION:	STANDARD:	P(X)\n";
	foreach my $s2 (sort keys %{$matrix{$s1}}) {
		print OUT "\t$s2:\t$matrix{$s1}{$s2}\n";

	}
	print OUT "END:\t1\n";
	print OUT "EMISSION:	SEQ:	CONTINUOUS\n";
	print OUT "\tPDF:	$pdf	PARAMETERS:	$model_idx{$pdf}\n";
}
print OUT "######################################################\n";
print OUT "//END\n";


###############################################################################
# Subroutines
###############################################################################

sub resample {
	my ($f1, $x) = @_;
	
	my @pr;
	my $total = 0;
	for (my $i = 0; $i < @$f1; $i++) {
		push @pr, $f1->[$i] + $total;
		$total += $f1->[$i];
	}
	
	my @data;
	for (my $i = 0; $i < 1e6; $i++) {
		my $val;
		if ($x == 0.5) {
			my $rnd = rand(1);
			for (my $k = 0; $k < @pr; $k++) {
				if ($rnd <= $pr[$k]) {
					$val = int $k * 0.5;
					last;
				}
			}
		} elsif ($x == 1.5) {
			my $rnd = rand(1);
			for (my $k = 0; $k < @pr; $k++) {
				if ($rnd <= $pr[$k]) {
					$val = int $k * 0.5;
					last;
				}
			}
			
			my $rnd1 = rand(1);
			my $v = 0;
			for (my $k = 0; $k < @pr; $k++) {
				if ($rnd <= $pr[$k]) {
					$v = $k;
					last;
				}
			}
			$val += $v;
		} elsif ($x == 2) { # works with any integer > 1
			for (my $j = 0; $j < $x; $j++) {
				my $rnd = rand(1);
				my $v = 0;
				for (my $k = 0; $k < @pr; $k++) {
					if ($rnd <= $pr[$k]) {
						$v = $k;
						last;
					}
				}
				$val += $v;
			}
		} else {
			die "unsupported value of X";
		}
		push @data, $val unless $val >= $LIMIT;
	}
	
	my @count;
	foreach my $val (@data) {$count[$val]++}
	return count2freq(\@count);	
}

sub count2freq {
	my ($count) = @_;
	
	# define undefined values with LaPlace normalization
	for (my $i = 0; $i < @$count; $i++) {$count->[$i]++}
	
	# smooth by averaging neighbors
	my @smooth;
	for (my $i = 0; $i < @$count -1; $i++) {
		$smooth[$i] = $count->[$i-1] + $count->[$i+1] + $count->[$i];
	}
	$smooth[0] = $smooth[1];
	$smooth[@$count-1] = $smooth[@$count-2];
		
	# convert to freq
	my $total;
	for (my $i = 0; $i < @smooth; $i++) {$total += $smooth[$i]}
	for (my $i = 0; $i < @smooth; $i++) {$smooth[$i] /= $total}
	
	return \@smooth;
}

sub run {
	my ($cmd) = @_;
	system($cmd) == 0 or die "$cmd failed\n";
}
