#!/usr/bin/perl
use strict;
use warnings;

die "usage: dosamatic.pl <directory from hmm_gen.pl>\n"
	unless @ARGV == 1;
my ($DIR) = @ARGV;

my $windowsize = `cat $DIR/windowsize.txt`;
chomp $windowsize;

my @file = `cat $DIR/filenames.txt`;
chomp @file;

open(OUT, ">$DIR/seq.fa") or die;
foreach my $file (@file) {
	my ($name) = $file =~ /(\w+)\.loc$/;
	print OUT ">$name\n";
	open(IN, "cpw $file $windowsize |") or die;
	while (<IN>) {print OUT}
	close IN;
}
close OUT;

system("dosamatic $DIR/cnv12.hmm $DIR/emission.txt $DIR/seq.fa > $DIR/seq.win")
	== 0 or die;
#system("gff_to_bed.py -d $DIR") == 0 or die;


# fix the GFF (window size and strip -S and -X)
my %gff;
open(IN, "$DIR/seq.win") or die;
while (<IN>) {
	next unless /\S/;
	my ($chrom, $tool, $state, $beg, $end) = split;
	push @{$gff{$chrom}}, {
		state => $state,
		beg => $beg,
		end => $end,
	};
}
close IN;

# create GFF output
open(OUT, ">$DIR/seq.gff") or die;
foreach my $chrom (sort keys %gff) {
	for (my $i = 0; $i < @{$gff{$chrom}}; $i++) {
		my $s1 = substr($gff{$chrom}[$i]{state}, 0, 2);
		for (my $j = $i +1; $j < @{$gff{$chrom}} -1; $j++) {
			my $s2 = substr($gff{$chrom}[$j]{state}, 0, 2);
			if ($s1 ne $s2) {
				print OUT join("\t", $chrom, 'dosamatic', $s1,
					$gff{$chrom}[$i  ]{beg} * $windowsize,
					$gff{$chrom}[$j-1]{end} * $windowsize,
					'.', '+', '.',
					$gff{$chrom}[$j-1]{end} * $windowsize - $gff{$chrom}[$i  ]{beg} * $windowsize), "\n";
				$i = $j - 1;
				last;
			}
		}
	}
}
