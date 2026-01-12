#!/usr/bin/perl

use strict;
use warnings;

my $file = "pepunit.lib";

`perl -pi -e 's/\r\n/\n/g' $file`;

my %H = (); my $head;
open IN, "$file";
while (<IN>){
	chomp;
	if (/>/){
		$_ =~ tr/\015//d;
		$head = $_;
		$H{$head} = "";
	}else{
		$_ =~ tr/\015//d;
		$_ =~ s/ //g;
		$H{$head} .= $_;
	}	
}
close IN;

open OUT, ">tmp.seq";
foreach my $key (sort keys %H){
	print OUT "$key\n$H{$key}\n";
}
close OUT;


`diamond makedb --in tmp.seq -d pepunit.db`;

`rm tmp.seq`;