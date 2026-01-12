#!/usr/bin/perl

use strict;
use warnings;

open IN, "ls *.txt |";
while (<IN>){
	chomp;
	`hmmpress $_`;
}
close IN;