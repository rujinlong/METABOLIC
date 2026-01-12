#!/usr/bin/env perl
# merge_hmm_profiles.pl
# Merges individual HMM files into a single database for efficient hmmpress/hmmsearch
# Usage: perl merge_hmm_profiles.pl <input_dir> <output_file> [filter_list]

use strict;
use warnings;

my $input_dir = $ARGV[0] or die "Usage: perl merge_hmm_profiles.pl <input_dir> <output_file> [filter_list]\n";
my $output_file = $ARGV[1] or die "Usage: perl merge_hmm_profiles.pl <input_dir> <output_file> [filter_list]\n";
my $filter_list = $ARGV[2]; # Optional: list of HMM files to include

# Load filter list if provided
my %include_hmm = ();
if ($filter_list && -e $filter_list) {
    open FILTER, "$filter_list" or die "Cannot open filter list: $filter_list\n";
    while (<FILTER>) {
        chomp;
        s/^\s+|\s+$//g;
        $include_hmm{$_} = 1 if $_;
    }
    close FILTER;
    print "Loaded " . scalar(keys %include_hmm) . " HMM profiles from filter list\n";
}

# Merge HMM files
my $count = 0;
open OUT, ">$output_file" or die "Cannot open output file: $output_file\n";
open IN, "ls $input_dir/*.hmm 2>/dev/null |";
while (<IN>) {
    chomp;
    my $hmm_file = $_;
    my ($basename) = $hmm_file =~ /([^\/]+)$/;
    
    # Skip if filter list provided and this file is not included
    if (%include_hmm && !$include_hmm{$basename}) {
        next;
    }
    
    open HMM, "$hmm_file" or next;
    while (<HMM>) {
        print OUT $_;
    }
    close HMM;
    $count++;
}
close IN;
close OUT;

print "Merged $count HMM profiles into $output_file\n";
