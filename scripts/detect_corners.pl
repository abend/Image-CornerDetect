#!/usr/bin/perl -w
use strict;
use warnings;

use FindBin qw($Bin);
use lib "$Bin/../lib";
use Image::CornerDetect;
use Data::Dumper;

my $file = $ARGV[0];
die "usage: $0 filename" unless $file;

my $detector = Image::CornerDetect->new({write_debug_image => 1});
print Dumper($detector->process($file));
