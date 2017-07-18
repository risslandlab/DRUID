#!/usr/bin/env perl

use v5.10.0;

use strict;
use warnings;

use autodie;

my $keep = 0;

while(<>){
    if( /^>/ ){
        $keep = ( /^>chr[0-9A-Z]{1,2}$/ ) ? 1 : 0;
    }
    if ($keep){
        print $_;
    }
}

           

