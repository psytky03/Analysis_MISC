use strict; 
use warnings;
my $listA= shift @ARGV;

open(FH, '<', $listA);
my @ids = <FH>;
close(FH);
chomp @ids;
my %id_list = map { $_ => undef } @ids;

