#!/usr/bin/perl
use strict;
use 5.010;
# perl hg19_50kb_bin.pl > hg19_50kb_bin.bed
# https://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes
my %chr_size = qw(
chr1	249250621
chr2	243199373
chr3	198022430
chr4	191154276
chr5	180915260
chr6	171115067
chr7	159138663
chrX	155270560
chr8	146364022
chr9	141213431
chr10	135534747
chr11	135006516
chr12	133851895
chr13	115169878
chr14	107349540
chr15	102531392
chr16	90354753
chr17	81195210
chr18	78077248
chr20	63025520
chrY	59373566
chr19	59128983
chr22	51304566
chr21	48129895
);

my $bin_size = 50000;

foreach my $chr (sort keys %chr_size){
    my $curr_size = $chr_size{$chr};
    my $chr_num = $chr;
    $chr_num =~ s/^chr//;
    for(my $i = 0; $i < $curr_size; $i += $bin_size+1){
        my $end = $i + $bin_size;
        if ($end > $curr_size){$end = $curr_size};      
        say "$chr_num\t$i\t$end";
    }
}

