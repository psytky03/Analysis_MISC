#!/usr/bin/perl
use strict;
use 5.010;

# Purpose
# calculate the expected missense/synonymous/non-sense mutation ratio 

# output (without mutation freq weight)
# perl calculate_mutation_possibility.pl | grep mis | awk '{sum1 += $9; sum2 += $10} END {print sum1, sum2}'
# perl calculate_mutation_possibility.pl | grep syn | awk '{sum1 += $9; sum2 += $10} END {print sum1, sum2}'
# perl calculate_mutation_possibility.pl | grep non | awk '{sum1 += $9; sum2 += $10} END {print sum1, sum2}'
# 6544.3 6560.52 (missense)
# 2058.9 2046.5 (synonymous)
# 368 378.67 (non sense)
# dN/dS = 6544.3/2058.9 = 3.18
# dS/dLoF = 2058.9/368 = 5.59

# output (with mutatioon freq weight)
# perl calculate_mutation_possibility.pl | grep mis | awk '{sum1 += $11; sum2 += $12} END {print sum1, sum2}'
# perl calculate_mutation_possibility.pl | grep syn | awk '{sum1 += $11; sum2 += $12} END {print sum1, sum2}'
# perl calculate_mutation_possibility.pl | grep non | awk '{sum1 += $11; sum2 += $12} END {print sum1, sum2}'
# 100176 100336 (missense)
# 44449.3 43890.6 (synonymous)
# 5891.38 5993.17 (non sense)
# dN/dS = 100176/44449.3 = 2.25
# dS/dLoF = 44449.3/5891 = 7.55


# 99517.3 99324.5 (missense)
# 46846.6 45207.1 (synonymous)
# 6801.97 6884.69 (non sense)
# 99517/46846
# 46846.6/6801
my %codon2aa = qw(
  TCA  S  TCC  S  TCG  S  TCT  S  TTC  F  TTT  F  TTA  L  TTG  L
  TAC  Y  TAT  Y  TAA  *  TAG  *  TGC  C  TGT  C  TGA  *  TGG  W
  CTA  L  CTC  L  CTG  L  CTT  L  CCA  P  CCC  P  CCG  P  CCT  P
  CAC  H  CAT  H  CAA  Q  CAG  Q  CGA  R  CGC  R  CGG  R  CGT  R
  ATA  I  ATC  I  ATT  I  ATG  M  ACA  T  ACC  T  ACG  T  ACT  T
  AAC  N  AAT  N  AAA  K  AAG  K  AGC  S  AGT  S  AGA  R  AGG  R
  GTA  V  GTC  V  GTG  V  GTT  V  GCA  A  GCC  A  GCG  A  GCT  A
  GAC  D  GAT  D  GAA  E  GAG  E  GGA  G  GGC  G  GGG  G  GGT  G
);

# Frequency of codon usage obtained from 
# https://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi?species=9606
my %codon2usage = qw(
TTT 17.6  TCT 15.2  TAT 12.2  TGT 10.6
TTC 20.3  TCC 17.7  TAC 15.3  TGC 12.6
TTA  7.7  TCA 12.2  TAA  1.0  TGA  1.6
TTG 12.9  TCG  4.4  TAG  0.8  TGG 13.2
CTT 13.2  CCT 17.5  CAT 10.9  CGT  4.5
CTC 19.6  CCC 19.8  CAC 15.1  CGC 10.4
CTA  7.2  CCA 16.9  CAA 12.3  CGA  6.2
CTG 39.6  CCG  6.9  CAG 34.2  CGG 11.4
ATT 16.0  ACT 13.1  AAT 17.0  AGT 12.1
ATC 20.8  ACC 18.9  AAC 19.1  AGC 19.5
ATA  7.5  ACA 15.1  AAA 24.4  AGA 12.2
ATG 22.0  ACG  6.1  AAG 31.9  AGG 12.0
GTT 11.0  GCT 18.4  GAT 21.8  GGT 10.8
GTC 14.5  GCC 27.7  GAC 25.1  GGC 22.2
GTA  7.1  GCA 15.8  GAA 29.0  GGA 16.5
GTG 28.1  GCG  7.4  GAG 39.6  GGG 16.5
);

# Frequency of codon usage obtained from 
# https://hive.biochemistry.gwu.edu/dna.cgi?cmd=tissue_codon_usage&id=586358&mode=cocoputs
my %codon2usage2 = qw(
TTT	17.22		TCT	16.96		TAT	12.15		TGT	10.47
TTC	17.51		TCC	17.29		TAC	13.47		TGC	10.83
TTA	 8.82		TCA	14.19		TAA	 0.44		TGA	 0.80
TTG	13.47		TCG	 4.05		TAG	 0.35		TGG	11.66
CTT	14.17		CCT	19.15		CAT	11.90		CGT	 4.56
CTC	17.81		CCC	18.86		CAC	14.62		CGC	 8.71
CTA	 7.48		CCA	18.79		CAA	14.17		CGA	 6.42
CTG	36.01		CCG	 6.15		CAG	35.28		CGG	10.62
ATT	16.58		ACT	14.30		AAT	18.52		AGT	14.06
ATC	18.68		ACC	17.77		AAC	18.27		AGC	19.66
ATA	 8.16		ACA	16.57		AAA	27.77		AGA	13.40
ATG	21.48		ACG	 5.59		AAG	31.79		AGG	12.17
GTT	11.80		GCT	18.89		GAT	24.09		GGT	10.82
GTC	13.44		GCC	25.64		GAC	24.25		GGC	19.73
GTA	 7.71		GCA	17.09		GAA	33.85		GGA	17.19
GTG	25.75		GCG	 5.94		GAG	39.40		GGG	15.26
);

# Human germline mutation rate obtained from Table 1 of 
# Milholland, Brandon, et al. "Differences between germline and somatic mutation rates in humans and mice." 
# Nature communications 8.1 (2017): 1-8.
my %codonfreq = qw(
    AC	5.87
    AT	9.68
    AG	27.5
    CA	9.78
    CT	41.21
    CG	5.97
    TA	9.68
    TC	27.5
    TG	5.87
    GA	41.21
    GC	5.97
    GT	9.78
);

# my %codonfreq = qw(
#     AC	3.97
#     AT	3.02
#     AG	17.30
#     CA	16.13
#     CT	51.47
#     CG	8.12
#     TA	3.02
#     TC	17.30
#     TG	3.97
#     GA	51.47
#     GC	8.12
#     GT	16.13
# );


my @dnas = qw(A C T G);

foreach my $codon (keys %codon2aa){
    my $aa = $codon2aa{$codon};
    my $aafreq1 = $codon2usage{$codon};
    my $aafreq2 = $codon2usage2{$codon};
    foreach my $i (0..2){
        my $pos = $i+1;
        foreach my $dna (@dnas){
            my $codon_mut = $codon;
            substr($codon_mut, $i, 1, $dna);
            next if $codon eq $codon_mut;
            my $aa_mut = $codon2aa{$codon_mut};

            my $tvti = substr($codon, $i, 1).substr($codon_mut, $i, 1);
            my $tvti_freq = $codonfreq{$tvti};

            my $aafreq1_wight = $aafreq1 * $tvti_freq;
            my $aafreq2_wight = $aafreq2 * $tvti_freq;

            my $mut_type;
            if ($aa eq "*"){
                $mut_type = "-"
            }elsif(
                $aa eq $aa_mut
            ){
                $mut_type = "syn"
            }elsif(
                ($aa ne $aa_mut) & ($aa_mut ne "*")
            ){
                $mut_type = "mis"
            }elsif(
                 ($aa ne $aa_mut) & ($aa_mut eq "*")
            ){
                $mut_type = "non"
            }
            say "$codon $aa $pos $dna $codon_mut $aa_mut $tvti $mut_type $aafreq1 $aafreq2 $aafreq1_wight $aafreq2_wight";
        }
    }
}
