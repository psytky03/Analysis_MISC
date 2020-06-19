#!/usr/bin/perl
use Modern::Perl;
use File::Path qw(make_path);

my $input = shift @ARGV; 
my $outdir = "/out";
make_path($outdir) unless -d $outdir;

my $outscript = "/out/scripts";
make_path($outscript) unless -d $outscript;

my $window = 2.5*1000*1000; #window size: 2.5Mb two sides  
my $threads = 10; 

my $header = <<EOF;
#!/bin/bash
#SBATCH --job-name=selscan
#SBATCH --partition=s5 --nodelist=pcom[21-22] -c 1
#SBATCH --qos=snp-normal
EOF

open (IN, "<", $input) or die "couldn't open the file $!\n";
while(<IN>){
    chomp;
    my ($ID, $gene) = (split /\t/, $_)[0,1];     
    my ($CHROM, $POS) = (split /_/, $ID)[0,1];
    my $start = $POS - $window;
    my $end = $POS + $window;
    my @commands;

    push (@commands, "bcftools view -r $CHROM:$start-$end -m2 -M2 -q '.01 minor' -Oz 1kgp.$CHROM.vcf.gz > $outdir/$gene.$ID.raw.vcf.gz");
    push (@commands, "~/tools/shapeit/shapeit2  --input-vcf $outdir/$gene.$ID.raw.vcf.gz -M ~/public/hg19physicalmap_selscan/genetic_map_chr$CHROM.txt -O $outdir/$gene.$ID --output-log $outdir/log.$gene.$ID.txt --thread $threads");
    push (@commands, "~/tools/shapeit/shapeit2  -convert --input-haps $outdir/$gene.$ID --output-vcf $outdir/$gene.$ID.vcf");
    push (@commands, "bgzip $outdir/$gene.$ID.vcf; tabix -p vcf $outdir/$gene.$ID.vcf.gz");
    push (@commands, "bcftools view $outdir/$gene.$ID.vcf.gz -H | awk '{print \$1, \$3, \"0\", \$2, \$4, \$5}' > $outdir/$gene.$ID.bim");
    push (@commands, "/opt/local/plink/1.90b3.38/plink --bim $outdir/$gene.$ID.bim --cm-map ~/public/hg19physicalmap_selscan/genetic_map_chr$CHROM.txt $CHROM --make-just-bim --out $outdir/$gene.$ID.gm");
    push (@commands, "cut -f 1-4 $outdir/$gene.$ID.gm.bim > $outdir/$gene.$ID.selscan.map");
    push (@commands, "~/tools/selscan/bin/linux/selscan --ihs --maf 0.01 --max-extend 0 --trunc-ok --vcf $outdir/$gene.$ID.vcf.gz --map $outdir/$gene.$ID.selscan.map --out $outdir/$gene.$ID.selscan");
    push (@commands, "~/tools/selscan/bin/linux/norm --ihs --files $outdir/$gene.$ID.selscan.ihs.out");
    
    open (OUT, ">", "$outscript/$gene.$ID.sh") or die "";
    print OUT $header;
    foreach (@commands){say OUT $_};

    system("chmod 755 $outscript/$gene.$ID.sh");
	system("sbatch -o $outdir/$gene.$ID.slurm.log $outscript/$gene.$ID.sh");

}
