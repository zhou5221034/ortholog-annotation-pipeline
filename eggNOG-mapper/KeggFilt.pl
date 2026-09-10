unless(@ARGV==2)
{
	die"Usage:perl $0 <KEGG.universal.tsv> <Animals/Plants/Fungi/Bacteria/Archaea/Protists>\n";
}
my $kegg_file=shift;
my $kegg_lib=shift;

$kegg_lib="/gpfs/users/liuzm/software/eggNOG-mapper/".$kegg_lib.".ko.txt";
open IN,$kegg_lib||die$!;
my %kegg_hash;
while(<IN>)
{
	my @arr=split(/\s+/);
	$kegg_hash{$arr[0]}=0;
}
close IN;

open IN,$kegg_file||die$!;
open OUT,">$kegg_file.filted"||die$!;
my $totals=0;
my $remains=0;
while(<IN>)
{
	my @arr=split(/\t/);
	if($arr[0] eq "Query")
	{
		print OUT $_;
	}
	else
	{
		$totals++;
		if(exists $kegg_hash{$arr[2]})
		{
			print OUT $_;
			$remains++;
		}
	}
}
close IN;
close OUT;
print "totals:\t$totals\n";
print "remains:\t$remains\n";
