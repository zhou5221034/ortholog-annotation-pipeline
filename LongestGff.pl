unless(@ARGV==1)
{
	die"Usage:perl $0 <gff> \n";
}

open IN,$ARGV[0]||die$!;
my %len_hash;
my %gene_hash;
while(<IN>)
{
	next if($_=~/#/);
	my @arr=split(/\s+/);
	if($arr[2] eq "mRNA" || $arr[2] eq "transcript")
	{
		my $mRNA=$1	if($_=~/ID=(\S+?)(;|$)/);
		my $gene=$1	if($_=~/Parent=(\S+?)(;|$)/);
		$gene_hash{$gene}{$mRNA}=0;
	}
	if($arr[2] eq "CDS")
	{
		my $mRNA=$1	if($_=~/Parent=(\S+?)(;|$)/);
		my $len=$arr[4]-$arr[3]+1;
		$len_hash{$mRNA}+=$len;
	}
}
close IN;

my %mRNA_hash;
my %longest_hash;
foreach my $site (sort keys %gene_hash)
{
	my $hash=$gene_hash{$site};
	foreach my $item (sort keys %$hash)
	{
		if(!exists $longest_hash{$site})
		{
			$mRNA_hash{$site}=$item;
			$longest_hash{$site}=$len_hash{$item};
		}
		else
		{
			if($longest_hash{$site} < $len_hash{$item})
			{
				$mRNA_hash{$site}=$item;
				$longest_hash{$site}=$len_hash{$item};
			}
		}
	}
}

open MAPID,">mapid"||die$!;

foreach my $site (keys %mRNA_hash)
{
	print MAPID $site,"\t",$mRNA_hash{$site},"\t",$longest_hash{$site},"\n";
}

close MAPID;
