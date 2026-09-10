unless(@ARGV==2)
{
	die"Usage:perl $0 <organism> <path>\n";
}
my $spc_file=shift;
my $path=shift;

open IN,$spc_file||die$!;
my %kingdom_hash;
while(<IN>)
{
	chomp;
	my @arr=split(/\t/);
	my $kingdom=(split(/;/,$arr[3]))[1];
	$kingdom_hash{$kingdom}{$arr[1]}=0;
}
close IN;

my %ko_hash;
foreach my $site (sort keys %kingdom_hash)
{
	my $hash=$kingdom_hash{$site};
	foreach my $item (keys %$hash)
	{
		my $kegg=$path."/".$item;
		if(-e $kegg)
		{
			open INN,$kegg||die$!;
			while(<INN>)
			{
				my $id=(split(/\s+/))[0];
				$ko_hash{$site}{$id}=0;
			}
		}
	}
}


foreach my $site (sort keys %ko_hash)
{
	my $out=$site.".ko.txt";
	open OUT,">$out"||die$!;
	my $hash=$ko_hash{$site};
	foreach my $item (keys %$hash)
	{
		print OUT $item,"\n";
	}
	close OUT;
}

