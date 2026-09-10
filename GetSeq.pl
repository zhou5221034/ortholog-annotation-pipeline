unless(@ARGV==2)
{
	die"Usage:perl $0 <id> <fasta>\n";
}
open IN,$ARGV[0]||die$!;
my %id_hash;
while(<IN>)
{
	my @arr=split(/\s+/);
	$id_hash{$arr[1]}=$arr[0];
}
close IN;

open IN,$ARGV[1]||die$!;
$/=">";<IN>;
while(<IN>)
{
	chomp;
	my $id=(split(/\s+/))[0];
	$_=~s/^.+?\n//g;
	$_=~s/\s+//g;
	if(exists $id_hash{$id})
	{
		print ">",$id_hash{$id},"\n",$_,"\n";
	}
}
close IN;
$/="\n";

