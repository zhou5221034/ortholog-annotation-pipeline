#!/usr/bin/perl -w
use strict;
use Cwd qw(abs_path);
use File::Basename qw(basename dirname);

unless(@ARGV>=3)
{
	die"Usage:perl $0 <reference> <query> <output> [coverage] [identity]\n";
}

my $reference=shift;
my $query=shift;
my $output=shift;
my $cov=shift;
my $iden=shift;
if(!$cov)
{
	$cov=0.3;
}
else
{
	if($cov>1 || $cov<=0)
	{
		die"coverage must be (0,1],explampe 0.3\n";
	}
}
if(!$iden)
{
	$iden=30;
}
else
{
	if($iden<=0 || $iden>100)
	{
		die"identity must be (0,100],explampe 30\n";
	}
}

#open OUT,">$output/tmp.fa"||die$!;

&DealFa($reference,"ref");
&DealFa($query,"tar");

#close OUT;

open CFG,">$output/config"||die$!;
print CFG "#species type\nref\tref\ntar\ttar\n#--END--\n\n#alignment cutoff\n#sp1 sp2 alignRate identity\nref\ttar\t$cov\t$iden\n#--END--\n";
close CFG;

open SH,">$output/runBest.sh"||die$!;
print SH "perl /gpfs/users/liuzm/software/Best2Best/runblastp.pl $output/ref.fa $output/tar.fa 10 $output\n";
print SH "perl /gpfs/users/liuzm/software/Best2Best/multi-process.pl -cpu 10 $output/blast_shell.sh\n";
print SH "cat $output/tar.fa.cut/*.m8 > $output/all_vs_all.1.m8\n";
print SH "cat $output/ref.fa $output/tar.fa > $output/tmp.fa\n";
print SH "perl /gpfs/users/liuzm/software/Best2Best/runblastp.pl $output/tar.fa $output/ref.fa 10 $output\n";
print SH "perl /gpfs/users/liuzm/software/Best2Best/multi-process.pl -cpu 10 $output/blast_shell.sh\n";
print SH "cat $output/ref.fa.cut/*.m8 > $output/all_vs_all.2.m8\n";
print SH "cat $output/all_vs_all.1.m8 $output/all_vs_all.2.m8 > $output/all_vs_all.m8\n";
print SH "perl /gpfs/users/liuzm/software/Best2Best/select_m8.pl $output/all_vs_all.m8 $output/config $output\n";
print SH "/gpfs/users/liuzm/software/Best2Best/solar.pl -a prot2prot -f m8 $output/ref_tar.m8 > $output/ref_tar.m8.solar\n";
print SH "perl /gpfs/users/liuzm/software/Best2Best/solar_add_realLen.pl $output/ref_tar.m8.solar $output/tmp.fa > $output/ref_tar.m8.solar.cor\n";
print SH "perl /gpfs/users/liuzm/software/Best2Best/solar_add_identity.pl --solar $output/ref_tar.m8.solar.cor --m8 $output/ref_tar.m8 > $output/ref_tar.m8.solar.cor.idAdd\n";
print SH "perl /gpfs/users/liuzm/software/Best2Best/ort1.pl $output/ref_tar.m8.solar.cor.idAdd ref tar > $output/ref_tar.m8.solar.cor.idAdd.rbh\n";
print SH "sed 's/_ref//g' $output/ref_tar.m8.solar.cor.idAdd.rbh |sed 's/_tar//g' |awk '{if(\$7==\"L1\") print}' > $output/result\n";
print SH "rm -r $output/all_vs_all.m8 $output/ref_tar.m8 $output/ref_tar.m8.solar $output/ref_tar.m8.solar.cor $output/ref_tar.m8.solar.cor.idAdd  $output/tar.fa* $output/ref.fa* $output/tmp.fa* $output/blast_shell.sh* $output/config $output/all_vs_all.1.m8 $output/all_vs_all.2.m8\n";
close SH;
#`perl /public/work/Personal/liuzhiming/software/Comp_Gene/best2best/new_bin/step1.pl $output/all_vs_all.m8 $output/tmp.fa config ./ 3`;


sub DealFa
{
	my ($file,$tag)=@_;
	open IN,$file||die$!;
	open OUT,">$output/$tag.fa"||die$!;
	$/=">";<IN>;
	while(my $seq=<IN>)
	{
		chomp $seq;
		my $id=(split(/\s+/,$seq))[0];
		$seq=~s/^.+?\n//g;
		$seq=~s/\s+//g;
		print OUT ">",$id,"_",$tag,"\n",$seq,"\n";
	}
	close IN;
	close OUT;
}
