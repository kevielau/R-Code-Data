
use strict;
use warnings;

my %hash=();

#��ȡ�ٴ��ļ��������浽hash����
open(RF,"clinical_TCGA.txt") or die $!;
while(my $line=<RF>){
	chomp($line);
	my @arr=split(/\t/,$line);
	my $sample=shift(@arr);
	if($.==1){
		$hash{"id"}=join("\t",@arr);
		next;
	}
	$hash{$sample}=join("\t",@arr);
}
close(RF);

#��ȡTMB�ļ����������ٴ���Ϣ����������"tmbClinical.txt"
open(RF,"singleGene_VDR_TCGA.txt") or die $!;
open(WF,">singleGeneClinical_VDR_TCGA.txt") or die $!;
while(my $line=<RF>){
	chomp($line);
	my @arr=split(/\t/,$line);
	my $sample=shift(@arr);
	my $tumor=pop(@arr);
	#my $Correlation=pop(@arr);
	#my $Pvalue=pop(@arr);
	my @samp1e=(localtime(time));
	if($.==1){
		print WF "id\t$hash{\"id\"}\t" . join("\t",@arr) . "\n";
		next;
	}
	my @sampleArr=split(/\-/,$sample);
	if($sampleArr[3]=~/^0/){
		my $sampleName="$sampleArr[0]-$sampleArr[1]-$sampleArr[2]";
		if(exists $hash{$sampleName}){
			print WF "$sample\t$hash{$sampleName}\t" . join("\t",@arr) . "\n";
			delete($hash{$sampleName});
		}
	}
}
close(WF);
close(RF);
