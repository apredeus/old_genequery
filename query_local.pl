#!/usr/bin/perl 

use warnings; 
use strict; 
use Time::HiRes qw (time);

my $t1=time;
my $dbase=shift @ARGV; 
my $query=shift @ARGV; 
my @genes;  
my $line="";

my $hohoa={}; 

open GL,"<",$query or die "$!"; 
open DB,"<",$dbase or die "$!"; 

while (<GL>) { 
  chomp; 
  push @genes,$_;
  $line=$line.$_;  
} 

## Pre-processing gene lists is complex and thus will be a separate module. Here we strictly require Entrez ID format. 
my $type=($line=~/^\d+$/)?"entrez":"not_entrez"; 
print STDERR "It seems like gene list is provided in Entrez format! Continuing..\n"; 
die "ERROR: this version only supports Entrez format query!" if ($type ne "entrez"); 

## Simple and elegant way to set up our "in-memory database" 
while (<DB>) { 
  chomp; 
  my @tt = split(/\s+/);
  push @{$hohoa->{$tt[0]}->{$tt[1]}},$tt[2]; 
}
my $t2=time;

printf STDERR "parsed module list, created database, time taken: %f\n",$t2-$t1;
open VALUES,">","$$.values.tmp" or die "$!";
open VALUES2,">","$$.values2.tmp" or die "$!";

## notice that we are using numerical sorting here (values descend) 
my @sortedgenes = sort {$b <=> $a} @genes; 
my $ngenes = $#genes; 

## The part that does the search 
foreach my $i (keys %{$hohoa}) {
  my $setoverlap=0; 
  my %overlaps; 
  my %modules;
  foreach my $j (keys %{$hohoa->{$i}}) {
    my @module=@{$hohoa->{$i}->{$j}};
    #my @sortedmodule = sort {$b <=> $a} @module; 
    #my $s1 = join ('',@sortedmodule); 
    #my $s2 = join ('',@module); 
    #die "Error: Failed the module sort test at set id $i, module number $j\n" if ($s1 ne $s2); 
    
    ## very efficient overlap calculation due to sorted arrays 
    my $m=0;
    my $n=0;
    my $overlap=0;
    while ($m<=$#module && $n<=$#sortedgenes) {
      if ($module[$m] == $sortedgenes[$n]) {
        $overlap++;
        $m++;
        $n++;
      } elsif ($module[$m] < $sortedgenes[$n]) {
        $n++;
      } else {
        $m++;
      }
    } 
   $overlaps{$j}=$overlap;
   $modules{$j}=$#module+1; 
   $setoverlap+=$overlap;
  }
  
  ## now have to print - loop again 
  foreach my $j (keys %{$hohoa->{$i}}) {  
    printf VALUES  "%s\t%d\t%d\t%d\t%d\n",$i,$overlaps{$j},$modules{$j}-$overlaps{$j},$setoverlap-$overlaps{$j},6000+$overlaps{$j}-$modules{$j}-$setoverlap;
    printf VALUES2 "%s\t%d\t%s\t%d\t%d\n",$i,$j,$setoverlap,$modules{$j},$overlaps{$j};
  } 
}
    
close VALUES; 
close VALUES2; 
my $t3=time;
printf STDERR "Writing the values file took %f s\n",$t3-$t2;   

## now use that optimized Fisher's exact test calculation program
system "fisher_test_right $$.values.tmp > $$.fisher.tmp"; 
my $t4=time; 
printf STDERR "Fisher's exact test calulations took %f s\n",$t4-$t3;   
system "paste $$.values2.tmp $$.fisher.tmp | sort -k7,7n"; 


close GL; 
close DB; 
unlink "$$.values.tmp" or die "$!"; 
unlink "$$.values2.tmp" or die "$!"; 
unlink "$$.fisher.tmp" or die "$!"; 
my $t5=time;
printf STDERR "Time used for script execution: %f\n",$t5-$t1; 
