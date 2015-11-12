#!/usr/bin/perl

## production version 4.0, last updated on 3.4.2014
## changes from v3: 
## -corrected the output for gene count (added whitespace trim)
## -corrected the output format to filter results by significance, and also to include "-inf" results correctly.
## -minor changes regarding output formatting.  


## set these directories to local values
my $bindir = "/home/apredeus/bin";
my $datadir = "/var/www/genequery/data";  
my $imgdir = "/genequery/images"; 

## CGI version of the query sript
## Uses temp file facility which stores files in /tmp 

use strict;
use warnings;
use CGI;
use CGI::Carp qw(fatalsToBrowser);
use File::Temp qw(tempfile); 
use HTML::Table; 

my $species=undef; 
my $db=undef;  
my $genes=undef;
my $q = CGI->new(\*STDIN);
print $q->header(-charset => 'UTF-8');
print $q->start_html;

my $values1 = File::Temp->new( UNLINK => 1, SUFFIX => '.tmp' );
my $values2 = File::Temp->new( UNLINK => 1, SUFFIX => '.tmp' );
my $fisher = File::Temp->new( UNLINK => 1, SUFFIX => '.tmp' );
my $output = File::Temp->new( UNLINK => 1, SUFFIX => '.tmp' );

$species = $q->param('species');
$db      = $q->param('db');
$genes   = $q->param('genes');
##print $q->b("Created the following temporary files: $values1 $values2 $fisher <br /><br />");
##print $q->b("Current directory is $dir <br /><br />");

if (defined $species && defined $db && defined $genes && !$genes =~ /^\s*$/) { 
  print $q->p("<b>processing data for the following species:</b><br /><br /> $species");
  print $q->b("processing the following list of IDs: <br /><br />");
  print $q->i("$genes<br /><br />");
  my $modlist="$datadir/Modules.$db.$species"; 
  my $names="$datadir/Names.$db.$species";
  my $reffile="$datadir/three_ids.$species"; 
  $genes =~ s/^\s+|\s+$//g;
  my @genes=split(/\s+/,$genes);  
  my @egenes1=();   ## array of gene IDs converted to Entrez IDs 
  my @egenes2=();  ## same as above but with all Entrez IDs unique - stats are also printed 
  my $qtype=undef; 
  
  my $hohoa={}; 
  my %names=();
  my %r2e = ();
  my %s2e = ();
  
  open MOD,$modlist or die "Error: Cannot open numerically sorted module list file!"; 
  open NAMES,$names or die "Error: Cannot open sample names file!"; 
  open REF,$reffile or die "Error: Cannot open reference file for species $species!"; 

  while (<NAMES>) {
    chomp;
    my @nn = split(/\t+/);
    $names{$nn[0]} = $nn[1];
  }
 
  while (<MOD>) { 
    chomp; 
    my @tt = split(/\s+/);
    push @{$hohoa->{$tt[0]}->{$tt[1]}},$tt[2]; 
  }

  while (<REF>) { 
    chomp; 
    my @tt = split (/\s+/); 
    $r2e{$tt[2]} = $tt[1]; 
    $s2e{uc($tt[3])} = $tt[1];
    ## uppercase just to be sure, symbols are often random and mixed   
  }

  if ($genes =~ m/^[0-9\s]*$/) { 
    print $q->p("<b>Genes are provided as <font color=\"FF0000\">Entrez IDs.</font></b><br />"); 
    $qtype = "entrez"; 
  } elsif ($genes =~ m/^[0-9\sNMR_]*$/) { 
    print $q->p("<b>Genes are provided as <font color=\"FF0000\">RefSeq IDs.</font></b><br />"); 
    $qtype = "refseq"; 
  } else { 
    print $q->p("<b>Genes are provided as <font color=\"FF0000\">gene symbols.</font></b><br />"); 
    $qtype = "symbol"; 
  }
  
  foreach my $gene (@genes) { 
    if ($qtype eq "entrez") { 
      push (@egenes1,$gene);  
    } elsif ($qtype eq "refseq" && $r2e{$gene} ne "") { 
      push (@egenes1,$r2e{$gene}); 
    } elsif ($qtype eq "symbol" && $s2e{uc($gene)} ne "") {
      push (@egenes1,$s2e{uc($gene)}); 
    }
  } 
  
  my %seen = ();
  foreach my $egene (@egenes1) {
    unless ($seen{$egene}) {
      push @egenes2,$egene;
      $seen{$egene} = 1;
    }
  }
  my $m1 = $#genes+1;
  my $m2 = $#egenes1+1;
  my $m3 = $#egenes2+1;
   
  print $q->p("<b>Query metrics:</b> sumbitted $m1 genes, found Entrez IDs for $m2, of them $m3 are unique.<br />"); 

  ## notice that we are using numerical sorting here (values must descend) 
  my @sortedgenes = sort {$b <=> $a} @egenes2; 
  my $ngenes = $#egenes2; 
  
  ## The part that does the search 
  foreach my $i (keys %{$hohoa}) {
    my $setoverlap=0; 
    my %overlaps=(); 
    my %modules=();
    foreach my $j (keys %{$hohoa->{$i}}) {
      my @module=@{$hohoa->{$i}->{$j}};
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
      printf $values1  "%s\t%d\t%d\t%d\t%d\n",$i,$overlaps{$j},$modules{$j}-$overlaps{$j},$setoverlap-$overlaps{$j},6000+$overlaps{$j}-$modules{$j}-$setoverlap;
      printf $values2  "%s\t%d\t%s\t%d\n",$i,$j,$modules{$j},$overlaps{$j};
    } 
  }
      
  ## now use that optimized Fisher's exact test calculation program
  my $dbid = undef; 
  if ($species eq "mm") {
    $dbid = ($db eq "2k") ? 1 : 2; 
  } else { 
    $dbid = ($db eq "2k") ? 3 : 4; 
  }
 ## print $q->p("<b>TEST: db is $db, dbid is $dbid, species is $species, m3 is $m3<br />");
  system "$bindir/fisher_test_right $values1 > $fisher"; 
  system "paste $values2 $fisher > $values1"; 
  system "$bindir/normal_pval_left $values1 $dbid $m3 > $values2";
  system "$bindir/process3.sh $values2 > $output"; 
 
  close MOD;
  close NAMES; 
  print $q->b("Following is the output of your analysis. <font color=\"FF0000\">Notice that \"-inf\" corresponds to a logp value of less than -325.</font> <br /><br /><br />"); 
  my $table = new HTML::Table();
  $table->setBorder(1); 
  $table->setWidth('85%'); 
  $table->addSectionRow ('thead',0,"Experiment ID","Platform","Experiment name","Module #","Genes in module","Overlap","Fisher p-value","log Fisher p-value","empirical p-value", "log empirical p-value");
  my $rc=2;
  if (-z $output) {
    print $q->b("<font color=\"FF0000\">No significant results were found for provided geneset.</font> <br /><br /><br />"); 
  } else { 
    while (<$output>) { 
      my @tt=split(/[_\s]+/);
      for (my $cc=0; $cc<=$#tt; $cc++) { 
        if ($cc==0) { 
          my $url="http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=".$tt[0]; 
          my $link="<a href=\"".$url."\">".$tt[0]."</a>";
          $table->setCell($rc,1,$link);
        } elsif ($cc==1) {
          my $url="http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=".$tt[1]; 
          my $link="<a href=\"".$url."\">".$tt[1]."</a>";
          $table->setCell($rc,2,$link);
        } elsif ($cc==2) {
          my $fullid=$tt[0]."_".$tt[1]; 
          my $url="$imgdir/$db/$fullid.png"; 
          my $link="<a href=\"".$url."\">".$tt[2]."</a>";
          $table->setCell($rc,3,$names{$fullid});
          $table->setCell($rc,4,$link);
        } else { 
          $table->setCell($rc,$cc+2,$tt[$cc]);
        } 
      } 
      $rc++;   
    } 
    $table->print;
  } 
} elsif (defined $species && defined $db && (!defined $genes || $genes =~ /^\s*$/)) { 
  die("ERROR: You have to provide a non-empty gene list!<br />");
} else { 
  die("ERROR: You have to specify the species and the database to proceed!<br />");
} 


print $q->end_html; 
