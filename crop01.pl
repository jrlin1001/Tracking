## Generate trace files for thumbnail 
## Jerry Lin 20140912

use warnings;
use strict;

my $well="002002";
if($ARGV[0]){$well=$ARGV[0];}

## READ POS FILE
open POS,"<output.pos-$well.csv";
my @poss;
my $i=1;
my $total =1;

while(<POS>){
  my $line =$_; chomp $line;
  if($line){$poss[$i]=$line;$i++;$total++;}
}
$total--;
close(POS);


## READ SLICE FILE
open SLICE,"<output.slice-$well.csv";
my @slices;

$i=1;
while(<SLICE>){
  my $line =$_; chomp $line;
  if($line){$slices[$i]=$line;$i++;}
}
close(SLICE);

## READ RATIO FILE
open RATIO,"<output.ratio-$well.csv";
my @ratios;

$i=1;
while(<RATIO>){
  my $line =$_; chomp $line;
  if($line){$ratios[$i]=$line;$i++;}
}
close(RATIO);

##WRITE OUTPUT FILES

system("mkdir $well");
system("cd $well");

for($i=1;$i<=$total;$i++){
  open TRACE,">.\\$well\\trace-$i.txt";
  if($slices[$i]){
  my @s=split ',',$slices[$i];
  my @p=split ',',$poss[$i];
  my @r=split ',',$ratios[$i];
  for(my $j=0;$j<@s;$j++){
    print TRACE "$s[$j]\t";
    my $tp = $p[$j];
    $tp =~s/-/\t/;
    print TRACE "$tp\t";
    print TRACE "$r[$j]\n";
  }#if(slices)
  }#for(j)
  close(TRACE);
}#for(i)








