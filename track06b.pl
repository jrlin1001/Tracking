## tracking script (Track06b.pl)
## Jerry Lin 2014/05/19

##
## for tracking two channel image
## nucleus and cytoplasmic ratio
## Usage: Perl trakc06.pl <begin site> <end site>
## Input files: Measure-mask-sXX.csv & Measure-Cytosol-sXX.csv
## Output files: Tracks-sXX.tsv & Tracking-sXX.csv
## additional output files: debug-sXX.log (XY position)  output.mean-sXX.csv (mean_nuc only output) output.ratio-sXX.csv (ratio only output)
## 06b update (05/19/2014)--> including area threshold to reject wrong tracks
##

use warnings;
use strict;

my $begin =1;
my $end =60;
if($ARGV[0]){$begin=$ARGV[0];}
if($ARGV[1]){$end=$ARGV[1];}

for(my $site=$begin;$site<$end+1;$site++)
{

##=======start tracking routine=========
my $filename='measure-mask-s'.$site.'.csv';
my $output = 'tracking-s'.$site.'.tsv';
my $trackf = 'tracks-s'.$site.'.csv';
my $cytof ='measure-cytosol-s'.$site.'.csv';



my $debugfile = 'debug-s'.$site.'.log';

open INFILE,"<$filename";
open DEBUG,">$debugfile";
open OUTFILE,">$output";
open TRACK,">$trackf";
open CYTO,"<$cytof";

my $debug = 1;
my %p2s;
my %p2i;

## set default positions
my $p_slice = 12;
my $p_x = 6;
my $p_y = 7;
my $p_area = 1;
my $p_mean =2;
my @cc;
my @matrix;
## ----------------------


##define parameter
my $limit_d =25;
##----------------

##Read Cytosol measurement
my $title = <CYTO>;
my @col = split ',',$title;
for(my $i=0; $i<@col; $i++){                   #analyze columns
  if($col[$i] eq "Area"){$p_area = $i;}
  if($col[$i] eq "Mean"){$p_mean = $i;}
}
my %measure;                       ## hash for cytosol measurement
my %mean_n;
my %ratio_nc;
my $index=1;

while(<CYTO>){
  my $line = $_;
  chomp $line;
  if($line){
    my @col = split ',',$line;
    my $area_nuc = $col[$p_area];
    my $mean_nuc = $col[$p_mean];
    $line = <CYTO>;
    chomp $line;
    @col = split ',',$line;
    my $area_cyto = $col[$p_area];
    my $mean_cyto = $col[$p_mean];

    $mean_cyto = ($mean_cyto*$area_cyto - $mean_nuc*$area_nuc)/($area_cyto-$area_nuc);
    $mean_n{$index}=$mean_nuc;

    $area_cyto =$area_cyto-$area_nuc;

    my $ratio = $mean_nuc/$mean_cyto;
    $ratio_nc{$index}=$ratio;

    $measure{$index}=$area_nuc.','.$mean_nuc.','.$area_cyto.','.$mean_cyto.','.$ratio;
    $index++;
  }#if
}#while

close(CYTO);
##------------------------------------------


## Read mask files
$title = <INFILE>;              #read title line
@col = split ',',$title;
for(my $i=0; $i<@col; $i++){                   #analyze columns
  if($col[$i] eq "Slice"){$p_slice=$i;}
  if($col[$i] eq "X"){$p_x = $i;}
  if($col[$i] eq "Y"){$p_y = $i;}
  if($col[$i] eq "Area"){$p_area = $i;}
  if($col[$i] eq "Mean"){$p_mean = $i;}
}
chomp $title;

$title = $title.",Area_Nuc,Mean_Nuc,Area_Cyto,Mean_Cyto,Ratio\n";

my $idx =0;
my $slice =1;

while(<INFILE>){
  my $line=$_;
  chomp $line;
  my @col = split ',',$line;

  if($col[$p_slice]==$slice){
     $matrix[$slice][$idx][0]=$col[0];
     $matrix[$slice][$idx][1]=$col[$p_x];
     $matrix[$slice][$idx][2]=$col[$p_y];
     $matrix[$slice][$idx][3]=$col[$p_mean];
     $matrix[$slice][$idx][4]=$col[$p_area];
     $matrix[$slice][$idx][5]=$line.','.$measure{$col[0]};                           ##including cytosol measurement
     $p2s{$col[0]}=$slice;
     $p2i{$col[0]}=$idx;
     $idx++;
  }else{
     print "$site site $slice slice has $idx cells\n";
     $cc[$slice]=$idx;
     $slice=$col[$p_slice];
     $idx=0;
     $matrix[$slice][$idx][0]=$col[0];
     $matrix[$slice][$idx][1]=$col[$p_x];
     $matrix[$slice][$idx][2]=$col[$p_y];
     $matrix[$slice][$idx][3]=$col[$p_mean];
     $matrix[$slice][$idx][4]=$col[$p_area];
     $matrix[$slice][$idx][5]=$line.','.$measure{$col[0]};                           #including cytosol measurement;
     $p2s{$col[0]}=$slice;
     $p2i{$col[0]}=$idx;
     $idx++;
  }#endif(col[p_slice])
}#while

$cc[$slice]=$idx;
print "$site site $slice slice has $idx cells\n";

close(INFILE);

my @track;
my @trackleng;


for(my $i=1; $i < $slice; $i++){
  print "Processing site $site slice $i ";

  for(my $j=0; $j<$cc[$i];$j++){
    print ".";
    my $p = $matrix[$i][$j][0];
    my $x = $matrix[$i][$j][1];
    my $y = $matrix[$i][$j][2];

    my $dist =0;
    my $next =0;
    my $min_x=0;
    my $min_y=0;

    my $mindist =100;
    for(my $k=0; $k < $cc[$i+1];$k++){                           #search next slice & find shortest path
       my $np = $matrix[$i+1][$k][0];
       my $nx = $matrix[$i+1][$k][1];
       my $ny = $matrix[$i+1][$k][2];
       $dist = sqrt(($x-$nx)*($x-$nx)+($y-$ny)*($y-$ny));
       if($dist<$mindist){$mindist = $dist; $next = $np;$min_x=$nx; $min_y=$ny;}
    }#for(k)
    if($debug){print DEBUG "$p\t$x\t$y\t$next\t$min_x\t$min_y\t$mindist\n";}

    my $flag=0;

    if($mindist < $limit_d){                                                #found next point

      my $node = 't'.$p.'t';
      for (my $l=0;$l< @track;$l++){                #append existed track
         if ($track[$l] =~m/($node)/){
           $track[$l]=$track[$l].$next.'t';
           $trackleng[$l]++;
           $flag=1;
         }#if
      }#for ($l)

      if($flag == 0){                          #add new track

        if ($i > 3){

          ## backtracking
          my $dist =0;
          my $previous =0;
          my $mindist =100;

          for(my $k=0; $k < $cc[$i-1];$k++){                           #search next slice & find shortest path
             my $pp = $matrix[$i-1][$k][0];
             my $px = $matrix[$i-1][$k][1];
             my $py = $matrix[$i-1][$k][2];
             $dist = sqrt(($x-$px)*($x-$px)+($y-$py)*($y-$py));
             if($dist<$mindist){$mindist = $dist; $previous = $pp;}
          }#for(k)

          if($mindist <$limit_d){
            my $node ='t'.$previous.'t';
            #if($debug){print "$previous $p $next \n";}
            my $found =0;

            for(my $l=0;$l<@track && !$found;$l++){
              if($track[$l] =~m/($node)/){                 #found previous track
                 my $oldtrack = $track[$l];
                 my $oldtrackleng = $trackleng[$l];
                 my $pos = index ($oldtrack,$node);
                 $oldtrack = substr($oldtrack,0,$pos);

                 my $link = $oldtrack.'t'.$p.'t'.$next.'t';
                 push @track,$link;
                 push @trackleng,$oldtrackleng;
                 $found = 1;
              }# if(track[l]=($node)
            }#for(my l)

           if($found ==0){                                        #no previous track but with previous point
                 my $link = 't'.$previous.'t'.$p.'t'.$next.'t';
                 push @track,$link;
                 push @trackleng,3;
           }#ifelse(track[i]=~m)

          }#if($midist)                                        #foudn previous point
        ##end backtracking

        }else{                                                 #first two slices
          my $link='t'.$p.'t'.$next.'t';
          push @track,$link;
          push @trackleng,2;
        }#if($i>2)

      }#if(flag = 0)                                           #add new track

    }#if(mindist)  --> found next point

  }#for (j) --> scan through each cell
  print"\n";
}#for(i)    --> scan through each slice

print OUTFILE "no\tTrack_length\tTrack\n";
open MEAN,">output.mean-s$site.csv";
open RATIO,">output.ratio-s$site.csv";



for(my $i=0;$i<@track;$i++){

  my $goodtrack =4;
  my $line = $track[$i];
  my @col = split(/t/,$line);
  my $average = 0;
  my $number =0;
  foreach my $key(@col){if($key){$average = $average + $matrix[$p2s{$key}][$p2i{$key}][4];$number++;}}
  $average = $average/$number;
  
  foreach my $key(@col){if($key){if($matrix[$p2s{$key}][$p2i{$key}][4] < 0.5*$average || $matrix[$p2s{$key}][$p2i{$key}][4] > 1.5*$average ){$goodtrack--;}}}

  if($goodtrack>0){print OUTFILE "$i\t$trackleng[$i]\t$track[$i]\n";}

  if($trackleng[$i]>$slice-2 && $goodtrack>0){
    print TRACK "****** Output track $i ** Length = $trackleng[$i]*******\n";
    print TRACK $title;
    my $line = $track[$i];
    my @col = split(/t/,$line);
    foreach my $key (@col){if($key){print TRACK "$matrix[$p2s{$key}][$p2i{$key}][5]\n";print MEAN "$mean_n{$key},"; print RATIO "$ratio_nc{$key},";}}#foreach
    print MEAN "\n";print RATIO "\n";
  }#if(trakcleng)
}#for

close(OUTFILE);
close(TRACK);
close(DEBUG);
close(MEAN);
close(RATIO);

}

