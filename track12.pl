## tracking script (Track11.pl)
## Jerry Lin 2014/09/01

##
## for tracking two channel image
## nucleus and cytoplasmic ratio
## Usage: Perl track10.pl <begin row> <end row> <begin col> <end col>
## Input files: Measure-mask-sXX.csv & Measure-Cytosol-sXX.csv
## Output files: Tracks-sXX.tsv & Tracking-sXX.csv
## additional output files: debug-sXX.log (XY position)  output.mean-sXX.csv (mean_nuc only output) output.ratio-sXX.csv (ratio only output)
## 06b update (05/19/2014)--> including area threshold to reject wrong tracks
## 07 update (06/30/2014) --> add area output files
## 08 update (07/11/2014) --> add slice/frame & position output files
## 09 update (2014/08/29) --> modify for operetta analysis
## 10 update (2014/08/31) --> dealing with skipping frames & add new output files (output.point & output.total)
## 11 update (2014/09/01) --> removing duplicates
## 12 update (2014/09/02) --> modify for processing both NIC and Operetta images
##


use warnings;
use strict;

my $rows = 2;
my $rowe = 7;
my $cols = 2;
my $cole = 11;

my @rowarray = ("000","001","002","003","004","005","006","007","008");
my @colarray = ("000","001","002","003","004","005","006","007","008","009","010","011","012");


if($ARGV[0]){$rows=$ARGV[0];}
if($ARGV[1]){$rowe=$ARGV[1];}
if($ARGV[2]){$cols=$ARGV[2];}
if($ARGV[3]){$cole=$ARGV[3];}

## if set cols or cole = 99, then activate NIC mode

for(my $row=$rows;$row<$rowe+1;$row++)
{
for(my $col=$cols;$col<$cole+1;$col++)
{

my $site;

if ($cols==99 || $cole==99){$site = "s$row";}else{$site = $rowarray[$row].$colarray[$col];}

print "Now processing $site\n";

##=======start tracking routine=========
my $filename='measure-mask-'.$site.'.csv';
my $output = 'tracking-'.$site.'.tsv';
my $trackf = 'tracks-'.$site.'.csv';
my $cytof ='measure-cytosol-'.$site.'.csv';



my $debugfile = 'debug-'.$site.'.log';

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


##define maximum displacement
my $limit_d=25;

##define maximun skipped frames
my $max_skip=10;

##define maximun error for mistracking
my $max_error=8;

## switch for insertion mode
my $ins_mode = 1;

## switch for removing duplicates
my $dup_off = 1;

## define replicate threshold
my $rep_thre = 0.9;


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
my %area_n;
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
    $area_n{$index}=$area_nuc;

    $area_cyto =$area_cyto-$area_nuc;

    my $ratio = $mean_nuc/$mean_cyto;
    $ratio_nc{$index}=$ratio;

    $measure{$index}=$area_nuc.','.$mean_nuc.','.$area_cyto.','.$mean_cyto.','.$ratio;
    $index++;
  }#if
}#while

close(CYTO);
##------------------------------------------

my %area;
my %slice;
my %maskmean;
my %pos;

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
  
  if($col[0]){
    $area{$col[0]}=$col[$p_area];                                #area hash
    $slice{$col[0]}=$col[$p_slice];                              #slice/frame hash
    $maskmean{$col[0]}=$col[$p_mean];                            #mask_mean hash
    $pos{$col[0]}=$col[$p_x].'-'.$col[$p_y];
  }

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
open MEAN,">output.mean-$site.csv";
open RATIO,">output.ratio-$site.csv";
open AREA,">output.area-$site.csv";
open POSI,">output.pos-$site.csv";
open SLICE,">output.slice-$site.csv";
open MASK,">output.mask-$site.csv";
open POINT,">output.point-$site.csv";
open TOTAL,">output.total-$site.csv";

my %total_m;
my %exist;

for(my $i=0;$i<@track;$i++){

  my $goodtrack =$max_error;
  my $line = $track[$i];
  my @col = split(/t/,$line);
  my $average = 0;
  my $number =0;
  my $avgarea = 0;
  my $avgtotal = 0;

  foreach my $key(@col){
    if($key){
      $average = $average + $maskmean{$key};
      $avgarea = $avgarea + $area{$key};
      $avgtotal = $avgtotal + $maskmean{$key}*$area{$key};
      $total_m{$key} = $maskmean{$key}*$area{$key};
      $number++;
      }#endif($key)
    }#foreach

  $average = $average/$number;
  $avgarea = $avgarea/$number;
  $avgtotal = $avgtotal/$number;

  foreach my $key(@col){                                       ## check unusual fluatuation on mask mean & area
    if($key){
      ## my $total = $maskmean{$key}*$area{$key};
      if($maskmean{$key} < 0.5*$average || $maskmean{$key} > 1.5*$average ){
        if($total_m{$key} < 0.75*$avgtotal || $total_m{$key} > 1.25*$avgtotal){$goodtrack--;}
        }
      if($area{$key} < 0.5*$avgarea || $area{$key} > 1.5*$avgarea){
        if($total_m{$key} < 0.75*$avgtotal || $total_m{$key} > 1.25*$avgtotal){$goodtrack--;}
        }
      if($total_m{$key} < 0.5*$avgtotal || $total_m{$key} > 1.5*$avgtotal){$goodtrack--;}
      }#endif($key)
    }#foreach

  if($goodtrack>0){print OUTFILE "$i\t$trackleng[$i]\t$track[$i]\n";}

  if($trackleng[$i]>$slice-$max_skip && $goodtrack>0){
    print TRACK "****** Output track $i ** Length = $trackleng[$i]*******\n";
    print TRACK $title;
    my $line = $track[$i];
    my @col = split(/t/,$line);
    my $number =0;
    my $p_mean;
    my $p_ratio;
    my $p_area;
    my $p_maskmean;
    my $p_slice;
    my $check = 1;
    my $leng = 0;
    my $rep =0;

    if($dup_off){
      foreach my $key (@col){$leng++;if($exist{$key}){$rep++;}}  ## read-through and check repeats
      if($rep/$leng > $rep_thre){$check =0;}
    }

   if($check){
    foreach my $key (@col){
      $number++;

      if($key){
        $exist{$key}++;
        if($number>2 && $ins_mode){
          if($slice{$key}-$p_slice > 1){
            my $ins_mean = ($p_mean+$mean_n{$key})/2;
            print MEAN "$ins_mean,";
            my $ins_ratio = ($p_ratio+$ratio_nc{$key})/2;
            print RATIO "$ins_ratio,";
            my $ins_area = ($p_area+$area{$key})/2;
            print AREA "$ins_area,";
            my $ins_mask = ($p_maskmean+$maskmean{$key})/2;
            print MASK "$ins_mask,";
          }
        }#endif($ins_mode)


        print TRACK "$matrix[$p2s{$key}][$p2i{$key}][5]\n";
        print MEAN "$mean_n{$key},";$p_mean = $mean_n{$key};
        print RATIO "$ratio_nc{$key},";$p_ratio = $ratio_nc{$key};
        print AREA "$area{$key},";$p_area = $area{$key};
        print POSI "$pos{$key},";
        print TOTAL "$total_m{$key},";
        print SLICE "$slice{$key},";$p_slice=$slice{$key};
        print MASK "$maskmean{$key},";$p_maskmean=$maskmean{$key};
        print POINT "$key,";
      }#endif
    }#foreach
    print MEAN "\n";print RATIO "\n";print AREA "\n";print POSI "\n";print SLICE "\n";print MASK "\n";print POINT "\n";print TOTAL "\n";
   }#if(check)
  }#if(trakcleng)
}#for

close(OUTFILE);
close(TRACK);
close(DEBUG);
close(MEAN);
close(RATIO);
close(AREA);
close(POSI);
close(SLICE);
close(MASK);
close(POINT);
close(TOTAL);

}#for (col)
}#for (row)
