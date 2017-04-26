#!/usr/bin/perl

use strict;
use warnings;

my ($open_exnode_file, $filename, $force, $force2, );

# knoten auf der rechten Seite
# my @biceps_bottom_nodes=(2,4,6,8,17,18,23,24,27);
# knoten Oben
#my @biceps_bottom_nodes = (5..8,11,12,14,16,18);

# nodes right
# my @biceps_bottom_nodes=(5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100,105,110,115,120,125);
# my @biceps_bottom_nodes=(3,6,9,12,15,18,21,24,27);
# my @biceps_bottom_nodes=(5,10,15,20,25,30,35,40,45);
# my @biceps_bottom_nodes=(25,50,75,100,125,150,175,200,225);
# my @biceps_bottom_nodes=(9,18,27,36,45,54,63,72,81,90,99,108,117,126,135,144,153,162,171,180,189,198,207,216,225);
# nodes left
# my @biceps_bottom_nodes=(1,6,11,16,21,26,31,36,41,46,51,56,61,66,71,76,81,86,91,96,101,106,111,116,121);

#my @biceps_bottom_nodes=(15,30,45,60,75,90,105,120,135);
#my @biceps_bottom_nodes=(17,34,51,68,85,102,119,136,151);
#my @biceps_bottom_nodes=(29,58,87,116,145,174,203,232,261);
my @biceps_bottom_nodes=(25,50,75,100,125,150,175,200,225);

print "Evaluating nodes: @biceps_bottom_nodes \n";
print "#######################################\n";

my @disps = (1..2000);

open WRITE_RESULTS, ">reac_forces.dat" or die "$!";
printf WRITE_RESULTS "# time   force\n";

  foreach my $disp (@disps)
  {
    $filename = sprintf("MainTime_1_%d", $disp);
    $open_exnode_file = "./$filename.part0.exnode";
    print "Reading: $open_exnode_file \n" ;
    $force = call_biceps_bottom(@biceps_bottom_nodes);

#    $open_exnode_file = "./$filename.part1.exnode";
#    print "Reading: $open_exnode_file \n" ;
#    $force2 = call_biceps_bottom(@biceps_bottom_nodes);
#    $force = $force + $force2;
#    
#    $open_exnode_file = "./$filename.part2.exnode";
#    print "Reading: $open_exnode_file \n" ;
#    $force2 = call_biceps_bottom(@biceps_bottom_nodes);
#    $force = $force + $force2;

#    $open_exnode_file = "./$filename.part3.exnode";
#    print "Reading: $open_exnode_file \n" ;
#    $force2 = call_biceps_bottom(@biceps_bottom_nodes);
#    $force = $force + $force2;

    printf WRITE_RESULTS "%3.7f  %2.6e \n", ($disp), $force;
  }

close WRITE;








system("rm all_fields field_* reac_forces");

sub call_biceps_bottom{

    my @biceps_bottom_nodes=@_;
	
    my $write_file="all_fields";
    my $line_no=-1;
    my $node_no;
    my $first=0;
    my $write4="field_4";
    my $write5="field_5";
    my $write6="field_6";
    my $sum_forces_x=0;
    my $sum_forces_y=0;
    my $sum_forces_z=0;
    my $magnitude=0;

# im prinzip nur die daten zu knoten aus den exnode files raussuchen und in die reaktionskräfte incl. ableitugen in all_fields sammeln 
    open READ, "$open_exnode_file" or die "$!"; #wurde von start_cmiss übergeben
    open WRITE, ">$write_file" or die "$!"; #all_fields
    while(my $line=<READ>){
      #m=The match operator; suche nach Node, \s - string \d - number      
	    if($line=~m/^ Node:[\s]+([\d]+)/) #wenn das zutrifft dann tue:
	      {$line_no=0;$node_no=$1;$first=1}
	    if(($line_no==11||$line_no==12||$line_no==13) && ($first==1))
	    {
	      print WRITE "#Node_x: ",$node_no," " if($line_no==11);
	      print WRITE "#Node_y: ",$node_no," " if($line_no==12);
	      print WRITE "#Node_z: ",$node_no," " if($line_no==13);
	      print WRITE $line;
	    }
	    $line_no++;
    }
    close READ;
    close WRITE;

# lese alle daten von all_fields und schreibe nur die von biceps_bottom_nodes in reac_forces
    my $open_file="all_fields";
    $write_file="reac_forces";
    open READ,   "$open_file" or die "$!";
    open WRITE, ">$write_file" or die "$!";
    open WRITE4, ">$write4" or die "$!";
    open WRITE5, ">$write5" or die "$!";
    open WRITE6, ">$write6" or die "$!";
    while(my $line=<READ>){
	    if($line=~m/^\#Node_([\w]): ([\d]+)/){
	      foreach my $node (@biceps_bottom_nodes){ 
		      if($2==$node && $1 eq "x"){
		        print WRITE4 "$line";
		        print WRITE "$line";
		      }
		      elsif($2==$node && $1 eq "y"){
		        print WRITE5 "$line";
		        print WRITE "$line";
		      }
		      elsif($2==$node && $1 eq "z"){
		        print WRITE6 "$line";
		        print WRITE "$line";
		      }
	      }
	    }
    }
    close READ;
    close WRITE;
    close WRITE4;
    close WRITE5;
    close WRITE6;

# lese field4,5,6 und addiere die reaktionskräfte in x,y & z
    open READ4, "$write4" or die "$!";
    open READ5, "$write5" or die "$!";
    open READ6, "$write6" or die "$!";

    while(my $line=<READ4>){
	    if ($line=~m/([+-]?[\d]+\.[\d]+E[+-]?[\d]+).*/){
	      $sum_forces_x+=$1;
	    }
    }

    while(my $line=<READ5>){
	    if ($line=~m/([+-]?[\d]+\.[\d]+E[+-]?[\d]+).*/){
	      $sum_forces_y+=$1;
	    }
    }

    while(my $line=<READ6>){
	    if ($line=~m/([+-]?[\d]+\.[\d]+E[+-]?[\d]+).*/){
	      $sum_forces_z+=$1;
	    }
    }
    close READ4;
    close READ5;
    close READ6;
    
#    $magnitude=sqrt($sum_forces_x**2+$sum_forces_y**2+$sum_forces_z**2);
    $magnitude=$sum_forces_x;

    return $magnitude;

} #ende subroutine
    
