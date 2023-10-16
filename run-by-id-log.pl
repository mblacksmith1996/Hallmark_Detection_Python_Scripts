#!/usr/bin/perl -w

use Fcntl qw(:flock); 

# Jeff Kidd
# simple perl script, for running SGE array jobs
# Hope it works!

if(@ARGV != 3)
{
  print "Three args, commands file, log file for completed commands, and the job ID to execute\n";
  exit(1)
}
($cf,$logf, $id)=@ARGV;



open IN, "< $cf " or die "cannot open in $cf files\n";
$lnum = 0;
while ($l = <IN>)
{
  chomp($l);
  $lnum++;
  if($lnum == $id)
  {
    system $l;
    open OUT, ">> $logf" or die "NO LOG FILE";
#    flock OUT, LOCK_EX;
    
    if ( $? == -1 )
    {
       print OUT "FAILED: $? $l\n";
    }
    else
    {
       print OUT "$l\n";

    }
    # flock OUT, $UNLOCK; #flock released when closed file
    close OUT;
  }
}
close IN;
