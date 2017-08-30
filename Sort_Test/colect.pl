#!/usr/bin/perl
use strict;
use warnings;
use utf8;


my ($i,$j,$k,$l,$str,@time);
my @types = (['-DMerge','Merge'],['-DRadix','Radix 8'],['-DHybrid','Hybrid'],['','qsort']);

for($i=0;$i<=22;$i++)
{
  $j = 2 ** $i;
  printf "%d\n" , $j;
  for($k=0;$k<=$#types;$k++)
  {
    $str = "";
    open(IN,"gcc sort_RN.c -O2 $types[$k][0] -DTimes=$j |");
    while(<IN>)
    {
      $str .= $_;
    }
    close(IN);
    $str = "";
    open(IN,"echo '' >  time.txt |");
    while(<IN>)
    {
      $str .= $_;
    }
    close(IN);
    @time = ();
    for($l=0;$l<5;$l++)
    {
      $str = "";
      open(IN,"time -o time.txt -a -p ./a.out |");
      while(<IN>)
      {
        $str .= $_;
      }
      close(IN);
    }
    $str = "";
    open(IN,"cat time.txt | grep user | sed 's|user ||g' | sort -n |");
    while(<IN>)
    {
      $str .= $_;
    }
    close(IN);
    $str =~ s"[a-zA-Z]+" "sg;
    $str =~ s"\n+" "sg;
    $str =~ s"\s+" "sg;
    $str =~ s"^ *""sg;
    $str =~ s" *$""sg;
    @time = split(/ /,$str);
    if($time[$#time/2] + 0.0 > 50.0)
    {
      printf "%-10s %s\n",$types[$k][1] , $time[$#time/2];
      splice(@types,$k,1);
      $k--;
      next;
    }
    elsif($time[$#time/2] + 0.0 > 10.0)
    {
      printf "%-10s %s\n",$types[$k][1] , $time[$#time/2];
      next;
    }
    for($l=5;$l<15;$l++)
    {
      $str = "";
      open(IN,"time -o time.txt -a -p ./a.out |");
      while(<IN>)
      {
        $str .= $_;
      }
      close(IN);
    }
    $str = "";
    open(IN,"cat time.txt | grep user | sed 's|user ||g' | sort -n |");
    while(<IN>)
    {
      $str .= $_;
    }
    close(IN);
    $str =~ s"[a-zA-Z]+" "sg;
    $str =~ s"\n+" "sg;
    $str =~ s"\s+" "sg;
    $str =~ s"^ *""sg;
    $str =~ s" *$""sg;
    @time = split(/ /,$str);
    if($time[$#time/2] + 0.0 > 5.0)
    {
      printf "%-10s %s\n",$types[$k][1] , $time[$#time/2];
      next;
    }
    for($l=15;$l<35;$l++)
    {
      $str = "";
      open(IN,"time -o time.txt -a -p ./a.out |");
      while(<IN>)
      {
        $str .= $_;
      }
      close(IN);
    }
    $str = "";
    open(IN,"cat time.txt | grep user | sed 's|user ||g' | sort -n |");
    while(<IN>)
    {
      $str .= $_;
    }
    close(IN);
    $str =~ s"[a-zA-Z]+" "sg;
    $str =~ s"\n+" "sg;
    $str =~ s"\s+" "sg;
    $str =~ s"^ *""sg;
    $str =~ s" *$""sg;
    @time = split(/ /,$str);
    printf "%-10s %s\n",$types[$k][1] , $time[$#time/2];
  }
  if($#types == -1)
  {
    last;
  }
}
