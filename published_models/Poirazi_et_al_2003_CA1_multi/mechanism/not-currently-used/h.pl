#!/usr/bin/perl

use PDL;

$vhalf=-90;
$K=8.5;

for $v (-70 .. 40) {
  
  $val= 1 - (1/(1+exp(($vhalf - $v)/$K)));
  push(@ninf, $val);
  push(@x,$v);

}

$ninf=pdl(@ninf);

line $ninf;

while (1==1) {;}

