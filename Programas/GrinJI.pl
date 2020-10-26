#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use List::Util;
use POSIX qw(DBL_MAX);
use PDL;
use PDL::NiceSlice;
use PDL::FFTW3;
use PDL::Complex;
use PDL::Graphics::Gnuplot;
use constant {PI=>4*atan2(1,1), LM=>log(DBL_MAX)};
use IO::Prompter;
#use Chart::Gnuplot;
#use GD::Graph;

#Calculate current density in an electrolytic cell with a 'point'
#electrode, to build GRIN structures.
#Notes GX
#Eq. Gx33

#Normalize to I and to the area.
my $a; #size of cell along x
my $b; #size of cell along y
my $C; #size of cell along z (water height)
my ($X0, $Y0, $Z0); #positioin of electrode
my $Z; #height to evaluate current.
my $N; #2N+1=number of points along each direction
my $I; #corriente
# $M; #2*$M+1=size of mask in reciprocal space
Getopt::Long::Configure('no_ignore_case');
GetOptions(
    'a=f'=>\$a,
    'b=f'=>\$b,
    'C=f'=>\$C,
    'X0=f'=>\$X0,
    'Y0=f'=>\$Y0,
    'Z0=f'=>\$Z0,
    'Z=f'=>\$Z,
    'N=f'=>\$N,
    'I=f'=>\$I,
    #'M=f'=>\$M
    )
    or usage();
usage("Undefined parameters") unless List::Util::all {defined $_}
($a, $b, $C, $X0, $Y0, $Z0, $Z, $I, $N); #, $M);
my $A=$a*$b; #area celda
my $u=sqrt($A); #unidades
my $a1=$a/$u; # largo adimensional
my $c=$C/$u;  # altura de liquido adimensional
my $x0=$X0/$u;
my $y0=$Y0/$u;
my $z0=$Z0/$u;
my $z=$Z/$u;
my $b1=1/$a1; #So area=1, this is our unit of area, and from here, distance
#Real space
my $x=zeroes(2*$N+1, 2*$N+1)->xvals*2*$a1/(2*$N+1);
my $y=zeroes(2*$N+1, 2*$N+1)->yvals*2*$b1/(2*$N+1);
#Reciprocal lattice shifted for Fourier transforms.
my $G=((zeroes(2*$N+1, 2*$N+1)->ndcoords-$N)*PI/pdl($a1, $b1))
    ->mv(1,0)->rotate(-$N)->mv(0,1)->mv(2,0)->rotate(-$N)->mv(0,2);
my $Gn=($G*$G)->sumover->sqrt; #norm;
my $sinh2Gc=sinh(mymin(2*$Gn*$c, LM));
$sinh2Gc->((0),(0)).=1; #to avoid division by 0 down
#my $jG=cosh(mymin($Gn*$z,LM))*
#    (exp(-mymin($Gn*$z0,LM))
#     -exp(-mymin(3*$Gn*$c,LM))*cosh(mymin($Gn*$c,LM))*sinh(mymin($Gn*$z0,LM))
#     /$sinh2Gc)
#    *cos($G->((0))*$x0)*cos($G->((1))*$y0); #*$mask;
#$jG->((0),(0)).=1; #Special value for G=0;
#my $jR=fft2($jG->r2C->real)->complex->re;
#my $j1=fft2($jG->r2C->real)->complex->re;
#my $jR=$j1/$A;
my $jG1=2*cosh(mymin($Gn*$z,LM))*
    sinh(mymin($Gn*$c, LM))*cosh(mymin($Gn*($c-$z0), LM))/$sinh2Gc
    *cos($G->((0))*$x0)*cos($G->((1))*$y0); #*$mask;
$jG1->((0),(0)).=1; #Special value for G=0;
my $jR1=fft2($jG1->r2C->real)->complex->re;
#my $jR=$jR1*$I/$A; #Se calcula j para cualquier I
my $X=$x*$u;
my $Y=$y*$u;

#figura vista desde arriba
my $ny=$N/2; #0 is y=0, N is y=b, 2N=2b
my $gw=PDL::Graphics::Gnuplot->new();
$gw->multiplot;
$gw->plot({
    #title=>["{/Times*2 a=$a cm, y=$b cm, C=$C, R0=($X0,$Y0,$Z0)cm, Z=$Z}",textcolor=>'"blue"'], #
    xrange=>[0,$a], yrange=>[0,$b],
    #view=>['equal','xy'], #parece no funcionar
    #size=>['ratio',-1],
    #mochan
    #size=>'1,.5',
    #origin=>[0,.5],
    #cris
    #size=>'1,.5',
    #origin=>[0,.5],
    justify=>1,
    label=>[1, "(a)", at=>"screen .8,.95", textcolor=>'"white"'],
    xlabel=>["x(cm)"] ,ylabel=>["y (cm)"], cblabel=>["j_⟂/I(1/cm²)"]
	  }
	  ,with=>'image', $X, $Y, $jR1);
# figura cortada en el  eje x ,
$gw->plot( {
    #title=>["{/Times*2 x=$a cm, y=$Y0 cm, C=$C cm, R_0=($X0,$Y0,$Z0)cm, Z=$Z}",textcolor=>'"black"'],
    xrange=>[0,$a], xlabel=>["x(cm)"], yrange=>[],
    ylabel=>["j_⟂/I(1/cm²)"], size=>'.875,.5', origin=>[0,0],
    label=>[2, "(b)", at=>"screen .8,.45"],
},
	   {with=>'lines'}, $X(:,($ny)),$jR1(:,($ny)));
$gw->end_multi;

prompt "Ya casi!", -single, -void;

sub usage {
    my $s=shift//"";
    print "$s\n";
    print <<'FIN';
#Calculate current density in an electrolytic cell with a 'point'
#electrode, to build GRIN structures.
#Notes GX
#Eq. Gx33

#Normalize to I and to the area.
my $a; #size of cell along x
my $b; #size of cell along y
my $C; #size of cell along z (water height)
my ($X0, $Y0, $Z0); #positioin of electrode
my $z; #height to evaluate current.
my $N; #2N+1=number of points along each direction
my $I; #corriente
#my $M; #2*$M+1=size of mask in reciprocal space
GetOptions(
    'a=f'=>\$a,
    'b=f'=>\$b,
    'C=f'=>\$C,
    'X0=f'=>\$X0,
    'Y0=f'=>\$Y0,
    'Z0=f'=>\$Z0,
    'z=f'=>\$z,
    'N=f'=>\$N,
    'I=f'=>\$I,
    #'M=f'=>\$M
    )
    or usage();
usage("Undefined parameters") unless List::Util::all {defined $_}
($a1, $b1, $C, $X0, $Y0, $Z0, $z, $I, $N); #, $M);
my $b1=1/$a1; #So area=1, this is our unit of area, and from here, distance

FIN
    exit(1);
}

sub mymax {
    my ($a1, $b1, $r)=@_;
    $r//=null;
    mymaxaux($a1, $b1, $r);
    return $r;
}
sub mymin {
    my ($a1, $b1, $r)=@_;
    $r//=null;
    myminaux($a1, $b1, $r);
    return $r;
}

BEGIN {
    thread_define 'mymaxaux(a1();b1();[o]r())', over {
	my $a1=shift;
	my $b1=shift;
	my $r=shift;
	$r.=$a1>=$b1?$a1:$b1;
    };
    thread_define 'myminaux(a1();b1();[o] r())', over {
	my $a1=shift;
	my $b1=shift;
	my $r=shift;
	$r.=$a1<=$b1?$a1:$b1;
    };
}

#perl GrinJI.pl --a=2.03 --b=1.48 --C=1 --X0=0.2 --Y0=0.74 --Z0=0.9 --Z=0 --I=5 --N=100
