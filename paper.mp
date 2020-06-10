g_use_svg = 0; %0 for postscript, 1 for svg
prologues:=3;
u=1cm;

outputtemplate := "%j-%c.eps";
if g_use_svg > 0:
	outputtemplate := "%j-%c.svg";
	outputformat := "svg";
fi

def refleja(expr p)=
  drawarrow p cutbefore (fullcircle scaled 10pt shifted point 0 of p) cutafter
    (fullcircle scaled 10pt shifted point 0 of reverse p) dashed evenly;
  draw ((0,-.1u)--(0,.1u)) rotated angle direction .5*length p of p
    shifted point .5*length p of p withcolor blue withpen pencircle
    scaled 3pt;
enddef;
interim ahlength:=1.5*ahlength;
a=2u;
b=1u;
c=1.5u;
N=1;
NM1=N+1;
NM12=N+.5;
NM112=N+1.5;

beginfig(1)
  z0=(.5u,.3u);
  z1=(-x0,y0);
  z2=(x0,-y0);
  z3=(-x0,-y0);
  for i:=-2*N-1 upto 2*NM1+1:
    draw (i*a,-2*NM12*b)--(i*a,2*NM112*b);
    draw (-2*NM12*a,i*b)--(2*NM112*a,i*b);
  endfor
  for i:=-N upto NM1:
    for j:= -N upto NM1:
      draw z0 shifted(2*i*a, 2*j*b) withpen pencircle scaled 5pt
	withcolor blue;
      draw z1 shifted(2*i*a, 2*j*b) withpen pencircle scaled 5pt
	withcolor blue;
      draw z2 shifted(2*i*a, 2*j*b) withpen pencircle scaled 5pt
	withcolor blue;
      draw z3 shifted(2*i*a, 2*j*b) withpen pencircle scaled 5pt
	withcolor blue;
    endfor
  endfor
  fill unitsquare xscaled a yscaled b withcolor .7white+.3blue;
  draw z0 withpen pencircle scaled 5pt withcolor red;
  refleja(z0{up}..z1{down});
  refleja(z0{right}..z2{left});
  refleja(z0{dir 25}..z1 shifted(2a,0){dir -25});
  refleja(z0{dir 45}..z2 shifted(0,2b){dir 135});
  refleja(z1{left}..z3{right});
  refleja(z1{dir 75}..z0 shifted(2*a,0){dir -75});
  refleja(z2{down}..z3{up});
  refleja(z2{dir -45}..z3 shifted(2*a,0){dir 45});
  refleja(z3{dir -80}..z2 shifted(2*a,0){dir 80});
  refleja(z2 shifted(2*a,0){dir 45}..z0 shifted(2*a,2*b){dir 135});
  refleja(z0 shifted(2*a,0){dir 100}..z1 shifted(-2a,0){dir -100});
  draw unitsquare shifted (-.5,-.5) xscaled 2a yscaled 2b shifted
    (4a,2b) dashed evenly scaled 2 withpen pencircle scaled 2pt;
  %interim ahangle:=180;
  drawdblarrow ((0,0)--(a,0)) shifted (0,-.2u) shifted (0,-2NM12*b);
  label.bot(btex $a$ etex scaled 2, (.5a,0) shifted (0,-.2u) shifted (0,-2NM12*b));
  drawdblarrow ((0,0)--(0,b)) shifted (-.2u,0) shifted (-2NM12*a,0);
  label.lft(btex $b$ etex scaled 2, (0,.5b) shifted (-.2u,0) shifted (-2NM12*a,0));
  drawdblarrow ((-a,0)--(a,0)) shifted (4a,0) shifted (0,-.2u) shifted (0,-2NM12*b);
  label.bot(btex $2a$ etex scaled 2, (4a,0) shifted (0,-.2u) shifted (0,-2NM12*b));
  drawdblarrow ((0,-b)--(0,b)) shifted (0,2b) shifted (-.2u,0) shifted (-2NM12*a,0);
  label.lft(btex $2b$ etex scaled 2, (0,2b) shifted (-.2u,0) shifted (-2NM12*a,0));
endfig;

beginfig(2)
  z0=(.5u,1.u);
  z1=(-x0,y0);
  z2=(x0,-y0);
  z3=(-x0,-y0);
  fill unitsquare xscaled a yscaled c withcolor .7white+.3blue;
  for i:=-2*N-1 upto 2*NM1+1:
    draw (i*a,-2*NM12*c)--(i*a,2*NM112*c);
    draw (-2*NM12*a,i*c)--(2*NM112*a,i*c);
  endfor
  for i:=-N upto NM1:
    for j:= -N upto NM1:
      if j mod 2=0:
	filldraw fullcircle scaled 4pt shifted z0 shifted(2*i*a, 2*j*c) withpen pencircle scaled 1.5pt
	  withcolor green;
	filldraw fullcircle scaled 4pt shifted z1 shifted(2*i*a, 2*j*c) withpen pencircle scaled 1.5pt
	  withcolor green;
	draw fullcircle scaled 4pt shifted z2 shifted(2*i*a, 2*j*c) withpen pencircle scaled 1.5pt
	  withcolor green;
	draw fullcircle scaled 4pt shifted z3 shifted(2*i*a, 2*j*c) withpen pencircle scaled 1.5pt
	  withcolor green;
      fi;
      if j mod 2=1:
	filldraw fullcircle scaled 4pt shifted z2 shifted(2*i*a, 2*j*c) withpen pencircle scaled 1.5pt
	  withcolor green;
	filldraw fullcircle scaled 4pt shifted z3 shifted(2*i*a, 2*j*c) withpen pencircle scaled 1.5pt
	  withcolor green;
	draw fullcircle scaled 4pt shifted z0 shifted(2*i*a, 2*j*c) withpen pencircle scaled 1.5pt
	  withcolor green;
	draw fullcircle scaled 4pt shifted z1 shifted(2*i*a, 2*j*c) withpen pencircle scaled 1.5pt
	  withcolor green;
      fi;
    endfor
  endfor
  for i:=-N upto NM1:
    draw z0 shifted(2*i*a, 0) withpen pencircle scaled 5pt
      withcolor blue;
    draw z1 shifted(2*i*a, 0) withpen pencircle scaled 5pt
      withcolor blue;
  endfor
  filldraw fullcircle scaled 4pt shifted z0 withpen pencircle scaled 1.5pt withcolor red;
  refleja(z0{dir -135}..z1{dir 135});
  refleja(z0{dir -45}..z2{dir -135});
  refleja(z0{dir 135}..z2 shifted(0,2c){dir 45});
  refleja(z2 shifted(0,2c){dir -45}..z0 shifted(0,-2c){dir -135});
  refleja(z2{dir 150}..z0 shifted(0,2c){dir 30});
  refleja(z0 shifted(0,-2c){dir 150}..z2 shifted(0,4c){dir 30});
  draw unitsquare shifted (-.5,-.5) xscaled 2a yscaled 4c shifted
    (4a,c) dashed evenly scaled 2 withpen pencircle scaled 2pt;
  %draw unitsquare shifted (-.5,-.5) xscaled 2a yscaled 2c shifted
  %  (-2a,0) dashed withdots scaled 2 withpen pencircle scaled 4pt;
  %draw subpath (1,4) of unitsquare shifted (-.5,-.5) xscaled 2a yscaled 2c shifted
  %  (-2a,2c) dashed withdots scaled 2 withpen pencircle scaled 4pt;
  %interim ahangle:=180;
  drawdblarrow ((0,0)--(a,0)) shifted (0,-.2u) shifted (0,-2NM12*c);
  label.bot(btex $a$ etex scaled 2, (.5a,0) shifted (0,-.2u) shifted (0,-2NM12*c));
  drawdblarrow ((0,0)--(0,c)) shifted (-.2u,0) shifted (-2NM12*a,0);
  label.lft(btex $c$ etex scaled 2, (0,.5c) shifted (-.2u,0) shifted (-2NM12*a,0));
  drawdblarrow ((-a,0)--(a,0)) shifted (4a,0) shifted (0,-.2u) shifted (0,-2NM12*c);
  label.bot(btex $2a$ etex scaled 2, (4a,0) shifted (0,-.2u) shifted (0,-2NM12*c));
  drawdblarrow ((0,-2c)--(0,2c)) shifted (0,c) shifted (.2u,0) shifted (2NM112*a,0);
  label.rt(btex $4c$ etex scaled 2, (0,c) shifted (.2u,0) shifted (2NM112*a,0));
  draw ((-a,y0)--(a,y0)) shifted (4a,0) dashed evenly;
  draw ((-a,-y0)--(a,-y0)) shifted (4a,0) dashed evenly;
  for i:=-2 upto 2:
    drawarrow ((i*a/3,.9y0)--(i*a/3,-.9y0)) shifted (4a,0);
  endfor
  draw ((-a,2c-y0)--(a,2c-y0)) shifted (4a,0) dashed evenly;
  draw ((-a,2c+y0)--(a,2c+y0)) shifted (4a,0) dashed evenly;
  for i:=-2 upto 2:
    drawarrow ((i*a/3,2c-.9y0)--(i*a/3,2c+.9y0)) shifted (4a,0);
  endfor
endfig;

def poro(expr n, size)=
  begingroup (for i=0 upto n-1:
      dir(i*360/(n)) scaled (.7+uniformdeviate(.3)) scaled (size)--
    endfor cycle)
  endgroup
enddef;

beginfig(3)
  picture pic[];
  fill unitsquare xscaled 30u yscaled 10u withcolor .9white;
  fill fullcircle scaled 6u shifted (5u,5u) withcolor .6white;
  fill fullcircle scaled 6u shifted (15u,5u) withcolor .6white;
  fill fullcircle scaled 6u shifted (25u,5u) withcolor .6white;
  for i=0 upto 30:
    begingroup
      save n;
      numeric n[];
      n0=10+i;
      n1=10/n0;
      n2=n1*.5u;
      n3=sqrt(n1);
      for j=0 upto n0:
	fill poro(20,n2) shifted (i*u+uniformdeviate(.7u), j*n1*u+uniformdeviate(.5*n1*u));
      endfor;
    endgroup;
  endfor;
  pic0=thelabel(btex A etex scaled 7, (5u,5u));
  draw pic0 withcolor white;
  pic1=thelabel(btex B etex scaled 7, (15u,5u));
  draw pic1 withcolor white;
  pic2=thelabel(btex C etex scaled 7, (25u,5u));
  draw pic2 withcolor white;
  clip currentpicture to unitsquare xscaled 30u yscaled 10u;
endfig;
end.
