set samples 40
set isosamples 41
set xlabel "E"
set ylabel "b"
set zlabel "Pcl"
set xrange [0:1]
set yrange [0:10]
Pcl(x,y,g) = y < (2*g/x)**0.25 ? 1 : 0
Psc(x,y,g) = (1 + exp(-3.14*(g**2/(x**0.5*y**3) - x**0.5*y*g/2. ) ) )**-1
