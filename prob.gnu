set samples 40
set isosamples 41
set xlabel "E"
set ylabel "b"
set zlabel "Pcl"
set xrange [0:1]
set yrange [0:10]
Pcl(x,y,g) = y < (2*g/x)**0.25 ? 1 : 0
Psc(x,y,g) = (1 + exp(22./7.*(0.5* x**0.5 * y - g/(y**3. * x**0.5) ) ) )**-1
