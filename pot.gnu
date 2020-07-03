set xlabel "r"
set ylabel "Veff"
set xrange [0:1]
set yrange [-1:10]
V(x,enr,b,g) = -0.5*g/x**4 + enr*b**2./x**2. 
r(x,enr,b,g) = x - sqrt(g/(enr * b**2.))
quad(enr,b,g) = -4. * enr**3. * b**6. / g**2.  
cubi(enr,b,g) = 36. * sqrt(enr**7. * b**14. / g**5.)
quar(enr,b,g) = -300. * enr**4. * b**8. / g**3.
Vts(enr,b,g) = 0.5* enr**2. * b**4. / g
Vhrm(x,enr,b,g) = Vts(enr,b,g) + 0.5*quad(enr,b,g)*r(x,enr,b,g)**2.  
Vvpt2(x,enr,b,g) = Vts(enr,b,g) + 0.5*quad(enr,b,g)*r(x,enr,b,g)**2. + 1./6.*cubi(enr,b,g)*r(x,enr,b,g)**3. + 1./24.*quar(enr,b,g)*r(x,enr,b,g)**4.
