FUNCTION fsch_par,M
  common sch_var, Mstar, phi_star, alpha

x = 10.d^(-0.4*(M - Mstar))
y=0.4 * alog(10) * phi_star * x^(alpha+1.d) * exp(-1.d*x)
return, y

END
