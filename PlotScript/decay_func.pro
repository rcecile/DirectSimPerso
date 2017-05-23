FUNCTION decay_func,k,ka

decay=1.4
h=0.679

P = -1. * (k/0.1/h) ^decay
f = 1. + ka[1] * k * exp(P) * sin (k * ka[0])

;stop
return,f

END

