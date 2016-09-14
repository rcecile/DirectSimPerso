FUNCTION decay_func,k,ka,h

amp=2.5
decay=1.4
P = -1. * exp(decay * alog(k)) / exp(decay * alog(0.1 * h))

f = 1. + amp * k * exp(P) * sin (2.*!pi*k / ka[0])
;stop
return,f

END

