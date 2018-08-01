FUNCTION sch_integ

Mmin = -24
Mmax = -13

;si = QROMO('fsch', Mmin, Mmax, /double,jmax=1000,K=3)
;print,'QROMO ',si

;si = QSIMP('fsch', Mmin, Mmax, /double,jmax=2000)
;print,'QSIMP ',si

si = QROMB('fsch', Mmin, Mmax, /double,jmax=1000,K=5)
return,  si

END
