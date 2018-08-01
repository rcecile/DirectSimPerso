FUNCTION fsch,M

range = 2

if (range eq 0) then begin ;; z < 0.5
   Mstar = -21.22d
   alpha = -1.37d
   phi_star = 0.00281
   
   EMstar = -20.99d
   Ealpha = -0.64d
   Ephi_star = 0.00152
   
   LMstar = -21.d
   Lalpha = -1.35d
   Lphi_star = 0.00212
   
   SMstar = -18.72d
   Salpha = -1.02d
   Sphi_star = 0.00426
endif


if (range eq 1) then begin ;; 0.5 <  z < 0.75
   Mstar = -21.79d
   alpha = -1.37d
   phi_star = 0.00202
   
   EMstar = -21.01d
   Ealpha = -0.39d
   Ephi_star = 0.00167
   
   LMstar = -21.29d
   Lalpha = -1.14d
   Lphi_star = 0.00234
   
   SMstar = -19.51d
   Salpha = -1.04d
   Sphi_star = 0.00278
endif

if (range eq 2) then begin ;; 0.75 <  z 
   Mstar = -21.62d
   alpha = -1.37d
   phi_star = 0.00195
   
   EMstar = -21.44d
   Ealpha = -0.72d
   Ephi_star = 0.00062
   
   LMstar = -21.15d
   Lalpha = -0.87d
   Lphi_star = 0.00222
   
   SMstar = -20.32d
   Salpha = -1.3d
   Sphi_star = 0.00177
endif


x = 10.d^(-0.4*(M - EMstar))
EphiM = 0.4 * alog(10) * Ephi_star * x^(Ealpha+1.d) * exp(-1.d*x)


x = 10.d^(-0.4*(M - LMstar))
LphiM = 0.4 * alog(10) * Lphi_star * x^(Lalpha+1.d) * exp(-1.d*x)


x = 10.d^(-0.4*(M - SMstar))
SphiM = 0.4 * alog(10) * Sphi_star * x^(Salpha+1.d) * exp(-1.d*x)


x = 10.d^(-0.4*(M - Mstar))
phiM = 0.4 * alog(10) * phi_star * x^(alpha+1.d) * exp(-1.d*x)

;return, phiM
return, EphiM+LphiM+SphiM

END
;0.11426614
;0.0034600820
; 0.071487924
; 0.021148100
;f=[0.0034600820,0.071487924,0.021148100]
;print,f/total(f)
