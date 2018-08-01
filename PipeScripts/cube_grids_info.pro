PRO cube_grids_info,shape

; shape = 0 = cube
; shape = 1 = shell

nslicez = 100
z_max = 3.

cell=8.
Nz=700.
h0=3200.
hmin = h0-Nz/2.*cell
hmax = h0+Nz/2.*cell

zslice = findgen(nslicez+1)*z_max/nslicez

h= dloscom(zslice)
hmin = Dloscom(zslice - (1.+zslice)*0.15)
hmax = Dloscom(zslice + (1.+zslice)*0.15)

;zref=[0.5,0.9,1.3,1.8];,0.5,0.9]
;thick = [140.,125.,75.,65.];,75,75] ;Nz du cube produit par cat_grid
;nxc=[350,640,900,1024,360,675]
;nxc=[160.,320.,480.,60.0]
;cellg=[8.,8.,8.,8.];,8,8]

;BAO
;zref=[0.62,0.95,1.45]
;thick = [125., 125.,175.]
;nxc=[120,225,320]

;BAO2
;zref=[0.61,0.96,1.55]
;thick = [81., 128.,175.]
;nxc=[75,225,320]

zref=[1.3]
thick = [125.]
nxc=[300]

cellg=[8.,8.,8.,8.]
;;;;;;;;;;;;;;;;;;;;;;;;;;;

err=[0., 0.03]
d = Dloscom(zref)
dmin = d - thick/2.*cellg
dmax = d + thick/2.*cellg
if (shape eq 0) then dmax = sqrt(dmax * dmax + 2. * (nxc * cellg * nxc * cellg))
zmin=zref
zmax=zmin
for i=0,n_elements(zref)-1 do zmin[i] = zfrlos(dmin[i],7)
for i=0,n_elements(zref)-1 do zmax[i] = zfrlos(dmax[i],7)

for i=0,n_elements(zref)-1 do begin
   print,'================================================================================================='

   for e=0,n_elements(err)-1 do begin
      z_min =  zmin[i] - (zmin[i]+1.)*err[e] *5. 
      z_max =  zmax[i] + (zmax[i]+1.)*err[e] *5. 
      d_min = Dloscom(z_min)
      d_max = Dloscom(z_max)


      i0 = max(where (h lt d[i]))
      i1 = max(where (h lt d_min))
      i2 = max(where (h lt d_max)) + 1
      i1 = max(where (hmax lt d_min))
      i2 = max(where (hmin lt d_max)) + 1

      derr = d_max - d[i]

      print,'z = ',zref[i], '  d = ', d[i], '  d err= ', derr, '  err = ', err[e], '  central slice = ',i0,'  range [',i1,', ',i2, $
            '] or in Mpc [',d_min,', ',d_max, '] for grids in the redshit range [',z_min,', ',z_max, ']',$
            " volume [Gpc3] = ",thick[i]*nxc[i]*nxc[i]*cellg[i]*cellg[i]*cellg[i]/1.e9
   endfor
   print,' '
endfor
print,'================================================================================================='
print,'Use these ranges to optimize the projection on grids (add cata file in case of errors on redshifts)'
print,'================================================================================================='

END
