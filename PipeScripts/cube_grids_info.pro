PRO cube_grids_info

nslicez = 70
cell=8.
Nz=700.
h0=3200.
hmin = h0-Nz/2.*cell
hmax = h0+Nz/2.*cell
h= h0 + (findgen(nslicez+1)-nslicez/2)*Nz*cell/nslicez

zref=[0.9,1.3,1.8,1.8]
thick = [125.,75,65,75] ;Nz du cube produit par cat_grid
nxc=[640,900,1024,500]
cellg=[8,8,8,16]
;;;;;;;;;;;;;;;;;;;;;;;;;;;

err=[0., 0.03]
d = Dloscom(zref)
dmin = d - thick/2.*cellg
dmax = d + thick/2.*cellg
zmin=zref
zmax=zmin
for i=0,n_elements(zref)-1 do zmin[i] = zfrlos(dmin[i],7)
for i=0,n_elements(zref)-1 do zmax[i] = zfrlos(dmax[i],7)


for i=0,n_elements(zref)-1 do begin
   print,'================================================================================================='

   for e=0,n_elements(err)-1 do begin
      z_min =  zmin[i] - (zmin[i]+1.)*err[e] *3. 
      z_max =  zmax[i] + (zmax[i]+1.)*err[e] *3. 
      d_min = Dloscom(z_min)
      d_max = Dloscom(z_max)


      i0 = max(where (h lt d[i]))
      i1 = max(where (h lt d_min))
      i2 = max(where (h lt d_max)) + 1

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
