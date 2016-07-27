PRO cube_grids_info


nslicez = 100
z_max = 3.
z_step = z_max / nslicez
zslice = findgen(nslicez+1)*z_step

zref=[0.7,1.4]
thick = [200.,200.] ;Nz du cube produit par cat_grid
nxc=[450,900]
cell=8.

;;;;;;;;;;;;;;;;;;;;;;;;;;;

err=[0., 0.03]
d = Dloscom(zref)
dmin = d - thick/2.*cell
dmax = d + thick/2.*cell
zmin=zref
zmax=zmin
for i=0,n_elements(zref)-1 do zmin[i] = zfrlos(dmin[i],7)
for i=0,n_elements(zref)-1 do zmax[i] = zfrlos(dmax[i],7)


for i=0,n_elements(zref)-1 do begin
   print,'================================================================================================='

   for e=0,n_elements(err)-1 do begin
      z_min =  zmin[i] - (zmin[i]+1.)*err[e] *3. 
      z_max =  zmax[i] + (zmax[i]+1.)*err[e] *3. 
      i0 = fix(zref[i]/z_step + 0.5)
      i1 = fix(z_min   /z_step + 0.5) - 1
      i2 = fix(z_max   /z_step + 0.5) + 1

      dmin = Dloscom(z_min)
      dmax = Dloscom(z_max)
      derr = dmax - d[i]

      print,'z = ',zref[i], '  d = ', d[i], '  d err= ', derr, '  err = ', err[e], '  central slice = ',i0,'  range [',i1,', ',i2, $
            '] or in Mpc [',dmin,', ',dmax, '] for grids in the redshit range [',z_min,', ',z_max, ']'
      print," volume [Gpc3] = ",thick[i]*nxc[i]*nxc[i]*cell*cell*cell/1.e9
   endfor
   print,' '
endfor
print,'================================================================================================='
print,'Use these ranges to optimize the projection on grids (add cata file in case of errors on redshifts)'
print,'================================================================================================='

END
