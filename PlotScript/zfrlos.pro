FUNCTION zfrlos,dlos,niter

;double CosmoCalc::ZFrLos(double loscom /* Mpc com */, int niter)
;// Recherche du redshift correspondant a une distance comobile
;// le long de la ligne de visee (radiale) egale a "loscom" Mpc
;// niter = nomber of iterations for precision measurement

if (niter lt 3) then niter = 5
dz = 0.1d
zmin=0.d 
zmax=0.d
while (dloscom(zmax) lt dlos) do zmax = zmax + dz

if (zmax eq 0.) then return,0.

for i=0,niter do begin
   zmin=zmax-dz
   if(zmin lt 0.) then zmin=0.  
   dz = dz / 10.d

   for z=zmin,zmax,dz do begin
      d = dloscom(z) 
      if (d lt dlos) then continue
      zmax = z
      break
   endfor
endfor

return, zmax

end
