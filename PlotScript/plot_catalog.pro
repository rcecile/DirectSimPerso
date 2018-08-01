PRO plot_catalog,minicat

;;les caractéristiques inf et sup 

lim_z_inf = 0.8
lim_z_sup = 1
lim_mag_sup = -20.5
lim_theta_sup = 6./180.*!pi


;;boucle sur les catalogues
test=0

for i=3,99 do begin
   cat="/sps/lsst/data/rcecile/Planck_BAO2/cat_lfZuccaAllFalse_zOrd_Slice"+strtrim(i,2)+".fits"
   hh = headfits(cat,ext=1)
   zmin = 1.d * sxpar(hh,'ZCAT_MIN')
   zmax = 1.d * sxpar(hh,'ZCAT_MAX')
   print,i, zmin, zmax

   if(zmin gt lim_z_sup OR zmax lt lim_z_inf) then continue  ;; on passe a la tranche suivante
   print,minicat
   if (minicat gt 0) then m=mrdfits(cat,1,h,range=[0,1e6]) else  m=mrdfits(cat,1,h)

;;coupure en magnitude et theta

   z0=(m.(3))
   Mag0=(m.(5))
   print,'avant coupure Ngal= ',n_elements(z0)
   ok=where( m.(2) le lim_theta_sup and Mag0 le lim_mag_sup  ,nok)
   print,'apres coupure Ngal= ',nok

;;coordonnés

   d0 = dloscom(z0[ok])
   xc0=d0*sin((m.(2))[ok])*cos((m.(1))[ok])
   yc0=d0*sin((m.(2))[ok])*sin((m.(1))[ok])
   zc0=d0*cos((m.(2))[ok])

   if (test eq 0 ) then begin
      xc = xc0
      yc = yc0
      zc = zc0
   endif else begin
      xc = [xc,xc0]
      yc = [yc,yc0]
      zc = [zc,zc0]
   end
   print,'             ', n_elements(zc)
   test=1
endfor

zc = zc - mean(zc)

h = hist_2d(xc,zc)
contour,h

stop 
ok2 = where(abs(xc) lt x_max and abs(yc) lt y_max and abs(zc) lt z_max)
print, minmax(xc),minmax(yc),minmax(zc)
print, min(xc[ok2]),max(xc[ok2]),min(yc[ok2]),max(yc[ok2]),min(zc[ok2]),max(zc[ok2])
xc = xc[ok2]
yc = yc[ok2]
zc = zc[ok2]

plot,xc,yc,/xs,/ys,psym=3,xtit='x Mpc/h',ytit='y Mpc/h',title='cata_Mag < -20.5   '



;;fichier
;openw,lun0, '/sps/lsst/data/benaadad/Planck_BAO/cat_z_08_1.txt', /get_lun
;printf, lun0, n_elements(zc)
;printf,  lun0, min(xc),max(xc),min(yc),max(yc),min(zc),max(zc), FORMAT = '(6(F10.2))'
;for i=0,n_elements(zc)-1 do  printf, lun0, xc[i],yc[i],zc[i]
;close, lun0
;free_lun, lun0
stop
END
