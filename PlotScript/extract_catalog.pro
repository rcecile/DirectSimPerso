PRO extract_catalog,minicat=1e4

;;les caractéristiques inf et sup 

lim_z_inf = 0.8
lim_z_sup = 1
lim_mag_sup = -20.5
lim_theta_sup = 7./180.*!pi
lim_cube=256


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
   if (minicat gt 0) then m=mrdfits(cat,1,h,range=[0,minicat]) else  m=mrdfits(cat,1,h)

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

h = hist_2d(xc[cube],zc[cube])
contour,h,/xs,/ys

cube = where(abs(xc) lt lim_cube and abs(yc) lt lim_cube and abs(zc) lt lim_cube)

; verification que le cube est bien rempli
h = hist_2d(xc[cube],zc[cube])
window,0,xs=500,ys=500
contour,h,/xs,/ys,xtit='X',ytit='Z'
h = hist_2d(xc[cube],yc[cube])
window,1,xs=500,ys=500
contour,h,/xs,/ys,xtit='X',ytit='Y'
h = hist_2d(zc[cube],yc[cube])
window,2,xs=500,ys=500
contour,h,/xs,/ys,xtit='Z',ytit='Y'

;;fichier
;openw,lun0, '/sps/lsst/data/benaadad/Planck_BAO/cat_z_08_1.txt', /get_lun
;printf, lun0, n_elements(zc)
;printf,  lun0, min(xc),max(xc),min(yc),max(yc),min(zc),max(zc), FORMAT = '(6(F10.2))'
;for i=0,n_elements(zc)-1 do  printf, lun0, xc[i],yc[i],zc[i]
;close, lun0
;free_lun, lun0
stop

END
