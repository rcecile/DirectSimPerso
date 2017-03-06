PRO plot_ps_noosc,saveplot
nsuff=['','_err0.03','_errP','_errPBDT8']
mycol=[80,205,240]

!p.charsize=2

loadct,39
restore,'temp_ascii.sav'        ; contient dir
dirn="/sps/lsst/data/rcecile/Planck_noBAO_PS/"

namez = ['0.5','0.9','1.5']
namezmean = ['0.51','0.93','1.57']
nx=['_120','_225','_320']
nz = n_elements(namez)
nk = n_elements(namek)
mytext=['spectroZ','Gauss 0.03','photoZ with podds cut']
mytext = ['theoretical spectrum wo BAO',mytext]

kBAO =2.*!pi/150.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
!p.multi=0

if (saveplot) then begin
; current plotting device.
   mydevice = !D.NAME
   SET_PLOT, 'PS'
   DEVICE, FILENAME='/sps/lsst/dev/rcecile/Fig/plot_ps_noosc_spec.eps', /PORTRAIT,/COLOR,XSIZE=8.8,YSIZE=5.5,FONT_SIZE=4
endif else window,0


suff = nsuff[1]
nref=0
for j=0,9 do begin
   for ig=0,4 do begin
      t  = dirn + 'PS_'+strtrim(j,2)+nx[iz]+'_z'+namez[iz]+suff+'_G0_wngal.txt' 
                                ; print,t
      check = FILE_TEST(t)
      if (check eq 0 ) then continue
                                ; print,'                read now'
      p = read_ascii(t, template =  TEMP_POW_SPEC_TXT) 
      if (j +ig eq 0) then   PSref=dblarr(n_elements(p.(0)))
      if (j +ig eq 0) then   SNref=dblarr(n_elements(p.(0)))
      PSref = PSref+p.(1)
      SNref = SNref+p.(6)
      nref = nref+1
   endfor
endfor
print,'REF computed from ',nref,' spectra'

plot,[0.02,0.1],[1e4,1e4],/xs,/ys,xra=[0.02,0.20],xtit='wavenumber [Mpc^-1]',ytit='Measured PS / theoretical PS',/nodata,$
     yra=[0.6,2],xma=[9,1],yma=[3,1],/xl
oplot,[kBAO,kBAO],[0,2],th=4,col=0,li=2

isuff = 0
for iz=0,nz-1 do begin
   
   xyouts,0.041,1.4,'150 Mpc',ori=90 

   t  = dirn + 'simu_ps'+nx[iz]+'_z'+namezmean[iz]+'_ntpk.txt'
   ps = read_ascii(t, template =  TEMP_POW_SPEC_TH)
   xt=ps.(0)
   pt=ps.(2)

   suff = nsuff[isuff]
   
   nsim=0
   for j=0,9 do begin
      t  = dirn + 'PS_'+strtrim(j,2)+nx[iz]+'_z'+namez[iz]+suff+'_G0_wngal.txt' 
                                ; print,t
      check = FILE_TEST(t)
      if (check eq 0 ) then continue
      nsim = nsim +  1
                                ; print,'                read now'
      p = read_ascii(t, template =  TEMP_POW_SPEC_TXT) 
      if (j eq 0) then   pref=dblarr(n_elements(p.(0)))
      xs = p.(0)
      ok = intarr(n_elements(xs))
      for i=0,n_elements(xs)-1 do ok[i] = where( abs(xs[i]-xt) eq min(abs(xs[i]-xt)))
      oplot,xs,(p.(1))/pt[ok],col=mycol[iz],th=1,li=1
      oplot,xs,(p.(1)-p.(6))/pt[ok],col=mycol[iz],th=1
   endfor
endfor
print,nsim

legend,'z = '+ namez,col=mycol,line=0,box=1,/fill,/right,/center,charsize=1.5,th=4
oplot,[0.02,0.02],[8e3,1e5],li=2,th=1
oplot,[0.12,0.12],[10e2,5e3],li=2,th=2
xyouts,.043,1.85,'Spectroscopic redshift',charsize=2


if (saveplot) then begin
   DEVICE, /CLOSE
   SET_PLOT, mydevice
endif

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

if (saveplot) then begin
; current plotting device.
   mydevice = !D.NAME
   SET_PLOT, 'PS'
   DEVICE, FILENAME='/sps/lsst/dev/rcecile/Fig/plot_ps_noosc_gauss.eps', /PORTRAIT,/COLOR,XSIZE=8.8,YSIZE=5.5,FONT_SIZE=4
endif else window,1


!p.multi=0
plot,[0.02,0.1],[1e4,1e4],/xs,/ys,xra=[0.02,0.20],xtit='wavenumber [Mpc^-1]',ytit='Measured PS / theoretical PS',/nodata,$
     yra=[0.,1.1],xma=[9,1],yma=[3,1],/xl
 oplot,[0.,0.2],[1,1],th=4,col=0,li=2
oplot,[kBAO,kBAO],[0,2],th=4,col=0,li=2

isuff = 1
for iz=0,nz-1 do begin   
   xyouts,0.041,.61,'150 Mpc',ori=90  

   t  = dirn + 'simu_ps'+nx[iz]+'_z'+namezmean[iz]+'_ntpk.txt'
   ps = read_ascii(t, template =  TEMP_POW_SPEC_TH)
   xt=ps.(0)
   pt=ps.(2)

   suff = nsuff[isuff]
   
   nsim=0
   for j=0,9 do begin
      t  = dirn + 'PS_'+strtrim(j,2)+nx[iz]+'_z'+namez[iz]+suff+'_G0_wngal.txt' 
                                ; print,t
      check = FILE_TEST(t)
      if (check eq 0 ) then continue
      nsim = nsim +  1
                                ; print,'                read now'
      p = read_ascii(t, template =  TEMP_POW_SPEC_TXT) 
      if (j eq 0) then   pref=dblarr(n_elements(p.(0)))
      xs = p.(0)
      ok = intarr(n_elements(xs))
      for i=0,n_elements(xs)-1 do ok[i] = where( abs(xs[i]-xt) eq min(abs(xs[i]-xt)))
         oplot,xs,(p.(1))/pt[ok],col=mycol[iz],th=1,li=1
         oplot,xs,(p.(1)-p.(6))/pt[ok],col=mycol[iz],th=1
      endfor
endfor
print,nsim

;legend,mytext[0:1],line=0,col=lcol[0:1],box=1,/fill,/left,/bottom,charsize=1.5,th=2
legend,'z = '+ namez,col=mycol,line=0,box=1,/fill,/right,/center,charsize=1.5,th=4
oplot,[0.02,0.02],[8e3,1e5],li=2,th=2
oplot,[0.12,0.12],[10e2,5e3],li=2,th=2
xyouts,.043,0.98,'Gaussian error 0.03(1+z)',charsize=2


if (saveplot) then begin
   DEVICE, /CLOSE
   SET_PLOT, mydevice
endif

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

if (saveplot) then begin
; current plotting device.
   mydevice = !D.NAME
   SET_PLOT, 'PS'
   DEVICE, FILENAME='/sps/lsst/dev/rcecile/Fig/plot_ps_noosc_errP.eps', /PORTRAIT,/COLOR,XSIZE=8.8,YSIZE=5.5,FONT_SIZE=4
   endif else window,2



plot,[0.02,0.1],[1e4,1e4],/xs,/ys,xra=[0.02,0.20],xtit='wavenumber [Mpc^-1]',ytit='Measured PS / theoretical PS',/nodata,$
     yra=[0.,1.1],xma=[9,1],yma=[3,1],/xl
oplot,[kBAO,kBAO],[0,2],th=4,col=0,li=2

isuff = 2
for iz=0,nz-1 do begin
   xyouts,0.041,.61,'150 Mpc',ori=90  

   
   t  = dirn + 'simu_ps'+nx[iz]+'_z'+namezmean[iz]+'_ntpk.txt'
   ps = read_ascii(t, template =  TEMP_POW_SPEC_TH)
   xt=ps.(0)
   pt=ps.(2)
   
   suff = nsuff[isuff]
   
   nsim=0
   for j=0,9 do begin
      t  = dirn + 'PS_'+strtrim(j,2)+nx[iz]+'_z'+namez[iz]+suff+'_G0_wngal.txt' 
                                ; print,t
      check = FILE_TEST(t)
      if (check eq 0 ) then continue
      nsim = nsim +  1
                                ; print,'                read now'
      p = read_ascii(t, template =  TEMP_POW_SPEC_TXT) 
      if (j eq 0) then   pref=dblarr(n_elements(p.(0)))
      xs = p.(0)
      ok = intarr(n_elements(xs))
      for i=0,n_elements(xs)-1 do ok[i] = where( abs(xs[i]-xt) eq min(abs(xs[i]-xt)))
      oplot,xs,(p.(1))/pt[ok],col=mycol[iz],th=1,li=1
      oplot,xs,(p.(1)-p.(6))/pt[ok],col=mycol[iz],th=1
   endfor
endfor
print,nsim

legend,'z = '+ namez,col=mycol,line=0,box=1,/fill,/right,/center,charsize=1.5,th=4
oplot,[0.02,0.02],[8e3,1e5],li=2,th=2
oplot,[0.12,0.12],[10e2,5e3],li=2,th=2
xyouts,.043,0.98,'photometric redshift',charsize=2

if (saveplot) then begin
   DEVICE, /CLOSE
   SET_PLOT, mydevice
endif


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

if (saveplot) then begin
; current plotting device.
   mydevice = !D.NAME
   SET_PLOT, 'PS'
   DEVICE, FILENAME='/sps/lsst/dev/rcecile/Fig/plot_ps_noosc_pbdt8.eps', /PORTRAIT,/COLOR,XSIZE=8.8,YSIZE=5.5,FONT_SIZE=4
   endif else window,3


plot,[0.02,0.1],[1e4,1e4],/xs,/ys,xra=[0.02,0.20],xtit='wavenumber [Mpc^-1]',ytit='Measured PS / theoretical PS',/nodata,$
     yra=[0.,1.1],xma=[9,1],yma=[3,1],/xl
oplot,[kBAO,kBAO],[0,2],th=4,col=0,li=2

isuff = 3
for iz=0,nz-1 do begin
   xyouts,0.041 ,.61,'150 Mpc',ori=90  

   
   t  = dirn + 'simu_ps'+nx[iz]+'_z'+namezmean[iz]+'_ntpk.txt'
   ps = read_ascii(t, template =  TEMP_POW_SPEC_TH)
   xt=ps.(0)
   pt=ps.(2)

   suff = nsuff[isuff]
   
   nsim=0
   for j=0,9 do begin
      t  = dirn + 'PS_'+strtrim(j,2)+nx[iz]+'_z'+namez[iz]+suff+'_G0_wngal.txt' 
                                ; print,t
      check = FILE_TEST(t)
      if (check eq 0 ) then continue
      nsim = nsim +  1
                                ; print,'                read now'
      p = read_ascii(t, template =  TEMP_POW_SPEC_TXT) 
      if (j eq 0) then   pref=dblarr(n_elements(p.(0)))
      xs = p.(0)
      ok = intarr(n_elements(xs))
      for i=0,n_elements(xs)-1 do ok[i] = where( abs(xs[i]-xt) eq min(abs(xs[i]-xt)))
      oplot,xs,(p.(1))/pt[ok],col=mycol[iz],th=1,li=1
      oplot,xs,(p.(1)-p.(6))/pt[ok],col=mycol[iz],th=1
   endfor
endfor
print,nsim

legend,'z = '+ namez,col=mycol,line=0,box=1,/fill,/right,/center,charsize=1.5,th=4
xyouts,.043,0.98,'photometric redshift',charsize=2
xyouts,.043,0.9,'after BDT 80% cut',charsize=2

if (saveplot) then begin
   DEVICE, /CLOSE
   SET_PLOT, mydevice
endif


end
