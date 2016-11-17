PRO plot_ps_noosc,saveplot
nsuff=['','_err0.03','_errPpodds','_errPBDT']

!p.charsize=2

loadct,39
restore,'temp_ascii.sav'        ; contient dir
dirn="/sps/lsst/data/rcecile/TJP_noBAO_PS/"

namez = ['0.5','0.9','1.3','1.8'];,'1.8']
nx= ['_350','_640', '_900', '_1024', '_500']
gsuff = ['', '', '', '', ' thin', ' thick']
nz = n_elements(namez)
nk = n_elements(namek)
mytext=['spectroZ','Gauss 0.03','photoZ with podds cut']
mytext = ['theoretical spectrum wo BAO',mytext]

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
!p.multi=0

if (saveplot) then begin
; current plotting device.
   mydevice = !D.NAME
   SET_PLOT, 'PS'
   DEVICE, FILENAME='/sps/lsst/dev/rcecile/Fig/plot_ps_noosc_spec.eps', /PORTRAIT,/COLOR,XSIZE=8.8,YSIZE=5.5,FONT_SIZE=4
endif else window,0


!p.multi=0
plot,[0.01,0.1],[1e4,1e4],/xs,/ys,xra=[0.01,0.15],xtit='wavenumber [Mpc^-1]',ytit='Measured PS / theoretical PS',/nodata,$
     yra=[0.6,2],xma=[9,1],yma=[3,1]
 oplot,[0.,0.2],[1,1],th=4,col=0,li=2

mycol=[50,95,210,250]
for iz=0,nz-1 do begin
;for iz=0,1 do begin
   
   t  = dirn + 'simu_ps'+nx[iz]+'_z'+namez[iz]+'_ntpk.txt'
   ps = read_ascii(t, template =  TEMP_POW_SPEC_TH)
   xt=ps.(0)
   pt=ps.(2)

;   for isuff=0,n_elements(nsuff)-1 do begin
   for isuff=0,0 do begin

      suff = nsuff[isuff]

      nsim=0
      for j=0,9 do begin
         t  = dirn + 'PS_G2_'+strtrim(j,2)+nx[iz]+'_z'+namez[iz]+suff+'_wngal.txt' 
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
endfor   

legend,'z = '+ namez+gsuff,col=mycol,line=0,box=1,/fill,/right,/top,charsize=1.5,th=4
oplot,[0.02,0.02],[8e3,1e5],li=2,th=1
oplot,[0.12,0.12],[10e2,5e3],li=2,th=2
xyouts,.045,0.75,'Spectroscopic redshift',charsize=2


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
plot,[0.01,0.1],[1e4,1e4],/xs,/ys,xra=[0.01,0.15],xtit='wavenumber [Mpc^-1]',ytit='Measured PS / theoretical PS',/nodata,$
     yra=[0.,1.1],xma=[9,1],yma=[3,1]
 oplot,[0.,0.2],[1,1],th=4,col=0,li=2

mycol=[50,95,210,250]
for iz=0,nz-1 do begin
;for iz=0,1 do begin
   
   t  = dirn + 'simu_ps'+nx[iz]+'_z'+namez[iz]+'_ntpk.txt'
   ps = read_ascii(t, template =  TEMP_POW_SPEC_TH)
   xt=ps.(0)
   pt=ps.(2)

;   for isuff=0,n_elements(nsuff)-1 do begin
   for isuff=1,1 do begin

      suff = nsuff[isuff]

      nsim=0
      for j=0,9 do begin
         t  = dirn + 'PS_G2_'+strtrim(j,2)+nx[iz]+'_z'+namez[iz]+suff+'_wngal.txt' 
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
endfor   

;legend,mytext[0:1],line=0,col=lcol[0:1],box=1,/fill,/left,/bottom,charsize=1.5,th=2
legend,'z = '+ namez+gsuff,col=mycol,line=0,box=1,/fill,/right,/center,charsize=1.5,th=4
oplot,[0.02,0.02],[8e3,1e5],li=2,th=2
oplot,[0.12,0.12],[10e2,5e3],li=2,th=2
xyouts,.05,0.4,'Gaussian error 0.03',charsize=2


if (saveplot) then begin
   DEVICE, /CLOSE
   SET_PLOT, mydevice
endif

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

if (saveplot) then begin
; current plotting device.
   mydevice = !D.NAME
   SET_PLOT, 'PS'
   DEVICE, FILENAME='/sps/lsst/dev/rcecile/Fig/plot_ps_noosc_podds.eps', /PORTRAIT,/COLOR,XSIZE=8.8,YSIZE=5.5,FONT_SIZE=4
   endif else window,2

lcol  = [0,80,150,210]


plot,[0.01,0.1],[1e4,1e4],/xs,/ys,xra=[0.01,0.15],xtit='wavenumber [Mpc^-1]',ytit='Measured PS / theoretical PS',/nodata,$
     yra=[0.,1.1],xma=[9,1],yma=[3,1]
 oplot,[0.0,0.2],[1,1],th=4,col=0,li=2

mycol=[50,95,210,250]
for iz=0,nz-1 do begin
;for iz=2,3 do begin
   
   t  = dirn + 'simu_ps'+nx[iz]+'_z'+namez[iz]+'_ntpk.txt'
   ps = read_ascii(t, template =  TEMP_POW_SPEC_TH)
   xt=ps.(0)
   pt=ps.(2)

;   for isuff=0,n_elements(nsuff)-1 do begin
   for isuff=2,2 do begin

      suff = nsuff[isuff]

      nsim=0
      for j=0,9 do begin
         t  = dirn + 'PS_G2_'+strtrim(j,2)+nx[iz]+'_z'+namez[iz]+suff+'_wngal.txt' 
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
endfor   

;legend,mytext[0:1],line=0,col=lcol[0:1],box=1,/fill,/left,/bottom,charsize=1.5,th=2
legend,'z = '+ namez+gsuff,col=mycol,line=0,box=1,/fill,/right,/center,charsize=1.5,th=4
oplot,[0.02,0.02],[8e3,1e5],li=2,th=2
oplot,[0.12,0.12],[10e2,5e3],li=2,th=2
xyouts,.035,0.35,'photometric redshift after ODDS cut',charsize=2

if (saveplot) then begin
   DEVICE, /CLOSE
   SET_PLOT, mydevice
endif


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

if (saveplot) then begin
; current plotting device.
   mydevice = !D.NAME
   SET_PLOT, 'PS'
   DEVICE, FILENAME='/sps/lsst/dev/rcecile/Fig/plot_ps_noosc_pbdt.eps', /PORTRAIT,/COLOR,XSIZE=8.8,YSIZE=5.5,FONT_SIZE=4
   endif else window,3

lcol  = [0,80,150,210]


plot,[0.01,0.1],[1e4,1e4],/xs,/ys,xra=[0.01,0.15],xtit='wavenumber [Mpc^-1]',ytit='Measured PS / theoretical PS',/nodata,$
     yra=[0.,1.1],xma=[9,1],yma=[3,1]
 oplot,[0.0,0.2],[1,1],th=4,col=0,li=2

mycol=[50,95,210,250]
for iz=0,nz-1 do begin
;for iz=2,3 do begin
   
   t  = dirn + 'simu_ps'+nx[iz]+'_z'+namez[iz]+'_ntpk.txt'
   ps = read_ascii(t, template =  TEMP_POW_SPEC_TH)
   xt=ps.(0)
   pt=ps.(2)

;   for isuff=0,n_elements(nsuff)-1 do begin
   for isuff=3,3 do begin

      suff = nsuff[isuff]

      nsim=0
      for j=0,9 do begin
         t  = dirn + 'PS_G2_'+strtrim(j,2)+nx[iz]+'_z'+namez[iz]+suff+'_wngal.txt' 
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
endfor   

;legend,mytext[0:1],line=0,col=lcol[0:1],box=1,/fill,/left,/bottom,charsize=1.5,th=2
legend,'z = '+ namez+gsuff,col=mycol,line=0,box=1,/fill,/right,/center,charsize=1.5,th=4
oplot,[0.02,0.02],[8e3,1e5],li=2,th=2
oplot,[0.12,0.12],[10e2,5e3],li=2,th=2
xyouts,.035,0.35,'photometric redshift after BDT cut',charsize=2

if (saveplot) then begin
   DEVICE, /CLOSE
   SET_PLOT, mydevice
endif

end
