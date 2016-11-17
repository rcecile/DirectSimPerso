PRO plot_ps_noosc,saveplot

nsuff=['','_err0.03','_errPpodds']

!p.charsize=2

loadct,39
restore,'temp_ascii.sav'        ; contient dir
dirn="/sps/lsst/data/rcecile/TJP_noBAO_PS/"

namez = ['0.9','1.3','1.8','1.8']
nx= ['_640', '_900', '_1024', '_500']
gsuff = ['', '', ' thin', ' thick']
nz = n_elements(namez)
nk = n_elements(namek)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
!p.multi=0

if (saveplot) then begin
; current plotting device.
   mydevice = !D.NAME
   SET_PLOT, 'PS'
   DEVICE, FILENAME='/sps/lsst/dev/rcecile/Fig/plot_ps_noosc_grids01.eps', /PORTRAIT,/COLOR,XSIZE=8.8,YSIZE=5.5,FONT_SIZE=4
endif else window,0

lpsym=[0,1,0,0]
lcol  = [0,80,150,210]
lpsym=[0,4,0,0]
lth=[2,4,4,4]


!p.multi=0
plot,[0.01,0.1],[1e4,1e4],/xs,/ys,xra=[0.01,0.15],/yl,xtit='wavenumber [Mpc^-1]',ytit='Power spectrum [(Mpc)^-3]',/nodata,$
     yra=[50,4e4],xma=[9,1],yma=[3,1]

for iz=0,1 do begin
   
   
   t  = dirn + 'simu_ps'+nx[iz]+'_z'+namez[iz]+'_ntpk.txt'
   ps = read_ascii(t, template =  TEMP_POW_SPEC_TH)
   xt=ps.(0)
   pt=ps.(2)
   oplot,xt,pt,li=iz*1,th=4


   for isuff=0,n_elements(nsuff)-1 do begin

      suff = nsuff[isuff]
      
      nsim=0
      for j=0,9 do begin
         t  = dirn + 'PS_G2_'+strtrim(j,2)+nx[iz]+'_z'+namez[iz]+suff+'_wngal.txt' 
         print,t
         check = FILE_TEST(t)
         if (check eq 0 ) then continue
         nsim = nsim +  1
         print,'                read now'
         p = read_ascii(t, template =  TEMP_POW_SPEC_TXT) 
         if (j eq 0) then   pref=dblarr(n_elements(p.(0)))
         xs = p.(0)
         oplot,xs,(p.(1)-p.(6)),col=lcol[isuff+1],th=1,li=iz*1
      endfor
 
     print,iz,isuff
   ;   read,xx

   endfor


endfor   

mytext=['spectroZ','Gauss 0.03','photoZ with podds cut']
mytext = ['theoretical spectrum wo BAO',mytext]
legend,mytext,line=0,col=lcol,box=1,/fill,/left,/bottom,charsize=1.5,th=2
legend,'z = '+ namez[0:1]+gsuff[0:1],line=indgen(2),box=1,/fill,/right,/top,charsize=1.5,th=4
oplot,[0.02,0.02],[2e3,1e5],li=2,th=2
oplot,[0.12,0.12],[.8e2,5e3],li=2,th=2

if (saveplot) then begin
   DEVICE, /CLOSE
   SET_PLOT, mydevice
endif

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

if (saveplot) then begin
; current plotting device.
   mydevice = !D.NAME
   SET_PLOT, 'PS'
   DEVICE, FILENAME='/sps/lsst/dev/rcecile/Fig/plot_ps_noosc_grids23.eps', /PORTRAIT,/COLOR,XSIZE=8.8,YSIZE=5.5,FONT_SIZE=4
   endif else window,1

lpsym=[0,1,0,0]
lcol  = [0,80,150,210]
lpsym=[0,4,0,0]
lth=[2,4,4,4]


!p.multi=0
plot,[0.01,0.1],[1e4,1e4],/xs,/ys,xra=[0.01,0.15],/yl,xtit='wavenumber [Mpc^-1]',ytit='Power spectrum [(Mpc)^-3]',/nodata,$
     yra=[10,4e4],xma=[9,1],yma=[3,1]

for iz=2,3 do begin
   
   
   t  = dirn + 'simu_ps'+nx[iz]+'_z'+namez[iz]+'_ntpk.txt'
   ps = read_ascii(t, template =  TEMP_POW_SPEC_TH)
   xt=ps.(0)
   pt=ps.(2)
   oplot,xt,pt,li=iz*1,th=4


   for isuff=0,n_elements(nsuff)-1 do begin

      suff = nsuff[isuff]
      
      nsim=0
      for j=0,9 do begin
         t  = dirn + 'PS_G2_'+strtrim(j,2)+nx[iz]+'_z'+namez[iz]+suff+'_wngal.txt' 
         print,t
         check = FILE_TEST(t)
         if (check eq 0 ) then continue
         nsim = nsim +  1
         print,'                read now'
         p = read_ascii(t, template =  TEMP_POW_SPEC_TXT) 
         if (j eq 0) then   pref=dblarr(n_elements(p.(0)))
         xs = p.(0)
         oplot,xs,(p.(1)-p.(6)),col=lcol[isuff+1],th=1,li=iz*1
      endfor
 
     print,iz,isuff
   ;   read,xx

   endfor


endfor   

mytext=['spectroZ','Gauss 0.03','photoZ with podds cut']
mytext = ['theoretical spectrum wo BAO',mytext]
legend,mytext,line=0,col=lcol,box=1,/fill,/left,/bottom,charsize=1.5,th=2
legend,'z = '+ namez[2:3]+gsuff[2:3],line=indgen(2),box=1,/fill,/right,/top,charsize=1.5,th=4
oplot,[0.02,0.02],[1e3,1e5],li=2,th=2
oplot,[0.06,0.06],[2e2,20e2],li=2,th=2
oplot,[0.12,0.12],[1e1,5e3],li=2,th=2

if (saveplot) then begin
   DEVICE, /CLOSE
   SET_PLOT, mydevice
endif


end
