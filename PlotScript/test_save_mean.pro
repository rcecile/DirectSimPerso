PRO test_save_mean,saveplot

namez = ['0.5','0.9','1.5']
namezmean=['0.525', '0.96', '1.54']
nx= ['_160','_300', '_225']
nameerr = ['', '_err0.03', '_errP','_errPBDT9','_errPBDT8']
;nameerr = ['', '_err0.03', '_errPBDT9','_errPBDT8']

loadct,12

lcol  = [95,210,35,135,  120,0]
;lcol  = [95,210,135,  120,0]
lerr=['spectroZ','Gaussian 0.03','photoZ','photoZ BDT 90%','photoZ BDT 80%','theoretical, no oscillation']
;lerr=['spectroZ','Gaussian 0.03','photoZ BDT 90%','photoZ BDT 80%','theoretical, no oscillation']
lpsym=[0,-5,-4]
llin = [0,3,2]

nz = n_elements(namez)
nerr = n_elements(nameerr)
nsim=10

restore,'temp_ascii.sav'   
dir="/sps/lsst/data/rcecile/Planck_BAO_PS/"
dirn="/sps/lsst/data/rcecile/Planck_noBAO_PS/"
myformat = '(F9.7," ",G14.6," ",G14.6," ",G14.6," ",G14.6," ",G14.6," ",G14.6," ",G14.6," ",I6," ",I6," ",G14.6)'
!p.thick=3
!p.charsize=2


if (saveplot ) then begin
; current plotting device.
      mydevice = !D.NAME
      SET_PLOT, 'PS'
      DEVICE, FILENAME='/sps/lsst/dev/rcecile/Fig/ps_undamped.eps', /PORTRAIT,/COLOR,XSIZE=8.8,YSIZE=5.5,FONT_SIZE=5
   endif


;plot,[0.01,0.1],[1,1],/xs,/ys,xra=[0.005,0.15],/nodata,yra=[1100.,.5e5],/yl,xtit='wavenumber [Mpc^-1]',ytit='Power spectrum',xma=[8.5,0.5],yma=[3.5,.5]
plot,[0.01,0.1],[1,1],/xs,/ys,xra=[0.005,0.15],/nodata,yra=[.5,1.5],xtit='wavenumber [Mpc^-1]',xma=[0.5,0.5],yma=[3.5,.5]

err=dblarr(nz,nerr,nsim)
for iz = 0,nz-1 do begin
   
   ; spectre th??orique no osc, no error
   ftheo  = dir + 'simu_ps'+nx[iz]+'_z'+namezmean[iz]+'_ntpk.txt'
   ptheo = read_ascii(ftheo, template =  TEMP_POW_SPEC_TH)     
   xt=ptheo.(0)
   
   for ir = 0, nerr-1 do begin
      suff = nameerr[ir]

      ; spectre observ??
      fobs  = dir + 'PS_cube'+nx[iz]+'_z'+namez[iz]+suff+'_wngal.txt' 
  ;    print,fobs
      pobs = read_ascii(fobs, template =  TEMP_POW_SPEC_TXT)     
      xs = pobs.(0)
      ok = intarr(n_elements(xs))
      for i=0,n_elements(xs)-1 do ok[i] = where( abs(xs[i]-xt) eq min(abs(xs[i]-xt)))

      ; < spectre simul?? no osc, with error >
      pnobao=dblarr(n_elements(xs))
      mynsim = 0
      for j=0,nsim-1 do begin
         fnobao  = dirn + 'PS_cube_'+strtrim(j,2)+nx[iz]+'_z'+namez[iz]+suff+'_wngal.txt' 
         check = FILE_TEST(fnobao)
         if (check eq 0) then continue
         mynsim = mynsim+1
         p = read_ascii(fnobao, template =  TEMP_POW_SPEC_TXT)   
         oplot,xs,(pobs.(1)) / (p.(1)) ,col=lcol[ir],th=1

         pnobao = pnobao + p.(1)
      endfor
      pnobao = pnobao / double(mynsim)
      print,suff, ' N sim OK =',mynsim

      oplot,xs,(pobs.(1)) / (pnobao) ,col=lcol[ir]
   endfor
endfor
;legend,'z = '+ namez+nsuff,col=lcol,box=1,line=0,/fill,/left,/bottom,charsize=1.5
;legend,'z = '+ namez,col=lcol,box=1,line=0,/fill,/left,/bottom,charsize=1.5
legend,lerr,li=[replicate(0,n_elements(nameerr)),2],col=lcol,box=1,/fill,/right,/top,charsize=1.5

oplot,[0,.2],[1,1],th=4
   if (saveplot) then begin
      DEVICE, /CLOSE
      SET_PLOT, mydevice
   endif

END
