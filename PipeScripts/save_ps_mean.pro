PRO save_ps_mean,saveplot

namez = ['0.5','0.9','1.5']
;namezmean=['0.525', '0.96', '1.54']
;nx= ['_160','_300', '_225']
nx= ['_120','_225', '_320']
namezmean=['0.51', '0.93', '1.57']
nameerr = ['', '_err0.03', '_errP','_errPBDT9','_errPBDT8']

loadct,12

lcol  = [95,210,35,135,  120,0]
;lcol  = [95,210,135,  120,0]
lerr=['spectroZ','Gaussian 0.03','photoZ','photoZ BDT 90%','photoZ BDT 80%','fiducial, no oscillation']
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


plot,[0.01,0.1],[1,1],/xs,/ys,xra=[0.005,0.15],/nodata,yra=[1100.,.5e5],/yl,xtit='wavenumber [Mpc^-1]',ytit='Undamped Power spectrum',xma=[4,0.5],yma=[3.,.5],ytickn=replicate(" " ,10)

err=dblarr(nz,nerr,nsim)
for iz = 0,nz-1 do begin
   
   ; spectre th??orique no osc, no error
   ftheo  = dirn + 'simu_ps'+nx[iz]+'_z'+namezmean[iz]+'_ntpk.txt'
   ptheo = read_ascii(ftheo, template =  TEMP_POW_SPEC_TH)     
   xt=ptheo.(0)
   
   fobs  = dir + 'PS_newErr'+nx[iz]+'_z'+namez[iz]+'_G0_wngal.txt' 
   pobs = read_ascii(fobs, template =  TEMP_POW_SPEC_TXT)   
   xs = pobs.(0)

   for ir = 0, nerr-1 do begin
      pnoise=dblarr(n_elements(xs))
      npnoise = 0
      for igd =0,4 do begin
         suff = nameerr[ir]
                                ; spectre observ??
         fobs  = dir + 'PS_newErr'+nx[iz]+'_z'+namez[iz]+suff+'_G'+strtrim(igd,2)+'_wngal.txt' 
                                ;    print,fobs
         pobs = read_ascii(fobs, template =  TEMP_POW_SPEC_TXT)     
         if (igd eq 0) then p2fit =pobs.(1) * 0.
         ok = intarr(n_elements(xs))
         for i=0,n_elements(xs)-1 do ok[i] = where( abs(xs[i]-xt) eq min(abs(xs[i]-xt)))
         
                                ; < spectre simul?? no osc, with error >
         pnobao=dblarr(n_elements(xs))
         pshot=dblarr(n_elements(xs),nsim)
         mynsim = 0
         for j=0,nsim-1 do begin
            fnobao  = dirn + 'PS_newErr_'+strtrim(j,2)+nx[iz]+'_z'+namez[iz]+suff+'_G'+strtrim(igd,2)+'_wngal.txt' 
            ;print,fnobao
            check = FILE_TEST(fnobao)
            if (check eq 0) then continue
            mynsim = mynsim+1
            p = read_ascii(fnobao, template =  TEMP_POW_SPEC_TXT)     
            pnobao = pnobao + p.(1) - p.(6)
            pnoise = pnoise +  p.(6)
            npnoise = npnoise+1
            pshot[*,j] = p.(6)
         endfor
         pnobao = pnobao / double(mynsim)
         print,suff, '   Grille n.',igd,' N sim OK =',mynsim
         std_shot = dblarr(n_elements(xs))
         for i=0,n_elements(xs)-1 do std_shot[i] = sqrt(mean(pshot[i,where(pshot[i,*] ne 0)]))
         okerr = where(p.(0) gt 0.02 and p.(0) lt 0.1)
         for j=0,nsim-1 do begin
            fnobao  = dirn + 'PS_newErr_'+strtrim(j,2)+nx[iz]+'_z'+namez[iz]+suff+'_G'+strtrim(igd,2)+'_wngal.txt' 
            check = FILE_TEST(fnobao)
            if (check eq 0) then continue
            p = read_ascii(fnobao, template =  TEMP_POW_SPEC_TXT)    
            rap = p.(1)/pnobao
            err[iz,ir,j] = stddev(rap[okerr])
                                ;   print,j
         endfor
                                ;    print,iz,ir,' ',mean(err[iz,ir,*])
         
          p2fit = p2fit + (pobs.(1)-pobs.(6)) / (pnobao) *(ptheo.(2)[ok])
         
        ; oplot,xs,(pobs.(1)-pobs.(6)) / (pnobao) *(ptheo.(2)[ok]),col=lcol[ir]
 
       endfor
      pnoise = pnoise / npnoise
      p2fit = p2fit / 5.
      oplot,xs,p2fit,col=lcol[ir]

      fname = dir + 'PS_2fit_newErr'+nx[iz]+'_z'+namez[iz]+suff+'_wngal.txt' 
  ;    print,fname
      openw,lun1, fname, /get_lun
 
      newp1 = p2fit + pnoise
      p7shot = pobs.(7)
      print,'SHOT  ',iz, ir, (pobs.(7))[25],p7shot[25]
      newp7 = (newp1-pobs.(6))/(pobs.(1)-pobs.(6))*p7shot ; erreur hors shot-noise"undampee", comme le signal
;      newp7 = (newp1-pnoise)/(pobs.(1)-pnoise)*p7shot ; erreur hors shot-noise"undampee", comme le signal
;      newp7 = newp7*sqrt(1.+1./10.)                     ; car 10
;      simus sans BAOs pour la base
     for i=0,n_elements(p2fit)-1 do  printf, lun1, (pobs.(0))[i], $
                                              newp1[i], $ ;spectre
                                              (pobs.(2))[i], (pobs.(3))[i], (pobs.(4))[i], (pobs.(5))[i], (pobs.(6))[i], $
                                              newp7[i], $ ; erreur
                                              (pobs.(8))[i], (pobs.(9))[i], (pobs.(10))[i],format=myformat
      close, lun1 
      free_lun, lun1
  ;    print,'  ',fname
   endfor
   oplot,xt,ptheo.(2),li=2
;legend,'z = '+ namez+nsuff,col=lcol,box=1,line=0,/fill,/left,/bottom,charsize=1.5
;legend,'z = '+ namez,col=lcol,box=1,line=0,/fill,/left,/bottom,charsize=1.5
endfor
legend,lerr,li=[replicate(0,n_elements(nameerr)),2],col=lcol,box=1,/fill,/right,/top,charsize=1.25

   if (saveplot) then begin
      DEVICE, /CLOSE
      SET_PLOT, mydevice
   endif

END
