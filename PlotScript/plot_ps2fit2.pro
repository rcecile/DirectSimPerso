PRO plot_ps2fit2,saveplot,iz

!p.charsize=2
!p.thick=2
!p.symsize=4
h= 0.679
h3 = h*h*h

loadct,12
restore,'temp_ascii.sav'        ; contient dir
dir="/sps/lsst/data/rcecile/Planck_BAO_PS/"


nsuff=['','_errP','_errPBDT8']
lerr=['spectroZ','photoZ','photoZ BDT 80%']
lcol  = [95,35,120]
lpsym=[1,6, 4]

z=[0.5,0.9, 1.5]
namez = ['0.5','0.9','1.5']
namezmean=['0.51', '0.93', '1.57']
nx= ['_120','_225', '_320']
kmax = intarr(n_elements(namez), n_elements(nsuff))
kmax[*,0] = [16, 16, 16]
kmax[*,1] = [14,14,14]
kmax[*,2] = [12,12,12]
lkmax = kmax/100.
skmax = strarr(n_elements(namez), n_elements(nsuff))
for i=0,n_elements(namez)-1 do for j=0,n_elements(nsuff)-1 do $
   if (kmax[i,j] ge 10) then skmax[i,j] = '_k0.'+strtrim(kmax[i,j],2) $
   else skmax[i,j] = '_k0.0'+strtrim(kmax[i,j],2) 
print,skmax
nz = n_elements(namez)


if (saveplot) then begin
; current plotting device.
   mydevice = !D.NAME
   SET_PLOT, 'PS'
   DEVICE, FILENAME='/sps/lsst/dev/rcecile/Fig/plot_ps2fit'+'_z'+strtrim(iz,2)+'.eps', /PORTRAIT,/COLOR,XSIZE=8.8,YSIZE=4.4,FONT_SIZE=4
endif else window,iz

if (saveplot) then PAT = bytarr(5,5) else PAT = bytarr(5,5)+255


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
!p.multi=0
;window,iz
myformat='(F10.1)'
ltext ='z = '+ namez
lline= 0

mxtit=[' ',' ','wavenumber [Mpc^-1]']
for ierr=0,2 do begin

   if (ierr eq 0) then plot,[0.02,0.1],[1e4,1e4],/xs,/ys,xra=[0.02,lkmax[iz]],ytit='Normalized PS [(Mpc)^-3]',/nodata,yra=[0.5,3.5],xtit=mxtit[iz],xma=[7,2],yma=[3,2]
   
   mysuff = nsuff[ierr]
 ;  t  =  dir + 'PS_2fit5grids'+nx[iz]+'_z'+namez[iz]+mysuff+'_wngal.txt'
   t  =  dir + 'PS_2fit_newErr'+nx[iz]+'_z'+namez[iz]+mysuff+'_wngal.txt'
   pcorr = read_ascii(t, template =  TEMP_2FIT)     
   xs = pcorr.(0)
                                ; spectre th??orique no osc, no error
   ftheo  = dir + 'simu_ps'+nx[iz]+'_z'+namezmean[iz]+'_ntpk.txt'
   ptheo = read_ascii(ftheo, template =  TEMP_POW_SPEC_TH)     
   xt=ptheo.(0)
   ok = intarr(n_elements(xs))
   for i=0,n_elements(xs)-1 do ok[i] = where( abs(xs[i]-xt) eq min(abs(xs[i]-xt)))
   
   ps = (pcorr.(1)-pcorr.(6))/(ptheo.(2))[ok]
   psm= (pcorr.(1)-pcorr.(6)-pcorr.(7))/(ptheo.(2))[ok] +2-ierr
   psp= (pcorr.(1)-pcorr.(6)+pcorr.(7))/(ptheo.(2))[ok] +2-ierr
   
   
   
   ;fit
 ;  ffit = dir + 'fit_SN5grids'+skmax[iz,ierr]+nx[iz]+'_z'+namez[iz]+mysuff+'_result.txt'
   ffit = dir + 'fit_newErr'+skmax[iz,ierr]+nx[iz]+'_z'+namez[iz]+mysuff+'_result.txt'
  print,ffit
; if (ir eq 3 and iz eq 2) then ffit = dir + 'fit_cube_Err'+nx[iz]+'_z'+namez[iz]+'_k0.03'+mysuff+'_result.txt'
  pfit = read_ascii(ffit, template =  TEMP_POW_SPEC_FIT)     

  oplot,[pfit.(1)-pfit.(4),pfit.(1)+pfit.(3),pfit.(1)+pfit.(3),pfit.(1)-pfit.(4),pfit.(1)-pfit.(4)],[.55,.55,1.45,1.45,.55]+2-ierr,col=lcol[ierr],th=5
  print,'FIT ',pfit.(1)-pfit.(4),pfit.(1)+pfit.(3)

  okk = where(xs ge 0.02 and xs le kmax[iz,ierr]/100.)

  oplot,xs[okk]+0.0003*ierr,ps[okk]+2-ierr,col=lcol[ierr],psym=lpsym[ierr],th=2
  errplot,xs[okk]+0.0003*ierr,psm[okk],psp[okk],col=lcol[ierr]
  oplot,xs,decay_func(xs,pfit.(1),h)+2-ierr,li=2,th=3 ,col=lcol[ierr]
   oplot,xt,ptheo.(1)/ptheo.(2)+2-ierr,li=2
  oplot,[1.,1.]*2.*!pi/151.,[0,4],th=2
  xyouts,0.0525,3-ierr-0.4,lerr[ierr],charsize=2
  xyouts,lkmax[iz,0]*.8,3.2,'z = '+namez[iz],charsize=3

endfor

if (saveplot) then begin
   DEVICE, /CLOSE
   SET_PLOT, mydevice
endif

end
;err_tot=[0.030515025 ,0.028397875, 0.037388461]
;err90_tot=[0.025927096 ,0.027347027, 0.032890236]
;err80_tot=[0.023046306 ,0.026009080, 0.030297966]

;print,err_tot/1.34900
;print,err90_tot/1.34900
;print,err80_tot/1.34900
