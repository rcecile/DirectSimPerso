PRO plot_ps2fit_bias_cata,saveplot,refsuff

!p.charsize=2
!p.thick=3
!p.symsize=4
h= 0.679
h3 = h*h*h

loadct,12
restore,'temp_ascii.sav'        ; contient dir
dir="/sps/lsst/data/rcecile/Planck_BAO_PS/"


nsuff=['_errP','_errPBDT9','_errPBDT8']
lerr=['photoZ','photoZ BDT 90%','photoZ BDT 80%']
lcol  = [35, 135,120]
lpsym=[1,6, 4]

z=[0.5,0.9, 1.5]
namez = ['0.5','0.9','1.5']
namezmean=['0.51', '0.93', '1.57']
nx= ['_120','_225', '_320']
kmax = intarr(n_elements(namez), n_elements(nsuff))
kmax_base = [20, 16, 12]
kmax_base = [12, 12, 12]
kmax[*,0] = kmax_base
kmax[*,1] = kmax_base;-2
kmax[*,2] = kmax_base;-4
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
   DEVICE, FILENAME='/sps/lsst/dev/rcecile/Fig/plot_ps2fit'+psuff+'_z'+strtrim(iz,2)+'.eps', /PORTRAIT,/COLOR,XSIZE=8.8,YSIZE=4.4,FONT_SIZE=4
endif 
if (saveplot) then PAT = bytarr(5,5) else PAT = bytarr(5,5)+255


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
!p.multi=0
;window,iz
myformat='(F10.1)'
ltext ='z = '+ namez
lline= 0

mxtit=[' ',' ','wavenumber [Mpc^-1]']
for iz=0,2  do begin
for ierr=0,2 do begin
   
   if (ierr+iz eq 0) then plot,[0.02,0.1],[1e4,1e4],/xs,/ys,xra=[0.02,lkmax[iz]],ytit='Normalized PS [(Mpc)^-3]',/nodata,yra=[0.5,3.5],xtit=mxtit[iz],xma=[7,2],yma=[3,2]
   
   mysuff = nsuff[ierr]
   t = dir + 'PS'+refsuff+'_2fit'+nx[iz]+'_z'+namez[iz]+mysuff+'_wngal.txt' 
   
   pcorr = read_ascii(t, template =  TEMP_2FIT)     
   xs = pcorr.(0)
                                ; spectre th??orique no osc, no error
   ftheo  = dir + 'simu_ps'+nx[iz]+'_z'+namezmean[iz]+'_ntpk.txt'
   ptheo = read_ascii(ftheo, template =  TEMP_POW_SPEC_TH)     
   xt=ptheo.(0)
   ok = intarr(n_elements(xs))
   for i=0,n_elements(xs)-1 do ok[i] = where( abs(xs[i]-xt) eq min(abs(xs[i]-xt)))
   pref = (ptheo.(2))[ok]
   okk = where(xs ge 0.02 and xs le kmax[iz,ierr]/100.)
   ps = ((pcorr.(1)-pcorr.(6))/pref)[okk]
   oplot,xs[okk],ps+2-iz,col=lcol[ierr],th=3
   print,minmax(xs[okk])
   oplot,xt,(ptheo.(1))/(ptheo.(2))+2-iz,li=2,th=2
   xyouts,0.025,2-iz+1.2,'z='+namez[iz],charsize=2,charthick=2
endfor
endfor
xyouts,0.09,2.65,lerr[0],col=lcol[0],charsize=2,charthick=2
xyouts,0.09,2.45,lerr[1],col=lcol[1],charsize=2,charthick=2
xyouts,0.09,2.25,lerr[2],col=lcol[2],charsize=2,charthick=2

;legend,lerr,line=0,col=lcol,box=1,/fill,/right,/top,charsize=1.5

if (saveplot) then begin
   DEVICE, /CLOSE
   SET_PLOT, mydevice
endif

end
