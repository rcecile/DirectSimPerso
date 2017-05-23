PRO plot_ps2fit2,saveplot

!p.charsize=2.5
!p.thick=3
!p.symsize=3
h= 0.679
h3 = h*h*h

loadct,12
restore,'temp_ascii.sav'        ; contient dir
dir="/sps/lsst/data/rcecile/Planck_BAO_PS/"


nsuff=['','_errP','_errPBDT8']
lerr=['spectroZ','photoZ','photoZ BDT 80%']
lcol  = [95,35,120]
lpsym=[1,6, 4]
lpsym=[4,4, 4]

z=[0.5,0.9, 1.5]
namez = ['0.5','0.9','1.5']
namezmean=['0.51', '0.93', '1.57']
nx= ['_120','_225', '_320']


kmax = intarr(n_elements(namez), n_elements(nsuff))
kmax[*,0] = [20,18,14]
kmax[*,1] = [18,16,10]
kmax[*,2] = [16,16,10]
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
   DEVICE, FILENAME='/sps/lsst/dev/rcecile/Fig/plot_ps2fit.eps', /PORTRAIT,/COLOR,XSIZE=17.6,YSIZE=12,FONT_SIZE=4
endif
if (saveplot) then PAT = bytarr(5,5) else PAT = bytarr(5,5)+255


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
!p.multi=[0,3,3]
;window,iz
myformat='(F10.1)'
ltext ='z = '+ namez
lline= 0

myxmargin1=[8,0,0]
myxmargin2=[0,0,2]
myymargin1=[0,0,4]
myymargin2=[.5,0,0]
myxtit=[' ',' ','wavenumber k [Mpc^-1]']
myytit=['oscillation PS [(Mpc)^-3]',' ',' ']
for iz=0,2 do begin
for ierr=0,2 do begin
   
   mysuff = nsuff[ierr]
   t  =  dir + 'PS_2fit_SN_Ref'+nx[iz]+'_z'+namez[iz]+mysuff+'_wngal.txt'
   pcorr = read_ascii(t, template =  TEMP_2FIT)     
   xs = pcorr.(0)
   okk = where(xs ge 0.02 and xs le kmax[iz,ierr]/100.)
                            ; spectre th??orique no osc, no error
   ftheo  = dir + 'simu_ps'+nx[iz]+'_z'+namezmean[iz]+'_ntpk.txt'
   ptheo = read_ascii(ftheo, template =  TEMP_POW_SPEC_TH)     
   xt=ptheo.(0)
   ok = intarr(n_elements(xs))
   for i=0,n_elements(xs)-1 do ok[i] = where( abs(xs[i]-xt) eq min(abs(xs[i]-xt)))
   ; print,'SHOT NOISE = ', minmax(pcorr.(6))
   ps = (pcorr.(1)-pcorr.(6))/(ptheo.(2))[ok]           
   psm= (pcorr.(1)-pcorr.(6)-pcorr.(7))/(ptheo.(2))[ok]
   psp= (pcorr.(1)-pcorr.(6)+pcorr.(7))/(ptheo.(2))[ok]
 
   plot,xs[okk],ps[okk],/xs,/ys,xra=[0.021,0.199],yra=[0.91,1.09],ytit=myytit[ierr],xtit=myxtit[iz],xma=[myxmargin1[ierr],myxmargin2[ierr]],yma=[myymargin1[iz],myymargin2[iz]],/nodata
   oplot,xs[okk],ps[okk],col=lcol[ierr],psym=lpsym[ierr]
   errplot,xs[okk],psm[okk],psp[okk],col=lcol[ierr]


  oplot,xt,ptheo.(1)/ptheo.(2),li=2
  oplot,[0,.2],[1,1],li=2,th=1
  ffitw = dir + 'fit_Ref'+skmax[iz,ierr]+nx[iz]+'_z'+namez[iz]+mysuff+'_result.txt'
  print,ffitw
  pfitw = read_ascii(ffitw, template =  TEMP_POW_SPEC_FIT_SA)     

  print,'        FIT with shot-noise sub',(pfitw.(1))[0],' - ',((pfitw.(4)))[0],' /+',((pfitw.(3)))[0],' (amp = ',(pfitw.(5))[0],')'
  oplot,xs,decay_func(xs,[pfitw.(1),pfitw.(5),pfitw.(6)]),col=lcol[ierr],th=2
;  xyouts,0.105,1.05,lerr[ierr],charsize=2
  legend,'z = '+namez[iz],charsize=2,/top,/right,/box

endfor
endfor
legend,lerr,charsize=1.5,li=1,psym=lpsym,col=lcol,/bottom,/right,/box

if (saveplot) then begin
   DEVICE, /CLOSE
   SET_PLOT, mydevice
endif

end
;  TEMP_POW_SPEC_FIT_SA = ASCII_TEMPLATE(ffitw)
; save,TEMP_POW_SPEC_FIT_SA,TEMP_SHOTNOISE,TEMP_SEL_FUNC,TEMP_POW_SPEC_ZXY,TEMP_POW_SPEC_TXT,TEMP_POW_SPEC_TH,TEMP_POW_SPEC_FITCHI2,TEMP_POW_SPEC_FIT,$
; TEMP_POW_SPEC,TEMP_PDF,TEMP_INFOS,TEMP_INFO,TEMP_2FIT,file='temp_ascii.sav' 
