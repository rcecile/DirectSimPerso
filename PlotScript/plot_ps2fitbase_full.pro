PRO plot_ps2fitbase,saveplot

!p.charsize=2.5
!p.thick=3
!p.symsize=3
h= 0.679
h3 = h*h*h

loadct,12
restore,'temp_ascii_new.sav'        ; contient dir
dir="/sps/lsst/data/rcecile/Planck_BAO2_PS/"


nsuff=['','_errP','_errPBDT8']
lerr=['spectroZ','photoZ','photoZ BDT 80%']
lcol  = [95,35,120]
lpsym=[1,6, 4]
lpsym=[4,4, 4]

z=[0.5,0.9, 1.5]
namez = ['0.5','0.9','1.3']
namezmean=['0.51', '0.93', '1.36']
nx= ['_120','_225', '_300']

kfid = 0.0417
kmin=[4,3,2]
kmax=[20,20,20,20,20,    20,15,15,15,15,    20,12,12,12,12]

nz = n_elements(namez)

if (saveplot) then begin
; current plotting device.
   mydevice = !D.NAME
   SET_PLOT, 'PS'
   DEVICE, FILENAME='/sps/lsst/users/rcecile/Fig/plot_ps2fitbase.eps', /PORTRAIT,/COLOR,XSIZE=17.6,YSIZE=12,FONT_SIZE=4
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
xt = findgen(1000)/5000.
for iz=0,2 do begin
print,'  '
for ierr=0,2 do begin
   print,'iz, ierr  ' ,namez[iz],nsuff[ierr]
   mysuff = nsuff[ierr]
   t  =  dir + 'PS_lfZuccaAllFalse'+nx[iz]+'_z'+namez[iz]+mysuff+'_wngal.txt'
   pcorr = read_ascii(t, template =  TEMP_POW_SPEC_4col)     
   xs = pcorr.(0)
   sk ='_k0.0' +strtrim(kmin[iz],2)+'_0.'+strtrim(kmax[iz*5+ierr],2)

   fbase  = dir + 'fit_lfZuccaAllFalse'+sk+nx[iz]+'_z'+namez[iz]+mysuff+'_baseline.txt'
   print,fbase
   pbase = read_ascii(fbase, template =  TEMP_POW_SPEC_BASE)     
   ps = (pcorr.(1)-pcorr.(2))/(pbase.(2))           
   psm= (pcorr.(1)-pcorr.(2)-pcorr.(3))/(pbase.(2))
   psp= (pcorr.(1)-pcorr.(2)+pcorr.(3))/(pbase.(2)) 
   ok = where(xs ge kmin[iz]/100. and xs le kmax[5*iz+ierr]/100.)
   plot,xs,ps,/xs,/ys,xra=[0.021,0.199],yra=[0.91,1.09],ytit=myytit[ierr],xtit=myxtit[iz],xma=[myxmargin1[ierr],myxmargin2[ierr]],yma=[myymargin1[iz],myymargin2[iz]],/nodata
   oplot,[kfid,kfid],[0,2],th=2
   oplot,[0,.2],[1,1],th=2
   oplot,xs[ok],ps[ok],col=lcol[ierr],psym=lpsym[ierr]
   errplot,xs[ok],psm[ok],psp[ok],col=lcol[ierr]

  ffitw = dir + 'fit_lfZuccaAllFalse'+sk+nx[iz]+'_z'+namez[iz]+mysuff+'_result.txt'
  pfitw = read_ascii(ffitw, template =  TEMP_POW_SPEC_FIT_RES)     

  print,'        FIT with shot-noise sub',(pfitw.(0))[0],' - ',((pfitw.(2)))[0],' /+',((pfitw.(1)))[0],' (amp = ',(pfitw.(3))[0],')'
  oplot,xt,decay_func(xt,[pfitw.(0)[0],pfitw.(3)[0]]),col=lcol[ierr],th=3
;  xyouts,0.105,1.05,lerr[ierr],charsize=2
  legend,'z = '+namez[iz],charsize=2,/top,/right,/box
;  oplot,2.*!pi/[pfitw.(0)[0],pfitw.(0)[0]],[0,2],li=1,th=1
  if (iz eq 2) then legend,lerr[ierr],charsize=2,li=0,psym=lpsym[ierr],col=lcol[ierr],/bottom,/center,/box
endfor
endfor
;legend,lerr,charsize=1.5,li=1,psym=lpsym,col=lcol,/bottom,/right,/box

if (saveplot) then begin
   DEVICE, /CLOSE
   SET_PLOT, mydevice
endif

end
;  TEMP_POW_SPEC_FIT_SA = ASCII_TEMPLATE(ffitw)
; save,TEMP_POW_SPEC_FIT_SA,TEMP_SHOTNOISE,TEMP_SEL_FUNC,TEMP_POW_SPEC_ZXY,TEMP_POW_SPEC_TXT,TEMP_POW_SPEC_TH,TEMP_POW_SPEC_FITCHI2,TEMP_POW_SPEC_FIT,$
; TEMP_POW_SPEC,TEMP_PDF,TEMP_INFOS,TEMP_INFO,TEMP_2FIT,file='temp_ascii.sav' 
