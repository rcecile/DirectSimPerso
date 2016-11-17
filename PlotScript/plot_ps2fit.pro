PRO plot_ps2fit,saveplot,option,iz

; restore,'temp_ascii.sav'        ; contient dir
; f='/sps/lsst/data/rcecile/TJP_BAO_PS/fit_G2_Err_1024_z1.8_chisq.txt'
; TEMP_POW_SPEC_FITCHI2 = ASCII_TEMPLATE(f)
; save,TEMP_2FIT,TEMP_INFO,TEMP_INFOS,TEMP_PDF,TEMP_POW_SPEC,TEMP_POW_SPEC_FIT,TEMP_POW_SPEC_FITCHI2,TEMP_POW_SPEC_TH,TEMP_POW_SPEC_TXT,TEMP_SEL_FUNC,TEMP_POW_SPEC_ZXY ,file='temp_ascii.sav'

!p.charsize=2
!p.thick=2
!p.symsize=4
h= 0.679
h3 = h*h*h

loadct,12
restore,'temp_ascii.sav'        ; contient dir
dir="/sps/lsst/data/rcecile/TJP_BAO_PS/"

z=[0.5,0.9, 1.3,1.8,1.8,1.8,1.8];,1.8]
namez = ['0.5','0.9','1.3','1.8','1.8','1.8','1.8'];,'1.8']
nx= ['_350','_640', '_900', '_1024', '_1024', '_1024', '_1024'];, '_500']
gsuff = ['', '', ' thin grid', ' thick grid']
gsuff = ['', '', '', '', '', '', '']
kmax = [12, 12, 12, 12,10,8,6]
lkmax = kmax/100.
skmax = strarr(n_elements(kmax))
for i=0,n_elements(kmax)-1 do $
   if (kmax[i] ge 10) then skmax[i] = '_k0.'+strtrim(kmax[i],2) $
   else skmax[i] = '_k0.0'+strtrim(kmax[i],2) 
print,skmax
nz = n_elements(namez)

if (option eq 0) then begin
   suff=['','_err0.03','_errPpodds']
   lerr=['spectroZ','Gaussian 0.03','photoZ podds']
   lpsym=[-1,-5,-4]
   lcol  = [80,35,  135]
   wtit='_SGP'
endif

if (option eq 1) then begin
   suff=['','_errPpodds']
   lerr=['spectroZ','photoZ podds']
   lpsym=[-1,-4]
   lcol  = [105, 135]
   wtit='_SP'
endif

if (option eq 2) then begin
   suff=['','_err0.03']
   lerr=['spectroZ','Gaussian 0.03']
   lpsym=[-1,-5]
   lcol  = [80,35]
   wtit='_SG'
endif

if (option eq 3) then begin
   suff=['','_errPpodds','_errPBDT']
   lerr=['spectroZ','photoZ podds','photoZ BDT']
   lpsym=[-1,-5,-3]
   lcol  = [80,135,125]
   wtit='_SPB'
endif


if (option eq 4) then begin
   suff=['','_err0.03','_errPpodds','_errPBDT']
   lerr=['spectroZ','Gaussian 0.03','photoZ podds','photoZ BDT']
   lpsym=[-1,-5,-3,-2]
   lcol  = [80,35,135,125]
   wtit='_SGPB'
endif


if (saveplot) then begin
; current plotting device.
   mydevice = !D.NAME
   SET_PLOT, 'PS'
   DEVICE, FILENAME='/sps/lsst/dev/rcecile/Fig/plot_ps2fit'+wtit+'_'+strtrim(fix(z[iz]*10),2)+nx[iz]+'.eps', /PORTRAIT,/COLOR,XSIZE=8.8,YSIZE=4.4,FONT_SIZE=4
endif


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
!p.multi=0
;window,iz
myformat='(F10.1)'
ltext ='z = '+ namez
lline= 0
plot,[0.02,0.1],[1e4,1e4],/xs,/ys,xra=[0.02,lkmax[iz]],ytit='Normalized PS [(Mpc)^-3]',/nodata,yra=[0.75,1.25],xtit='wavenumber [Mpc^-1]',xma=[7,2],yma=[3,2]

for ir=0,n_elements(suff)-1 do begin

   mysuff = suff[ir]
  ; spectre observe corrige du damping de l'erreur 
   t  =  dir + 'PS_G2_2fitErr'+nx[iz]+'_z'+namez[iz]+mysuff+'_wngal.txt'
;   if (ir eq 2 and iz eq 2) then t  =  dir + 'PS_G2_2fitErr'+nx[iz]+'_z'+namez[iz]+'_k0.06'+mysuff+'_wngal.txt'
   pcorr = read_ascii(t, template =  TEMP_2FIT)     
   print,t
   xs = pcorr.(0)

   ; spectre th??orique no osc, no error
   ftheo  = dir + 'simu_ps'+nx[iz]+'_z'+namez[iz]+'_ntpk.txt'
   ptheo = read_ascii(ftheo, template =  TEMP_POW_SPEC_TH)     
   xt=ptheo.(0)
   ok = intarr(n_elements(xs))
   for i=0,n_elements(xs)-1 do ok[i] = where( abs(xs[i]-xt) eq min(abs(xs[i]-xt)))

   ps = (pcorr.(1)-pcorr.(6))/(ptheo.(2))[ok]
   psm= (pcorr.(1)-pcorr.(6)-pcorr.(7))/(ptheo.(2))[ok]
   psp= (pcorr.(1)-pcorr.(6)+pcorr.(7))/(ptheo.(2))[ok]

   oplot,xs+0.0003*ir,ps,col=lcol[ir],psym=lpsym[ir],th=2

   print,ir
   errplot,xs+0.0003*ir,psm,psp,col=lcol[ir]

   ;fit
   ffit = dir + 'fit_G2_Err'+skmax[iz]+nx[iz]+'_z'+namez[iz]+mysuff+'_result.txt'
   print,ffit
  ; if (ir eq 3 and iz eq 2) then ffit = dir + 'fit_G2_Err'+nx[iz]+'_z'+namez[iz]+'_k0.03'+mysuff+'_result.txt'
   pfit = read_ascii(ffit, template =  TEMP_POW_SPEC_FIT)     
   oplot,xs,decay_func(xs,pfit.(1),h),col=lcol[ir],li=2,th=5

   oplot,[pfit.(1)-pfit.(4),pfit.(1)-pfit.(4)],[0,2],col=lcol[ir],th=5
   oplot,[pfit.(1)+pfit.(3),pfit.(1)+pfit.(3)],[0,2],col=lcol[ir],th=5
   oplot,xt,ptheo.(1)/ptheo.(2),li=2

   kpt = string(pfit.(1)*1000.,format=myformat)+' +/-'+string(pfit.(3)*1000.,format=myformat)
   spt = string(2.*!pi/pfit.(1),format=myformat)+' +/-'+string(2.*!pi*pfit.(3)/pfit.(1)/pfit.(1),format=myformat)

   chi2fit = dir + 'fit_G2_Err'+skmax[iz]+nx[iz]+'_z'+namez[iz]+mysuff+'_chisq.txt'
   c2fit = read_ascii(chi2fit, template =  TEMP_POW_SPEC_FITCHI2)     
   if (ir eq 0 ) then kerr =lerr[ir] + '   k_a ='+kpt + ' chi2 = ' + strtrim(min(c2fit.(1)),2)  else $
      kerr = [kerr,' ======== ' + lerr[ir]+'   k_a ='+kpt +' chi2 = ' +strtrim( min(c2fit.(1)),2) ]
   if (ir eq 0 ) then serr =lerr[ir] + '   s_a ='+spt   else $
      serr = [serr,' ======== ' + lerr[ir]+'   s_a ='+spt   ]
 ;  read,xx
endfor
legend,lerr,col=lcol,li=0,box=1,/fill,/right,/top,charsize=1.5
legend,'z = '+namez[iz]+gsuff[iz],box=1,/fill,/right,/bottom,charsize=2
print,"========================================================= "
print,serr
print,"========================================================= "
print,kerr
print,"========================================================= "

oplot,[0,2],[1,1],th=2

if (saveplot) then begin
   DEVICE, /CLOSE
   SET_PLOT, mydevice
endif

end
