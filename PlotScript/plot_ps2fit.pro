PRO plot_ps2fit,saveplot,option,iz

!p.charsize=2
!p.thick=2
!p.symsize=4
h= 0.679
h3 = h*h*h

loadct,39
restore,'temp_ascii.sav'        ; contient dir
dir="/sps/lsst/data/rcecile/TJP_BAO_PS/"


namez = ['0.7','1.4']
z=[0.7, 1.4]
nx= ['_450', '_875']

z=[0.9, 1.3,1.8]
namez = ['0.9','1.3','1.8']
nx= ['_640', '_900', '_1000']

nz = n_elements(namez)

if (option eq 0) then begin
   suff=['','_err0.03','_errPpodds']
   lerr=['spectroZ','Gaussian 0.03','photoZ podds']
   lpsym=[-1,-5,-4]
   lcol  = [80, 150, 210]
   wtit='_SGP'
endif

if (option eq 1) then begin
   suff=['','_errPpodds']
   lerr=['spectroZ','photoZ podds']
   lpsym=[-1,-4]
   lcol  = [80, 210]
   wtit='_SP'
endif

if (option eq 2) then begin
   suff=['','_err0.03']
   lerr=['spectroZ','Gaussian 0.03']
   lpsym=[-1,-5]
   lcol  = [80, 150]
   wtit='_SG'
endif


if (saveplot) then begin
; current plotting device.
   mydevice = !D.NAME
   SET_PLOT, 'PS'
   DEVICE, FILENAME='/sps/lsst/dev/rcecile/Fig/plot_ps2fit'+wtit+'_'+strtrim(fix(z[iz]*10),2)+'.eps', /PORTRAIT,/COLOR,XSIZE=8.8,YSIZE=4.4,FONT_SIZE=4
endif


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
!p.multi=0
myformat='(F10.0)'
ltext ='z = '+ namez
lline= 0
;plot,[0.02,0.1],[1e4,1e4],/xs,/ys,xra=[0.02,0.1],ytit='[(Mpc)^-3]',/nodata,yra=[0.9,1.2],xtit='wavenumber [Mpc^-1]',tit='Normalized PS @ z = '+namez[iz],xma=[7,2],yma=[3,2]
plot,[0.02,0.1],[1e4,1e4],/xs,/ys,xra=[0.02,0.1],ytit='[(Mpc)^-3]',/nodata,yra=[0.75,1.25],xtit='wavenumber [Mpc^-1]',tit='Normalized PS @ z = '+namez[iz],xma=[7,2],yma=[3,2]

for ir=0,n_elements(suff)-1 do begin

   mysuff = suff[ir]
  ; spectre observ?? corrig?? du damping de l'erreur 
   t  =  dir + 'PS_G2_2fitErr'+nx[iz]+'_z'+namez[iz]+mysuff+'_wngal.txt'
;   t  =  dir + 'PS_G2_2fitErr'+nx[iz]+'_z'+namez[iz]+mysuff+'_wngal.txt'
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
   ffit = dir + 'fit_G2_Err_k0.15'+nx[iz]+'_z'+namez[iz]+mysuff+'_result.txt'
  ; if (ir eq 3 and iz eq 2) then ffit = dir + 'fit_G2_Err'+nx[iz]+'_z'+namez[iz]+'_k0.03'+mysuff+'_result.txt'
   pfit = read_ascii(ffit, template =  TEMP_POW_SPEC_FIT)     
   oplot,xs,decay_func(xs,pfit.(1),h),col=lcol[ir],li=2,th=5

   oplot,[pfit.(1)-pfit.(4),pfit.(1)-pfit.(4)],[0,2],col=lcol[ir],th=5
   oplot,[pfit.(1)+pfit.(3),pfit.(1)+pfit.(3)],[0,2],col=lcol[ir],th=5
   oplot,xt,ptheo.(1)/ptheo.(2),li=2

   kpt = string(pfit.(1)*1000.,format=myformat)+' +/-'+string(pfit.(3)*1000.,format=myformat)
   spt = string(2.*!pi/pfit.(1),format=myformat)+' +/-'+string(2.*!pi*pfit.(3)/pfit.(1)/pfit.(1),format=myformat)

   if (ir eq 0 ) then kerr =lerr[ir] + '   k_a ='+kpt else kerr = [kerr,lerr[ir]+'   k_a ='+kpt]
   if (ir eq 0 ) then serr =lerr[ir] + '   s_a ='+spt else serr = [serr,lerr[ir]+'   s_a ='+spt]
 ;  read,xx
endfor
legend,lerr,col=lcol,li=0,box=1,/fill,/right,/top,charsize=1.5

print,"====================================================== "
print,serr
print,"====================================================== "
print,kerr
print,"====================================================== "

oplot,[0,2],[1,1],th=2

if (saveplot) then begin
   DEVICE, /CLOSE
   SET_PLOT, mydevice
endif

end
