PRO check_ps_err,iz,in

; check_ps_err,iz,0 : BAO, eisenstein
; check_ps_err,iz,1 : BAO, Julien

!p.charsize=2
!p.thick=2
!p.symsize=4
h= 0.679
h3 = h*h*h

loadct,12
restore,'temp_ascii.sav'        ; contient dir
dir="/sps/lsst/data/rcecile/Planck_BAO_PS/"
dirJ="/sps/lsst/data/jsouchar/PS_fit/"
dirJfit="/sps/lsst/data/jsouchar/Planck_BAO_PS/"

if (in eq 1) then dirfit = dirJfit else dirfit = dir
if (in eq 1) then dirps = dirJ     else dirps = dir

nsuff=['','_errP','_errPBDT8']
lerr=['spectroZ','photoZ','photoZ BDT 80%']
lcol  = [95,35,120]
lpsym=[1,6, 4]

z=[0.5,0.9, 1.5]
namez = ['0.5','0.9','1.5']
namezmean=['0.51', '0.93', '1.57']
;nx= ['_160','_300', '_225']
nx= ['_120','_225', '_320']
nz = n_elements(namez)
kmax = intarr(n_elements(namez), n_elements(nsuff))
kmax[*,0] = [20,20,16]
kmax[*,1] = [20,20,14]
kmax[*,2] = [16,16,12]
skmax = strarr(n_elements(namez), n_elements(nsuff))
for i=0,n_elements(namez)-1 do for j=0,n_elements(nsuff)-1 do $
   if (kmax[i,j] ge 10) then skmax[i,j] = '_k0.'+strtrim(kmax[i,j],2) $
   else skmax[i,j] = '_k0.0'+strtrim(kmax[i,j],2) 
print,skmax


window,10*in+iz



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
!p.multi=0
;window,iz
myformat='(F10.1)'
ltext ='z = '+ namez
lline= 0

mxtit=[' ',' ','wavenumber [Mpc^-1]']
!p.multi=[0,2,3]
for ierr=0,2 do begin

   plot,[0.02,0.1],[1e4,1e4],/xs,/ys,xra=[0.02,0.20],ytit='Normalized PS [(Mpc)^-3]',/nodata,yra=[0.91,1.09],xtit=mxtit[iz],xma=[7,2],yma=[3,2]
   
   mysuff = nsuff[ierr]
   t  =  dir + 'PS_2fit_SN_Ref'+nx[iz]+'_z'+namez[iz]+mysuff+'_wngal.txt'
   if (in ge 1) then t = dir + 'PS_2fitJulienSN'+nx[iz]+'_z'+namez[iz]+mysuff+'_wngal.txt' 
   pcorr = read_ascii(t, template =  TEMP_2FIT)     
   xs = pcorr.(0)
                                ; spectre th??orique no osc, no error
   ftheo  = dirps + 'simu_ps'+nx[iz]+'_z'+namezmean[iz]+'_ntpk.txt'
   if (in ge 1) then ftheo  = dirps + 'PS_2fitJulienSN'+nx[iz]+'_z'+namez[iz]+mysuff+'_nosc_out.txt'
   ptheo = read_ascii(ftheo, template =  TEMP_POW_SPEC_TH)     
   xt=ptheo.(0)
   base = interpol((ptheo.(2)),(ptheo.(0)),pcorr.(0))

   ps = (pcorr.(1)-pcorr.(6))/base
   psm= (pcorr.(1)-pcorr.(6)-pcorr.(7))/base
   psp= (pcorr.(1)-pcorr.(6)+pcorr.(7))/base 

 ;  oplot,xs,ps,col=lcol[ierr],psym=lpsym[ierr],th=2
   errplot,xs,psm,psp,col=lcol[ierr]

   ffitw = dirfit +'fit_Ref'+skmax[iz,ierr]+nx[iz]+'_z'+namez[iz]+mysuff+'_result.txt'
   if (in ge 1) then ffitw = dirfit + 'fit_data_smoothSN'+skmax[iz,ierr]+nx[iz]+'_z'+namez[iz]+mysuff+'_result.txt'
   print,ffitw
;   TEMP_POW_SPEC_FIT_SA = ASCII_TEMPLATE(ffitw)
;   save,TEMP_2FIT,TEMP_INFO,TEMP_INFOS,TEMP_PDF,TEMP_POW_SPEC,TEMP_POW_SPEC_FIT,TEMP_POW_SPEC_FITCHI2,TEMP_POW_SPEC_TH,TEMP_POW_SPEC_TXT,TEMP_SEL_FUNC,TEMP_POW_SPEC_ZXY,TEMP_POW_SPEC_PS ,TEMP_POW_SPEC_FIT_SA,file='temp_ascii.sav'
   pfitw = read_ascii(ffitw, template =  TEMP_POW_SPEC_FIT_SA)     
   print,'        FIT with shot-noise sub',(pfitw.(1)),' - ',((pfitw.(4))),' /+',((pfitw.(3))),' (amp = ',(pfitw.(5)),')'
   fit = decay_func(xs,[pfitw.(1),pfitw.(5),pfitw.(6)])
   oplot,xs,fit,col=210,th=3 

   diff_err =  ((pcorr.(1)-pcorr.(6)) - fit*base )/ pcorr.(7)
   dok = where(pcorr.(0) ge 0.02 and pcorr.(0) le kmax[iz,ierr]/100.)
   hd = histogram(diff_err[dok],min=-5,max=5,bin=0.5)
   xhd = findgen(n_elements(hd))*0.5-5
   plot,xhd,hd,th=3
   res = gaussfit(xhd,hd,a,nterms=3,est=[5.,0.,1.])
   print,'amplitude = ',a[0],', central value = ',a[1],', sigma= ',a[2]
   gth = a[0]* exp(-0.5*xhd*xhd)
   oplot,xhd,gth,col=210,th=4
   oplot,xhd,res,col=190,li=2,th=4

   legend,'sigma = '+strtrim(a[2],2),col=190,box=1,/fill,/right,/top,charsize=1.5
endfor

;write_jpeg,'/sps/lsst/dev/rcecile/BAO_InOut/check_ps_err_z'+namez[iz]+'.jpg' ,tvrd(true=3),true=3
stop
end
;err_tot=[0.030515025 ,0.028397875, 0.037388461]
;err90_tot=[0.025927096 ,0.027347027, 0.032890236]
;err80_tot=[0.023046306 ,0.026009080, 0.030297966]

;print,err_tot/1.34900
;print,err90_tot/1.34900
;print,err80_tot/1.34900
