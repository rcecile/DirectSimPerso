PRO check_ps_err,iz

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
;nx= ['_160','_300', '_225']
nx= ['_120','_225', '_320']
nz = n_elements(namez)

window,iz



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
!p.multi=0
;window,iz
myformat='(F10.1)'
ltext ='z = '+ namez
lline= 0

mxtit=[' ',' ','wavenumber [Mpc^-1]']
!p.multi=[0,2,3]
for ierr=0,2 do begin

   plot,[0.02,0.1],[1e4,1e4],/xs,/ys,xra=[0.02,0.20],ytit='Normalized PS [(Mpc)^-3]',/nodata,yra=[0.5,1.5],xtit=mxtit[iz],xma=[7,2],yma=[3,2]
   
   mysuff = nsuff[ierr]
   t  =  dir + 'PS_2fit5grids'+nx[iz]+'_z'+namez[iz]+mysuff+'_wngal.txt'
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
   psm= (pcorr.(1)-pcorr.(6)-pcorr.(7))/(ptheo.(2))[ok]
   psp= (pcorr.(1)-pcorr.(6)+pcorr.(7))/(ptheo.(2))[ok] 
   
   oplot,xs,ps,col=lcol[ierr],psym=lpsym[ierr],th=2
   errplot,xs,psm,psp,col=lcol[ierr]
   oplot,xs,ptheo.(1)[ok]/ptheo.(2)[ok],col=210,th=3

   diff_err =  (pcorr.(1)-pcorr.(6) -ptheo.(1)[ok] )/pcorr.(7)
   dok = where(pcorr.(0) ge 0.02 and pcorr.(0) le 0.2)
   hd = histogram(diff_err[dok],min=-5,max=5,bin=0.1)
   xhd = findgen(n_elements(hd))*0.1-5
   plot,xhd,hd,th=3
   res = gaussfit(xhd,hd,a,nterms=3,est=[5.,0.,1.])
   print,'amplitude = ',a[0],', central value = ',a[1],', sigma= ',a[2]

   gth = a[0]* exp(-0.5*xhd*xhd)
   oplot,xhd,gth,col=210,th=4
   oplot,xhd,res,col=190,li=2,th=4
   legend,'sigma = '+strtrim(a[2],2),col=190,box=1,/fill,/right,/top,charsize=1.5
endfor

write_jpeg,'/sps/lsst/dev/rcecile/BAO_InOut/check_ps_err_z'+namez[iz]+'.jpg' ,tvrd(true=3),true=3
stop
end
;err_tot=[0.030515025 ,0.028397875, 0.037388461]
;err90_tot=[0.025927096 ,0.027347027, 0.032890236]
;err80_tot=[0.023046306 ,0.026009080, 0.030297966]

;print,err_tot/1.34900
;print,err90_tot/1.34900
;print,err80_tot/1.34900
