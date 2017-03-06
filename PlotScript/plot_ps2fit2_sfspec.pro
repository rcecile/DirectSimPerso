PRO plot_ps2fit2_sfspec,saveplot

!p.charsize=2
!p.thick=3
!p.symsize=4
h= 0.679
h3 = h*h*h

loadct,39
restore,'temp_ascii.sav'        ; contient dir
dir="/sps/lsst/data/rcecile/Planck_BAO_PS/"


nsuff=['_errP','_errPBDT8']
datasuff=['','_SFspec']
lerr=['photoZ','photoZ BDT 80%']
ldata=['noBAO with appropriate SF','noBAO with SF from z_s']
ltyp=[0,2]

z=[0.5,1.5]
namez = ['0.5','1.5']
namezmean=['0.51', '1.57']
nx= ['_120','_320']
nz = n_elements(namez)


if (saveplot) then begin
; current plotting device.
   mydevice = !D.NAME
   SET_PLOT, 'PS'
   DEVICE, FILENAME='/sps/lsst/dev/rcecile/Fig/plot_ps2fit_SFspec.eps', /PORTRAIT,/COLOR,XSIZE=8.8,YSIZE=6.6,FONT_SIZE=4
endif 
if (saveplot) then PAT = bytarr(5,5) else PAT = bytarr(5,5)+255


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
!p.multi=0
;window,iz
myformat='(F10.1)'
ltext ='z = '+ namez
lline= 0

for iz=0,1  do begin
   ftheo  = dir + 'simu_ps'+nx[iz]+'_z'+namezmean[iz]+'_ntpk.txt'
   pth = read_ascii(ftheo, template =  TEMP_POW_SPEC_TH)     
   xt=pth.(0)
   
   for ierr=0,1 do begin
      mysuff = nsuff[ierr]

      if (ierr+iz eq 0) then plot,[0.02,0.1],[1e4,1e4],/xs,/ys,xra=[0.02,0.12],ytit='Normalized PS [(Mpc)^-3]',/nodata,yra=[0.5,5.],xtit='wavenumber [Mpc^-1]',xma=[4.5,2],yma=[3,.5]
      
      for id=0,1 do begin
         t = dir + 'PS'+datasuff[id]+'_2fit'+nx[iz]+'_z'+namez[iz]+mysuff+'_wngal.txt' 
         pcorr = read_ascii(t, template =  TEMP_2FIT)     
         xs = pcorr.(0)
         ok = intarr(n_elements(xs))
         for i=0,n_elements(xs)-1 do ok[i] = where( abs(xs[i]-xt) eq min(abs(xs[i]-xt)))
         ptheo = (pth.(2))[ok]
         if (id eq 0) then oplot,xt,(pth.(1))/(pth.(2))+3-iz*2-ierr,col=190,th=5



         ps = (pcorr.(1)-pcorr.(6))/ ptheo
         oplot,xs,ps+3-iz*2-ierr,li=ltyp[id],th=3

         xyouts,0.095-ierr*0.01,3-iz*2-ierr+1.2,'z='+namez[iz]+', '+lerr[ierr],charsize=1.7
      endfor
   endfor
endfor
legend,ldata,line=ltyp,box=1,/fill,/left,/top,charsize=1.5

print,ldata
print,ltyp

if (saveplot) then begin
   DEVICE, /CLOSE
   SET_PLOT, mydevice
endif

end
