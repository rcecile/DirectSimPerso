PRO plot_ps_osc,saveplot

nsuff=['','_err0.03','_errP','_errPBDT9','_errPBDT8']
prod = '_nZ'

!p.charsize=2
!p.thick=3
loadct,12
restore,'temp_ascii.sav'        ; contient dir
dir="/sps/lsst/data/rcecile/Planck_BAO2_PS/"

namez = ['0.5','0.9','1.5']
namezmean=['0.51', '0.93', '1.57']
nx=['_120','_225','_320']

nline=[2,2,0,0];,0,0]
nz = n_elements(namez)
nk = n_elements(namek)
ngal_mean=[15.66,4.65,0.84,  15.53,4.66,0.84,  16.31,4.82,0.85,  14.20,4.59,0.73,   12.23,4.29,0.61]
SN = 1./ngal_mean*8.*8.*8.

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
!p.multi=0

if (saveplot eq 1) then begin
; current plotting device.
   mydevice = !D.NAME
   SET_PLOT, 'PS'
   DEVICE, FILENAME='/sps/lsst/dev/rcecile/Fig/plot_ps_wosc'+prod+'.eps', /PORTRAIT,/COLOR,XSIZE=8.8,YSIZE=5.5,FONT_SIZE=4
endif

lcol  = [0,95,35,  135, 120]
lcol  = [0,95,210,35,135,  120]


!p.multi=0
plot,[0.01,0.1],[1e4,1e4],/xs,/ys,xra=[0.0015,0.19],/yl,xtit='wavenumber [Mpc^-1]',ytit='Power spectrum [(Mpc)^-3]',/nodata,yra=[2,5e4],xma=[9,1],yma=[3,1]
plot,[0.01,0.1],[1e4,1e4],/xs,/ys,xra=[0.0015,0.19],/yl,xtit='wavenumber [Mpc^-1]',ytit='Power spectrum [(Mpc)^-3]',/nodata,yra=[800,5e4],xma=[9,1],yma=[3,1]

for iz=0,nz-1 do begin


   for isuff=0,n_elements(nsuff)-1 do begin

      suff = nsuff[isuff]

      t  = dir + 'PS'+prod+nx[iz]+'_z'+namez[iz]+suff+'_G0_wngal.txt' 
      print,t
      check = FILE_TEST(t)
      if (check eq 0 ) then continue
;     TEMP_POW_SPEC_TXT2 = ASCII_TEMPLATE(t)
 ;    save,TEMP_2FIT,TEMP_INFO,TEMP_INFOS,TEMP_PDF,TEMP_POW_SPEC,TEMP_POW_SPEC_FIT,TEMP_POW_SPEC_FITCHI2,TEMP_POW_SPEC_TH,TEMP_POW_SPEC_TXT,TEMP_SEL_FUNC,TEMP_POW_SPEC_ZXY,TEMP_POW_SPEC_FIT_SA,TEMP_POW_SPEC_TXT2,file='temp_ascii.sav'
      p = read_ascii(t, template =  TEMP_POW_SPEC_TXT2) 
      xs = p.(0)
      pref = p.(1)-mean((p.(3))[where(p.(0) gt 0.2)])
      print,suff,xs[10],pref[12]

      

      oplot,xs,pref,col=lcol[isuff+1]
      print,iz,mean((p.(3))[where(p.(0) gt 0.2)]),mean((p.(3))[where(p.(0) gt 0.2)])/SN[isuff*3 +iz]
   ;   errplot,xs,pref-p.(4),pref+p.(4),col=lcol[isuff+1]
   ;   oplot,xs,p.(1)-SN[isuff*3 +iz],col=lcol[isuff+1],li=2
   endfor
   
   
   t  = dir + 'simu_ps'+nx[iz]+'_z'+namezmean[iz]+'_ntpk.txt'
   ps = read_ascii(t, template =  TEMP_POW_SPEC_TH)
   print,(ps.(1))[10:12]
   xt=ps.(0)
   pt=ps.(1)
   oplot,xt,pt,li=2


endfor   

mytext=['spectroZ','Gauss 0.03','photoZ','photoZ with BDT 90% cut','photoZ with BDT 80% cut']
mytext = ['fiducial model with BAO',mytext]
legend,mytext,line=[2,0,0,0,0,0],col=lcol,box=1,/fill,/left,/bottom,charsize=1.5
;legend,'z = '+ namez,line=indgen(nz)*2,box=1,/fill,/right,/top,charsize=1.5,th=4

if (saveplot eq 1) then begin
   DEVICE, /CLOSE
   SET_PLOT, mydevice
   print,'plot saved'
endif


end
