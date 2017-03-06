PRO plot_ps_osc_ref,saveplot

nsuff=['','_err0.03']
proj = '_cube'
;proj = '_shell'
proj = '' ; cube en fait

!p.charsize=2
!p.thick=3
loadct,12
restore,'temp_ascii.sav'        ; contient dir
dir="/sps/lsst/data/rcecile/Planck_BAO_PS/"

namez = ['0.5','0.9','1.5']
namezmean = ['0.51','0.93','1.57']
nx=['_120','_225','_320']

nline=[2,2,0,0];,0,0]
nz = n_elements(namez)
nk = n_elements(namek)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
!p.multi=0

if (saveplot eq 1) then begin
; current plotting device.
   mydevice = !D.NAME
   SET_PLOT, 'PS'
   DEVICE, FILENAME='/sps/lsst/dev/rcecile/Fig/plot_ps_wosc_ref'+proj+'.eps', /PORTRAIT,/COLOR,XSIZE=8.8,YSIZE=5.5,FONT_SIZE=4
endif

lcol  = [0,95,35,  135, 120]
lcol  = [0,95,210,35,135,  120]


!p.multi=0
plot,[0.01,0.1],[1e4,1e4],/xs,/ys,xra=[0.01,0.25],/yl,xtit='wavenumber [Mpc^-1]',ytit='Power spectrum [(Mpc)^-3]',/nodata,yra=[20,7e4],xma=[9,1],yma=[3,.5],/xl

for iz=0,nz-1 do begin


   for isuff=0,n_elements(nsuff)-1 do begin

      suff = nsuff[isuff]
      for igd =0,4 do begin
         t  = dir + 'PS'+proj+nx[iz]+'_z'+namez[iz]+suff+'_G'+strtrim(igd,2)+'_wngal.txt' 
         print,t
         check = FILE_TEST(t)
         if (check eq 0 ) then continue
         p = read_ascii(t, template =  TEMP_POW_SPEC_TXT) 
         xs = p.(0)
         pref = (p.(1)-p.(6))
         print,suff,xs[10],pref[10]
         oplot,xs,pref,col=lcol[isuff+1] ;,li=nline[iz];,th=lth[isuff];,psym=lpsym[isuff+1],symsize=2
      endfor
   endfor
   
   ; simu avec zmean indiscernable de celle avec z
   t  = dir + 'simu_ps'+nx[iz]+'_z'+namezmean[iz]+'_ntpk.txt'
   ps = read_ascii(t, template =  TEMP_POW_SPEC_TH)
   print,(ps.(1))[10:12]
   xt=ps.(0)
   pt=ps.(1)
   oplot,xt,pt,li=2

endfor   

mytext=['spectroscopic redshift','Gaussian 0.03 error model']
mytext = ['fiducial model with BAO',mytext]
legend,mytext,line=[2,0,0],col=lcol,box=1,/fill,/left,/bottom,charsize=1.5
;legend,'z = '+ namez,line=indgen(nz)*2,box=1,/fill,/right,/top,charsize=1.5,th=4

if (saveplot eq 1) then begin
   DEVICE, /CLOSE
   SET_PLOT, mydevice
   print,'plot saved'
endif


end
