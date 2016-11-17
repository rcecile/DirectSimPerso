PRO test_nocata,doplot

!p.charsize=2
!p.thick=3
!p.symsize=4
h= 0.679
h3 = h*h*h

loadct,39
mycol=[50,95,210,250]

restore,'temp_ascii.sav'        ; contient dir

err=['spectroZ','gauss 0.03','Photo-z','Photo-z PODDS']
err='_errPpodds'
print,err

dir='/sps/lsst/data/rcecile/TJP_BAO_PS/'
suff='G2_'

D=['0.5','0.9','1.3','1.8'];,'1.8']
snx=['350','640','900','1024'];,'500']

if (doplot) then begin
; current plotting device.
   mydevice = !D.NAME
   SET_PLOT, 'PS'
   DEVICE, FILENAME='/sps/lsst/dev/rcecile/Fig/plot_PS_gaussSF.eps',/PORTRAIT,/COLOR,XSIZE=8.8,YSIZE=4.4,FONT_SIZE=4
endif 

plot,[0.04,0.2],[0.,20],/nodata,/xs,/ys,xra=[0.01,0.15],yra=[0.5,7.5],xtit='Wavenumber [Mpc^-1]',ytit='Power spectra ratio',xma=[5,.5],yma=[3.5,.5]
for iz=0,n_elements(D)-1 do begin
   name=dir+"PS_gauss_"+snx[iz]+"_z"+D[iz]+err+"_wngal.txt"
   print,name
   result = FILE_TEST(name)
   if (result eq 0) then continue
   p1 = read_ascii(name, template =  TEMP_2FIT)  

   name=dir+"PS_specz_"+snx[iz]+"_z"+D[iz]+err+"_wngal.txt"
   print,name
   result = FILE_TEST(name)
   if (result eq 0) then continue

   p2 = read_ascii(name, template =  TEMP_2FIT)  
    
   name=dir+"PS_G2_"+snx[iz]+"_z"+D[iz]+err+"_wngal.txt"
   print,(p1.(1))[5:10]
   print,name
   p0 = read_ascii(name, template =  TEMP_2FIT)     

   xs = p0.(0)
   oplot,xs,(p1.(1)-p1.(6))/(p0.(1)-p0.(6)),col=mycol[iz]
;   oplot,xs,(p2.(1)-p2.(6))/(p0.(1)-p0.(6)),col=mycol[iz],li=2
;   read,xx
endfor
legend,'z = '+D,box=1,col=mycol,line=0,/fill,/right,/top,charsize=1.5

if (doplot) then begin
   DEVICE, /CLOSE
   SET_PLOT, mydevice
endif
;read,xx

if (doplot) then begin
; current plotting device.
   mydevice = !D.NAME
   SET_PLOT, 'PS'
   DEVICE, FILENAME='/sps/lsst/dev/rcecile/Fig/plot_PS_noCata.eps',/PORTRAIT,/COLOR,XSIZE=8.8,YSIZE=4.4,FONT_SIZE=4
endif 

plot,[0.04,0.2],[0.,20],/nodata,/xs,/ys,xra=[0.01,0.15],yra=[0.5,99],xtit='Wavenumber [Mpc^-1]',ytit='Power spectra ratio',/yl,xma=[6,2],yma=[3.5,.5]
for iz=0,n_elements(D)-1 do begin

   name=dir+"PS_G2_cata015_"+snx[iz]+"_z"+D[iz]+err+"_wngal.txt"
   print,name
   result = FILE_TEST(name)
   if (result eq 0) then continue

   p1 = read_ascii(name, template =  TEMP_2FIT)
     
   name=dir+"PS_G2_"+snx[iz]+"_z"+D[iz]+err+"_wngal.txt"
   print,(p1.(1))[5:10]
   print,name
   p0 = read_ascii(name, template =  TEMP_2FIT)     
   xs = p0.(0)
   
   oplot,xs,(p1.(1)-p1.(6))/(p0.(1)-p0.(6)),col=mycol[iz]
endfor
legend,'z = '+D,box=1,col=mycol,line=0,/fill,/left,/top,charsize=1.5

if (doplot) then begin
   DEVICE, /CLOSE
   SET_PLOT, mydevice
endif

END
