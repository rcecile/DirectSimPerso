PRO test_sf,doplot

!p.charsize=2
!p.thick=3
!p.symsize=4
h= 0.679
h3 = h*h*h


loadct,12
lcol  = [35, 135,120]

restore,'temp_ascii.sav'        ; contient dir

nameerr=['Photo-z','Photo-z BDT 90%','Photo-z BDT 80%']
errlist=['_errP','_errPBDT9','_errPBDT8']

dir='/sps/lsst/data/rcecile/Planck_BAO_PS/'
suff='noCata_'

D=['0.5','0.9'];,'1.5']
snx=['160','300','225']

!p.multi=[0,1,2]

if (doplot) then begin
; current plotting device.
   mydevice = !D.NAME
   SET_PLOT, 'PS'
   DEVICE, FILENAME='/sps/lsst/dev/rcecile/Fig/plot_PS_SFspec.eps',/PORTRAIT,/COLOR,XSIZE=8.8,YSIZE=6.6,FONT_SIZE=4
endif 

for iz=0,n_elements(D)-1 do begin

   if (iz eq 0) then plot,[0.04,0.2],[0.,20],/nodata,/xs,/ys,xra=[0.01,0.15],xticknam=replicate(" ",10),ytit='PS excess [%]',xma=[8,.5],yma=[0,.5],yra=[-20,200] $
   else plot,[0.04,0.2],[0.,20],/nodata,/xs,/ys,xra=[0.01,0.15],yra=[-20,200],xtit='Wavenumber [Mpc^-1]',ytit='PS excess [%]',xma=[8,.5],yma=[3,0]

   for ir =0,n_elements(errlist)-1 do begin

      name=dir+"PS_spec_"+snx[iz]+"_z"+D[iz]+errlist[ir]+"_wngal.txt"
      p1 = read_ascii(name, template =  TEMP_2FIT)
      
      name=dir+"PS_cube_"+snx[iz]+"_z"+D[iz]+errlist[ir]+"_wngal.txt"
      p0 = read_ascii(name, template =  TEMP_2FIT)     
      xs = p0.(0)
      
      oplot,xs,(p1.(1)-p0.(1))/(p0.(1)-p0.(6))*100.,col=lcol[ir]
      ok = where(xs ge 0.02 and xs le 0.15)
      print,ir,mean(((p1.(1)-p0.(1))/(p0.(1)-p0.(6))*100.)[ok])
   endfor

   xyouts,0.05,9,'z = '+D[iz],charsize=2
   print,'grille suivante'
endfor
mytext=['photoZ','photoZ with BDT 90% cut','photoZ with BDT 80% cut']
legend,mytext,line=0,col=lcol,box=1,/fill,/right,/top,charsize=1.5
;
if (doplot) then begin
   DEVICE, /CLOSE
   SET_PLOT, mydevice
endif


END
