PRO plot_ps_osc2,saveplot

nsuff=['','_err0.03','_errP','_errPBDT9','_errPBDT8']
proj = '_cube'
;proj = '_shell'
proj = ''

!p.charsize=4
!p.thick=3
loadct,12
restore,'temp_ascii.sav'        ; contient dir
dir="/sps/lsst/data/rcecile/Planck_BAO_PS/"

namez = ['0.5','0.9','1.5']
;namezmean = ['0.525','0.96','1.54']
;nx=['_160','_300','_225']
namezmean = ['0.51','0.93','1.57']
nx=['_120','_225','_320']

nline=[2,2,0,0];,0,0]
nz = n_elements(namez)
nk = n_elements(namek)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

if (saveplot eq 1) then begin
; current plotting device.
   mydevice = !D.NAME
   SET_PLOT, 'PS'
   DEVICE, FILENAME='/sps/lsst/dev/rcecile/Fig/plot_ps_wosc2'+proj+'.eps', /PORTRAIT,/COLOR,XSIZE=8.8,YSIZE=5.5,FONT_SIZE=4
endif

lcol  = [0,95,35,  135, 120]
lcol  = [0,95,210,35,135,  120]

ymar0=[0,0,3]
ymar1=[0.5,0,0]
xtitall=[' ',' ','wavenumber [Mpc^-1]']
mytext=['photoZ','photoZ with BDT 90% cut','photoZ with BDT 80% cut']

!p.multi=[0,1,nz]

for iz=0,nz-1 do begin

   plot,[0.01,0.1],[1e4,1e4],/xs,/ys,xra=[0.015,0.2],xtit=xtitall[iz],ytit='<PS ratio>',/nodata,yra=[0,1.9],xma=[9,1],yma=[ymar0[iz],ymar1[iz]],/xl
   
   
   for isuff=2,n_elements(nsuff)-1 do begin
      
      suff = nsuff[isuff]
      for igd=0,4 do begin

         t  = dir + 'PS'+proj+nx[iz]+'_z'+namez[iz]+nsuff[1]+'_G'+strtrim(igd,2)+'_wngal.txt' 
         pref = read_ascii(t, template =  TEMP_POW_SPEC_TXT) 
         xs = pref.(0)
         
         t  = dir + 'PS'+proj+nx[iz]+'_z'+namez[iz]+suff+'_G'+strtrim(igd,2)+'_wngal.txt' 
         p = read_ascii(t, template =  TEMP_POW_SPEC_TXT) 
         p2 = (p.(1)-p.(6))/(pref.(1)-pref.(6))
         if (igd eq 0) then pmean = p2 else pmean = pmean + p2
       ;  oplot,xs,p2,col=lcol[isuff+1] ;,li=nline[iz];,th=lth[isuff];,psym=lpsym[isuff+1],symsize=2
      endfor
      oplot,xs,pmean/5.,col=lcol[isuff+1] ;,li=nline[iz];,th=lth[isuff];,psym=lpsym[isuff+1],symsize=2
   endfor

   oplot,[0.01,0.5],[1,1],li=2
   xyouts,0.0182,0.2,'grids @ z='+namez[iz],charsize=2
   if (iz eq 1) then legend,mytext,line=0,col=lcol[3:*],box=1,/fill,/right,/bottom,charsize=1.5
endfor   


if (saveplot eq 1) then begin
   DEVICE, /CLOSE
   SET_PLOT, mydevice
   print,'plot saved'
endif


end
