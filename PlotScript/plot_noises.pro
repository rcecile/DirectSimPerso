PRO plot_noises,saveplot
mycol=[80,205,240]
z = [0.5,0.9,1.3]

nx_cell= [120.,225., 320., 300.]
nz_cell= [125.,125.,175.,125.]
cell=8.
!p.charsize=1.5



restore,'temp_ascii_new.sav'        ; contient dir
dir="/sps/lsst/data/rcecile/Planck_BAO2_PS/"

namez = ['0.5','0.9','1.3']
lnamez = ['1.3','0.9','0.5']
namezmean = ['0.51','0.93','1.36']
nx=['_120','_225','_300']
lth=[2,6,10]
nz = n_elements(namez)
nk = n_elements(namek)
nameerr = ['', '_err0.03', '_errP','_errPBDT9','_errPBDT8']
lerr=['spectroZ','Gaussian 0.03','photoZ','photoZ BDT 90%','photoZ BDT 80%']
loadct,39
lcol  = [95,210,35,135,  120,0]
loadct,12

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


minmy=[0.,4]
maxmy=[0.5,0.]
nameerr = ['', '_err0.03', '_errP','_errPBDT9','_errPBDT8']
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
!p.multi=0
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
if (saveplot) then begin
; current plotting device.
   mydevice = !D.NAME
   SET_PLOT, 'PS'
   DEVICE, FILENAME='/sps/lsst/users/rcecile/Fig/plot_noises.eps', /PORTRAIT,/COLOR,XSIZE=8.8,YSIZE=5.5,FONT_SIZE=4
endif else window,0

 xbao = 2*!pi/150

!p.multi=0
for isuff=0,4 do begin
 if (isuff eq 0) then   plot,[0.01,0.1],[1e4,1e4],/xs,/ys,xra=[0.01,0.2],xtit='Wavenumber k [Mpc^-1]',ytit='Error on PS [%]',/nodata,$
        yra=[.3,12],xma=[5.5,.5],yma=[3.5,.5],/yl,/xl


 if (isuff eq 0) then oplot,[xbao,xbao],[.1,5],th=20,col=230
 xyouts,xbao*0.95,0.35,'BAO scale',ori=90

   suff = nameerr[isuff]
   for iz=2,0,-1 do begin
      t = dir + 'PS_lfZuccaAllFalse'+nx[iz]+'_z'+namez[iz]+suff+'_wngal.txt' 
      print,'READ ',t 
      p = read_ascii(t, template = TEMP_POW_SPEC_4col ) 
      xs = p.(0)
      sig = p.(3)/(p.(1)-p.(2))*100.
      print,sig
;         if (iz eq 0) then oplot, xs,sig,col=lcol[isuff],th=3,li=2 else $
;            if (iz eq 2) then oplot, xs,sig,col=lcol[isuff],th=3 else $
;               oplot, xs,sig,col=lcol[isuff],th=2,psym=4,symsi=2
      oplot, xs,sig,col=lcol[isuff],th=lth[iz]
   endfor
  ; read,x
endfor
;legend,'z='+lnamez,li=[0,0,2],psym=[-3,4,-3],box=1,col=lcol[4],/fill,/right,/top,charsize=1.5,th=3
legend,'z='+lnamez,li=[0,0,0],th=[10,6,2],box=1,col=lcol[4],/fill,/right,/top,charsize=1.5

legend,lerr,li=0,col=lcol,box=1,/fill,/left,/bottom,charsize=1.5,th=3
if (saveplot) then begin
   DEVICE, /CLOSE
   SET_PLOT, mydevice
endif

END
