PRO plot_grid_slice,saveplot
nsuff=['','_G0.03','_errP','_errPBDT9','_errPBDT8']

!p.charsize=1.75

loadct,39
restore,'temp_ascii.sav'        ; contient dir
dir="/sps/lsst/data/rcecile/Planck_BAO2_grids/"

namez = ['0.5','0.9','1.5']
nx=[120,225,320]
nz=[125,125,175]

mylev=[0,0.4,.8,1.,3]
mylev=[0,0.25,.5,1,2]
mycol = findgen(5)*50+30
for iz=1,1 do begin
   if (saveplot) then begin
; current plotting device.
      mydevice = !D.NAME
      SET_PLOT, 'PS'
      DEVICE, FILENAME='/sps/lsst/dev/rcecile/Fig/plot_grid_slice'+strtrim(iz,2)+'.eps', /PORTRAIT,/COLOR,XSIZE=8.8,YSIZE=7.5,FONT_SIZE=4
   endif else window,iz+1,XSIZE=880,YSIZE=600

   g1 = mrdfits(dir+"grids_SN_"+strtrim(nx[iz],2)+"_z"+namez[iz]+"_5cubes.fits",1,h)
   g2 = mrdfits(dir+"grids_SN_"+strtrim(nx[iz],2)+"_z"+namez[iz]+nsuff[4]+"_5cubes.fits",1,h)
   !p.multi=[0,2,2]
   contour,g1[*,*,nz[iz]/2],/xs,/ys,lev=mylev,c_col=mycol,xma=[7,0.4],yma=[0.25,.5],xtickname=replicate(" ",10),ytit='x-axis [cell]',/fill
   contour,g2[*,*,nz[iz]/2],/xs,/ys,lev=mylev,c_col=mycol,xma=[0.4,7],ytickname=replicate(" ",10),yma=[0.25,.5],xtickname=replicate(" ",10),/fill
   xyouts,335,50,'transverse plane',ori=90;,charsize=1
   contour,g1[nx[iz]/2,*,*],/xs,/ys,lev=mylev,c_col=mycol,xma=[7,0.4],yma=[8,0.25],ytit='z-axis [cell]',xtit='y-axis [cell]',/fill
   xyouts,1,-80,'with spectroscopic redshift';,charsize=1
   contour,g2[nx[iz]/2,*,*],/xs,/ys,lev=mylev,c_col=mycol,xma=[0.4,7],ytickname=replicate(" ",10),yma=[8,0.25],xtit='y-axis [cell]',/fill
   xyouts,1,-80,'with photometric BDT 80% redshift';,charsize=1
   xyouts,335,1,'line-of-sight',ori=90;,charsize=1


   if (saveplot) then begin
      DEVICE, /CLOSE
      SET_PLOT, mydevice
   endif else stop

endfor

END
