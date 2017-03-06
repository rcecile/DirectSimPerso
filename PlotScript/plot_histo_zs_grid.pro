PRO plot_histo_zs_grid,doplot

if (doplot) then begin
; current plotting device.
   mydevice = !D.NAME
   SET_PLOT, 'PS'
   DEVICE, FILENAME='/sps/lsst/dev/rcecile/Fig/plot_histo_zs_grid.eps',/PORTRAIT,/COLOR,XSIZE=8.8,YSIZE=6.6,FONT_SIZE=4
endif 

dir='/sps/lsst/data/rcecile/Planck_BAO/'
!p.thick=3
!p.charsize=3
lcol  = [95,35,135, 120,210]
mycol=[80,205,240]

suff = ['','_errP','_errPBDT9','_errPBDT8','_err0.03']
z_min=[0.34,0.72,1.21]
z_max=[0.73,1.27,1.95]
isuff = 1
!p.multi=[0,1,3]

mxma = [7,2]
myti = 'Ngal'

mytext=['photoZ','photoZ with BDT 90% cut','photoZ with BDT 80% cut','Gaussian 0.03']
for ig=0,n_elements(z_min)-1 do begin
   if (ig eq 0) then myma = [0,1.5]
   if (ig eq 1) then myma = [0,0]
   if (ig eq 2) then myma = [3,0]
   if (ig lt 2) then mxti = '' else mxti = 'Spectroscopic redshift'
   loadct,12
   plot,[0.2,2.45],[1.1,70e6],/xs,/ys,/yl,ymar=myma,xma=mxma,yra=[1.1,70e6],xtit=mxti,ytit=myti,/nodata

   for isuff=1,n_elements(suff)-1 do begin
      restore,dir+'histo_zs_statZp'+suff[isuff]+'.sav'
      all_histg = all_hist
      hg = z*0.
      hp = z*0.
      ok = where(z ge z_min[ig] and z le z_max[ig],nok,comp=okzs)
      for iok=0,nok-1 do begin
         hg = hg +  all_histg[*,ok[iok]]
      endfor
      oplot,z,hg,col=lcol[isuff]
   endfor
   if (ig eq 0) then legend,mytext,line=0,col=lcol[1:*],box=1,/fill,/right,/top,charsize=1.5
   loadct,39
   oplot,[z_min[ig],z_min[ig]],[0.1,1e8],th=6,col=mycol[ig]
   oplot,[z_max[ig],z_max[ig]],[0.1,1e8],th=6,col=mycol[ig]

endfor

if (doplot) then begin
   DEVICE, /CLOSE
   SET_PLOT, mydevice
endif

END
