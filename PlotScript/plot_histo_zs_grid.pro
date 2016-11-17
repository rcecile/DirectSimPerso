PRO plot_histo_zs_grid,doplot

if (doplot) then begin
; current plotting device.
   mydevice = !D.NAME
   SET_PLOT, 'PS'
   DEVICE, FILENAME='/sps/lsst/dev/rcecile/Fig/plot_histo_zs_grid.eps',/PORTRAIT,/COLOR,XSIZE=8.8,YSIZE=6.6,FONT_SIZE=4
endif 

dir='/sps/lsst/data/rcecile/TJP_BAO/'
loadct,39
!p.thick=3
!p.charsize=1.5
mycol=[50,95,210,250]

suff = ['','_errPpodds','_errPBDT','_err0.03']
z_min=[0.35,0.7,1.15,1.6];,1.45]
z_max=[0.68,1.1,1.45,2.];,2.2]
isuff = 1
!p.multi=[0,2,2]
restore,dir+'histo_zs_zp'+suff[1]+'.sav'
all_histp = all_hist
restore,dir+'histo_zs_zp'+suff[3]+'.sav'
all_histg = all_hist
for ig=0,n_elements(z_min)-1 do begin
   hg = z*0.
   hp = z*0.
   ok = where(z ge z_min[ig] and z le z_max[ig],nok,comp=okzs)
   for iok=0,nok-1 do begin
      hg = hg +  all_histg[*,ok[iok]]
      hp = hp +  all_histp[*,ok[iok]]
   endfor
   if (ig lt 2) then myma = [1,3] else myma=[3,1]
   if (ig lt 2) then mxti = '' else mxti = 'Spectroscopic redshift'
   if (ig mod 2 eq 0) then mxma = [7,1] else mxma=[3,4]
   if (ig mod 2 eq 0) then myti = 'Ngal' else myti = ''
   plot,z,hp,/xs,/ys,/yl,ymar=myma,xma=mxma,yra=[1.1,70e6],xtit=mxti,ytit=myti
   oplot,z,hg,li=2
   oplot,[z_min[ig],z_min[ig]],[0.1,1e8],col=mycol[ig],th=4
   oplot,[z_max[ig],z_max[ig]],[0.1,1e8],col=mycol[ig],th=4
endfor

if (doplot) then begin
   DEVICE, /CLOSE
   SET_PLOT, mydevice
endif

END
