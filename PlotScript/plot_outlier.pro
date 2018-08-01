PRO plot_outlier,saveplot

dir='/sps/lsst/data/rcecile/Planck_BAO2/'
loadct,12
lcol  = [0,210,105,  135, 120]
!p.thick=3
!p.charsize=3
lcol  = [35,105,210,  135, 120]
lcol  = [35,135, 120]

!x.margin=[8.5,0.5]

!p.multi=[0,1,3]

if (saveplot) then begin  
; current plotting device.
   mydevice = !D.NAME
   SET_PLOT, 'PS'
   DEVICE, FILENAME='/sps/lsst/users/rcecile/Fig/plot_outlier.eps', /PORTRAIT,/COLOR,XSIZE=6.6,YSIZE=6.6,FONT_SIZE=4
endif


sxtit = "Photometric redshift" 
restore,dir+'statZp_lfZuccaAllFalse_errP.sav'
p_out15 = poutlier
p_out08 = poutlier08
p_out10 = poutlier10
p_out12 = poutlier12
restore,dir+'statZp_lfZuccaAllFalse_errPBDT9.sav'
d_out15 = poutlier
d_out08 = poutlier08
d_out10 = poutlier10
d_out12 = poutlier12
restore,dir+'statZp_lfZuccaAllFalse_errPBDT8.sav'
b_out15 = poutlier
b_out08 = poutlier08
b_out10 = poutlier10
b_out12 = poutlier12
zstat = zslice

ref=['photoZ','photoZ, BDT 90% cut','photoZ, BDT 80% cut']
mtext='e_z > ' + ['0.08','0.10','0.12','0.15']

plot,zstat,p_out08,/xs,/ys,yma=[.5,0.5],yra=[0.0,29],/nodata,ytit="f_out [%]",xra=[0.2,2.45],th=2,li=2,xma=[5.5,0.2]
oplot,zstat,p_out08,col=lcol[0],th=2,li=2
oplot,zstat,p_out10,col=lcol[0],th=2
oplot,zstat,p_out12,col=lcol[0],th=2,li=2
oplot,zstat,p_out15,col=lcol[0]
oplot,[0,3],[10.,10.],li=2

legend,mtext,lin=[2,0,2,0],thick=[2,2,2,4],box=1,/fill,/left,/top,charsize=1.
xyouts,0.9,15,ref[0],charsize=2

plot,zstat,d_out08,/xs,/ys,yma=[1.5,-.5],yra=[0.,29],/nodata,ytit="f_out [%]",xra=[0.2,2.45],th=2,li=2,xma=[5.5,0.2]
oplot,zstat,d_out08,col=lcol[1],th=2,li=2
oplot,zstat,d_out10,col=lcol[1],th=2
oplot,zstat,d_out12,col=lcol[1],th=2,li=2
oplot,zstat,d_out15,col=lcol[1]
oplot,[0,3],[10.,10.],li=2

xyouts,0.9,15,ref[1],charsize=2

plot,zstat,b_out08,/xs,/ys,yma=[3,-1.5],yra=[0.0,29],/nodata,ytit="f_out [%]",xra=[0.2,2.45],th=2,li=2,xma=[5.5,0.2],xtit='Photometric redshift'
oplot,zstat,b_out08,col=lcol[2],th=2,li=2
oplot,zstat,b_out10,col=lcol[2],th=2
oplot,zstat,b_out12,col=lcol[2],th=2,li=2
oplot,zstat,b_out15,col=lcol[2]
oplot,[0,3],[10.,10.],li=2

xyouts,0.9,15,ref[2],charsize=2



if (saveplot) then begin
   DEVICE, /CLOSE
   SET_PLOT, mydevice
endif

zref=[0.9,1.3,1.8,1.8]
zmin=[0.72,1.16,1.64,1.45]
zmax=[1.1,1.45,1.97,2.21]




END
