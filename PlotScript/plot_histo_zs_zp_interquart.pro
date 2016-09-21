PRO plot_histo_zs_zp_interquart,saveplot

dir='/sps/lsst/data/rcecile/TJP_BAO/'
restore,dir+'histo_zs_zp_requirement.sav'
loadct,39
!p.thick=3
!p.charsize=3

if (saveplot) then begin
; current plotting device.
   mydevice = !D.NAME
   SET_PLOT, 'PS'
   DEVICE, FILENAME='/sps/lsst/dev/rcecile/Fig/plot_stat_interquartile.eps', /PORTRAIT,/COLOR,XSIZE=8.8,YSIZE=4.4,FONT_SIZE=4
endif


!p.multi=[0,1,3]
plot,zstat,pinterquart*100,/xs,/ys,yma=[0,0.2],yra=[0,0.1]*100,/nodata,ytit="Sigma [x100]",xma=[6,0.5]
oplot,zstat,pinterquart*100,col=50
oplot,zstat,dinterquart*100,col=123
oplot,[0,3],[0.02,0.02]*100,li=2
oplot,[0,3],[0.05,0.05]*100,li=2

plot,zstat,poutlier,/xs,/ys,yma=[0,0],yra=[0,20],/nodata,ytit="Outliers [%]",xma=[6,0.5]
oplot,zstat,poutlier,col=50
oplot,zstat,doutlier,col=123
oplot,[0,3],[10.,10.],li=2


plot,zstat,pbias*1000,/xs,/ys,yma=[3,0],yra=[0,0.015]*1000,/nodata,ytit="Bias [x1000]",xtit="Spectroscopic redshift z_s",xma=[6,0.5]
oplot,zstat,pbias*1000,col=50
oplot,zstat,dbias*1000,col=123
oplot,[0,3],[0.003,0.003]*1000,li=2



if (saveplot) then begin
   DEVICE, /CLOSE
   SET_PLOT, mydevice
endif


END