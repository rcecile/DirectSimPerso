PRO plot_hist_zs_odds_bdt,saveplot

dir='/sps/lsst/data/rcecile/TJP_BAO/'
restore,dir+'histo_zs_zp_allstat.sav'
loadct,39
!p.thick=3
!p.charsize=3

if (saveplot) then begin
; current plotting device.
   mydevice = !D.NAME
   SET_PLOT, 'PS'
   DEVICE, FILENAME='/sps/lsst/dev/rcecile/Fig/plot_stat_interquartile.eps', /PORTRAIT,/COLOR,XSIZE=8.8,YSIZE=4.4,FONT_SIZE=4
endif

!x.margin=[6,0.5]

!p.multi=[0,1,3]
plot,zstat,pinterquart[*,0]*100,/xs,/ys,yma=[0,0.2],yra=[0,0.1]*100,/nodata,ytit="Sigma [x100]"
oplot,zstat,pinterquart[*,0]*100,col=50
;oplot,zstat,pinterquart[*,4]*100,col=123
;oplot,zstat,pinterquart[*,2]*100,col=225
oplot,[0,3],[0.02,0.02]*100,li=2
oplot,[0,3],[0.05,0.05]*100,li=2

plot,zstat,poutlier[*,0],/xs,/ys,yma=[0,0],yra=[0,20],/nodata,ytit="Outliers [%]"
oplot,zstat,poutlier[*,0],col=50
;oplot,zstat,poutlier[*,4],col=123
;oplot,zstat,poutlier[*,2],col=225
oplot,[0,3],[10.,10.],li=2


plot,zstat,pbias[*,0]*1000,/xs,/ys,yma=[3,0],yra=[0,0.015]*1000,/nodata,ytit="Bias [x1000]"$
     ,xtit="Spectroscopic redshift z_s" ,yticks=3,yminor=5
oplot,zstat,pbias[*,0]*1000,col=50
oplot,zstat,pbias[*,0]*(-1000.),col=50,li=2
;oplot,zstat,pbias[*,4]*1000,col=123
;oplot,zstat,pbias[*,4]*(-1000.),col=123,li=2
;oplot,zstat,pbias[*,2]*1000,col=225
;oplot,zstat,pbias[*,2]*(-1000.),col=225,li=2
oplot,[0,3],[0.003,0.003]*1000,li=2
oplot,[0,3],[-0.003,-0.003]*1000,li=2



if (saveplot) then begin
   DEVICE, /CLOSE
   SET_PLOT, mydevice
endif


END
