PRO plot_histo_zs_zp_interquart,saveplot

dir='/sps/lsst/data/rcecile/TJP_BAO/'
loadct,12
lcol  = [0,105,35,  135, 120]
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
restore,dir+'histo_zs_zp_errP.sav'
p_sig = pinterquart
p_out = poutlier
p_bias = pbias
restore,dir+'histo_zs_zp_errPpodds.sav'
d_sig = pinterquart
d_out = poutlier
d_bias = pbias
restore,dir+'histo_zs_zp_errPBDT.sav'
b_sig = pinterquart
b_out = poutlier
b_bias = pbias


;   save,zstat,z,pinterquart,pbias,poutlier,all_hist,file=dir+'histo_zs_zp'+suff[isuff]+'.sav'
plot,zstat,p_sig*100,/xs,/ys,yma=[0,0.2],yra=[0,0.1]*100,/nodata,ytit="Sigma [x100]"
oplot,zstat,p_sig*100,col=lcol[2]
oplot,zstat,d_sig*100,col=lcol[3]
oplot,zstat,b_sig*100,col=lcol[4]
oplot,[0,3],[0.02,0.02]*100,li=2
oplot,[0,3],[0.05,0.05]*100,li=2

plot,zstat,p_out,/xs,/ys,yma=[0,0],yra=[0,20],/nodata,ytit="Outliers [%]"
oplot,zstat,p_out,col=lcol[2]
oplot,zstat,d_out,col=lcol[3]
oplot,zstat,b_out,col=lcol[4]
oplot,[0,3],[10.,10.],li=2


plot,zstat,p_bias*1000,/xs,/ys,yma=[3,0],yra=[0,0.015]*1000,/nodata,ytit="Bias [x1000]"$
     ,xtit="Spectroscopic redshift z_s" ,yticks=3,yminor=5
oplot,zstat,p_bias*1000,col=lcol[2]
oplot,zstat,d_bias*1000,col=lcol[3]
oplot,zstat,b_bias*1000,col=lcol[4]
oplot,zstat,p_bias*(-1000.),col=lcol[2],li=2
oplot,zstat,d_bias*(-1000.),col=lcol[3],li=2
oplot,zstat,b_bias*(-1000.),col=lcol[4],li=2
oplot,[0,3],[0.003,0.003]*1000,li=2
oplot,[0,3],[-0.003,-0.003]*1000,li=2



if (saveplot) then begin
   DEVICE, /CLOSE
   SET_PLOT, mydevice
endif

zref=[0.9,1.3,1.8,1.8]
zmin=[0.72,1.16,1.64,1.45]
zmax=[1.1,1.45,1.97,2.21]

print,'        z      biasx1000       sigmax100      %outlier'
for iz=0,n_elements(zref)-1 do begin
   ok = where(zstat ge zmin[iz] AND zstat le zmax[iz])
   zb = mean(d_bias[ok]*1000)
   zs = mean(d_sig[ok]*100)
   zo = mean(d_out[ok])
   print,zref[iz],zb,zs,zo

endfor

stop
END
