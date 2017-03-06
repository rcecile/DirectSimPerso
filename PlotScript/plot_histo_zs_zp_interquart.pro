PRO plot_histo_zs_zp_interquart,saveplot,mode
; plot_histo_zs_zp_interquart,0,0 : plot as a function of z_s
; plot_histo_zs_zp_interquart,0,1 : plot as a function of z_p

dir='/sps/lsst/data/rcecile/Planck_BAO/'
loadct,12
lcol  = [0,210,105,  135, 120]
!p.thick=3
!p.charsize=4
lcol  = [35,105,210,  135, 120]
lcol  = [35,135, 120]

!x.margin=[8.5,0.5]

!p.multi=[0,1,3]

if (saveplot) then begin  
; current plotting device.
   mydevice = !D.NAME
   SET_PLOT, 'PS'
   if (mode eq 0) then DEVICE, FILENAME='/sps/lsst/dev/rcecile/Fig/plot_stat_interquartile.eps',    /PORTRAIT,/COLOR,XSIZE=8.8,YSIZE=4.4,FONT_SIZE=4
   if (mode eq 1) then DEVICE, FILENAME='/sps/lsst/dev/rcecile/Fig/plot_stat_interquartile_zp.eps', /PORTRAIT,/COLOR,XSIZE=6.6,YSIZE=6.6,FONT_SIZE=4
endif

if (mode eq 0) then begin
   sxtit = "Spectroscopic redshift" 
   restore,dir+'histo_zs_zp_errP.sav'
   p_sig = pinterquart
   p_out = poutlier
   p_bias = pbias
   
   restore,dir+'histo_zs_zp_errPBDT9.sav'
   d_sig = pinterquart
   d_out = poutlier
   d_bias = pbias
   restore,dir+'histo_zs_zp_errPBDT8.sav'
   b_sig = pinterquart
   b_out = poutlier
   b_bias = pbias
endif
   restore,dir+'histo_zs_zp_errPBDT8.sav'

if (mode eq 1) then begin
   sxtit = "Photometric redshift" 
   restore,dir+'stat_zp_1pzs_errP.sav'
   p_sig = pinterquart
   p_out = poutlier
   p_bias = pbias
   restore,dir+'stat_zp_1pzs_errPBDT9.sav'
   d_sig = pinterquart
   d_out = poutlier
   d_bias = pbias
   restore,dir+'stat_zp_1pzs_errPBDT8.sav'
   b_sig = pinterquart
   b_out = poutlier
   b_bias = pbias
endif

plot,zstat,p_sig*100,/xs,/ys,yma=[0,0.2],yra=[0,0.115]*100,/nodata,ytit="IQR x100",xra=[0.2,2.45]
oplot,zstat,p_sig*100,col=lcol[0]
oplot,zstat,d_sig*100,col=lcol[1]
oplot,zstat,b_sig*100,col=lcol[2]
oplot,[0,3],[0.02,0.02]*100,li=2
oplot,[0,3],[0.05,0.05]*100,li=2

mtext=['photoZ','photoZ, BDT 90% cut','photoZ, BDT 80% cut']
legend,mtext,lin=0,col=lcol,box=1,/fill,/left,/top,charsize=1.5


plot,zstat,p_out,/xs,/ys,yma=[0,0],yra=[0.01,80],/nodata,ytit="f_out [%]",xra=[0.2,2.45],/yl
oplot,zstat,p_out,col=lcol[0]
oplot,zstat,d_out,col=lcol[1]
oplot,zstat,b_out,col=lcol[2]
oplot,[0,3],[10.,10.],li=2


plot,zstat,(p_bias)*1000,/xs,/ys,yma=[3,0],yra=[-6.5,14.5],/nodata,ytit="b x1000",xra=[0.2,2.45]$
     ,xtit=sxtit
oplot,[0,3],[0,0],col=210
oplot,zstat,(p_bias)*1000,col=lcol[0]
oplot,zstat,(d_bias)*1000,col=lcol[1]
oplot,zstat,(b_bias)*1000,col=lcol[2]
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

openw,lun0, '/sps/lsst/data/rcecile/Planck_BAO_grids/Bias_photoZ.txt', /get_lun
for i=0,n_elements(z)-1 do  printf, lun0, zstat[i], p_bias[i]
close, lun0
free_lun, lun0


openw,lun1, '/sps/lsst/data/rcecile/Planck_BAO_grids/Bias_photoZ_BDT90.txt', /get_lun
for i=0,n_elements(z)-1 do  printf, lun1, zstat[i], d_bias[i]
close, lun1
free_lun, lun1


openw,lun2, '/sps/lsst/data/rcecile/Planck_BAO_grids/Bias_photoZ_BDT80.txt', /get_lun
for i=0,n_elements(z)-1 do  printf, lun2, zstat[i], b_bias[i]
close, lun2
free_lun, lun2


openw,lun0, '/sps/lsst/data/rcecile/Planck_BAO_grids/Sigma_photoZ.txt', /get_lun
for i=0,n_elements(z)-1 do  printf, lun0, zstat[i], p_sig[i]
close, lun0
free_lun, lun0


openw,lun1, '/sps/lsst/data/rcecile/Planck_BAO_grids/Sigma_photoZ_BDT90.txt', /get_lun
for i=0,n_elements(z)-1 do  printf, lun1, zstat[i], d_sig[i]
close, lun1
free_lun, lun1


openw,lun2, '/sps/lsst/data/rcecile/Planck_BAO_grids/Sigma_photoZ_BDT80.txt', /get_lun
for i=0,n_elements(z)-1 do  printf, lun2, zstat[i], b_sig[i]
close, lun2
free_lun, lun2

openw,lun0, '/sps/lsst/data/rcecile/Planck_BAO_grids/Outliers_photoZ.txt', /get_lun
for i=0,n_elements(z)-1 do  printf, lun0, zstat[i], p_out[i]
close, lun0
free_lun, lun0


openw,lun1, '/sps/lsst/data/rcecile/Planck_BAO_grids/Outliers_photoZ_BDT90.txt', /get_lun
for i=0,n_elements(z)-1 do  printf, lun1, zstat[i], d_out[i]
close, lun1
free_lun, lun1


openw,lun2, '/sps/lsst/data/rcecile/Planck_BAO_grids/Outliers_photoZ_BDT80.txt', /get_lun
for i=0,n_elements(z)-1 do  printf, lun2, zstat[i], b_out[i]
close, lun2
free_lun, lun2



END
