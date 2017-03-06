PRO plot_stat_zs_zp,doplot

dir='/sps/lsst/data/rcecile/Planck_BAO/'
suff = ['_errP','_errPBDT9','_errPBDT8']

loadct,12
lcol  = [35,135, 120]
!p.charsize=1.5
!p.thick=4

if (doplot eq 1) then begin
   mydevice = !D.NAME
   SET_PLOT, 'PS'
   DEVICE, FILENAME='/sps/lsst/dev/rcecile/Fig/plot_good_cata.eps', /PORTRAIT,/COLOR,XSIZE=8.8,YSIZE=5.5,FONT_SIZE=4
endif


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

restore,dir+'histo_zs_zp'+suff[0]+'.sav'
hp = all_hist
print,'ok'
restore,dir+'histo_zs_zp'+suff[1]+'.sav'
hd = all_hist
restore,dir+'histo_zs_zp'+suff[2]+'.sav'
hb = all_hist

htot = findgen(n_elements(z))
fp_1 = htot
fp_5 = htot
fp_15 = htot
fpodds_1 = htot
fpodds_5 = htot
fpodds_15 = htot
fbdt_1 = htot
fbdt_5 = htot
fbdt_15 = htot
for i=0,n_elements(z)-1 do htot[i] = total(hp[*,i])

for i=0,n_elements(z)-1 do fp_5[i] = total(hp[where( abs(z-z[i])/(1.+z[i]) le 0.05),i])
for i=0,n_elements(z)-1 do fp_15[i] = total(hp[where( abs(z-z[i])/(1.+z[i]) gt 0.15),i])
for i=0,n_elements(z)-1 do fpodds_5[i] = total(hd[where( abs(z-z[i])/(1.+z[i]) le 0.05),i])
for i=0,n_elements(z)-1 do fpodds_15[i] = total(hd[where( abs(z-z[i])/(1.+z[i]) gt 0.15),i])
for i=0,n_elements(z)-1 do fbdt_5[i] = total(hb[where( abs(z-z[i])/(1.+z[i]) le 0.05),i])
for i=0,n_elements(z)-1 do fbdt_15[i] = total(hb[where( abs(z-z[i])/(1.+z[i]) gt 0.15),i])

!p.multi=0
;plot,z,[0,100],/xs,/ys,xtit='Redshift',ytit='Fraction of galaxies [%]',yra=[0.1,100],/nodata,xra=[0.2,2.45],xma=[7,1],yma=[3,1],/ylog
;oplot,z,fp_5/htot*100,col=lcol[0],li=0
;oplot,z,fp_15/htot*100,col=lcol[0],li=2
;oplot,z,fpodds_5/htot*100,col=lcol[1],li=0
;oplot,z,fpodds_15/htot*100,col=lcol[1],li=2
;oplot,z,fbdt_5/htot*100,col=lcol[2],li=0
;oplot,z,fbdt_15/htot*100,col=lcol[2],li=2
;mtext=['"good" photoZ','"good" photoZ, BDT 90% cut','"good" photoZ, BDT 80% cut','"catastrophic" photoZ','"catastrophic" photoZ, BDT 90% cut','"catastrophic" photoZ, BDT 80% cut']
;legend,mtext,lin=[0,0,0,2,2,2],col=[lcol,lcol],box=1,/fill,/left,/bottom
mtext=['photoZ','photoZ, BDT 90% cut','photoZ, BDT 80% cut']
;legend,mtext,col=[lcol],box=1,/fill,/left,/center

!p.multi=[0,1,2]
plot,z,[0,100],/xs,/ys,ytit='"good"/total [%]',yra=[10,100],/nodata,xra=[0.2,2.45],xma=[8.5,1],yma=[0,0.5],xtickname=replicate(" ",10)
oplot,z,fp_5/htot*100,col=lcol[0],li=0
oplot,z,fpodds_5/htot*100,col=lcol[1],li=0
oplot,z,fbdt_5/htot*100,col=lcol[2],li=0
legend,mtext,lin=0,col=[lcol],box=1,/fill,/left,/bottom

plot,z,[0,100],/xs,/ys,xtit='Redshift',ytit='"catastrophic"/total [%]',yra=[0.01,40],/nodata,xra=[0.2,2.45],xma=[8.5,1],yma=[3,0],/yl
oplot,z,fp_15/htot*100,col=lcol[0];,li=2
oplot,z,fpodds_15/htot*100,col=lcol[1];,li=2
oplot,z,fbdt_15/htot*100,col=lcol[2];,li=2
;mtext=['"good" photoZ','"good" photoZ, BDT 90% cut','"good" photoZ, BDT 80% cut','"catastrophic" photoZ','"catastrophic" photoZ, BDT 90% cut','"catastrophic" photoZ, BDT 80% cut']
;legend,mtext,lin=[0,0,0,2,2,2],col=[lcol,lcol],box=1,/fill,/left,/bottom
mtext=['photoZ','photoZ, BDT 90% cut','photoZ, BDT 80% cut']
;legend,mtext,col=[lcol],box=1,/fill,/left,/center

if (doplot eq 1) then begin
   DEVICE, /CLOSE
   SET_PLOT, mydevice
endif

END

