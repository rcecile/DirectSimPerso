PRO plot_stat_zs_zp,doplot

dir='/sps/lsst/data/rcecile/TJP_BAO/'
loadct,39
restore,dir+'/histo_zs_zp_0.01.sav'

if (doplot eq 1) then begin
   mydevice = !D.NAME
   SET_PLOT, 'PS'
   DEVICE, FILENAME='/sps/lsst/dev/rcecile/Fig/plot_good_cata.eps', /PORTRAIT,/COLOR,XSIZE=8.8,YSIZE=5.5,FONT_SIZE=4
endif


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

!p.thick=4
!p.charsize=2
diff_bin=0.001
diff_max=3.
idx=findgen(diff_max/diff_bin*2+1)
x = -diff_max + idx*diff_bin
idx1 = where( abs(x) lt 0.01)
idx3 = where( abs(x) lt 0.03)
idx5 = where( abs(x) lt 0.05)
idx9 = where( abs(x) lt 0.09)
tot=total(hpdiff)

print,' '
print,'[%] within 0.01, 0.03, 0.05, 0.09 for z_p'
print,total(hpdiff[idx1])/tot*100,total(hpdiff[idx3])/tot*100,total(hpdiff[idx5])/tot*100,total(hpdiff[idx9])/tot*100,format='(f4.0)'
print,' '
print,'[%] within 0.01, 0.03, 0.05, 0.09 for z_d compared to zp'
print,total(hddiff[idx1])/tot*100,total(hddiff[idx3])/tot*100,total(hddiff[idx5])/tot*100,total(hddiff[idx9])/tot*100,format='(f4.0)'
tot=total(hddiff)
print,' '
print,'[%] within 0.01, 0.03, 0.05, 0.09 for z_d compared to zd'
print,total(hddiff[idx1])/tot*100,total(hddiff[idx3])/tot*100,total(hddiff[idx5])/tot*100,total(hddiff[idx9])/tot*100,format='(f4.0)'

htot = findgen(n_elements(z))
fp_1 = htot
fp_3 = htot
fp_9 = htot
fpodds_1 = htot
fpodds_3 = htot
fpodds_9 = htot
for i=0,n_elements(z)-1 do htot[i] = total(hp[*,i])

for i=0,n_elements(z)-1 do fp_3[i] = total(hp[where( abs(z-z[i])/(1.+z[i]) le 0.03),i])
for i=0,n_elements(z)-1 do fp_9[i] = total(hp[where( abs(z-z[i])/(1.+z[i]) gt 0.09),i])
for i=0,n_elements(z)-1 do fpodds_3[i] = total(hd[where( abs(z-z[i])/(1.+z[i]) le 0.03),i])
for i=0,n_elements(z)-1 do fpodds_9[i] = total(hd[where( abs(z-z[i])/(1.+z[i]) gt 0.09),i])

!p.multi=0
plot,z,[0,100],/xs,/ys,xtit='Redshift',ytit='Fraction of galaxies [%]',yra=[0,100],/nodata,xra=[0.1,2.6],xma=[7,1],yma=[3,1]
oplot,z,fp_3/htot*100,col=198,li=0
oplot,z,fp_9/htot*100,col=198,li=2
oplot,z,fpodds_3/htot*100,col=240,li=0
oplot,z,fpodds_9/htot*100,col=240,li=2
mtext=['"good" photoZ','"good" photoZ, ODDS cut','"catastrophic" photoZ','"catastrophic" photoZ, ODDS cut']
legend,mtext,lin=[0,0,2,2],col=[198,240,198,240],box=1,/fill,/center,/top,charsize=1.5,th=4
oplot,[0,3],[68.27,68.27],th=1

if (doplot eq 1) then begin
   DEVICE, /CLOSE
   SET_PLOT, mydevice
endif

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
!p.thick=4
!p.charsize=2
;plot,[-.05,.05],[0,23],/xs,/ys,xtit='Delta Redshift',ytit='Fraction of galaxies [%]',/nodata,xma=[7,1],yma=[3,1]
gsigma=findgen(200)
gcentre=gsigma

for i=20,220-1,1 do begin
   x= (z-z[i])/(1.+z[i])
   y= hp[*,i]/total(hd[*,i])*100.
   ok = where(abs(x) lt 0.05)
 ;  plot,x,y,col=i,/xs,/ys,xra=[-0.05,0.05]
   yfit = GAUSSFIT(x[ok], y[ok], coeff, NTERMS=3)
   gsigma[i-20] = coeff[2]
   gcentre[i-20] = coeff[1]
  ; oplot,x,coeff[0]*exp(-0.5*(x-coeff[1])*(x-coeff[1])/coeff[2]/coeff[2]),li=1
endfor


sigz=findgen(200) ; stdev dispersion in photoz
eta = sigz
ez = eta
sigz2=findgen(200) ; stdev dispersion in photoz
eta2 = sigz
ez2 = eta
ez3 = eta
for i=20,220-1,1 do begin
   t = hd[i,*]
   t = reform(t)
   x= (z-z[i])/(1.+z[i])
   sigz[i-20] = sqrt(total(t*x*x)/total(t))
   cata = where(abs(x) gt 3.*sigz[i-20])
   eta[i-20] = total(t[cata])/total(t)*100.
   ez[i-20] = abs(total (t* x )/total(t))
   ok = where(abs(x) lt 0.15)
   sigz2[i-20] = sqrt(total(t[ok]*x[ok]*x[ok])/total(t[ok]))
   cata = where(abs(x) gt 3.*sigz2[i-20])
   eta2[i-20] = total(t[cata])/total(t)*100.
   ez2[i-20] = abs(total (t[ok]* x[ok] )/total(t[ok]))
   def = where(t gt 0)
   tt = replicate(x[def[0]],t[def[0]])
   for j=def[1],n_elements(x)-1 do if (t[j] gt 0) then tt = [tt,replicate(x[j],t[j])]
   ez3[i-20] = median(tt)
   print,i,median(tt),mean(tt)
endfor

if (doplot eq 1) then begin
   mydevice = !D.NAME
   SET_PLOT, 'PS'
   DEVICE, FILENAME='/sps/lsst/dev/rcecile/Fig/plot_zp_requirement.eps', /PORTRAIT,/COLOR,XSIZE=8.8,YSIZE=5.5,FONT_SIZE=4
endif

!p.multi=[0,1,3]
plot,z[20:219],sigz,/xs,/ys,ytit="Dispersion of z_p",xma=[10,1],yma=[0,1],yra=[0,0.15],/nodata
oplot,z[20:219],sigz,col=70
oplot,z[20:219],sigz2,col=70,li=2
oplot,[0,3],[0.05,0.05],th=2
oplot,[0,3],[0.02,0.02],li=2,th=2
plot,z[20:219],eta,/xs,/ys,ytit="Catastrophic [%]",xma=[10,1],yma=[0,0],yra=[0,12],/nodata
oplot,z[20:219],eta,col=75
oplot,z[20:219],eta2,col=75,li=2
oplot,[0,3],[10,10],th=2
plot,z[20:219],ez,/xs,/ys,ytit="Bias of z_p",xma=[10,1],yma=[4,0],yra=[0.00005,0.05],/nodata,/yl,xtit='Spectroscopic redshift'
oplot,z[20:219],ez,col=80
oplot,z[20:219],ez2,col=80,li=2
oplot,z[20:219],ez3,col=0,li=2
oplot,[0,3],[0.003,0.003],th=2

if (doplot eq 1) then begin
   DEVICE, /CLOSE
   SET_PLOT, mydevice
endif


END

