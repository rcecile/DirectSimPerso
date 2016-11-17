PRO plot_stat_zs_zp_mpc3,doplot

dir='/sps/lsst/data/rcecile/TJP_BAO/'
loadct,39

suff = ['_errPpodds','_errPBDT','_errP','_err0.03']
nsuff = ['_errPpodds','_errPBDT','_errP','_errG']
mytit = ['',' with ODDS cut',' with BDT cut','with Gaussian error']
 ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


loadct,12
lcol  = [0,135,120,35]

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

print,'PODDS'
print,' '
print,'[%] within 0.01, 0.03, 0.05, 0.09 for z_d compared to zp'
print,total(hddiff[idx1])/tot*100,total(hddiff[idx3])/tot*100,total(hddiff[idx5])/tot*100,total(hddiff[idx9])/tot*100,format='(f4.0)'
tot=total(hddiff)
print,' '
print,'[%] within 0.01, 0.03, 0.05, 0.09 for z_d compared to zd'
print,total(hddiff[idx1])/tot*100,total(hddiff[idx3])/tot*100,total(hddiff[idx5])/tot*100,total(hddiff[idx9])/tot*100,format='(f4.0)'

print,'BDT'
print,' '
print,'[%] within 0.01, 0.03, 0.05, 0.09 for z_b compared to zp'
print,total(hbdiff[idx1])/tot*100,total(hbdiff[idx3])/tot*100,total(hbdiff[idx5])/tot*100,total(hbdiff[idx9])/tot*100,format='(f4.0)'
tot=total(hbdiff)
print,' '
print,'[%] within 0.01, 0.03, 0.05, 0.09 for z_b compared to zb'
print,total(hbdiff[idx1])/tot*100,total(hbdiff[idx3])/tot*100,total(hbdiff[idx5])/tot*100,total(hbdiff[idx9])/tot*100,format='(f4.0)'

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
oplot,z,fp_3/htot*100,col=lcol[2],li=0
oplot,z,fp_9/htot*100,col=lcol[2],li=2
oplot,z,fpodds_3/htot*100,col=lcol[3],li=0
oplot,z,fpodds_9/htot*100,col=lcol[3],li=2
mtext=['"good" photoZ','"good" photoZ, ODDS cut','"catastrophic" photoZ','"catastrophic" photoZ, ODDS cut']
legend,mtext,lin=[0,0,2,2],col=[lcol[2],lcol[3],lcol[2],lcol[3]],box=1,/fill,/center,/top,charsize=1.5,th=4
oplot,[0,3],[68.27,68.27],th=1,col=lcol[1]

if (doplot eq 1) then begin
   DEVICE, /CLOSE
   SET_PLOT, mydevice
endif
read,xx
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;   save,hp,hd,hpdiff,hddiff,z,binz,zmax,diff_bin,diff_max,i,file=dir+'histo_zs_zp_0.01.sav'
nz = n_elements(z)
n = z
nd = z
n_mpc = z
nd_mpc = z
s_cube = 1600.*1600.*8.*8.
l_cube = 1600.*8.
diag_cube = l_cube / sqrt(2.)

theta = 1.04720
for i=0,nz-1 do begin
   n[i] = total(hp[i,*])
   nd[i] = total(hd[i,*])

   d= dloscom(z[i])
   s_disque = 2. * !pi *d*d * tan(theta)* tan(theta)
   s_coq = 2. * !pi *d*d * (1.-cos(theta))
   l_disque = d * tan(theta)
   l_coq = d

   if(l_disque le diag_cube) then surf = s_disque else begin
      alpha = acos(diag_cube/l_disque)
      s_segment = d*d/2. * (alpha- sin (alpha))
      surf = s_disque - 4.*s_segment
   endelse
;   surf = s_coq

   n_mpc[i] = n[i] / surf /  80.
   nd_mpc[i] = nd[i] / surf /  80.
   
endfor
print,'<Ngal /  arcmin2> all   = ',total(n) / (!pi * tan(1.04720)* tan(1.04720) *180./!pi *180./!pi *60.*60.)
print,'<Ngal /  arcmin2> podds   = ',total(nd) / (!pi * tan(1.04720)* tan(1.04720) *180./!pi *180./!pi *60.*60.)
;print,'<Ngal /  arcmin2> all   = ',total(n) / (!pi * 2. * (1.-cos(1.04720)) *180./!pi *180./!pi *60.*60.)
;print,'<Ngal /  arcmin2> podds   = ',total(nd) / (!pi * 2. * (1.-cos(1.04720)) *180./!pi *180./!pi *60.*60.)
print,' '
print,' '

for i=0,5 do begin
   ok = where(z gt 0.5*(i) and z lt 0.5*(i+1))
   print,'z range ',0.5*(i),'-',0.5*(i+1)
   print,'<Ngal /  arcmin2> all   = ',total(n[ok]) / (!pi * tan(1.04720)* tan(1.04720) *180./!pi *180./!pi *60.*60.)
   print,'<Ngal /  arcmin2> podds   = ',total(nd[ok]) / (!pi * tan(1.04720)* tan(1.04720) *180./!pi *180./!pi *60.*60.)
;   print,'<Ngal /  arcmin2> all   = ',total(n[ok]) / (!pi * 2. * (1.-cos(1.04720)) *180./!pi *180./!pi *60.*60.)
;   print,'<Ngal /  arcmin2> podds   = ',total(nd[ok]) / (!pi * 2. * (1.-cos(1.04720)) *180./!pi *180./!pi *60.*60.)
endfor


if (doplot eq 1) then begin
   mydevice = !D.NAME
   SET_PLOT, 'PS'
   DEVICE, FILENAME='/sps/lsst/dev/rcecile/Fig/plot_stat_mpc3.eps', /PORTRAIT,/COLOR,XSIZE=8.8,YSIZE=5.5,FONT_SIZE=4
endif

 plot,z,n_mpc,/xs,/ys,/yl,yra=[1e-6,0.5]/3.,xtit='Redshift',ytit='Ngal [/ Mpc^3]',/nodata,xra=[0.1,2.6],xma=[7,1],yma=[3,1]
oplot,z,n_mpc,col=lcol[2]
oplot,z,nd_mpc,col=lcol[3]
legend,['photoZ','photoZ with ODDS cut'],lin=[0,0],col=[lcol[2],lcol[3]],box=1,/fill,/right,/top,charsize=1.5,th=4

if (doplot eq 1) then begin
   DEVICE, /CLOSE
   SET_PLOT, mydevice
endif

END
