PRO self2,doplot
; doplot = 1 : plot on screen, = 2 : save .eps file

restore,'temp_ascii.sav'
!p.charsize=2.
!p.multi=0
!p.thick=3
loadct,12
lcol  = [95,210,35,135, 120]

dirg='/sps/lsst/data/rcecile/Planck_BAO2_grids/'

st = dirg+'SelFunc_All__specz_nofz.txt'
sg = dirg+'SelFunc_All_err0.03_nofz.txt'
sp = dirg+'SelFunc_All_errP_nofz.txt'
sb8 = dirg+'SelFunc_All_errPBDT8_nofz.txt'
sb9 = dirg+'SelFunc_All_errPBDT9_nofz.txt'


if (doplot eq 2) then begin
; current plotting device.
   mydevice = !D.NAME
   SET_PLOT, 'PS'
   DEVICE, FILENAME='/sps/lsst/dev/rcecile/Fig/self_G1.eps', /PORTRAIT,/COLOR,XSIZE=8.8,YSIZE=8.8,FONT_SIZE=4
endif

print,'DOPLOT ',doplot

nsmooth = 9

!p.multi=[0,2,3]
s = read_ascii(st, template =  TEMP_SEL_FUNC)
z = s.(0)
z1 = where(s.(0) le 0.5)
z2 = where(s.(0) gt 0.5 and s.(0) le 0.75)
z3 = where(s.(0) gt 0.75)
specz = s.(1)
specz[z1] = smooth((s.(1))[z1],nsmooth)
specz[z2] = smooth((s.(1))[z2],nsmooth)
specz[z3] = smooth((s.(1))[z3],nsmooth)
oplot,z,specz
if (doplot eq 1) then begin
   plot,z,s.(1)/specz,/xs,/ys,xra=[0.1,2.65],yra=[0.9,1.1]
   xyouts,0.5,1.05,'spectroZ'
endif

s = read_ascii(sg, template =  TEMP_SEL_FUNC)
gauss = s.(1)
gauss[z1] = smooth((s.(1))[z1],nsmooth)
gauss[z2] = smooth((s.(1))[z2],nsmooth)
gauss[z3] = smooth((s.(1))[z3],nsmooth)
if (doplot eq 1) then begin
   plot,z,s.(1)/gauss,/xs,/ys,xra=[0.1,2.65],yra=[0.9,1.1]
   xyouts,0.5,1.05,'Gaussian error: 0.03*(1+z)'
endif

s = read_ascii(sp, template =  TEMP_SEL_FUNC)
photz = s.(1)
photz[z1] = smooth((s.(1))[z1],nsmooth)
photz[z2] = smooth((s.(1))[z2],nsmooth)
photz[z3] = smooth((s.(1))[z3],nsmooth)
if (doplot eq 1) then begin
   plot,z,s.(1)/photz,/xs,/ys,xra=[0.1,2.65],yra=[0.9,1.1]
   xyouts,0.5,1.05,'Photometric error'
endif

s = read_ascii(sb9, template =  TEMP_SEL_FUNC)
bdt9 = s.(1)
bdt9[z1] = smooth((s.(1))[z1],nsmooth)
bdt9[z2] = smooth((s.(1))[z2],nsmooth)
bdt9[z3] = smooth((s.(1))[z3],nsmooth)
if (doplot eq 1) then begin
   plot,z,s.(1)/bdt9,/xs,/ys,xra=[0.1,2.65],yra=[0.9,1.1]
   xyouts,0.5,1.05,'Photometric error + BDT90 cut'
endif

s = read_ascii(sb8, template =  TEMP_SEL_FUNC)
bdt8 = s.(1)
bdt8[z1] = smooth((s.(1))[z1],nsmooth)
bdt8[z2] = smooth((s.(1))[z2],nsmooth)
bdt8[z3] = smooth((s.(1))[z3],nsmooth)
if (doplot eq 1) then begin
   plot,z,s.(1)/bdt8,/xs,/ys,xra=[0.1,2.65],yra=[0.9,1.1]
   xyouts,0.5,1.05,'Photometric error + BDT80 cut'
endif


!p.multi=[0,1,2]
plot,[0.1,2.65],[0.00002,1],/xs,ys=8,ytit='log10(SF)',xra=[0.2,2.45],yra=[-4,0],yma=[0,.5],/nodata,xmar=[6,5.5]
oplot,z,alog10(specz),col=lcol[0]
oplot,z,alog10(gauss),col=lcol[1],li=2
oplot,z,alog10(photz),col=lcol[2]
oplot,z,alog10(bdt9),col=lcol[3]
oplot,z,alog10(bdt8),col=lcol[4]

what=['spectroZ','Gauss 0.03','photoZ','photoZ + 90% cut','photoZ + 80% cut']
legend,what,col=lcol,line=[0,2,0,0,0],box=1,psym=0,/fill,/left,/bottom,charsize=1.5

; voir coef dans cat_mpc3.pro
AXIS, YAXIS=1, YRANGE = (!Y.CRANGE)+ alog10(0.10397683), ys=1,  YTITLE = 'log10(d) [Mpc^-3]'


plot,[0.1,2.65],[0,2],/xs,/ys,ytit='SF ratio',xtit='Spectroscopic redshift z_s',xra=[0.2,2.45],yma=[3,0],yra=[0.1,1.75],/nodata,xmar=[6,5.5]
oplot,z,smooth(specz/specz,7),col=lcol[0]
oplot,z,smooth(gauss/specz,7),col=lcol[1],li=2
oplot,z,smooth(photz/specz,7),col=lcol[2]
oplot,z,smooth(bdt9/specz,7),col=lcol[3]
oplot,z,smooth(bdt8/specz,7),col=lcol[4]

if (doplot eq 2) then begin
   DEVICE, /CLOSE
   SET_PLOT, mydevice
endif


;stop
   
openw,lun0, '/sps/lsst/data/rcecile/Planck_BAO2_grids/SelFunc_speczmean_nofz.txt', /get_lun
for i=0,n_elements(z)-1 do  printf, lun0, z[i], specz[i]
close, lun0
free_lun, lun0
   
openw,lun1, '/sps/lsst/data/rcecile/Planck_BAO2_grids/SelFunc_gaussmean_nofz.txt', /get_lun
for i=0,n_elements(z)-1 do  printf, lun1, z[i], gauss[i]
close, lun1
free_lun, lun1
   
openw,lun2, '/sps/lsst/data/rcecile/Planck_BAO2_grids/SelFunc_photzmean_nofz.txt', /get_lun
for i=0,n_elements(z)-1 do  printf, lun2, z[i], photz[i]
close, lun2
free_lun, lun2
   
openw,lun3, '/sps/lsst/data/rcecile/Planck_BAO2_grids/SelFunc_pbdt9mean_nofz.txt', /get_lun
for i=0,n_elements(z)-1 do  printf, lun3, z[i], bdt9[i]
close, lun3
free_lun, lun3

   
openw,lun4, '/sps/lsst/data/rcecile/Planck_BAO2_grids/SelFunc_pbdt8mean_nofz.txt', /get_lun
for i=0,n_elements(z)-1 do  printf, lun4, z[i], bdt8[i]
close, lun4
free_lun, lun4

if (doplot eq 2) then stop

print,'1 if you want to check the written selection functions, q if not'
read,xx

restore,'temp_ascii.sav'
!p.multi=0
s0 = read_ascii('/sps/lsst/data/rcecile/Planck_BAO2_grids/SelFunc_speczmean_nofz.txt', template =  TEMP_SEL_FUNC)
s1 = read_ascii('/sps/lsst/data/rcecile/Planck_BAO2_grids/SelFunc_gaussmean_nofz.txt', template =  TEMP_SEL_FUNC)
s2 = read_ascii('/sps/lsst/data/rcecile/Planck_BAO2_grids/SelFunc_photzmean_nofz.txt', template =  TEMP_SEL_FUNC)
s3 = read_ascii('/sps/lsst/data/rcecile/Planck_BAO2_grids/SelFunc_pbdt9mean_nofz.txt', template =  TEMP_SEL_FUNC)
s4 = read_ascii('/sps/lsst/data/rcecile/Planck_BAO2_grids/SelFunc_pbdt8mean_nofz.txt', template =  TEMP_SEL_FUNC)
plot,[0.1,2.65],[0.0002,1],/xs,/ys,ytit='Mean selection function',/yl,xra=[0.1,3.],yra=[0.00002,1],/nodata
oplot,s0.(0),s0.(1),col=lcol[0]
oplot,s1.(0),s1.(1),col=lcol[1]
oplot,s2.(0),s2.(1),col=lcol[2]
oplot,s3.(0),s3.(1),col=lcol[3]
oplot,s4.(0),s4.(1),col=lcol[4]

s = read_ascii('/sps/lsst/data/rcecile/Planck_BAO2_grids/SelFunc_All__ZONLY.txt', template =  TEMP_SEL_FUNC)
plot,s.(0),s.(1),/xs,/ys,/yl,xra=[0.2,2.45],yra=[1e4,3e7]
oplot,s.(0),s.(1)*s0.(1),col=123
ok = where(s.(0) ge 0.2 and s.(0) lt 2.45)
ntot= s.(1)
ngold = s.(1)*s0.(1)
print,total(ntot[ok]),total(ngold[ok])
stop

end
