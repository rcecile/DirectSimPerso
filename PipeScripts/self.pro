PRO self,doplot,isuff
; doplot = 1 : plot on screen, = 2 : save .eps file

restore,'../PlotScript/temp_ascii_new.sav'
!p.charsize=2.
!p.multi=0
!p.thick=3
loadct,12
lcol  = [95,210,35,135, 120]

dirg='/sps/lsst/data/rcecile/Planck_BAO2_grids/'
suff=['All','lfZuccaAllFalse']


sz = dirg+'SelFunc_'+suff[isuff]+'_ZONLY.txt'
st = dirg+'SelFunc_'+suff[isuff]+'_specz_nofz.txt'
sg = dirg+'SelFunc_'+suff[isuff]+'_err0.03_nofz.txt'
sp = dirg+'SelFunc_'+suff[isuff]+'_errP_nofz.txt'
sb8 = dirg+'SelFunc_'+suff[isuff]+'_errPBDT8_nofz.txt'
sb9 = dirg+'SelFunc_'+suff[isuff]+'_errPBDT9_nofz.txt'


fst = read_ascii(st, template =  TEMP_SEL_FUNC)
fsz = read_ascii(sz, template =  TEMP_SEL_FUNC)
z = fst.(0)

plot,z,fst.(1),/xs,/ys,xra=[0.1,2.65],/yl
oplot,z,smooth(fst.(1),27),col=123
; 27 correspond a 0.03*(1.08)/dz

if (doplot eq 2) then begin
; current plotting device.
   mydevice = !D.NAME
   SET_PLOT, 'PS'
   DEVICE, FILENAME='/sps/lsst/users/rcecile/Fig/self_G1.eps', /PORTRAIT,/COLOR,XSIZE=8.8,YSIZE=8.8,FONT_SIZE=4
endif

print,'DOPLOT ',doplot

s = read_ascii(st, template =  TEMP_SEL_FUNC)
z = s.(0)
specz = s.(1)

s = read_ascii(sg, template =  TEMP_SEL_FUNC)
gauss = s.(1)

s = read_ascii(sp, template =  TEMP_SEL_FUNC)
photz = s.(1)

s = read_ascii(sb9, template =  TEMP_SEL_FUNC)
bdt9 = s.(1)

s = read_ascii(sb8, template =  TEMP_SEL_FUNC)
bdt8 = s.(1)


!p.multi=[0,1,2]
plot,[0.1,2.65],[0.00002,1],/xs,ys=8,ytit='log10(SF)',xra=[0.2,2.45],yra=[-4,0],yma=[1,.5],/nodata,xmar=[6,5.5]
oplot,z,alog10(smooth(specz,9)),col=lcol[0]
oplot,z,alog10(smooth(gauss,9)),col=lcol[1],li=2
oplot,z,alog10(smooth(photz,9)),col=lcol[2]
oplot,z,alog10(smooth(bdt9,9)),col=lcol[3]
oplot,z,alog10(smooth(bdt8,9)),col=lcol[4]

what=['spectroZ','Gauss 0.03','photoZ','photoZ + 90% cut','photoZ + 80% cut']
legend,what,col=lcol,line=[0,2,0,0,0],box=1,psym=0,/fill,/left,/bottom,charsize=1.5

; voir coef dans cat_mpc3.pro
AXIS, YAXIS=1, YRANGE = (!Y.CRANGE)+ alog10(0.10397683), ys=1,  YTITLE = 'log10(d) [Mpc^-3]'



plot,[0.1,2.65],[0,2],/xs,/ys,ytit='SF ratio',xtit='Spectroscopic redshift z_s',xra=[0.2,2.45],yma=[3,-1],yra=[0.1,1.75],/nodata,xmar=[6,5.5],xtickname=replicate('  ',10)
oplot,z,smooth(specz/specz,9),col=lcol[0]
oplot,z,smooth(gauss/specz,9),col=lcol[1],li=2
oplot,z,smooth(photz/specz,9),col=lcol[2]
oplot,z,smooth(bdt9/specz,9),col=lcol[3]
oplot,z,smooth(bdt8/specz,9),col=lcol[4]

if (doplot eq 2) then begin
   DEVICE, /CLOSE
   SET_PLOT, mydevice
endif


   
openw,lun0, '/sps/lsst/data/rcecile/Planck_BAO2_grids/SelFunc_'+suff[isuff]+'_speczmean_nofz.txt', /get_lun
for i=0,n_elements(z)-1 do  printf, lun0, z[i], specz[i]
close, lun0
free_lun, lun0
   
openw,lun1, '/sps/lsst/data/rcecile/Planck_BAO2_grids/SelFunc_'+suff[isuff]+'_gaussmean_nofz.txt', /get_lun
for i=0,n_elements(z)-1 do  printf, lun1, z[i], gauss[i]
close, lun1
free_lun, lun1
   
openw,lun2, '/sps/lsst/data/rcecile/Planck_BAO2_grids/SelFunc_'+suff[isuff]+'_photzmean_nofz.txt', /get_lun
for i=0,n_elements(z)-1 do  printf, lun2, z[i], photz[i]
close, lun2
free_lun, lun2
   
openw,lun3, '/sps/lsst/data/rcecile/Planck_BAO2_grids/SelFunc_'+suff[isuff]+'_pbdt9mean_nofz.txt', /get_lun
for i=0,n_elements(z)-1 do  printf, lun3, z[i], bdt9[i]
close, lun3
free_lun, lun3

   
openw,lun4, '/sps/lsst/data/rcecile/Planck_BAO2_grids/SelFunc_'+suff[isuff]+'_pbdt8mean_nofz.txt', /get_lun
for i=0,n_elements(z)-1 do  printf, lun4, z[i], bdt8[i]
close, lun4
free_lun, lun4

if (doplot eq 2) then stop

print,'1 if you want to check the written selection functions, q if not'
read,xx

restore,'temp_ascii_new.sav'
!p.multi=0
s0 = read_ascii('/sps/lsst/data/rcecile/Planck_BAO2_grids/SelFunc_'+suff[isuff]+'_speczmean_nofz.txt', template =  TEMP_SEL_FUNC)
s1 = read_ascii('/sps/lsst/data/rcecile/Planck_BAO2_grids/SelFunc_'+suff[isuff]+'_gaussmean_nofz.txt', template =  TEMP_SEL_FUNC)
s2 = read_ascii('/sps/lsst/data/rcecile/Planck_BAO2_grids/SelFunc_'+suff[isuff]+'_photzmean_nofz.txt', template =  TEMP_SEL_FUNC)
s3 = read_ascii('/sps/lsst/data/rcecile/Planck_BAO2_grids/SelFunc_'+suff[isuff]+'_pbdt9mean_nofz.txt', template =  TEMP_SEL_FUNC)
s4 = read_ascii('/sps/lsst/data/rcecile/Planck_BAO2_grids/SelFunc_'+suff[isuff]+'_pbdt8mean_nofz.txt', template =  TEMP_SEL_FUNC)
plot,[0.1,2.65],[0.0002,1],/xs,/ys,ytit='Mean selection function',/yl,xra=[0.1,3.],yra=[0.00002,1],/nodata
oplot,s0.(0),s0.(1),col=lcol[0]
oplot,s1.(0),s1.(1),col=lcol[1]
oplot,s2.(0),s2.(1),col=lcol[2]
oplot,s3.(0),s3.(1),col=lcol[3]
oplot,s4.(0),s4.(1),col=lcol[4]

read,xx

s = read_ascii('/sps/lsst/data/rcecile/Planck_BAO2_grids/SelFunc_'+suff[isuff]+'_ZONLY.txt', template =  TEMP_SEL_FUNC)
plot,s.(0),s.(1),/xs,/ys,/yl,xra=[0.2,2.45],yra=[1e4,3e7]
oplot,s.(0),s.(1)*s0.(1),col=123
ok = where(s.(0) ge 0.2 and s.(0) lt 2.45)
ntot= s.(1)
ngold = s.(1)*s0.(1)
print,total(ntot[ok]),total(ngold[ok])
stop

end
