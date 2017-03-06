PRO cat_mpc3,doplot

loadct,12
lcol  = [95,35,135, 120]
!p.charsize=1.5
!p.thick=4
suff = ['','_errP','_errPBDT9','_errPBDT8']
nsuff = n_elements(suff)

if (doplot eq 1) then begin
   mydevice = !D.NAME
   SET_PLOT, 'PS'
   DEVICE, FILENAME='/sps/lsst/dev/rcecile/Fig/plot_stat_mpc3.eps', /PORTRAIT,/COLOR,XSIZE=8.8,YSIZE=5.5,FONT_SIZE=4
endif

nslice = 100
ntot = dblarr(nsuff)
hp = dblarr(nslice,nsuff)
z_max = 3
zslice = findgen(nslice+1)*z_max/nslice

dslice = dloscom( (shift(zslice,-1)+zslice)/2. )
dslice = dslice[0:nslice-1]
thick_slice = dloscom( shift(zslice,-1)) - dloscom(zslice)
thick_slice = thick_slice[0:nslice-1]
surf_slice = 2. * !pi * dslice * dslice * (1. - cos(60./180.*!pi))
vslice = surf_slice  * thick_slice

dir='/sps/lsst/data/rcecile/Planck_BAO/'
for isuff = 0,nsuff-1 do begin
   for i=3,nslice -1 do begin
      if (isuff eq 0) then file  = "cat_zOrd_Slice"+strtrim(i,2)+".fits" else file  = "cat_zOrdZP_Slice"+strtrim(i,2)+suff[isuff]+".fits"
      hh=headfits(dir+file,ext=1)
      n = 1.d * sxpar(hh,'NAXIS2')
      hp[i,isuff] = n / vslice[i]
      ntot[isuff] = ntot[isuff] + n
   endfor
endfor

plot,zslice,hp[*,0],/xs,/ys,xtit='Redshift',ytit='Galaxy density [Mpc^-3]',/nodata,xra=[0.2,2.45],/yl,yra=[1e-5,.1],xmar=[7,1],ymar=[3,.5]
for i=0,nsuff-1 do oplot,zslice,hp[*,i],col=lcol[i]
mytext=['spectroZ','photoZ','photoZ with BDT 90% cut','photoZ with BDT 80% cut']
legend,mytext,line=0,col=lcol,box=1,/fill,/right,/top,charsize=1.5

print,'N galaxies / arc-min 2 =', ntot / (!pi * 180. /!pi * 180. /!pi *60.*60.) 

if (doplot eq 1) then begin
   DEVICE, /CLOSE
   SET_PLOT, mydevice
endif

dirg='/sps/lsst/data/rcecile/Planck_BAO_grids/'
st = [dirg+'SelFunc_gold__specz_nofz.txt']
s = read_ascii(st, template =  TEMP_SEL_FUNC)
z = s.(0)
specz = s.(1)

!p.multi=[0,1,2]
plot,z,specz,/xs,/ys
plot,zslice,hp[*,0]


stop
end

