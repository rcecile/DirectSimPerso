loadct,39
restore,'temp_ascii_new.sav'        ; contient dir
dir="/sps/lsst/data/rcecile/Planck_BAO2_PS/"
g = mrdfits(dir+"simu_ps_120_z0.51_r.fits",0,h)
print,h

mylev=[0,0.25,.5,1,2]*50
mycol = findgen(5)*50+30

plot,g
plot,g[0:1e4],/xs,/ys,psym=3

!p.multi=[0,1,2]
contour,g[*,*,50],/xs,/ys,lev=mylev,c_col=mycol,ytit='x-axis [cell]',/fill
contour,g[50,*,*],/xs,/ys,lev=mylev,c_col=mycol,ytit='x-axis [cell]',/fill
