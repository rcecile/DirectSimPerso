PRO plot_cat

m = mrdfits("/sps/lsst/data/rcecile/Planck_BAO/cat_zOrd_Slice29.fits",1,h)

ok = where(m.(3) gt 0.88 and m.(3) lt 0.8805)
theta = (m.(2))[ok]
phi = (m.(1))[ok]
zs = (m.(3))[ok]

dc = dloscom(mean(zs))

x=dc*cos(phi)*sin(theta)            ;
y=dc*sin(phi)*sin(theta)            ;
z=dc*cos(theta)                    ;
window,1
plot,x,y,/xs,/ys,psym=3


m0 = mrdfits("/sps/lsst/data/rcecile/TJP_BAO/cat_G2_Slice32.fits",1,h)
ok = where(m0.(3) gt 0.88 and m0.(3) lt 0.8805)
help,ok
theta0 = (m0.(2))[ok]
phi0 = (m0.(1))[ok]
zs0 = (m0.(3))[ok]

dc0 = dloscom(mean(zs0))

x0=dc0*cos(phi0)*sin(theta0)            ;
y0=dc0*sin(phi0)*sin(theta0)            ;
z0=dc0*cos(theta0)                    ;
window,0
plot,x0,y0,/xs,/ys,psym=3



END

