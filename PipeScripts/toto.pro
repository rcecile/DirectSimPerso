g=mrdfits("/sps/lsst/data/rcecile/TJP_BAO/cat_zOrd_Slice55.fits",1,h)

theta=g.(2)
phi=g.(1)
z=g.(3)

dc = dloscom(z)
r = dc*tan(theta)
x = r*cos(phi)   
y = r*sin(phi)   

h=hist_2d(x,y,bin1=100,bin2=100)
contour,h,/xs,/ys
         100013286438          100003631382          100018275236          100000492116          100018367934          100004375554
          100007676121          100008697253          100018606976          100005940330          100006095782


g0=mrdfits("/sps/lsst/data/rcecile/TJP_BAO/cat_Slice10.fits",1,h,range=[0,10])

h = headfits("/sps/lsst/data/rcecile/TJP_BAO/cat_zOrd_Slice22.fits",ext=1)
zmin = sxpar(h,'ZCAT_MIN')
zmax = sxpar(h,'ZCAT_MAX')
 print,'z range  ', zmin,zmax
g=mrdfits("/sps/lsst/data/rcecile/TJP_BAO/cat_zOrd_Slice22.fits",1,h)

i=8
ok = where(g.(0) eq (g0.(0))[i])
print,okra
print,g0[i]
print,g[ok]


dir='/sps/lsst/data/rcecile/TJP_BAO/'
nslice=0ll
.r
for i=0,70-1 do begin
   name=dir+'cat_Slice'+strtrim(i,2)+'.fits'
   h = headfits(name,ext=1)
   n = sxpar(h,'NAXIS2')
   nslice = nslice + n
endfor
end
nsliceord=0ll
.r
for i=3,100-1 do begin
   name=dir+'cat_zOrd_Slice'+strtrim(i,2)+'.fits'
   h = headfits(name,ext=1)
   n = sxpar(h,'NAXIS2')
   nsliceord = nsliceord + n
endfor
end

print,nslice
print,nsliceord
