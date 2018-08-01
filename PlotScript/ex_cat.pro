PRO ex_cat

m=mrdfits("/sps/lsst/data/rcecile/Planck_BAO2/cat_lfZuccaAllFalse_zOrd_Slice15.fits",1,h,range=[1e4,2e4])
help,m,/str

z=m.(3)
t=m.(2)
p=m.(1)

ok = where(t lt 15./180.*!pi)

m=mrdfits("/sps/lsst/data/rcecile/Planck_BAO2/cat_lfZuccaAllFalse_zOrd_Slice16.fits",1,h,range=[1e4,2e4])
z=[z0,m.(3)]
t=[t0,m.(2)]
p=[p0,m.(1)]

ok = where(m.(2) lt 15./180.*!pi)
z=[z,(m.(3))[ok]]
t=[t,(m.(2))[ok]]
p=[p,(m.(1))[ok]]

openw,lun0, '/sps/lsst/data/rcecile/Planck_BAO2/cat_mag165_theta15.txt', /get_lun
for i=0,n_elements(z)-1 do  printf, lun0, p[ok],t[ok],z[ok]
close, lun0
free_lun, lun0
   
END
