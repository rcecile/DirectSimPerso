PRO stat_basique

dir='/sps/lsst/data/rcecile/Planck_BAO/'

n= 0ll
for i=3,99 do begin
   h = headfits(dir+'cat_zOrd_Slice'+strtrim(i,2)+'.fits',ext=1)
   ngal  = sxpar(h,'NAXIS2')
   print,i,ngal
   n = n+ ngal
endfor

print,n, n/!pi /2 / (1+cps(1.04))

END
