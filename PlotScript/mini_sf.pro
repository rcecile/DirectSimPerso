PRO mini_sf,isuff,nslice
dir='/sps/lsst/data/rcecile/Planck_BAO2/'
if nslice le 0 then nslice=70
zstat = dblarr(nslice)
nONLY=dblarr(nslice)
nGOLD=dblarr(nslice)
z=nGOLD

lsuff=['lfZuccaAllFalse','Zucca']
suff=lsuff[isuff]

for i = 0,nslice-1 do begin
   file='cat_'+suff+'_Slice'+strtrim(i,2)+'.fits'
   hh=headfits(dir+file,ext=1)
   n = sxpar(hh,'NAXIS2')
   nGOLD[i] = 1.*n
   zi = sxpar(hh,'ZCAT_MIN')
	z[i] = zi
   file='cat_'+suff+'_Slice'+strtrim(i,2)+'_ZONLY.fits'
   hh=headfits(dir+file,ext=1)
   n = sxpar(hh,'NAXIS2')
   nONLY[i] = 1.*n
   print,i,nGOLD[i],nONLY[i],nGOLD[i]/nONLY[i]*100.
endfor

plot,z,nGold,/xs,/ys
stop

end
