PRO print_stat_gal

nx=['_120','_225','_300']
ncx = [120,225,300]
ncz = [125,125,125]
namez = ['0.5','0.9','1.3']


myformat ='(F7.1, F7.2)'
print,'Table 3'
suff = ['','_errPBDT9','_errPBDT8']
next = [12,11,11,11,11]

for iz=0,2 do begin
   print,' '
   for i=0,n_elements(suff)-1 do begin
     h=headfits("/sps/lsst/data/rcecile/Planck_BAO2_grids/grids_lfZuccaAllFalse"+nx[iz]+"_z"+namez[iz]+suff[i]+"_5cubes.fits",ext=next[i])
     ngal  = sxpar(h,'NGRID',/silent)
     mys = string(ngal/1e6,1.*ngal/ncx[iz]/ncx[iz]/ncz[iz],format=myformat)
     print,'iz = ', namez[iz],'  ', suff[i],mys

   endfor
endfor


suff = ['','_err0.03','_errP','_errPBDT9','_errPBDT8']
print,' '
print,'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ '
print,'pour comp_sn'

for i=0,n_elements(suff)-1 do begin
   print,' '
   for iz=0,2 do begin
     h=headfits("/sps/lsst/data/rcecile/Planck_BAO2_grids/grids_lfZuccaAllFalse"+nx[iz]+"_z"+namez[iz]+suff[i]+"_5cubes.fits",ext=next[i])
     ngal  = sxpar(h,'NGRID',/silent)

     mys = string(ngal/1e6,1.*ngal/ncx[iz]/ncx[iz]/ncz[iz],format=myformat)
     print,'iz = ', namez[iz],'  ', suff[i],mys

   endfor
endfor

END
