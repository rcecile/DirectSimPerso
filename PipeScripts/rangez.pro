PRO rangez
; Catalog must be totaly produced before use of this routine

dir='/sps/lsst/data/rcecile/Planck_noBAO/'
nslice = 70
nslicez = 100
z_max = 3.

zslice = findgen(nslicez+1)*z_max/nslicez

for isim=4,4 do begin
for k=51,51 do begin
   nok_tot = 0
   for i=0,nslice-1 do begin
      name=dir+'cat_gold_'+strtrim(isim,2)+'_Slice'+strtrim(i,2)+'.fits'
      h = headfits(name,ext=1)
      n = sxpar(h,'NAXIS2')
      zmin = sxpar(h,'ZCAT_MIN')
      zmax = sxpar(h,'ZCAT_MAX')
      print,'z range  ',isim, k,i,zmin,zmax,' for range ' ,zslice[k],'-',zslice[k+1]

      if (zmin ge zslice[k+1]  OR zmax lt zslice[k] ) then continue
      m=mrdfits(name ,1)
      ok = where(m.(3) ge zslice[k] and m.(3) lt zslice[k+1],nok)
      print,'TOTAL   ',nok_tot,k,i
      if (nok_tot eq 0) then mz = m[ok] else mz = [mz,m[ok]]
      nok_tot = nok_tot + nok
      print,k,i,zslice[k],nok
      
   endfor
   m = 0
   if (nok_tot eq 0) then continue

   print,'Ready to sort the catalog'
   msz = mz[sort(mz.(3))]
   names=dir+'cat_zOrdCLOSE_'+strtrim(isim,2)+'_Slice'+strtrim(k,2)+'.fits'
   MWRFITS, msz, names, h
   print,'############# ', names,'  written'
   sxaddpar,h,'ZCAT_MIN',zslice[k]
   sxaddpar,h,'ZCAT_MAX',zslice[k+1]
   modfits,names,0,h,exten_no=1    ;Update header of ext 1, let data unchaged 
   print,'   header updated'
   msz=0
   fits_close,0

endfor
endfor

END
