PRO rangez_zp,isuff
; Catalog must be totaly produced before use of this routine

suff = ['_errPBDT9','_errPBDT8','_errP','_err0.03']

dir='/sps/lsst/data/rcecile/Planck_BAO2/'
nslice = 100
nslicez = 100
z_max = 3.
pinterquart = dblarr(nslicez)
pbias = dblarr(nslicez)
poutlier = dblarr(nslicez)

zslice = findgen(nslicez+1)*z_max/nslicez

for k=3,nslice-1 do begin
   names=dir+'cat_AllzOrdZP_Slice'+strtrim(k,2)+suff[isuff]+'.fits'
   nok_tot = 0

   for i=3,nslice-1 do begin

      name=dir+'cat_AllzOrd_Slice'+strtrim(i,2)+suff[isuff]+'.fits'
      print,k,i,name
      zpmin = zslice[i]   - 0.15*(1.+zslice[i])
      zpmax = zslice[i+1] + 0.15*(1.+zslice[i+1])
      print,'zp cat range  ', k,i,zpmin,zpmax,' for slice range ' ,zslice[k],'-',zslice[k+1]
      if (zpmin ge zslice[k+1]  OR zpmax lt zslice[k] ) then continue

      m=mrdfits(name ,1)
      ok = where(m.(6) ge zslice[k] and m.(6) lt zslice[k+1],nok)
      print,'TOTAL   ',nok_tot,k,i
      if (nok_tot eq 0) then mz = m[ok] else mz = [mz,m[ok]]
      nok_tot = nok_tot + nok
      print,k,i,zslice[k],nok
      
   endfor
   name=dir+'cat_AllSlice'+suff[isuff]+'_cataLT10.fits'
   check = FILE_TEST(name)
   if (check eq 1) then begin
      m=mrdfits(name ,1)
      ok = where(m.(6) ge zslice[k] and m.(6) lt zslice[k+1],nok)
      print,'TOTAL   with cata',nok_tot,k,i
      if (nok_tot eq 0) then mz = m[ok] else mz = [mz,m[ok]]
      nok_tot = nok_tot + nok
      print,k,i,zslice[k],nok
   endif

   m = 0
   if (nok_tot eq 0) then continue

   print,'Ready to sort the catalog'
   msz = mz[sort(mz.(3))]
   MWRFITS, msz, names, h
   print,'############# ', names,'  written'
   sxaddpar,h,'ZCAT_MIN',zslice[k]
   sxaddpar,h,'ZCAT_MAX',zslice[k+1]
   modfits,names,0,h,exten_no=1    ;Update header of ext 1, let data unchaged 
   print,'   header updated'
   
   t = (msz.(3)-msz.(6))/(1.+msz.(6))
   msz=0

   RSTAT,t, Med, Hinge1, Hinge2
   pbias[k] = Med
   pinterquart[k] = Hinge2-Hinge1
   out = where(abs(t) gt 0.15,nout)
   poutlier[k] = double(nout)/double(n_elements(t))*100.

   print,'As a function of z_s --------------------------------------------------'
   print,'-------------------------------------------------------------------------------------------------'
   print,'stat ',k,'       bias ',pbias[k],'    sigma ',pinterquart[k],'   outlier ',poutlier[k]
   print,'-------------------------------------------------------------------------------------------------'
   save,z,pinterquart,pbias,poutlier,file=dir+'stat_zp_zs'+suff[isuff]+'.sav'


endfor

END

