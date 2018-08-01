PRO rangez_zp,isuff,mode

; mode = 0 : create zOrdZP slices + compute stat, histo
; mode = 1 : read zOrdZP slices + compute stat, histo
; mode = 2 : read zOrdZP slices + compute stat only
; Catalog must be totaly produced before use of this routine

suff = ['_errP','_errPBDT9','_errPBDT8','_err0.03']

dir='/sps/lsst/data/rcecile/Planck_BAO2/'
statname=dir+'statZp_lfZuccaAllFalse'+suff[isuff]+'.sav'
histname=dir+'histo_zs_zp_lfZuccaAllFalse'+suff[isuff]+'.sav'

nslice = 100
z_max = 3.
pinterquart = dblarr(nslice)
pbias = dblarr(nslice)
poutlier = dblarr(nslice)
poutlier08 = dblarr(nslice)
poutlier10 = dblarr(nslice)
poutlier12 = dblarr(nslice)
binz = 1./nslice

zslice = findgen(nslice+1)*z_max/nslice

binz=0.01
zmax=3.
z=findgen(zmax/binz+1)*binz
nz = n_elements(z)
all_hist = dblarr(nz,nz)

for k=3,nslice-1 do begin
   names=dir+'cat_lfZuccaAllFalse_zOrdZP_Slice'+strtrim(k,2)+suff[isuff]+'.fits'
   nok_tot = 0

   if (mode eq 0 ) then begin
      for i=3,nslice-1 do begin 
         
         name=dir+'cat_lfZuccaAllFalse_zOrd_Slice'+strtrim(i,2)+suff[isuff]+'.fits'
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
      name=dir+'cat_lfZuccaAllFalse_zOrd_AllSlice'+suff[isuff]+'_cataLT10.fits'
      m=mrdfits(name ,1)
      ok = where(m.(6) ge zslice[k] and m.(6) lt zslice[k+1],nok)
      print,'TOTAL   with cata',nok_tot,k,i
      if (nok_tot eq 0) then mz = m[ok] else mz = [mz,m[ok]]
      nok_tot = nok_tot + nok
      print,k,i,zslice[k],nok
      
      m = 0
      if (nok_tot eq 0) then continue
      
      print,'Ready to sort the catalog'
      msz = mz[sort(mz.(3))]
      MWRFITS, msz, names, h
      print,'############# ', names,'  written'
      sxaddpar,h,'ZCAT_MIN',zslice[k]
      sxaddpar,h,'ZCAT_MAX',zslice[k+1]
      modfits,names,0,h,exten_no=1 ;Update header of ext 1, let data unchaged 
      print,'   header updated'
      
   endif

   if (mode gt 0) then msz = mrdfits(names,1)

   print,'Complete histogram and stats'
   t = (msz.(3)-msz.(6))/(1.+msz.(6))

   RSTAT,t, Med, Hinge1, Hinge2
   pbias[k] = Med
   pinterquart[k] = Hinge2-Hinge1
   out = where(abs(t) gt 0.15,nout)
   poutlier[k] = double(nout)/double(n_elements(t))*100.


   out = where(abs(t) gt 0.08,nout)
   poutlier08[k] = double(nout)/double(n_elements(t))*100.
   out = where(abs(t) gt 0.10,nout)
   poutlier10[k] = double(nout)/double(n_elements(t))*100.
   out = where(abs(t) gt 0.12,nout)
   poutlier12[k] = double(nout)/double(n_elements(t))*100.



   print,'As a function of z_s --------------------------------------------------'
   print,'-------------------------------------------------------------------------------------------------'
   print,'stat ',k,'       bias ',pbias[k],'    sigma ',pinterquart[k],'   outlier ',poutlier[k]
   print,'-------------------------------------------------------------------------------------------------'
   save,zslice,pinterquart,pbias,poutlier,poutlier08,poutlier10,poutlier12,file=statname

   if (mode le 1) then begin 
      hh=hist_2d(msz.(3),msz.(6),bin1=binz,bin2=binz,max1=z_max,max2=z_max,min1=0,min2=0)
      all_hist = all_hist+hh
      contour,all_hist,/xs,/ys,lev=1000
      save,z,all_hist,file=histname
   endif
   msz=0

endfor

END

