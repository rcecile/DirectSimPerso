PRO prepa_hist_zs_zp_all_type,isuff

dir='/sps/lsst/data/rcecile/Planck_BAO2/'
nslice =100

suff = ['_errP','_errPBDT9','_errPBDT8']
pinterquart = dblarr(nslice)
pbias = dblarr(nslice)
poutlier = dblarr(nslice)

histname=dir+'histo_zs_zp_lfZuccaAllFalse_type'+suff[isuff]+'.sav'

; 2D histo
binz=0.01
zmax=3.
z=findgen(zmax/binz+1)*binz
nz = n_elements(z)
all_hist0 = dblarr(nz,nz)
all_hist1 = dblarr(nz,nz)
all_hist2 = dblarr(nz,nz)

; better to do this part in qlogin mode
for i = 3,nslice-1 do begin
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   file='cat_lfZuccaAllFalse_zOrd_Slice'+strtrim(i,2)+suff[isuff]+'.fits' 
   print,file
   m=mrdfits(dir+file,1,h,col=['ZS','ZP','TYPE'])

   ok0 = where(fix(m.(2)) eq 1,nok0)
   ok1 = where(fix(m.(2)) eq 2,nok1)
   ok2 = where(fix(m.(2)) eq 3,nok2)

   if (nok0 gt 0) then begin
      hh=hist_2d((m.(0))[ok0],(m.(1))[ok0],bin1=binz,bin2=binz,max1=zmax,max2=zmax,min1=0,min2=0)
      all_hist0 = all_hist0+hh
   endif
   
    if (nok1 gt 0) then begin
       hh=hist_2d((m.(0))[ok1],(m.(1))[ok1],bin1=binz,bin2=binz,max1=zmax,max2=zmax,min1=0,min2=0)
       all_hist1 = all_hist1+hh
    endif
    
    if (nok2 gt 0) then begin
       hh=hist_2d((m.(0))[ok2],(m.(1))[ok2],bin1=binz,bin2=binz,max1=zmax,max2=zmax,min1=0,min2=0)
       all_hist2 = all_hist2+hh
    endif
   
   file='cat_lfZuccaAllFalse_zOrd_Slice'+strtrim(i,2)+suff[isuff]+'_cata.fits'
   print,dir+file
   m=mrdfits(dir+file,1,h,col=['ZS','ZP','TYPE'])
   print,'MIN/MAX ZS ' ,minmax(m.(0))	
   print,'  MIN/MAX photoZ ', minmax(m.(1))
   if (n_elements(m.(0)) gt 1) then begin
      
      ok0 = where(fix(m.(2)) eq 1,nok0)
      ok1 = where(fix(m.(2)) eq 2,nok1)
      ok2 = where(fix(m.(2)) eq 3,nok2)
     
      if (nok0 gt 0) then begin
         hh=hist_2d((m.(0))[ok0],(m.(1))[ok0],bin1=binz,bin2=binz,max1=zmax,max2=zmax,min1=0,min2=0)
         all_hist0 = all_hist0+hh
      endif
      
      if (nok1 gt 0) then begin
         hh=hist_2d((m.(0))[ok1],(m.(1))[ok1],bin1=binz,bin2=binz,max1=zmax,max2=zmax,min1=0,min2=0)
         all_hist1 = all_hist1+hh
      endif
      
      if (nok2 gt 0) then begin
         hh=hist_2d((m.(0))[ok2],(m.(1))[ok2],bin1=binz,bin2=binz,max1=zmax,max2=zmax,min1=0,min2=0)
         all_hist2 = all_hist2+hh
      endif
   endif
   window,0,xs=400,ys=400
   contour,all_hist0,/xs,/ys,lev=1000
   window,1,xs=400,ys=400
   contour,all_hist1,/xs,/ys,lev=1000
   window,2,xs=400,ys=400
   contour,all_hist2,/xs,/ys,lev=1000

   save,z,all_hist0,all_hist1,all_hist2,file=histname
endfor


END
