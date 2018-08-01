PRO prepa_hist_zs_zp_all,isuff

dir='/sps/lsst/data/rcecile/Planck_BAO2/'
nslice =100
statname=dir+'statZs_lfZuccaAllFalse'+suff[isuff]+'.sav'
histname=dir+'histo_zs_zp_lfZuccaAllFalse'+suff[isuff]+'.sav'

suff = ['_errP','_errPBDT9','_errPBDT8','_err0.03']
pinterquart = dblarr(nslice)
pbias = dblarr(nslice)
poutlier = dblarr(nslice)

; 2D histo
binz=0.01
zmax=3.
z=findgen(zmax/binz+1)*binz
nz = n_elements(z)
all_hist = dblarr(nz,nz)

; better to do this part in qlogin mode
for i = 3,nslice-1 do begin
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   file='cat_lfZuccaAllFalse_zOrd_Slice'+strtrim(i,2)+suff[isuff]+'.fits' 
   print,file
   if(isuff le 2) then m=mrdfits(dir+file,1,hh,col=['ZS','ZP']) else m=mrdfits(dir+file,1,h,col=['ZS','ZG'])  
   zmin = sxpar(h,'ZCAT_MIN')
   zmax = sxpar(h,'ZCAT_MAX')
   zstat[i] = (zmin+zmax)/2.

   hh=hist_2d(m.(0),m.(1),bin1=binz,bin2=binz,max1=zmax,max2=zmax,min1=0,min2=0)
   all_hist = all_hist+hh
   
   ok = where(m.(1) lt 9)
   t = (m.(0)-m.(1))/(1.+m.(1))

   file='cat_lfZuccaAllFalse_zOrd_Slice'+strtrim(i,2)+suff[isuff]+'_cata.fits'
   print,dir+file
   check = FILE_TEST(dir+file)
   if (check eq 1) then begin
      if(isuff le 2) then m=mrdfits(dir+file,1,h,col=['ZS','ZP']) else m=mrdfits(dir+file,1,h,col=['ZS','ZG'])  
      print,'MIN/MAX ZS ' ,minmax(m.(0))	
      print,'  MIN/MAX photoZ ', minmax(m.(1))
      if (n_elements(m.(0)) gt 1) then begin
         hh=hist_2d(m.(0),m.(1),bin1=binz,bin2=binz,max1=zmax,max2=zmax,min1=0,min2=0)
         all_hist = all_hist+hh
      endif
   endif
   contour,all_hist,/xs,/ys,lev=1000
   
   ok = where(m.(1) lt 9)
   t = [t,((m.(0)-m.(1))/(1.+m.(1)))[ok]]
   
   RSTAT,t, Med, Hinge1, Hinge2
   pbias[i] = Med
   pinterquart[i] = Hinge2-Hinge1
   out = where(abs(t) gt 0.15,nout)
   poutlier[i] = double(nout)/double(n_elements(t))*100.

   print,'As a function of z_s --------------------------------------------------'
   print,'-------------------------------------------------------------------------------------------------'
   print,'stat ',i,'       bias ',pbias[i],'    sigma ',pinterquart[i],'   outlier ',poutlier[i]
   print,'-------------------------------------------------------------------------------------------------'

   save,z,all_hist,file=histname
   save,zstat,pinterquart,pbias,poutlier,file=statname
endfor


END
