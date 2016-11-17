PRO prepa_hist_zs_zp_all,isuff

dir='/sps/lsst/data/rcecile/Planck_BAO/'
dir='/sps/lsst/data/rcecile/TJP_BAO/'
nslice =70

suff = ['_errP','_errPpodds','_errPBDT','_err0.03']
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
for i = 1000,nslice-1 do begin
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   file='cat_G2_Slice'+strtrim(i,2)+suff[isuff]+'.fits'
   print,file
   if(isuff le 2) then m=mrdfits(dir+file,1,h,col=['ZS','ZP']) else m=mrdfits(dir+file,1,h,col=['ZS','ZG'])  
 

   hh=hist_2d(m.(0),m.(1),bin1=binz,bin2=binz,max1=zmax,max2=zmax,min1=0,min2=0)
   all_hist = all_hist+hh
   ok = where(m.(1) lt 9)
   t = (m.(1)-m.(0))/(1.+m.(0))


   file='cat_G2_Slice'+strtrim(i,2)+suff[isuff]+'_cata.fits'
   print,file
   if(isuff le 2) then m=mrdfits(dir+file,1,h,col=['ZS','ZP']) else m=mrdfits(dir+file,1,h,col=['ZS','ZG'])  
print,'MIN/MAX ZS ' ,minmax(m.(0))	
print,'  MIN/MAX photoZ ', minmax(m.(1))
   hh=hist_2d(m.(0),m.(1),bin1=binz,bin2=binz,max1=zmax,max2=zmax,min1=0,min2=0)
   all_hist = all_hist+hh
   contour,all_hist,/xs,/ys,lev=1000

   ok = where(m.(1) lt 9)
   t = [t,((m.(1)-m.(0))/(1.+m.(0)))[ok]]

   RSTAT,t, Med, Hinge1, Hinge2
   pbias[i] = Med
   pinterquart[i] = Hinge2-Hinge1
   out = where(abs(t) gt 0.15,nout)
   poutlier[i] = double(nout)/double(n_elements(t))*100.

   print,'-------------------------------------------------------------------------------------------------'
   print,'stat ',i,'       bias ',pbias[i],'    sigma ',pinterquart[i],'   outlier ',poutlier[i]
   print,'-------------------------------------------------------------------------------------------------'
   save,z,pinterquart,pbias,poutlier,all_hist,file=dir+'histo_zs_zp'+suff[isuff]+'.sav'
endfor

zstat = dblarr(nslice)
restore,dir+'histo_zs_zp'+suff[isuff]+'.sav'
for i = 0,nslice-1 do begin
   file='cat_G2_Slice'+strtrim(i,2)+'.fits'
   print,file
   hh=headfits(dir+file,ext=1)
   zmin = sxpar(hh,'ZCAT_MIN')
   zmax = sxpar(hh,'ZCAT_MAX')
   zstat[i] = (zmin+zmax)/2.
   print,zmin,zmax,zstat[i]
endfor
save,zstat,z,pinterquart,pbias,poutlier,all_hist,file=dir+'histo_zs_zp'+suff[isuff]+'.sav'

END
