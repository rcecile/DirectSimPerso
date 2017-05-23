PRO prepa_hist_zs_zp_all,isuff,option,icase

; option = 0 : histo +  stat
; option 1 : stat
; icase=-1 : avec BAO, 0--9 : no BAO

if (icase lt 0) then dir='/sps/lsst/data/rcecile/Planck_BAO2/' $
else dir='/sps/lsst/data/rcecile/Planck_noBAO/' 
nslice =100

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
;for i = 4,nslice-1 do begin
for i=4,99 do begin
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   if (icase lt 0) then file='cat_AllzOrd_Slice'+strtrim(i,2)+suff[isuff]+'.fits' $
   else  file='cat_AllzOrd_'+strtrim(icase,2)+'_Slice'+strtrim(i,2)+suff[isuff]+'.fits' 
   print,file
   if(isuff le 2) then m=mrdfits(dir+file,1,h,col=['ZS','ZP']) else m=mrdfits(dir+file,1,h,col=['ZS','ZG'])  
 

   if (option eq 0) then begin
      hh=hist_2d(m.(0),m.(1),bin1=binz,bin2=binz,max1=zmax,max2=zmax,min1=0,min2=0)
      all_hist = all_hist+hh
   endif
   
   ok = where(m.(1) lt 9)
   t = (m.(0)-m.(1))/(1.+m.(1))


  if (icase lt 0) then  file='cat_AllzOrd_Slice'+strtrim(i,2)+suff[isuff]+'_cata.fits' $
  else  file='cat_AllzOrd_'+strtrim(icase,2)+'_Slice'+strtrim(i,2)+suff[isuff]+'_cata.fits'
   print,dir+file
   check = FILE_TEST(dir+file)
   if (check eq 1) then begin
      if(isuff le 2) then m=mrdfits(dir+file,1,h,col=['ZS','ZP']) else m=mrdfits(dir+file,1,h,col=['ZS','ZG'])  
      print,'MIN/MAX ZS ' ,minmax(m.(0))	
      print,'  MIN/MAX photoZ ', minmax(m.(1))
      if (n_elements(m.(0)) gt 1 and option eq 0) then begin
         hh=hist_2d(m.(0),m.(1),bin1=binz,bin2=binz,max1=zmax,max2=zmax,min1=0,min2=0)
         all_hist = all_hist+hh
      endif
   endif
   if (option eq 0) then contour,all_hist,/xs,/ys,lev=1000
   
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
   if (icase lt 0) then if (option eq 0) then histname=dir+'histo_zs_statZp'+suff[isuff]+'.sav' $
   else histname=dir+'statZp'+suff[isuff]+'.sav'
   if (icase ge 0) then if (option eq 0) then histname=dir+'histo_zs_statZp_'+strtrim(icase,2)+suff[isuff]+'.sav' $
   else histname=dir+'statZp_'+strtrim(icase,2)+suff[isuff]+'.sav'



   if (option eq 0) then save,z,pinterquart,pbias,poutlier,all_hist,file=histname $
   else save,z,pinterquart,pbias,poutlier,file=histname
endfor

zstat = dblarr(nslice)
if (option eq 0) then restore,histname
for i = 3,nslice-1 do begin
   if (icase lt 0) then  file='cat_AllzOrd_Slice'+strtrim(i,2)+'.fits' $
   else file='cat_AllzOrd_'+strtrim(icase,2)+'_Slice'+strtrim(i,2)+'.fits' 
   print,file
   hh=headfits(dir+file,ext=1)
   zmin = sxpar(hh,'ZCAT_MIN')
   zmax = sxpar(hh,'ZCAT_MAX')
   zstat[i] = (zmin+zmax)/2.
   print,zmin,zmax,zstat[i]
endfor
if (icase lt 0) then if (option eq 0) then filename=dir+'histo_zs_statZp'+suff[isuff]+'.sav' $
else filename=dir+'statZp'+suff[isuff]+'.sav'
if (icase GE 0) then if (option eq 0) then filename=dir+'histo_zs_statZp_'+strtrim(icase,2)+suff[isuff]+'.sav' $
else filename=dir+'statZp_'+strtrim(icase,2)+suff[isuff]+'.sav'

if (option eq 0) then save,zstat,z,pinterquart,pbias,poutlier,all_hist,filename $
                           else save,zstat,z,pinterquart,pbias,poutlier,filename

END
